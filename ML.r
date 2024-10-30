# Importing and loading necessary packages
library(tidyverse)
library(caret)
library(glmnet)
library(data.table)

# Load the data
datak <- fread('/scratch/projects/truonglab/Kallisto_pseudoaligns/Kall_filtered_tfs.csv')
datazm <- fread('/scratch/sp7512/ZhenMa/Kall_filtered_tfs_ZhenMa.csv')
datat <- fread('/scratch/sp7512/01.RawData/Kall_filtered_tfs_Tomoki.csv')

data <- rbindlist(list(datak, datazm, datat), use.names = TRUE, fill = TRUE)

# Combine SRA ids and cell identity
data[, SRA_state := paste(SRA, cell_identity, sep = "_")]

# Attach external gene name to enst id
data[, gene_isoforms := paste(external_gene_name, target_id, sep = "_")]

# Create a smaller table with just the 3 columns of interest
data_filt <- data[, .(gene_isoforms, tpm, SRA_state)]

# Pivot the smaller table to get SRA_state as column names and genes+enst ids as the columns
pivot_df <- dcast(data_filt, SRA_state ~ gene_isoforms, value.var = "tpm", fill = 0)

# Filter to include only cols where at least one value is 0
'''
this fuck is removing the sra_state column entirely...
filtered_pivot_df <- pivot_df[, lapply(.SD, function(x) if (any(x == 0)) x else NULL), .SDcols = -1]
'''

# Append the cell identity into a new column separate from the index
filtered_pivot_df[, cell_identity := sapply(strsplit(SRA_state, "_"), function(x) if (length(x) == 2) x[2] else paste(x[-1], collapse = "_"))]

# Count occurrences of each cell_identity
cell_identity_counts <- table(filtered_pivot_df$cell_identity)

# Create a boolean mask to filter rows with cell_identity occurring more than once
mask <- filtered_pivot_df$cell_identity %in% names(cell_identity_counts[cell_identity_counts > 1])

# Filter the DataFrame using the mask
filtered_pivot_df <- filtered_pivot_df[mask]

# Separate the ensmbl ids from the cell_identity column
TFs <- setdiff(names(filtered_pivot_df), c("SRA_state", "cell_identity"))

# Set predictors and predicted
X <- as.matrix(filtered_pivot_df[, ..TFs])
y <- filtered_pivot_df$cell_identity

# Train test split
set.seed(42)
trainIndex <- createDataPartition(y, p = .75, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

print('y_test:')
print(y_test)
print('y_train:')
print(y_train)

# Train logistic regression model
clf_logistic <- glmnet(X, y, family = "multinomial", alpha = 0, lambda = 0)

# Fit the logistic regression model using statsmodels to get p-values
X_with_intercept <- cbind(1, X)  # Add intercept term
model <- multinom(y ~ ., data = as.data.frame(X_with_intercept))
summary_model <- summary(model)

# Get the p-values for the coefficients
p_values <- summary_model$coefficients[, -1, drop = FALSE]

all_coefs <- list()

feature_names <- colnames(X)
class_labels <- levels(y)

# Looping through all classes to get a coefficient and p-value for each TF per class
for (index in seq_along(class_labels)) {
    class_name <- class_labels[index]
    class_coefs <- coef(model)[index, -1]
    coef_df <- data.frame(
        TF = feature_names,
        Coefficient = class_coefs,
        P_value = p_values[, index],
        cell_identity = class_name
    )
    all_coefs[[index]] <- coef_df
}

# Concatenate
all_coefs_df <- do.call(rbind, all_coefs)

# Write to CSV
fwrite(all_coefs_df, '/scratch/sp7512/thesis/TF_coefs_and_pvalues.csv')
