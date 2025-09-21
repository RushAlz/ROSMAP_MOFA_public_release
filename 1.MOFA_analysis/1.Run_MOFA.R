############
# RUN MOFA #
############
library(matsbyname)
library(basilisk)
library(MOFA2)

# Retriving omics and cell type data
rna_AC = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_AC.rds")
rna_DLPF = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_DLPFC.rds")
rna_PCG = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_PCG.rds")
tmt = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/dlpfc_tmt.rds")
acetil = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/histoneAcetylation_H3K9ac.rds")
metabol = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/brain_metabolomics.rds")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/sn_RNA_norm.RDATA") #cell_types (sn_RNA_norm)


var_names <- c("rna_AC", "rna_DLPF", "rna_PCG", "sn_RNA_norm", "tmt","metabol", "acetil")

lst <- list()
for (i in var_names) {
  var_value <- get(i)
  lst[[i]] <- var_value
}

# Creating matrices
for (i in seq_along(lst)){
  matrixx <- data.matrix(lst[[i]])
  lst[[i]] <- matrixx
}

# Extracting column names from the matrices
unique_cols <- unique(unlist(lapply(lst, colnames)))

# Adding missing columns to the matrices
for (i in seq_along(lst)) {
  matrixx <- lst[[i]]
  
  missing_cols <- setdiff(unique_cols, colnames(matrixx))
  
  for (y in missing_cols) {
    matrixx <- cbind(matrixx, NA)
    colnames(matrixx)[ncol(matrixx)] <- y
  }
  
  lst[[i]] <- matrixx
}

# Rearranging column orders
for (i in seq_along(lst)){
  matrixx <- sort_rows_cols(lst[[i]])
  lst[[i]] <- matrixx
}

# Create MOFA object
mofa <- create_mofa(lst)

# Data opts
data_opts <- get_default_data_options(mofa)
data_opts$scale_views <- TRUE

# Model opts
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 50

# Training opts
train_opts <- get_default_training_options(mofa)

# Prepare MOFA
mofa <- prepare_mofa(mofa, data_options = data_opts, model_options = model_opts, training_options = train_opts)

# RUN MOFA
mofa <- run_mofa(mofa, use_basilisk = TRUE)






