####################################################
# Plotting  the correlation matrix between factors #
####################################################
library(MOFA2)
library(ggplot2)
library(dplyr)
library(circlize)
library(reshape2)

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# Retrieving participants phenotypes data
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/meta.RDATA")
meta <- meta[(meta$projid %in% mofa@samples_metadata$sample) , ]
samples_metadata(mofa) <- meta

# MOFA function to plot the correlation matrix between factors
plot_factor_cor(mofa)


