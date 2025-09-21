library(MOFA2)
library(reshape2)
library(stringr)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyverse)

# Function to add asterisks to significant terms 
my_cell_fun = function(j, i, x, y, width, height, fill, table_1) {
  if(!is.na(table_1[i, j]) & abs(table_1[i, j]) > -log10(0.001)) {
    grid.text("***", x, y, gp = gpar(fontsize = 12), vjust = 0.77)
  } else if(!is.na(table_1[i, j]) & abs(table_1[i, j]) > -log10(0.01)) {
    grid.text("**", x, y, gp = gpar(fontsize = 12), vjust = 0.77)
  } else if(!is.na(table_1[i, j]) & abs(table_1[i, j]) > -log10(0.05)) {
    grid.text("*", x, y, gp = gpar(fontsize = 12), vjust = 0.77)
  }
}

# Retriving omics and cell type data
rna_AC = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_AC.rds")
rna_DLPF = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_DLPFC.rds")
rna_PCG = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_PCG.rds")
tmt = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/dlpfc_tmt.rds")
acetil = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/histoneAcetylation_H3K9ac.rds")
metabol = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/brain_metabolomics.rds")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/sn_RNA_norm.RDATA") #cell_types (sn_RNA_norm)

# Utils
source("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Figuras_MOFA/Scripts/covars/util_covar.R")

# Retrieving participants phenotypes data
phenotypes = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/basic_Apr2022_selected_list_Jul2024.rds")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/pheno_list_Jul2024.RData")

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# Checking if samples are in the same order
phenotype_dt = phenotypes[match(mofa@samples_metadata$sample, rownames(phenotypes)), ]
identical(rownames(phenotype_dt), mofa@samples_metadata$sample)
samples_metadata(mofa) <- samples_metadata(mofa) %>% left_join(phenotype_dt %>% rownames_to_column("projid"),by = c("sample"="projid"))
rownames(samples_metadata(mofa)) = samples_metadata(mofa)$sample

# MOFAs Z matrix
mofa_factors = mofa@expectations$Z$group1 %>% as.data.frame()

####################################################################################################################################################
all_samples <- unique(c(colnames(rna_AC),colnames(rna_DLPF),colnames(rna_PCG),colnames(acetil),colnames(tmt),colnames(metabol),colnames(sn_RNA_norm)))

prot_samples <- intersect(colnames(tmt),all_samples)

except_prot_samples <- setdiff(all_samples,prot_samples)

set.seed(123)
subset_except_prot <- sample(except_prot_samples, length(prot_samples))

phenotype_prot <- phenotype_dt[prot_samples,]
phenotype_except_prot <- phenotype_dt[subset_except_prot,]

selected_factors = paste0("Factor",c(14))
mofa_F14 = mofa_factors[,selected_factors, drop = FALSE]

mofa_prot <- mofa_F14[prot_samples,,drop=FALSE]
mofa_except_prot <- mofa_F14[subset_except_prot,,drop = FALSE]
################################################################################################################################################

# Performing linear regressions (adjusted by age of death, sex and years of education)
res_test = run_module_trait_association(data4linear_reg = mofa_prot,
                                        phenotype_dt = phenotype_prot,
                                        pheno_list = pheno_list,
                                        covariates = c("age_death","msex","educ"),
                                        verbose = F)

res_beta <- res_test$matrix_beta

# Adjusting P-values by all phenotypes together
matrix_pvalue = res_test$matrix_pvalue
adj_matrix_pvalue = matrix(p.adjust(as.vector(as.matrix(matrix_pvalue)), method='fdr'),ncol=ncol(matrix_pvalue))
dimnames(adj_matrix_pvalue) = dimnames(matrix_pvalue)

# Converting adjusted p-values to signed -log10 scale and transposing the matrix
res_pvalues <- -log10(adj_matrix_pvalue)
res_pvalues <- res_pvalues * sign(res_beta)
res_pvalues <- t(res_pvalues)


# Adjusting heatmap colors
abs_max_logP = max(abs(res_pvalues), na.rm = T)
coll <- colorRamp2(quantile(c(abs_max_logP, -abs_max_logP)), 
                   c("#0571B0", "#92C5DE", "#F7F7F7", "#F4A582", "#CA0020"))

# Creating column annotations for heatmap
column_annotation = data.frame(pheno_label = colnames(res_pvalues))
column_annotation %>% left_join(pheno_data %>% rownames_to_column("pheno_label")) %>%
  dplyr::select(category, pheno_label, description) %>%
  column_to_rownames("pheno_label") -> column_annotation

col_pal = RColorBrewer::brewer.pal(n = length(unique(column_annotation$category)), name ="Dark2")
names(col_pal) = sort(unique(column_annotation$category))
col_pal = c(Cognition="#66A61E",
            Pathology="#E6AB02",
            `Vascular Pathology`="#A6761D",
            Motor="#7570B3",
            Disabilities="#1B9E77",
            Depression="#D95F02",
            Genetics="#E7298A")


column_order = order(column_annotation$description)
column_annotation = column_annotation[column_order,,drop=F]
col_ha = columnAnnotation(`Category` = column_annotation$category,
                          col = list(`Category` = col_pal),
                          height = unit(0.3,"mm"),
                          annotation_height = unit(0.3, "mm"),
                          simple_anno_size = unit(2, "mm"),
                          annotation_label = "  "
)

# Rearanging column and row orders for 'res_pvalues'
res_pvalues = res_pvalues[,column_order, drop = FALSE]
res_pvalues = res_pvalues[gtools::mixedorder(rownames(res_pvalues), decreasing = T),,drop = FALSE]
colnames(res_pvalues) = column_annotation$description

col_split = SetNames(column_annotation$category, column_annotation$category)

# Abbreviating Depression to Depres.
col_split_char <- as.character(col_split)
col_split_char <- sub("Depression", "Depres.", col_split_char)
col_split <- factor(col_split_char, levels = c("Cognition","Pathology","Vascular Pathology", "Motor","Disabilities", "Depres.","Genetics"))

# Renaming some covariates
colnames(res_pvalues)[colnames(res_pvalues) == "AD NIA-Reagan score"] <- "NIA-Reagan score"
colnames(res_pvalues)[colnames(res_pvalues) == "Depressive symptoms (CES-D)"] <- "Depressive symptoms (mCES-D)"
colnames(res_pvalues)[colnames(res_pvalues) == "Basic activities of daily living (ADL) "] <- "Basic activities of daily living"
colnames(res_pvalues)[colnames(res_pvalues) == "Instrumental activities of daily living (IADL)"] <- "Instrumental activities of daily living"
colnames(res_pvalues)[colnames(res_pvalues) == 'Cognitive decline'] <- 'Cognitive decline (slope)'
colnames(res_pvalues)[colnames(res_pvalues) == 'Cognitive resilience'] <- 'Cognitive resilience (slope)'
colnames(res_pvalues)[colnames(res_pvalues) == "Alzheimer's Dementia"] <- "Alzheimer's dementia"
colnames(res_pvalues)[colnames(res_pvalues) == "Beta-Amyloid"] <- "Beta-amyloid"
colnames(res_pvalues)[colnames(res_pvalues) == "Lewy Body disease"] <- "Lewy body disease"
colnames(res_pvalues)[colnames(res_pvalues) == "Cerebral Atherosclerosis"] <- "Cerebral atherosclerosis"
colnames(res_pvalues)[colnames(res_pvalues) == "Cerebral Infarctions (Gross)"] <- "Cerebral infarctions (Gross)"
colnames(res_pvalues)[colnames(res_pvalues) == "Cerebral Infarctions (Micro)"] <- "Cerebral infarctions (Micro)"
colnames(res_pvalues)[colnames(res_pvalues) == "Major Depressive Disorder"] <- "Major depressive disorder"

ht_opt$TITLE_PADDING = unit(c(4,4), "points")

rownames(res_pvalues) <- 'Factor 14 \nn = 596'

# Heatmap
prot_heat <- Heatmap(
  res_pvalues, 
  cell_fun = function(j, i, x, y, width, height, fill) {
    my_cell_fun(j, i, x, y, width, height, fill, res_pvalues) 
  },
  heatmap_legend_param = list(
    col_fun = coll,
    title = bquote(-log[10](adj.~italic(P))~(signed)),
    title_position = "leftcenter-rot",
    legend_height = unit(4, "cm") ),
  #column_title = "Association of Factors and Phenotypes",
  col = coll,
  column_split = col_split,
  column_title_gp = gpar(
    fill = SetNames(rep("white",length(column_annotation$category)), column_annotation$category),
    alpha = 1,
    fontsize = 8
  ),
  row_names_side = "left", show_row_names = T,
  cluster_rows = F, cluster_columns = F,
  column_names_rot = 55,
  column_names_gp = gpar(fontsize = 8.5),
  row_names_gp = gpar(fontsize = 9),
  bottom_annotation = col_ha,
  show_row_dend = F, show_column_dend = F, 
  rect_gp = gpar(col = "black", lwd = 1))
###################################################################################################################################################

# Performing linear regressions (adjusted by age of death, sex and years of education)
res_test = run_module_trait_association(data4linear_reg = mofa_except_prot,
                                        phenotype_dt = phenotype_except_prot,
                                        pheno_list = pheno_list,
                                        covariates = c("age_death","msex","educ"),
                                        verbose = F)

res_beta <- res_test$matrix_beta

# Adjusting P-values by all phenotypes together
matrix_pvalue = res_test$matrix_pvalue
adj_matrix_pvalue = matrix(p.adjust(as.vector(as.matrix(matrix_pvalue)), method='fdr'),ncol=ncol(matrix_pvalue))
dimnames(adj_matrix_pvalue) = dimnames(matrix_pvalue)

# Converting adjusted p-values to signed -log10 scale and transposing the matrix
res_pvalues <- -log10(adj_matrix_pvalue)
res_pvalues <- res_pvalues * sign(res_beta)
res_pvalues <- t(res_pvalues)

# Creating column annotations for heatmap
column_annotation = data.frame(pheno_label = colnames(res_pvalues))
column_annotation %>% left_join(pheno_data %>% rownames_to_column("pheno_label")) %>%
  dplyr::select(category, pheno_label, description) %>%
  column_to_rownames("pheno_label") -> column_annotation

col_pal = RColorBrewer::brewer.pal(n = length(unique(column_annotation$category)), name ="Dark2")
names(col_pal) = sort(unique(column_annotation$category))
col_pal = c(Cognition="#66A61E",
            Pathology="#E6AB02",
            `Vascular Pathology`="#A6761D",
            Motor="#7570B3",
            Disabilities="#1B9E77",
            Depression="#D95F02",
            Genetics="#E7298A")


column_order = order(column_annotation$description)
column_annotation = column_annotation[column_order,,drop=F]
col_ha = columnAnnotation(`Category` = column_annotation$category,
                          col = list(`Category` = col_pal),
                          height = unit(0.3,"mm"),
                          annotation_height = unit(0.3, "mm"),
                          simple_anno_size = unit(2, "mm"),
                          annotation_label = "  "
)

# Rearanging column and row orders for 'res_pvalues'
res_pvalues = res_pvalues[,column_order, drop = FALSE]
res_pvalues = res_pvalues[gtools::mixedorder(rownames(res_pvalues), decreasing = T),,drop = FALSE]
colnames(res_pvalues) = column_annotation$description

col_split = SetNames(column_annotation$category, column_annotation$category)

# Abbreviating Depression to Depres.
col_split_char <- as.character(col_split)
col_split_char <- sub("Depression", "Depres.", col_split_char)
col_split <- factor(col_split_char, levels = c("Cognition","Pathology","Vascular Pathology", "Motor","Disabilities", "Depres.","Genetics"))

# Renaming some covariates
colnames(res_pvalues)[colnames(res_pvalues) == "AD NIA-Reagan score"] <- "NIA-Reagan score"
colnames(res_pvalues)[colnames(res_pvalues) == "Depressive symptoms (CES-D)"] <- "Depressive symptoms (mCES-D)"
colnames(res_pvalues)[colnames(res_pvalues) == "Basic activities of daily living (ADL) "] <- "Basic activities of daily living"
colnames(res_pvalues)[colnames(res_pvalues) == "Instrumental activities of daily living (IADL)"] <- "Instrumental activities of daily living"
colnames(res_pvalues)[colnames(res_pvalues) == 'Cognitive decline'] <- 'Cognitive decline (slope)'
colnames(res_pvalues)[colnames(res_pvalues) == 'Cognitive resilience'] <- 'Cognitive resilience (slope)'
colnames(res_pvalues)[colnames(res_pvalues) == "Alzheimer's Dementia"] <- "Alzheimer's dementia"
colnames(res_pvalues)[colnames(res_pvalues) == "Beta-Amyloid"] <- "Beta-amyloid"
colnames(res_pvalues)[colnames(res_pvalues) == "Lewy Body disease"] <- "Lewy body disease"
colnames(res_pvalues)[colnames(res_pvalues) == "Cerebral Atherosclerosis"] <- "Cerebral atherosclerosis"
colnames(res_pvalues)[colnames(res_pvalues) == "Cerebral Infarctions (Gross)"] <- "Cerebral infarctions (Gross)"
colnames(res_pvalues)[colnames(res_pvalues) == "Cerebral Infarctions (Micro)"] <- "Cerebral infarctions (Micro)"
colnames(res_pvalues)[colnames(res_pvalues) == "Major Depressive Disorder"] <- "Major depressive disorder"

ht_opt$TITLE_PADDING = unit(c(4,4), "points")

rownames(res_pvalues) <- 'Factor 14 \nn = 596'

# Heatmap
except_heat <- Heatmap(
  res_pvalues, 
  cell_fun = function(j, i, x, y, width, height, fill) {
    my_cell_fun(j, i, x, y, width, height, fill, res_pvalues) 
  },
  heatmap_legend_param = list(
    col_fun = coll,
    title = bquote(-log[10](adj.~italic(P))~(signed)),
    title_position = "leftcenter-rot",
    legend_height = unit(4, "cm") ),
  col = coll,
  column_split = col_split,
  column_title_gp = gpar(
    fill = SetNames(rep("white",length(column_annotation$category)), column_annotation$category),
    alpha = 1,
    fontsize = 8
  ),
  row_names_side = "left", show_row_names = T,
  cluster_rows = F, cluster_columns = F,
  column_names_rot = 55,
  column_names_gp = gpar(fontsize = 8.5),
  row_names_gp = gpar(fontsize = 9),
  bottom_annotation = col_ha,
  show_row_dend = F, show_column_dend = F, 
  rect_gp = gpar(col = "black", lwd = 1))

prot_heat %v% except_heat


pdf(file = "C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Figuras_MOFA/pre_inkscape/geral/covars_with_and_without_proteomics.pdf",
    width =11.5,
    height = 4)
prot_heat %v% except_heat
dev.off()
