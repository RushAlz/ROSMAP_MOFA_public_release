#############################################################################################
# PLOT THE RELATIONSHIP WITH FACTORS AND COVARIATES AND THE VARIANCE EXPLAINED BY PHENOTYPE #
#############################################################################################
library(MOFA2)
library(reshape2)
library(stringr)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyverse)

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

# Selecting only the AD-related factors
selected_factors = paste0("Factor",c(1, 2, 3, 4, 8, 12, 14, 26, 42))
mofa_factors = mofa_factors[,selected_factors]

# Performing linear regressions (adjusted by age of death, sex and years of education)
res_test = run_module_trait_association(data4linear_reg = mofa_factors,
                                        phenotype_dt = phenotype_dt,
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
res_pvalues = res_pvalues[,column_order]
res_pvalues = res_pvalues[gtools::mixedorder(rownames(res_pvalues), decreasing = T),]
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

# Heatmap
Heatmap(
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



###### variancePartition
library(variancePartition)

# 1) Perform variance partitioning for each variable separately (including sex, age and educ in the model)
vp_ind_pheno_r2 = list()
for(pheno_i in names(pheno_list)){
  if(pheno_list[pheno_i]=="binomial"){
    pheno_predictor = paste0("~ (1|",pheno_i,")")
  }else{
    pheno_predictor = paste0("~ ",pheno_i)
  }
  pheno_form <- as.formula(paste0(pheno_predictor,  " + (1|msex) + age_death + educ"))
  varPart_tx <- suppressWarnings( fitExtractVarPartModel(t(mofa_factors), pheno_form, phenotype_dt) )
  vp_ind_pheno_r2[[pheno_i]] = SetNames(varPart_tx[[pheno_i]],rownames(varPart_tx))
}
vp_ind_pheno_r2_df = as.data.frame(vp_ind_pheno_r2)

# 2) Perform variance partitioning for all variables together (removing variables with 0 variance explained)
vp_df = data.frame()
for(factor_i in colnames(mofa_factors)){
  non_zero_pheno = colnames(vp_ind_pheno_r2_df)[vp_ind_pheno_r2_df[factor_i,]>0]
  pheno_list_non_zero = pheno_list[non_zero_pheno]
  binary_pheno = c(names(pheno_list_non_zero)[pheno_list_non_zero=="binomial"])
  binary_form = paste0(paste(paste0("(1|",binary_pheno,")"), collapse = " + "))
  other_pheno = colnames(phenotype_dt)[!colnames(phenotype_dt)%in%binary_pheno]
  other_pheno = other_pheno[other_pheno %in% non_zero_pheno]
  formula_string = as.formula(paste0("~ ", 
                                     binary_form," + ", 
                                     paste(other_pheno, collapse = " + "),
                                     " + (1|msex) + age_death"))
  
  varPart_tx <- suppressWarnings( fitExtractVarPartModel(t(mofa_factors)[factor_i,,drop=F], formula_string, phenotype_dt, useWeights = F) )
  vp <- sortCols( varPart_tx )
  vp_df = bind_rows(vp_df, as.data.frame(vp))
}

column_annotation_VP = column_annotation
column_annotation_VP = na.omit(column_annotation_VP[colnames(vp_df),])
column_annotation_VP$color = col_pal[column_annotation_VP$category]
column_annotation_VP = column_annotation_VP %>%
  bind_rows(data.frame(category = "Sex", description = "Sex", color =  '#808080', row.names = "msex")) %>%
  bind_rows(data.frame(category = "Age at death", description = "Age at death", color = '#555555', row.names = "age_death")) %>%
  bind_rows(data.frame(category = "Residuals", description = "Residuals", color = '#AAAAAA', row.names = "Residuals"))
column_annotation_VP$variable = rownames(column_annotation_VP)

vp2plot = vp_df %>% as.data.frame() %>% rownames_to_column("factor") %>% 
  pivot_longer(cols = -factor, names_to = "variable", values_to = "Variance Explained") %>%
  mutate("Variance Explained" = `Variance Explained`) %>%
  left_join(column_annotation_VP, by = c("variable")) %>%
  mutate(category = factor(category,
                           levels = c(levels(column_annotation$category),"Age at death","Sex","Residuals")),
         factor = factor(factor, levels = gtools::mixedsort(unique(factor))))

vp2plot_color_map = unique(vp2plot[,c("color","category")])
vp2plot_color_pal = SetNames(vp2plot_color_map$color,vp2plot_color_map$category)[levels(vp2plot$category)]


ggplot(vp2plot %>%
         filter(factor %in% selected_factors), 
       aes(x = `Variance Explained`, y = factor, fill = category)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  theme_minimal(base_family = "Arial", base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.margin = margin(0, 0, 0, 0),
        plot.title = element_text(face = "plain", size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  ) +
  scale_fill_manual(values = vp2plot_color_pal) +
  coord_cartesian(xlim=c(0, 0.20)) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion(0,0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y= "Factors", fill = "", title = 'Variance Explained by Phenotype Category')



