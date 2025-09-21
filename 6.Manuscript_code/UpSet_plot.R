################################
# Upsetplot of multiomics data #
################################
library(MOFA2)

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# Retriving omics and cell type data
rna_AC = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_AC.rds")
rna_DLPF = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_DLPFC.rds")
rna_PCG = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_PCG.rds")
tmt = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/dlpfc_tmt.rds")
acetil = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/histoneAcetylation_H3K9ac.rds")
metabol = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/brain_metabolomics.rds")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/sn_RNA_norm.RDATA")

var_names <- c("rna_AC", "rna_DLPF", "rna_PCG", "sn_RNA_norm", "tmt","metabol", "acetil")

# Retrieving participants IDs
IDs_vector <- c()
for(var in var_names){
  var_matrix_colnames <- colnames(get(var))
  IDs_vector <- c(IDs_vector,var_matrix_colnames)
}
IDs_vector <- unique(IDs_vector)

# Creating presence data.frame
presence_df <- data.frame(
  ID = IDs_vector,
  H3K9ac = as.numeric(IDs_vector %in% colnames(acetil)),
  `RNA (AC)` = as.numeric(IDs_vector %in% colnames(rna_AC)),
  `RNA (PCG)` = as.numeric(IDs_vector %in% colnames(rna_PCG)),
  `RNA (DLPFC)` = as.numeric(IDs_vector %in% colnames(rna_DLPF)),
  Proteomics = as.numeric(IDs_vector %in% colnames(tmt)),
  Metabolomics = as.numeric(IDs_vector %in% colnames(metabol)),
  `Cell Types` = as.numeric(IDs_vector %in% colnames(sn_RNA_norm)),
  check.names = FALSE
)
rownames(presence_df) <- presence_df$ID
presence_df$ID <- NULL

library(ComplexUpset)
library(ggplot2)
# Upsetplot
xx = ComplexUpset::upset(
  presence_df,
  colnames(presence_df),
  name = 'Samples',
  min_size = 1,
  matrix = 
    (
      intersection_matrix(geom=geom_point(shape='circle filled', size=8)) 
      + scale_color_manual(
        values = c('H3K9ac' = "#00468BFF",
                   'RNA (AC)' =  "#ED0000FF",
                   'RNA (DLPFC)' = "#42B540FF",
                   'RNA (PCG)' = "#0099B4FF",
                   'Proteomics' = "#925E9FFF",
                   'Metabolomics' =  "#FDAF91FF",
                   'Cell Types' = "#AD002AFF"),
        guide=guide_legend(override.aes=list(shape='circle'))
      )
    ),
  queries = list
  (
    upset_query(set='H3K9ac', fill = "#00468BFF"),
    upset_query(set='RNA (AC)', fill = "#ED0000FF"),
    upset_query(set='RNA (DLPFC)', fill = "#42B540FF"),
    upset_query(set='RNA (PCG)', fill = "#0099B4FF"),
    upset_query(set='Proteomics', fill = "#925E9FFF"),
    upset_query(set='Metabolomics', fill = "#FDAF91FF"),
    upset_query(set='Cell Types', fill = "#AD002AFF")
  ), 
  base_annotations = list
  (
    'Intersection size' = intersection_size(
      counts = TRUE,
      mapping = aes(fill='bars_color'),
      text = list(size = 12) 
    ) + scale_fill_manual(values=c('bars_color'='#DC9445FF'), guide='none') +
      theme(
        axis.text.y = element_text(size=26), 
        axis.title.y = element_text(size = 36, face = 'plain') 
      )
  ),
  set_sizes = (
    upset_set_size() +
      theme(
        axis.text.x = element_text(size = 26), 
        axis.title.x = element_text(size = 36) 
      )
  )
) + theme(
  text = element_text(size = 36), 
  axis.text.y = element_text(size = 32)
)
