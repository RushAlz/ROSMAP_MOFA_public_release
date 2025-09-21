############################
# List of modules for GSEA #
############################

# AD GWAS
bellengez <- read.delim("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/new_genes_list/gene_lists/AD_GWAS/gwas_ad_bellengez2022.txt")
AD_GWAS <- as.list(bellengez)
names(AD_GWAS) <- 'AD_GWAS'

# DAM (Damage-associated microglia)
DAM_files <- list.files(path = "C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/new_genes_list/gene_lists/DAM_lists", pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
DAM_files <- DAM_files[-3]

DAM_files <- lapply(DAM_files, function(x) read.delim(x))

DAM_down <- DAM_files[[1]]$Human.gene.name
DAM_down <- DAM_down[DAM_down != ""]

DAM_up <- DAM_files[[2]]$Human_gene_name
DAM_up <- DAM_up[DAM_up != ""]

DAM <- list()
DAM[['DAM_down']] <- DAM_down
DAM[['DAM_up']] <- DAM_up

# HAM (human AD microglia)
HAM <- list.files(path = "C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/new_genes_list/gene_lists/HAM_lists", pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
HAM <- lapply(HAM, function(x) read.delim(x))
HAM <- lapply(HAM, function(x) unlist(x, use.names = FALSE))
names(HAM) <- c('HAM_down','HAM_up')

# Modules proteomics
modules_johnson_etal <- list.files(path = "C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/new_genes_list/gene_lists/modules_johnson_etal", pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
modules_johnson <- lapply(modules_johnson_etal, function(x) read.delim(x))
modules_johnson <- lapply(modules_johnson, function(x) unlist(x, use.names = FALSE))
names(modules_johnson) <- sub(".*/(.*?)\\.txt$", "\\1",modules_johnson_etal)

modules_jhonson_exclude <- c("m10_ambiguous", "m14_protein_folding", "m15_ambiguous", "m16_rna_binding", "m17_transcription", "m19_axonogenesis",
                             "m21_mhc_complex_immune", "m23_ambiguous", "m26_complement_acute_phase", "m27_extracellular_matrix", "m28_ribosome_translation",
                             "m30_proteasome", "m31_axon_node_ion_channel", "m32_ambiguous", "m33_ambiguous", "m34_ambiguous", "m35_ambiguous", 
                             "m36_neurotransmitter_regulation", "m37_endosome", "m38_heat_shock_folding", "m39_translation_initiation", "m40_ambiguous", "m41_ambiguous",
                             "m43_ribonucleoprotein_binding", "m44_ribosome_translation", "m6_ribosome", "m9_golgi")

modules_johnson <- modules_johnson[!names(modules_johnson) %in% modules_jhonson_exclude]

# Modules transcriptomics
modules_mostafavi_2018 <- list.files(path = "C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/new_genes_list/gene_lists/modules_mostafavi_2018", pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
modules_mostafavi <- lapply(modules_mostafavi_2018, function(x) read.delim(x))
modules_mostafavi <- lapply(modules_mostafavi, function(x) unlist(x, use.names = FALSE))
names(modules_mostafavi) <- sub(".*/(.*?)\\.txt$", "\\1",modules_mostafavi_2018)

mostafavi_to_exclude <- c("Sara_M1", "Sara_M10", "Sara_M106", "Sara_M107", "Sara_M108", "Sara_M11", "Sara_M110", "Sara_M112", "Sara_M113", "Sara_M115",
                          "Sara_M116", "Sara_M117", "Sara_M118", "Sara_M119", "Sara_M121", "Sara_M122", "Sara_M123", "Sara_M125", "Sara_M126", "Sara_M127",
                          "Sara_M128", "Sara_M16", "Sara_M17", "Sara_M18", "Sara_M187", "Sara_M19", "Sara_M2", "Sara_M21", "Sara_M22", "Sara_M23", "Sara_M233",
                          "Sara_M234", "Sara_M257", "Sara_M3", "Sara_M4", "Sara_M5", "Sara_M6", "Sara_M8", "Sara_M9")

modules_mostafavi <- modules_mostafavi[!names(modules_mostafavi) %in% mostafavi_to_exclude]


# PIG (plaque-induced genes)
PIG <- read.delim("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/new_genes_list/gene_lists/PIG_lists/PIG_orthologs.txt")

PIG <- as.list(PIG)
PIG <- PIG[-2]
names(PIG) <- 'PIG'

# Cell types
library(openxlsx)
cell_types <- read.xlsx('C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/table_2_cell_type_DEGs_v2.xlsx', sheet = 2)

sub_cell_types_list <- list()
for(states in unique(cell_types$state)){
  filtering_genes <- filter(cell_types, state == states & p_val_adj < 0.05 & p_val_adj != 0)$gene
  sub_cell_types_list[[states]] <- filtering_genes
}
sub_cell_types_list <- sub_cell_types_list[!is.na(names(sub_cell_types_list))]

#####################################################################################################################################
all_lists <- list(AD_GWAS = AD_GWAS, HAM = HAM, DAM = DAM, modules_johnson = modules_johnson,  modules_mostafavi = modules_mostafavi, 
                  PIG = PIG,sub_cell_types_list = sub_cell_types_list)

MODULES_LIST <- list()
for (list_name in names(all_lists)) {
  for (vec_name in names(all_lists[[list_name]])) {
    MODULES_LIST[[vec_name]] <- all_lists[[list_name]][[vec_name]]
  }
}
