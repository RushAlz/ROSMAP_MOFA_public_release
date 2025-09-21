##########################
# List of genes for GSEA #
##########################

library(msigdbr)
library(purrr)

# Retrieving term lists from KEGG, Reactome and Gene Ontology

## KEGG
KEGG_msigdbr_df <- msigdbr(species = "human", category = "C2", subcategory = "KEGG")
KEGG_pathways = split(x = KEGG_msigdbr_df$gene_symbol, f = KEGG_msigdbr_df$gs_name)

## REACTOME
REACTOME_msigdbr_df <- msigdbr(species = "human", category = "C2", subcategory = "REACTOME")
REACTOME_pathways = split(x = REACTOME_msigdbr_df$gene_symbol, f = REACTOME_msigdbr_df$gs_name)

## GENE ONTOLOGY (MOLECULAR FUNCTION)
GO_MF_msigdbr_df <- msigdbr(species = "human", category = "C5", subcategory = "GO:MF")
GO_MF_terms = split(x = GO_MF_msigdbr_df$gene_symbol, f = GO_MF_msigdbr_df$gs_name)

## GENE ONTOLOGY (BILOGICAL PROCESS)
GO_BP_msigdbr_df <- msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
GO_BP_terms = split(x = GO_BP_msigdbr_df$gene_symbol, f = GO_BP_msigdbr_df$gs_name)

## GENE ONTOLOGY (CELLULAR COMPONENT)
GO_CC_tmsigdbr_df <- msigdbr(species = "human", category = "C5", subcategory = "GO:CC")
GO_CC_terms = split(x = GO_CC_tmsigdbr_df$gene_symbol, f = GO_CC_tmsigdbr_df$gs_name)

# Grouping the terms
ALL_DATA_BASE <-  unlist(list(GO_BP_terms, GO_CC_terms, GO_MF_terms, KEGG_pathways, REACTOME_pathways), recursive = FALSE)


# Retrieving term list for proteomics

GO_FEA <- unlist(list(GO_BP_terms, GO_CC_terms, GO_MF_terms), recursive = FALSE)

# Proteomics data
prot_domain = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/dlpfc_tmt.rds")

# Clearing up the names
prot_domain = rownames(prot_domain)
prot_domain = gsub('.*_','',prot_domain)

# Keep only elements of 'GO_FEA' where all values are present in 'prot_domain'
prot_domain <- keep(GO_FEA, function(x) all(x %in% prot_domain))
