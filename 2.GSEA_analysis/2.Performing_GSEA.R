##########################
# SCRIPT TO PERFORM GSEA #
##########################
library(fgsea)
library(MOFA2)
library(foreach)
library(reshape2)
library(matsbyname)
library(dplyr)


# Retriving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# Retriving enrichment lists
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/ALL_DATA_BASE.RDATA")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/new_genes_list/MODULES_LIST.RDATA")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/prot_domain.RDATA")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/new_metabol_list.RDATA")

# Excluding non relevant terms from the enrichment lists
ALL_DATA_BASE <- ALL_DATA_BASE[!grepl("^KEGG", names(ALL_DATA_BASE))]

ALL_DATA_BASE <- ALL_DATA_BASE[!names(ALL_DATA_BASE) %in% c('cognitive decline', 'relevant_metabol', 'np diagnosis', 'clinical diagnosis', 
                                                            'cognition', 'nia-reagan score', 'tangles', 'global np', 'proved_AD', 'amyloid',
                                                            'unknown','MODULO_109')]


# Appending enrichment lists
ALL_DATA_BASE <- append(ALL_DATA_BASE, MODULES_LIST, after = length(ALL_DATA_BASE))
prot_domain <- append(prot_domain, MODULES_LIST, after = length(prot_domain))
#################################################################################################################################################

# Func. to perform GSEA
run_GSEA <- function(my_list_of_genes, pathways_list = GO_FEA, threads = 2){
  
  if(!is.list(my_list_of_genes)){
    my_list_of_genes = list(GeneList = my_list_of_genes)
  }
  enrich_res_df = foreach(i = 1:length(my_list_of_genes), .combine=rbind) %dopar% {
    GeneSet = names(my_list_of_genes)[i]
    fgsea.ovca = fgseaMultilevel(pathways = pathways_list, 
                                 stats = my_list_of_genes[[GeneSet]], 
                                 nproc = threads, 
                                 minSize = 10, 
                                 maxSize = Inf, 
                                 nPermSimple = 10000,
                                 eps = 0)
    fgsea.ovca$leadingEdge = NULL
    fgsea.ovca$GeneSet = GeneSet
    fgsea.ovca
  }
  enrich_res_df = enrich_res_df[order(enrich_res_df$padj),]
  return(enrich_res_df)  
}

# Enrichment functions, tailored for each data type
rna_GSEA_enrichment <- function(Views){
  
  # Selecting only the AD-related factors
  facs <- c(1,2,3,4,8,12,14,26,42)
  
  factor_list <- list()
  
  # Run GSEA for every factor
  for (i in facs){
    
    # Grabbing all features weights
    enrich = plot_top_weights(mofa, view = Views, factor = i, nfeatures = nrow(mofa@data[[Views]]$group1))
    
    # Changing the genes names, for better readability. 'r' is now a vector with the new genes names
    df <- enrich$data[c('feature_id', 'value', 'sign')]
    r <- df
    r<- gsub(".*_(.*?)(\\..)?_.*_.*", "\\1", df$feature_id)
    r = unique(r)
    
    # Adding the sign for each loading. 'final_values' is now a named vector with features values
    values <- c()
    for (t in 1:length(r)) {
      val <- df$value[t]
      sign <- df$sign[t]
      if (sign == '-') {
        val <- -val}
      values <- append(values, val)
    }
    final_values <- setNames(values, r)
    
    # Performing GSEA - 'final_values' is now a data.table with the enrichment information
    final_values <- run_GSEA(final_values, ALL_DATA_BASE)
    
    # Adding columns (fator, var_explained, view, mofa)
    final_values$factor <- i
    factr_name <- paste0('Factor',i)
    all_var <- plot_variance_explained(mofa, max_r2=15)$data
    var_exp <- subset(all_var, factor == factr_name & view == Views)$value
    final_values$var_exp <- var_exp
    final_values$view <- Views
    final_values$mofa <- 'standart_50f'
    
    factor_list[[paste0('factor',i)]] <- final_values
  }
  
  final_df <- do.call(rbind, factor_list)
  
  if(Views == 'rna_AC'){
    
    assign('AC_table', final_df, envir = .GlobalEnv)
  }
  else if(Views == 'rna_PCG'){
    
    assign('PCG_table', final_df, envir = .GlobalEnv)
  }
  else{
    
    assign('DLPF_table', final_df, envir = .GlobalEnv)
  }
  
}
acetil_methyl_enrichment <- function(Views){
  
  # Selecting only the AD-related factors
  facs <- c(1,2,3,4,8,12,14,26,42)
  
  factor_list <- list()
  
  # Dowloading metadata for histone acetylation or DNA methylation
  if (Views == 'acetil') {
    metadata <- readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/histoneAcetylation_H3K9ac.metadata.rds")
  } else {
    metadata <- readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/DNAmethylation.metadata.rds")
  }
  
  # Run GSEA for each factor
  for (i in facs){
    
    # Grabbing all features weights
    enrich = plot_top_weights(mofa, view = Views, factor = i, nfeatures = nrow(mofa@data[[Views]]$group1))
    df <- enrich$data[c('feature_id', 'value', 'sign')]
    
    # Preparing 'metadata' to merge with 'df'. The goal is to map 'df' peaks to the closest genes in 'metadata'
    temp_metadata <- metadata[,c('probeid','closest_gene')]
    temp_metadata <- subset(temp_metadata, probeid %in% df$feature_id)
    colnames(temp_metadata)[colnames(temp_metadata) == "probeid"] <- "feature_id"
    
    #Merging
    merge_metadata <- merge(temp_metadata, df, by = 'feature_id')
    merge_metadata <-merge_metadata[,-1]
    
    # Some peaks share the same closest gene; so here we are keeping only one entry per gene (the one with higher absolute loading)
    for (j in unique(merge_metadata$closest_gene)){
      
      same_names = subset(merge_metadata, closest_gene == j)
      if (any(same_names$sign == '+')){
        
        max_value <- unique(max(same_names$value))
        exclude_names <- subset(same_names, value !=max_value)
        merge_metadata <- anti_join(merge_metadata, exclude_names, by = c("closest_gene", 'value', 'sign'))
      }
      
      else{
        
        max_value <- unique(min(same_names$value))
        exclude_names <- subset(same_names, value !=max_value)
        merge_metadata <- anti_join(merge_metadata, exclude_names, by = c("closest_gene", 'value', 'sign'))
      }
    }
    
    # Adding the sign for each loading 'final_values' is now a named vector with features values
    values <- c()
    for (t in 1:length(merge_metadata$closest_gene)) {
      val <- merge_metadata$value[t]
      sign <- merge_metadata$sign[t]
      if (sign == '-') {
        val <- -val}
      values <- append(values, val)
    }
    final_values <- setNames(values, merge_metadata$closest_gene)
    
    # Performing GSEA - 'final_values' is now a data.table with the enrichment information
    final_values <- run_GSEA(final_values, ALL_DATA_BASE)
    
    # Adding columns (fator, var_explained, view, mofa)
    final_values$factor <- i
    factr_name <- paste0('Factor',i)
    all_var <- plot_variance_explained(mofa, max_r2=15)$data
    var_exp <- subset(all_var, factor == factr_name & view == Views)$value
    final_values$var_exp <- var_exp
    
    final_values$view <- Views
    final_values$mofa <- 'standart_50f'
    
    factor_list[[paste0('factor',i)]] <- final_values
    
  }
  
  final_df <- do.call(rbind, factor_list)
  
  if(Views == 'metil'){
    
    assign('metil_table', final_df, envir = .GlobalEnv)
  }
  
  else{
    
    assign('acetyl_table', final_df, envir = .GlobalEnv)
  }
}
prot_enrichment <- function(){
  
  # Selecting only the AD-related factors
  facs <- c(1,2,3,4,8,12,14,26,42)
  factor_list <- list()
  
  # Run GSEA for each factor
  for (i in facs){
    
    # Grabbing all features weights and cleaning up the names
    PROT_enrichment = plot_top_weights(mofa, view = "tmt", factor = i, nfeatures = nrow(mofa@data$tmt$group1))
    PROT_enrichment <- PROT_enrichment$data[c('feature_id', 'value', 'sign')]
    prot_names <- gsub('.*_','', PROT_enrichment$feature_id)
    prot_names <- unique(prot_names)
    
    # Adding the sign for each loading 'PROT_enrichment' is now a named vector with features values
    prot_values <- c()
    for (j in 1:length(prot_names)) {
      val <- PROT_enrichment$value[j]
      sign <- PROT_enrichment$sign[j]
      if (sign == '-') {
        val <- -val}
      prot_values <- append(prot_values, val)
    }
    PROT_enrichment <-  setNames(prot_values, prot_names)
    
    # Performing GSEA - 'final_values' is now a data.table with the enrichment information
    final_values = run_GSEA(PROT_enrichment, prot_domain)
    
    
    # Adding columns (fator, var_explained, view, mofa)
    final_values$factor <- i
    factr_name <- paste0('Factor',i)
    all_var <- plot_variance_explained(mofa, max_r2=15)$data
    var_exp <- subset(all_var, factor == factr_name & view == 'tmt')$value
    final_values$var_exp <- var_exp
    
    final_values$view <- 'tmt'
    final_values$mofa <- 'standart_50f'
    
    factor_list[[paste0('factor',i)]] <- final_values
    
  }
  
  final_df <- do.call(rbind, factor_list)
  assign('tmt_table', final_df, envir = .GlobalEnv)
}
metabol_enrichment <- function(){
  
  # Selecting only the AD-related factors
  facs <- c(1,2,3,4,8,12,14,26,42)
  factor_list <- list()
  
  # Run GSEA for every factor
  for (i in facs){
    
    # Grabbing all features weights
    enrich = plot_top_weights(mofa, view = 'metabol', factor = as.numeric(i), nfeatures = nrow(mofa@data$metabol$group1))
    
    # Changing the metabolites names, for better readability. 'r' is now a vector with the new metabolites names
    df <- enrich$data[c('feature_id', 'value', 'sign')]
    df$feature_id <- gsub('^X', "", df$feature_id)
    r <- df$feature_id
    
    # Adding the sign for each loading. 'final_values' is now a named vector with features values
    values <- c()
    for (t in 1:length(r)) {
      val <- df$value[t]
      sign <- df$sign[t]
      if (sign == '-') {
        val <- -val}
      values <- append(values, val)
    }
    
    final_values <- setNames(values, r)
    
    # Performing GSEA - 'final_values' is now a data.table with the enrichment information
    final_values <- run_GSEA(final_values, new_metabol_list)
    
    # Adding columns (fator, var_explained, view, mofa)
    final_values$factor <- i
    factr_name <- paste0('Factor',i)
    all_var <- plot_variance_explained(mofa, max_r2=15)$data
    var_exp <- subset(all_var, factor == factr_name & view == 'metabol')$value
    final_values$var_exp <- var_exp
    
    final_values$view <- 'metabol'
    final_values$mofa <- 'standart_50f'
    
    factor_list[[paste0('factor',i)]] <- final_values
  }
  
  final_df <- do.call(rbind, factor_list)
  
  assign('metabol_table', final_df, envir = .GlobalEnv)
}

# Executing functions
rna_GSEA_enrichment(Views = 'rna_AC')
rna_GSEA_enrichment(Views = 'rna_PCG')
rna_GSEA_enrichment(Views = 'rna_DLPF')
acetil_methyl_enrichment(Views = 'acetil')
prot_enrichment()
metabol_enrichment()

# Binding the results
final_table_v3 <- rbind(AC_table,PCG_table,DLPF_table,acetyl_table,tmt_table, metabol_table)
