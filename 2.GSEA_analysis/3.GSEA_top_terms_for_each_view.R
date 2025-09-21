#############################################################################################
# SCRIPT TO GENERATE A HEATMAP WITH THE TOP ENRICHED TERMS (BASED ON P-VALUE) FOR EACH VIEW #
#############################################################################################
library(MOFA2)
library(reshape2)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggsci)

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")


##################################################################################################################################
# Retrieving GSEA data
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Best_mofas_analysis/standart_50f/final_table_v3.RDATA")

# Renaming protein modules
final_table_v3$pathway <- gsub("m(\\d+)_(.*)", "Prot-M\\1 (\\2)", final_table_v3$pathway)

# Renaming RNA modules
final_table_v3$pathway <- gsub('Sara_(M\\d+)','RNA-\\1', final_table_v3$pathway)

# Renaming metabolite
final_table_v3$pathway <- gsub('bioenergetic_pathway', 'Bioenergetic Pathway', final_table_v3$pathway)

# Renaming views
final_table_v3 <- final_table_v3  %>%
  mutate(view = gsub("rna_AC", "RNA (AC)", view) %>%
           gsub('rna_PCG', 'RNA (PCG)', .) %>%
           gsub('rna_DLPF', 'RNA (DLPFC)', .) %>%
           gsub('acetil', 'H3K9ac', .) %>%
           gsub('tmt', 'Proteomics', .) %>%
           gsub('metabol', 'Metabolomics', .))

# Rearranging views order
custom_order <- c('H3K9ac', 'RNA (AC)', 'RNA (DLPFC)', 'RNA (PCG)', 'Proteomics','Metabolomics')
final_table_v3$view <- factor(final_table_v3$view, levels = custom_order)
final_table_v3 <- final_table_v3[order(final_table_v3$view), ]

final_table_v3 <- as.data.frame(final_table_v3)
final_table <- final_table_v3
##################################################################################################################################

# Retrieving GO data so it can be used to reduce terms redundancy later on
library(rrvgo)
library(org.Hs.eg.db)
library(msigdbr)
library(tidyverse)
library(purrr)

GO_msigdbr = list()
## GENE ONTOLOGY (MOLECULAR FUNCTION)
GO_msigdbr[["MF"]] <- msigdbr(species = "human", category = "C5", subcategory = "GO:MF")
## GENE ONTOLOGY (BIOLOGICAL PROCESS)
GO_msigdbr[["BP"]] <- msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
## GENE ONTOLOGY (CELLULAR COMPONENT)
GO_msigdbr[["CC"]] <- msigdbr(species = "human", category = "C5", subcategory = "GO:CC")
#########################################################################################################

# Function to add asterisks to significant terms
cell_fun = function(j, i, x, y, width, height, fill, table_1) {
  if(!is.na(table_1[i, j]) & abs(table_1[i, j]) > -log10(0.001)) {
    grid.text("***", x, y, gp = gpar(fontsize = 12), vjust = 0.77)
  } else if(!is.na(table_1[i, j]) & abs(table_1[i, j]) > -log10(0.01)) {
    grid.text("**", x, y, gp = gpar(fontsize = 12), vjust = 0.77)
  } else if(!is.na(table_1[i, j]) & abs(table_1[i, j]) > -log10(0.05)) {
    grid.text("*", x, y, gp = gpar(fontsize = 12), vjust = 0.77)
  }
}

# Main function - produces a ComplexHeatmap of the best terms for each view
arrange_terms <- function(chosen_factor){
  
  ######################################################################################################################## 
  expected_columns <- c('H3K9ac', 'RNA (AC)', 'RNA (DLPFC)', 'RNA (PCG)', 'Proteomics')
  brewer_colors <- brewer.pal(5, 'RdBu')
  coll <- colorRamp2(c(-4, 0, 4), c(brewer_colors[5], "white", brewer_colors[1]))
  heatmap_list <- list()
  NES_list <- list()
  pval_list <- list()
  ############################################## Others ##################################################################
  
  # Keeps track of how many terms are being retrieved
  val_new_lists = 0
  
  # Filtering by factor, padj and by term names (here we retrieve the terms named as 'Others')
  temp_final_table_new_lists <- final_table
  temp_final_table_new_lists <- filter(temp_final_table_new_lists, padj <0.05 & factor == as.character(chosen_factor))
  temp_final_table_new_lists <- filter(temp_final_table_new_lists, pathway %in% c("AD GWAS", "PIG","HAM down","HAM up","DAM down",
                                                                                  "DAM up"))
  if(nrow(temp_final_table_new_lists) > 0){
    
    # Saving the terms to be used on the heatmap
    save_names_new_lists = unique(temp_final_table_new_lists$pathway)
    
    # Mapping the saved terms into a matrix with padj values (transformed to log10) as values
    temp_table_new_lists_v3 <- final_table[final_table$pathway %in% save_names_new_lists,]
    temp_table_new_lists_v3 <- filter(temp_table_new_lists_v3, padj <0.05 & factor == as.character(chosen_factor))
    temp_table_new_lists_v3$padj <- -log10(temp_table_new_lists_v3$padj)*sign(temp_table_new_lists_v3$NES)
    temp_table_new_lists_v3 <- dcast(temp_table_new_lists_v3, view ~ pathway, value.var = 'padj')
    temp_table_new_lists_v3 <- replace(temp_table_new_lists_v3, is.na(temp_table_new_lists_v3), 0)
    rownames(temp_table_new_lists_v3) <- temp_table_new_lists_v3$view
    temp_table_new_lists_v3 <- temp_table_new_lists_v3[, c(-1), drop = FALSE]
    temp_table_new_lists_v3 <- t(temp_table_new_lists_v3)
    
    # Mapping the saved terms into a matrix with NES scores as values
    temp_table_new_lists_v2 <- final_table[final_table$pathway %in% rownames(temp_table_new_lists_v3),]
    temp_table_new_lists_v2 <- filter(temp_table_new_lists_v2, factor == as.character(chosen_factor))
    temp_table_new_lists_v2 <- dcast(temp_table_new_lists_v2, view ~ pathway, value.var = 'NES')
    temp_table_new_lists_v2 <- replace(temp_table_new_lists_v2, is.na(temp_table_new_lists_v2), 0)
    rownames(temp_table_new_lists_v2) <- temp_table_new_lists_v2$view
    temp_table_new_lists_v2 <- temp_table_new_lists_v2[, c(-1), drop = FALSE]
    temp_table_new_lists_v2 <- t(temp_table_new_lists_v2)
    
    
    # Adding missing columns to the matrices (if any)
    expected_cols_lists <- setdiff(expected_columns, colnames(temp_table_new_lists_v2))
    for(cols in expected_cols_lists){
      temp_table_new_lists_v2 <- cbind(temp_table_new_lists_v2, cols = 0)
      colnames(temp_table_new_lists_v2)[ncol(temp_table_new_lists_v2)] <- cols
    }
    expected_cols_lists_v3 <- setdiff(expected_columns, colnames(temp_table_new_lists_v3))
    for(cols in expected_cols_lists_v3){
      temp_table_new_lists_v3 <- cbind(temp_table_new_lists_v3, cols = 0)
      colnames(temp_table_new_lists_v3)[ncol(temp_table_new_lists_v3)] <- cols
    }
    
    # Ordering columns
    temp_table_new_lists_v2 <- temp_table_new_lists_v2[,expected_columns, drop = FALSE]
    temp_table_new_lists_v3 <- temp_table_new_lists_v3[,expected_columns, drop = FALSE]
    
    # Ordering rows
    table_1_new_lists <- as.matrix(temp_table_new_lists_v2[rownames(temp_table_new_lists_v3), , drop=FALSE])
    
    val_new_lists = nrow(table_1_new_lists)
    
    NES_list[['lists']] <- table_1_new_lists
    pval_list[['lists']] <- temp_table_new_lists_v3
    
    # HEATMAP
    row_ha_lists = rowAnnotation(
      Category = as.factor(rep('Others', nrow(table_1_new_lists))),
      col = list(Category = c(`Others` = "#E18727FF")),  
      annotation_height = unit(nrow(table_1_new_lists), "mm"), 
      width = unit(1, "mm"), 
      annotation_label = "  ")
    
    heat_lists=Heatmap(table_1_new_lists,
                       name = 'NES',
                       column_title = paste('Best Terms For Each View'),
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 8)),
                       col =coll,
                       rect_gp = gpar(col = 'black', lwd = 1, lty = 1),
                       border_gp = gpar(col = 'black', lwd      =1),
                       show_row_dend = FALSE,
                       show_column_dend = FALSE,
                       row_names_side = 'left',
                       row_names_gp = gpar(fontsize = 10),
                       show_column_names = TRUE,
                       column_names_rot = 45,
                       column_names_gp = gpar(fontsize = 9),
                       column_order = colnames(table_1_new_lists)[order(as.integer(sub("Factor", "", colnames(table_1_new_lists))),  decreasing = TRUE)],
                       show_row_names = TRUE,
                       cluster_rows = TRUE,
                       right_annotation = row_ha_lists,
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         cell_fun(j, i, x, y, width, height, fill, abs(temp_table_new_lists_v3))
                       })
    heatmap_list <- list()
    heatmap_list[['lists']] <- heat_lists
    
  }
  
  ############################################ CELL TYPES #################################################################
  
  # Keeps track of how many terms are being retrieved
  val_new_cels = 0
  
  # Filtering by factor, padj and by term names (here we retrieve the terms related to Cell types)
  temp_final_table_new_cels <- final_table
  temp_final_table_new_cels <- filter(temp_final_table_new_cels, padj <0.05 & factor == as.character(chosen_factor))
  green_cell_types <- c(
    "Mic.1", "Mic.2", "Mic.3", "Mic.4", "Mic.5", "Mic.6", "Mic.7", "Mic.8", "Mic.9", "Mic.10",
    "Mic.11", "Mic.12", "Mic.13", "Mic.14", "Mic.15", "Mic.16", "Macrophages", "Monocytes", 
    "Ast.1", "Ast.2", "Ast.3", "Ast.4", "Ast.5", "Ast.6", "Ast.7", "Ast.8", "Ast.9", "Ast.10",
    "Oli.1", "Oli.2", "Oli.3", "Oli.4", "Oli.5", "Oli.6", "Oli.7", "Oli.8", "Oli.9", "Oli.10", 
    "Oli.11", "Oli.12", "Oli.13", "OPC.1", "OPC.2", "OPC.3", "NFOL/MOL", "COP", 
    "End.1", "End.2", "End.3", "End.4", "End.5", "Arteriole", "Venule", 
    "Peri.1", "Peri.2", "SMC.1", "SMC.2", "Fib.1", "Fib.2", "Fib.3", 
    "Exc.1", "Exc.2", "Exc.3", "Exc.4", "Exc.5", "Exc.6", "Exc.7", "Exc.8", "Exc.9", "Exc.10", 
    "Exc.11", "Exc.12", "Exc.13", "Exc.14", "Exc.15", "Exc.16", 
    "Inh.1", "Inh.2", "Inh.3", "Inh.4", "Inh.5", "Inh.6", "Inh.7", "Inh.8", "Inh.9", "Inh.10", 
    "Inh.11", "Inh.12", "Inh.13", "Inh.14", "Inh.15", "Inh.16"
  )
  temp_final_table_new_cels <- filter(temp_final_table_new_cels, pathway %in% green_cell_types)
  
  
  save_names_cell_types <- c()
  
  # Saving the most significant terms for each view to be used on the heatmap (limiting to 5 terms)
  while(length(save_names_cell_types)<5 && nrow(temp_final_table_new_cels) != 0){
    
    for(views in unique(as.character(temp_final_table_new_cels$view))){
      n_views = length(unique(as.character(temp_final_table_new_cels$view)))
      
      ordered_terms <- filter(temp_final_table_new_cels, view == views)
      if(nrow(ordered_terms) == 0){next}
      ordered_terms <- ordered_terms[order(ordered_terms$padj),]
      ordered_terms <- ordered_terms$pathway[1]
      
      save_names_cell_types <- c(save_names_cell_types, ordered_terms)
      temp_final_table_new_cels <- temp_final_table_new_cels[temp_final_table_new_cels$pathway != ordered_terms,]
      
      if(length(save_names_cell_types)>=5 ||
         nrow(temp_final_table_new_cels) == 0||
         length(unique(as.character(temp_final_table_new_cels$view))) < n_views){
        break
      }
    }
  }
  
  if(length(save_names_cell_types) > 0){
    
    # Mapping the saved terms into a matrix with padj values (transformed to log10) as values
    temp_table_new_cels_v3 <- final_table[final_table$pathway %in% save_names_cell_types,]
    temp_table_new_cels_v3 <- filter(temp_table_new_cels_v3, padj <0.05 & factor == as.character(chosen_factor))
    temp_table_new_cels_v3$padj <- -log10(temp_table_new_cels_v3$padj)*sign(temp_table_new_cels_v3$NES)
    temp_table_new_cels_v3 <- dcast(temp_table_new_cels_v3, view ~ pathway, value.var = 'padj')
    temp_table_new_cels_v3 <- replace(temp_table_new_cels_v3, is.na(temp_table_new_cels_v3), 0)
    rownames(temp_table_new_cels_v3) <- temp_table_new_cels_v3$view
    temp_table_new_cels_v3 <- temp_table_new_cels_v3[, c(-1), drop = FALSE]
    temp_table_new_cels_v3 <- t(temp_table_new_cels_v3)
    
    
    # Mapping the saved terms into a matrix with NES scores as values
    temp_table_new_cels_v2 <- final_table[final_table$pathway %in% rownames(temp_table_new_cels_v3),]
    temp_table_new_cels_v2 <- filter(temp_table_new_cels_v2, factor == as.character(chosen_factor))
    temp_table_new_cels_v2 <- dcast(temp_table_new_cels_v2, view ~ pathway, value.var = 'NES')
    temp_table_new_cels_v2 <- replace(temp_table_new_cels_v2, is.na(temp_table_new_cels_v2), 0)
    rownames(temp_table_new_cels_v2) <- temp_table_new_cels_v2$view
    temp_table_new_cels_v2 <- temp_table_new_cels_v2[, c(-1), drop = FALSE]
    temp_table_new_cels_v2 <- t(temp_table_new_cels_v2)
    
    # Adding missing columns to the matrices (if any)
    expected_cols_lists <- setdiff(expected_columns, colnames(temp_table_new_cels_v2))
    for(cols in expected_cols_lists){
      temp_table_new_cels_v2 <- cbind(temp_table_new_cels_v2, cols = 0)
      colnames(temp_table_new_cels_v2)[ncol(temp_table_new_cels_v2)] <- cols
    }
    expected_cols_lists_v3 <- setdiff(expected_columns, colnames(temp_table_new_cels_v3))
    for(cols in expected_cols_lists_v3){
      temp_table_new_cels_v3 <- cbind(temp_table_new_cels_v3, cols = 0)
      colnames(temp_table_new_cels_v3)[ncol(temp_table_new_cels_v3)] <- cols
    }
    
    # Ordering columns
    temp_table_new_cels_v2 <- temp_table_new_cels_v2[,expected_columns, drop = FALSE]
    temp_table_new_cels_v3 <- temp_table_new_cels_v3[,expected_columns, drop = FALSE]
    
    # Ordering rows
    table_1_new_cels <- as.matrix(temp_table_new_cels_v2[rownames(temp_table_new_cels_v3), , drop=FALSE])
    
    val_new_cels = nrow(table_1_new_cels)
    
    NES_list[['cells']] <- table_1_new_cels
    pval_list[['cells']] <- temp_table_new_cels_v3
    
    # HEATMAP
    row_ha_cels = rowAnnotation(
      Category = as.factor(rep('Cell Types', nrow(table_1_new_cels))),
      col = list(Category = c(`Cell Types` = "#7876B1FF")),  
      annotation_height = unit(nrow(table_1_new_cels), "mm"),  
      annotation_label = "  ")
    
    heat_cels=Heatmap(table_1_new_cels,
                      name = 'NES',
                      column_title = paste('Best Terms For Each View'),
                      heatmap_legend_param = list(title_gp = gpar(fontsize = 8)),
                      col =coll,
                      rect_gp = gpar(col = 'black', lwd = 1, lty = 1),
                      border_gp = gpar(col = 'black', lwd      =1),
                      show_row_dend = FALSE,
                      show_column_dend = FALSE,
                      row_names_side = 'left',
                      row_names_gp = gpar(fontsize = 10),
                      show_column_names = TRUE,
                      column_names_rot = 45,
                      column_names_gp = gpar(fontsize = 9),
                      column_order = colnames(table_1_new_cels)[order(as.integer(sub("Factor", "", colnames(table_1_new_cels))),  decreasing = TRUE)],
                      show_row_names = TRUE,
                      cluster_rows = TRUE,
                      right_annotation = row_ha_cels,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        cell_fun(j, i, x, y, width, height, fill, abs(temp_table_new_cels_v3))
                      })
    
    heatmap_list[['cells']] <- heat_cels
  }
  
  
  ################################### GO and REACTOME #######################################################################
  
  # Filtering by factor and padj, and by term names — retrieving GO, Reactome, and metabolites terms (excluded later)
  temp_final_table <- final_table
  temp_final_table <- filter(temp_final_table, padj <0.05 & factor == as.character(chosen_factor))
  temp_final_table <- temp_final_table[grep('GOBP|GOMF|GOCC|REACTOME|Urea cycle; 
                                            Arginine and Proline Metabolism|Phospholipid Metabolism|Histidine Metabolism|
                                            Phosphatidylcholine (PC)|Sphingomyelins|Endocannabinoid|Leucine, Isoleucine and Valine Metabolism|
                                            Lysophospholipid|Glutathione Metabolism|Glycine, Serine and Threonine Metabolism|Glycolysis, 
                                            Gluconeogenesis, and Pyruvate Metabolism|Hexosylceramides (HCER)|
                                            Phosphatidylethanolamine (PE)|Purine Metabolism, Adenine containing|
                                            Fatty Acid, Monohydroxy|Lysine Metabolism|Bioenergetic Pathway|
                                            Glutamate Metabolism|Nicotinate and Nicotinamide Metabolism|
                                            Long Chain Polyunsaturated Fatty Acid (n3 and n6)|
                                            Methionine, Cysteine, SAM and Taurine Metabolism',temp_final_table$pathway),]
  
  # Retrieving Reactome terms (KEGG not available)
  kegg_reactome_df <- temp_final_table[grepl('KEGG|REACTOME',temp_final_table$pathway),]
  
  # Retrieving GO terms
  go_df <- temp_final_table[grepl('GOBP|GOMF|GOCC',temp_final_table$pathway),]
  
  # Reducing term redundancy
  gsea_results_list = list()
  simMatrix = list()
  reducedTerms = list()
  top_reducedTerms = list()
  scatterplot_terms = list()
  for(GO_class in c("MF","BP","CC")){
    gsea_results_list[[GO_class]] = go_df %>% 
      inner_join(unique(GO_msigdbr[[GO_class]][,c("gs_exact_source","gs_name")]), by = c("pathway"="gs_name"))
    
    simMatrix[[GO_class]] <- calculateSimMatrix(gsea_results_list[[GO_class]]$gs_exact_source,
                                                orgdb="org.Hs.eg.db",
                                                ont=GO_class,
                                                method="Rel")
    scores <- setNames(-log10(gsea_results_list[[GO_class]]$padj), gsea_results_list[[GO_class]]$gs_exact_source)
    
    if(chosen_factor == '26'){
      reducedTerms[[GO_class]] <- reduceSimMatrix(simMatrix[[GO_class]],
                                                  scores,
                                                  threshold=0.7,
                                                  orgdb="org.Hs.eg.db")
    }else{
      reducedTerms[[GO_class]] <- reduceSimMatrix(simMatrix[[GO_class]],
                                                  scores,
                                                  threshold=0.9,
                                                  orgdb="org.Hs.eg.db")
    }
    
    
    top_reducedTerms[[GO_class]]  = reducedTerms[[GO_class]] %>% 
      group_by(cluster) %>%
      summarise(n = n(), 
                top_parentTerm = parentTerm[which.max(score)], 
                top_parentGO = parent[which.max(score)],
                top_score = max(score)) 
    
    scatterplot_terms[[GO_class]]  = scatterPlot(simMatrix[[GO_class]], reducedTerms[[GO_class]])
  }
  top_reducedTerms_df = purrr::reduce(top_reducedTerms,bind_rows)
  
  GO_msigdbr_map = bind_rows(
    cbind(GO = "MF", unique(GO_msigdbr[["MF"]][,c("gs_exact_source","gs_name")])),
    cbind(GO = "BP", unique(GO_msigdbr[["BP"]][,c("gs_exact_source","gs_name")])),
    cbind(GO = "CC", unique(GO_msigdbr[["CC"]][,c("gs_exact_source","gs_name")])))
  
  gsea_results_filt = go_df %>% left_join(GO_msigdbr_map, by = c("pathway"="gs_name")) %>%
    inner_join(top_reducedTerms_df, by = c("gs_exact_source"="top_parentGO")) %>%
    dplyr::select(-c(top_score,n,cluster,gs_exact_source,GO,top_parentTerm))
  
  temp_final_table <- rbind(gsea_results_filt, kegg_reactome_df)
  
  save_names <- c()
  
  # Saving the most significant terms for each view to be used on the heatmap (limiting to 28 terms)
  while(length(save_names)<28 && nrow(temp_final_table) != 0){
    
    for(views in unique(as.character(temp_final_table$view))){
      n_views = length(unique(as.character(temp_final_table$view)))
      
      ordered_terms <- filter(temp_final_table, view == views)
      if(nrow(ordered_terms) == 0){next}
      ordered_terms <- ordered_terms[order(ordered_terms$padj),]
      ordered_terms <- ordered_terms$pathway[1]
      
      save_names <- c(save_names, ordered_terms)
      temp_final_table <- temp_final_table[temp_final_table$pathway != ordered_terms,]
      
      if(length(save_names)>=28 ||
         nrow(temp_final_table) == 0||
         length(unique(as.character(temp_final_table$view))) < n_views){
        break
      }
    }
  }
  
  # Mapping the saved terms into a matrix with padj values (transformed to log10) as values
  temp_table_v3 <- final_table[final_table$pathway %in% save_names,]
  temp_table_v3 <- filter(temp_table_v3, padj <0.05 & factor == as.character(chosen_factor))
  temp_table_v3$padj <- -log10(temp_table_v3$padj)*sign(temp_table_v3$NES)
  temp_table_v3 <- dcast(temp_table_v3, view ~ pathway, value.var = 'padj')
  temp_table_v3 <- replace(temp_table_v3, is.na(temp_table_v3), 0)
  rownames(temp_table_v3) <- temp_table_v3$view
  temp_table_v3 <- temp_table_v3[, c(-1), drop = FALSE]
  temp_table_v3 <- t(temp_table_v3)
  
  # Mapping the saved terms into a matrix with NES scores as values
  temp_table_v2 <- final_table[final_table$pathway %in% rownames(temp_table_v3),]
  temp_table_v2 <- filter(temp_table_v2, factor == as.character(chosen_factor))
  temp_table_v2 <- dcast(temp_table_v2, view ~ pathway, value.var = 'NES')
  temp_table_v2 <- replace(temp_table_v2, is.na(temp_table_v2), 0)
  rownames(temp_table_v2) <- temp_table_v2$view
  temp_table_v2 <- temp_table_v2[, c(-1), drop = FALSE]
  temp_table_v2 <- t(temp_table_v2)
  
  
  # Adding missing columns to the matrices (if any)
  expected_cols_lists <- setdiff(expected_columns, colnames(temp_table_v2))
  for(cols in expected_cols_lists){
    temp_table_v2 <- cbind(temp_table_v2, cols = 0)
    colnames(temp_table_v2)[ncol(temp_table_v2)] <- cols
  }
  expected_cols_lists_v3 <- setdiff(expected_columns, colnames(temp_table_v3))
  for(cols in expected_cols_lists_v3){
    temp_table_v3 <- cbind(temp_table_v3, cols = 0)
    colnames(temp_table_v3)[ncol(temp_table_v3)] <- cols
  }
  
  # Ordering columns
  temp_table_v2 <- temp_table_v2[,expected_columns, drop = FALSE]
  temp_table_v3 <- temp_table_v3[,expected_columns, drop = FALSE]
  
  # Ordering rows
  table_1 <- as.matrix(temp_table_v2[rownames(temp_table_v3), , drop=FALSE])
  
  assign(paste0('F',chosen_factor,'names'),table_1,envir = .GlobalEnv)
  
  library(stringr)
  
  # Limiting name lenght
  rownames(table_1) <- str_trunc(rownames(table_1), 37, side = 'left')
  
  val_db = nrow(table_1)
  
  # Clearing column names
  cols <- gsub("_", " ", tolower(rownames(table_1)))
  cols <- sub("^gobp ", "", cols)
  cols <- sub("^gocc ", "", cols)
  cols <- sub("^gomf ", "", cols)
  cols <- sub("^kegg ", "", cols)
  cols <- sub("^reactome ", "", cols)
  rownames(table_1) <- cols
  
  NES_list[['DB']] <- table_1
  pval_list[['DB']] <- temp_table_v3
  
  
  # HEATMAP
  row_ha_DB = rowAnnotation(
    Category = as.factor(rep('GO\nReactome', nrow(table_1))),
    col = list(Category = c(`GO\nReactome` = "#20854EFF")),  
    annotation_height = unit(nrow(table_1), "mm"),  
    annotation_label = "  ")
  
  heat_DB=Heatmap(table_1,
                  name = 'NES',
                  column_title = paste('Best Terms For Each View'),
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 8)),
                  col =coll,
                  rect_gp = gpar(col = 'black', lwd = 1, lty = 1),
                  border_gp = gpar(col = 'black', lwd      =1),
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  row_names_side = 'left',
                  row_names_gp = gpar(fontsize = 10),
                  show_column_names = TRUE,
                  column_names_rot = 45,
                  column_names_gp = gpar(fontsize = 9),
                  column_order = colnames(table_1)[order(as.integer(sub("Factor", "", colnames(table_1))),  decreasing = TRUE)],
                  show_row_names = TRUE,
                  cluster_rows = TRUE,
                  right_annotation = row_ha_DB,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    cell_fun(j, i, x, y, width, height, fill, abs(temp_table_v3))
                  })
  
  heatmap_list[['DB']] <- heat_DB
  
  
  
  ################################## MODULES ###########################################################################
  
  # Filtering by factor, padj, and term names — here we retrieve terms related to protein and RNA modules
  temp_final_table_modules <- final_table
  temp_final_table_modules <- filter(temp_final_table_modules, padj <0.05 & factor == as.character(chosen_factor))
  temp_final_table_modules <- temp_final_table_modules[grep('Prot-|RNA-',temp_final_table_modules$pathway),]
  
  save_names_modules <- c()
  
  # Limiting the number of modules terms based on the maximum allowed in the heatmap (40)
  lim_terms_modules = 40 - (val_new_lists + val_new_cels + val_db)
  
  # Saving the most significant terms for each view to be used on the heatmap
  while(length(save_names_modules)<lim_terms_modules && nrow(temp_final_table_modules) != 0){
    
    for(views in unique(as.character(temp_final_table_modules$view))){
      n_views = length(unique(as.character(temp_final_table_modules$view)))
      
      ordered_terms <- filter(temp_final_table_modules, view == views)
      ordered_terms <- ordered_terms[order(ordered_terms$padj),]
      ordered_terms <- ordered_terms$pathway[1]
      
      save_names_modules <- c(save_names_modules, ordered_terms)
      temp_final_table_modules <- temp_final_table_modules[temp_final_table_modules$pathway != ordered_terms,]
      
      if(length(save_names_modules)>=lim_terms_modules ||
         nrow(temp_final_table_modules) == 0 || 
         length(unique(as.character(temp_final_table_modules$view))) < n_views){
        break
      }
    }
  }
  
  # Mapping the saved terms into a matrix with padj values (transformed to log10) as values
  temp_table_modules_v3 <- final_table[final_table$pathway %in% save_names_modules,]
  temp_table_modules_v3 <- filter(temp_table_modules_v3, padj <0.05 & factor == as.character(chosen_factor))
  temp_table_modules_v3$padj <- -log10(temp_table_modules_v3$padj)*sign(temp_table_modules_v3$NES)
  temp_table_modules_v3 <- dcast(temp_table_modules_v3, view ~ pathway, value.var = 'padj')
  temp_table_modules_v3 <- replace(temp_table_modules_v3, is.na(temp_table_modules_v3), 0)
  rownames(temp_table_modules_v3) <- temp_table_modules_v3$view
  temp_table_modules_v3 <- temp_table_modules_v3[, c(-1), drop = FALSE]
  temp_table_modules_v3 <- t(temp_table_modules_v3)
  
  # Mapping the saved terms into a matrix with NES scores as values
  temp_table_modules_v2 <- final_table[final_table$pathway %in% rownames(temp_table_modules_v3),]
  temp_table_modules_v2 <- filter(temp_table_modules_v2, factor == as.character(chosen_factor))
  temp_table_modules_v2 <- dcast(temp_table_modules_v2, view ~ pathway, value.var = 'NES')
  temp_table_modules_v2 <- replace(temp_table_modules_v2, is.na(temp_table_modules_v2), 0)
  rownames(temp_table_modules_v2) <- temp_table_modules_v2$view
  temp_table_modules_v2 <- temp_table_modules_v2[, c(-1), drop = FALSE]
  temp_table_modules_v2 <- t(temp_table_modules_v2)
  
  # Adding missing columns to the matrices (if any)
  expected_cols_lists <- setdiff(expected_columns, colnames(temp_table_modules_v2))
  for(cols in expected_cols_lists){
    temp_table_modules_v2 <- cbind(temp_table_modules_v2, cols = 0)
    colnames(temp_table_modules_v2)[ncol(temp_table_modules_v2)] <- cols
  }
  expected_cols_lists_v3 <- setdiff(expected_columns, colnames(temp_table_modules_v3))
  for(cols in expected_cols_lists_v3){
    temp_table_modules_v3 <- cbind(temp_table_modules_v3, cols = 0)
    colnames(temp_table_modules_v3)[ncol(temp_table_modules_v3)] <- cols
  }
  
  # Ordering columns
  temp_table_modules_v2 <- temp_table_modules_v2[,expected_columns, drop = FALSE]
  temp_table_modules_v3 <- temp_table_modules_v3[,expected_columns, drop = FALSE]
  
  # Ordering rows
  table_1_modules <- as.matrix(temp_table_modules_v2[rownames(temp_table_modules_v3), , drop=FALSE])
  
  # Clearing row names
  rownames(table_1_modules) <- gsub('_',' ',rownames(table_1_modules))
  
  NES_list[['modules']] <- table_1_modules
  pval_list[['modules']] <- temp_table_modules_v3
  
  row_ha_modules = rowAnnotation(
    Category = as.factor(rep('Modules', nrow(table_1_modules))),
    col = list(Category = c(`Modules` = "#6F99ADFF")),  
    annotation_height = unit(nrow(table_1_modules), "mm"),  
    annotation_label = "  ")
  
  # HEATMAP
  heat_modules=Heatmap(table_1_modules,
                       name = 'NES',
                       column_title = paste('Best Terms For Each View'),
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 8)),
                       col =coll,
                       rect_gp = gpar(col = 'black', lwd = 1, lty = 1),
                       border_gp = gpar(col = 'black', lwd      =1),
                       show_row_dend = FALSE,
                       show_column_dend = FALSE,
                       row_names_side = 'left',
                       row_names_gp = gpar(fontsize = 10),
                       show_column_names = TRUE,
                       column_names_rot = 45,
                       column_names_gp = gpar(fontsize = 9),
                       column_order = colnames(table_1_modules)[order(as.integer(sub("Factor", "", colnames(table_1_modules))),  decreasing = TRUE)],
                       show_row_names = TRUE,
                       cluster_rows = TRUE,
                       right_annotation = row_ha_modules,
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         cell_fun(j, i, x, y, width, height, fill, abs(temp_table_modules_v3))
                       })
  
  heatmap_list[['modules']] <- heat_modules
  
  
  ######################################## VARIANCE EXPLAINED #####################################################################
  
  # MOFA function that retrieves the % variance explained for all factors
  heat_1 <- plot_variance_explained(mofa, max_r2=15)$data
  
  # Filtering by factor and excluding the views 'metabolites' and 'cell types'
  heat_1 <- filter(heat_1, factor == paste0('Factor',chosen_factor), !(view %in% c('sn_RNA_norm','metabol')))
  
  # Reshaping the matrix
  heat_1 <- dcast(heat_1, factor ~ view, value.var = "value")
  rownames(heat_1) <- heat_1$factor
  heat_1 <- heat_1[,-1]
  heat_1 <- heat_1[,c(5,1,2,3,4)]
  colnames(heat_1) <- c('H3K9ac', 'RNA (AC)', 'RNA (DLPFC)', 'RNA (PCG)', 'Proteomics')
  heat_1 <- as.matrix(heat_1)
  
  col_fun = colorRamp2(
    seq(min(heat_1), max(heat_1), length.out = 5),  
    pal_material("purple")(5)
  )
  
  annotation = rowAnnotation(
    explanation = anno_text(
      rep("% variance explained"),
      gp = gpar(fontsize = 8, fontface = "bold")
    )
  )
  
  # HEATMAP
  heat_var_xp <- Heatmap(heat_1,
                         column_title = paste('Best Terms For Each View'),
                         show_heatmap_legend = FALSE,
                         col=col_fun,
                         show_row_names = FALSE,
                         show_column_names = FALSE,
                         show_row_dend = FALSE,
                         show_column_dend = FALSE,
                         column_order = colnames(heat_1),
                         right_annotation = annotation,
                         cell_fun = function(j, i, x, y, width, height, fill) {  
                           grid.text(sprintf("%.2f", heat_1[i, j]),
                                     x, y, gp = gpar(col = "black", fontsize = 8))}
  )
  heatmap_list[['var_xp']] <- heat_var_xp
  ###########################################################################################################
  
  all_matrix_NES<- do.call(rbind,NES_list)
  all_matrix_pval <- do.call(rbind, pval_list)
  
  heatmap_order <- c('var_xp','DB','lists','modules','cells')
  heatmap_list <- heatmap_list[intersect(heatmap_order, names(heatmap_list))]
  combined_heatmap <- Reduce(`%v%`, heatmap_list)
  draw(combined_heatmap, gap = unit(1, 'mm'))
  
  
}

#ex:
#arrange_terms(4)
