#######################################
# Factors enrichment characterization #
#######################################
library(dplyr)
library(reshape2)
library(ComplexHeatmap)

# Function to identify terms through keywords
find_matches <- function(pathway, keywords) {
  matches <- keywords[sapply(keywords, function(k) grepl(k, pathway, ignore.case = TRUE))]
  paste(matches, collapse = ", ")  # Return matched keywords as a comma-separated string
}

# Retrieving GSEA data
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Best_mofas_analysis/standart_50f/final_table_v3.RDATA")
final_table_v3 <- filter(final_table_v3, padj < 0.05)

##################################################################################################################################################

# Immune keywords
Immune <- c('IMMUNE','INFLAMMATORY','LEUKOCYTE','IMMUNITY','LYMPHOCYTE','HUMORAL','IMMUNOGLOBULIN',
            'DEFENSE', 'COMPLEMENT', 'VIRUS', 'BACTERIUM', 'FUNGUS','_B_CELL','_T_CELL','CYTOKINE',
            'TUMOR_NECROSIS','INTERFERON','INTERLEUKIN','MHC','T_HELPER','GRANULOCYTE'
            )

# Applying function 'find_matches'
final_table_v3$matched_keyword <- sapply(final_table_v3$pathway, find_matches, keywords = Immune)

# Filtering to include only terms where a keyword was found
final_table_v3_matched <- final_table_v3[final_table_v3$matched_keyword != "", ]

# Keep only terms where a match was found
Immune = final_table_v3_matched[, c("pathway", "matched_keyword")]

# Excluding terms that overlap with other categories
Immune = filter(Immune, !(pathway %in% c("GOBP_LEUKOCYTE_ADHESION_TO_VASCULAR_ENDOTHELIAL_CELL", 
                                         "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_ADHESION_TO_VASCULAR_ENDOTHELIAL_CELL")))

##################################################################################################################################################

# Mitochondria keywords
Mitochondria <- c('MITOCHONDRIAL','MITOCHONDRIA','MITOCHONDRION','OXIDATIVE_PHOSPHORYLATION','ELECTRON_TRANSPORT',
                  'RESPIRATION','REACTIVE_OXYGEN_SPECIES', 'OXIDOREDUTASE_COMPLEX','RESPIROSSOME')

# Applying function 'find_matches'
final_table_v3$matched_keyword <- sapply(final_table_v3$pathway, find_matches, keywords = Mitochondria)

# Filtering to include only terms where a keyword was found
final_table_v3_matched <- final_table_v3[final_table_v3$matched_keyword != "", ]

# Keep only terms where a match was found
Mitochondria = final_table_v3_matched[, c("pathway", "matched_keyword")]

# Excluding terms that overlap with other categories
Mitochondria = filter(Mitochondria, !(pathway %in% c("GOCC_MITOCHONDRIAL_LARGE_RIBOSOMAL_SUBUNIT", "GOCC_MITOCHONDRIAL_SMALL_RIBOSOMAL_SUBUNIT")))

##################################################################################################################################################

# Ribosome kewwords
Ribosome <- c('RIBOSOME','RIBOSOMAL')

# Applying function 'find_matches'
final_table_v3$matched_keyword <- sapply(final_table_v3$pathway, find_matches, keywords = Ribosome)

# Filtering to include only terms where a keyword was found
final_table_v3_matched <- final_table_v3[final_table_v3$matched_keyword != "", ]

# Keep only terms where a match was found
Ribosome = final_table_v3_matched[, c("pathway", "matched_keyword")]

##################################################################################################################################################

# Proteostasis kewwords
Proteostasis <-c('_UNFOLDED_PROTEIN','HEAT_SHOCK','PROTEOLYSIS','CHAPERONE','TOPOLOGICALLY_INCORRECT_PROTEIN','HSF1',
                 'FOLDING')

# Applying function 'find_matches'
final_table_v3$matched_keyword <- sapply(final_table_v3$pathway, find_matches, keywords = Proteostasis)

# Filtering to include only terms where a keyword was found
final_table_v3_matched <- final_table_v3[final_table_v3$matched_keyword != "", ]

# Keep only terms where a match was found
Proteostasis = final_table_v3_matched[, c("pathway", "matched_keyword")]

##################################################################################################################################################

# Cytoskeleton kewwords
Cytoskeleton <- c('MICROTUBULE','CYTOSKELETON','CYTOSKELETAL','_ACTIN_',
                  'STRESS_FIBER','SUPRAMOLECULAR_FIBER', 'ACTOMYOSIN','_RHO_')

# Applying function 'find_matches'
final_table_v3$matched_keyword <- sapply(final_table_v3$pathway, find_matches, keywords = Cytoskeleton)

# Filtering to include only terms where a keyword was found
final_table_v3_matched <- final_table_v3[final_table_v3$matched_keyword != "", ]

# Keep only terms where a match was found
Cytoskeleton = final_table_v3_matched[, c("pathway", "matched_keyword")]

##################################################################################################################################################

# Vasculature keywords
Vasculature <- c('ANGIOGENESIS', 'VASCULAR', 'ENDOTHELIUM', 'ENDOTHELIAL', 'BLOOD_VESSEL',  
                 'VASOCONSTRICTION',  'ARTERY', 'VENULE', 'ARTERIOLE', 'CIRCULATION','CYRCULATORY', 'ISCHEMIA','VEGF')  

# Applying function 'find_matches'
final_table_v3$matched_keyword <- sapply(final_table_v3$pathway, find_matches, keywords = Vasculature)

# Filtering to include only terms where a keyword was found
final_table_v3_matched <- final_table_v3[final_table_v3$matched_keyword != "", ]

# Keep only terms where a match was found
Vasculature = final_table_v3_matched[, c("pathway", "matched_keyword")]

##################################################################################################################################################

# RNA processing keywords
rna_processing <- c('RNA_PROCESSING','SPLICING','SPLICEOSOME','SNRNP','INTRON','EXON',
                    'RNA_MODIFICATION','RNA_BINDING','RIBONUCLEASE',
                    'RNA_TRANSPORT')

# Applying function 'find_matches'
final_table_v3$matched_keyword <- sapply(final_table_v3$pathway, find_matches, keywords = rna_processing)

# Filtering to include only terms where a keyword was found
final_table_v3_matched <- final_table_v3[final_table_v3$matched_keyword != "", ]

# Keep only terms where a match was found
rna_processing = final_table_v3_matched[, c("pathway", "matched_keyword")]

##################################################################################################################################################

# Synapse keywords
Synapse <- c('NEURON_PROJECTION','AXON','SYNAPSE', 'SYNAPTIC',
             'DENDRITE','ACTION_POTENTIAL',
             'NEUROTRANSMITTER')

# Applying function 'find_matches'
final_table_v3$matched_keyword <- sapply(final_table_v3$pathway, find_matches, keywords = Synapse)

# Filtering to include only terms where a keyword was found
final_table_v3_matched <- final_table_v3[final_table_v3$matched_keyword != "", ]

# Keep only terms where a match was found
Synapse = final_table_v3_matched[, c("pathway", "matched_keyword")]

##################################################################################################################################################

# Heatmaps
heat_map_list <- list()
for(categories in c('Immune',
                    'Mitochondria',
                    'Ribosome',
                    'Synapse',
                    'Proteostasis',
                    'Cytoskeleton',
                    'Vasculature',
                    'rna_processing')){
  
  #' Some factors might have the same term in more then one view, 
  #' here we keep only the term for the view with the lowest p-adj for each factor (excluding duplicates)
  filter_pathways = filter(final_table_v3,pathway %in% get(categories)$pathway)
  filter_padj <- filter_pathways %>%
    group_by(pathway,factor) %>%
    slice_min(order_by = padj, n=1) %>%
    ungroup()
  
  # Select the top 30 most significant terms per factor
  top_terms <- filter_padj %>%
    group_by(factor) %>%
    slice_min(order_by = padj, n = 30) %>%
    ungroup()
  
  # Log transformation
  top_terms$padj <- -log10(top_terms$padj)
  
  # Reshape data, handle missing values, set row names, and adjust columns
  to_matrix <- dcast(top_terms, factor ~ pathway, value.var = 'padj')
  to_matrix <- replace(to_matrix, is.na(to_matrix), 0)
  rownames(to_matrix) <- paste('Factor', to_matrix$factor)
  to_matrix <- to_matrix[, c(-1), drop = FALSE]
  
  
  expected_rows_lists <- setdiff(c("Factor 1",
                                   "Factor 2",
                                   "Factor 3",
                                   "Factor 4",
                                   "Factor 8",
                                   "Factor 12",
                                   "Factor 14",
                                   "Factor 26",
                                   "Factor 42"), rownames(to_matrix))
  
  # Add missing rows with zeros and convert to matrix
  for(rowss in expected_rows_lists){
    to_matrix <- rbind(to_matrix, rowss = 0)
    rownames(to_matrix)[nrow(to_matrix)] <- rowss
  }
  to_matrix <- as.matrix(to_matrix)
  
  # Order rows numerically in decreasing order based on row names
  row_order <- order(as.numeric(gsub("Factor ", "", rownames(to_matrix))), decreasing = TRUE)
  to_matrix <- to_matrix[row_order, , drop = FALSE]
  
  if(categories == 'rna_processing'){
    categories  = 'RNA processing'
  }
  
  library(circlize)
  library(RColorBrewer)
  
  # Setting colors and intervals
  color_palette <- brewer.pal(9, "YlOrRd")
  col_H4 = colorRamp2(c(0,5,10,15,30,50,70), c('white',color_palette[4:9]))
  
  # Heatmap
  heat_map <- Heatmap(to_matrix,
                      name = '-log(padj)',
                      column_title_gp = gpar(fontsize = 28),
                      column_title = paste(categories),
                      heatmap_legend_param = list(title_gp = gpar(fontsize = 22),
                                                  at = c(0,5,10,15,30,50,70),
                                                  labels = c(0,5,10,15,30,50,70),
                                                  break_dist = 1.1,
                                                  legend_height = unit(6,'cm'),
                                                  grid_width = unit(1, "cm"),
                                                  labels_gp = gpar(fontsize = 18),
                                                  legend_gap = unit(1.5, "cm")
                      ),
                      col =col_H4,
                      rect_gp = gpar(col = 'black', lwd = 1, lty = 0.1),
                      border_gp = gpar(col = 'black', lwd      =1),
                      show_row_dend = FALSE,
                      show_column_dend = FALSE,
                      row_names_side = 'left',
                      row_names_gp = gpar(fontsize = 28),
                      show_column_names = FALSE,
                      column_title_rot = 45,
                      show_row_names = TRUE,
                      cluster_rows = FALSE,
                      row_order = order(as.numeric(rownames(to_matrix))),
                      column_title_side = "bottom"
  )
  
  heat_map_list[[categories]] <- heat_map
  
}

# Ploting final heatmap
combined_heatmap <- Reduce(`+`, heat_map_list)
draw(combined_heatmap, gap = unit(1, 'mm'))