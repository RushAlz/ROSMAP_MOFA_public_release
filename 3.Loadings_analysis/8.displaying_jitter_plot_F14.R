library(dplyr)
library(MOFA2)
library(stringr)

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# MOFA function that retrieves all weights
all_weights=get_weights(mofa,factors = 14,scale = TRUE, as.data.frame = TRUE)

# Filtering acetylation peaks
filterd_for_acetil <- filter(all_weights, view == 'acetil')

# Retrieving metadata for histone acetylation
metadata <- readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/histoneAcetylation_H3K9ac.metadata.rds")
metadata <- metadata[,c('probeid','closest_gene')]

# Preparing 'metadata' to merge with 'filterd_for_acetil'. The goal is to map 'filterd_for_acetil' peaks to the closest genes in 'metadata'
metadata <- subset(metadata, probeid %in% filterd_for_acetil$feature)
colnames(metadata)[colnames(metadata) == "probeid"] <- "feature"

# Merging
merge_metadata <- merge(metadata, filterd_for_acetil, by = 'feature')
merge_metadata <-merge_metadata[,-1]
colnames(merge_metadata)[1] <- 'feature'

# Clearing up gene names
merge_metadata$feature <- gsub('\\..*','',merge_metadata$feature)
merge_metadata$feature <- gsub('-.*','',merge_metadata$feature)

# Some peaks share the same closest gene; so here we are keeping only one entry per gene (the one with higher absolute loading)
acetil_df <- merge_metadata %>%
  group_by(factor, feature) %>%
  filter(
    if (any(value >= 0)) value == max(value) else value == min(value)
  ) %>%
  ungroup()

# Mergin all features weights with the new acetil names
all_weights <- subset(all_weights, all_weights$view != 'acetil')
all_weights <- rbind(all_weights, acetil_df)

# Clean and standardize feature names, extract sign, and convert values to absolute
cleared_df = all_weights %>%
  mutate(
    sign = case_when(
      value >= 0 ~ '+',
      value < 0 ~ '-'),
    value = abs(value),
    feature_id = feature,
    feature = case_when(
      view %in% c('rna_AC','rna_DLPF', 'rna_PCG') ~ str_replace_all(str_replace(feature, '.*_(.*)_(.*)_(.*)', '\\1'), '\\.',''),
      view == 'tmt' ~ str_replace_all(str_replace(feature, '.*-',''), '.*_',''),
      TRUE ~ feature
    )
    
  )

##########################################################
# DISPLAY JITTER PLOT WITH FEATURE WEIGHTS FOR FACTOR 14 #
##########################################################
library(readxl)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(patchwork)

# Function to display a jitter plot with feature weights
jitter_ploting_weights <- function(df,categories,filtered_val){
  
  # Box.padding
  bb=1
  # Point.padding
  pp= 0.2
  # Font size
  siz = 6
  
  # Function to select which features are going to be displayed
  choosing_names <- function(View){
    
    # Filtering 'df' by view (df contains all feature weights for the factor)
    df_view <- filter(df, view == View)
    df_view$feature_adjusted <- as.numeric(as.factor(df_view$feature))
    df_view$labels <- 'nshow'
    
    # Number of positive and negative features per view to be displayed
    num_show = 6
    
    #' Filtering 'df_view' by:
    #' combined_genes: genes that belong to selected pathways;
    #' repeated_genes: genes appearing in multiple views;
    #' g_was: genes from a genome-wide association study;
    #' value > filtered_val: feature weights with an absolute value greater than filtered_val
    temp_view <- filter(df_view, feature %in% unique(c(combined_genes,repeated_genes,g_was)) & value >filtered_val)
    
    #' Further filtering 'temp_view' by:
    #' g_was;
    #' repeated_genes;
    #' store_genes: custom vector of genes related to Alzheimer's Disease;
    #' combined_genes;
    #' real_values > 0: features with positive weights
    pos_gwas <- filter(temp_view, feature %in% g_was & real_values > 0)$feature
    pos_rep_terms <- filter(temp_view, feature %in% repeated_genes & real_values >0)$feature
    pos_weights_list <- filter(temp_view, feature %in% store_genes & real_values > 0)$feature
    pos_weights_geral <- filter(temp_view, feature %in% combined_genes & real_values > 0)$feature
    
    # Selecting genes in a priority order, until it reaches num_show
    pos_genes_to_show <- c()
    for(genes in c(pos_gwas,pos_rep_terms,pos_weights_list,pos_weights_geral)){
      pos_genes_to_show <- unique(c(pos_genes_to_show,genes))
      pos_genes_to_show <- head(pos_genes_to_show, n = num_show)
      if(length(pos_genes_to_show) == num_show){
        break
      }
    }
    
    # Marking the selected features to be displayed
    df_view$labels <- ifelse(df_view$feature %in% pos_genes_to_show, 'yshow',df_view$labels)
    
    # Repeating the process for features with negative weights
    neg_gwas <- filter(temp_view, feature %in% g_was & real_values < 0)$feature
    neg_rep_terms <- filter(temp_view, feature %in% repeated_genes & real_values <0)$feature
    neg_weights_list <- filter(temp_view, feature %in% store_genes & real_values < 0)$feature
    neg_weights_geral <- filter(temp_view, feature %in% combined_genes & real_values < 0)$feature
    
    neg_genes_to_show <- c()
    for(genes in c(neg_gwas,neg_rep_terms,neg_weights_list,neg_weights_geral)){
      neg_genes_to_show <- unique(c(neg_genes_to_show,genes))
      neg_genes_to_show <- head(neg_genes_to_show, n = num_show)
      if(length(neg_genes_to_show) == num_show){
        break
      }
    }
    df_view$labels <- ifelse(df_view$feature %in% neg_genes_to_show, 'yshow',df_view$labels)
    return(df_view)
  }
  
  ############################################################################################################################
  
  df_acetil <- choosing_names('acetil')
  
  # Jitter plot
  acetil_plot <- ggplot(df_acetil, aes(x = feature, y = real_values, color = category)) +
    geom_jitter(aes(alpha = ifelse(category == "genes", 0.4, 1)), 
                width = 0.2, 
                size = 3.5, 
                shape = 20) +
    geom_text_repel(aes(
      label = ifelse(labels == 'yshow', feature, ""),
      fontface = ifelse(feature %in% g_was, 'bold','plain')),
      size = siz,
      box.padding = bb,
      point.padding = pp,
      max.overlaps = Inf,
      show.legend = FALSE)+
    scale_color_manual(values = categories) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 24),
      panel.grid.major = element_blank(),        
      panel.grid.minor = element_blank(),        
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 17),
      axis.text.y = element_text(size = 15)) +
    labs(title = "H3K9ac", x = NULL, y = "Weights") +
    scale_alpha_identity() +
    scale_y_continuous(limits = c(-1.02, 1.02), 
                       expand = expansion(0, 0),
                       breaks = c(-1.0,filtered_val*-1,0.0,filtered_val,1.0))
  
  ############################################################################################################################
  
  df_AC <- choosing_names('rna_AC')
  
  # Jitter plot
  AC_plot <- ggplot(df_AC, aes(x = feature, y = real_values, color = category)) +
    geom_jitter(aes(alpha = ifelse(category == "genes", 0.4, 1)), 
                width = 0.2, 
                size = 3.5, 
                shape = 20) +
    geom_text_repel(aes(
      label = ifelse(labels == 'yshow', feature, ""),
      fontface = ifelse(feature %in% g_was, 'bold','plain')),
      size = siz,
      box.padding = bb,
      point.padding = pp,
      max.overlaps = Inf,
      show.legend = FALSE)+
    scale_color_manual(values = categories) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 24),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),        
      panel.grid.minor = element_blank(),        
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()) +
    labs(title = "RNA (AC)", x = NULL, y = NULL) +
    scale_alpha_identity() +
    scale_y_continuous(limits = c(-1.02, 1.02), 
                       expand = expansion(0, 0))
  
  ############################################################################################################################
  
  df_DLPFC <- choosing_names('rna_DLPF')
  
  # Jitter plot
  DLPFC_plot <-ggplot(df_DLPFC, aes(x = feature_adjusted, y = real_values, color = category)) +
    geom_jitter(aes(alpha = ifelse(category == "genes", 0.4, 1)), 
                width = 0.2, 
                size = 3.5, 
                shape = 20) +
    geom_text_repel(aes(
      label = ifelse(labels == 'yshow', feature, ""),
      fontface = ifelse(feature %in% g_was, 'bold','plain')),
      size = siz,
      box.padding = bb,
      point.padding = pp,
      max.overlaps = Inf,
      show.legend = FALSE)+
    scale_color_manual(values = categories) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 24),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),        
      panel.grid.minor = element_blank(),        
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()) +
    labs(title = "RNA (DLPFC)", x = NULL, y = NULL) +
    scale_alpha_identity() +
    scale_y_continuous(limits = c(-1.02, 1.02), 
                       expand = expansion(0, 0))+
    scale_x_continuous(expand = c(0, 0))
  
  ###########################################################################################################
  
  df_PCG <- choosing_names('rna_PCG')
  
  # Jitter plot
  PCG_plot <- ggplot(df_PCG, aes(x = feature_adjusted, y = real_values, color = category)) +
    geom_jitter(aes(alpha = ifelse(category == "genes", 0.4, 1)), 
                width = 0.2, 
                size = 3.5, 
                shape = 20) +
    geom_text_repel(aes(
      label = ifelse(labels == 'yshow', feature, ""),
      fontface = ifelse(feature %in% g_was, 'bold','plain')),
      size = siz,
      box.padding = bb,
      point.padding = pp,
      max.overlaps = Inf,
      show.legend = FALSE)+
    scale_color_manual(values = categories) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 24),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),        
      panel.grid.minor = element_blank(),        
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()) +
    labs(title = "RNA (PCG)", x = NULL, y = NULL) +
    scale_alpha_identity() +
    scale_y_continuous(limits = c(-1.02, 1.02), 
                       expand = expansion(0, 0))+
    scale_x_continuous(expand = c(0, 0))
  
  ###########################################################################################################################
  
  df_tmt <- choosing_names('tmt')
  
  df_tmt <- df_tmt %>%
    mutate(feature_adjusted = case_when(
      feature == 'MGST' ~ 1,
      TRUE ~ feature_adjusted
    ))
  
  # Jitter plot
  tmt_plot <- ggplot(df_tmt, aes(x = feature_adjusted, y = real_values, color = category)) +
    geom_jitter(aes(alpha = ifelse(category == "genes", 0.4, 1)),
                width = 0.5,
                size = 3.5,
                shape = 20) +
    geom_text_repel(aes(
      label = ifelse(labels == 'yshow', feature, ""),
      fontface = ifelse(feature %in% g_was, 'bold','plain')),
      size = siz,
      box.padding = bb,
      point.padding = pp,
      max.overlaps = Inf,
      show.legend = FALSE,
      parse = FALSE)+
    scale_color_manual(values = categories, drop = FALSE) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 24),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 18)
    ) +
    labs(title = "Proteomics", x = NULL, y = NULL) +
    scale_alpha_identity() +
    scale_y_continuous(limits = c(-1.02, 1.02),
                       expand = expansion(0, 0)) +
    scale_x_continuous(expand = c(0, 0))  
  
  # Adding dotted line on weight = 0.4
  acetil_plot <- acetil_plot + geom_hline(yintercept = filtered_val, linetype = "dotted", color = "black") + guides(color = "none")
  AC_plot <- AC_plot + geom_hline(yintercept = filtered_val, linetype = "dotted", color = "black") + guides(color = "none")
  DLPFC_plot <- DLPFC_plot + geom_hline(yintercept = filtered_val, linetype = "dotted", color = "black") + guides(color = "none")
  PCG_plot <- PCG_plot + geom_hline(yintercept = filtered_val, linetype = "dotted", color = "black") + guides(color = "none")
  tmt_plot <- tmt_plot + geom_hline(yintercept = filtered_val, linetype = "dotted", color = "black") + guides(color = "none")
  
  acetil_plot <- acetil_plot + geom_hline(yintercept = -filtered_val, linetype = "dotted", color = "black") + guides(color = "none")
  AC_plot <- AC_plot + geom_hline(yintercept = -filtered_val, linetype = "dotted", color = "black") + guides(color = "none")
  DLPFC_plot <- DLPFC_plot + geom_hline(yintercept = -filtered_val, linetype = "dotted", color = "black") + guides(color = "none")
  PCG_plot <- PCG_plot + geom_hline(yintercept = -filtered_val, linetype = "dotted", color = "black") + guides(color = "none")
  tmt_plot <- tmt_plot + geom_hline(yintercept = -filtered_val, linetype = "dotted", color = "black") + guides(color = "none")
  
  # Adding legend
  acetil_plot <- acetil_plot + guides(color = "none")
  AC_plot <- AC_plot + guides(color = 'none')
  DLPFC_plot <- DLPFC_plot + guides(color = 'none')
  PCG_plot <- PCG_plot + guides(color = "none")
  tmt_plot <- tmt_plot + guides(color = guide_legend(override.aes = list(size = 7)))
  
  # Setting plot margins
  combined_plot <- (acetil_plot & theme(plot.margin = unit(c(1, 1, 1, 1), "mm"))) + 
    (AC_plot & theme(plot.margin = unit(c(1, 0.2, 1, 0.2), "mm"))) + 
    (DLPFC_plot & theme(plot.margin = unit(c(1, 0.2, 1, 0.2), "mm"))) + 
    (PCG_plot & theme(plot.margin = unit(c(1, 0.2, 1, 0.2), "mm"))) + 
    (tmt_plot & theme(plot.margin = unit(c(1, 0.2,1, 1, 1), "mm"))) + 
    plot_layout(guides = "collect", ncol = 5) & 
    theme(legend.position = "right")
  
  combined_plot
  
  
}

# Retrieving lists of pathways/modules with their corresponding genes
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/ALL_DATA_BASE.RDATA")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/new_genes_list/MODULES_LIST.RDATA")

df_weights_F14 <- cleared_df

# A custom vector of genes related to AD
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/F14_gene_vec")

# Adding the gene APOE to the g_was vector
g_was <- MODULES_LIST[1]
g_was=unlist(g_was)
g_was <- c(g_was,'APOE')
names(g_was)[74] <- 'AD_GWAS74'

# Genes that are often repeated among views
repeated_genes <- c("DNAJC6", "VGF")

# Threshold for filtering feature weights by absolute value
val_filter = 0.5
df_filtered <- filter(df_weights_F14, value > val_filter)

##################################################### Chosen pathways #############################################################################

# Retrieving the genes from the pathways lists

MAPK <- MODULES_LIST[grep("m7_mapk_metabolism",names(MODULES_LIST))]
MAPK = unlist(MAPK)
MAPK = unique(MAPK)
MAPK = intersect(df_filtered$feature,MAPK)

synapse <- MODULES_LIST[grep("m5_post_synaptic_density",names(MODULES_LIST))]
synapse = unlist(synapse)
synapse = unique(synapse)
synapse = intersect(df_filtered$feature,synapse)

prot_transport <-  MODULES_LIST['m8_protein_transport']
prot_transport <- unlist(prot_transport)
prot_transport <- unique(prot_transport)
prot_transport <- intersect(df_filtered$feature,prot_transport)

# Concatenating genes in a vector
combined_genes  = c(MAPK,prot_transport,synapse,g_was)

# Checking which genes are in more than one of the combined pathways
gene_counts <- table(unlist(combined_genes))
genes_in_more_than_one <- names(gene_counts[gene_counts > 1])
##################################################################################################################################################

# Adding the sign to the absolute weight values
df_weights_F14$real_values <- ifelse(df_weights_F14$sign == '-', df_weights_F14$value*-1, df_weights_F14$value*1)
df_weights_F14$category <- 'genes'

# Filtering by genes present in 'combined_genes' and by 'val_filter'
df_weights_F14x <- filter(df_weights_F14, feature %in% combined_genes)
df_weights_F14x <- filter(df_weights_F14x, value > val_filter)

# Updating the 'category' column values according to each pathway
df_weights_F14x <- df_weights_F14x %>%
  mutate(
    category = case_when(
      feature %in% g_was ~ 'AD GWAS',
      feature %in% synapse ~ 'Post-synaptic density',
      feature %in% prot_transport ~ 'Protein transport',
      feature %in% MAPK ~ 'MAPK metabolism',
      TRUE ~ 'genes'
    )
  )

# Filtering by genes not present in 'combined_genes' and by 'val_filter'
df_weights_F14 <- filter(df_weights_F14, !(feature %in% combined_genes) | value <= val_filter)

#' Recombining 'df_weights_F14x' and 'df_weights_F14' to restore the original data frame structure, 
#' but now with different values on the column 'category'
df_weights_F14 <- bind_rows(df_weights_F14, df_weights_F14x)

# Rearranging order
df_weights_F14$category <- factor(df_weights_F14$category, levels = c(
  'AD GWAS',
  'MAPK metabolism',
  'Post-synaptic density',
  'Protein transport',
  "genes"
))

# Display plot
jitter_plot_F14 <- jitter_ploting_weights(df_weights_F14,filtered_val = val_filter, 
                                          categories = c('AD GWAS' = 'black',
                                                         'MAPK metabolism' = "#99991EFF",
                                                         'Post-synaptic density' = "#CC0000FF",
                                                         'Protein transport' = "#631879FF"
                                          ))

ggsave("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Figuras_MOFA/pre_inkscape/weight/jitter/new/new_jitter_F14_v6.pdf",
       plot = jitter_plot_F14 , width = 500, height = 255, units = 'mm')
