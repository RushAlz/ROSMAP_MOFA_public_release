#########################################################
# DISPLAY JITTER PLOT WITH FEATURE WEIGHTS FOR FACTOR 8 #
#########################################################
library(MOFA2)        
library(dplyr)        
library(ggplot2)      
library(ggrepel)      
library(scales)       
library(ggthemes)     
library(ggeasy)       
library(ggsci)        

# Set working directory
setwd("/pastel/Github_scripts/OtherProjects/MOFA/scripts")
source("map_gene_to_function.R")

phenotypes = readRDS("/pastel/resources/data2share/basic_Apr2022_customJul2024.rds") 
rownames(phenotypes) = phenotypes$projid 

acetil_meta = readRDS("/pastel/resources/data2share/histoneAcetylation_H3K9ac.metadata.rds")
acetil_meta$gene = acetil_meta$closest_gene

# Load your existing MOFA data
load("/pastel/projects/MOFA/model/standart_50f.RDATA")
MOFAobject = mofa

###############################################################################
factor_i = 8
val_filter = 0.3

create_jitter <- function(factor_i, val_filter){
  
  loadings_Fi = get_weights(MOFAobject, scale = T, abs = F, factors = factor_i, as.data.frame = T)
  
  df_weights_Fi = loadings_Fi %>%
    filter(view %in% c("rna_AC","rna_DLPF","rna_PCG","tmt","acetil")) %>%
    mutate(gene = feature) %>%
    mutate(gene = case_when(
      view == "rna_AC" ~ gsub("(.*?)_(.*)_rna_AC","\\2",gene),
      view == "rna_DLPF" ~ gsub("(.*?)_(.*)_rna_DLPF","\\2",gene),
      view == "rna_PCG" ~ gsub("(.*?)_(.*)_rna_PCG","\\2",gene),
      view == "tmt" ~ gsub("(.*);(.*)","\\1",gsub("(.*?)_(.*)","\\2",gene)),
      view == "acetil" ~ acetil_meta[,c("probeid","gene")][match(gene, acetil_meta$probeid),"gene"])
    ) %>% 
    mutate(absolute_score = abs(value)) %>%
    mutate(view = gsub("rna_AC", "RNA (AC)", view) %>%
             gsub('rna_PCG', 'RNA (PCG)', .) %>%
             gsub('rna_DLPF', 'RNA (DLPFC)', .) %>%
             gsub('acetil', 'H3K9ac', .) %>%
             gsub('tmt', 'Protein', .) %>%
             gsub('metabol', 'Metabolites', .) %>%
             gsub('sn_RNA_norm', 'Cell propotions', .)) %>%
    arrange(-absolute_score)
  df_weights_Fi$view = factor(df_weights_Fi$view, levels = c("H3K9ac","RNA (AC)","RNA (DLPFC)","RNA (PCG)","Protein")) 
  
  # Threshold for filtering feature weights by absolute value
  final_Fi_df = df_weights_Fi %>% mutate(highlight = ifelse(absolute_score > val_filter, T, F))
  final_Fi_df_top_up = final_Fi_df %>% group_by(view) %>% top_n(5, wt = value) %>% pull(gene)
  final_Fi_df_top_down = final_Fi_df %>% group_by(view) %>% top_n(5, wt = -value) %>% pull(gene)
  final_Fi_df_top = unique(c(final_Fi_df_top_up, final_Fi_df_top_down))
  final_Fi_df$top_weights <- ifelse(final_Fi_df$gene %in% final_Fi_df_top & final_Fi_df$highlight, final_Fi_df$gene, NA)
  
  genes_to_plot = unique(final_Fi_df$top_weights[!is.na(final_Fi_df$top_weights)])
  
  gene_function = classify_gene_list(genes_to_plot, method = "go_db", verbose = TRUE)
  
  category_colors = setNames(ggsci::pal_d3()(length(get_category_definitions())), 
                             names(get_category_definitions()))
  category_colors["Other"] = "grey50"
  
  gene_categories = split(gene_function$Gene, gene_function$Category)
  
  # Create a function to map genes to categories, handling multiple category assignments
  get_gene_category <- function(gene_name) {
    for (category in names(gene_categories)) {
      if (gene_name %in% gene_categories[[category]]) {
        return(category)
      }
    }
    return("Other") # Default for genes not in the defined lists
  }
  
  library(scales)
  # 1a) Define a cubic transform
  power_trans <- function(p) {
    trans_new(
      name      = paste0("power_", p),
      transform = function(x) sign(x) * abs(x)^p,
      inverse   = function(u) sign(u) * abs(u)^(1/p)
    )
  }
  
  # Add category information to the data frame
  final_Fi_df$gene_category <- sapply(final_Fi_df$gene, get_gene_category)
  final_Fi_df$gene_category = factor(final_Fi_df$gene_category, levels = c(names(get_category_definitions()), "Other"))
  final_Fi_df$label_color <- category_colors[final_Fi_df$gene_category]
  final_Fi_df$alpha_val <- ifelse(final_Fi_df$gene %in% genes_to_plot, 1, 0.4)
  final_Fi_df[!final_Fi_df$highlight, "alpha_val"] = 0.4
  final_Fi_df[!final_Fi_df$highlight, "gene_category"] = "Other"
  final_Fi_df[!final_Fi_df$highlight, "label_color"] = "grey50"
  
  
  plott = ggplot(final_Fi_df, aes(x = view, y = value)) +
    geom_hline(yintercept = c(-val_filter, val_filter), linetype = "dashed", color = "grey") +
    geom_jitter(aes(color = gene_category, alpha = alpha_val),
                position = position_jitter(seed = 1), show.legend = T) +
    ggthemes::theme_tufte(base_family = "Helvetica") +
    scale_y_continuous(
      trans  = power_trans(p = 1.8),
      breaks = seq(-1, 1, by = 0.2),
      limits = c(-1.01, 1.01)
    ) +
    facet_wrap(~view, scales = "free_x", nrow = 1) +
    ggrepel::geom_text_repel(
      aes(label = top_weights, color = gene_category),
      size = 3, fontface = "italic",
      max.overlaps = 10, show.legend = F,
      position = position_jitter(seed = 1)
    ) +
    scale_color_manual(values = category_colors, name = "Gene Category", drop = FALSE) +
    theme(legend.position = "right",
          strip.text = element_text(size = 14),
          legend.key.size = unit(0.5, "cm"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8)) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
    ggeasy::easy_remove_axes(what = c("tick","text"), which = "x") +
    guides(alpha = "none") +
    labs(x = "Views", y = "Weights")
  
  return(plott)
}

fct_2_plot = c(1,2,3,4,8,12,14,26,42)
for(i in fct_2_plot){
  plott = create_jitter(factor_i = i, val_filter = 0.3)
  pdf(paste0("/pastel/Github_scripts/ROSMAP-MOFA/figures/F",i,"_jitter.pdf"), width = 10, height = 5)
  print(plott)
  dev.off()
}
