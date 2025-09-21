##############################################
# SCRIPT TO DISPLAY THE TOP FEATURE WEIGHTS #
##############################################
library(MOFA2)
library(readxl)
library(dplyr)
library(patchwork)
library(ggplot2)
library(stringr)

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# Choosing the factor and the number of features to be displayed
factors = 4
nFeatures = 5

#################################### ACETYL PLOTS ###########################################

# Retrieving metadata for histone acetylation
metadata <- readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/histoneAcetylation_H3K9ac.metadata.rds")

# MOFA function that retrieves all acetylations peaks weights
enrich_H3K9ac = plot_top_weights(mofa, view = 'acetil', factor = factors,
                                 nfeatures = nrow(mofa@data[['acetil']]$group1))$data


# Preparing 'metadata' to merge with 'enrich_H3K9ac'. The goal is to map 'enrich_H3K9ac' peaks to the closest genes in 'metadata'
metadata <- metadata[,c('probeid','closest_gene')]
metadata <- subset(metadata, probeid %in% enrich_H3K9ac$feature)
colnames(metadata)[colnames(metadata) == "probeid"] <- "feature"

# Merging
merge_metadata <- merge(metadata, enrich_H3K9ac, by = 'feature')
merge_metadata <-merge_metadata[,-1]
colnames(merge_metadata)[1] <- 'feature'

# Selecting the 100 highest absolute weights
merge_metadata <- merge_metadata %>% 
  arrange(desc(value)) %>% 
  head(100)

# Clearing up gene names
merge_metadata$feature <- gsub('\\..*','',merge_metadata$feature)
merge_metadata$feature <- gsub('-.*','',merge_metadata$feature)

# Some peaks share the same closest gene; so here we are keeping only one entry per gene (the one with higher absolute loading)
for (i in unique(merge_metadata$feature)){
  
  same_names = subset(merge_metadata, feature == i)
  if (any(same_names$sign == '+')){
    
    max_value <- unique(max(same_names$value))
    exclude_names <- subset(same_names, value !=max_value)
    merge_metadata <- anti_join(merge_metadata, exclude_names, by = c("feature", 'value', 'sign'))
  }
  
  else{
    
    max_value <- unique(min(same_names$value))
    exclude_names <- subset(same_names, value !=max_value)
    merge_metadata <- anti_join(merge_metadata, exclude_names, by = c("feature", 'value', 'sign'))
  }
}

# Ordering the n top positive features to be displayed 
filter_pos <- filter(merge_metadata, sign == '+')
filter_pos <- filter_pos[1:nFeatures,]
filter_pos$feature <- factor(filter_pos$feature, levels = rev(filter_pos$feature))
filter_pos <- filter_pos[order(filter_pos$value, decreasing = TRUE), ]

# Ordering the n top positive features to be displayed 
filter_neg <- filter(merge_metadata, sign == '-')
filter_neg <- filter_neg[1:nFeatures,]
filter_neg$feature <- factor(filter_neg$feature, levels = rev(filter_neg$feature))
filter_neg <- filter_neg[order(filter_neg$value, decreasing = TRUE), ]
####################################### RNA AC PLOTS ############################################

# MOFA function that retrieves the positive RNA (AC) features weights
enrich_pos_AC = plot_top_weights(mofa, view = 'rna_AC', factor = factors, nfeatures = nFeatures, sign = 'positive')
enrich_pos_AC <- enrich_pos_AC$data

# Preparing features names
enrich_pos_AC$feature <- gsub('.*_(.*)_(.*)_(.*)', '\\1', enrich_pos_AC$feature)
enrich_pos_AC$feature <- gsub('\\.','',enrich_pos_AC$feature)
enrich_pos_AC$feature <- paste(".", enrich_pos_AC$feature, sep = "")

# Ordering the top positive features to be displayed 
enrich_pos_AC$feature <- factor(enrich_pos_AC$feature, levels = rev(enrich_pos_AC$feature))
enrich_pos_AC <- enrich_pos_AC[order(enrich_pos_AC$value, decreasing = TRUE), ]

# MOFA function that retrieves the negative RNA (AC) features weights
enrich_neg_AC = plot_top_weights(mofa, view = 'rna_AC', factor = factors, nfeatures = nFeatures, sign = 'negative')
enrich_neg_AC <- enrich_neg_AC$data

# Preparing features names
enrich_neg_AC$feature <- gsub('.*_(.*)_(.*)_(.*)', '\\1', enrich_neg_AC$feature)
enrich_neg_AC$feature <- gsub('\\.','',enrich_neg_AC$feature)
enrich_neg_AC$feature <- paste(".", enrich_neg_AC$feature, sep = "")

# Ordering the top negative features to be displayed
enrich_neg_AC$feature <- factor(enrich_neg_AC$feature, levels = rev(enrich_neg_AC$feature))
enrich_neg_AC <- enrich_neg_AC[order(enrich_neg_AC$value, decreasing = T), ]
####################################### RNA DLPFC PLOTS ############################################

# MOFA function that retrieves the positive RNA (DLPFC) features weights
enrich_pos_DLPFC = plot_top_weights(mofa, view = 'rna_DLPF', factor = factors, nfeatures = nFeatures, sign = 'positive')
enrich_pos_DLPFC <- enrich_pos_DLPFC$data

# Preparing features names
enrich_pos_DLPFC$feature <- gsub('.*_(.*)_(.*)_(.*)', '\\1', enrich_pos_DLPFC$feature)
enrich_pos_DLPFC$feature <- gsub('\\.','',enrich_pos_DLPFC$feature)
enrich_pos_DLPFC$feature <- paste(";", enrich_pos_DLPFC$feature, sep = "")

# Ordering the top positive features to be displayed
enrich_pos_DLPFC$feature <- factor(enrich_pos_DLPFC$feature, levels = rev(enrich_pos_DLPFC$feature))
enrich_pos_DLPFC <- enrich_pos_DLPFC[order(enrich_pos_DLPFC$value, decreasing = TRUE), ]

# MOFA function that retrieves the negative RNA (DLPFC) features weights
enrich_neg_DLPFC = plot_top_weights(mofa, view = 'rna_DLPF', factor = factors, nfeatures = nFeatures, sign = 'negative')
enrich_neg_DLPFC <- enrich_neg_DLPFC$data

# Preparing features names
enrich_neg_DLPFC$feature <- gsub('.*_(.*)_(.*)_(.*)', '\\1', enrich_neg_DLPFC$feature)
enrich_neg_DLPFC$feature <- gsub('\\.','',enrich_neg_DLPFC$feature)
enrich_neg_DLPFC$feature <- paste(";", enrich_neg_DLPFC$feature, sep = "")

# Ordering the top negative features to be displayed
enrich_neg_DLPFC$feature <- factor(enrich_neg_DLPFC$feature, levels = rev(enrich_neg_DLPFC$feature))
enrich_neg_DLPFC <- enrich_neg_DLPFC[order(enrich_neg_DLPFC$value, decreasing = T), ]
####################################### RNA PCG PLOTS ############################################

# MOFA function that retrieves the positive RNA (PCG) features weights
enrich_pos_PCG = plot_top_weights(mofa, view = 'rna_PCG', factor = factors, nfeatures = nFeatures, sign = 'positive')
enrich_pos_PCG <- enrich_pos_PCG$data

# Preparing features names
enrich_pos_PCG$feature <- gsub('.*_(.*)_(.*)_(.*)', '\\1', enrich_pos_PCG$feature)
enrich_pos_PCG$feature <- gsub('\\.','',enrich_pos_PCG$feature)
enrich_pos_PCG$feature <- paste(",", enrich_pos_PCG$feature, sep = "")

# Ordering the top positive features to be displayed
enrich_pos_PCG$feature <- factor(enrich_pos_PCG$feature, levels = rev(enrich_pos_PCG$feature))
enrich_pos_PCG <- enrich_pos_PCG[order(enrich_pos_PCG$value, decreasing = TRUE), ]

# MOFA function that retrieves the negative RNA (PCG) features weights
enrich_neg_PCG = plot_top_weights(mofa, view = 'rna_PCG', factor = factors, nfeatures = nFeatures, sign = 'negative')
enrich_neg_PCG <- enrich_neg_PCG$data

# Preparing features names
enrich_neg_PCG$feature <- gsub('.*_(.*)_(.*)_(.*)', '\\1', enrich_neg_PCG$feature)
enrich_neg_PCG$feature <- gsub('\\.','',enrich_neg_PCG$feature)
enrich_neg_PCG$feature <- paste(",", enrich_neg_PCG$feature, sep = "")

# Ordering the top negative features to be displayed
enrich_neg_PCG$feature <- factor(enrich_neg_PCG$feature, levels = rev(enrich_neg_PCG$feature))
enrich_neg_PCG <- enrich_neg_PCG[order(enrich_neg_PCG$value, decreasing = T), ]
###################################### TMT PLOTS ####################################################

# MOFA function that retrieves the positive proteomics features weights
enrich_pos_TMT = plot_top_weights(mofa, view = 'tmt', factor = factors, nfeatures = nFeatures, sign = 'positive')
enrich_pos_TMT <- enrich_pos_TMT$data

# Preparing features names
enrich_pos_TMT$feature <- gsub('.*-','', enrich_pos_TMT$feature)
enrich_pos_TMT$feature <- gsub('.*_','',enrich_pos_TMT$feature )
enrich_pos_TMT$feature <- paste("*", enrich_pos_TMT$feature, sep = "")

# Ordering the top positive features to be displayed
enrich_pos_TMT$feature <- factor(enrich_pos_TMT$feature, levels = rev(enrich_pos_TMT$feature))
enrich_pos_TMT <- enrich_pos_TMT[order(enrich_pos_TMT$value, decreasing = TRUE), ]

# MOFA function that retrieves the negative proteomics features weights
enrich_neg_TMT = plot_top_weights(mofa, view = 'tmt', factor = factors, nfeatures = nFeatures, sign = 'negative')
enrich_neg_TMT <- enrich_neg_TMT$data

# Preparing features names
enrich_neg_TMT$feature <- gsub('.*-','', enrich_neg_TMT$feature)
enrich_neg_TMT$feature <- gsub('.*_','',enrich_neg_TMT$feature)
enrich_neg_TMT$feature <- paste("*", enrich_neg_TMT$feature, sep = "")

# Ordering the top negative features to be displayed
enrich_neg_TMT$feature <- factor(enrich_neg_TMT$feature, levels = rev(enrich_neg_TMT$feature))
enrich_neg_TMT <- enrich_neg_TMT[order(enrich_neg_TMT$value, decreasing = TRUE), ]
###################################### METABOLITES ##################################################################

# MOFA function that retrieves all metabolomics features weights
enrich_all_met = plot_top_weights(mofa, view = 'metabol', factor = factors, nfeatures = nrow(mofa@data$metabol$group1))
enrich_all_met <- enrich_all_met$data

# Cleaning up the names
enrich_all_met$feature <- gsub('^X', '', enrich_all_met$feature)

# Preparing a data frame with more readable names to merge with 'enrich_all_met'
correct_names <- read_excel("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/Metabolites/s1_nomes_IDs.xlsx")
colnames(correct_names) <- c('correct_name', 'super_pathway', 'sub_pathway', 'comp_id', 'platform', 'chemical_id', 'RI', 'mass', 'pubchem', 
                             'cas', 'KEGG', 'HMDB', 'comp_IDstr')
correct_names <- correct_names[-1,]
correct_names$feature <- correct_names$correct_name
correct_names <- correct_names[,c(1,14)]
correct_names$feature <- gsub('[^A-Za-z0-9]','.',correct_names$feature)
correct_names$feature <- gsub('^X','',correct_names$feature)

# Merging
merge_metabol <- merge(correct_names, enrich_all_met, by = 'feature')
merge_metabol <- merge_metabol[,-1]
colnames(merge_metabol)[1] <- 'feature'

# Ordering the n top positive features to be displayed 
enrich_pos_metabol <- filter(merge_metabol,sign == '+')
enrich_pos_metabol <- enrich_pos_metabol[order(enrich_pos_metabol$value, decreasing = T), ]
enrich_pos_metabol <- enrich_pos_metabol[(1:nFeatures) ,]
enrich_pos_metabol$feature <- factor(enrich_pos_metabol$feature, levels = rev(enrich_pos_metabol$feature))
enrich_pos_metabol <- enrich_pos_metabol[order(enrich_pos_metabol$value, decreasing = T), ]

# Ordering the n top negative features to be displayed
enrich_neg_metabol <- filter(merge_metabol,sign == '-')
enrich_neg_metabol <- enrich_neg_metabol[order(enrich_neg_metabol$value, decreasing = T), ]
enrich_neg_metabol <- enrich_neg_metabol[(1:nFeatures) ,]
enrich_neg_metabol$feature <- factor(enrich_neg_metabol$feature, levels = rev(enrich_neg_metabol$feature))
enrich_neg_metabol <- enrich_neg_metabol[order(enrich_neg_metabol$value, decreasing = T), ]
######################################## CELL TYPE #######################################################################

# MOFA function that retrieves the positive cell types features weights
enrich_pos_cell_type = plot_top_weights(mofa, view = 'sn_RNA_norm', factor = factors, nfeatures = nFeatures, sign = 'positive')
enrich_pos_cell_type <- enrich_pos_cell_type$data

# Preparing features names
enrich_pos_cell_type$feature <- toupper(enrich_pos_cell_type$feature)

# Ordering the top positive features to be displayed
enrich_pos_cell_type$feature <- factor(enrich_pos_cell_type$feature, levels = rev(enrich_pos_cell_type$feature))
enrich_pos_cell_type <- enrich_pos_cell_type[order(enrich_pos_cell_type$value, decreasing = TRUE), ]

# MOFA function that retrieves the negative cell types features weights
enrich_neg_cell_type = plot_top_weights(mofa, view = 'sn_RNA_norm', factor = factors, nfeatures = nFeatures, sign = 'negative')
enrich_neg_cell_type <- enrich_neg_cell_type$data

# Preparing features names
enrich_neg_cell_type$feature <- toupper(enrich_neg_cell_type$feature)

# Ordering the top negative features to be displayed
enrich_neg_cell_type$feature <- factor(enrich_neg_cell_type$feature, levels = rev(enrich_neg_cell_type$feature))
enrich_neg_cell_type <- enrich_neg_cell_type[order(enrich_neg_cell_type$value, decreasing = T), ]

####################################################################################################

# Binding all data.frames
all_df_weight <- rbind(filter_pos, enrich_pos_AC,enrich_pos_DLPFC,enrich_pos_PCG,enrich_pos_TMT,
                       enrich_pos_metabol,enrich_pos_cell_type,filter_neg,enrich_neg_AC,
                       enrich_neg_DLPFC,enrich_neg_PCG,enrich_neg_TMT,enrich_neg_metabol,
                       enrich_neg_cell_type)


all_df_weight$view <- factor(all_df_weight$view, levels = c("acetil",'rna_AC','rna_DLPF', 'rna_PCG',
                                                            'tmt', 'metabol', 'sn_RNA_norm'))
all_df_weight$sign <- factor(all_df_weight$sign, levels = c("+","-"))

all_df_weight$feature <- factor(all_df_weight$feature, ordered = TRUE)


Y_lab_ftsize = 10
Y_txt_ftsize = 7
X_lab_ftsize = 10
X_txt_ftsize = 8
pos_neg_size = 10

expnd = c(0,0)
marg = c(0,0,0,0)

# Ploting
all_df_weight_plot <- ggplot(all_df_weight, aes(x = value, y = feature)) +
  geom_segment(aes(xend = 0, yend = feature), color = "black", stat = 'identity', linewidth =0.5) +
  geom_point(aes(x = value, y = feature), color = "black", size = 0.5, stat = 'identity') +
  labs(title = 'Top Features For Each View', x = "Weight", y = "") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.margin = margin(marg),
        axis.title.y = element_text(size = Y_lab_ftsize, colour = 'black'),
        axis.text.y = element_text(size = Y_txt_ftsize, colour = 'black'),
        axis.title.x = element_text(size = X_lab_ftsize, colour = 'black'),
        axis.text.x = element_text(size = X_txt_ftsize, colour = 'black'),
        panel.grid.major.y = element_line(color = "gray", size = 0.1, linetype = "solid"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.02),
                     expand = expansion(add =expnd)) +
  facet_wrap(~view ~ sign, scales = 'free_y', ncol = 2, strip.position = 'left')

all_df_weight_plot