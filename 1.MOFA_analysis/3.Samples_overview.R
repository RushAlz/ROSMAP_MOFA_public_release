###########################
# Ploting Samples oveview #
###########################
library(MOFA2)
library(ggplot2)
library(dplyr)
library(ggsci)

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# Retrieving participants phenotypes data
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/meta.RDATA")
meta <- meta[(meta$projid %in% mofa@samples_metadata$sample) , ]
samples_metadata(mofa) <- meta

# MOFA function to do a tile plot showing the missing value structure of the input data
xx <- plot_data_overview(mofa)

# Renaming views
xx$data$view_label <- sub("rna_AC\nD=17309", 'RNA (AC)\n(n=721)', xx$data$view_label)
xx$data$view_label <- sub("rna_DLPF\nD=17309", 'RNA (DLPFC)\n(n=1210)', xx$data$view_label)
xx$data$view_label <- sub("rna_PCG\nD=17309", 'RNA (PCG)\n(n=656)', xx$data$view_label)
xx$data$view_label <- sub("sn_RNA_norm\nD=96", 'Cell Types\n(n=424)', xx$data$view_label)
xx$data$view_label <- sub("tmt\nD=8425", 'Proteomics\n(n=596)', xx$data$view_label)
xx$data$view_label <- sub("metabol\nD=667", 'Metabolomics\n(n=500)', xx$data$view_label)
xx$data$view_label <- sub("acetil\nD=26384", 'H3K9ac\n(n=632)', xx$data$view_label)

# Reordeing views
xx$data$view_label <- factor(xx$data$view_label, levels = rev(c('H3K9ac\n(n=632)','RNA (AC)\n(n=721)',
                                                                'RNA (DLPFC)\n(n=1210)','RNA (PCG)\n(n=656)',
                                                                'Proteomics\n(n=596)','Metabolomics\n(n=500)',
                                                                'Cell Types\n(n=424)')))

# Setting the colors for each view
s <- pal_lancet(alpha = 1)(7)
names(s) <- c('acetil','rna_AC','rna_DLPF','rna_PCG','tmt','metabol','sn_RNA_norm')
xx$plot_env$colors <- s
xx <- xx +
  scale_fill_manual(values = xx$plot_env$colors)

# Renaming title
xx$data$group_label <- 'Samples Overview'

# Renaming X-axis lable
xx$labels$x <- 'Samples (n=1358)'

# Adjusting text size
xx$theme$text$size <- 14
xx <- xx +  theme(axis.text.y = element_text(size = 8))
xx <- xx +  theme(axis.title.x = element_text(size = 8))
xx <- xx + labs(x = "Samples")

xx
