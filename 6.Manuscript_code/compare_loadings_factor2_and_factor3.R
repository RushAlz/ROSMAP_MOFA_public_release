library(ggpubr)
library(MOFAdata)
library(MOFA2)
library(tidyverse)
setwd("/pastel/Github_scripts/ROSMAP-MOFA")

# C2: curated gene sets from online pathway databases, publications in PubMed, and knowledge of domain experts.
data("MSigDB_v6.0_C2_human") 
# C5: extracted from the Gene Ontology data.base
data("MSigDB_v6.0_C5_human") 

load("/pastel/projects/MOFA/model/standart_50f.RDATA")
phenotypes = readRDS("/pastel/resources/data2share/basic_Apr2022_selected_list_Jul2024.rds") 
phenotypes$projid = rownames(phenotypes)
load("/pastel/resources/data2share/pheno_list_Jul2024.RData")

sample_metadata <- data.frame(
  projid = samples_names(mofa)[[1]]
) %>% left_join(phenotypes)
sample_metadata$sample = sample_metadata$projid
samples_metadata(mofa) <- sample_metadata

loadings = bind_rows(
  get_weights(mofa, as.data.frame = T, scale = T, factors = c(2)),
  get_weights(mofa, as.data.frame = T, scale = T, factors = c(3))
)
# summary(loadings$value)
# loadings %>% arrange(-value) %>% head()

l = loadings %>% filter(factor %in% c("Factor2","Factor3")) %>%
  # filter(abs(value) > 0.001) %>%
  mutate(feature = paste0(view, "|", feature)) %>% 
  pivot_wider(id_cols = feature, names_from = c(factor), values_from = value) %>%
  mutate(view = gsub("(.*?)\\|(.*)","\\1",feature)) %>%
  group_by(view) %>% 
  na.omit() %>% arrange(-Factor2) %>%
  mutate(view = gsub("rna_AC", "RNA (AC)", view) %>%
           gsub('rna_PCG', 'RNA (PCG)', .) %>%
           gsub('rna_DLPF', 'RNA (DLPFC)', .) %>%
           gsub('acetil', 'H3K9ac', .) %>%
           gsub('tmt', 'Protein', .) %>%
           gsub('metabol', 'Metabolites', .) %>%
           gsub('sn_RNA_norm', 'Cell propotions', .))
  

s <- ggsci::pal_lancet(alpha = 1)(7)
names(s) <- c("H3K9ac", "RNA (AC)", "RNA (DLPFC)", "RNA (PCG)", "Protein", "Metabolites", "Cell propotions")
# names(s) <- c('acetil','rna_AC','rna_DLPF','rna_PCG','tmt','metabol','sn_RNA_norm')

l %>% reframe(cor = cor(Factor2,Factor3)) 

plot_a = ggplot(l, aes(x = Factor2, y = Factor3)) + 
  geom_point(aes(color = view), show.legend = F) + 
  geom_smooth(method = "lm", se = FALSE, show.legend = F) +
  stat_cor(label.y.npc="top", label.x.npc = "left", method = "pearson") +
  theme_minimal() + 
  scale_color_manual(values = s) +
  facet_wrap(~view, ncol = 1) +
  labs(title = "Feature weights comparison between Factors 2 and 3\n(scaled by factor across views)", 
       x = "Factor 2", y = "Factor 3", color = "View")
plot_a

df = l %>% filter(view == "RNA (DLPFC)")

library(mclust)
# Fit Gaussian Mixture Model (GMM) with 2 clusters
gmm_result <- Mclust(df[,c("Factor2","Factor3")], G = 2)
# Get cluster probabilities
df$prob1 <- gmm_result$z[,1]  # Probability of belonging to cluster 1
df$prob2 <- gmm_result$z[,2]  # Probability of belonging to cluster 2
# Assign cluster with the highest probability (optional)
df$cluster <- as.factor(gmm_result$classification)

plot_b = ggplot(df, aes(x = Factor2, y = Factor3)) + 
  geom_point(color = "grey20", size = 1, alpha = 0.8) +  # All dots are black
  stat_ellipse(data = df, aes(x = Factor2, y = Factor3, group = cluster, fill = cluster, color = cluster, alpha = cluster),
               geom = "polygon",  linewidth = 1, level = 0.95) +  # Confidence ellipses
  scale_color_manual(values = c("blue", "red")) +  # Adjust colors
  scale_fill_manual(values = c("blue", "red")) +  # Adjust colors
  scale_alpha_manual(values = c("1" = 0.1, "2" = 0.3)) +  # Adjust transparency
  theme_minimal() +
  labs(title = "Features weights between Factor2 and Factor3 in DLPFC RNA",
    subtitle = "Gaussian Mixture Model Clustering with Confidence Ellipses")
plot_b

df$gene_symbol = gsub("(.*?)\\|(.*?)\\.(.*)_(.*?)_(.*)","\\4",df$feature)

# Define probability threshold for overlapping assignment (e.g., > 0.3)
threshold <- 0.05
# Create subsets for overlapping clusters
df_cluster1 <- df %>% filter(prob1 > threshold)
df_cluster2 <- df %>% filter(prob2 > threshold)

my_list_of_genes = list(cluster1_factor2 = setNames(df_cluster1$Factor2, df_cluster1$gene_symbol),
     cluster1_factor3 = setNames(df_cluster1$Factor3, df_cluster1$gene_symbol),
     cluster2_factor2 = setNames(df_cluster2$Factor2, df_cluster2$gene_symbol),
     cluster2_factor3 = setNames(df_cluster2$Factor3, df_cluster2$gene_symbol))

run_GSEA2 <- function(my_list_of_genes, pathways_list = GO_FEA, threads = 2){
  # Libraries
  if (!require("clusterProfiler")) BiocManager::install("clusterProfiler"); library(clusterProfiler)
  if (!require("foreach")) install.packages("foreach"); library(foreach)
  if (!require("doParallel")) install.packages("doParallel"); library(doParallel)
  doParallel::registerDoParallel(cores = threads)
  
  if(!is.list(my_list_of_genes)){
    my_list_of_genes = list(GeneList = my_list_of_genes)
  }
  
  enrich_res_df = foreach(i = 1:length(my_list_of_genes), .combine=rbind) %do% {
    GeneSet = names(my_list_of_genes)[i]
    set.seed(2020)
    gsea_res <- GSEA(
      geneList = sort(my_list_of_genes[[GeneSet]],decreasing = T), # Ordered ranked gene list
      minGSSize = 25, # Minimum gene set size
      maxGSSize = Inf, # Maximum gene set set
      pvalueCutoff = 0.05, # p-value cutoff
      eps = 0, # Boundary for calculating the p value
      seed = TRUE, # Set seed to make results reproducible
      pAdjustMethod = "BH", # Benjamini-Hochberg correction
      nPermSimple = 10000, # Number of permutations
      TERM2GENE = dplyr::select(
        pathways_list,
        gs_name,
        gene_symbol
      )
    )
    gsea_res@result$core_enrichment = NULL
    gsea_res@result$GeneSet = GeneSet
    gsea_res@result
  }
  return(enrich_res_df)  
}

library(rrvgo)
library(org.Hs.eg.db)
library(msigdbr)
library(tidyverse)
library(sigPathway)
library(purrr)

GO_msigdbr = list()
## GENE ONTOLOGY (MOLECULAR FUNCTION)
GO_msigdbr[["MF"]] <- msigdbr(species = "human", category = "C5", subcategory = "GO:MF")
## GENE ONTOLOGY (BILOGICAL PROCESS)
GO_msigdbr[["BP"]] <- msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
## GENE ONTOLOGY (CELLULAR COMPONENT)
GO_msigdbr[["CC"]] <- msigdbr(species = "human", category = "C5", subcategory = "GO:CC")
ALL_DATA_BASE = map_df(GO_msigdbr,bind_rows)

gsea_results = run_GSEA2(my_list_of_genes, pathways_list = ALL_DATA_BASE, threads = 6)
top_by_geneSet = gsea_results %>% group_by(GeneSet) %>% top_n(10, -log10(p.adjust)) 
# Filter by top terms
gsea_results_filtered = gsea_results %>% filter(ID %in% unique(top_by_geneSet$ID)) %>%
  dplyr::select(GeneSet, Description, NES, p.adjust) %>%
  mutate(Description = stringr::str_trunc(Description,40))

plot_c = ggplot(gsea_results_filtered, aes(x = GeneSet, y = Description)) +
  geom_point(aes(size = -log10(p.adjust), color = NES)) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +  
  theme_minimal() +
  # scale_y_discrete(labels = stringr::str_trunc(n = 40)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.95)) +
  labs(title = "GSEA Results: Top 10 Terms per Query",
       x = "Query",
       y = "Enriched Term",
       color = "NES",
       size = "-log10(P-value)") 
plot_c

pdf("figures/GSEA_Factor3_DLPFC.pdf", width = 12, height = 12)
ggarrange(plot_a, 
          ggarrange(plot_b, plot_c, ncol = 1, nrow = 2, heights = c(1,2)), ncol = 2, nrow = 1, widths = c(0.7,1))
dev.off()
png("figures/GSEA_Factor3_DLPFC.png", width = 12, height = 12, res = 300, units = "in")
ggarrange(plot_a, 
          ggarrange(plot_b, plot_c, ncol = 1, nrow = 2, heights = c(1,2)), ncol = 2, nrow = 1, widths = c(0.7,1))
dev.off()
