#########################################
# SCRIPT TO DISPLAY SAMPLES CLUSTERING #
########################################
library(MOFA2)
library(reshape2)
library(stringr)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyverse)
library(data.table)
library(doParallel)
library(foreach)

get_overlap_stats = function(df, column_i, column_j){
  library(GeneOverlap)
  go.obj.df = data.frame()
  for(mod_i in unique(speakeasy_cluster_df_combined[[column_i]])){
    for(mod_j in unique(speakeasy_cluster_df_combined[[column_j]])){
      donors_mod_i = speakeasy_cluster_df_combined$projid[speakeasy_cluster_df_combined[[column_i]] == mod_i]
      donors_mod_j = speakeasy_cluster_df_combined$projid[speakeasy_cluster_df_combined[[column_j]] == mod_j]
      go.obj <- newGeneOverlap(donors_mod_i,donors_mod_j, genome.size=length(unique(speakeasy_cluster_df_combined$projid)))
      go.obj <- testGeneOverlap(go.obj)
      jaccard_df = data.frame(mod_i = mod_i, mod_j = mod_j, jaccard = getJaccard(go.obj), pval = getPval(go.obj))
      go.obj.df = bind_rows(go.obj.df, jaccard_df)
    }
  }
  return(go.obj.df)
}

# Function to calculate 1000 permutations of the cluster labels in each phenotypes to get the null distribution and then calculate pvalues
calc_prs_expression_permutation <- function(se_cluster_pheno, 
                                            columns_to_analyze, 
                                            stats_to_summarise = "mean", 
                                            nperm = 1000, 
                                            make_plot = F){
  
  if(stats_to_summarise == "sum"){
    se_cluster_pheno[,columns_to_analyze,drop=F] = se_cluster_pheno[,columns_to_analyze,drop=F] %>% mutate_if(is.character,as.numeric) %>% mutate_if(is.factor,as.numeric) %>% mutate_if(is.logical,as.numeric)
  }
  
  # Transform the data into long format. Columns projid, cluster, pheno, value 
  se_cluster_pheno_long = se_cluster_pheno %>%
    select(projid, cluster, all_of(columns_to_analyze)) %>%
    pivot_longer(cols = all_of(columns_to_analyze),
                 names_to = "pheno",
                 values_to = "value",
                 values_transform = list(value = as.numeric)) %>% distinct()
  
  observed_stats = se_cluster_pheno %>%
    group_by(cluster) %>%
    summarise(across(all_of(columns_to_analyze), (!!sym(stats_to_summarise)), na.rm = T)) %>%
    pivot_longer(cols = all_of(columns_to_analyze), names_to = "pheno", values_to = "value")
  
  # Manually calculate the permutations
  permuted_stats = data.frame()
  
  n.cores <- 4
  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  #check cluster definition (optional)
  print(my.cluster)
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  #how many workers are available? (optional)
  foreach::getDoParWorkers()
  
  gc()
  set.seed(123)
  permuted_stats <- foreach(i = 1:nperm,
                            .packages = 'tidyverse',
                            .combine = bind_rows) %dopar% {
                              # Randomly permute cluster labels for each phenotype
                              se_cluster_perm = se_cluster_pheno %>% 
                                mutate(cluster = sample(cluster)) %>% ungroup()
                              
                              # calculate stats for permuted data
                              se_cluster_stats = se_cluster_perm %>% 
                                group_by(cluster) %>% 
                                summarise(across(all_of(columns_to_analyze), (!!sym(stats_to_summarise)), na.rm = T)) %>%
                                pivot_longer(cols = all_of(columns_to_analyze), names_to = "pheno", values_to = "value") %>%
                                mutate(permutation = as.character(i))
                              
                              # Append permuted stats to the permuted_stats dataframe
                              return(se_cluster_stats)
                            }
  
  parallel::stopCluster(cl = my.cluster)
  gc()
  
  # Get a P-value for the observed vs permuted stats
  stats_centered = observed_stats %>% mutate(permutation = "observed") %>%
    bind_rows(permuted_stats) %>%
    group_by(pheno, cluster) %>%
    mutate(value_centered = as.numeric(scale(value, center = T, scale = T)))
  
  perm_pvalues = stats_centered %>% 
    filter(permutation != "observed") %>%
    rename(null_value = value_centered) %>%
    left_join(stats_centered %>%
                filter(permutation == "observed") %>%
                rename(obs_value = value_centered),
              by = c("cluster","pheno"), suffix = c("_null","_obs")) %>%
    group_by(pheno,cluster) %>%
    summarise(estimate = mean(null_value),
              num_over_obs = sum( abs(null_value)>=abs(obs_value) ), 
              pval = sum( abs(null_value)>=abs(obs_value) ) / nperm)
  
  if(make_plot){
    null_dist_plot = permuted_stats %>% filter(pheno %in% columns_to_analyze) %>%
      ggplot(aes(value)) +
      geom_histogram() +
      geom_vline(data = observed_stats %>% filter(pheno %in% columns_to_analyze), 
                 aes(xintercept = value), color = "red") +
      facet_wrap(pheno~cluster, scales = "free", ncol = length(unique(permuted_stats$cluster))) + theme_classic()
  }else{
    null_dist_plot = NULL
  }
  
  return(list(perm_stats = perm_pvalues, 
              plot = null_dist_plot))
}

make_phenotype_plot = function(speakeasy_cluster_df_filtered, selected_factors = paste0("Factor",c(1, 2, 3, 4, 8, 12, 14, 26, 42))){
  # Add phenotypes
  se_cluster_pheno = speakeasy_cluster_df_filtered %>% 
    left_join(phenotype_dt %>% rownames_to_column("projid")) %>%
    left_join(mofa_factors %>% rownames_to_column("projid"))
  
  # Order clusters by mean age_death
  cluster_mean_age_death = se_cluster_pheno %>% group_by(cluster) %>% reframe(mean_age_death = mean(age_death,na.rm=T))
  se_cluster_pheno$cluster = factor(as.character(se_cluster_pheno$cluster),
                                    levels = cluster_mean_age_death$cluster[order(cluster_mean_age_death$mean_age_death)])
  
  ####################################################################################
  se_cluster_pheno$fsex = ifelse(as.numeric(se_cluster_pheno$msex)==1,1,0)
  
  df_mean_factor = se_cluster_pheno %>% 
    pivot_longer(cols = colnames(mofa_factors),
                 names_to = "pheno", values_to = "value") %>% 
    group_by(pheno, cluster) %>% 
    summarise(mean_pheno = mean(value, na.rm = T), n = n()) %>% 
    ungroup() %>% na.omit() %>% filter(pheno %in% selected_factors)
  df_mean_factor_perm = calc_prs_expression_permutation(se_cluster_pheno, 
                                                        selected_factors, 
                                                        nperm = 1000, 
                                                        stats_to_summarise = "mean", 
                                                        make_plot = T)
  df_mean_factor = df_mean_factor %>% 
    left_join(df_mean_factor_perm$perm_stats, by = c("pheno","cluster"))
  
  df_pheno_mean = se_cluster_pheno %>%
    pivot_longer(cols = c("age_death","educ",names(pheno_list)[pheno_list=="ordinal"]),
                 names_to = "pheno", values_to = "value") %>%
    group_by(pheno, cluster) %>%
    summarise(mean_pheno = mean(value, na.rm = T), n = n()) %>%
    ungroup()
  df_pheno_mean_perm = calc_prs_expression_permutation(se_cluster_pheno, 
                                                       c("age_death","educ",names(pheno_list)[pheno_list=="ordinal"]), 
                                                       nperm = 1000, 
                                                       stats_to_summarise = "mean", 
                                                       make_plot = T)
  df_pheno_mean = df_pheno_mean %>% 
    left_join(df_pheno_mean_perm$perm_stats,by = c("pheno","cluster"))
  
  df_pheno_median = se_cluster_pheno %>% 
    pivot_longer(cols = c(names(pheno_list)[pheno_list=="gaussian"]), 
                 names_to = "pheno", values_to = "value") %>% 
    group_by(pheno, cluster) %>% 
    summarise(median_pheno = median(value, na.rm = T), n = n()) %>% 
    ungroup() 
  df_pheno_median_perm = calc_prs_expression_permutation(se_cluster_pheno,
                                                         names(pheno_list)[pheno_list=="gaussian"], 
                                                         nperm = 1000, 
                                                         stats_to_summarise = "median", 
                                                         make_plot = T)
  df_pheno_median = df_pheno_median %>% 
    left_join(df_pheno_median_perm$perm_stats, by = c("pheno","cluster"))
  
  df_pheno_count = se_cluster_pheno %>% 
    pivot_longer(cols = c("fsex",names(pheno_list)[pheno_list=="binomial"]), 
                 names_to = "pheno", values_to = "value", values_transform = as.numeric) %>% 
    dplyr::select(pheno,cluster,value) %>% 
    group_by(pheno, cluster) %>%
    summarise(sum_pheno = sum(value, na.rm = T), n = n()) %>% 
    mutate(frac_pheno = sum_pheno/n) %>%
    ungroup()
  df_pheno_count_perm = calc_prs_expression_permutation(se_cluster_pheno,
                                                        columns_to_analyze = c("fsex",names(pheno_list)[pheno_list=="binomial"]), 
                                                        nperm = 1000, 
                                                        stats_to_summarise = "sum", 
                                                        make_plot = T)
  df_pheno_count = df_pheno_count %>% 
    left_join(df_pheno_count_perm$perm_stats, by = c("pheno","cluster"))
  
  df_pheno_all = bind_rows(
    df_mean_factor %>% dplyr::rename("value"="mean_pheno") %>% dplyr::select(pheno,cluster,value,estimate,pval),
    df_pheno_mean %>% dplyr::rename("value"="mean_pheno") %>% dplyr::select(pheno,cluster,value,estimate,pval),
    df_pheno_count %>% dplyr::rename("value"="frac_pheno") %>% dplyr::select(pheno,cluster,value,estimate,pval),
    df_pheno_median %>% dplyr::rename("value"="median_pheno") %>% dplyr::select(pheno,cluster,value,estimate,pval))
  df_pheno_all$padj = p.adjust(df_pheno_all$pval, method = "fdr")
  
  ################################################################################################
  df_pheno_mat = df_pheno_all %>% dplyr::select(pheno,cluster,value) %>%
    pivot_wider(names_from = "pheno", values_from = "value") %>%
    as.data.frame()
  rownames(df_pheno_mat) = df_pheno_mat$cluster
  df_pheno_mat$cluster = NULL
  
  df_pheno_mat_t = t(df_pheno_mat) %>% as.data.frame()
  
  df_pheno_pval = df_pheno_all %>% dplyr::select(pheno,cluster,padj ) %>%
    pivot_wider(names_from = "pheno", values_from = "padj") %>%
    as.data.frame()
  rownames(df_pheno_pval) = df_pheno_pval$cluster
  df_pheno_pval$cluster = NULL
  
  df_pheno_pval_t = t(df_pheno_pval) %>% as.data.frame()
  
  column_ord = c('M1','M2','M3','M11','M12','M13','M14','M6','M7','M8','M9')
  #column_ord = levels(se_cluster_pheno$cluster)
  
  pheno_data_ord = na.omit(pheno_data[rownames(df_pheno_mat_t),]) %>% arrange(description)
  row_ord = c("age_death","educ","fsex",rownames(pheno_data_ord),selected_factors)
  
  df_pheno_mat_t[row_ord,] -> df_pheno_mat_t
  rownames(df_pheno_mat_t) = c("Age at death", "Years of education", "Sex (female)", 
                               as.character(pheno_data_ord$description), 
                               selected_factors)
  
  df_pheno_pval_t[row_ord,] -> df_pheno_pval_t
  rownames(df_pheno_pval_t) = rownames(df_pheno_mat_t)
  
  categories_ = c("Demographics", "Demographics", "Demographics", 
                  as.character(pheno_data_ord$category), rep("MOFA Factors",length(selected_factors)))
  
  df_pheno_mat_t[,column_ord] -> df_pheno_mat_t
  df_pheno_pval_t[,column_ord] -> df_pheno_pval_t
  
  to_keep = apply(df_pheno_mat_t, 1, var) != 0
  categories_ = categories_[to_keep]
  
  df_pheno_mat_t = df_pheno_mat_t[to_keep, ]
  df_pheno_pval_t = df_pheno_pval_t[to_keep, ]
  
  df_pheno_mat_t_scaled = as.matrix(t(scale(t(df_pheno_mat_t))))
  
  row_annotation = data.frame(description = rownames(df_pheno_mat_t_scaled))
  row_annotation %>% left_join(pheno_data) %>%
    dplyr::select(category, description) %>%
    column_to_rownames("description") -> row_annotation
  catergory_order = levels(row_annotation$category)
  
  row_annotation$category = factor(categories_, 
                                   levels = c("Demographics",catergory_order,"MOFA Factors"))
  
  row_pal = RColorBrewer::brewer.pal(n = min(8,length(unique(row_annotation$category))), name ="Dark2")
  if(length(unique(row_annotation$category)) > 8) 
    row_pal = c(row_pal, RColorBrewer::brewer.pal(n = length(unique(row_annotation$category)), name ="Set3"))
  
  names(row_pal) = sort(unique(row_annotation$category))
  row_pal = c(Demographics="#1F78B4",
              Cognition="#66A61E",
              Pathology="#E6AB02",
              `Vascular Pathology`="#A6761D",
              Motor="#7570B3",
              Disabilities="#1B9E77",
              Depression="#D95F02",
              Genetics="#E7298A",
              `MOFA Factors`="grey50")
  
  row_ha = rowAnnotation(`Category` = row_annotation$category,
                         col = list(`Category` = row_pal),
                         show_legend = T,
                         simple_anno_size = unit(2, "mm"),
                         annotation_label = "  "
  )
  
  row_split = setNames(row_annotation$category, "")
  
  df_pheno_pval_t.signif <- symnum(as.matrix(df_pheno_pval_t), corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.05, 1), 
                                   symbols = c(quote("*")," "))
  
  ht_opt$TITLE_PADDING = unit(c(4,4), "points")
  
  col_ha_prep <- speakeasy_cluster_df_filtered %>%
    select(cluster, cluster_lv1) %>%
    distinct()  %>%
    slice(match(colnames(df_pheno_mat_t_scaled), cluster))
  
  cluster_ha_col <- HeatmapAnnotation(
    cluster_lv1 = factor(col_ha_prep$cluster_lv1),
    col = list(cluster_lv1 = c('1' = "#F39B7FFF"  , '2' = "#99991EFF", '3' = 'brown')),
    annotation_legend_param = list(cluster_lv1 = list(title = "Cluster")),
    show_annotation_name = FALSE
  )
  
  colnames(df_pheno_mat_t_scaled) <-  c('M1    \n(n=180)', 'M2    \n(n=170)', 'M3    \n(n=153)', 'M11   \n(n=123)', 'M12   \n(n=120)', 
                                        'M13  \n(n=96)', 'M14  \n(n=44)', 'M6    \n(n=156)', 'M7    \n(n=144)', 'M8   \n(n=79)', 'M9   \n(n=74)')
  
  htmap = Heatmap(
    df_pheno_mat_t_scaled, 
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text( df_pheno_pval_t.signif[i,j], x, y, gp = gpar(fontsize = 12))
    },
    cluster_rows = F, cluster_columns = F,
    col = colorRamp2(range(df_pheno_mat_t_scaled), hcl_palette = "RdBu", reverse = TRUE),
    row_split = row_split, row_title = " ",
    row_names_side = "right", show_row_names = T,
    right_annotation = row_ha,
    bottom_annotation = cluster_ha_col,
    column_names_rot = 55,
    heatmap_legend_param = list(
      title = "Scaled aggr. values",
      legend_width = unit(4, "cm"), 
      direction = "vertical"),
    rect_gp = gpar(col = "black", lwd = 1))
  
  assign('clust_heatmap',htmap, envir = .GlobalEnv)
  draw(htmap, merge_legend = TRUE)
  
  ################################################################################################
  return(list(
    df_pheno_all = df_pheno_all,
    df_pheno_mat_t = df_pheno_mat_t,
    df_pheno_mat_t_scaled = df_pheno_mat_t_scaled,
    df_pheno_pval_t.signif = df_pheno_pval_t.signif,
    df_mean_factor_perm = df_mean_factor_perm,
    df_pheno_mean_perm = df_pheno_mean_perm,
    df_pheno_median_perm = df_pheno_median_perm,
    df_pheno_count_perm = df_pheno_count_perm))
}

################################################################################
# Utils
source("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Figuras_MOFA/Scripts/covars/util_covar.R")

# Participants phenotypes data
phenotypes = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/basic_Apr2022_selected_list_Jul2024.rds")

# A data.frame containing the covariates description
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/pheno_list_Jul2024.RData")

# Renaming some covariates
pheno_data <- pheno_data %>%
  mutate(description = gsub('Cognitive decline','Cognitive decline (slope)',description)) %>%
  mutate(description = gsub('Cognitive resilience', 'Cognitive resilience (slope)',description)) %>%
  mutate(description = gsub('AD NIA-Reagan score','NIA-Reagan score', description))

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# Check if MOFA and participants phenotypes samples are in the same order
phenotype_dt = phenotypes[match(mofa@samples_metadata$sample, rownames(phenotypes)), ]
identical(rownames(phenotype_dt), mofa@samples_metadata$sample)
samples_metadata(mofa) <- samples_metadata(mofa) %>% left_join(phenotype_dt %>% rownames_to_column("projid"),by = c("sample"="projid"))
rownames(samples_metadata(mofa)) = samples_metadata(mofa)$sample

# MOFAs Z matrix
mofa_factors = mofa@expectations$Z$group1 %>% as.data.frame()

################################################################################################

# Speakeasy consensus clustering results
speakeasy1_cluster_df = fread("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/projidBycluster.txt")

# Formating 'speakeasy1_cluster_df'
speakeasy1_cluster_df$projid = sprintf("%08d",speakeasy1_cluster_df$projid)
speakeasy1_cluster_df$cluster = paste0("M",speakeasy1_cluster_df$cluster_lv2)
speakeasy1_cluster_df_filtered = speakeasy1_cluster_df %>% 
  group_by(cluster) %>% mutate(n = n()) %>% ungroup() %>%
  filter(n > 15) %>% distinct()
unique(speakeasy1_cluster_df_filtered$cluster)

# Changing the order of lvl1 clusters
speakeasy1_cluster_df_filtered <- speakeasy1_cluster_df_filtered %>%
  mutate(
    cluster_lv1 = case_when(
      cluster_lv1 == 2 ~ 3,
      cluster_lv1 == 3 ~ 2,
      TRUE ~ cluster_lv1
    )
  )

speakeasy1_results = make_phenotype_plot(speakeasy_cluster_df_filtered = speakeasy1_cluster_df_filtered)

###################################################################################################
library(dplyr)
library(ggstats)
library(ggplot2)
library(rstatix)
library(ggsci)
library(stringr)
library(ggpubr)
library(rlang)

title_size = 14
title_xy = 11
text_xy = 9

# Joining phenotypes and clusterization results 
phenotype_dt$projid <- rownames(phenotype_dt)
for_bloxplots = phenotype_dt %>% inner_join(speakeasy1_cluster_df_filtered, by = "projid")


# Boxplot of global cognition function (cluster lvl 1)
for_bloxplots_cogn <- for_bloxplots
for_bloxplots_cogn$cluster_lv1 <- as.factor(for_bloxplots_cogn$cluster_lv1)

# Performing t-test
stat.testt <- t_test(for_bloxplots_cogn, cogn_global_lv ~ cluster_lv1, p.adjust.method = 'fdr')
stat.testt <- stat.testt[stat.testt$p.adj <0.05,]
stat.testt$cluster_lv1 <- ''

# Setting brackets positions
first_bracket <- max(na.omit(for_bloxplots_cogn$cogn_global_lv))
increments <- seq(first_bracket, by = ( 0.11*first_bracket), length.out = nrow(stat.testt))
stat.testt$y.position <- increments

# Boxplot
cogn_plot <- ggplot(for_bloxplots_cogn, aes(x = cluster_lv1, y = cogn_global_lv, fill = cluster_lv1)) +
  geom_boxplot() +
  labs(x = "clusters", y = NULL, title = 'Global Cogn.', fill = NULL) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = title_size, face = 'plain'),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.grid.major.y= element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(size = title_xy),
    axis.title.y = element_text(size = title_xy),
    axis.text.x = element_text(size = text_xy, color = 'black'),
    axis.text.y = element_blank()
  )+
  scale_fill_manual(values =c("#F39B7FFF", "#99991EFF",'brown','gray')) +
  stat_pvalue_manual(stat.testt, 
                     label = "p.adj.signif", 
                     bracket.size = 0.5, 
                     bracket.shorten	= 0.4,
                     tip.length = 0.01,
                     bracket.nudge.y= 0.2
  )

# Selecting AD-related clusters (clusters lvl 2)
for_bloxplots <- filter(for_bloxplots, cluster %in% c('M8','M7','M6','M9'))
for_bloxplots$cluster <- as.factor(for_bloxplots$cluster)

# Boxplot of neurofibrillary tangles

# Performing t-test
stat.testt <- t_test(for_bloxplots, nft_sqrt ~ cluster,p.adjust.method = 'fdr')
stat.testt <- stat.testt[stat.testt$p.adj <0.05,]
stat.testt$cluster <- ''

# Setting brackets positions
first_bracket <- max(na.omit(for_bloxplots$nft_sqrt))
increments <- seq(first_bracket, by = ( 0.11*first_bracket), length.out = nrow(stat.testt))
stat.testt$y.position <- increments

# Boxplot
nft_plot <- ggplot(for_bloxplots, aes(x = cluster, y = nft_sqrt, fill = cluster)) +
  geom_boxplot() +
  labs(x = "clusters", y = NULL, title = 'Neurofibrilary Tangles', fill = NULL) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = title_size, face = 'plain'),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.grid.major.y= element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(size = title_xy),
    axis.title.y = element_text(size = title_xy),
    axis.text.x = element_text(size = text_xy, color = 'black'),
    axis.text.y = element_blank()
  )+
  scale_fill_manual(values = pal_ucscgb(alpha = 1)(length(unique(for_bloxplots$nft_sqrt)))) +
  stat_pvalue_manual(stat.testt, 
                     label = "p.adj.signif", 
                     bracket.size = 0.5, 
                     bracket.shorten	= 0.4,
                     tip.length = 0.01,
                     bracket.nudge.y= 0.2
  )

# Boxplot of global AD pathology

# Performing t-test
stat.testt <- t_test(for_bloxplots, gpath ~ cluster,p.adjust.method = 'fdr')
stat.testt <- stat.testt[stat.testt$p.adj <0.05,]
stat.testt$cluster <- ''

# Setting brackets positions
first_bracket <- max(na.omit(for_bloxplots$gpath))
increments <- seq(first_bracket, by = ( 0.1*first_bracket), length.out = nrow(stat.testt))
stat.testt$y.position <- increments

# Boxplot
gpath_plot <- ggplot(for_bloxplots, aes(x = cluster, y = gpath, fill = cluster)) +
  geom_boxplot() +
  labs(x = "clusters", y = NULL, title = 'Global AD pathology', fill = NULL) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = title_size, face = 'plain'),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.grid.major.y= element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(size = title_xy),
    axis.title.y = element_text(size = title_xy),
    axis.text.x = element_text(size = text_xy, color = 'black'),
    axis.text.y = element_blank()
  )+
  scale_fill_manual(values = pal_ucscgb(alpha = 1)(length(unique(for_bloxplots$gpath)))) +
  stat_pvalue_manual(stat.testt, 
                     label = "p.adj.signif", 
                     bracket.size = 0.5, 
                     bracket.shorten	= 0.4,
                     tip.length = 0.01,
                     bracket.nudge.y= 0
  )

# Boxplot of motor dexterity

# Performing t-test
stat.testt <- t_test(for_bloxplots, motor_dexterity_lv ~ cluster,p.adjust.method = 'fdr')
stat.testt <- stat.testt[stat.testt$p.adj <0.05,]
stat.testt$cluster <- ''

# Setting brackets positions
first_bracket <- max(na.omit(for_bloxplots$motor_dexterity_lv))
increments <- seq(first_bracket, by = ( 0.1*first_bracket), length.out = nrow(stat.testt))
stat.testt$y.position <- increments

# Boxplot
motor_dex_plot <- ggplot(for_bloxplots, aes(x = cluster, y = motor_dexterity_lv, fill = cluster)) +
  geom_boxplot() +
  labs(x = "clusters", y = NULL, title = 'Motor dexterity', fill = NULL) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = title_size, face = 'plain'),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.grid.major.y= element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(size = title_xy),
    axis.title.y = element_text(size = title_xy),
    axis.text.x = element_text(size = text_xy, color = 'black'),
    axis.text.y = element_blank()
  )+
  scale_fill_manual(values = pal_ucscgb(alpha = 1)(length(unique(for_bloxplots$motor_dexterity_lv)))) +
  stat_pvalue_manual(stat.testt, 
                     label = "p.adj.signif", 
                     bracket.size = 0.5, 
                     bracket.shorten	= 0.4,
                     tip.length = 0.01,
                     bracket.nudge.y= 0
  )

# Displaying all plots togheter
library(patchwork)
library(grid)
library(ggplotify)
heatmap_grob <- grid.grabExpr(draw(clust_heatmap, merge_legend = TRUE,annotation_legend_list = list(), 
                                   legend_gap = unit(60, "mm")))
heatmap_grob <- as.ggplot(heatmap_grob)
heatmap_grob <- heatmap_grob + theme(aspect.ratio = 1.2)
boxplots <- cogn_plot/gpath_plot/nft_plot/motor_dex_plot
(heatmap_grob | boxplots) + plot_layout(widths = c(1,0.2))
