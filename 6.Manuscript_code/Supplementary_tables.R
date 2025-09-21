#########################
# Supplementary tables #
########################

# GSEA table

library(dplyr)

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

# Filtering significant padj
final_table <- filter(final_table, padj < 0.05)

library(writexl)

# Creating a list of factors with their respective data
xx = split(final_table, final_table$factor)
names(xx) <- c('Factor 1','Factor 2','Factor 3','Factor 4','Factor 8','Factor 12','Factor 14', 'Factor 26','Factor 42')

empty_df <- data.frame()
sheets <- c(list("Blank" = empty_df), xx)
#######################################################################################################################

# Factor loadings table (only loadings > 0.4)

# Factor 1 loadings
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/df_weights_F1")
df_filtered1 <- filter(df_weights_F1, value > 0.4)

# Factor 2 loadings
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/df_weights_F2")
df_filtered2 <- filter(df_weights_F2, value > 0.4)

# Factor 3 loadings
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/df_weights_F3")
df_filtered3 <- filter(df_weights_F3, value > 0.4)

# Factor 4 loadings
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/df_weights_F4")
df_filtered4 <- filter(df_weights_F4, value > 0.4)

# Factor 8 loadings
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/df_weights_F8")
df_filtered8 <- filter(df_weights_F8, value > 0.4)

# Factor 12 loadings
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/df_weights_F12")
df_filtered12 <- filter(df_weights_F12, value > 0.4)

# Factor 14 loadings
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/df_weights_F14")
df_filtered14 <- filter(df_weights_F14, value > 0.4)

# Factor 26 loadings
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/df_weights_F26")
df_filtered26 <- filter(df_weights_F26, value > 0.4)

# Factor 42 loadings
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/df_weights_F42")
df_filtered42 <- filter(df_weights_F42, value > 0.4)

# Grouping the factor loadings on a list
yy = list(df_filtered1,df_filtered2,df_filtered3,df_filtered4,df_filtered8,df_filtered12,df_filtered14,df_filtered26,df_filtered42)
names(yy) <- c('Factor 1','Factor 2','Factor 3','Factor 4','Factor 8','Factor 12','Factor 14', 'Factor 26','Factor 42')


empty_df <- data.frame()
sheets <- c(list("Blank" = empty_df), yy)
#############################################################################################################################

# Participants characteristics table

library(reporter)
library(magrittr)
library(tidyverse)
library(moments)
library(gt)
library(writexl)

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# Meta 1
phenotypes = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/basic_Apr2022_selected_list_Jul2024.rds")
phenotype_dt = phenotypes[match(mofa@samples_metadata$sample, rownames(phenotypes)), ]
identical(rownames(phenotype_dt), mofa@samples_metadata$sample)

# Metadata 2
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/meta.RDATA")
phenotype2 = meta[match(mofa@samples_metadata$sample, meta$projid), ]
rownames(phenotype2) <- phenotype2$projid
identical(rownames(phenotype2), mofa@samples_metadata$sample)


# Removing overlapping columns from metadata 1 and 2
common_cols <- intersect(names(phenotype_dt), names(phenotype2))
phenotype2 <- phenotype2[,!(names(phenotype2) %in% common_cols)]

# Combining the two metadada tables
phenotype3 <- cbind(phenotype_dt,phenotype2)

# Retriving omics and cell type data
rna_AC = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_AC.rds")
rna_DLPF = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_DLPFC.rds")
rna_PCG = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/rnaseq_bulk_PCG.rds")
tmt = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/dlpfc_tmt.rds")
acetil = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/histoneAcetylation_H3K9ac.rds")
metabol = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/brain_metabolomics.rds")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/sn_RNA_norm.RDATA")#sn_RNA_norm

# Subset the merged metadata to create sample-specific metadata tables for each omic dataset.
`RNA (AC)` <- phenotype3[match(colnames(rna_AC), rownames(phenotype3)),]
`RNA (DLPFC)` <- phenotype3[match(colnames(rna_DLPF), rownames(phenotype3)),]
`RNA (PCG)` <- phenotype3[match(colnames(rna_PCG), rownames(phenotype3)),]
Proteomics <- phenotype3[match(colnames(tmt), rownames(phenotype3)),]
H3K9ac <- phenotype3[match(colnames(acetil), rownames(phenotype3)),]
Metabolites <- phenotype3[match(colnames(metabol), rownames(phenotype3)),]
`Cell types` <- phenotype3[match(colnames(sn_RNA_norm), rownames(phenotype3)),]

# Creating data.frame with participants data
participants_summary = phenotype3 %>% 
  summarise(n = n(),
            `Age at death, mean (SD), y` = sprintf("%.1f (%.1f)", 
                                                   mean(age_death), sd(age_death)),
            
            `Male sex, No. (%)` = sprintf("%.d (%.1f)", 
                                          sum(msex==1), 100*sum(msex==1)/length(msex)),
            
            `Educational level, mean (SD), y` = sprintf("%.1f (%.1f)",
                                                        mean(educ, na.rm = TRUE), sd(educ, na.rm = TRUE)),
            
            `Alzheimer's dementia, No. (%)` =  sprintf("%.d (%.1f)", 
                                                       sum(ad_dementia_status==TRUE, na.rm = T), 
                                                       100*sum(ad_dementia_status==TRUE, na.rm = T)/length(ad_dementia_status)),
            
            `Cognitive impairment, No. (%)` =  sprintf("%.d (%.1f)", 
                                                       sum(cog_impairment_status==TRUE, na.rm = T), 
                                                       100*sum(cog_impairment_status==TRUE, na.rm = T)/length(cog_impairment_status)),
            
            `Global cognition function, mean (SD) - at death` = sprintf("%.2f (%.2f)", 
                                                                        mean(cogn_global_lv, na.rm = T),
                                                                        sd(cogn_global_lv, na.rm = T)),
            
            `Cognitive decline, mean (SD) - at death` = sprintf("%.2f (%.2f)", 
                                                                mean(cogng_demog_slope, na.rm = T),
                                                                sd(cogng_demog_slope, na.rm = T)),
            
            `Cognitive resilience, mean (SD) - at death` = sprintf("%.2f (%.2f)", 
                                                                   mean(cogng_path_slope, na.rm = T),
                                                                   sd(cogng_path_slope, na.rm = T)),
            
            `NIA-Reagan AD, No. (%)ᵃ` = sprintf("%d (%.1f)", sum(niareagansc>2, na.rm = T),
                                                100*sum(niareagansc>2, na.rm = T)/sum(!is.na(niareagansc))),
            
            `Global AD pathologic score, median (IQR)` = sprintf("%.2f (%.2f-%.2f)", 
                                                                 summary(gpath, na.rm = T)[3],
                                                                 summary(gpath, na.rm = T)[2],
                                                                 summary(gpath, na.rm = T)[5]),
            
            
            `Beta-amyloid, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                mean(amyloid, na.rm = T),
                                                sd(amyloid, na.rm = T)),
            
            `Tangle density, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                  mean(tangles, na.rm = T),
                                                  sd(tangles, na.rm = T)),
            
            `Diffuse plaques score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                         mean(plaq_d, na.rm = T),
                                                         sd(plaq_d, na.rm = T)),
            
            `Neuritic plaques score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                          mean(plaq_n, na.rm = T),
                                                          sd(plaq_n, na.rm = T)),
            
            `Neurofibrillary tangles score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                                 mean(nft, na.rm = T),
                                                                 sd(nft, na.rm = T)),
            
            `TDP-43, No. (%)ᵇ` = sprintf("%d (%.1f)", 
                                         sum(tdp_43_binary == TRUE, na.rm = T), 
                                         100*sum(tdp_43_binary == TRUE, na.rm = T)/sum(!is.na(tdp_43_binary))),
            
            `Neocortical Lewy bodies, No. (%)*` = sprintf("%d (%.1f)",
                                                          sum(dlbany==1, na.rm = T), 
                                                          100*sum(dlbany==1, na.rm = T)/sum(!is.na(dlbany))),
            
            `PD pathology, No. (%)*` = sprintf("%d (%.1f)",
                                               sum(PD_pathology_YN==TRUE,na.rm = T),
                                               100*sum(PD_pathology_YN == TRUE, na.rm = T)/sum(!is.na(PD_pathology_YN))),
            
            `Hippocampal sclerosis, No. (%)*` = sprintf("%d (%.1f)",
                                                        sum(hipscl_binary==TRUE, na.rm = T), 
                                                        100*sum(hipscl_binary==TRUE, na.rm = T)/sum(!is.na(hipscl_binary))),
            
            `Arteriolosclerosis, No. (%)ᶜ` = sprintf("%d (%.1f)",
                                                     sum(arteriol_scler>1, na.rm = T), 
                                                     100*sum(arteriol_scler>1, na.rm = T)/sum(!is.na(arteriol_scler))),
            
            `Cerebral atherosclerosis, No. (%)ᶜ` = sprintf("%d (%.1f)",
                                                           sum(cvda_4gp2>1, na.rm = T), 
                                                           100*sum(cvda_4gp2>1, na.rm = T)/sum(!is.na(cvda_4gp2))),
            
            `Cerebral amyloid angiopathy, No. (%)ᶜ` = sprintf("%d (%.1f)", sum(caa_4gp>1, na.rm = T), 
                                                              100*sum(caa_4gp>1, na.rm = T)/sum(!is.na(caa_4gp))),
            
            `Macroscopic infarcts, No. (%)` = sprintf("%d (%.1f)",
                                                      sum(ci_num2_gct==1),
                                                      100*sum(ci_num2_gct==1)/length(ci_num2_gct)),
            
            `Microinfarcts, No. (%)` = sprintf("%d (%.1f)",
                                               sum(ci_num2_mct==1),
                                               100*sum(ci_num2_mct==1)/length(ci_num2_mct)),
            
            `Frailty score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                 mean(frailty_z_lv, na.rm = T),
                                                 sd(frailty_z_lv, na.rm = T)),
            
            `Bradykinesia score , median (IQR)` = sprintf("%.2f (%.2f-%.2f)", 
                                                          summary(bradysc_lv, na.rm = T)[3],
                                                          summary(bradysc_lv, na.rm = T)[2],
                                                          summary(bradysc_lv, na.rm = T)[5]),
            
            `Gait score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                              mean(gaitsc_lv, na.rm = T),
                                              sd(gaitsc_lv, na.rm = T)),
            
            `Rigidity score , median (IQR)` = sprintf("%.2f (%.2f-%.2f)", 
                                                      summary(rigidsc_lv, na.rm = T)[3],
                                                      summary(rigidsc_lv, na.rm = T)[2],
                                                      summary(rigidsc_lv, na.rm = T)[5]),
            
            `Tremor score , median (IQR)` = sprintf("%.2f (%.2f-%.2f)", 
                                                    summary(tremsc_lv, na.rm = T)[3],
                                                    summary(tremsc_lv, na.rm = T)[2],
                                                    summary(tremsc_lv, na.rm = T)[5]),
            
            `Global parkinsonian score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                             mean(sqrt_parksc_demog_slope, na.rm = T),
                                                             sd(sqrt_parksc_demog_slope, na.rm = T)),
            
            `Motor dexterity, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                   mean(motor_dexterity_lv, na.rm = T),
                                                   sd(motor_dexterity_lv, na.rm = T)),
            
            `Motor gait, mean (SD)` = sprintf("%.2f (%.2f)", 
                                              mean(motor_gait_lv, na.rm = T),
                                              sd(motor_gait_lv, na.rm = T)),
            
            `Motor hand strength, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                       mean(motor_handstreng_lv, na.rm = T),
                                                       sd(motor_handstreng_lv, na.rm = T)),
            
            `Motor functions, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                   mean(motor10_demog_slope, na.rm = T),
                                                   sd(motor10_demog_slope, na.rm = T)),
            
            `Basic activities of daily living, mean (SD)` = sprintf("%.2f (%.2f-%.2f)", 
                                                                    summary(katzsum_lv, na.rm = T)[3],
                                                                    summary(katzsum_lv, na.rm = T)[2],
                                                                    summary(katzsum_lv, na.rm = T)[5]),
            
            `Instrumental activities of daily living, mean (SD)` = sprintf("%.2f (%.2f-%.2f)", 
                                                                           summary(iadlsum_lv, na.rm = T)[3],
                                                                           summary(iadlsum_lv, na.rm = T)[2],
                                                                           summary(iadlsum_lv, na.rm = T)[5]),
            
            `Mobility disability, mean (SD)` = sprintf("%.2f (%.2f-%.2f)", 
                                                       summary(rosbsum_lv, na.rm = T)[3],
                                                       summary(rosbsum_lv, na.rm = T)[2],
                                                       summary(rosbsum_lv, na.rm = T)[5]),
            
            `Depressive symptoms (CES-D), median (IQR)` = sprintf("%.2f (%.2f-%.2f)", 
                                                                  summary(cesdsum_lv, na.rm = T)[3],
                                                                  summary(cesdsum_lv, na.rm = T)[2],
                                                                  summary(cesdsum_lv, na.rm = T)[5]),
            
            `Major Depressive Disorder, NO. (%)ᵈ` = sprintf("%d (%.1f)",
                                                            sum(r_depres_status==TRUE, na.rm = T), 
                                                            100*sum(r_depres_status==TRUE, na.rm = T)/sum(!is.na(r_depres_status))),
            
            `APOE - e4, NO. (%)` = sprintf("%d (%.1f)",
                                           sum(apoe_any4==T, na.rm = T), 
                                           100*sum(apoe_any4==T, na.rm = T)/sum(!is.na(apoe_any4))),
            
            `TOMM40'523-L, NO. (%)` = sprintf("%d (%.1f)",
                                              sum(tomm40_long==T, na.rm = T), 
                                              100*sum(tomm40_long==T, na.rm = T)/sum(!is.na(tomm40_long))),
            
  ) %>%
  t() %>% as.data.frame()

colnames(participants_summary) = 'Overall Data'
participants_summary = participants_summary %>% rownames_to_column("Description")

# Creating the same data.frame but for individual omics
for(omics in c('RNA (AC)','RNA (DLPFC)','RNA (PCG)','Proteomics','H3K9ac','Metabolites','Cell types')){
  
  participants_summary_2 = get(omics) %>% 
    summarise(n = n(),
              `Age at death, mean (SD), y` = sprintf("%.1f (%.1f)", 
                                                     mean(age_death), sd(age_death)),
              
              `Male sex, No. (%)` = sprintf("%.d (%.1f)", 
                                            sum(msex==1), 100*sum(msex==1)/length(msex)),
              
              `Educational level, mean (SD), y` = sprintf("%.1f (%.1f)",
                                                          mean(educ, na.rm = TRUE), sd(educ, na.rm = TRUE)),
              
              `Alzheimer's dementia, No. (%)` =  sprintf("%.d (%.1f)", 
                                                         sum(ad_dementia_status==TRUE, na.rm = T), 
                                                         100*sum(ad_dementia_status==TRUE, na.rm = T)/length(ad_dementia_status)),
              
              `Cognitive impairment, No. (%)` =  sprintf("%.d (%.1f)", 
                                                         sum(cog_impairment_status==TRUE, na.rm = T), 
                                                         100*sum(cog_impairment_status==TRUE, na.rm = T)/length(cog_impairment_status)),
              
              `Global cognition function, mean (SD) - at death` = sprintf("%.2f (%.2f)", 
                                                                          mean(cogn_global_lv, na.rm = T),
                                                                          sd(cogn_global_lv, na.rm = T)),
              
              `Cognitive decline, mean (SD) - at death` = sprintf("%.2f (%.2f)", 
                                                                  mean(cogng_demog_slope, na.rm = T),
                                                                  sd(cogng_demog_slope, na.rm = T)),
              
              `Cognitive resilience, mean (SD) - at death` = sprintf("%.2f (%.2f)", 
                                                                     mean(cogng_path_slope, na.rm = T),
                                                                     sd(cogng_path_slope, na.rm = T)),
              
              `NIA-Reagan AD, No. (%)ᵃ` = sprintf("%d (%.1f)", sum(niareagansc>2, na.rm = T),
                                                  100*sum(niareagansc>2, na.rm = T)/sum(!is.na(niareagansc))),
              
              `Global AD pathologic score, median (IQR)` = sprintf("%.2f (%.2f-%.2f)", 
                                                                   summary(gpath, na.rm = T)[3],
                                                                   summary(gpath, na.rm = T)[2],
                                                                   summary(gpath, na.rm = T)[5]),
              
              
              `Beta-amyloid, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                  mean(amyloid, na.rm = T),
                                                  sd(amyloid, na.rm = T)),
              
              `Tangle density, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                    mean(tangles, na.rm = T),
                                                    sd(tangles, na.rm = T)),
              
              `Diffuse plaques score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                           mean(plaq_d, na.rm = T),
                                                           sd(plaq_d, na.rm = T)),
              
              `Neuritic plaques score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                            mean(plaq_n, na.rm = T),
                                                            sd(plaq_n, na.rm = T)),
              
              `Neurofibrillary tangles score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                                   mean(nft, na.rm = T),
                                                                   sd(nft, na.rm = T)),
              
              `TDP-43, No. (%)ᵇ` = sprintf("%d (%.1f)", 
                                           sum(tdp_43_binary == TRUE, na.rm = T), 
                                           100*sum(tdp_43_binary == TRUE, na.rm = T)/sum(!is.na(tdp_43_binary))),
              
              `Neocortical Lewy bodies, No. (%)*` = sprintf("%d (%.1f)",
                                                            sum(dlbany==1, na.rm = T), 
                                                            100*sum(dlbany==1, na.rm = T)/sum(!is.na(dlbany))),
              
              `PD pathology, No. (%)*` = sprintf("%d (%.1f)",
                                                 sum(PD_pathology_YN==TRUE,na.rm = T),
                                                 100*sum(PD_pathology_YN == TRUE, na.rm = T)/sum(!is.na(PD_pathology_YN))),
              
              `Hippocampal sclerosis, No. (%)*` = sprintf("%d (%.1f)",
                                                          sum(hipscl_binary==TRUE, na.rm = T), 
                                                          100*sum(hipscl_binary==TRUE, na.rm = T)/sum(!is.na(hipscl_binary))),
              
              `Arteriolosclerosis, No. (%)ᶜ` = sprintf("%d (%.1f)",
                                                       sum(arteriol_scler>1, na.rm = T), 
                                                       100*sum(arteriol_scler>1, na.rm = T)/sum(!is.na(arteriol_scler))),
              
              `Cerebral atherosclerosis, No. (%)ᶜ` = sprintf("%d (%.1f)",
                                                             sum(cvda_4gp2>1, na.rm = T), 
                                                             100*sum(cvda_4gp2>1, na.rm = T)/sum(!is.na(cvda_4gp2))),
              
              `Cerebral amyloid angiopathy, No. (%)ᶜ` = sprintf("%d (%.1f)", sum(caa_4gp>1, na.rm = T), 
                                                                100*sum(caa_4gp>1, na.rm = T)/sum(!is.na(caa_4gp))),
              
              `Macroscopic infarcts, No. (%)` = sprintf("%d (%.1f)",
                                                        sum(ci_num2_gct==1),
                                                        100*sum(ci_num2_gct==1)/length(ci_num2_gct)),
              
              `Microinfarcts, No. (%)` = sprintf("%d (%.1f)",
                                                 sum(ci_num2_mct==1),
                                                 100*sum(ci_num2_mct==1)/length(ci_num2_mct)),
              
              `Frailty score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                   mean(frailty_z_lv, na.rm = T),
                                                   sd(frailty_z_lv, na.rm = T)),
              
              `Bradykinesia score , median (IQR)` = sprintf("%.2f (%.2f-%.2f)", 
                                                            summary(bradysc_lv, na.rm = T)[3],
                                                            summary(bradysc_lv, na.rm = T)[2],
                                                            summary(bradysc_lv, na.rm = T)[5]),
              
              `Gait score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                mean(gaitsc_lv, na.rm = T),
                                                sd(gaitsc_lv, na.rm = T)),
              
              `Rigidity score , median (IQR)` = sprintf("%.2f (%.2f-%.2f)", 
                                                        summary(rigidsc_lv, na.rm = T)[3],
                                                        summary(rigidsc_lv, na.rm = T)[2],
                                                        summary(rigidsc_lv, na.rm = T)[5]),
              
              `Tremor score , median (IQR)` = sprintf("%.2f (%.2f-%.2f)", 
                                                      summary(tremsc_lv, na.rm = T)[3],
                                                      summary(tremsc_lv, na.rm = T)[2],
                                                      summary(tremsc_lv, na.rm = T)[5]),
              
              `Global parkinsonian score, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                               mean(sqrt_parksc_demog_slope, na.rm = T),
                                                               sd(sqrt_parksc_demog_slope, na.rm = T)),
              
              `Motor dexterity, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                     mean(motor_dexterity_lv, na.rm = T),
                                                     sd(motor_dexterity_lv, na.rm = T)),
              
              `Motor gait, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                mean(motor_gait_lv, na.rm = T),
                                                sd(motor_gait_lv, na.rm = T)),
              
              `Motor hand strength, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                         mean(motor_handstreng_lv, na.rm = T),
                                                         sd(motor_handstreng_lv, na.rm = T)),
              
              `Motor functions, mean (SD)` = sprintf("%.2f (%.2f)", 
                                                     mean(motor10_demog_slope, na.rm = T),
                                                     sd(motor10_demog_slope, na.rm = T)),
              
              `Basic activities of daily living, mean (SD)` = sprintf("%.2f (%.2f-%.2f)", 
                                                                      summary(katzsum_lv, na.rm = T)[3],
                                                                      summary(katzsum_lv, na.rm = T)[2],
                                                                      summary(katzsum_lv, na.rm = T)[5]),
              
              `Instrumental activities of daily living, mean (SD)` = sprintf("%.2f (%.2f-%.2f)", 
                                                                             summary(iadlsum_lv, na.rm = T)[3],
                                                                             summary(iadlsum_lv, na.rm = T)[2],
                                                                             summary(iadlsum_lv, na.rm = T)[5]),
              
              `Mobility disability, mean (SD)` = sprintf("%.2f (%.2f-%.2f)", 
                                                         summary(rosbsum_lv, na.rm = T)[3],
                                                         summary(rosbsum_lv, na.rm = T)[2],
                                                         summary(rosbsum_lv, na.rm = T)[5]),
              
              `Depressive symptoms (CES-D), median (IQR)` = sprintf("%.2f (%.2f-%.2f)", 
                                                                    summary(cesdsum_lv, na.rm = T)[3],
                                                                    summary(cesdsum_lv, na.rm = T)[2],
                                                                    summary(cesdsum_lv, na.rm = T)[5]),
              
              `Major Depressive Disorder, NO. (%)ᵈ` = sprintf("%d (%.1f)",
                                                              sum(r_depres_status==TRUE, na.rm = T), 
                                                              100*sum(r_depres_status==TRUE, na.rm = T)/sum(!is.na(r_depres_status))),
              
              `APOE - e4, NO. (%)` = sprintf("%d (%.1f)",
                                             sum(apoe_any4==T, na.rm = T), 
                                             100*sum(apoe_any4==T, na.rm = T)/sum(!is.na(apoe_any4))),
              
              `TOMM40'523-L, NO. (%)` = sprintf("%d (%.1f)",
                                                sum(tomm40_long==T, na.rm = T), 
                                                100*sum(tomm40_long==T, na.rm = T)/sum(!is.na(tomm40_long))),
              
    )%>% 
    t() %>% as.data.frame()
  
  colnames(participants_summary_2) = omics
  participants_summary_2 = participants_summary_2 %>% rownames_to_column("Description")
  
  # Joining data.frames
  participants_summary=left_join(participants_summary,participants_summary_2)
}

# Adding footnotes
participants_summary_gt = participants_summary %>%
  gt() %>% 
  tab_footnote("ᵃ Intermediate or high likelihood.") %>%
  tab_footnote("ᵇ Inclusion beyond the amygdala.") %>%
  tab_footnote("ᶜ Moderate or severe.") %>%
  tab_footnote("ᵈ Highly probable.") %>%
  tab_footnote("* Check variables definitions.")

participants_summary_gt

#gtsave(participants_summary_gt, "C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Figuras_MOFA/pre_inkscape/geral/table1.pdf")


############################################################################################################################################

# A .txt document with all factor loadings

load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/ALL_df_weights_F1")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/ALL_df_weights_F2")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/ALL_df_weights_F3")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/ALL_df_weights_F4")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/ALL_df_weights_F8")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/ALL_df_weights_F12")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/ALL_df_weights_F14")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/ALL_df_weights_F26")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/weights/ALL_df_weights_F42")

yy = list(all_df_weight_F1, all_df_weight_F2, all_df_weight_F3, all_df_weight_F4, all_df_weight_F8, all_df_weight_F12, all_df_weight_F14,
          all_df_weight_F26, all_df_weight_F42)
combined_dfs <- do.call(rbind,yy)
##############################################################################################################################################

library(MOFA2)
library(reshape2)
library(stringr)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyverse)
library(writexl)

# Utils
source("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Figuras_MOFA/Scripts/covars/util_covar.R")

# Retrieving participants phenotypes data
phenotypes = readRDS("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/basic_Apr2022_selected_list_Jul2024.rds")
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/pheno_list_Jul2024.RData")

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# Checking if samples are in the same order
phenotype_dt = phenotypes[match(mofa@samples_metadata$sample, rownames(phenotypes)), ]
identical(rownames(phenotype_dt), mofa@samples_metadata$sample)
samples_metadata(mofa) <- samples_metadata(mofa) %>% left_join(phenotype_dt %>% rownames_to_column("projid"),by = c("sample"="projid"))
rownames(samples_metadata(mofa)) = samples_metadata(mofa)$sample

# MOFAs Z matrix
mofa_factors = mofa@expectations$Z$group1 %>% as.data.frame()

# Selecting only the AD-related factors
selected_factors = paste0("Factor",c(1, 2, 3, 4, 8, 12, 14, 26, 42))
mofa_factors = mofa_factors[,selected_factors]

# Performing linear regressions (adjusted by age of death, sex and years of education)
res_test = run_module_trait_association(data4linear_reg = mofa_factors,
                                        phenotype_dt = phenotype_dt,
                                        pheno_list = pheno_list,
                                        covariates = c("age_death","msex","educ"),
                                        verbose = F)

# Retrieving beta values
res_beta <- res_test$matrix_beta
res_beta <- t(res_beta)

# Adjusting P-values by all phenotypes together.
matrix_pvalue = res_test$matrix_pvalue
adj_matrix_pvalue = matrix(p.adjust(as.vector(as.matrix(matrix_pvalue)), method='fdr'),ncol=ncol(matrix_pvalue))
dimnames(adj_matrix_pvalue) = dimnames(matrix_pvalue)
res_pvalues <- t(adj_matrix_pvalue)


# Reshaping adj p-value matrix to apropriate format
pval_long <- res_pvalues %>%
  as.data.frame() %>%
  rownames_to_column(var = "Factor") %>%
  pivot_longer(-Factor, names_to = "phenotype", values_to = "adj.pvalue")

# Reshaping beta matrix to apropriate format
beta_long <- res_beta %>%
  as.data.frame() %>%
  rownames_to_column(var = "Factor") %>%
  pivot_longer(-Factor, names_to = "phenotype", values_to = "estimate")

# Merge the two long tables by Factor + phenotype
merged_df <- left_join(beta_long, pval_long, by = c("Factor", "phenotype"))


# The same process but for all the 50 factors

mofa_factors = mofa@expectations$Z$group1 %>% as.data.frame()
selected_factors = paste0("Factor",seq(1,50))
mofa_factors = mofa_factors[,selected_factors]

res_test = run_module_trait_association(data4linear_reg = mofa_factors,
                                        phenotype_dt = phenotype_dt,
                                        pheno_list = pheno_list,
                                        covariates = c("age_death","msex","educ"),
                                        verbose = F)

res_beta <- res_test$matrix_beta
res_beta <- t(res_beta)

matrix_pvalue = res_test$matrix_pvalue
adj_matrix_pvalue = matrix(p.adjust(as.vector(as.matrix(matrix_pvalue)), method='fdr'),ncol=ncol(matrix_pvalue))
dimnames(adj_matrix_pvalue) = dimnames(matrix_pvalue)
res_pvalues <- t(adj_matrix_pvalue)

pval_long <- res_pvalues %>%
  as.data.frame() %>%
  rownames_to_column(var = "Factor") %>%
  pivot_longer(-Factor, names_to = "phenotype", values_to = "adj.pvalue")

beta_long <- res_beta %>%
  as.data.frame() %>%
  rownames_to_column(var = "Factor") %>%
  pivot_longer(-Factor, names_to = "phenotype", values_to = "estimate")

merged_df <- left_join(beta_long, pval_long, by = c("Factor", "phenotype"))
