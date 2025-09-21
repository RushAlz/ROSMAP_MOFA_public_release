###################################
#  List of metabolites for GSEA   #
###################################
library(readxl)
library(dplyr)

# Retrieving the correctly formated metabolites names
names_ids <- read_excel("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/Metabolites/s1_nomes_IDs.xlsx")
colnames(names_ids) <- c('name', 'super_pathway', 'sub_pathway', 'comp_id', 'platform', 'chemical_id', 'RI', 'mass', 'pubchem', 'cas', 'KEGG', 'HMDB', 'comp_IDstr')
names_ids <- names_ids[-1,]

# Unformating names (so they match MOFAs metabolites names)
names_ids$name <- gsub("[^A-Za-z0-9]", '.',names_ids$name)
names_ids$name <- gsub("^X", '',names_ids$name)

# Creating a list of metabolites and their respective sub pathways
new_metabol_list <- list()
for(i in unique(names_ids$sub_pathway)){
  
  x <- filter(names_ids, sub_pathway == i)
  new_metabol_list[i] <- x
}

# Excluding empty entry
new_metabol_list <- new_metabol_list[-102]

# Adding unknown metabolites
new_metabol_list$unknown <- names_ids$name[602:667]


# Retrieving metabolites related to AD traits
s5 = "C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/Metabolites/s5_all_results.xlsx"

# Storing relevant metabolites on a vector
relevant_metabol <- c()
for(i in excel_sheets(s5)[-1]){
  
  var_s5 <- read_excel(s5, sheet = i)
  x <- filter(var_s5, adj_p<0.05)
  relevant_metabol <- c(relevant_metabol, x$name)
}

# Unformating names (so they match MOFAs metabolites names) 
relevant_metabol <- unique(relevant_metabol)
relevant_metabol <- gsub("[^A-Za-z0-9]", '.', relevant_metabol)
relevant_metabol <- gsub("^X", "", relevant_metabol)


# Assigning the metabolites from 'relevant_metabol' into a list with their respective AD traits
relevant_trait <- list()
for(i in excel_sheets(s5)[-1]){
  
  var_s5 <- read_excel(s5, sheet = i)
  x <- filter(var_s5, adj_p<0.05)
  relevant_trait[i] <- x
}

# Unformating names (so they match MOFAs metabolites names)
relevant_trait <- lapply(relevant_trait, function(x) gsub("[^A-Za-z0-9]", '.', x))
relevant_trait <- lapply(relevant_trait, function(x) gsub("^X", '.', x))

#Lista com os metabólitos comprovados

# Metabolites found to be associated with AD in ROSMAP cohort and in BLSA cohort
BLSA_names <- c('gamma-aminobutyrate (GABA)', 'choline', 'N-acetylglutamate', 'S-adenosylhomocysteine (SAH)', 'S-adenosylmethionine (SAM)')

# Metabolites found to be associated with AD in ROSMAP cohort and in Mayo cohort
mayo = read_excel("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/Metabolites/s8_mayo_results.xlsx")
mayo_names<- mayo$`Supplementary Table 11: Metabolic associations with AD in Mayo cohort.`[2:31]

# Merging 'mayo_names' and 'BLSA_names', and unformating names (so they match MOFAs metabolites names)
proved_AD <- c(mayo_names, BLSA_names)
proved_AD <- gsub("[^A-Za-z0-9]", '.', proved_AD)
proved_AD <- gsub("^X", "", proved_AD)


# Metabolites related to the bioenergetic pathway
carboxyethyl_conjugates <- c('1-carboxyethylisoleucine', '1-carboxyethylleucine', '1-carboxyethylphenylalanine', '1-carboxyethyltyrosine', 
                             '1-carboxyethylvaline')
others<- c('glucose', 'glycerate', 'glucose 6-phosphate', "1,5-anhydroglucitol (1,5-AG)", 'valine', 'isoleucine',
           'leucine', "beta-hydroxyisovalerate", '3-hydroxyisobutyrate', 'isobutyrylcarnitine (C4)',
           'tiglyl carnitine (C5)', '2-methylbutyrylcarnitine (C5)', 'glutarylcarnitine (C5-DC)',
           '5-dodecenoylcarnitine (C12:1)')

# Unformating names (so they match MOFAs metabolites names)
bioenergetic_pathway <- c(others, carboxyethyl_conjugates)
bioenergetic_pathway <- gsub("[^A-Za-z0-9]", '.', bioenergetic_pathway)
bioenergetic_pathway <- gsub("^X", "", bioenergetic_pathway)

# Metabolites related to cholesterol metabolism and sterol pathway
sterol_pathway = c("7-HOCA", '4-cholesten-3-one', 
                   '7-hydroxycholesterol (alpha or beta)')

# Unformating names (so they match MOFAs metabolites names)
sterol_pathway <- gsub("[^A-Za-z0-9]", '.', sterol_pathway)
sterol_pathway <- gsub("^X", "", sterol_pathway)

# Metabolites related to neuroinflammation 

neuroinflammation_pathway <- c('4-hydroxy-nonenal-glutathione', 'cysteinylglycine disulfide*', 
                               'ophthalmate', '15-KETE', '12-HHTrE', "eicosapentaenoate (EPA; 20:5n3)",
                               "docosahexaenoate (DHA; 22:6n3)")

# Unformating names (so they match MOFAs metabolites names)
neuroinflammation_pathway <- gsub("[^A-Za-z0-9]", '.', neuroinflammation_pathway)
neuroinflammation_pathway <- gsub("^X", "", neuroinflammation_pathway)

# Metabolites related to osmoregulation

osmoregulation_pathway <- c('2-aminoadipate', 'arginine', 'glycerophosphorylcholine (GPC)', 
                            'myo-inositol', 'serine', 'urea','betaine')

# Unformating names (so they match MOFAs metabolites names)
osmoregulation_pathway <- gsub("[^A-Za-z0-9]", '.', osmoregulation_pathway)
osmoregulation_pathway <- gsub("^X", "", osmoregulation_pathway)

# Combining all elements into one list
new_metabol_list$relevant_metabol <- relevant_metabol
new_metabol_list <- c(relevant_trait, new_metabol_list)
new_metabol_list$proved_AD <- proved_AD

new_metabol_list$bioenergetic_pathway <- bioenergetic_pathway
new_metabol_list$sterol_pathway <- sterol_pathway
new_metabol_list$neuroinflammation_pathway <- neuroinflammation_pathway
new_metabol_list$osmoregulation_pathway <- osmoregulation_pathway

