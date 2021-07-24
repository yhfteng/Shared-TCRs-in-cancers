library(readr)
library(dplyr)
library(tidyr)
library(reshape2)
library(tibble)
library(stringr)

args=commandArgs(trailingOnly=TRUE)
savefilepath <- args[1]
wd1 <- args[2] # filepath containing both HNN and NPC freq_min2reads 

########3 Merge all freq_min2reads_File#################################################
setwd(wd1)

filelist <- list.files(path = ".", pattern="*.tsv", full.names=TRUE)
all_files <- lapply(filelist, function(x) {
    file <- read_tsv(x)
    # Get the start of filename 
    prefix = sub("_.*", "", x)
    small_file <- file %>% select(CDR3, V, J, CDR3_count, CDR3_freq_min2reads)
    colnames(small_file) <- c('CDR3', 'V', 'J', paste0(prefix, '_count'), paste0(prefix, '_freq'))
    return(small_file)
})

merged_freq_min2reads <- all_files %>% Reduce(function(df1, df2) full_join(df1, df2, by=c('CDR3', 'V','J')), .)

###### Merge all binary.df files from both cancers#######################################
HNN <- read_tsv('/usr/HNN/HNN_binaryDF.tsv')

NPC <- read_tsv('/usr/NPC/NPC_binaryDF.tsv')

merged.binary.df <- full_join(HNN,NPC, by='rowname')
merged.binary.df[is.na(merged.binary.df)] <- 0
write.table(merged.binary.df, paste0(savefilepath, '2cancers_merged.binary.df.tsv'), sep='\t', row.names=F)

########calculate Jaccard index ###################
dat_mat <- merged.binary.df %>% select( contains('_counts'))
rownames(dat_mat) <- merged.binary.df$rowname

#obtain a pairwise list of different samples
pairwise.list <- t(combn(colnames(dat_mat) , 2))
pairwise.list <- as.data.frame(pairwise.list)

#create function for jaccard's similarity index with own formula, based on binary (present/absent) data
jaccard.method1 <- function(col1, col2) {
  sums <- rowSums(dat_mat[,c(col1, col2)])
  similarity <- length(sums[sums==2])
  total <- length(sums[sums==1]) + similarity
  jaccard <- similarity/total
  return(jaccard)
}

#reiterate jaccard.function for each pair in pairwise.list
jaccard.vector <- apply(pairwise.list, 1, function(x) jaccard.method1(x[1], x[2]))
jaccard.index.df <- cbind(pairwise.list, jaccard.vector)
jaccard.index.df <- jaccard.index.df %>% mutate(V1 = str_sub(V1, 3,8),
						V2 = str_sub(V2, 3,8))

write.table(jaccard.index.df, paste0(savefilepath,'2cancers_Jaccard_index_df.tsv'), sep='\t', row.names=F)

jaccard.table <- spread(jaccard.index.df, V2, jaccard.vector)
write.table(jaccard.table, paste0(savefilepath, '2cancers_Jaccard_table.tsv'), sep='\t', row.names=F)

############## check if clonotypes match to healthy donors and or 1 or 2 cancers#########
bmc_normal <- read_tsv('/usr/TCRs_Healthy.tsv')
normal <- bmc_normal %>% mutate(V=str_replace(vGeneName, "TCRBV", "TRBV"),
                      J=str_replace(jGeneName, "TCRBJ", "TRBJ"),
                      V=str_replace_all(V, "0", ""),
                      J=str_replace_all(J, "0", ""),
                      Clonotype = paste0(aminoAcid, '_', V, '_', J)) %>% 
                      select(Clonotype, aminoAcid, V, J, tcr, n_cmv_public)

shared_clonotype_HNN <- read_tsv('/usr/HNN/HNN_shared30_ClonotypeList.tsv')
shared_clonotype_NPC <- read_tsv('/usr/NPC/NPC_shared30_ClonotypeList.tsv')

anno <- merged.binary.df %>% 
select(rowname, contains('Num_Samples_Present'), contains('Freq_Samples_Present'))
colnames(anno) <- c('Clonotype', 'Num_Samples_Present_HNN', 'Num_Samples_Present_NPC', 'Freq_Samples_Present_HNN', 'Freq_Samples_Present_NPC')

shared_tidyFreq <- function(shared, CancerType){
  s <- unlist(shared)
  dat <-  merged_freq_min2reads %>%
          mutate(Clonotype = paste0(CDR3, '_',V, '_', J)) %>%
          select(Clonotype, contains('_freq')) %>%
          filter(Clonotype %in% s) %>%
	  gather(Sample,Freq, 2:30) %>%
          filter(!is.na(Freq)) %>% 
	  mutate(Patient=str_sub(Sample,3,8))%>%
	  mutate(Cancer = str_sub(Patient,1,3))%>%
	  mutate(Cancer = case_when(Cancer == 'HNN' ~ 'HNN', Cancer == 'NCC' ~ 'NPC')) %>%
          filter(Cancer == CancerType) %>%
          left_join(anno, by='Clonotype')
  return(dat)
}

HNN_shared_tidy <- shared_tidyFreq(shared_clonotype_HNN, 'HNN')
NPC_shared_tidy <- shared_tidyFreq(shared_clonotype_NPC, 'NPC')

dat <- rbind(HNN_shared_tidy, NPC_shared_tidy)
dat <- left_join(dat,normal,by='Clonotype')

dat <- dat %>%
  mutate(Color_breaks = case_when(Cancer == 'HNN' & Num_Samples_Present_HNN == 1 | Num_Samples_Present_HNN == 2  ~ '1-2',
                                  Cancer == 'HNN' & Num_Samples_Present_HNN == 3 | Num_Samples_Present_HNN == 4  ~ '3-4',
                                  Cancer == 'HNN' & Num_Samples_Present_HNN == 5 | Num_Samples_Present_HNN == 6 ~ '5-6',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC == 1 | Num_Samples_Present_NPC == 2  ~ '1-2',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC  == 3 | Num_Samples_Present_NPC == 4  ~ '3-4',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC  == 5 | Num_Samples_Present_NPC == 6  ~ '5-6',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC  == 7 | Num_Samples_Present_NPC == 8  ~ '7-8',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC == 9 | Num_Samples_Present_NPC == 10 ~ '9-10',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC >= 11 ~ '11-12')) %>%
  mutate(Type= ifelse(!is.na(n_cmv_public), 'Healthy',
                                ifelse(Freq_Samples_Present_HNN >= 30 & Freq_Samples_Present_NPC >= 30, 'Both', 'Single'))) %>%
  mutate(Sample = str_replace(Sample, 'NCC', 'NPC'))

write.table(dat, '/usr/sharedClonotypes_Freq_min2reads_withAnnotations.tsv', sep='\t', row.names=F)

######### prepare clonally expanded freq file#############

clonalExpanded_HNN <- read_tsv('/usr/HNN/HNN_clonallyExpanded_df.tsv')
clonalExpanded_NPC <- read_tsv('/usr/NPC/NPC_clonallyExpanded_df.tsv')

CE_tidyFreq <- function(df){
  dat <-  df %>%
          mutate(Clonotype = paste0(CDR3, '_',V, '_', J)) %>%
          mutate(Patient=str_sub(Sample,3,8))%>%
          mutate(Cancer = str_sub(Patient,1,3))%>%
          mutate(Cancer = case_when(Cancer == 'HNN' ~ 'HNN', Cancer == 'NCC' ~ 'NPC')) %>%
          left_join(anno, by='Clonotype')
  return(dat)
}

HNN_CE <- CE_tidyFreq(clonalExpanded_HNN)
NPC_CE <- CE_tidyFreq(clonalExpanded_NPC)

dat_CE <- rbind(HNN_CE,NPC_CE)
dat_CE <- left_join(dat_CE,normal,by='Clonotype')

dat_CE <- dat_CE %>%
  mutate(Color_breaks = case_when(Cancer == 'HNN' & Num_Samples_Present_HNN == 1 | Num_Samples_Present_HNN == 2  ~ '1-2',
                                  Cancer == 'HNN' & Num_Samples_Present_HNN == 3 | Num_Samples_Present_HNN == 4  ~ '3-4',
                                  Cancer == 'HNN' & Num_Samples_Present_HNN == 5 | Num_Samples_Present_HNN == 6 ~ '5-6',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC == 1 | Num_Samples_Present_NPC == 2  ~ '1-2',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC  == 3 | Num_Samples_Present_NPC == 4  ~ '3-4',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC  == 5 | Num_Samples_Present_NPC == 6  ~ '5-6',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC  == 7 | Num_Samples_Present_NPC == 8  ~ '7-8',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC == 9 | Num_Samples_Present_NPC == 10 ~ '9-10',
                                  Cancer == 'NPC' & Num_Samples_Present_NPC >= 11 ~ '11-12')) %>%
  mutate(Type= ifelse(!is.na(n_cmv_public), 'Healthy',
                                ifelse(Freq_Samples_Present_HNN >= 30 & Freq_Samples_Present_NPC >= 30, 'Both', 'Single'))) %>%
  mutate(Sample = str_replace(Sample, 'NCC', 'NPC'))

write.table(dat, '/usr/clonallyExpanded_Clonotypes_Freq_min2reads_withAnnotations.tsv', sep='\t', row.names=F)

############get summary table of shared freq##########
dat_summary <- dat %>% group_by(Freq <0.1) %>% summarise(Count= n())

##proportions graph#######
shared_Proportions <- dat %>%
              group_by(Cancer, Clonotype, Type) %>%
              summarise(Num_pat_perClonotype_perPanCancer=n()) %>%
              group_by(Cancer,Type) %>%
              summarise(Num_Group= n())

CE_annotations <- dat_CE %>%
  mutate(Type_HNN = case_when(!is.na(n_cmv_public) ~ 'Healthy',
                              Num_Samples_Present_HNN == 1 ~ 'Private',
                              Num_Samples_Present_HNN > 1 & Freq_Samples_Present_HNN <30 ~ 'Less_Shared',
                              Freq_Samples_Present_HNN >= 30 ~ 'Shared'),
         Type_NPC = case_when(!is.na(n_cmv_public) ~ 'Healthy',
                              Num_Samples_Present_NPC == 1 ~ 'Private',
                              Num_Samples_Present_NPC > 1 & Freq_Samples_Present_NPC <
                                30 ~ 'Less_Shared',
                              Freq_Samples_Present_NPC >= 30 ~ 'Shared'),
         SpecificGroup = case_when(Type_HNN == 'Healthy' | Type_NPC == 'Healthy' ~ 'Healthy',
                                   Type_HNN == 'Shared' & Type_NPC == 'Shared' ~ 'Shared_Both_Cancers',
                                   Type_HNN == 'Shared' | Type_NPC == 'Shared' ~ 'Shared_Single_Cancer',
                                   Type_HNN != 'Shared' & Cancer == 'HNN' ~ Type_HNN,
                                   Type_NPC != 'Shared' & Cancer == 'NPC' ~ Type_NPC))

CE_proportions <- CE_annotations %>%
  group_by(Cancer, Clonotype,SpecificGroup) %>%
  summarise(Num_clonotype_perCancer_perSpecificGroup = n()) %>%
  group_by(Cancer, SpecificGroup) %>%
  summarise(Num_perSpecificGroup = n())

write_tsv(shared_Proportions,'/usr/shared_Proportions.tsv')
write_tsv(CE_annotations, '/usr/clonallyExpanded_annotations.tsv')
write_tsv(CE_proportions , '/usr/clonallyExpanded_Proportions.tsv')


