library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(purrr)
library(tibble)

args=commandArgs(trailingOnly=TRUE)
wd1 <- args[1]
savefilepath <- args[2]
cancertype <- args[3]
wd2 <- args[4]
numSamples <- args[5]

###### Merge all count files############################################################
setwd(wd1)

filelist <- list.files(path = ".", pattern="*.tsv", full.names=TRUE) 
all_files <- lapply(filelist, function(x) {
    file <- read_tsv(x)
    # Get the start of filename 
    prefix = sub("_.*", "", x)
    small_file <- file %>% select(CDR3, V, J, CDR3_count)
    colnames(small_file) <- c('CDR3', 'V', 'J', paste0(prefix, '_counts'))
    return(small_file)
})

merged_counts <- all_files %>% Reduce(function(df1, df2) full_join(df1, df2, by=c('CDR3', 'V','J')), .)

########3 Merge all freq_min2reads_File#################################################
setwd(wd2)

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


######### Get binary presence/absence dataframe for merged counts########################

get.binary.df <- function(dat, numSamples){
  numSamples <- as.integer(numSamples)
  dat[is.na(dat)] <- 0
  dat <- dat %>% mutate(Clonotype = paste0(CDR3, '_', V, '_', J))
  mat <- dat %>% select(contains('counts')) 
  mat[mat < 2] <- 0
  mat[mat > 1 ] <- 1
  rownames(mat) <- dat$Clonotype
  df <- rownames_to_column(mat) %>%
    mutate(Num_Samples_Present = rowSums(.[2:(numSamples+1)])) %>%
    mutate(Freq_Samples_Present = Num_Samples_Present / numSamples * 100)
  return(df)
}

binary.df <- get.binary.df(merged_counts, numSamples)
write_tsv(binary.df, paste0(savefilepath,cancertype, '_binaryDF.tsv'))

############### get counts of shared / private #########################################

get.count <- function(dat){
  tmp <- dat %>% group_by(Num_Samples_Present) %>% summarise(CDR3_Count=n())
  return(tmp)
} #get counts of shared and private CDR3s

counts_histo <- get.count(binary.df)
counts_distribution <- counts_histo %>% filter ( Num_Samples_Present > 0) %>%
  mutate(log2_CDR3_Count = log2(CDR3_Count)) %>%
  mutate(log3_CDR3_Count = log2(CDR3_Count + 1)) %>%
  mutate ( Freq_CDR3_Count = CDR3_Count / sum(CDR3_Count) *100)

write_tsv(counts_distribution, paste0(savefilepath,cancertype, '_sharedCounts_histo.tsv'))

##################get list of private, less shared and shared TCRs#############################
bmc_normal <- read_tsv('/usr/TCRS_healthy.tsv')
normal <- bmc_normal %>% mutate(V=str_replace(vGeneName, "TCRBV", "TRBV"),
                      J=str_replace(jGeneName, "TCRBJ", "TRBJ"),
                      V=str_replace_all(V, "0", ""),
                      J=str_replace_all(J, "0", ""),
                      Clonotype = paste0(aminoAcid, '_', V, '_', J)) %>%
                      select(Clonotype, aminoAcid, V, J, tcr, n_cmv_public)

####annotate shared/private list with normal healthy
shared30_removeHealthy <- binary.df %>% filter(Freq_Samples_Present >= 30) %>% select(rowname) %>% left_join(normal, by=c('rowname' = 'Clonotype')) %>% filter(is.na(n_cmv_public)) %>% select(rowname)

private_removeHealthy <- binary.df %>%  filter(Num_Samples_Present == 1 ) %>% select(rowname) %>% left_join(normal, by=c('rowname' = 'Clonotype')) %>% filter(is.na(n_cmv_public)) %>%
select(rowname)

lessShared_removeHealthy <- binary.df %>%  filter(Num_Samples_Present > 1 & Freq_Samples_Present < 30 ) %>% select(rowname) %>% left_join(normal, by=c('rowname' = 'Clonotype')) %>% filter(is.na(n_cmv_public)) %>% select(rowname)

write_tsv(shared30_removeHealthy, paste0(savefilepath,cancertype, '_shared30_ClonotypeList_removeHealthy.tsv'))
write_tsv(private_removeHealthy, paste0(savefilepath,cancertype, '_private_ClonotypeList_removeHealthy.tsv'))
write_tsv(lessShared_removeHealthy, paste0(savefilepath,cancertype, '_lessShared_ClonotypeList_removeHealthy.tsv'))

############get length of shared, less shared and private clonotypes#####################

CDR3_df <- binary.df %>%
                separate(rowname, c('CDR3','V','J'), sep='_')%>%
                mutate(CDR3_length=str_length(CDR3)) %>%
                mutate(Clonotype = paste0(CDR3,'_',V, '_', J)) %>%
    left_join(normal, by = 'Clonotype')

removeHealthy.binarydf <- CDR3_df %>% filter( is.na(n_cmv_public))

write_tsv(removeHealthy.binarydf, paste0(savefilepath,cancertype, 'removeHealthy_binaryDF.tsv'))

private_CDR3_length <- removeHealthy.binarydf %>% filter(Num_Samples_Present == 1) %>%
  group_by(CDR3_length) %>%
  summarise(Private_CDR3Length_Count=n())

shared30_CDR3_length <- removeHealthy.binarydf %>%
  filter(Freq_Samples_Present >= 30) %>%
  group_by(CDR3_length) %>%
  summarise(Shared30_CDR3Length_Count=n())

less_shared_length<- removeHealthy.binarydf %>%
  filter(Num_Samples_Present > 1 & Freq_Samples_Present < 30) %>%
  group_by(CDR3_length) %>%
  summarise(lessShared_CDR3Length_Count=n())

#merge private and shared CDR3 together
merged <- full_join(private_CDR3_length, shared30_CDR3_length, by='CDR3_length')
merged_length <- full_join(merged, less_shared_length, by='CDR3_length')
merged_length[is.na(merged_length)] <- 0
merged_length <- merged_length %>% mutate(Private_Freq = Private_CDR3Length_Count/sum(Private_CDR3Length_Count) *100,
                                          Shared30_Freq= Shared30_CDR3Length_Count/sum(Shared30_CDR3Length_Count)*100,
                                          LessShared_Freq = lessShared_CDR3Length_Count/sum(lessShared_CDR3Length_Count)*100)

write_tsv(merged_length, paste0(savefilepath,cancertype, '_CDR3length.tsv'))

##KS test to determine difference between lengths#

kstest1 <- ks.test(unlist(merged_length$Private_CDR3Length_Count), unlist(merged_length$Shared30_CDR3Length_Count))
kstest2 <- ks.test(unlist(merged_length$Private_CDR3Length_Count), unlist(merged_length$lessShared_CDR3Length_Count))
kstest3 <- ks.test(unlist(merged_length$lessShared_CDR3Length_Count), unlist(merged_length$Shared30_CDR3Length_Count))
ks.results1 <- as.data.frame(flatten(kstest1))
ks.results2 <- as.data.frame(flatten(kstest2))
ks.results3 <- as.data.frame(flatten(kstest3))
ks.results <- rbind(ks.results1, ks.results2, ks.results3)

write_tsv(ks.results, paste0(savefilepath,cancertype, '_CDR3length_KStest.tsv'))

private_CDR3_VJusage <- removeHealthy.binarydf %>%
                filter(Num_Samples_Present == 1) %>%
                group_by(V.x,J.x) %>%
                summarise (VJusages_privateCDR3=n()) 
sum <- sum(private_CDR3_VJusage$VJusages_privateCDR3)
private_CDR3_VJusage <- private_CDR3_VJusage %>%
                mutate(Freq_VJgene_privateCDR3= VJusages_privateCDR3/sum*100)

####analyse VJ usage#######
shared_CDR3_VJusage <- removeHealthy.binarydf %>%
                filter(Freq_Samples_Present >= 30) %>%
                group_by(V.x,J.x) %>%
                summarise (VJusages_sharedCDR3=n()) 
sum <- sum(shared_CDR3_VJusage$VJusages_sharedCDR3)
shared_CDR3_VJusage <- shared_CDR3_VJusage %>%
                mutate(Freq_VJgene_sharedCDR3= VJusages_sharedCDR3/sum*100)

VJusage <- full_join(private_CDR3_VJusage, shared_CDR3_VJusage, by=c('V.x','J.x')) 
VJusage[is.na(VJusage)] <- 0
VJusage <- VJusage %>%
    mutate(Fold_Change_SharedFreq_divide_PrivateFreq = Freq_VJgene_sharedCDR3/Freq_VJgene_privateCDR3)

write_tsv(VJusage, paste0(savefilepath,cancertype, '_VJusage.tsv'))

############get list of clonally expanded Clonotypes ####################################

clonaldf <-  clonal %>% mutate(Clonotype = paste0(CDR3, '_', V, '_', J)) %>%
    left_join(normal, by = 'Clonotype')
clonalList_removeHealthy <- clonaldf %>% filter(is.na(n_cmv_public)) %>% select(Clonotype)
 
write_tsv(clonalList_removeHealthy, paste0(savefilepath,cancertype, '_clonallyExpanded_removeHealthy.tsv'))

############## directly matched the shared clonotype sequences to public DB################
vdjdb <- read_tsv('/usr/vdjdb_full.txt')

### extract out HomoSapiens related sequences and rows with cdr3b sequences
vdjdb_HS <- vdjdb %>% filter(species == 'HomoSapiens') %>% filter(! is.na(cdr3.beta))

shared_withVDJDB <- shared30_removeHealthy %>% separate(rowname, c('CDR3','V','J'), sep='_') %>% left_join(vdjdb_HS, by=c('CDR3' = 'cdr3.beta'))

write_tsv(shared_withVDJDB, paste0(savefilepath,cancertype, '_directClonotypeMatch_VDJDB_removeHealthy.tsv'))

############ directly matched the shared clonotype sequences to McPAS###################
McPAS <- read_csv('/usr/McPAS-TCR.csv')

McPAS_HS <- McPAS %>% filter(Species == 'Human') %>% filter(! is.na(CDR3.beta.aa))

shared_withMcPAS <- shared30_removeHealthy %>% separate(rowname, c('CDR3','V','J'), sep='_') %>% left_join(McPAS_HS, by=c('CDR3' = 'CDR3.beta.aa')) 

write_tsv(shared_withMcPAS, paste0(savefilepath,cancertype, '_directClonotypeMatch_McPAS_removeHealthy.tsv'))
