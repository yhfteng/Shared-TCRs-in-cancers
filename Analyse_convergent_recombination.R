library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)

args=commandArgs(trailingOnly=TRUE)
filepath <- args[1]
setwd(filepath)
binary.df <- read_tsv(args[2]) #input remove healthy binary df
shared <- read_tsv(args[3]) # input shared30 remove healthy list
savefilepath <- args[4]

################put all the files of interest into a list###########################
filenames <- list.files(path=filepath, pattern ="_Clonotype_df.tsv")
fileID <- substring(filenames,1,6)
filelist <- lapply(filenames, read_tsv)
names(filelist) <- fileID

#############################extract shared CDR3#########################################
shared_clonotype <- unlist(shared)

####### create 'clonotype column in individual clonotype_df file and count number of nt##
get.ntCount.eachCDR3 <- function(dat){
  selected_df <- dat %>%
		mutate(Clonotype= paste0('C',`CDR3-IMGT.x`,'F', '_', V, '_', J)) %>%  
                  filter(Clonotype %in% shared_clonotype) %>%
                  group_by(Clonotype, `CDR3-IMGT.y`) %>%
                  summarise(Counts_of_nt_per_Clonotype=n())
  return (selected_df)
}

#########execute function on all files in filelist###################################
shared_clonotype_nt <- lapply(filelist, get.ntCount.eachCDR3)
ans <- map_df(shared_clonotype_nt, ~as.data.frame(.x), id='FileNum')

write.table(ans, paste0(savefilepath, '_sharedClonotype_NtList_Counts.tsv'), sep='\t', row.names=F)

num_pat_withsameNt_perCDR3 <- ans %>%
    group_by(Clonotype,`CDR3-IMGT.y`) %>%
    summarise(Num_Pat_sameNt_perClonotype=n())

write.table(num_pat_withsameNt_perCDR3, paste0(savefilepath,'_num_Pat_withsameNt_persharedClonotype.tsv'), sep='\t', row.names=F)

num_Nt_perCDR3 <- num_pat_withsameNt_perCDR3 %>% group_by(Clonotype) %>% summarise(num_Nt_persharedClonotype = n())
colnames(num_Nt_perCDR3) <- c('Clonotype', 'num_Nt_persharedClonotype')

write.table(num_Nt_perCDR3, paste0(savefilepath,'_num_NT_persharedClonotype.tsv'), sep='\t', row.names=F)

############join NTList_Counts with annotated clonotype list (num_of_samples_present)to plot####

sharedCDR3_numSamples <- binary.df %>% select (Clonotype, Num_Samples_Present, Freq_Samples_Present) %>% filter(Clonotype %in% shared_clonotype)

join<-left_join(num_Nt_perCDR3, sharedCDR3_numSamples, by=c('Clonotype'='Clonotype'))

write.table(join, paste0(savefilepath,'_df_to_plot_convergentRecombination_removeHealthy.tsv'), sep='\t', row.names=F)
