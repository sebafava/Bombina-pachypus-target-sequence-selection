# new general run

rm(list=ls()) #clears all variables

objects() # clear all objects


## load libs

lib<-c("tidyverse","DESeq2","edgeR","GenomicFeatures", 'scales', 'dplyr')
lapply(lib,library,character.only=T)

setwd("/home/sebastiano/Scrivania/DOTTORATO\ +\ TESI/Progetto_dottorato_Bombina-pachypus/script e file exon selection")



# miscellaneous #########################################

TPM_Norm_exon_length <- read.table("TPM_Norm_1_reads_filt_averaged_exon_lengths>150_by_quantile_miscellaneous.csv",h=T, sep=",")


# Ordinare il dataframe in ordine decrescente di average_TPM, calcolare la somma cumulativa della colonna exon_length, filtrare le righe fino a che la somma cumulativa non raggiunge 7 milioni
data_selected <- TPM_Norm_exon_length %>%
  arrange(desc(average_TPM)) %>%
  mutate(cumulative_length = cumsum(exon_length)) %>%
  filter(cumulative_length <= 7000000)

write.table(data_selected,"TPM_Norm_1_reads_filt_averaged_exon_lengths>150_HIGH_EXP_7Mb_spam_miscellaneous.csv", row.names = FALSE, sep=",")



# Ordinare il dataframe in ordine crescente di average_TPM
data_sorted <- TPM_Norm_exon_length %>%
  arrange(average_TPM)

# Calcolare la somma cumulativa della colonna exon_length
data_sorted <- data_sorted %>%
  mutate(cumulative_length = cumsum(exon_length))

# Filtrare le righe fino a che la somma cumulativa non raggiunge 7 milioni
data_filtered <- data_sorted %>%
  filter(cumulative_length <= 7000000)


write.table(data_filtered,"TPM_Norm_1_reads_filt_averaged_exon_lengths>150_LOW_EXP_7Mb_spam_miscellaneous.csv", row.names = FALSE, sep=",")




# brain #########################################

TPM_Norm_exon_length <- read.table("TPM_Norm_1_reads_filt_averaged_exon_lengths>150_by_quantile_brain.csv", h=T, sep=",")


# Ordinare il dataframe in ordine decrescente di average_TPM, calcolare la somma cumulativa della colonna exon_length, filtrare le righe fino a che la somma cumulativa non raggiunge 7 milioni
data_selected <- TPM_Norm_exon_length %>%
  arrange(desc(average_TPM)) %>%
  mutate(cumulative_length = cumsum(exon_length)) %>%
  filter(cumulative_length <= 7000000)

write.table(data_selected,"TPM_Norm_1_reads_filt_averaged_exon_lengths>150_HIGH_EXP_7Mb_spam_brain.csv", row.names = FALSE, sep=",")



# Ordinare il dataframe in ordine crescente di average_TPM
data_sorted <- TPM_Norm_exon_length %>%
  arrange(average_TPM)

# Calcolare la somma cumulativa della colonna exon_length
data_sorted <- data_sorted %>%
  mutate(cumulative_length = cumsum(exon_length))

# Filtrare le righe fino a che la somma cumulativa non raggiunge 7 milioni
data_filtered <- data_sorted %>%
  filter(cumulative_length <= 7000000)


write.table(data_filtered,"TPM_Norm_1_reads_filt_averaged_exon_lengths>150_LOW_EXP_7Mb_spam_brain.csv", row.names = FALSE, sep=",")




# gonads #########################################


TPM_Norm_exon_length <- read.table("TPM_Norm_1_reads_filt_averaged_exon_lengths>150_by_quantile_gonads.csv", h=T, sep=",")


# Ordinare il dataframe in ordine decrescente di average_TPM, calcolare la somma cumulativa della colonna exon_length, filtrare le righe fino a che la somma cumulativa non raggiunge 7 milioni
data_selected <- TPM_Norm_exon_length %>%
  arrange(desc(average_TPM)) %>%
  mutate(cumulative_length = cumsum(exon_length)) %>%
  filter(cumulative_length <= 7000000)

write.table(data_selected,"TPM_Norm_1_reads_filt_averaged_exon_lengths>150_HIGH_EXP_7Mb_spam_gonads.csv", row.names = FALSE, sep=",")



# Ordinare il dataframe in ordine crescente di average_TPM
data_sorted <- TPM_Norm_exon_length %>%
  arrange(average_TPM)

# Calcolare la somma cumulativa della colonna exon_length
data_sorted <- data_sorted %>%
  mutate(cumulative_length = cumsum(exon_length))

# Filtrare le righe fino a che la somma cumulativa non raggiunge 7 milioni
data_filtered <- data_sorted %>%
  filter(cumulative_length <= 7000000)


write.table(data_filtered,"TPM_Norm_1_reads_filt_averaged_exon_lengths>150_LOW_EXP_7Mb_spam_gonads.csv", row.names = FALSE, sep=",")
