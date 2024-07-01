# new general run

rm(list=ls()) #clears all variables

objects() # clear all objects


## load libs

lib<-c("tidyverse","DESeq2","edgeR","GenomicFeatures", 'scales', 'dplyr')
lapply(lib,library,character.only=T)


setwd("/home/sebastiano/Scrivania/DOTTORATO\ +\ TESI/Progetto_dottorato_Bombina-pachypus/script e file exon selection")


# miscellaneous #########################################

all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths<- read.table("all_quant_htseq_miscellaneous_TPM_Norm_1_reads_filt_averaged_exon_lengths.csv",h=T,sep=",")

filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths <- all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths[all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$exon_length > 150, ]

quantiles <- quantile(filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$average_TPM, probs = c(0.2, 0.4, 0.6, 0.8))

labels <- c("1st quantile", "2nd quantile", "3rd quantile", "4th quantile", "5th quantile")

filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$quantiles <- cut(filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$average_TPM, breaks = c(-Inf, quantiles, Inf), labels = labels, include.lowest = TRUE)

write.table(filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths,"TPM_Norm_1_reads_filt_averaged_exon_lengths>150_by_quantile_miscellaneous.csv", row.names = FALSE, sep=",")



# brain #########################################

all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths<- read.table("all_quant_htseq_brain_TPM_Norm_1_reads_filt_averaged_exon_lengths.csv",h=T,sep=",")

filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths <- all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths[all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$exon_length > 150, ]

quantiles <- quantile(filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$average_TPM, probs = c(0.2, 0.4, 0.6, 0.8))

labels <- c("1st quantile", "2nd quantile", "3rd quantile", "4th quantile", "5th quantile")

filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$quantiles <- cut(filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$average_TPM, breaks = c(-Inf, quantiles, Inf), labels = labels, include.lowest = TRUE)

write.table(filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths,"TPM_Norm_1_reads_filt_averaged_exon_lengths>150_by_quantile_brain.csv", row.names = FALSE, sep=",")



# gonads #########################################

all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths<- read.table("all_quant_htseq_gonads_TPM_Norm_1_reads_filt_averaged_exon_lengths.csv",h=T,sep=",")

filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths <- all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths[all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$exon_length > 150, ]

quantiles <- quantile(filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$average_TPM, probs = c(0.2, 0.4, 0.6, 0.8))

labels <- c("1st quantile", "2nd quantile", "3rd quantile", "4th quantile", "5th quantile")

filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$quantiles <- cut(filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths$average_TPM, breaks = c(-Inf, quantiles, Inf), labels = labels, include.lowest = TRUE)

write.table(filtered_all_quant_htseq_TPM_Norm_reads_filt_averaged_exon_lengths,"TPM_Norm_1_reads_filt_averaged_exon_lengths>150_by_quantile_gonads.csv", row.names = FALSE, sep=",")

