# new general run

rm(list=ls()) #clears all variables
objects() # clear all objects


## load libs

lib<-c("tidyverse","DESeq2","edgeR","GenomicFeatures", 'scales', 'dplyr')
lapply(lib,library,character.only=T)



setwd("/home/sebastiano/Scrivania/DOTTORATO\ +\ TESI/Progetto_dottorato_Bombina-pachypus/script e file exon selection")


gene_counts_exon_lenght<- read.table("all_quant_htseq_exon_length.csv", h=T, sep=",")


# subsampling gonads #########################################

gonads_samples <- c("exon", "BP87GD", "BP91GD", "BP92GD", "BP95GD", "BP100GD", "BP103GD", "BP402GD", "scaffold", "exon_length")
gene_counts_exon_lenght_gonads <- gene_counts_exon_lenght[, gonads_samples]

write.table(gene_counts_exon_lenght_gonads,"all_quant_htseq_exon_length--nonunique all_gonads.csv", col.names = TRUE, row.names = FALSE, sep="," )


# subsampling brain #########################################

brain_samples <- c("exon", "BP87BR", "BP91BR", "BP92BR", "BP95BR", "BP100BR", "BP103BR", "BP402BR", "scaffold", "exon_length")
gene_counts_exon_lenght_brain <- gene_counts_exon_lenght[, brain_samples]

write.table(gene_counts_exon_lenght_brain,"all_quant_htseq_exon_length--nonunique all_brain.csv", col.names = TRUE, row.names = FALSE, sep="," )


# subsampling miscellaneous #########################################

miscellaneous_samples <- c("exon", "BP87BR", "BP92GD", "BP92_pool_L2", "BP92_pool_L4", "BP103BR", "BP402GD", "scaffold", "exon_length")
gene_counts_exon_lenght_miscellaneous <- gene_counts_exon_lenght[, miscellaneous_samples]

write.table(gene_counts_exon_lenght_miscellaneous,"all_quant_htseq_exon_length--nonunique all_miscellaneous.csv", col.names = TRUE, row.names = FALSE, sep="," )





# miscellaneous #########################################

gene_counts_miscellaneous <- read.table("all_quant_htseq_exon_length--nonunique all_miscellaneous.csv", h=T, sep=",")


#filtering at least 1 reads per exon in all samples
filtered_data <- gene_counts_miscellaneous %>% 
  filter_all(all_vars(. >= 1))
write.table(filtered_data,"all_quant_htseq_exon_length_1_reads_filt_miscellaneous.csv", col.names = TRUE, row.names = FALSE, sep=",")




# gonads #########################################

gene_counts_gonads <- read.table("all_quant_htseq_exon_length--nonunique all_gonads.csv", h=T, sep=",")


#filtering at least 1 reads per exon in all samples
filtered_data <- gene_counts_gonads %>% 
  filter_all(all_vars(. >= 1))
write.table(filtered_data,"all_quant_htseq_exon_length_1_reads_filt_gonads.csv", col.names = TRUE, row.names = FALSE, sep="," )



# brain #########################################

gene_counts_brain <- read.table("all_quant_htseq_exon_length--nonunique all_brain.csv", h=T, sep=",")


#filtering at least 1 reads per exon in all samples
filtered_data <- gene_counts_brain %>% 
  filter_all(all_vars(. >= 1))
write.table(filtered_data,"all_quant_htseq_exon_length_1_reads_filt_brain.csv", col.names = TRUE, row.names = FALSE, sep="," )









