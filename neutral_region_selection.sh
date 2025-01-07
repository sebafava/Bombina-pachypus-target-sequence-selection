#!/bin/bash

bedtools=bedtools-2.31

MASTER=/cluster_data/home/genomic/bombina/B.pachypus
genome_bed=$MASTER/genomeDEF/BoPac.JBAT.review.FINAL.bed
genome_index=$MASTER/genomeDEF/BoPac.JBAT.review.FINAL.fa.fai
genes_bed=$MASTER/neutral_region_selection/bed_files/BoPac.genes.bed
TE_bed=$MASTER/neutral_region_selection/bed_files/BoPac.TEannotation_merged.bed
peaks_bed=$MASTER/neutral_region_selection/bed_files/bo_pac_peaks.bed
TE_offset=3000

# Inizializza il file di report
report_file=offset_report.txt
echo "Offset(kb)	Total Length(bp)" > $report_file

# Itera sui valori di k_offset
for k_offset in 5 10 30 50; do
  offset=${k_offset}000

  # Espande le regioni dei geni, dei picchi e degli elementi trasponibili
  $bedtools slop -i $genes_bed -g $genome_index -b $offset > BoPac_expanded_genes_${offset}.bed
  $bedtools slop -i $peaks_bed -g $genome_index -b $offset > BoPac_expanded_peaks_${offset}.bed
  $bedtools slop -i $TE_bed -g $genome_index -b $TE_offset > BoPac_expanded_TE_${TE_offset}.bed

  # Combina ed esclude regioni
  cat BoPac_expanded_genes_${offset}.bed BoPac_expanded_peaks_${offset}.bed BoPac_expanded_TE_${TE_offset}.bed \
    | sort -k1,1 -k2,2n \
    | $bedtools merge > TE_GENE_PEAKS_${offset}_offset.bed

  $bedtools subtract -a $genome_bed -b TE_GENE_PEAKS_${offset}_offset.bed > no-TE_no-GENE_no-PEAKS_${offset}_offset.bed

  # Merging delle regioni con coordinate sovrapposte
  $bedtools merge -i no-TE_no-GENE_no-PEAKS_${offset}_offset.bed -c 2 -o collapse > no-TE_no-GENE_no-PEAKS_${offset}_offset_merged.bed

  # Filtra per lunghezza minima di 150 bp
  awk '($3 - $2) >= 150' no-TE_no-GENE_no-PEAKS_${offset}_offset_merged.bed > no-TE_no-GENE_no-PEAKS_${offset}_offset_merged_filtered.bed

  # Calcola l'estensione totale delle regioni filtrate
  total_length=$(awk '{sum += $3 - $2} END {print sum}' no-TE_no-GENE_no-PEAKS_${offset}_offset_merged_filtered.bed)

  # Aggiungi il risultato al report
  echo -e "${k_offset}K\t$total_length" >> $report_file

  # Crea una directory per il valore corrente di k_offset e sposta i file
  mkdir -p ${k_offset}K
  mv *$offset* ${k_offset}K

done

# Output del report finale
echo "Report generato in $report_file"





############################# produrre un file fasta sulla base delle coordinate delle regioni neutrali selezionate (TE_offset=3000, GENE and PEAKS offset= 30000)

genome=/cluster_data/home/genomic/bombina/B.pachypus/genomeDEF/BoPac.JBAT.review.FINAL.fa
bedtools=bedtools-2.31

$bedtools getfasta -fi $genome -bed ./no-TE_no-GENE_no-PEAKS_30000_offset_merged_filtered.bed -fo no-TE_no-GENE_no-PEAKS_30000_offset_merged_filtered.fa -name -s



