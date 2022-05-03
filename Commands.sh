##################################
## Supplementary figures 1 to 7 ##
##################################

##########################
## Mutation frequencies
rm -rf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Mutation_frequency/*
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Mutation_frequency_analysis/RData/* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Mutation_frequency


####################
## Mutation rates
rm -rf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Mutation_rate/*
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Mutation_rate_analysis/RData/* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Mutation_rate





###################################
## Supplementary figures 8 and 9 ##
###################################

## These figures were manually created:
##
## Fig S8: /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/GeneHancer/TFBS_CRE_association.odp
##
## Fig S9: /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Model/Xseq_explained.odp





#####################################
## Supplementary figures 10 and 11 ##
#####################################

###############################################################
## Fig S10: Highlighted genes: BRCA-US, HNSC-US, LIHC-US, and LUSC-US
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/plots/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P1.pdf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P1.pdf

convert -density 500 -trim /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P1.pdf -quality 350 /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P1.jpeg


###############################################################
## Fig S11: Highlighted genes: LUAD-US, STAD-US, and UCEC-US
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/plots/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P2.pdf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P2.pdf

convert -density 500 -trim /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P2.pdf -quality 350 /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P2.jpeg





#############################
## Supplementary figure 12 ##
#############################

###########################################################################################
## Fig S12: number of mutation in TFBS and Exons across highlighted genes in all cohorts
mkdir -p /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Nb_TFBS_vs_Exonic_mutations ;
rm -rf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Nb_TFBS_vs_Exonic_mutations/* ;
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/protein_coding_genes_data/Nb_TFBS_vs_Exonic_mutations/rdata/Analysis_samples_w_TFBS_Exonic_mutations.rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Nb_TFBS_vs_Exonic_mutations/Analysis_samples_w_TFBS_Exonic_mutations.rdata





#############################
## Supplementary figure 13 ##
#############################

###########################################################
## Fig S13: Aggregated networks with all predicted genes 

for COHORT in BRCA HNSC LIHC LUAD LUSC STAD UCEC
do
	echo $COHORT
	mkdir -p /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Aggregated_networks/$COHORT ;
	rm -rf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Aggregated_networks/$COHORT/* ;
	rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/$COHORT-US/premiRNA/xseq/trans/Rdata/Aggregated_network_$COHORT-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Aggregated_networks/$COHORT/Aggregated_network_$COHORT-US.Rdata
done





####################################
## Supplementary figures 14 to 17 ##
####################################

## S14: GO Biological Process 2018 enrichR heatmaps
## S15: KEGG 2019 enrichR heatmaps
## S16: WikiPathways human 2019 enrichR heatmaps
## S17: Panther 2016 enrichR heatmaps

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/enrichR_analysis/RData/*.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/FEA_heatmaps/





#############################
## Supplementary figure 18 ##
#############################
 
####################################################
## Fig S18: pan-cancer highlighted miRNAS heatmap 
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/plots/xseq_trans_pancancer_heatmap_premiRNA_analysis_mutation_type_TFBS_results_summary_plot.pdf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_premiRNA_analysis_mutation_type_TFBS_results_summary_plot.pdf 

convert -density 500 -trim /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_premiRNA_analysis_mutation_type_TFBS_results_summary_plot.pdf  -quality 350 /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_premiRNA_analysis_mutation_type_TFBS_results_summary_plot.jpeg





#############################
## Supplementary figure 19 ##
#############################

#############################################
## Fig S19: Dysregulation network fraction

rm -rf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/BRCA-US/
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/BRCA-US/premiRNA/xseq/trans/Rdata/* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/BRCA-US/

rm -rf /home/jiamondra/Documents/PostDoc/Mathelier_lab/Manuscript/Supp_Material/xseq_ranked_genes/HNSC-US/
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/HNSC-US/premiRNA/xseq/trans/Rdata/* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/HNSC-US/

rm -rf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LIHC-US/
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LIHC-US/premiRNA/xseq/trans/Rdata/* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LIHC-US/

rm -rf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LUAD-US/
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LUAD-US/premiRNA/xseq/trans/Rdata/* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LUAD-US/

rm -rf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LUSC-US/
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LUSC-US/premiRNA/xseq/trans/Rdata/* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LUSC-US/

rm -rf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/STAD-US/
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/STAD-US/premiRNA/xseq/trans/Rdata/* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/STAD-US/

rm -rf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/UCEC-US/
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/UCEC-US/premiRNA/xseq/trans/Rdata/* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/UCEC-US/





#####################################
## Supplementary figures 20 and 21 ##
#####################################

## S20: BASIS dysregulated genes functional enrichment: KEGG + Wikipathways
## S21: BASIS dysregulated genes functional enrichment: GO + Panther

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/All/miRNA/xseq/trans/targetScan_clean_nb_conn_100_minW_0.8_minscore_0.8/enrichR_analysis/RData/Enriched_terms_plots_BASIS.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/FEA_heatmaps/Enriched_terms_plots_BASIS.RData


#####################################
## Supplementary figures 22 and 23 ##
#####################################

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/All/miRNA/xseq/trans/targetScan_clean_nb_conn_100_minW_0.8_minscore_0.8/enrichR_analysis/RData/Enriched_terms_plots_BASIS.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/FEA_heatmaps/Enriched_terms_plots_BASIS.RData




###################
##  ##
###########
mkdir -p /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_results/miRNA ;

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/*-US/premiRNA/xseq/trans/RData/xseq_trans_premiRNA_all_mutation_types_all_effects_*-US_mirna_network_targetScanClean.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_results/miRNA/


rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/All/miRNA/xseq/trans/targetScan_clean_nb_conn_100_minW_0.8_minscore_0.8/tables/xseq_trans_miRNA_all_mutation_types_all_effects_BASIS_network_targetScan_clean.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_results/miRNA/BASIS_All_miRNA.Rdata


rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/ER+/miRNA/xseq/trans/targetScan_clean_nb_conn_100_minW_0.8_minscore_0.8/tables/xseq_trans_miRNA_all_mutation_types_all_effects_BASIS_network_targetScan_clean.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_results/miRNA/BASIS_ER+_miRNA.Rdata


rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/ER-/miRNA/xseq/trans/targetScan_clean_nb_conn_100_minW_0.8_minscore_0.8/tables/xseq_trans_miRNA_all_mutation_types_all_effects_BASIS_network_targetScan_clean.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_results/miRNA/BASIS_ER-_miRNA.Rdata


mkdir -p /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_results/PCG ;

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/*-US/protein_coding_genes_data/xseq/trans/RData/xseq_trans_protein_coding_all_mutation_types_all_effects_*-US_PPI_network_xseq.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_results/PCG/

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/All/PCG/xseq/trans/tables/xseq_trans_protein_coding_all_mutation_types_all_effects_BASIS_PPI_network_xseq_nb_conn_1000_minW_0.8_min_connW_0.8.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_results/PCG/BASIS_All_PCG.Rdata

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/ER+/PCG/xseq/trans/tables/xseq_trans_protein_coding_all_mutation_types_all_effects_BASIS_PPI_network_xseq_nb_conn_1000_minW_0.8_min_connW_0.8.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_results/PCG/BASIS_ER+_PCG.Rdata

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/ER-/PCG/xseq/trans/tables/xseq_trans_protein_coding_all_mutation_types_all_effects_BASIS_PPI_network_xseq_nb_conn_1000_minW_0.8_min_connW_0.8.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_results/PCG/BASIS_ER-_PCG.Rdata



###############
## Figure 2A ##
###############
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/plots/xseq_trans_pancancer_heatmap_PCG_all_mutation_types_triangles.pdf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Main_manuscript/Fig2/xseq_trans_pancancer_heatmap_PCG_all_mutation_types_triangles.pdf

convert -density 500 -trim /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Main_manuscript/Fig2/xseq_trans_pancancer_heatmap_PCG_all_mutation_types_triangles.pdf -quality 350 /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Main_manuscript/Fig2/xseq_trans_pancancer_heatmap_PCG_all_mutation_types_triangles.jpeg


######################
## Figures 5A and B ##
######################

## 5A
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/plots/xseq_trans_pancancer_heatmap_premiRNA_analysis_mutation_type_TFBS_results_summary_plot1.pdf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Main_manuscript/Fig5/Fig5A.pdf

convert -density 500 -trim /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Main_manuscript/Fig5/Fig5A.pdf -quality 350 /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Main_manuscript/Fig5/Fig5A.jpeg

## 5B
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LIHC-US/premiRNA/xseq/trans/dysregulation_heatmaps/hsa-mir-20a-5p::hsa-mir-20a_TFBS_dysregulation_heatmap.pdf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA/hsa-mir-20a-5p::hsa-mir-20a_TFBS_dysregulation_heatmap.pdf

convert -density 500 -trim /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA/hsa-mir-20a-5p::hsa-mir-20a_TFBS_dysregulation_heatmap.pdf -quality 350 /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Main_manuscript/Fig5/dysreg_heatmap_example.jpeg




###########################
## Numbers in manuscript ##
###########################

## Number of mature miRNAs (pancancer)
more  /storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/filtered_networks/All_geneList_Pancancer_analysis.txt  | cut -f2,4 | grep '::' | sed 's/::/\t/g' | sort | uniq | cut -f1 | sort | uniq | wc -l

## Number of miRNA precursors (pancancer)
more  /storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/filtered_networks/All_geneList_Pancancer_analysis.txt  | cut -f2,4 | grep '::' | sed 's/::/\t/g' | sort | uniq | cut -f2 | sort | uniq | wc -l

## Mature miRNA - Nb of cohorts
more  /storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/filtered_networks/All_geneList_Pancancer_analysis.txt  | cut -f2,4 | grep '::' | sed 's/::/\t/g' | sort | uniq | cut -f1,3 | cut -f1 | uniq -c | sort -h -r 

## Precursosr in more than 5 cohorts
more  /storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/filtered_networks/All_geneList_Pancancer_analysis.txt  | cut -f2,4 | grep '::' | sed 's/::/\t/g' | sort | uniq | cut -f1,2,3 | cut -f1,2 | uniq -c | sort -h -r | cut -d '' -f4




















more /storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/protein_coding_genes_data/xseq/trans/tables/xseq_trans_protein_coding_results_with_annotations.tab | grep -v LoF | cut -f1 | grep -v gene | sort -u | wc -l

## Count number of mature miRNAs predicted in TCGA:
more more /storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/tables/xseq_trans_table_all_analysis_all_cancer_projects_cat_premiRNA.tab | cut -f 2 | perl -lne  '$_ =~ s/::.+$//gi; print $_' | sort | uniq | grep -v gene | wc -l

## Count number of pre-miRNAs predicted in TCGA: 
more /storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/tables/xseq_trans_premiRNA_results_with_annotations.tab | cut -f1 | grep -v gene | sort -u | wc -l

## Count number of highlighted PCGs in TCGA cohorts
more /storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/protein_coding_genes_data/xseq/trans/tables/xseq_trans_protein_coding_results_with_annotations.tab | grep -v LoF | cut -f1 | grep -v gene | sort -u | wc -l


#############
## Table 1 ##
#############
ls /storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/*-US/tables/specimens/Specimen_Donor_information_Tumour_*-US_WGS.tab | grep -v Blood | xargs cat | cut -f1,2,3,4,7 > /storage/mathelierarea/processed/jamondra/Projects/dysmir/tmp/Supp_table_1_selected_TCGA_samples.tab


#############
## Table 3 ##
#############
#ls /storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/*-US/premiRNA/xseq/trans/filtered_networks/Dysregulated_network_fraction_*-US.tab | xargs cat > /storage/mathelierarea/processed/jamondra/Projects/dysmir/tmp/Supp_table_3_drivers_and_dysregulated_genes.tab

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/tmp/Supp_table_3_drivers_and_dysregulated_genes.tab /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscript/Tables/Supp_table_3_predicted_genes_and_dysregulated_targets.tab


#############
## Table 5 ##
#############
#ls /storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/*-US/premiRNA/xseq/trans/filtered_networks/Driver_and_dysregulated_genes_summary_*-US.tab | xargs cat > /storage/mathelierarea/processed/jamondra/Projects/dysmir/tmp/Supp_table_5_summary_number_highlighted_drivers_and_dysreg_targets.tab

rsync -rptuvl jamondra@biotin3.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/tmp/Supp_table_5_summary_number_highlighted_drivers_and_dysreg_targets.tab /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscript/Tables/Supp_table_5_summary_number_highlighted_drivers_and_dysreg_targets.tab


###################################
## Pancancer landscape paragraph ##
###################################

## Nb LoF (protein coding genes)
more /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscript/Tables/Additional_file_3_drivers_and_dysregulated_genes.tab | grep LoF | cut -f2 | grep -v '::' | sort | uniq -c |wc -l

## Nb TFBS (protein coding genes)
more /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscript/Tables/Additional_file_3_drivers_and_dysregulated_genes.tab | grep TFBS | cut -f2 | grep -v '::' | sort | uniq -c |wc -l

## Genes predicted by LoF and TFBSs
more /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscript/Tables/Additional_file_3_drivers_and_dysregulated_genes.tab  | cut -f2,3 | grep -v '::' | sort | uniq -c

## Total number pf predicted protein coding genes
more /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscript/Tables/Additional_file_3_drivers_and_dysregulated_genes.tab |  cut -f2 | grep -v '::' | sort | uniq | wc -l

###############################
## Pancancer miRNAs analysis ##
###############################

## Number of mature miRNAs
more /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscript/Tables/Additional_file_3_drivers_and_dysregulated_genes.tab |  cut -f2 | grep '::' | perl -lne ' $_ =~ s/^.+:://gi; print $_;' | sort | uniq | wc -l

## Number of precursors







##################################################
## Predicted miRNA target distribution
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/RData/Nb_targets_distribution.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Nb_targets_miRNA/Nb_targets_distribution.RData 

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/RData/Nb_targets_vs_Nb_cohorts.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Nb_targets_miRNA/Nb_targets_vs_Nb_cohorts.RData





#####################
## Gene rank plots ##
#####################
rm -rf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/BRCA-US/
rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/BRCA-US/premiRNA/xseq/trans/Rdata/Rank* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/BRCA-US/

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/HNSC-US/premiRNA/xseq/trans/Rdata/Rank* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/HNSC-US/

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LIHC-US/premiRNA/xseq/trans/Rdata/Rank* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LIHC-US/

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LUAD-US/premiRNA/xseq/trans/Rdata/Rank* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LUAD-US/

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LUSC-US/premiRNA/xseq/trans/Rdata/Rank* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LUSC-US/

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/STAD-US/premiRNA/xseq/trans/Rdata/Rank* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/STAD-US/

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/UCEC-US/premiRNA/xseq/trans/Rdata/Rank* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/UCEC-US/

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/All/miRNA/xseq/trans/targetScan_clean_nb_conn_10000_minW_0.8_minscore_0.8/RData/Rank* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/BASIS/All/

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/ER+/miRNA/xseq/trans/targetScan_clean_nb_conn_10000_minW_0.8_minscore_0.8/RData/Rank* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/BASIS/ER+/

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/ER-/miRNA/xseq/trans/targetScan_clean_nb_conn_10000_minW_0.8_minscore_0.8/RData/Rank* /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/BASIS/ER-/





############################################
## Dysregulation heatmap examples - miRNA ##
############################################
#rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/BRCA-US/premiRNA/xseq/trans/Rdata/Dysregulation_heatmaps_BRCA-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA

#rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/HNSC-US/premiRNA/xseq/trans/Rdata/Dysregulation_heatmaps_HNSC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA

#rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LIHC-US/premiRNA/xseq/trans/Rdata/Dysregulation_heatmaps_LIHC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA

#rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LUAD-US/premiRNA/xseq/trans/Rdata/Dysregulation_heatmaps_LUAD-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA

#rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LUSC-US/premiRNA/xseq/trans/Rdata/Dysregulation_heatmaps_LUSC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA

#rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/STAD-US/premiRNA/xseq/trans/Rdata/Dysregulation_heatmaps_STAD-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA

#rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/UCEC-US/premiRNA/xseq/trans/Rdata/Dysregulation_heatmaps_UCEC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA




############################################
## Examples dysregulation heatmaps - PCGs ##
############################################

#mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/BRCA-US/Dysregulation_heatmaps_BRCA-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/c/PCG/Dysregulation_heatmaps_BRCA-US.Rdata

#mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LUAD-US/Dysregulation_heatmaps_LUAD-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/PCG/Dysregulation_heatmaps_LUAD-US.Rdata

#mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LUSC-US/Dysregulation_heatmaps_LUSC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/PCG/Dysregulation_heatmaps_LUSC-US.Rdata


##############################################
## Examples dysregulation heatmaps - miRNAs ##
##############################################

#mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/PCG/Dysregulation_heatmaps_BRCA-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA/Dysregulation_heatmaps_BRCA-US.Rdata

#mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LIHC-US/Dysregulation_heatmaps_LUSC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA/Dysregulation_heatmaps_LUSC-US.Rdata

#mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/UCEC-US/Dysregulation_heatmaps_UCEC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_heatmaps_examples/miRNA/Dysregulation_heatmaps_UCEC-US.Rdata


####################################
## Fraction dysregulated networks ##
####################################

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/BRCA-US/premiRNA/xseq/trans/Rdata/Dysregulated_network_fraction_BRCA-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Fraction_dysreg_net/Dysregulated_network_fraction_BRCA-US.Rdata

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/HNSC-US/premiRNA/xseq/trans/Rdata/Dysregulated_network_fraction_HNSC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Fraction_dysreg_net/Dysregulated_network_fraction_HNSC-US.Rdata

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LIHC-US/premiRNA/xseq/trans/Rdata/Dysregulated_network_fraction_LIHC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Fraction_dysreg_net/Dysregulated_network_fraction_LIHC-US.Rdata

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LUAD-US/premiRNA/xseq/trans/Rdata/Dysregulated_network_fraction_LUAD-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Fraction_dysreg_net/Dysregulated_network_fraction_LUAD-US.Rdata

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/LUSC-US/premiRNA/xseq/trans/Rdata/Dysregulated_network_fraction_LUSC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Fraction_dysreg_net/Dysregulated_network_fraction_LUSC-US.Rdata

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/STAD-US/premiRNA/xseq/trans/Rdata/Dysregulated_network_fraction_STAD-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Fraction_dysreg_net/Dysregulated_network_fraction_STAD-US.Rdata

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/UCEC-US/premiRNA/xseq/trans/Rdata/Dysregulated_network_fraction_UCEC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Fraction_dysreg_net/Dysregulated_network_fraction_UCEC-US.Rdata


######################################
## Number of targets miRNA networks ##
######################################

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/RData/Nb_targets_distribution.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Nb_targets_miRNA/Nb_targets_distribution.RData

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/RData/Nb_targets_vs_Nb_cohorts.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Nb_targets_miRNA/Nb_targets_vs_Nb_cohorts.RData


##############################################################
## dysregulated networks intersection among miRNAs and PCGs ##
##############################################################

mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/BRCA-US/Dysreg_networks_intersection_BRCA-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_net_intersection/Dysreg_networks_intersection_BRCA-US.Rdata

mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/HNSC-US/Dysreg_networks_intersection_HNSC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_net_intersection/Dysreg_networks_intersection_HNSC-US.Rdata

mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LIHC-US/Dysreg_networks_intersection_LIHC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_net_intersection/Dysreg_networks_intersection_LIHC-US.Rdata

mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LUAD-US/Dysreg_networks_intersection_LUAD-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_net_intersection/Dysreg_networks_intersection_LUAD-US.Rdata

mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/LUSC-US/Dysreg_networks_intersection_LUSC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_net_intersection/Dysreg_networks_intersection_LUSC-US.Rdata

mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/STAD-US/Dysreg_networks_intersection_STAD-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_net_intersection/Dysreg_networks_intersection_STAD-US.Rdata

mv /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/xseq_ranked_genes/UCEC-US/Dysreg_networks_intersection_UCEC-US.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Dysreg_net_intersection/Dysreg_networks_intersection_UCEC-US.Rdata



###################################
## ER composition: BASIS vs TCGA ##
###################################

rsync -rptuvl jamondra@biotin3.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/TEST/ER_MATCH/tables/BASIS_ER_status_tab.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/ER_status/BASIS_ER_status_tab.RData

rsync -rptuvl jamondra@biotin3.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/TEST/ER_MATCH/tables/TCGA_ER_status_tab.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/ER_status/TCGA_ER_status_tab.RData


















<!-- ```{r Nb_TFBS_per_gene_PCG, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="", fig.align = "center"} -->
  <!-- # Distribution of associated TFBSs among predicted, non-predicted, TFs and cancer-related protein-coding genes represented as A) Cumulative Distribution Function, and B) Histogram of frequencies. -->
  <!-- cohorts.all <- c(cohorts, "BASIS_All", "BASIS_ER+", "BASIS_ER-") -->
    <!-- selected.genes.list  <- vector(mode = "list", length = length(cohorts.all)) -->
      <!-- genes.in.cohort.list <- vector(mode = "list", length = length(cohorts.all)) -->
        <!-- Nb.mut.samples.pcg   <- list() -->
          <!-- for (co in cohorts.all) { -->
              
              <!--   if (!grepl(co, pattern = "BASIS")) { -->
                  
                  <!--     co <- paste0(co, "-US") -->
                    
                    <!--     load(paste0("xseq_results/PCG/xseq_trans_protein_coding_all_mutation_types_all_effects_", co, "_PPI_network_xseq.Rdata")) ## xseq.pred.all.concat -->
                  
                  <!--   } else { -->
                      
                      <!--     load(paste0("xseq_results/PCG/", co, "_PCG.Rdata")) -->
                      
                      <!--   } -->
              
              <!--   ## Parse xseq output -->
              <!--   xseq.pred.pcg <- xseq.pred.all.concat %>%  -->
                <!--                     filter(mut_effect == "TFBS") %>%  -->
                <!--                     select(sample, gene, prob_mut, prob_gene) %>%  -->
                <!--                     distinct() %>%  -->
                <!--                     arrange(desc(prob_gene))  -->
                
                
                <!--   Nb.mut.samples.pcg[[co]] <- xseq.pred.pcg %>%  -->
                  <!--                                 group_by(gene) %>%  -->
                  <!--                                 summarise(Nb_samples = n()) %>%  -->
                  <!--                                 arrange(desc(Nb_samples)) %>%  -->
                  <!--                                 rename(Gene = gene) %>%  -->
                  <!--                                 mutate(Cohort = co) -->
                  
                  
                  <!--   genes.in.cohort.list[[co]] <- unique(xseq.pred.pcg$gene) -->
                    <!--   selected.genes.list[[co]]  <- unique(subset(xseq.pred.pcg, prob_gene >= xseq.th.pcg[[co]])$gene) -->
                      <!-- } -->
          
          <!-- selected.pcg.list  <- selected.genes.list -->
            <!-- pcg.in.cohort.list <- genes.in.cohort.list -->
              
              <!-- predicted.genes.pancancer <- unique(as.vector(unlist(selected.genes.list))) -->
                <!-- xseq.analyzed.genes       <- unique(as.vector(unlist(genes.in.cohort.list))) -->
                  <!-- reported.cancer.genes <- xseq.sup.mat.tab$inToGen_names -->
                    <!-- TF.genes              <- xseq.sup.mat.tab$Human_TF_names -->
                      <!-- ## A subset with only the TFBS associated to the genes analyzed by xseq across the cohorts -->
                      <!-- TFBS.xseq.analyzed.genes.tfbs <- all.associated.features.tab %>%  -->
                        <!--                               filter(Associated_gene %in% xseq.analyzed.genes) %>%  -->
                        <!--                               distinct() -->
                        <!-- # Nb.TFBSs.per.gene.per.TF <- TFBS.xseq.analyzed.genes %>%  -->
                        <!-- #                               distinct() %>%  -->
                        <!-- #                               group_by(Associated_gene, TF) %>%  -->
                        <!-- #                               summarise(Nb_TFBS = n(), .groups = "drop") %>%  -->
                        <!-- #                               arrange(desc(Nb_TFBS)) %>%  -->
                        <!-- #                               mutate(Predicted_gene = ifelse(Associated_gene %in% predicted.genes.pancancer, yes = "Predicted", no = "Non-predicted")) -->
                        <!-- Nb.TFBSs.per.pcg <- TFBS.xseq.analyzed.genes.tfbs %>%  -->
                          <!--                         distinct() %>%  -->
                          <!--                         group_by(Associated_gene) %>%  -->
                          <!--                         summarise(Nb_TFBS = n(), .groups = "drop") %>%  -->
                          <!--                         arrange(desc(Nb_TFBS)) -->
                          <!-- Nb.TFBSs.per.pcg$Predicted     <- ifelse(Nb.TFBSs.per.pcg$Associated_gene %in% predicted.genes.pancancer, yes = "Predicted", no = "") -->
                            <!-- Nb.TFBSs.per.pcg$Non_Predicted <- ifelse(!Nb.TFBSs.per.pcg$Associated_gene %in% predicted.genes.pancancer, yes = "Non-Predicted", no = "") -->
                              <!-- Nb.TFBSs.per.pcg$Cancer <- ifelse(Nb.TFBSs.per.pcg$Associated_gene %in% reported.cancer.genes, yes = "Cancer gene", no = "") -->
                                <!-- Nb.TFBSs.per.pcg$TF <- ifelse(Nb.TFBSs.per.pcg$Associated_gene %in% TF.genes, yes = "TF", no = "") -->
                                  <!-- Nb.TFBSs.per.pcg.classes <- reshape2::melt(Nb.TFBSs.per.pcg, -->
                                                                                    <!--                                id.vars = c("Nb_TFBS", "Associated_gene"), -->
                                                                                    <!--                                measure.vars = c("Predicted", "Non_Predicted", "Cancer", "TF")) %>%  -->
                                    <!--                       filter(value != "") %>%  -->
                                    <!--                       within(rm(variable)) %>%  -->
                                    <!--                       rename(Class = value) -->
                                    <!-- ## Plot 1: ECDF -->
                                    <!-- TFBS.distr.ecdf.plot <- ggplot(Nb.TFBSs.per.pcg.classes, aes(x = Nb_TFBS, colour = Class)) + -->
                                      <!--         stat_ecdf(size = 1.15, pad = FALSE) + -->
                                      <!--         theme_bw() + -->
                                      <!--         scale_color_manual(values = c("#E69F00", "#56B4E9", "#1b9e77", "#d95f02")) +  -->
                                      <!--         labs(x = "Nb of TFBSs assigned to genes", y = "ECDF") + -->
                                      <!--         theme(legend.position = "none")  -->
                                      <!-- ## Plot 2: histogram with real values -->
                                      <!-- TFBS.distr.hist.real.plot <- ggplot(Nb.TFBSs.per.pcg.classes, aes(x = Nb_TFBS, fill = Class)) + -->
                                        <!--         geom_histogram(binwidth = 50, position = "identity") + -->
                                        <!--         scale_fill_manual(values = c("#E69F00", "#56B4E9", "#1b9e77", "#d95f02")) + -->
                                        <!--         theme_bw() + -->
                                        <!--         labs(x = "Nb of TFBSs assigned to genes", y = "Frequency") + -->
                                        <!--         facet_grid(Class ~ .) + -->
                                        <!--         theme(legend.position = "right")  -->
                                        <!-- ## Plot 3: histogram with log2 transformed values -->
                                        <!-- TFBS.distr.hist.log2.plot <- ggplot(Nb.TFBSs.per.pcg.classes, aes(x = log2(Nb_TFBS), fill = Class)) + -->
                                          <!--         geom_histogram(binwidth = 0.1, position = "identity") + -->
                                          <!--         scale_fill_manual(values = c("#E69F00", "#56B4E9", "#1b9e77", "#d95f02")) + -->
                                          <!--         theme_bw() + -->
                                          <!--         facet_grid(Class ~ .) -->
                                          <!-- ## Plot 4: histogram with log10 transformed values -->
                                          <!-- TFBS.distr.hist.log10.plot <- ggplot(Nb.TFBSs.per.pcg.classes, aes(x = log10(Nb_TFBS), fill = Class)) + -->
                                            <!--         geom_histogram(binwidth = 0.1, position = "identity") + -->
                                            <!--         scale_fill_manual(values = c("#E69F00", "#56B4E9", "#1b9e77", "#d95f02")) + -->
                                            <!--         theme_bw() + -->
                                            <!--         facet_grid(Class ~ .) -->
                                            <!-- # predicted.TFBS.values <- subset(Nb.TFBSs.per.pcg, Predicted_gene == "Predicted")$Nb_TFBS -->
                                            <!-- # non.predicted.TFBS.values <- subset(Nb.TFBSs.per.pcg, Predicted_gene == "Non-predicted")$Nb_TFBS -->
                                            <!-- # predicted.TFBS.values     <- log10(predicted.TFBS.values) -->
                                            <!-- # non.predicted.TFBS.values <- log10(non.predicted.TFBS.values) -->
                                            <!-- # fligner.test(x = list(predicted.TFBS.values, non.predicted.TFBS.values)) -->
                                            <!-- #  -->
                                            <!-- #  -->
                                            <!-- # wilcox.test(x = predicted.TFBS.values, y = non.predicted.TFBS.values, alternative = "two.sided", paired = FALSE, conf.int = 0.95) -->
                                            <!-- # enriched.terms.gg <- plot_grid(TFBS.distr.ecdf.plot, -->
                                            <!-- #                                TFBS.distr.hist.real.plot, -->
                                            <!-- #                                ncol = 1, align = "v", labels = c('A)', 'B)'), label_size = 13, -->
                                            <!-- #                                label_colour = "black", hjust = -0.25, vjust = 1.15, -->
                                            <!-- #                                rel_heights = c(1, 1.5), rel_widths = c(1, 1), -->
                                            <!-- #                                axis = "right") -->
                                            <!-- #  -->
                                            <!-- # ## Combine plots -->
                                            <!-- # enriched.terms.gg -->
                                            <!-- ``` -->
                                            
                                            
                                            
                                            
                                            <!-- ```{r Nb_connections_per_PCG, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Distribution of number of associated TFBSs (A and B) and number of targets in the regulatory network (C and D) among predicted, non-predicted, TFs and cancer-related protein-coding genes (aggregated from all evaluated cohorts) represented as Cumulative Distribution Function and Histogram of frequencies. The maximum number of targets of each protein-coding gene in the regulatory network was restricted to 1000.", fig.height = 10, fig.width = 10, fig.align = "center"} -->
                                            <!-- ## Genes in the network that were analyzed by xseq -->
                                            <!-- Nb.targets.per.pcg <- xseq.net.ppi %>%  -->
                                              <!--                                 filter(Gene %in% xseq.analyzed.genes) %>%  -->
                                              <!--                                 group_by(Gene) %>%  -->
                                              <!--                                 summarise(Nb_target = n(), .groups = "drop") %>%  -->
                                              <!--                                 arrange(desc(Nb_target)) -->
                                              <!-- Nb.targets.per.pcg$Predicted     <- ifelse(Nb.targets.per.pcg$Gene %in% predicted.genes.pancancer, yes = "Predicted", no = "") -->
                                                <!-- Nb.targets.per.pcg$Non_Predicted <- ifelse(!Nb.targets.per.pcg$Gene %in% predicted.genes.pancancer, yes = "Non-Predicted", no = "") -->
                                                  <!-- Nb.targets.per.pcg$Cancer        <- ifelse(Nb.targets.per.pcg$Gene %in% reported.cancer.genes, yes = "Cancer gene", no = "") -->
                                                    <!-- Nb.targets.per.pcg$TF            <- ifelse(Nb.targets.per.pcg$Gene %in% TF.genes, yes = "TF", no = "") -->
                                                      <!-- Nb.targets.per.pcg.classes <- reshape2::melt(Nb.targets.per.pcg, -->
                                                                                                          <!--                                         id.vars = c("Gene", "Nb_target"), -->
                                                                                                          <!--                                         measure.vars = c("Predicted", "Non_Predicted", "Cancer", "TF")) %>%  -->
                                                        <!--                                         filter(value != "") %>%  -->
                                                        <!--                                         within(rm(variable)) %>%  -->
                                                        <!--                                         rename(Class = value) -->
                                                        <!-- ## Plot 1: ECDF -->
                                                        <!-- conn.distr.ecdf.plot <- ggplot(Nb.targets.per.pcg.classes, aes(x = Nb_target, colour = Class)) + -->
                                                          <!--         stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                          <!--         theme_bw() + -->
                                                          <!--         scale_color_manual(values = c("#E69F00", "#56B4E9", "#1b9e77", "#d95f02")) +  -->
                                                          <!--         labs(x = "Nb of connections in network - Protein-coding genes", y = "ECDF") + -->
                                                          <!--         theme(legend.position = "none")  -->
                                                          <!-- ## Plot 2: histogram with real values -->
                                                          <!-- conn.distr.hist.real.plot <- ggplot(Nb.targets.per.pcg.classes, aes(x = Nb_target, fill = Class)) + -->
                                                            <!--         geom_histogram(binwidth = 25, position = "identity") + -->
                                                            <!--         scale_fill_manual(values = c("#E69F00", "#56B4E9", "#1b9e77", "#d95f02")) + -->
                                                            <!--         theme_bw() + -->
                                                            <!--         labs(x = "Nb of TFBSs assigned to genes", y = "Frequency") + -->
                                                            <!--         facet_grid(Class ~ .) + -->
                                                            <!--         theme(legend.position = "right")  -->
                                                            <!-- enriched.terms.gg <- plot_grid(TFBS.distr.ecdf.plot, -->
                                                                                                  <!--                                TFBS.distr.hist.real.plot, -->
                                                                                                  <!--                                conn.distr.ecdf.plot, -->
                                                                                                  <!--                                conn.distr.hist.real.plot, -->
                                                                                                  <!--                                ncol = 2, align = "hv", labels = c('A)', 'B)', 'C)', 'D)'), label_size = 13, -->
                                                                                                  <!--                                label_colour = "black", hjust = -0.25, vjust = 1.15, -->
                                                                                                  <!--                                rel_heights = c(1, 1, 1, 1), rel_widths = c(1, 1, 1, 1), -->
                                                                                                  <!--                                axis = "right") -->
                                                              <!-- ## Combine plots -->
                                                              <!-- enriched.terms.gg -->
                                                              <!-- ``` -->
                                                              
                                                              
                                                              
                                                              <!-- ```{r eCDF_plot1_pcg, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Empirical Cumulative Distribution Functions of number of TFBS, number of targets and number of mutations in protein-coding genes, separated by cohorts (part 1).", fig.height = 10, fig.width = 8, fig.align = "center"} -->
                                                              
                                                              <!-- analyzed.pcg.summary <- Nb.targets.per.pcg %>%  -->
                                                                <!--                           select(Gene, Nb_target, Predicted) %>%  -->
                                                                <!--                           rename(Associated_gene = Gene) %>%  -->
                                                                <!--                           left_join(Nb.TFBSs.per.pcg, by = "Associated_gene") %>%  -->
                                                                <!--                           rename(Predicted = Predicted.x, -->
                                                                                                        <!--                                  Gene      = Associated_gene) %>%  -->
                                                                <!--                           select(Gene, Nb_target, Nb_TFBS, Predicted, Cancer, TF) %>%  -->
                                                                <!--                           mutate(Predicted = ifelse(Predicted == "", yes = "Non-Predicted", no = Predicted)) -->
                                                                
                                                                <!-- analyzed.pcg.cohorts.summary <- lapply(Nb.mut.samples.pcg, function(l){ -->
                                                                    <!--                                   merge(l, analyzed.pcg.summary, by = "Gene") %>%  -->
                                                                    <!--                                   reshape2::melt(id.vars = c("Gene", "Nb_target", "Nb_TFBS", "Nb_samples", "Cohort"), -->
                                                                                                                            <!--                                                  measure.vars = c("Predicted", "Cancer", "TF")) %>%  -->
                                                                    <!--                                                  filter(value != "") %>%  -->
                                                                    <!--                                                  within(rm(variable)) %>%  -->
                                                                    <!--                                                  rename(Class = value) -->
                                                                    <!--                                 }) -->
                                                                  <!-- ############### -->
                                                                  <!-- ## Example 1 ## -->
                                                                  <!-- ############### -->
                                                                  <!-- # aa <- rbindlist(analyzed.pcg.cohorts.summary) %>%  -->
                                                                  <!-- #   ggplot(aes(x = Nb_TFBS, colour = Class)) + -->
                                                                  <!-- #                     stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                  <!-- #                     theme_bw() + -->
                                                                  <!-- #                     scale_color_manual(values = cols) +  -->
                                                                  <!-- #                     labs(x = "Nb of TFBSs assigned to protein-coding genes", y = "ECDF") + -->
                                                                  <!-- #                     theme(legend.position = "none") + -->
                                                                  <!-- #                     facet_wrap(. ~ Cohort, ncol = 1) -->
                                                                  <!-- #  -->
                                                                  <!-- #  -->
                                                                  <!-- # bb <- rbindlist(analyzed.pcg.cohorts.summary) %>%  -->
                                                                  <!-- #   ggplot(aes(x = Nb_target, colour = Class)) + -->
                                                                  <!-- #                     stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                  <!-- #                     theme_bw() + -->
                                                                  <!-- #                     scale_color_manual(values = cols) +  -->
                                                                  <!-- #                     labs(x = "Nb of TFBSs assigned to protein-coding genes", y = "ECDF") + -->
                                                                  <!-- #                     theme(legend.position = "none") + -->
                                                                  <!-- #                     facet_wrap(. ~ Cohort, ncol = 1) -->
                                                                  <!-- #  -->
                                                                  <!-- #  -->
                                                                  <!-- # cc <- rbindlist(analyzed.pcg.cohorts.summary) %>%  -->
                                                                  <!-- #   ggplot(aes(x = Nb_samples, colour = Class)) + -->
                                                                  <!-- #                     stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                  <!-- #                     theme_bw() + -->
                                                                  <!-- #                     scale_color_manual(values = cols) +  -->
                                                                  <!-- #                     labs(x = "Nb of TFBSs assigned to protein-coding genes", y = "ECDF") + -->
                                                                  <!-- #                     theme(legend.position = "none") + -->
                                                                  <!-- #                     facet_wrap(. ~ Cohort, ncol = 1) -->
                                                                  <!-- #  -->
                                                                  <!-- # enriched.terms.gg <- plot_grid(aa, -->
                                                                  <!-- #                                bb, -->
                                                                  <!-- #                                cc, -->
                                                                  <!-- #                                ncol = 3, align = "hv", hjust = -0.25, vjust = 1.15, -->
                                                                  <!-- #                                rel_heights = c(1, 1), rel_widths = c(1, 1), -->
                                                                  <!-- #                                axis = "right") -->
                                                                  <!-- # enriched.terms.gg -->
                                                                  
                                                                  <!-- ############### -->
                                                                  <!-- ## Example 2 ## -->
                                                                  <!-- ############### -->
                                                                  <!-- cols <- c("Cancer gene" = "#E69F00", "Non-Predicted" = "#56B4E9", "Predicted" = "#1b9e77", "TF" = "#d95f02") -->
                                                                    <!-- nn <- -->
                                                                      <!-- reshape2::melt(rbindlist(analyzed.pcg.cohorts.summary), -->
                                                                                            <!--                id.vars = c("Gene", "Cohort", "Class"), -->
                                                                                            <!--                measure.vars = c("Nb_target", "Nb_TFBS", "Nb_samples")) %>% -->
                                                                      <!--                filter(value != "") -->
                                                                      <!-- cols.subset <- c("BASIS_All", "BASIS_ER+", "BASIS_ER-", "BRCA-US", "HNSC-US") -->
                                                                        <!-- nn %>% -->
                                                                        <!--   filter(Cohort %in% cols.subset) %>%  -->
                                                                        <!--   ggplot(aes(x = value, colour = Class)) + -->
                                                                        <!--                     stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                        <!--                     theme_bw() + -->
                                                                        <!--                     scale_color_manual(values = cols) + -->
                                                                        <!--                     labs(x = "Nb of TFBSs assigned to protein-coding genes", y = "ECDF") + -->
                                                                        <!--                     theme(legend.position = "bottom") + -->
                                                                        <!--                     facet_grid(vars(Cohort), vars(variable), scales = "free") -->
                                                                        <!-- # View(analyzed.pcg.cohorts.summary[[10]]) -->
                                                                        <!-- # cols.point <- c("Non-Predicted" = "#E69F00", "Predicted" = "#1b9e77") -->
                                                                        <!-- #  -->
                                                                        <!-- #  -->
                                                                        <!-- # cols <- c("Cancer gene" = "#E69F00", "Non-Predicted" = "#56B4E9", "Predicted" = "#1b9e77", "TF" = "#d95f02") -->
                                                                        <!-- #  -->
                                                                        <!-- # ecdf.tfbs <- ggplot(analyzed.pcg.cohorts.summary[[10]], aes(x = Nb_TFBS, colour = Class)) + -->
                                                                        <!-- #                     stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                        <!-- #                     theme_bw() + -->
                                                                        <!-- #                     scale_color_manual(values = cols) +  -->
                                                                        <!-- #                     labs(x = "Nb of TFBSs assigned to protein-coding genes", y = "ECDF") + -->
                                                                        <!-- #                     theme(legend.position = "none")  -->
                                                                        <!-- #  -->
                                                                        <!-- # ecdf.targets <- ggplot(analyzed.pcg.cohorts.summary[[10]], aes(x = Nb_target, colour = Class)) + -->
                                                                        <!-- #                        stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                        <!-- #                        theme_bw() + -->
                                                                        <!-- #                        scale_color_manual(values = cols) +  -->
                                                                        <!-- #                        labs(x = "Nb of targets of protein-coding genes", y = "ECDF") + -->
                                                                        <!-- #                        theme(legend.position = "none")  -->
                                                                        <!-- #  -->
                                                                        <!-- # ecdf.samples <- ggplot(analyzed.pcg.cohorts.summary[[10]], aes(x = Nb_samples, colour = Class)) + -->
                                                                        <!-- #                        stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                        <!-- #                        theme_bw() + -->
                                                                        <!-- #                        scale_color_manual(values = cols) +  -->
                                                                        <!-- #                        labs(x = "Nb of mutated samples", y = "ECDF") + -->
                                                                        <!-- #                        theme(legend.position = "right")  -->
                                                                        <!-- #  -->
                                                                        <!-- #  -->
                                                                        <!-- # enriched.terms.gg <- plot_grid(ecdf.tfbs, -->
                                                                        <!-- #                                ecdf.targets, -->
                                                                        <!-- #                                ecdf.samples, -->
                                                                        <!-- #                                ncol = 3, align = "hv", labels = c('A)', 'B)', 'C)'), label_size = 20, -->
                                                                        <!-- #                                label_colour = "black", hjust = -0.25, vjust = 1.15, -->
                                                                        <!-- #                                rel_heights = c(1, 1, 1), rel_widths = c(1, 1, 1), -->
                                                                        <!-- #                                axis = "right") -->
                                                                        <!-- #  -->
                                                                        <!-- # ## Combine plots -->
                                                                        <!-- # enriched.terms.gg -->
                                                                        <!-- ``` -->
                                                                        
                                                                        
                                                                        
                                                                        <!-- ```{r eCDF_plot2, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Empirical Cumulative Distribution Functions of number of TFBS, number of targets and number of mutations in protein-coding genes, separated by cohorts (part 1)", fig.height = 10, fig.width = 8, fig.align = "center"} -->
                                                                        <!-- cols <- c("Cancer gene" = "#E69F00", "Non-Predicted" = "#56B4E9", "Predicted" = "#1b9e77", "TF" = "#d95f02") -->
                                                                          <!-- nn <- -->
                                                                            <!-- reshape2::melt(rbindlist(analyzed.pcg.cohorts.summary), -->
                                                                                                  <!--                id.vars = c("Gene", "Cohort", "Class"), -->
                                                                                                  <!--                measure.vars = c("Nb_target", "Nb_TFBS", "Nb_samples")) %>% -->
                                                                            <!--                filter(value != "") -->
                                                                            <!-- cols.subset <- c("LIHC-US", "LUAD-US", "LUSC-US", "STAD-US", "UCEC-US") -->
                                                                              <!-- nn %>% -->
                                                                              <!--   filter(Cohort %in% cols.subset) %>%  -->
                                                                              <!--   ggplot(aes(x = value, colour = Class)) + -->
                                                                              <!--                     stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                              <!--                     theme_bw() + -->
                                                                              <!--                     scale_color_manual(values = cols) + -->
                                                                              <!--                     labs(x = "Nb of TFBSs assigned to protein-coding genes", y = "ECDF") + -->
                                                                              <!--                     theme(legend.position = "bottom") + -->
                                                                              <!--                     facet_grid(vars(Cohort), vars(variable), scales = "free") -->
                                                                              <!-- ``` -->
                                                                              
                                                                              
                                                                              
                                                                              <!-- ```{r Corr_PCG, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Number of target genes and TFBS associated with the predicted and non-predicted protein-coding genes (all cohorts).", fig.height = 10, fig.width = 8, fig.align = "center"} -->
                                                                              <!-- analyzed.genes.cohort.list <- sapply(genes.in.cohort.list, data.frame) -->
                                                                                <!-- names(analyzed.genes.cohort.list) <- names(genes.in.cohort.list) -->
                                                                                  <!-- analyzed.genes.cohort.list <- analyzed.genes.cohort.list[11:20] ## 1-10 entries are NULL -->
                                                                                  
                                                                                  
                                                                                  
                                                                                  <!-- cols.point <- c("Non-Predicted" = "#E69F00", "Predicted" = "#1b9e77") -->
                                                                                    
                                                                                    <!-- p2 <- ggplot(analyzed.pcg.summary, aes(x = Nb_target, y = Nb_TFBS, colour = Predicted)) + -->
                                                                                      <!--   geom_point() + -->
                                                                                      <!--   theme_bw() + -->
                                                                                      <!--   scale_color_manual(values = cols.point) + -->
                                                                                      <!--   theme(legend.position = "bottom")  -->
                                                                                      <!-- ggMarginal(p2, type = "boxplot", groupColour = TRUE, groupFill = TRUE) -->
                                                                                      <!-- ``` -->
                                                                                      
                                                                                      
                                                                                      
                                                                                      <!-- ```{r Nb_TFBS_per_gene_miRNA, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="", fig.height = 10, fig.width = 8, fig.align = "center"} -->
                                                                                      
                                                                                      <!-- cohorts.all <- c(cohorts, "BASIS_All", "BASIS_ER+", "BASIS_ER-") -->
                                                                                        <!-- selected.miRNA.list  <- vector(mode = "list", length = length(cohorts.all)) -->
                                                                                          <!-- miRNA.in.cohort.list <- vector(mode = "list", length = length(cohorts.all)) -->
                                                                                            <!-- Nb.mut.samples.mirna <- list() -->
                                                                                              <!-- for (co in cohorts.all) {  -->
                                                                                                  <!--    if (!grepl(co, pattern = "BASIS")) {  -->
                                                                                                      <!--      co <- paste0(co, "-US")  -->
                                                                                                        <!--      load(paste0("xseq_results/miRNA/xseq_trans_premiRNA_all_mutation_types_all_effects_", co, "_mirna_network_targetScanClean.Rdata")) ## xseq.pred.all.concat  -->
                                                                                                      <!--    } else {  -->
                                                                                                          <!--      load(paste0("xseq_results/miRNA/", co, "_miRNA.Rdata"))  -->
                                                                                                          <!--    }  -->
                                                                                                  
                                                                                                  <!--   ## Parse xseq output -->
                                                                                                  <!--   xseq.pred.mirna <- xseq.pred.all.concat %>%  -->
                                                                                                    <!--                       filter(mut_effect == "TFBS") %>%  -->
                                                                                                    <!--                       select(sample, gene, prob_mut, prob_gene) %>%  -->
                                                                                                    <!--                       distinct() %>%  -->
                                                                                                    <!--                       arrange(desc(prob_gene))  -->
                                                                                                    
                                                                                                    
                                                                                                    
                                                                                                    <!--       Nb.mut.samples.mirna[[co]] <- xseq.pred.mirna %>%  -->
                                                                                                      <!--                                 group_by(gene) %>%  -->
                                                                                                      <!--                                 summarise(Nb_samples = n()) %>%  -->
                                                                                                      <!--                                 arrange(desc(Nb_samples)) %>%  -->
                                                                                                      <!--                                 rename(Gene = gene) %>%  -->
                                                                                                      <!--                                 mutate(Cohort = co) -->
                                                                                                      
                                                                                                      <!--   if (grepl(co, pattern = "BASIS")) { -->
                                                                                                          <!--       Nb.mut.samples.mirna[[co]] <- merge(Nb.mut.samples.mirna[[co]], mirna.premirna.dict, -->
                                                                                                                                                           <!--                                           by.x = "Gene", -->
                                                                                                                                                           <!--                                           by.y = "miR", -->
                                                                                                                                                           <!--                                           allow.cartesian = TRUE) %>%  -->
                                                                                                            <!--                                     select(mature_ID, Nb_samples, Cohort) %>%  -->
                                                                                                            <!--                                     rename(Gene = mature_ID) %>%  -->
                                                                                                            <!--                                     tibble() -->
                                                                                                            
                                                                                                            <!--   } -->
                                                                                                      
                                                                                                      
                                                                                                      <!--   miRNA.in.cohort.list[[co]] <- unique(xseq.pred.mirna$gene) -->
                                                                                                        <!--   selected.miRNA.list[[co]]  <- unique(subset(xseq.pred.mirna, prob_gene >= xseq.th.mirna[[co]])$gene) -->
                                                                                                          <!-- } -->
                                                                                              <!-- selected.mirna.list  <- selected.miRNA.list -->
                                                                                                <!-- mirna.in.cohort.list <- miRNA.in.cohort.list -->
                                                                                                  
                                                                                                  <!-- miRNA.names.list <- list(Predicted = unique(as.vector(unlist(selected.miRNA.list))), -->
                                                                                                                                  <!--                          Analyzed  = unique(as.vector(unlist(miRNA.in.cohort.list)))) -->
                                                                                                    
                                                                                                    <!-- reported.cancer.miRNA <- xseq.sup.mat.tab$mirna_cancer_genes %>%  -->
                                                                                                      <!--                           mutate(ID = paste(miRNA, premiRNA, sep = "::")) %>%  -->
                                                                                                      <!--                           select(ID) %>%  -->
                                                                                                      <!--                           distinct() %>%  -->
                                                                                                      <!--                           unlist() %>%  -->
                                                                                                      <!--                           as.vector() %>%  -->
                                                                                                      <!--                           tolower() -->
                                                                                                      
                                                                                                      <!-- miRNA.genes           <- xseq.sup.mat.tab$mirna_history_mirbase %>%  -->
                                                                                                        <!--                           mutate(ID = paste(miRNA, premiRNA, sep = "::")) %>%  -->
                                                                                                        <!--                           select(ID) %>%  -->
                                                                                                        <!--                           distinct() %>%  -->
                                                                                                        <!--                           unlist() %>%  -->
                                                                                                        <!--                           as.vector() %>%  -->
                                                                                                        <!--                           tolower() -->
                                                                                                        
                                                                                                        <!-- ## The TFBSs are associated to premiRNA, not to mature miRNAs, convert the premiRNA  -->
                                                                                                        <!-- ## names to mature miRNA before matching the miRNA analyzed by xseq -->
                                                                                                        <!-- for (l in 1:length(miRNA.names.list)) { -->
                                                                                                            
                                                                                                            <!--   l.tmp <- miRNA.names.list[[l]][which(!grepl(miRNA.names.list[[l]], pattern = "::"))] -->
                                                                                                              <!--   l.tmp <- data.frame(miRNA_names = tolower(l.tmp)) -->
                                                                                                                <!--   l.tmp <- merge(l.tmp, mirna.premirna.dict, -->
                                                                                                                                        <!--                  by.x = "miRNA_names", -->
                                                                                                                                        <!--                  by.y = "miR", -->
                                                                                                                                        <!--                  allow.cartesian = TRUE) -->
                                                                                                                  
                                                                                                                  <!--   miRNA.names.list[[l]] <- c(miRNA.names.list[[l]], l.tmp$mature_ID) -->
                                                                                                                    <!-- } -->
                                                                                                        <!-- predicted.genes.pancancer <- miRNA.names.list[[1]][grepl(miRNA.names.list[[1]],pattern = "::")] -->
                                                                                                          <!-- xseq.analyzed.genes       <- miRNA.names.list[[2]][grepl(miRNA.names.list[[2]],pattern = "::")] -->
                                                                                                            <!-- all.associated.features.tab.mir <- all.associated.features.tab %>%  -->
                                                                                                              <!--                                     filter(grepl(Associated_gene, pattern = "^hsa-")) -->
                                                                                                              <!-- all.associated.features.tab.mir <- merge(all.associated.features.tab.mir, mirna.premirna.dict, -->
                                                                                                                                                              <!--                                       by.x = "Associated_gene", -->
                                                                                                                                                              <!--                                       by.y = "premiR", -->
                                                                                                                                                              <!--                                       allow.cartesian = TRUE) %>%  -->
                                                                                                                <!--                                     within(rm(Associated_gene, miR)) %>%  -->
                                                                                                                <!--                                     rename(Associated_gene = mature_ID) -->
                                                                                                                <!-- ## A subset with only the TFBS associated to the genes analyzed by xseq across the cohorts -->
                                                                                                                <!-- # TFBS.xseq.analyzed.miRNA <- all.associated.features.tab %>%  -->
                                                                                                                <!-- #                               filter(Associated_gene %in% xseq.analyzed.genes) %>%  -->
                                                                                                                <!-- #                               distinct() -->
                                                                                                                <!-- Nb.TFBSs.per.miRNA <- all.associated.features.tab.mir %>%  -->
                                                                                                                  <!--                         distinct() %>%  -->
                                                                                                                  <!--                         group_by(Associated_gene) %>%  -->
                                                                                                                  <!--                         summarise(Nb_TFBS = n(), .groups = "drop") %>%  -->
                                                                                                                  <!--                         arrange(desc(Nb_TFBS)) -->
                                                                                                                  <!-- Nb.TFBSs.per.miRNA$Predicted     <- ifelse(Nb.TFBSs.per.miRNA$Associated_gene %in% predicted.genes.pancancer, yes = "Predicted", no = "") -->
                                                                                                                    <!-- Nb.TFBSs.per.miRNA$Non_Predicted <- ifelse(!Nb.TFBSs.per.miRNA$Associated_gene %in% predicted.genes.pancancer, yes = "Non-Predicted", no = "") -->
                                                                                                                      <!-- Nb.TFBSs.per.miRNA$Cancer <- ifelse(Nb.TFBSs.per.miRNA$Associated_gene %in% reported.cancer.miRNA, yes = "Cancer miRNA", no = "") -->
                                                                                                                        <!-- # Nb.TFBSs.per.miRNA$miRNA <- ifelse(Nb.TFBSs.per.miRNA$Associated_gene %in% miRNA.genes, yes = "miRNA", no = "") -->
                                                                                                                        
                                                                                                                        
                                                                                                                        <!-- Nb.TFBSs.per.miRNA.classes <- reshape2::melt(Nb.TFBSs.per.miRNA, -->
                                                                                                                                                                            <!--                                 id.vars = c("Nb_TFBS", "Associated_gene"), -->
                                                                                                                                                                            <!--                                 measure.vars = c("Predicted", "Non_Predicted", "Cancer")) %>%  -->
                                                                                                                          <!--                                 filter(value != "") %>%  -->
                                                                                                                          <!--                                 within(rm(variable)) %>%  -->
                                                                                                                          <!--                                 rename(Class = value) -->
                                                                                                                          <!-- cols <- c("Cancer miRNA" = "#E69F00", "Non-Predicted" = "#56B4E9", "Predicted" = "#1b9e77") -->
                                                                                                                            <!-- Nb.TFBSs.per.miRNA.classes$Class <- factor(Nb.TFBSs.per.miRNA.classes$Class, levels = c("Cancer miRNA", "Non-Predicted", "Predicted")) -->
                                                                                                                              <!-- ## Plot 1: ECDF -->
                                                                                                                              <!-- TFBS.distr.ecdf.plot <- ggplot(Nb.TFBSs.per.miRNA.classes, aes(x = Nb_TFBS, colour = Class)) + -->
                                                                                                                                <!--         stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                                                                                <!--         theme_bw() + -->
                                                                                                                                <!--         scale_color_manual(values = cols) +  -->
                                                                                                                                <!--         labs(x = "Nb of TFBSs assigned to miRNA genes", y = "ECDF") + -->
                                                                                                                                <!--         theme(legend.position = "none")  -->
                                                                                                                                <!-- ## Plot 2: histogram with real values -->
                                                                                                                                <!-- TFBS.distr.hist.real.plot <- ggplot(Nb.TFBSs.per.miRNA.classes, aes(x = Nb_TFBS, fill = Class)) + -->
                                                                                                                                  <!--         geom_histogram(binwidth = 50, position = "identity") + -->
                                                                                                                                  <!--         scale_fill_manual(values = cols) + -->
                                                                                                                                  <!--         theme_bw() + -->
                                                                                                                                  <!--         labs(x = "Nb of TFBSs assigned to miRNA genes", y = "Frequency") + -->
                                                                                                                                  <!--         facet_grid(Class ~ .) + -->
                                                                                                                                  <!--         theme(legend.position = "right") -->
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  <!-- ############################## -->
                                                                                                                                  <!-- ## Rename BASIS miRNA names ## -->
                                                                                                                                  <!-- ############################## -->
                                                                                                                                  <!-- # selected.mirna.list -->
                                                                                                                                  <!-- # mirna.in.cohort.list -->
                                                                                                                                  <!-- for (l in 1:length(selected.mirna.list)) { -->
                                                                                                                                      
                                                                                                                                      <!--   l.tmp <- selected.mirna.list[[l]][which(!grepl(selected.mirna.list[[l]], pattern = "::"))] -->
                                                                                                                                        <!--   l.tmp <- data.frame(miRNA_names = tolower(l.tmp)) -->
                                                                                                                                          <!--   l.tmp <- merge(l.tmp, mirna.premirna.dict, -->
                                                                                                                                                                  <!--                  by.x = "miRNA_names", -->
                                                                                                                                                                  <!--                  by.y = "miR", -->
                                                                                                                                                                  <!--                  allow.cartesian = TRUE) -->
                                                                                                                                            
                                                                                                                                            <!--   selected.mirna.list[[l]] <- c(selected.mirna.list[[l]], l.tmp$mature_ID) -->
                                                                                                                                              <!-- } -->
                                                                                                                                  
                                                                                                                                  <!-- for (l in 1:length(mirna.in.cohort.list)) { -->
                                                                                                                                      
                                                                                                                                      <!--   l.tmp <- mirna.in.cohort.list[[l]][which(!grepl(mirna.in.cohort.list[[l]], pattern = "::"))] -->
                                                                                                                                        <!--   l.tmp <- data.frame(miRNA_names = tolower(l.tmp)) -->
                                                                                                                                          <!--   l.tmp <- merge(l.tmp, mirna.premirna.dict, -->
                                                                                                                                                                  <!--                  by.x = "miRNA_names", -->
                                                                                                                                                                  <!--                  by.y = "miR", -->
                                                                                                                                                                  <!--                  allow.cartesian = TRUE) -->
                                                                                                                                            
                                                                                                                                            <!--   mirna.in.cohort.list[[l]] <- c(mirna.in.cohort.list[[l]], l.tmp$mature_ID) -->
                                                                                                                                              <!-- } -->
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  <!-- selected.mirna.list  <- sapply(selected.mirna.list, function(x){ -->
                                                                                                                                      <!--                           x[grepl(x, pattern = "::")] -->
                                                                                                                                      <!--                         }) -->
                                                                                                                                    <!-- mirna.in.cohort.list <- sapply(mirna.in.cohort.list, function(x){ -->
                                                                                                                                        <!--                           x[grepl(x, pattern = "::")] -->
                                                                                                                                        <!--                         }) -->
                                                                                                                                      <!-- ``` -->
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      <!-- ```{r Nb_connections_per_miRNA, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Distribution of number of associated TFBSs (A and B) and number of targets in the regulatory network (C and D) among predicted, non-predicted, and cancer-related miRNA genes (aggregated from all evaluated cohorts) represented as Cumulative Distribution Function and Histogram of frequencies. The maximum number of targets of each miRNA in the regulatory network was restricted to 100.", fig.height = 8, fig.width = 9, fig.align = "center"} -->
                                                                                                                                      <!-- # predicted.genes.pancancer -->
                                                                                                                                      <!-- # xseq.analyzed.genes -->
                                                                                                                                      <!-- # reported.cancer.miRNA -->
                                                                                                                                      <!-- # miRNA.genes -->
                                                                                                                                      <!-- #  -->
                                                                                                                                      <!-- # xseq.net.mirna -->
                                                                                                                                      <!-- xseq.net.mirna <- merge(xseq.net.mirna, mirna.premirna.dict, -->
                                                                                                                                                                     <!--                          by.x = "Gene", -->
                                                                                                                                                                     <!--                          by.y = "miR", -->
                                                                                                                                                                     <!--                          allow.cartesian = TRUE) -->
                                                                                                                                        <!-- ## Genes in the network that were analyzed by xseq -->
                                                                                                                                        <!-- xseq.net.mirna.analyzed.genes <- xseq.net.mirna %>%  -->
                                                                                                                                          <!--                                 filter(mature_ID %in% xseq.analyzed.genes) %>%  -->
                                                                                                                                          <!--                                 group_by(mature_ID) %>%  -->
                                                                                                                                          <!--                                 summarise(Nb_target = n(), .groups = "drop") %>%  -->
                                                                                                                                          <!--                                 arrange(desc(Nb_target)) %>%  -->
                                                                                                                                          <!--                                 rename(Gene = mature_ID) -->
                                                                                                                                          <!-- xseq.net.mirna.analyzed.genes$Predicted     <- ifelse(xseq.net.mirna.analyzed.genes$Gene %in% predicted.genes.pancancer, yes = "Predicted", no = "") -->
                                                                                                                                            <!-- xseq.net.mirna.analyzed.genes$Non_Predicted <- ifelse(!xseq.net.mirna.analyzed.genes$Gene %in% predicted.genes.pancancer, yes = "Non-Predicted", no = "") -->
                                                                                                                                              <!-- xseq.net.mirna.analyzed.genes$Cancer        <- ifelse(xseq.net.mirna.analyzed.genes$Gene %in% reported.cancer.miRNA, yes = "Cancer miRNA", no = "") -->
                                                                                                                                                <!-- xseq.net.mirna.analyzed.genes$miRNA            <- ifelse(xseq.net.mirna.analyzed.genes$Gene %in% miRNA.genes, yes = "miRNA", no = "") -->
                                                                                                                                                  <!-- xseq.net.mirna.analyzed.genes.classes <- reshape2::melt(xseq.net.mirna.analyzed.genes, -->
                                                                                                                                                                                                                 <!--                                         id.vars = c("Gene", "Nb_target"), -->
                                                                                                                                                                                                                 <!--                                         measure.vars = c("Predicted", "Non_Predicted", "Cancer", "miRNA")) %>%  -->
                                                                                                                                                    <!--                                         filter(value != "") %>%  -->
                                                                                                                                                    <!--                                         within(rm(variable)) %>%  -->
                                                                                                                                                    <!--                                         rename(Class = value) -->
                                                                                                                                                    <!-- ## Plot 1: ECDF -->
                                                                                                                                                    <!-- conn.distr.ecdf.plot <- ggplot(xseq.net.mirna.analyzed.genes.classes, aes(x = Nb_target, colour = Class)) + -->
                                                                                                                                                      <!--         stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                                                                                                      <!--         theme_bw() + -->
                                                                                                                                                      <!--         scale_color_manual(values = cols) +  -->
                                                                                                                                                      <!--         labs(x = "Nb of targets - miRNA genes", y = "ECDF") + -->
                                                                                                                                                      <!--         theme(legend.position = "none")  -->
                                                                                                                                                      <!-- ## Plot 2: histogram with real values -->
                                                                                                                                                      <!-- conn.distr.hist.real.plot <- ggplot(xseq.net.mirna.analyzed.genes.classes, aes(x = Nb_target, fill = Class)) + -->
                                                                                                                                                        <!--         geom_histogram(binwidth = 5, position = "identity") + -->
                                                                                                                                                        <!--         scale_fill_manual(values = cols) + -->
                                                                                                                                                        <!--         theme_bw() + -->
                                                                                                                                                        <!--         labs(x = "Nb of TFBSs assigned to miRNA genes", y = "Frequency") + -->
                                                                                                                                                        <!--         facet_grid(Class ~ .) + -->
                                                                                                                                                        <!--         theme(legend.position = "right")  -->
                                                                                                                                                        <!-- enriched.terms.gg <- plot_grid(TFBS.distr.ecdf.plot, -->
                                                                                                                                                                                              <!--                                TFBS.distr.hist.real.plot, -->
                                                                                                                                                                                              <!--                                conn.distr.ecdf.plot, -->
                                                                                                                                                                                              <!--                                conn.distr.hist.real.plot, -->
                                                                                                                                                                                              <!--                                ncol = 2, align = "hv", labels = c('A)', 'B)', 'C)', 'D)'), label_size = 13, -->
                                                                                                                                                                                              <!--                                label_colour = "black", hjust = -0.25, vjust = 1.15, -->
                                                                                                                                                                                              <!--                                rel_heights = c(1, 1, 1, 1), rel_widths = c(1.5, 1, 1.5, 1), -->
                                                                                                                                                                                              <!--                                axis = "right") -->
                                                                                                                                                          <!-- ## Combine plots -->
                                                                                                                                                          <!-- enriched.terms.gg -->
                                                                                                                                                          <!-- ``` -->
                                                                                                                                                          
                                                                                                                                                          
                                                                                                                                                          <!-- ```{r export_Rdata_selected_genes, cache=FALSE, include=FALSE, echo=FALSE, eval=TRUE} -->
                                                                                                                                                          
                                                                                                                                                          <!-- all.predicted.genes <- list(Predicted_PCG   = selected.pcg.list, -->
                                                                                                                                                                                             <!--                             Analyzed_PCG    = pcg.in.cohort.list, -->
                                                                                                                                                                                             <!--                             Predicted_miRNA = selected.mirna.list, -->
                                                                                                                                                                                             <!--                             Analyzed_miRNA  = mirna.in.cohort.list) -->
                                                                                                                                                            
                                                                                                                                                            <!-- save(all.predicted.genes, file = "Predicted_genes_all_cohorts.RData") -->
                                                                                                                                                            
                                                                                                                                                            <!-- ``` -->
                                                                                                                                                            
                                                                                                                                                            
                                                                                                                                                            
                                                                                                                                                            <!-- ```{r eCDF_plot1_miRNA, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Empirical Cumulative Distribution Functions of number of TFBS, number of targets and number of mutations in miRNA genes, separated by cohorts (part 1).", fig.height = 10, fig.width = 8, fig.align = "center"} -->
                                                                                                                                                            
                                                                                                                                                            <!-- analyzed.mirnas.cohort.list <- sapply(miRNA.in.cohort.list, data.frame) -->
                                                                                                                                                              <!-- names(analyzed.mirnas.cohort.list) <- names(miRNA.in.cohort.list) -->
                                                                                                                                                                <!-- analyzed.mirnas.cohort.list <- analyzed.mirnas.cohort.list[11:20] ## 1-10 entries are NULL -->
                                                                                                                                                                <!-- Nb.targets.per.mirna <- xseq.net.mirna.analyzed.genes -->
                                                                                                                                                                  <!-- analyzed.mirna.summary <- Nb.targets.per.mirna %>%  -->
                                                                                                                                                                    <!--                           select(Gene, Nb_target, Predicted) %>%  -->
                                                                                                                                                                    <!--                           rename(Associated_gene = Gene) %>%  -->
                                                                                                                                                                    <!--                           left_join(Nb.TFBSs.per.miRNA, by = "Associated_gene") %>%  -->
                                                                                                                                                                    <!--                           rename(Predicted = Predicted.x, -->
                                                                                                                                                                                                            <!--                                  Gene      = Associated_gene) %>%  -->
                                                                                                                                                                    <!--                           select(Gene, Nb_target, Nb_TFBS, Predicted, Cancer) %>%  -->
                                                                                                                                                                    <!--                           mutate(Predicted = ifelse(Predicted == "", yes = "Non-Predicted", no = Predicted)) -->
                                                                                                                                                                    
                                                                                                                                                                    <!-- analyzed.mirna.cohorts.summary <- lapply(Nb.mut.samples.mirna, function(l){ -->
                                                                                                                                                                        <!--                                   merge(l, analyzed.mirna.summary, by = "Gene") %>%  -->
                                                                                                                                                                        <!--                                   reshape2::melt(id.vars = c("Gene", "Nb_target", "Nb_TFBS", "Nb_samples", "Cohort"), -->
                                                                                                                                                                                                                                <!--                                                  measure.vars = c("Predicted", "Cancer")) %>%  -->
                                                                                                                                                                        <!--                                                  filter(value != "") %>%  -->
                                                                                                                                                                        <!--                                                  within(rm(variable)) %>%  -->
                                                                                                                                                                        <!--                                                  rename(Class = value) -->
                                                                                                                                                                        <!--                                 }) -->
                                                                                                                                                                      <!-- cols <- c("Cancer miRNA" = "#E69F00", "Non-Predicted" = "#56B4E9", "Predicted" = "#1b9e77") -->
                                                                                                                                                                        <!-- nn <- -->
                                                                                                                                                                          <!-- reshape2::melt(rbindlist(analyzed.mirna.cohorts.summary), -->
                                                                                                                                                                                                <!--                id.vars = c("Gene", "Cohort", "Class"), -->
                                                                                                                                                                                                <!--                measure.vars = c("Nb_target", "Nb_TFBS", "Nb_samples")) %>% -->
                                                                                                                                                                          <!--                filter(value != "") -->
                                                                                                                                                                          <!-- cols.subset <- c("BASIS_All", "BASIS_ER+", "BASIS_ER-", "BRCA-US", "HNSC-US") -->
                                                                                                                                                                            <!-- nn %>% -->
                                                                                                                                                                            <!--   filter(Cohort %in% cols.subset) %>%  -->
                                                                                                                                                                            <!--   ggplot(aes(x = value, colour = Class)) + -->
                                                                                                                                                                            <!--                     stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                                                                                                                            <!--                     theme_bw() + -->
                                                                                                                                                                            <!--                     scale_color_manual(values = cols) + -->
                                                                                                                                                                            <!--                     labs(x = "Nb of TFBSs assigned to protein-coding genes", y = "ECDF") + -->
                                                                                                                                                                            <!--                     theme(legend.position = "bottom") + -->
                                                                                                                                                                            <!--                     facet_grid(vars(Cohort), vars(variable), scales = "free") -->
                                                                                                                                                                            <!-- ``` -->
                                                                                                                                                                            
                                                                                                                                                                            
                                                                                                                                                                            
                                                                                                                                                                            <!-- ```{r eCDF_plot2_miRNA, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Empirical Cumulative Distribution Functions of number of TFBS, number of targets and number of mutations in miRNA genes, separated by cohorts (part 1).", fig.height = 10, fig.width = 8, fig.align = "center"} -->
                                                                                                                                                                            
                                                                                                                                                                            <!-- analyzed.mirna.cohorts.summary <- lapply(Nb.mut.samples.mirna, function(l){ -->
                                                                                                                                                                                <!--                                   merge(l, analyzed.mirna.summary, by = "Gene") %>%  -->
                                                                                                                                                                                <!--                                   reshape2::melt(id.vars = c("Gene", "Nb_target", "Nb_TFBS", "Nb_samples", "Cohort"), -->
                                                                                                                                                                                                                                        <!--                                                  measure.vars = c("Predicted", "Cancer")) %>%  -->
                                                                                                                                                                                <!--                                                  filter(value != "") %>%  -->
                                                                                                                                                                                <!--                                                  within(rm(variable)) %>%  -->
                                                                                                                                                                                <!--                                                  rename(Class = value) -->
                                                                                                                                                                                <!--                                 }) -->
                                                                                                                                                                              <!-- cols <- c("Cancer miRNA" = "#E69F00", "Non-Predicted" = "#56B4E9", "Predicted" = "#1b9e77") -->
                                                                                                                                                                                <!-- nn <- -->
                                                                                                                                                                                  <!-- reshape2::melt(rbindlist(analyzed.mirna.cohorts.summary), -->
                                                                                                                                                                                                        <!--                id.vars = c("Gene", "Cohort", "Class"), -->
                                                                                                                                                                                                        <!--                measure.vars = c("Nb_target", "Nb_TFBS", "Nb_samples")) %>% -->
                                                                                                                                                                                  <!--                filter(value != "") -->
                                                                                                                                                                                  <!-- cols.subset <- c("LIHC-US", "LUAD-US", "LUSC-US", "STAD-US", "UCEC-US") -->
                                                                                                                                                                                    <!-- nn %>% -->
                                                                                                                                                                                    <!--   filter(Cohort %in% cols.subset) %>%  -->
                                                                                                                                                                                    <!--   ggplot(aes(x = value, colour = Class)) + -->
                                                                                                                                                                                    <!--                     stat_ecdf(size = 1.15, pad = FALSE) + -->
                                                                                                                                                                                    <!--                     theme_bw() + -->
                                                                                                                                                                                    <!--                     scale_color_manual(values = cols) + -->
                                                                                                                                                                                    <!--                     labs(x = "Nb of TFBSs assigned to protein-coding genes", y = "ECDF") + -->
                                                                                                                                                                                    <!--                     theme(legend.position = "bottom") + -->
                                                                                                                                                                                    <!--                     facet_grid(vars(Cohort), vars(variable), scales = "free") -->
                                                                                                                                                                                    <!-- ``` -->
                                                                                                                                                                                    
                                                                                                                                                                                    
                                                                                                                                                                                    
                                                                                                                                                                                    <!-- ```{r Corr_miRNA, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Number of target genes and TFBS associated with the predicted and non-predicted miRNA genes (all cohorts).", fig.height = 10, fig.width = 8, fig.align = "center", message=FALSE, warning=FALSE} -->
                                                                                                                                                                                    
                                                                                                                                                                                    <!-- cols.point <- c("Non-Predicted" = "#E69F00", "Predicted" = "#1b9e77") -->
                                                                                                                                                                                      <!-- p2 <- ggplot(analyzed.mirna.summary, aes(x = Nb_target, y = Nb_TFBS, colour = Predicted)) + -->
                                                                                                                                                                                        <!--   geom_point() + -->
                                                                                                                                                                                        <!--   theme_bw() + -->
                                                                                                                                                                                        <!--   scale_color_manual(values = cols.point) + -->
                                                                                                                                                                                        <!--   theme(legend.position = "bottom")  -->
                                                                                                                                                                                        <!-- ggMarginal(p2, type = "boxplot", groupColour = TRUE, groupFill = TRUE) -->
                                                                                                                                                                                        <!-- ``` -->
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        ```{r TFBS_mut_rdata, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, }
                                                                                                                                                                                      # library(dplyr)
                                                                                                                                                                                      # library(data.table)
                                                                                                                                                                                      # 
                                                                                                                                                                                      # 
                                                                                                                                                                                      # cohorts       <- c("BRCA-US", "HNSC-US", "LIHC-US", "LUAD-US", "LUSC-US", "STAD-US", "UCEC-US")
                                                                                                                                                                                      # TFBS.mut.list <- vector(mode = "list", length = length(cohorts))
                                                                                                                                                                                      # for (co in cohorts){
                                                                                                                                                                                      # 
                                                                                                                                                                                      #   mutations.file <- paste0("/storage/scratch/TCGA/ICGC/Variant_calling_download/", co, "/VCF_muse/VCF/post_process/ICGC_mutations_WXS_and_TFBS_annotated_", co, ".tab")
                                                                                                                                                                                      # 
                                                                                                                                                                                      # 
                                                                                                                                                                                      #   mutations <- fread(mutations.file)
                                                                                                                                                                                      # 
                                                                                                                                                                                      # 
                                                                                                                                                                                      #   ## Group the TFs matching the same mutation (chr,start, end)
                                                                                                                                                                                      #   TFBS.mut.list[[co]] <- mutations %>%
                                                                                                                                                                                      #                           filter(effect == "TFBS") %>%
                                                                                                                                                                                      #                           within(rm(variant_type, effect)) %>%
                                                                                                                                                                                      #                           mutate(ID = paste(chr, start, end, sep = "_")) %>% 
                                                                                                                                                                                      #                           # filter(sample == "SP89443") %>% 
                                                                                                                                                                                      #                           group_by(ID) %>%
                                                                                                                                                                                      #                           mutate(TFs = paste(unique(TF), collapse = ",")) %>%
                                                                                                                                                                                      #                           ungroup() %>%
                                                                                                                                                                                      #                           within(rm(TF)) %>%
                                                                                                                                                                                      #                           distinct() %>%
                                                                                                                                                                                      #                           select(sample, hgnc_symbol, TFs, Cohort)
                                                                                                                                                                                      # 
                                                                                                                                                                                      # }
                                                                                                                                                                                      # 
                                                                                                                                                                                      # save(TFBS.mut.list, file = "/storage/mathelierarea/processed/jamondra/Projects/dysmir/TFBS_mutations_per_cohort.Rdata")
                                                                                                                                                                                      ```
                                                                                                                                                                                      
                                                                                                                                                                                      
                                                                                                                                                                                      <!-- \pagebreak -->
                                                                                                                                                                                        <!-- \begin{landscape} -->
                                                                                                                                                                                        <!-- ```{r FEA_BASIS_1, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Functional enrichment analysis considering dysregulated genes in the networks of predicted miRNA and protein-coding genes in the BASIS cohort. The barplots represent the significance (-log$_{10}$(P-value)) of the top-10 most enriched terms (rows) associated to the dysregulated genes when considering protein-coding genes with LoF mutations (PCGs::LoF; A-C), protein-coding genes with cis-regulatory mutations (PCG::TFBS; B-D), and miRNAs with cis-regulatory mutations (miRNAs::TFBS; E). Enriched terms from KEGG 2021 Human (A-B) and WikiPathways (C-E) are provided.", fig.height = 100, fig.width = 180, fig.align = "center"} -->
                                                                                                                                                                                        <!-- load("FEA_heatmaps/Enriched_terms_plots_BASIS.RData")  ## GSEA.enrich.plots -->
                                                                                                                                                                                      <!-- ########## -->
                                                                                                                                                                                        <!-- ## KEGG ## -->
                                                                                                                                                                                        <!-- ########## -->
                                                                                                                                                                                        <!-- ## KEGG TFBS PCG -->
                                                                                                                                                                                        <!-- gene.names <- GSEA.enrich.plots[["KEGG_2021_Human"]][["PCG.TFBS"]][["gene_names"]] -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.kegg <- GSEA.enrich.plots[["KEGG_2021_Human"]][["PCG.TFBS"]][["plot"]] -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.kegg.plot <- -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.kegg + -->
                                                                                                                                                                                        <!--   scale_x_discrete(limits = gene.names) + -->
                                                                                                                                                                                        <!--   scale_fill_distiller(palette = "Blues", direction = +1) + -->
                                                                                                                                                                                        <!--   labs(title = " PCG.TFBS - KEGG") + -->
                                                                                                                                                                                        <!--   theme(text = element_text(size = 90), -->
                                                                                                                                                                                                       <!--         axis.text.x = element_text(hjust=1, size = 100), -->
                                                                                                                                                                                                       <!--         axis.text.y = element_text(hjust=1, size = 100)) + -->
                                                                                                                                                                                        <!--   guides(fill = guide_colourbar(barwidth = 9, barheight = 30)) + -->
                                                                                                                                                                                        <!--         theme( -->
                                                                                                                                                                                                              <!--         panel.background = element_rect(fill = "white"), -->
                                                                                                                                                                                                              <!--         plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"), -->
                                                                                                                                                                                                              <!--         plot.background = element_rect( -->
                                                                                                                                                                                                                                                             <!--           fill = "#fefefe", -->
                                                                                                                                                                                                                                                             <!--           colour = "white", -->
                                                                                                                                                                                                                                                             <!--           size = 1)) -->
                                                                                                                                                                                        <!-- ## KEGG LoF -->
                                                                                                                                                                                        <!-- # gene.names <- GSEA.enrich.plots[["KEGG_2019_Human"]][["PCG.LoF"]][["gene_names"]] -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.kegg <- GSEA.enrich.plots[["KEGG_2019_Human"]][["PCG.LoF"]][["plot"]] -->
                                                                                                                                                                                        <!-- #  -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.kegg.plot <- -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.kegg + -->
                                                                                                                                                                                        <!-- #   scale_x_discrete(limits = gene.names) + -->
                                                                                                                                                                                        <!-- #   scale_fill_distiller(palette = "Blues", direction = +1) + -->
                                                                                                                                                                                        <!-- #   guides(fill = guide_colourbar(barwidth = 9, barheight = 30)) + -->
                                                                                                                                                                                        <!-- #   labs(title = " PCG.LoF - KEGG") + -->
                                                                                                                                                                                        <!-- #   theme(text = element_text(size = 100), -->
                                                                                                                                                                                        <!-- #           axis.text.x = element_text(hjust=1, size = 100), -->
                                                                                                                                                                                        <!-- #           axis.text.y = element_text(hjust=1, size = 100)) + -->
                                                                                                                                                                                        <!-- #         theme( -->
                                                                                                                                                                                        <!-- #         panel.background = element_rect(fill = "white"), -->
                                                                                                                                                                                        <!-- #         plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"), -->
                                                                                                                                                                                        <!-- #         plot.background = element_rect( -->
                                                                                                                                                                                        <!-- #           fill = "#fefefe", -->
                                                                                                                                                                                        <!-- #           colour = "white", -->
                                                                                                                                                                                        <!-- #           size = 1)) -->
                                                                                                                                                                                        <!-- ################## -->
                                                                                                                                                                                        <!-- ## WikiPathways ## -->
                                                                                                                                                                                        <!-- ################## -->
                                                                                                                                                                                        <!-- ## WikiPathways TFBS PCG -->
                                                                                                                                                                                        <!-- gene.names <- GSEA.enrich.plots[["WikiPathways_2019_Human"]][["PCG.TFBS"]][["gene_names"]] -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.wiki <- GSEA.enrich.plots[["WikiPathways_2019_Human"]][["PCG.TFBS"]][["plot"]] -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.wiki.plot <- -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.wiki + -->
                                                                                                                                                                                        <!--   scale_x_discrete(limits = gene.names) + -->
                                                                                                                                                                                        <!-- scale_fill_distiller(palette = "Blues", direction = +1) + -->
                                                                                                                                                                                        <!-- labs(title = "PCG.TFBS - WikiPathways") + -->
                                                                                                                                                                                        <!-- guides(fill = guide_colourbar(barwidth = 9, barheight = 30)) + -->
                                                                                                                                                                                        <!-- theme(text = element_text(size = 100), -->
                                                                                                                                                                                                     <!--           axis.text.x = element_text(hjust=1, size = 100), -->
                                                                                                                                                                                                     <!--           axis.text.y = element_text(hjust=1, size = 100)) + -->
                                                                                                                                                                                        <!--         theme( -->
                                                                                                                                                                                                              <!--         panel.background = element_rect(fill = "white"), -->
                                                                                                                                                                                                              <!--         plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"), -->
                                                                                                                                                                                                              <!--         plot.background = element_rect( -->
                                                                                                                                                                                                                                                             <!--           fill = "#fefefe", -->
                                                                                                                                                                                                                                                             <!--           colour = "white", -->
                                                                                                                                                                                                                                                             <!--           size = 1)) -->
                                                                                                                                                                                        <!-- ## WikiPathways LoF -->
                                                                                                                                                                                        <!-- # gene.names <- GSEA.enrich.plots[["WikiPathways_2019_Human"]][["PCG.LoF"]][["gene_names"]] -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.wiki <- GSEA.enrich.plots[["WikiPathways_2019_Human"]][["PCG.LoF"]][["plot"]] -->
                                                                                                                                                                                        <!-- #  -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.wiki.plot <- -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.wiki + -->
                                                                                                                                                                                        <!-- #   scale_x_discrete(limits = gene.names) + -->
                                                                                                                                                                                        <!-- # scale_fill_distiller(palette = "Blues", direction = +1) + -->
                                                                                                                                                                                        <!-- # guides(fill = guide_colourbar(barwidth = 9, barheight = 30)) + -->
                                                                                                                                                                                        <!-- # labs(title = "PCG.LoF - WikiPathways") + -->
                                                                                                                                                                                        <!-- # theme(text = element_text(size = 100), -->
                                                                                                                                                                                        <!-- #           axis.text.x = element_text(hjust=1, size = 100), -->
                                                                                                                                                                                        <!-- #           axis.text.y = element_text(hjust=1, size = 100)) + -->
                                                                                                                                                                                        <!-- #         theme( -->
                                                                                                                                                                                        <!-- #         panel.background = element_rect(fill = "white"), -->
                                                                                                                                                                                        <!-- #         plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"), -->
                                                                                                                                                                                        <!-- #         plot.background = element_rect( -->
                                                                                                                                                                                        <!-- #           fill = "#fefefe", -->
                                                                                                                                                                                        <!-- #           colour = "white", -->
                                                                                                                                                                                        <!-- #           size = 1)) -->
                                                                                                                                                                                        <!-- ## WikiPathways TFBS miRNA -->
                                                                                                                                                                                        <!-- # gene.names <- GSEA.enrich.plots[["WikiPathways_2019_Human"]][["premiRNA.TFBS"]][["gene_names"]] -->
                                                                                                                                                                                        <!-- # basis.premiRNA.tfbs.wiki <- GSEA.enrich.plots[["WikiPathways_2019_Human"]][["premiRNA.TFBS"]][["plot"]] -->
                                                                                                                                                                                        <!-- #  -->
                                                                                                                                                                                        <!-- # basis.premiRNA.tfbs.wiki.plot <- -->
                                                                                                                                                                                        <!-- # basis.premiRNA.tfbs.wiki + -->
                                                                                                                                                                                        <!-- #   scale_x_discrete(limits = gene.names) + -->
                                                                                                                                                                                        <!-- # scale_fill_distiller(palette = "Blues", direction = +1) + -->
                                                                                                                                                                                        <!-- # guides(fill = guide_colourbar(barwidth = 9, barheight = 30)) + -->
                                                                                                                                                                                        <!-- # labs(title = "miRNA.TFBS - WikiPathways") + -->
                                                                                                                                                                                        <!-- # theme(text = element_text(size = 100), -->
                                                                                                                                                                                        <!-- #           axis.text.x = element_text(hjust=1, size = 100), -->
                                                                                                                                                                                        <!-- #           axis.text.y = element_text(hjust=1, size = 100)) #+ -->
                                                                                                                                                                                        <!-- # #         theme( -->
                                                                                                                                                                                        <!-- # #         panel.background = element_rect(fill = "white"), -->
                                                                                                                                                                                        <!-- # #         plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"), -->
                                                                                                                                                                                        <!-- # #         plot.background = element_rect( -->
                                                                                                                                                                                        <!-- # #           fill = "#fefefe", -->
                                                                                                                                                                                        <!-- # #           colour = "white", -->
                                                                                                                                                                                        <!-- # #           size = 1)) -->
                                                                                                                                                                                        <!-- # enriched.terms.gg <- plot_grid(basis.pcg.tfbs.kegg.plot, -->
                                                                                                                                                                                        <!-- #                                basis.pcg.tfbs.wiki.plot, -->
                                                                                                                                                                                        <!-- #                                ncol = 2, align = "hv", labels = c('A)', 'B)'), label_size = 150, -->
                                                                                                                                                                                        <!-- #                                label_colour = "black", hjust = -2, vjust = 3, -->
                                                                                                                                                                                        <!-- #                                rel_heights = c(1, 1)) -->
                                                                                                                                                                                        <!-- #  -->
                                                                                                                                                                                        <!-- # ## Combine plots -->
                                                                                                                                                                                        <!-- # enriched.terms.gg -->
                                                                                                                                                                                        <!-- #  -->
                                                                                                                                                                                        <!-- # rm(GSEA.enrich.plots) -->
                                                                                                                                                                                        <!-- # rm(enriched.terms.gg) -->
                                                                                                                                                                                        <!-- ``` -->
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        <!-- ```{r FEA_BASIS_2, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Functional enrichment analysis considering dysregulated genes in the networks of predicted miRNA and protein-coding genes in the BASIS cohort. The heatmaps represent the significance (-log$_{10}$(P-value)) of the top-10 most enriched terms (rows) associated to the dysregulated genes when considering protein-coding genes with LoF mutations (PCGs::LoF; A-C), protein-coding genes with cis-regulatory mutations (PCG::TFBS; B-D), and miRNAs with cis-regulatory mutations (miRNAs::TFBS; E). Enriched terms from GO biological processes (GO-BP; A-B) and Panther 2016 (C-E) are provided.", fig.height = 100, fig.width = 200, fig.align = "center"} -->
                                                                                                                                                                                        <!-- load("FEA_heatmaps/Enriched_terms_plots_BASIS.RData")  ## GSEA.enrich.plots -->
                                                                                                                                                                                      <!-- ########################### -->
                                                                                                                                                                                        <!-- ## GO Biological process ## -->
                                                                                                                                                                                        <!-- ########################### -->
                                                                                                                                                                                        <!-- ## GO-BP TFBS PCG -->
                                                                                                                                                                                        <!-- gene.names <- GSEA.enrich.plots[["GO_Biological_Process_2021"]][["PCG.TFBS"]][["gene_names"]] -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.go.bp <- GSEA.enrich.plots[["GO_Biological_Process_2021"]][["PCG.TFBS"]][["plot"]] -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.go.bp.plot <- -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.go.bp + -->
                                                                                                                                                                                        <!--   scale_x_discrete(limits = gene.names) + -->
                                                                                                                                                                                        <!-- scale_fill_distiller(palette = "Blues", direction = +1) + -->
                                                                                                                                                                                        <!-- guides(fill = guide_colourbar(barwidth = 9, barheight = 30)) + -->
                                                                                                                                                                                        <!-- labs(title = "PCG.TFBS - GO-BP") + -->
                                                                                                                                                                                        <!-- theme(text = element_text(size = 110), -->
                                                                                                                                                                                                     <!--           axis.text.x = element_text(hjust=1, size = 110), -->
                                                                                                                                                                                                     <!--           axis.text.y = element_text(hjust=1, size = 110)) + -->
                                                                                                                                                                                        <!--         theme( -->
                                                                                                                                                                                                              <!--         panel.background = element_rect(fill = "white"), -->
                                                                                                                                                                                                              <!--         plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"), -->
                                                                                                                                                                                                              <!--         plot.background = element_rect( -->
                                                                                                                                                                                                                                                             <!--           fill = "#fefefe", -->
                                                                                                                                                                                                                                                             <!--           colour = "white", -->
                                                                                                                                                                                                                                                             <!--           size = 1)) -->
                                                                                                                                                                                        <!-- ## GO-BP LoF -->
                                                                                                                                                                                        <!-- # gene.names <- GSEA.enrich.plots[["GO_Biological_Process_2018"]][["PCG.LoF"]][["gene_names"]] -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.go.bp <- GSEA.enrich.plots[["GO_Biological_Process_2018"]][["PCG.LoF"]][["plot"]] -->
                                                                                                                                                                                        <!-- #  -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.go.bp.plot <- -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.go.bp + -->
                                                                                                                                                                                        <!-- #   scale_x_discrete(limits = gene.names) + -->
                                                                                                                                                                                        <!-- # scale_fill_distiller(palette = "Blues", direction = +1) + -->
                                                                                                                                                                                        <!-- # guides(fill = guide_colourbar(barwidth = 9, barheight = 30)) + -->
                                                                                                                                                                                        <!-- # labs(title = "PCG.LoF - GO-BP") + -->
                                                                                                                                                                                        <!-- # theme(text = element_text(size = 110), -->
                                                                                                                                                                                        <!-- #           axis.text.x = element_text(hjust=1, size = 110), -->
                                                                                                                                                                                        <!-- #           axis.text.y = element_text(hjust=1, size = 110)) + -->
                                                                                                                                                                                        <!-- #         theme( -->
                                                                                                                                                                                        <!-- #         panel.background = element_rect(fill = "white"), -->
                                                                                                                                                                                        <!-- #         plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"), -->
                                                                                                                                                                                        <!-- #         plot.background = element_rect( -->
                                                                                                                                                                                        <!-- #           fill = "#fefefe", -->
                                                                                                                                                                                        <!-- #           colour = "white", -->
                                                                                                                                                                                        <!-- #           size = 1)) -->
                                                                                                                                                                                        <!-- ################## -->
                                                                                                                                                                                        <!-- ## Panther 2016 ## -->
                                                                                                                                                                                        <!-- ################## -->
                                                                                                                                                                                        <!-- ## Panther TFBS PCG -->
                                                                                                                                                                                        <!-- gene.names <- GSEA.enrich.plots[["Panther_2016"]][["PCG.TFBS"]][["gene_names"]] -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.panther <- GSEA.enrich.plots[["Panther_2016"]][["PCG.TFBS"]][["plot"]] -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.panther.plot <- -->
                                                                                                                                                                                        <!-- basis.pcg.tfbs.panther + -->
                                                                                                                                                                                        <!--   scale_x_discrete(limits = gene.names) + -->
                                                                                                                                                                                        <!-- scale_fill_distiller(palette = "Blues", direction = +1) + -->
                                                                                                                                                                                        <!-- guides(fill = guide_colourbar(barwidth = 9, barheight = 30)) + -->
                                                                                                                                                                                        <!-- labs(title = "PCG.TFBS - Panther") + -->
                                                                                                                                                                                        <!-- theme(text = element_text(size = 110), -->
                                                                                                                                                                                                     <!--           axis.text.x = element_text(hjust=1, size = 110), -->
                                                                                                                                                                                                     <!--           axis.text.y = element_text(hjust=1, size = 110)) + -->
                                                                                                                                                                                        <!--         theme( -->
                                                                                                                                                                                                              <!--         panel.background = element_rect(fill = "white"), -->
                                                                                                                                                                                                              <!--         plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"), -->
                                                                                                                                                                                                              <!--         plot.background = element_rect( -->
                                                                                                                                                                                                                                                             <!--           fill = "#fefefe", -->
                                                                                                                                                                                                                                                             <!--           colour = "white", -->
                                                                                                                                                                                                                                                             <!--           size = 1)) -->
                                                                                                                                                                                        <!-- ## Panther LoF -->
                                                                                                                                                                                        <!-- # gene.names <- GSEA.enrich.plots[["Panther_2016"]][["PCG.LoF"]][["gene_names"]] -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.panther <- GSEA.enrich.plots[["Panther_2016"]][["PCG.LoF"]][["plot"]] -->
                                                                                                                                                                                        <!-- #  -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.panther.plot <- -->
                                                                                                                                                                                        <!-- # basis.pcg.lof.panther + -->
                                                                                                                                                                                        <!-- #   scale_x_discrete(limits = gene.names) + -->
                                                                                                                                                                                        <!-- # scale_fill_distiller(palette = "Blues", direction = +1) + -->
                                                                                                                                                                                        <!-- # guides(fill = guide_colourbar(barwidth = 9, barheight = 30)) + -->
                                                                                                                                                                                        <!-- # labs(title = "PCG.LoF - Panther") + -->
                                                                                                                                                                                        <!-- # theme(text = element_text(size = 110), -->
                                                                                                                                                                                        <!-- #           axis.text.x = element_text(hjust=1, size = 110), -->
                                                                                                                                                                                        <!-- #           axis.text.y = element_text(hjust=1, size = 110)) + -->
                                                                                                                                                                                        <!-- #         theme( -->
                                                                                                                                                                                        <!-- #         panel.background = element_rect(fill = "white"), -->
                                                                                                                                                                                        <!-- #         plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"), -->
                                                                                                                                                                                        <!-- #         plot.background = element_rect( -->
                                                                                                                                                                                        <!-- #           fill = "#fefefe", -->
                                                                                                                                                                                        <!-- #           colour = "white", -->
                                                                                                                                                                                        <!-- #           size = 1)) -->
                                                                                                                                                                                        <!-- ## Panther TFBS miRNA -->
                                                                                                                                                                                        <!-- # gene.names <- GSEA.enrich.plots[["Panther_2016"]][["premiRNA.TFBS"]][["gene_names"]] -->
                                                                                                                                                                                        <!-- # basis.premiRNA.tfbs.panther <- GSEA.enrich.plots[["Panther_2016"]][["premiRNA.TFBS"]][["plot"]] -->
                                                                                                                                                                                        <!-- #  -->
                                                                                                                                                                                        <!-- # basis.premiRNA.tfbs.panther.plot <- -->
                                                                                                                                                                                        <!-- # basis.premiRNA.tfbs.panther + -->
                                                                                                                                                                                        <!-- #   scale_x_discrete(limits = gene.names) + -->
                                                                                                                                                                                        <!-- # scale_fill_distiller(palette = "Blues", direction = +1) + -->
                                                                                                                                                                                        <!-- # guides(fill = guide_colourbar(barwidth = 9, barheight = 30)) + -->
                                                                                                                                                                                        <!-- # labs(title = "miRNA.TFBS- Panther") + -->
                                                                                                                                                                                        <!-- # theme(text = element_text(size = 110), -->
                                                                                                                                                                                        <!-- #           axis.text.x = element_text(hjust=1, size = 110), -->
                                                                                                                                                                                        <!-- #           axis.text.y = element_text(hjust=1, size = 110)) + -->
                                                                                                                                                                                        <!-- #         theme( -->
                                                                                                                                                                                        <!-- #         panel.background = element_rect(fill = "white"), -->
                                                                                                                                                                                        <!-- #         plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"), -->
                                                                                                                                                                                        <!-- #         plot.background = element_rect( -->
                                                                                                                                                                                        <!-- #           fill = "#fefefe", -->
                                                                                                                                                                                        <!-- #           colour = "white", -->
                                                                                                                                                                                        <!-- #           size = 1)) -->
                                                                                                                                                                                        <!-- enriched.terms.gg <- plot_grid(basis.pcg.tfbs.go.bp.plot, -->
                                                                                                                                                                                                                              <!--                                basis.pcg.tfbs.kegg.plot, -->
                                                                                                                                                                                                                              <!--                                # NULL,NULL, -->
                                                                                                                                                                                                                              <!--                                basis.pcg.tfbs.wiki.plot, -->
                                                                                                                                                                                                                              <!--                                basis.pcg.tfbs.panther.plot, -->
                                                                                                                                                                                                                              <!--                                ncol = 2, align = "hv", labels = c('A)', 'B)', "C)", "D)"), label_size = 150, -->
                                                                                                                                                                                                                              <!--                                label_colour = "black", hjust = -5, vjust = 1, -->
                                                                                                                                                                                                                              <!--                                rel_heights = c(1, 1, 1, 1)) -->
                                                                                                                                                                                        <!-- ## Combine plots -->
                                                                                                                                                                                        <!-- enriched.terms.gg -->
                                                                                                                                                                                        <!-- rm(GSEA.enrich.plots) -->
                                                                                                                                                                                        <!-- rm(enriched.terms.gg) -->
                                                                                                                                                                                        <!-- ``` -->
                                                                                                                                                                                        <!-- \end{landscape} -->
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        <!-- \pagebreak -->
                                                                                                                                                                                        <!-- ```{r Ratio_mutated_vs_selected_miRNA, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Ratio of dysregulated vs mutated samples in the miRNA analysis. Each panel represents a cohort with all the evaluated miRNAs ranked by the DAC probability. The color of each dot corresponds to the ratio between the number of samples with SSD >= 0.5 (i.e., with high association of dysregulation) divided by the number of mutated samples. The black line corresponds to the cohort-specific threshold defined by FDR of 0.05.", fig.height = 32, fig.width = 25, fig.align = "center", warning = FALSE} -->
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        <!-- cohorts.all <- c(cohorts, "BASIS_All", "BASIS_ER+", "BASIS_ER-") -->
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        <!-- rank.plot.list <- vector(mode = "list", length = length(cohorts)) -->
                                                                                                                                                                                        <!-- for (co in cohorts.all) { -->
                                                                                                                                                                                            
                                                                                                                                                                                            <!--   if (!grepl(co, pattern = "BASIS")) { -->
                                                                                                                                                                                                <!--     co <- paste0(co, "-US") -->
                                                                                                                                                                                                  
                                                                                                                                                                                                  <!--     load(paste0("xseq_results/miRNA/xseq_trans_premiRNA_all_mutation_types_all_effects_", co, "_mirna_network_targetScanClean.Rdata")) ## xseq.pred.all.concat -->
                                                                                                                                                                                                
                                                                                                                                                                                                <!--   } else { -->
                                                                                                                                                                                                    
                                                                                                                                                                                                    <!--     load(paste0("xseq_results/miRNA/", co, "_miRNA.Rdata")) -->
                                                                                                                                                                                                    
                                                                                                                                                                                                    <!--   } -->
                                                                                                                                                                                            
                                                                                                                                                                                            
                                                                                                                                                                                            
                                                                                                                                                                                            
                                                                                                                                                                                            <!--   xseq.pred.mirna <- xseq.pred.all.concat %>%  -->
                                                                                                                                                                                              <!--                     select(sample, gene, prob_mut, prob_gene) %>%  -->
                                                                                                                                                                                              <!--                     distinct() %>%  -->
                                                                                                                                                                                              <!--                     group_by(gene, prob_gene) %>%  -->
                                                                                                                                                                                              <!--                     summarise(N_mutated_samples  = n(), -->
                                                                                                                                                                                                                                   <!--                               N_selected_samples = sum(prob_mut >= 0.5), .groups = "drop") %>%  -->
                                                                                                                                                                                              <!--                     mutate(Ratio = N_selected_samples/N_mutated_samples, -->
                                                                                                                                                                                                                                <!--                            Lab   = ifelse(prob_gene >= xseq.th.mirna[[co]], yes = gene, no = "")) %>%  -->
                                                                                                                                                                                              <!--                     arrange(desc(prob_gene)) -->
                                                                                                                                                                                              <!-- xseq.pred.mirna$Rank <- seq_len(nrow(xseq.pred.mirna)) -->
                                                                                                                                                                                                <!-- xseq.pred.mirna$Lab  <- gsub(xseq.pred.mirna$Lab, pattern = "::.+$", replacement = "", perl = T) -->
                                                                                                                                                                                                  
                                                                                                                                                                                                  
                                                                                                                                                                                                  <!-- rank.plot <- ggplot(xseq.pred.mirna, aes(y = prob_gene, x = Rank, colour = Ratio, size = N_mutated_samples, label = gene)) + -->
                                                                                                                                                                                                    <!--                 geom_line(data = xseq.pred.mirna, aes(x = Rank, y = prob_gene), inherit.aes = F, color = '#9ecae1', alpha = 0.35) +  -->
                                                                                                                                                                                                    <!--                 geom_hline(yintercept = xseq.th.mirna[[co]]) + -->
                                                                                                                                                                                                    <!--                 geom_point(shape = 20) + -->
                                                                                                                                                                                                    <!--                 theme_bw() + -->
                                                                                                                                                                                                    <!--                 ylim(c(0, 1)) + -->
                                                                                                                                                                                                    <!--                 labs(x = "Rank", y = "DAC probability", title = co) + -->
                                                                                                                                                                                                    <!--                 geom_label_repel(aes(label = Lab), -->
                                                                                                                                                                                                                                            <!--                                  box.padding = 0.1, -->
                                                                                                                                                                                                                                            <!--                                  direction = "both",  -->
                                                                                                                                                                                                                                            <!--                                  force = 40, -->
                                                                                                                                                                                                                                            <!--                                  max.overlaps = 30, -->
                                                                                                                                                                                                                                            <!--                                  size = 2.5, -->
                                                                                                                                                                                                                                            <!--                                  color = "black", -->
                                                                                                                                                                                                                                            <!--                                  max.iter = 20000, -->
                                                                                                                                                                                                                                            <!--                                  segment.color = '#d9d9d9') -->
                                                                                                                                                                                                    
                                                                                                                                                                                                    <!-- rank.plot.list[[co]] <- rank.plot -->
                                                                                                                                                                                                      
                                                                                                                                                                                                      <!-- } -->
                                                                                                                                                                                        
                                                                                                                                                                                        <!-- ranks.w.ratio.selected.mut <- plot_grid(rank.plot.list[["BASIS_All"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["BASIS_ER+"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["BASIS_ER-"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["BRCA-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["HNSC-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["LIHC-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["LUAD-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["LUSC-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["STAD-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["UCEC-US"]], -->
                                                                                                                                                                                                                                       <!--                                         ncol = 2, align = "hv", labels = c('A)', 'B)', 'C)', 'D)', 'E)', 'F)', 'G)', "H)", "I)", "J)"), label_size = 30, -->
                                                                                                                                                                                                                                       <!--                                         label_colour = "black", hjust = 0, vjust = 1.3, -->
                                                                                                                                                                                                                                       <!--                                         rel_heights = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1), rel_widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)) -->
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        <!-- ranks.w.ratio.selected.mut -->
                                                                                                                                                                                        <!-- ``` -->
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        <!-- \pagebreak -->
                                                                                                                                                                                        <!-- ```{r Ratio_mutated_vs_selected_PCG, cache=TRUE, include=TRUE, echo=FALSE, eval=TRUE, fig.cap="Ratio of dysregulated vs mutated samples in the protein-coding gene analysis. Each panel represents a cohort with all the evaluated genes ranked by the DAC probability. The color of each dot corresponds to the ratio between the number of samples with SSD >= 0.5 (i.e., with high association of dysregulation) divided by the number of mutated samples. The black line corresponds to the cohort-specific threshold defined by FDR of 0.05.", fig.height = 32, fig.width = 25, fig.align = "center", warning = FALSE} -->
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        <!-- cohorts.all <- c(cohorts, "BASIS_All", "BASIS_ER+", "BASIS_ER-") -->
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        <!-- rank.plot.list <- vector(mode = "list", length = length(cohorts)) -->
                                                                                                                                                                                        <!-- for (co in cohorts.all) { -->
                                                                                                                                                                                            
                                                                                                                                                                                            <!--   if (!grepl(co, pattern = "BASIS")) { -->
                                                                                                                                                                                                <!--     co <- paste0(co, "-US") -->
                                                                                                                                                                                                  
                                                                                                                                                                                                  <!--     load(paste0("xseq_results/PCG/xseq_trans_protein_coding_all_mutation_types_all_effects_", co, "_PPI_network_xseq.Rdata")) ## xseq.pred.all.concat -->
                                                                                                                                                                                                
                                                                                                                                                                                                <!--   } else { -->
                                                                                                                                                                                                    
                                                                                                                                                                                                    <!--     load(paste0("xseq_results/PCG/", co, "_PCG.Rdata")) -->
                                                                                                                                                                                                    
                                                                                                                                                                                                    <!--   } -->
                                                                                                                                                                                            
                                                                                                                                                                                            
                                                                                                                                                                                            <!--   xseq.pred.pcg <- xseq.pred.all.concat %>%  -->
                                                                                                                                                                                              <!--                     filter(mut_effect == "TFBS") %>%  -->
                                                                                                                                                                                              <!--                     select(sample, gene, prob_mut, prob_gene) %>%  -->
                                                                                                                                                                                              <!--                     distinct() %>%  -->
                                                                                                                                                                                              <!--                     group_by(gene, prob_gene) %>%  -->
                                                                                                                                                                                              <!--                     summarise(N_mutated_samples  = n(), -->
                                                                                                                                                                                                                                   <!--                               N_selected_samples = sum(prob_mut >= 0.5), -->
                                                                                                                                                                                                                                   <!--                               .groups = "drop") %>%  -->
                                                                                                                                                                                              <!--                     mutate(Ratio = N_selected_samples/N_mutated_samples, -->
                                                                                                                                                                                                                                <!--                            Lab   = ifelse(prob_gene >= xseq.th.pcg[[co]], yes = gene, no = "")) %>%  -->
                                                                                                                                                                                              <!--                     arrange(desc(prob_gene)) -->
                                                                                                                                                                                              <!--   xseq.pred.pcg$Rank <- seq_len(nrow(xseq.pred.pcg)) -->
                                                                                                                                                                                                
                                                                                                                                                                                                <!--   rank.plot <- ggplot(xseq.pred.pcg, aes(y = prob_gene, x = Rank, colour = Ratio, size = N_mutated_samples, label = gene)) + -->
                                                                                                                                                                                                  <!--                 geom_line(data = xseq.pred.pcg, aes(x = Rank, y = prob_gene), inherit.aes = F, color = '#9ecae1', alpha = 0.35) +  -->
                                                                                                                                                                                                  <!--                 geom_hline(yintercept = xseq.th.pcg[[co]]) + -->
                                                                                                                                                                                                  <!--                 geom_point(shape = 20) + -->
                                                                                                                                                                                                  <!--                 theme_bw() + -->
                                                                                                                                                                                                  <!--                 ylim(c(0, 1)) + -->
                                                                                                                                                                                                  <!--                 labs(x = "Rank", y = "DAC probability", title = co) + -->
                                                                                                                                                                                                  <!--                 geom_label_repel(aes(label = Lab), -->
                                                                                                                                                                                                                                          <!--                                  box.padding = 0.1, -->
                                                                                                                                                                                                                                          <!--                                  direction = "both",  -->
                                                                                                                                                                                                                                          <!--                                  force = 40, -->
                                                                                                                                                                                                                                          <!--                                  max.overlaps = 30, -->
                                                                                                                                                                                                                                          <!--                                  size = 2.5, -->
                                                                                                                                                                                                                                          <!--                                  color = "black", -->
                                                                                                                                                                                                                                          <!--                                  max.iter = 20000, -->
                                                                                                                                                                                                                                          <!--                                  segment.color = '#d9d9d9') -->
                                                                                                                                                                                                  
                                                                                                                                                                                                  <!-- rank.plot.list[[co]] <- rank.plot -->
                                                                                                                                                                                                    
                                                                                                                                                                                                    <!-- } -->
                                                                                                                                                                                        
                                                                                                                                                                                        <!-- ranks.w.ratio.selected.mut <- plot_grid(rank.plot.list[["BASIS_All"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["BASIS_ER+"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["BASIS_ER-"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["BRCA-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["HNSC-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["LIHC-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["LUAD-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["LUSC-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["STAD-US"]], -->
                                                                                                                                                                                                                                       <!--                                         rank.plot.list[["UCEC-US"]], -->
                                                                                                                                                                                                                                       <!--                                         ncol = 2, align = "hv", labels = c('A)', 'B)', 'C)', 'D)', 'E)', 'F)', 'G)', "H)", "I)", "J)"), label_size = 30, -->
                                                                                                                                                                                                                                       <!--                                         label_colour = "black", hjust = 0, vjust = 1.3, -->
                                                                                                                                                                                                                                       <!--                                         rel_heights = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1), rel_widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)) -->
                                                                                                                                                                                        
                                                                                                                                                                                        
                                                                                                                                                                                        <!-- ranks.w.ratio.selected.mut -->
                                                                                                                                                                                        <!-- ``` -->


