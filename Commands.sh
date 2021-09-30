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
#rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/plots/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P1.pdf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P1.pdf

convert -density 500 -trim /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P1.pdf -quality 350 /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P1.jpeg


###############################################################
## Fig S11: Highlighted genes: LUAD-US, STAD-US, and UCEC-US
#rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/plots/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P2.pdf /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_PCG_sep_by_cohort_P2.pdf

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



#####################################
## Dysregulation network fractions ##
#####################################

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


######################
## enrichR heatmaps ##
######################

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/enrichR_analysis/RData/*.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/FEA_heatmaps/

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/results_436_samples/All/miRNA/xseq/trans/targetScan_clean_nb_conn_10000_minW_0.8_minscore_0.8/enrichR_analysis/RData/Enriched_terms_plots_BASIS.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/FEA_heatmaps/


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


########################################################
## Heatmap highlighted PCGs (in at least two cohorts) ##
########################################################

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/RData/xseq_trans_pancancer_heatmap_premiRNA_analysis_mutation_type_TFBS_results_summary_plot.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_premiRNA_analysis_mutation_type_TFBS_results_summary_plot.RData

rsync -rptuvl jamondra@biotin2.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/ICGC/results_ICGC/Pancancer_analysis/premiRNA/xseq/trans/RData/xseq_trans_pancancer_heatmap_PCG_all_mutation_types.Rdata /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/Highlighted_genes_pancancer/xseq_trans_pancancer_heatmap_PCG_all_mutation_types.Rdata


###################################
## ER composition: BASIS vs TCGA ##
###################################

rsync -rptuvl jamondra@biotin3.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/TEST/ER_MATCH/tables/BASIS_ER_status_tab.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/ER_status/BASIS_ER_status_tab.RData

rsync -rptuvl jamondra@biotin3.hpc.uio.no:/storage/mathelierarea/processed/jamondra/Projects/dysmir/BASIS/TEST/ER_MATCH/tables/TCGA_ER_status_tab.RData /home/jamondra/Documents/PostDoc/Mathelier_lab/Manuscripts/dysmiR/dysmiR_supp_notes/ER_status/TCGA_ER_status_tab.RData







