# Predicting core p53 target lncRNAs and their functions
This page provides the following analysis sections used for the paper: **Decoding the lncRNAome across diverse cellular stresses reveals core p53-effector pan-cancer suppressive lncRNAs.**
 
 **1. Meta_analysis** - Meta-analysis was conducted based on 14 RNA-seq datasets where cells were treated with either p53-activating stimuli or vehicle control.
 
 &emsp;**R script** 
 	
	Meta_analysis/meta_analysis.R - The R script contains the analysis code to predict core p53-effector lncRNAs. 
	
 
 &emsp;**Data files (Input)** 
 	
	a) Meta_analysis/input_RNAseq - Contains gene expression fold-change matrices for 14 RNA-seq datasets.
	b) Meta_analysis/gencode_v27.gtf - List of lncRNAs present in Gencode v27 database (PMID:30357393).
	c) Meta_analysis/p53-target_lnc_chipbase_atleast_5_data.txt - List of potential p53-regulated lncRNAs that have the p53-binding site(s) in five or more ChIP-seq datasets. The data were obtained from ChIPBase v2.0 database (PMID:27924033).
 
 &emsp;**Results (Output)**
 
 	a) Meta_analysis/res_rra.csv - Meta-analysis results showing potential p53-effector lncRNAs, their expression levels across 14 RNA-seq datasets, meta-analysis p-values, and the number of p53 ChIP-seq datasets where the reported lncRNAs were directly bound by p53. 
	b) The script generates Figure 2E and source data for Supplementary Table S4.
	

**2. TCGA_pathway** - KEGG pathway and Hallmark geneset enrichment analysis across ten TCGA cancer types to predict p53-effector lncRNA functions. 

&emsp;**R script** 

	TCGA_pathway/TCGA_hallmark_kegg_enrichment.R - The R script contains the analysis code to identify p53-effector lncRNA-pathway associations.

&emsp;**Data files (Input)** 

	TCGA_pathway/*/lnc_mRNA_[negative/positive]_regression.GMT - Negatively or positively associated lncRNA-mRNA pairs in indicated TCGA cancer type.
	hallmark_categories.txt - Grouping of Hallmark genesets into eight broad categories.
	
	*: TCGA cancer type (e.g. 'blca')

**Results (Output)**

	Figure 4B, Supplementary Tables S4-S8, source data for Figure 4a and supplementary Figures S3-S5. 



**3. CRISPR_RNAi** - Analysis of CRISPR and RNAi screening data across cancer cell lines of ten cancer types to predict cancer cell survival/growth associated lncRNAs. 
	
&emsp;**R scripts**

	Following R scripts were used in the generation of Figures and Tables related to RNAi and CRISPR screening data analyses-
		a) CRISPR_RNAi/RNAi/lncRNA_RNAi.R  - Analyzes RNAi screening data
		b) CRISPR_RNAi/CRISPR/lncRNA_CRISPR.R - Analyzes CRISPR screening data 
		c) CRISPR_RNAi/CRISPR_RNAi_TCGA_correlation/correlation.R - Computes correlations between the results obtained from the RNAi and CRISPR screening data and TCGA patient data.
			
&emsp;**Data files (Input)** 
		
	CRISPR_RNAi/RNAi - 
		a) Download 'D2_combined_gene_dep_scores.csv' file from depmap DEMETER2 Data v6 (https://depmap.org/portal/download/all/)  
		b) depmap/*_mapped_cell_lines.csv - list of cell lines
		c) */lnc_mRNA_[negative/positive]_regression.GMT - Negatively or positively associated lncRNA-mRNA pairs in indicated TCGA cancer type

	CRISPR_RNAi/CRISPR - 
		a) Download 'Achilles_gene_effect.csv' file from depmap version 19Q3 (https://depmap.org/portal/download/all/)
		b) depmap/*_mapped_cell_lines.csv - list of cell lines
		c) cellline_id_to_name.txt - Mapping of cell line id to cell line name
		d) */lnc_mRNA_[negative/positive]_regression.GMT - Negatively or positively associated lncRNA-mRNA pairs in indicated TCGA cancer type

	CRISPR_RNAi/CRISPR_RNAi_TCGA_correlation -
		a) crispr_res_p53_effector_lnc_suppress_growth_fdr_0.001.csv - p53 effector lncRNAs potentially suppress cancer cell survival/growth (CRISPR data)
		b) rnai_res_p53_effector_lnc_suppress_growth_fdr_0.001.csv - p53 effector lncRNAs potentially suppress cancer cell survival/growth (RNAi data)
		c) lncRNA_alter_proliferation_cancer_number.csv - p53 effector lncRNAs potentially suppress/induce proliferation across TCGA cancer types
	
	*: TCGA cancer type (e.g. 'blca')

**Results (Output)**
		
		a) CRISPR_RNAi/[CRISPR/RNAi]/*/[negative/positive]_regression/[negative/positive]_regression_res_gsea_*.csv - Gene Set Enrichment Analysis (GSEA) results show which lncRNAs potentially suppress (negative normalized enrichment score (NES)) or induce (positive NES) survival/growth of specific cancer celllines with statistical significance. 	
		b) Figures 2C-2H and source data for Supplementary Tables 9-12.
		
	*: TCGA cancer type (e.g. 'blca')


**4. PTSL_RNAseq** - RNA-sequencing was performed to evaluate gene expression changes following ectopic expression of PTSL or vector control in Lung Adenocarcinoma (LUAD) cells.

&emsp;**R script** 

	PTSL_RNAseq/PTSL_RNAseq.R - The R script contains the analysis code to determine altered genes and pathways in A549 LUAD cells overexpressing PTSL compared to controls.

&emsp;**Data files (Input)** 

	a) RNA-seq data can be downloaed from Gene Expression Omnibus (Accession ID is provided in the article).
	b) PTSL_RNAseq/regression_tcga_luad/lnc_mRNA_[negative/positive]_regression.GMT - Negatively or positively associated lncRNA-mRNA pairs in TCGA LUAD.
	c) PTSL_RNAseq/ts_genes.txt - List of known tumor-suppressor genes significantly (atleast 2 fold-change with FDR < 0.05) reduced in TCGA LUAD compared to normal samples.
	d) PTSL_RNAseq/onco_genes.txt - List of known oncogenes significantly (atleast 2 fold-change with FDR < 0.05) elevated in TCGA LUAD compared to normal samples.
	e) tcga_luad_favorable_genes.txt - List of genes whose elevated levels increase LUAD patient's overall survival.
	f) tcga_luad_unfavorable_genes.txt - List of genes whose elevated levels reduce LUAD patient's overall survival.

**Results (Output)**
	Differentially expressed genes in A549 LUAD cells overexpressing PTSL compared to controls, Figures 5C-E, source data for Figures 5B and 7D.


TCGA : The Cancer Genome Atlas; * : TCGA cancer type (e.g. 'blca'); FDR: False Discovery Rate
