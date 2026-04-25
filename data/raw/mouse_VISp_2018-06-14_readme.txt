=================================
Mouse VISp
=================================

RNA sequencing data of single cells isolated from mouse primary visual cortex (VISp).
The data set includes 15413 single cells collected from various cortical layers.
The sequencing results were aligned to exons and introns in the GRCm38.p3 reference genome using the STAR algorithm, and aggregated counts at the gene level were calculated.
For more details, please see the Documentation tab in the Cell Types web application.

Gene expression data matrices
	mouse_VISp_2018-06-14_exon-matrix.csv
		Contains the (row, column) matrix of read counts obtained for each (gene, sample) based on alignment to exons.
		The first row contains the unique identifiers of the samples (cells)
		The first column contains the unique identifiers of the genes
	mouse_VISp_2018-06-14_intron-matrix.csv
		Contains the (row, column) matrix of read counts obtained for each (gene, sample) based on alignment to introns.
		The first row contains the unique identifiers of the samples (cells)
		The first column contains the unique identifiers of the genes

		
Sample information (mouse_VISp_2018-06-14_samples-columns.csv)
	sample_name
		Unique sample identifier
	sample_id
		Unique numeric sample identifier
	sample_type
		Indicates that the samples are cells
	organism
		Donor species
	donor_id
		Unique donor identifier
	sex
		Sex of donor
	age_days
		Age of donor in days
	eye_condition
		Annotation of eye conditions identified in each animal, followed by the side on which the eye was affected. L, left; R, right
	genotype
		Full genotype of donor
	driver_lines
		Cre and/or Flp driver line(s) for each donor
	reporter_lines
		Cre and/or Flp-dependent reporter(s) for each donor
	brain_hemisphere
		Brain hemisphere of dissected cells
	brain_region
		Brain region targeted for sampling
	brain_subregion
		Brain sub-region targeted for sampling	
	injection_label_direction
		Type of injection experiment performed (or "No Injection" if no injection was used)
	injection_primary
		Center of injection site, as verified by inspection of coronal target site slices. Abbreviations match the Allen Mouse Brain Atlas.
	injection_secondary
		All labeled regions, as verified by inspection of coronal target site slices. Abbreviations match the Allen Mouse Brain Atlas.
	injection_tract
		Regions that show labeled cells along the needle insertion injection tract. Abbreviations match the Allen Mouse Brain Atlas.
	injection_material
		Virus injected
	injection_exclusion_criterion
		Reason for excluding injected samples (OK in injection experiment passes)
	facs_date
		The date on which cells were collected by FACS. All dissections were performed on the same date as the FACS date for each sample.
	facs_container
		FACS container unique identifier
	facs_sort_criteria
		Gating criteria used to select cells for sorting
	rna_amplification_set
		Amplification plate
	library_prep_set
		Library plate
	library_prep_avg_size_bp
		Average bp size of Library (Fragment Analyzer™ Automated CE)
	seq_name
		Unique identifier for sequencing
	seq_tube
		Sequencing Lane
	seq_batch
		Sequencing Batch
	total_reads
		Total number of sequencing reads
	percent_exon_reads
		% unique genomic reads aligned to exons (STAR)
	percent_intron_reads
		% unique genomic reads aligned to introns (STAR)
	percent_intergenic_reads
		% unique genomic reads aligned to intergenic sequence (STAR) [corrected 10/2018]
	percent_rrna_reads
		% total reads aligned to rRNA (STAR)
	percent_unique_reads
		% total reads aligned to genome and unique (STAR)
	percent_synth_reads
		% total reads aligned to ERCC synthetic mRNA (STAR)      
	percent_ecoli_reads
		% total reads aligned to E. coli (STAR)
	percent_aligned_reads_total
		% total reads aligned (STAR)
	complexity_cg
		Dinucleotide odds ratio (PRINSEQ)
	genes_detected_cpm_criterion
		# of genes with CPM values greater than 0, intron and exon counts included
	genes_detected_fpkm_criterion
		# of genes with FPKM values greater than 0, only exon counts included
	tdt_cpm
		# of reads mapping to tdTamato transgenic insert sequence, normalized to counts per million
	gfp_cpm
		# of reads mapping to gfp transgenic insert sequence, normalized to counts per million
	class
		Broad cell type class (e.g. GABAergic, Glutamatergic or Non-neuronal), or "No Class" if identified as low-quality or a doublet, or if originally assigned to an outlier or donor-specific cluster 
	subclass
		Cell type sub-class annotation
	cluster
		Cell type cluster name
	confusion_score
		The ratio of the probabilities that each cell was co-clustered with the cells from its second-best cluster and with the cells from its assigned cluster, based on co-clustering scores from 100 bootstrapped rounds of clustering. Confusion score is a measure of the confidence of cell type assignment, the lower the value the higher the confidence.
	cluster_correlation
		Pearson correlation of each cell with the centroid of gene expression for its assigned cluster.
	core_intermediate_call
		Core or Intermediate cell type calls based on centroid classifier validation. Cells that were assigned to their cluster in >90/100 trials are called “Core”. All others are “Intermediate”.

		
Gene information (mouse_VISp_2018-06-14_genes-rows.csv)
	gene
		Gene symbol
	chromosome
		Chromosome location of gene
	entrez_id
		NCBI Entrez ID
	gene_name
		Gene name