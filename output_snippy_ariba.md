# Output description for SRVCNM169 analysis.

## Pipeline overview

* [FastQC](#fastqc) v0.11.8 - read quality control.
* [Trimmomatic](#trimming) v.0.33 - adapter and low quality trimming.
* [Unicycler](#unicycler) v.0.4.6 - prokaryote assembler.
* [Quast](#quast) v4.1 - assemblies quality control
* [Kmerfinder](#kmerfinder) v.3.1 - species and contamination determination.
* [Snippy](#snippy) v.4.4.0 - haploid variant calling and core genome alignment
* [BamUtil](#bamutil) v.1.0.13 - bam statistics.
* [Picard WGSMetrics](#picard) v1.140 - bam statistics.
* [RAxML](#raxml) v.8.2.9 - maximum likelihood phylogeny.
* [canSNPer](#cansnper) v.1.0.9 - A hierarchical genotype classifier of clonal pathogens

**Note:**

Depending on the analysis, we could have some ANALYSIS_IDs. This ANALYSIS_IDs are going to be composed of the date of the analysis, and an analysis identification. You can find a README in de ANALYSIS folder with a brief description of the different analysis.

## Preprocessing
### FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output directory: `01-fastqc`**

* `{sample_id}/{sample_id}_R[12]_fastqc.html`
  * html report. This file can be opened in your favourite web browser (Firefox/chrome preferable) and it contains the different graphs that fastqc calculates for QC.
* `{sample_id}/{sample_id}_R[12]_fastqc`
  * older with fastqc output in plain text.
* `{sample_id}/{sample_id}_R[12]_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

### Trimming
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (1) is used for removal of adapter contamination and trimming of low quality regions.
Parameters included for trimming are:
-  Nucleotides with phred quality < 10 in 3'end.
-  Mean phred quality < 20 in a 4 nucleotide window.
-  Read lenght < 50

**Results directory: `02-preprocessing`**
- Files:
   - `{sample_id}/{sample_id}_R[12]_filtered.fastq.gz`: contains high quality reads with both forward and reverse tags surviving.
   - `{sample_id}/{sample_id}_R[12]_unpaired.fastq.gz`: contains high quality reads with only forward or reverse tags surviving.

 **Note**:To see how your reads look after trimming, look at the FastQC reports in the 03-preprocQC directory

### Kmerfinder
[Kmerfinder](https://cge.cbs.dtu.dk/services/KmerFinder) (2) is a software used for species indentification and the determination of possible contamination in the sample. We use this software using the bacterial database provided by the developers, and with the "winner takes it all" algorithm. You can check [here](https://cge.cbs.dtu.dk/services/KmerFinder/output.php) for a description of the columns in the output.

**Output directory:** `04-kmerfinder/{sample_id}/`

* `data.json`
  * results in json format.
* `results.spa`
  * results in spa format.
* `results.txt`
  * results in txt format. **This is the format you have to use if you are going to open it with excel.**

**NOTE**: You can also find in `99-stats` a summary of all samples results (`kmerfinder.csv`).

## Phylogenetic analysis
### Snippy
[Snippy](https://github.com/tseemann/snippy) Snippy finds SNPs between a haploid reference genome and your NGS sequence reads. It will find both substitutions (snps) and insertions/deletions (indels). It can then take a set of Snippy results using the same reference and generate a core SNP alignment (and ultimately a phylogenomic tree).

**Output directory:** `05-snippy`

* `{sample_id}`

	* `snps.tab`: A simple tab-separated summary of all the variants
	* `snps.csv`: A comma-separated version of the .tab file
	* `snps.html`: A HTML version of the .tab file
	* `snps.vcf`: The final annotated variants in VCF format
	* `snps.bed`: The variants in BED format
	* `snps.gff`: The variants in GFF3 format
	* `snps.bam`: The alignments in BAM format. Includes unmapped, multimapping reads. Excludes duplicates.
	* `snps.bam.bai`: Index for the .bam file
	* `snps.log`: A log file with the commands run and their outputs
	* `snps.aligned.fa`: A version of the reference but with - at position with depth=0 and N for 0 < depth < --mincov (does not have variants)
	* `snps.consensus.fa`: A version of the reference genome with all variants instantiated
	* `snps.consensus.subs.fa`: A version of the reference genome with only substitution variants instantiated
	* `snps.raw.vcf`: The unfiltered variant calls from Freebayes
	* `snps.filt.vcf`: The filtered variant calls from Freebayes
	* `snps.vcf.gz`: Compressed .vcf file via BGZIP
* `core.aln`: A core SNP alignment in the --aformat format (default FASTA)
* `core.full.aln`: A whole genome SNP alignment (includes invariant sites)
* `core.tab`: Tab-separated columnar list of core SNP sites with alleles but NO annotations
* `core.vcf`: Multi-sample VCF file with genotype GT tags for all discovered alleles
* `core.txt`: Tab-separated columnar list of alignment/core-size statistics
* `core.ref.fa`: FASTA version/copy of the --ref
* `core.self_mask.bed`: BED file generated if --mask auto is used.

### BamUtil
[BamUtil](https://genome.sph.umich.edu/wiki/BamUtil) provides useful stats for the mapping step.

**Output directory:** `99-stats`
* `bamstat.csv`: summary file with all bamUtil information in one file for all samples.

### Picard WGSMetrics
[Picard WGSMetrics](https://broadinstitute.github.io/picard/picard-metric-definitions.html) provides useful stats for the mapping step.

**Output directory:** `99-stats`
* `wgsmetrics.csv`: summary file with all Picard WGSMetrics information in one file for all samples.

### RAxML
[RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) (6) performs a maximum likelihood phylogeny from a multiple sequence alignment.

**Output directory:** `07-raxml`

* `RAxML_bipartitions.RAXML*`: phylogenetic tree in newick format with bootstrap values.

## Strain characterization
### Ariba
[Ariba](https://github.com/sanger-pathogens/ariba) (7) is a tool that identifies antibiotic resistance genes by running local assemblies. It can also be used for MLST calling.

**Output directory:** `08-ariba`

* `get_prep_ref`: databases references for resistance, virulence, etc.
* `db_*_mlst`: mlst reference database.
* `run_mlst`: mlst results.
	* `{sample_id}/{database_run_folder}`
		* `{sample_id}report.tsv`: report with detected genes (field explanation in [ANEXI](#anexi))
		* `assembled_genes.fa.gz`: sequence for the assembled resistance/virulence/etc genes.
* `run_dbs`: databases results.
	* `{sample_id}/{database_run_folder}`
		* `{sample_id}report.tsv`: report with detected genes (field explanation in [ANEXI](#anexi))
		* `assembled_genes.fa.gz`: sequence for the assembled resistance/virulence/etc genes.
* `summary`: summary files for all samples.

## Custom analysis.

## Final results.
We provide a folder with the most relevant results, in most of the cases these are the results you are interested in, and the ones you need to check and interpret.

**Output directory:** `RESULTS`
* `snps_raxml_100bootstrap.nwk`: Phylogenetic tree obtained with RAxML and SNIPPY+Gubbins SNP matrix.
* `summary_ariba`: folder with summary data from ariba (resistance, virulence, mlst, plasmids).
* `wgsmetrics_all.csv`: mapping summary statistics for all samples.
* `kmerfinder.csv`: identification and contamination detection for all samples.


## REFERENCES
1. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
2. Rapid and precise alignment of raw reads against redundant databases with KMA Philip T.L.C. Clausen, Frank M. Aarestrup, Ole Lund.
3. Wick RR, Judd LM, Gorrie CL, Holt KE. Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 2017.
4. Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics 2014 Jul 15;30(14):2068-9.
5. Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi and Glenn Tesler. QUAST: quality assessment tool for genome assemblies, Bioinformatics (2013) 29 (8): 1072-1075.
6. A. Stamatakis: "RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies". In Bioinformatics, 2014
7. Hunt et al. ARIBA: rapid antimicrobial resistance genotyping directly from sequencing reads. Microbiology society.2017.

## ANEXES
### ANEX I

| Column                 | Description |
|------------------------|-------------|
| 1.  ariba_ref_name     | ariba name of reference sequence chosen from cluster (needs to rename to stop some tools breaking) |
| 2.  ref_name           | original name of reference sequence chosen from cluster, before renaming |
| 3.  gene               | 1=gene, 0=non-coding (same as metadata column 2) |
| 4.  var_only           | 1=variant only, 0=presence/absence (same as metadata column 3) |
| 5.  flag               | cluster flag |
| 6.  reads              | number of reads in this cluster |
| 7.  cluster            | name of cluster |
| 8.  ref_len            | length of reference sequence |
| 9.  ref_base_assembled | number of reference nucleotides assembled by this contig |
| 10.  pc_ident           | %identity between reference sequence and contig |
| 11. ctg                | name of contig matching reference |
| 12. ctg_len            | length of contig |
| 13. ctg_cov            | mean mapped read depth of this contig |
| 14. known_var          | is this a known SNP from reference metadata? 1 or 0 |
| 15. var_type           | The type of variant. Currently only SNP supported |
| 16. var_seq_type       | Variant sequence type. if known_var=1, n or p for nucleotide or protein |
| 17. known_var_change   | if known_var=1, the wild/variant change, eg I42L |
| 18. has_known_var      | if known_var=1, 1 or 0 for whether or not the assembly has the variant |
| 19. ref_ctg_change     | amino acid or nucleotide change between reference and contig, eg I42L |
| 20. ref_ctg_effect     | effect of change between reference and contig, eg SYS, NONSYN (amino acid changes only) |
| 21. ref_start          | start position of variant in reference |
| 22. ref_end            | end position of variant in reference |
| 23. ref_nt             | nucleotide(s) in reference at variant position |
| 24. ctg_start          | start position of variant in contig |
| 25. ctg_end            | end position of variant in contig |
| 26. ctg_nt             | nucleotide(s) in contig at variant position |
| 27. smtls_total_depth  | total read depth at variant start position in contig, reported by mpileup |
| 28. smtls_nts          | nucleotides on contig, as reported by mpileup. The first is the contig nucleotide |
| 29. smtls_nts_depth    | depths on contig, as reported by mpileup. One number per nucleotide in the previous column |
| 30. var_description    | description of variant from reference metdata |
| 31. free_text          | other free text about reference sequence, from reference metadata |

