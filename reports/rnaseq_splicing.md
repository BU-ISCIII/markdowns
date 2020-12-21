# Output description for rnaseq-nf pipeline
**rnaseq-nf** is a bioinformatics best-practice analysis pipeline used for RNA-seq data. The pipeline is focused in counts and differential expression analyses.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* <a href="#fastqc">FastQC</a> v0.11.8 - read quality control
* <a href="#trimming">Trimmomatic</a> v.0.38 - adapter and low quality trimming.
* <a href="#star">STAR</a> v2.6.1d - alignment
* <a href="#rseqc">RSeQC</a> v3.0.0 - RNA quality control metrics
  * <a href="#bamstat">Bam stat</a>
  * <a href="#clippingprofile">Clipping profile</a>
  * <a href="#genebodycoverage">Gene body coverage</a>
  * <a href="#inferexperiment">Infer Experiment</a>
  * <a href="#innerdistance">Inner distance</a>
  * <a href="#junctionannotation">Junction annotation</a>
  * <a href="#junctionsaturation">Junction saturation</a>
  * <a href="#readdistribution">Read distribution</a>
  * <a href="#readduplication">Read duplication</a>
* <a href="#preseq">Preseq</a> v2.0.3 - library complexity
* <a href="#picard">Picard</a> v2.18.27 - Identify duplicate reads
* <a href="#dupradar">dupRadar</a> v1.12.1 - technical / biological read duplication
* <a href="#featurecounts">Subread(Featurecounts)</a> v1.6.4 - gene counts, biotype counts, rRNA estimation.
* <a href="#stringtie">StringTie</a> v1.3.5 - FPKMs for genes and transcripts
* <a href="#edger">edgeR</a> v3.24.1 - create MDS plot and sample pairwise distance heatmap / dendrogram
* <a href="#multiqc">MultiQC</a> v1.7 - aggregate report, describing results of the whole pipeline
* Custom analysis:
  * <a href="#gffcompare">GFF Compare</a> v0.11.2 - Compare reference GFF file and Stringtie's de novo generated GFF files.
  * <a href="#deseq2">DESeq2</a> v1.18.1 - Differential expression analysis and plots

## Preprocessing
### FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)<a name="fastqc"> </a><a href="#fastqc_reference">[1]</a> gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other over-represented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

  **Note**:The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the 02-preprocessing/FastQC directory.

**Output directory:** `01-fastqc/``
* `logs/`
  * Folder with the log files per sample of the FastQC.
* `zips/`
  * `{sample_id}_R[12]_fastqc.zip`
    * zip file containing the FastQC report, tab-delimited data file and plot images

### Trimming
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)<a name="trimming"> </a><a href="#trimming_reference">[2]</a> is used for removal of adapter contamination and trimming of low quality regions.
Parameters included for trimming are:
-  Nucleotides with phred quality < 10 in 3'end.
-  Mean phred quality < 20 in a 4 nucleotide window.
-  Read lenght < 50

MultiQC reports the percentage of bases removed by trimming in bar plot showing percentage or reads trimmed in forward and reverse.

**Results directory:** `02-preprocessing/`
* `FastQC/`
  * `{sample_id}_filtered_R[12].fastqc.html`: html report of the trimmed reads.
  * `{sample_id}_filtered_R[12].fastqc.html.zip`: zip compression of above file.
* `logs/`
  * Folder with the log files per sample of the FastQC.

**NOTE:** Trimmed reads are not delivered to the researcher by default due to disk space issues. If you are interested in using them, please contact us and we will add them to your delivery.

## Alignment
### STAR
[STAR](https://github.com/alexdobin/STAR)<a name="star"> </a><a href="#STAR_reference">[3]</a> is a read aligner designed for RNA sequencing.  STAR stands for Spliced Transcripts Alignment to a Reference, it produces results comparable to TopHat (the aligned previously used by NGI for RNA alignments) but is much faster.

The STAR section of the MultiQC report shows a bar plot with alignment rates: good samples should have most reads as _Uniquely mapped_ and few _Unmapped_ reads.

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/star_alignment_plot.png" width="600" height="350" align="middle" alt="STAR alignment plot"/>

**Output directory:** [03-alignment/](../ANALYSIS/03-alignment/)

* `{sample_id}_filteredAligned.sortedByCoord.out.bam`
  * The aligned BAM file
* `{sample_id}_filteredAligned.sortedByCoord.out.bam.bai`
  * The aligned BAM file index
* `logs/`
  * `{sample_id}_filteredLog.final.out`
    * The STAR alignment report, contains mapping results summary
  * `{sample_id}_filteredLog.out` and `{sample_id}filteredLog.progress.out`
    * STAR log files, containing a lot of detailed information about the run. Typically only useful for debugging purposes.
  * `{sample_id}_filteredSJ.out.tab`
    * Filtered splice junctions detected in the mapping

## Alignment Quality Control
### RSeQC
<a name="rseqc">RSeQC</a> <a href="#RSeQC_reference">[4]</a> is a package of scripts designed to evaluate the quality of RNA seq data. You can find out more about the package at the [RSeQC website](http://rseqc.sourceforge.net/).

This pipeline runs several, but not all RSeQC scripts. All of these results are summarized within the MultiQC report and described below.

**Output directory:** `04-rseqc/`

These are all quality metrics files and contains the raw data used for the plots in the MultiQC report. In general, the `.r` files are R scripts for generating the figures, the `.txt` are summary files, the `.xls` are data tables and the `.pdf` files are summary figures.

#### Bam stat
**Output:** <a name="bamstat">bam_stat/</a>

This script gives numerous statistics about the aligned BAM files produced by STAR. A typical output looks as follows:

**From the `{sample_id}.bam_stat.txt` file:**

```txt
#Output (all numbers are read count)
#==================================================
Total records:                                 41465027
QC failed:                                     0
Optical/PCR duplicate:                         0
Non Primary Hits                               8720455
Unmapped reads:                                0

mapq < mapq_cut (non-unique):                  3127757
mapq >= mapq_cut (unique):                     29616815
Read-1:                                        14841738
Read-2:                                        14775077
Reads map to '+':                              14805391
Reads map to '-':                              14811424
Non-splice reads:                              25455360
Splice reads:                                  4161455
Reads mapped in proper pairs:                  21856264
Proper-paired reads map to different chrom:    7648
```

MultiQC plots each of these statistics in a dot plot. Each sample in the project is a dot - hover to see the sample highlighted across all fields.

RSeQC documentation: [bam_stat.py](http://rseqc.sourceforge.net/#bam-stat-py)

#### Clipping Profile
**Output directory:** <a name="clippingprofile">clipping_profile/</a>

This program is used to estimate clipping profile of RNA-seq reads from BAM or SAM file. Note that to use this function, CIGAR strings within SAM/BAM file should have ‘S’ operation (This means your reads aligner should support clipped mapping).

* `data/`
  * `{sample_id}.clipping_profile.xls`
    * Contains 3 columns: the first column is position (starting from 0) of read in 5’->3’ direction; the second column is the number of reads clipped at this position; the third column is the number of reads non-clipped at this position.
* `plots/`
  * `{sample_id}.clipping_profile.pdf`
    * PDF with the following graph(s)

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/clipping_profile.png" width="550" height="600" align="middle" alt="Clipping profile"/>

RSeQC documentation: [clipping_profile.py] http://rseqc.sourceforge.net/#clipping-profile-py

#### Infer Experiment
**Output:** <a name="inferexperiment">infer_experiment/</a>

This script predicts the mode of library preparation (sense-stranded or antisense-stranded) according to how aligned reads overlay gene features in the reference genome.
Example output from an unstranded (~50% sense/antisense) library of paired end data:

**From MultiQC report:**
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/rseqc_infer_experiment_plot.png" width="600" height="400" align="middle" alt="Infer experiment"/>

**From the `{sample_id}.infer_experiment.txt` file:**

```txt
This is PairEnd Data
Fraction of reads failed to determine: 0.0409
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4839
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4752
```

RSeQC documentation: [infer_experiment.py](http://rseqc.sourceforge.net/#infer-experiment-py)

#### Inner distance
Only for paired-end data. The <a name="innerdistance">inner distance</a> script tries to calculate the inner distance between two paired RNA reads. It is the distance between the end of read 1 to the start of read 2, and it is sometimes confused with the insert size (see [this blog post](http://thegenomefactory.blogspot.com.au/2013/08/paired-end-read-confusion-library.html)).

**Output:** `inner_distance/`

* `{sample_id}.inner_distance.txt`
  * First column is read ID
  * Second column is inner distance. Could be negative value if PE reads were overlapped or mapping error (e.g. Read1_start < Read2_start, while Read1_end >> Read2_end due to spliced mapping of read1).
  * Third column indicates how paired reads were mapped: PE_within_same_exon, PE_within_diff_exon,PE_reads_overlap.
* `data/`
  * `{sample_id}.inner_distance_freq.txt`
    * Inner distance starts
    * Inner distance ends
    * Number of read pairs
    * Note the first 2 columns are left side half open interval
* `plots/`
  * `{sample_id}.inner_distance.pdf`
    * Histogram representing the previous txt file's information.

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/inner_distance_concept.png" width="500" height="250" align="middle" alt="Inner distance concept"/>

> _Credit: modified from RSeQC documentation._

Note that values can be negative if the reads overlap. A typical set of samples may look like this:
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/rseqc_inner_distance_plot.png" width="600" height="450" align="middle" alt="Inner distance plot"/>

This plot will not be generated for single-end data. Very short inner distances are often seen in old or degraded samples (_eg._ FFPE).

RSeQC documentation: [inner_distance.py](http://rseqc.sourceforge.net/#inner-distance-py)

#### Junction annotation
<a name="junctionannotation">Junction annotation</a> compares detected splice junctions to a reference gene model. An RNA read can be spliced 2 or more times, each time is called a splicing event.

**Output:** `junction_annotation/`

* `logs/`
* `{sample_id}.junction_annotation_log.txt`
  * Log files of the analysis.
* `data/`
* `{sample_id}.junction.xls`
  * Data tables used to generate the PDF plots.
  * First column is chromosome ID.
  * Start position of junction (0 based)
  * End position of junction (1 based)
  * Number of splice events supporting this junction
  * 'Annotated', 'complete_novel' or 'partial_novel'
* `events/`
  * `{sample_id}.splice_events.pdf`
    * Shows the distribution of the splicing events found.
* `junctions/`
  * `{sample_id}.splice_junction.pdf`
    * Shows the distribution of the splicing junctions found.

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/rseqc_junction_annotation_junctions_plot.png" width="600" height="400" align="middle" alt="Junction annotation"/>

RSeQC documentation: [junction_annotation.py](http://rseqc.sourceforge.net/#junction-annotation-py)


#### Junction saturation
**Output:** <a name="junctionannotation">junction_saturation/</a>

* `{sample_id}.junctionSaturation_plot.pdf`
  * Shows the number of splice sites detected at the data at various levels of subsampling. A sample that reaches a plateau before getting to 100% data indicates that all junctions in the library have been detected, and that further sequencing will not yield more observations. A good sample should approach such a plateau of _Known junctions_, very deep sequencing is typically requires to saturate all _Novel Junctions_ in a sample.

None of the lines in this example have plateaued and thus these samples could reveal more alternative splicing information if they were sequenced deeper.

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/rseqc_junction_saturation_plot.png" width="600" height="400" align="middle" alt="Junction saturation"/>

RSeQC documentation: [junction_saturation.py](http://rseqc.sourceforge.net/#junction-saturation-py)


#### Read Distribution
**Output:**<a name="readdistribution">read_distribution/</a>

This tool calculates how mapped reads are distributed over genomic features. A good result for a standard RNA seq experiments is generally to have as many exonic reads as possible (`CDS_Exons`). A large amount of intronic reads could be indicative of DNA contamination in your sample or some other problem.

* `{sample_id}.read_distribution.txt`: The ouput is a txt table with four columns:
    * Total_bases: This does NOT include those QC fail,duplicate and non-primary hit reads
    * Tag_count: Number of tags that can be unambiguously assigned the 10 groups.
    * Tags/Kb: Tags per kilobase

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/rseqc_read_distribution_plot.png" width="600" height="400" align="middle" alt="Read distribution"/>

RSeQC documentation: [read_distribution.py](http://rseqc.sourceforge.net/#read-distribution-py)


#### Read duplication
**Output:**<a name="readduplication">read_duplication/</a>

* `{sample_id}.read_duplication.DupRate_plot.pdf`
  * This plot shows the number of reads (y-axis) with a given number of exact duplicates (x-axis). Most reads in an RNA-seq library should have a low number of exact duplicates. Samples which have many reads with many duplicates (a large area under the curve) may be suffering excessive technical duplication.
* `dup_pos/`
  * `dup_pos/{sample_id}.read_duplication.pos.DupRate.xls`
    * Table with the read duplication rate determined from mapping position of read. First column is “occurrence” or duplication times, second column is number of uniquely mapped reads.
* `dup_seq/`
  * `dup_seq/{sample_id}.read_duplication.seq.DupRate.xls`
    * Table with th read duplication rate determined from sequence of read. First column is “occurrence” or duplication times, second column is number of uniquely mapped reads.

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/rseqc_read_dups_plot.png" width="600" height="450" align="middle" alt="Read duplication"/>

RSeQC documentation: [read_duplication.py](http://rseqc.sourceforge.net/#read-duplication-py)


## Counts
### Preseq
[Preseq](http://smithlabresearch.org/software/preseq/)<a name="preseq"> </a><a href="#preseq_reference">[5]</a> estimates the complexity of a library, showing how many additional unique reads are sequenced for increasing the total read count. A shallow curve indicates that the library has reached complexity saturation and further sequencing would likely not add further unique reads. The dashed line shows a perfectly complex library where total reads = unique reads.

Note that these are predictive numbers only, not absolute. The MultiQC plot can sometimes give extreme sequencing depth on the X axis - click and drag from the left side of the plot to zoom in on more realistic numbers.

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/preseq_plot.png" width="600" height="450" align="middle" alt="Preseq complexity curve"/>

**Output directory:** `05-preseq/`

* `{sample_id}.ccurve.txt`
  * This file contains plot values for the complexity curve, plotted in the MultiQC report.

### Picard
[Picard](https://broadinstitute.github.io/picard/index.html)<a name="picard"> </a><a href="#picard_reference">[6]</a> is a set of command line tools for manipulating high-throughput sequencing (HTS) data. In this case we used it to locate and tag duplicate reads in BAM files.

**Output directory:** `06-removeDuplicates/picard/`

* `metrics/`
  * `metrics/{sample_id}.markDups_metrics.txt`
    * Read duplication metrics.

Picard documentation: [Picard docs](https://broadinstitute.github.io/picard/command-line-overview.html)

### dupRadar
[dupRadar](https://www.bioconductor.org/packages/release/bioc/html/dupRadar.html)<a name="duparadar"> </a><a href="#dupradar_reference">[7]</a> is a Bioconductor library for R. It plots the duplication rate against expression (RPKM) for every gene. A good sample with little technical duplication will only show high numbers of duplicates for highly expressed genes. Samples with technical duplication will have high duplication for all genes, irrespective of transcription level.

> _Credit: [dupRadar documentation](https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html)_

**Output directory:** `06-removeDuplicates/dupRadar/`

* `{sample_id}.markDups_dup_intercept_mqc.txt`
  * General stats file.
* `box_plots/`
  * `box_plots/{sample_id}.markDups_duprateExpBoxplot.pdf`
    * These box plots show the percentage of read duplication vs the mean expression. The boxplot in the left side represents a good experiment.
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/dupRadar_boxplot.png" width="750" height="400" align="middle" alt="dupRadar boxplot"/>
* `gene_data/`
  * `gene_data/{sample_id}.markDups_dupMatrix.txt`
    * The duplication matrix contains read counts in different scenarios and RPK and RPKM values for every gene.
* `histograms/`
  * `histograms/{sample_id}.markDups_expressionHist.pdf`
    * This histogram represents the distribution of RPK values per genes. This would help in identifying skewed distributions with unusual amount of lowly expressed genes, or to detect no consensus between replicates.
* `intercepts_slopes/`
  * `intercepts_slopes/{sample_id}.markDups_intercept_slope.txt`
* `scatter_curve_data/`
  * `scatter_curve_data/{sample_id}.markDups_duprateExpDensCurve_mqc.tx`
    * Provides duplication rate quality control for RNA-Seq datasets. Highly expressed genes can be expected to have a lot of duplicate reads, but high numbers of duplicates at low read counts can indicate low library complexity with technical duplication.
* `scatter_plots/`
  * `scatter_plots/{sample_id}.markDups.bam_duprateExpDens.pdf`
    * This plot relates the normalized number of counts per gene (RPK, as a quantification of the gene expression) and the fraction represented by duplicated reads. A good experiment is shown in the left pic.
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/dupRadar_plot.png" width="700" height="400" align="middle" alt="dupRadar scatter plot"/>

DupRadar documentation: [dupRadar docs](https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html)

### featureCounts
[featureCounts](http://bioinf.wehi.edu.au/featureCounts/)<a name="featurecounts"> </a><a href="#featurecounts_reference">[8]</a> from the subread package summarizes the read distribution over genomic features such as genes, exons, promotors, gene bodies, genomic bins and chromosomal locations.
RNA reads should mostly overlap genes, so be assigned.

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/featureCounts_assignment_plot.png" width="600" height="400" align="middle" alt="featureCounts assignment plot"/>

We also use featureCounts to count overlaps with different classes of features. This gives a good idea of where aligned reads are ending up and can show potential problems such as rRNA contamination.
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/featureCounts_biotype_plot.png" width="600" height="400" align="middle" alt="featureCounts biotype"/>

**Output directory:** `07-featureCounts/`

* `biotype_counts/{sample_id}_biotype_counts_mqc.txt`
  * Read counts for the different gene biotypes that featureCounts distinguishes.
* `gene_counts/{sample_id}_gene.featureCounts.txt`
  * Read the counts for each gene provided in the reference `gtf` file
* `gene_counts_summaries/{sample_id}_gene.featureCounts.txt.summary`
  * Summary file, containing statistics about the counts
* `merged_gene_counts.txt`
  * File with gene counts merged between all samples.

### StringTie
[StringTie](https://ccb.jhu.edu/software/stringtie/)<a name="stringtie"> </a><a href="#StringTie_reference">[9]</a> assembles RNA-Seq alignments into potential transcripts. It assembles and quantitates full-length transcripts representing multiple splice variants for each gene locus.

StringTie outputs FPKM metrics for genes and transcripts as well as the transcript features that it generates.

**Output directory:** `08-stringtieFPKM/`

* `{sample_id}.gene_abund.txt`
  * Gene abundances, FPKM values.
* `ballgown/`
  * Folder with ballgown data for each sample. For more information about these files and its usage read [this page](https://github.com/alyssafrazee/ballgown)
* `transcipts/`
  * `transcripts/{sample_id}_transcripts.gtf`
    * This `.gtf` file contains all of the assembled transcipts from StringTie
* `cov_refs/`
  * `{sample_id}.cov_refs.gtf`
    * This `.gtf` file contains the transcripts that are fully covered by reads.

## Sample correlation and Differential Expression (DE)

### edgeR
Sample correlation with <a name="edger">edgeR</a>. [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) <a href="#edgeR_reference">[10]</a> is a Bioconductor package for R used for RNA-seq data analysis. The script included in the pipeline uses edgeR to normalize read counts and create a heatmap / dendrogram showing pairwise euclidean distance (sample similarity). It also creates a 2D MDS scatter plot showing sample grouping. These help to show sample similarity and can reveal batch effects and sample groupings.

**Heatmap:**

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/mqc_hcplot_hocmzpdjsq.png" width="600" height="600" align="middle" alt="edgeR heatmap"/>

**MDS plot:**

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/mqc_hcplot_ltqchiyxfz.png" width="600" height="400" align="middle" alt="edgeR MDS plot"/>

**Output directory:** `09-sample_correlation/`

* `edgeR_MDS_plot.pdf`
  * MDS scatter plot, showing sample similarity
* `edgeR_MDS_distance_matrix.txt`
  * Distance matrix containing raw data from MDS analysis
* `edgeR_MDS_plot_coordinates_mqc.txt`
  * Scatter plot coordinates from MDS plot, used for MultiQC report
* `log2CPM_sample_distances_dendrogram.pdf`
  * Dendrogram plot showing the euclidian distance between your samples
* `log2CPM_sample_distances_heatmap.pdf`
  * Heatmap plot showing the euclidian distance between your samples
* `log2CPM_sample_distances_mqc.csv`
  * Raw data used for heatmap and dendrogram plots.

## Final reports
### MultiQC
[MultiQC](http://multiqc.info)<a name="multiqc"> </a><a href="#MultiQC_reference">[11]</a> is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualized in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory:** `99-stats/MultiQC/`

* `multiqc_report.html`
  * A standalone HTML file that can be viewed in your web browser
* `multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)

## Custom analysis
This part contains the customized analyses realized for the researcher specifically for this analysis.

### New splicing events.

To predict new splicing events, we will first use <a name="stringtie">StringTie</a> in a different mode to allow it to predict new transcript that are not contained in the reference GTF file.

**Output directory:** `10-stringtieSplice`
* `{sample_id}.gene_abund.txt`
  * Gene abundances, FPKM values.
* `ballgown/`
  * Folder with ballgown data for each sample. For more information about these files and its usage read [this page](https://github.com/alyssafrazee/ballgown)
* `transcipts/`
  * `transcripts/{sample_id}_transcripts.gtf`
    * This `.gtf` file contains all of the assembled transcripts from StringTie predicting de novo transcripts.
* `cov_refs/`
  * `{sample_id}.cov_refs.gtf`
    * This `.gtf` file contains the transcripts that are fully covered by reads, including de novo transcripts.

These new transcript files generated by StringTie will be further used in <a href="#gffcompare">GFF Compare</a> <a href="#gffcompare_reference">[12]</a>. GFF compare is a program can be used to compare, merge, annotate and estimate accuracy of one or more GFF files (the “query” files), when compared with a reference annotation (also provided as GFF).

**Output directories:**
* `10-stringtieSplice`
  * `transcripts/{sample_id}_transcripts.gtf.refmap`: This tab-delimited file lists, for each reference transcript, which query transcript either fully or partially matches that reference transcript.
  * `transcripts/{sample_id}_transcripts.gtf.refmap`: This tab delimited file lists the most closely matching reference transcript for each query transcript. Classification of the transcripts based in the following schema:

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/gffcompare_codes.png" width="600" height="400" align="middle" alt="gff compare classification"/>

* `12-gffcompare/`
  * `{sample_id}.annotation.gtf`: Annotation GFT file of the transcripts
  * `{sample_id}.loci`: File with the LOCI of the transcripts.
  * `{sample_id}.stats`: various statistics related to the “accuracy” (or a measure of agreement) of the input transcripts when compared to reference annotation data. Sensitivity and Precision values are estimated at various levels, which are largely an increasingly stringent way of evaluating the accuracy/correctness of a set of predicted (assembled) transcripts (transfrags), when compared to the reference annotation provided with the -r option. The levels are as follows:
  * `{sample_id}.trackings`: This file matches transcripts up between samples. Each row represents a transcript structure that is preserved (structurally equivalent) across all the input GTF files. GffCompare considers transcripts "matching" (i.e. structurally equivalent) if all their introns are identical. Note that "matching" transcripts are allowed to differ on the length of the first and last exons, since these lengths can usually vary across samples for the same biological transcript. *** In our case each sample was analyzed once agains the reference file so this file it's not sample relative relevant.**

> _Credit: For more information visit [GFF compare documentation](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)_

### DESeq2
Differential expression analysis with <a name="deseq2">DESeq2</a>. [DESeq2](https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf) <a href="#DESeq2_reference">[13]</a> is a Bioconductor package for R used for RNA-seq data analysis. The script included in the pipeline uses DESeq2 to normalize read counts and create a heatmap / dendrogram showing pairwise euclidean distance (sample similarity). It also creates other plots to evaluate the sample dispersion. It also provides PCA plots to evaluate sample grouping.

**MA plot**
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/ma_plot.png" width="400" height="400" align="middle" alt="DESeq2 MA plot"/>

**Sample to sample heatmap**
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/sample_to_sample.png" width="500" height="300" align="middle" alt="DESeq2 sample to sample heatmap"/>

**PCA plot**
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/PCA_plot.png" width="500" height="300" align="middle" alt="DESeq2 PCA plot"/>

**Normalized Boxplot**
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/boxplot.jpg" width="400" height="400" align="middle" alt="DESeq2 normalized boxplot"/>

**Cook Boxplot**
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/cooks_boxplot.png" width="400" height="400" align="middle" alt="DESeq2 Cook boxplot"/>

**Dispersion Estimate**
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/disp_calc.png" width="400" height="400" align="middle" alt="DESeq2 dispersion estimate"/>

**Pvalue test histogram**

<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/pvalue_hist.png" width="300" height="300" align="middle" alt="DESeq2 pvalue histogram"/>

**Top20 genes heatmap**
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/heatmap_top20.png" width="500" height="500" align="middle" alt="DESeq2 top20 heatmap"/>

**Hierarchical clustering**
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/hclust.png" width="700" height="500" align="middle" alt="DESeq2 hierarchical clustering"/>

**Standard desviation curve**
<img src="/processing_Data/bioinformatics/pipelines/rnaseq-nf/docs/images/standard_deviation.png" width="500" height="500" align="middle" alt="DESeq2 SD curve"/>

**Output directory:** `11-DESeq2/{comaprison_folder}/`
* `DE_matrix.txt`
  * Comparative table with the differential expression of two conditions.
* `clinical_data.txt`
  * Table with the metadata used to compare the samples.
* `maPlot_all.pdf`
  * MA plot of the DESeq analysis results for all the samples
* `heatmap_sample_to_sample.pdf`
  * Heatmap with the euclidean distance between samples.
* `plotPCA.pdf`
  * PCA plot of the samples for the rlog and the vsd.
      * rlog refers to the regularized log transformation, which transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size.
      * vsd refers to variance stabilizing transformation (VST), which calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) and then transforms the count data (normalized by division by the size factors or normalization factors), yielding a matrix of values which are now approximately homoskedastic (having constant variance along the range of mean values). The transformation also normalizes with respect to library size.
* `boxplot.pdf`
  * PDF file with the box_plots
      * Box plot of the normalized Counts
      * Box plot of the counts cook distances to see if one sample is consistently higher than others.
* `plotDispersions.pdf`
  * PDF file with plots to analyze the dispersion of the samples
      * Dispersion calc is the per-gene dispersion estimate together with the fitted mean-dispersion relationship.
      * Histogram with the test of the differential expression pvalues
* `cluster_dendrogram.pdf`
  * PDF file with the hierarchical clustering of the samples. The input data comes from the normalization of the counts. For the normalization DESeq uses the normalization of the ratios where the counts are divided by sample-specific size factors determined by median ratio of gene counts relative ro geometric mean per gene.
* `heatmapCount_top20.pdf`
  * Heat map of the top 20 genes with the higher normalized mean count. The normalization is the same that the one of the hierarchical clustering.
* `heatmapCount_all_genes.pdf`
  * Heatmap of the normalized counts of all the genes.
* `plotSD.pdf`
  * Standard deviation of the transformed data, across samples, against the mean, using the shifted logarithm transformation.--->

## References
1. <a name="fastqc_reference"></a> <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">web FASTQC</a> P. Babraham Bioinformatics - FastQC. A Quality Control tool for High Throughput Sequence Data. 2012
trimming_reference
* <a name="trimming_reference"></a> Bolger AM, Lohse M, Usadel B. <a href="https://www.ncbi.nlm.nih.gov/pubmed/24695404">Trimmomatic: a flexible trimmer for Illumina sequence data.</a> Bioinformatics. 2014 Apr 28;30(15):2114–20.
* <a name="STAR_reference"></a> Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. <a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=STAR%3A+ultrafast+universal+RNA-seq+aligner">STAR: ultrafast universal RNA-seq aligner.</a> Bioinformatics. 2013 Jan 1;29(1):15-21.
* <a name="RSeQC_reference"></a> Wang, L., Wang, S., Li, W. <a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=RSeQC%3A+quality+control+of+RNA-seq+experiments">RSeQC: quality control of RNA-seq experiments.</a> Bioinformatics. 2012 Aug 15;28(16):2184-5.
* <a name="preseq_reference"></a> <a href="http://smithlabresearch.org/software/preseq/">web Preseq</a> The Smith Lab. Computational Genomics Research.
* <a name="picard_reference"></a> <a href="https://github.com/broadinstitute/picard">web Picard toolkit</a> Broad Institute.
* <a name="dupradar_reference"></a> Sayols S, Scherzinger D, Klein H.<a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1276-2">dupRadar: a Bioconductor package for theassessment of PCR artifacts in RNA-Seqdata</a> BMC Bioinformatics (2016) 17:428.
* <a name="featurecounts_reference"></a> Liao Y, Smyth GK, Shi W. <a href="http://www.ncbi.nlm.nih.gov/pubmed/24227677">featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features.</a> Bioinformatics, 30(7):923-30, 2014
* <a name="StringTie_reference"></a> Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL. <a href="https://www.ncbi.nlm.nih.gov/pubmed/25690850">StringTie enables improved reconstruction of a transcriptome from RNA-seq reads.</a> Nat Biotechnol. 2015 Mar;33(3):290-5.
* <a name="edgeR_reference"></a> Robinson MD, McCarthy DJ, Smyth GK. <a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=edgeR%3A+a+Bioconductor+package+for+differential+expression+analysis+of+digital+gene+expression+data.">edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.</a> Bioinformatics. 2010 Jan 1;26(1):139-40.
* <a name="MultiQC_reference"></a> Ewels P, Magnusson M, Lundin S, Käller M. <a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=MultiQC%3A+Summarize+analysis+results+for+multiple+tools+and+samples+in+a+single+report">MultiQC: summarize analysis results for multiple tools and samples in a single report.</a> Bioinformatics. 2016 Oct 1;32(19):3047-8.
* <a name="gffcompare_reference"></a> Pertea G and Pertea M. <a href="https://doi.org/10.12688/f1000research.23297.1"> GFF Utilities: GffRead and GffCompare.</a> F1000Research 2020, 9:304.
* <a name="DESeq2_reference"></a> Love MI, Huber W, Anders S. <a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=Moderated+estimation+of+fold+change+and+dispersion+for+RNA-seq+data+with+DESeq2.">Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.</a> Genome Biol. 2014;15(12):550.
