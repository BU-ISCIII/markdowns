<p>
<img src="../assets/BU_ISCIII_logo.png" alt="logo" width="200" align="left"/>
</p>
<br>
<br>

## Pipeline description
Pipeline for gene panel variant calling and annotation.

## Pipeline overview

* [FastQC](#fastqc) v0.11.8 - read quality control.
* [Trimmomatic](#trimming) v0.33 - adapter and low quality trimming.
* [BWA](#bwa) v0.7.12 -  mapping against reference genome.
* [SAMtools](#samtools) v1.2 - alignment result processing
* [VarScan2](#varscan2) v2.3.9 - variant calling
* [KGGSeq](#kggseq) v1.2 - variants annotation
* [Picard](#picard) v1.140 - Mapping statistics
* [MultiQC](#multiqc) v1.9 - Present QC for raw reads, alignment, assembly and variant calling

## Preprocessing
### FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output directory: `01-fastQC`**

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

#### Trimming FastQC
After the trimming steps, FastQC for new quality control is performed.

**Output directory: `03-preprocQC`**

* `{sample_id}/{sample_id}_R[12]_filtered_fastqc.html`
  * html report. This file can be opened in your favourite web browser (Firefox/chrome preferable) and it contains the different graphs that fastqc calculates for QC.
* `{sample_id}/{sample_id}_R[12]_filtered_fastqc`
  * older with fastqc output in plain text.
* `{sample_id}/{sample_id}_R[12]_filtered_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## Mapping
### BWA
[BWA](http://bio-bwa.sourceforge.net/) or Burrows-Wheeler Aligner, is designed for mapping low-divergent sequence reads against reference genomes. The result alignment files are further processed with [SAMtools](http://samtools.sourceforge.net/), sam format is converted to bam, sorted and an index .bai is generated.

**Output directory: `04-mapping`**

* `{sample_id}_sorted.bam`
  * Sorted aligned bam file.
* `{sample_id}_sorted.bam.bai`
  * Index file for soreted aligned bam.

### SAMtools
[Samtools](http://www.htslib.org/) is a suite of programs for interacting with high-throughput sequencing data. We use the [mpiluep](http://www.htslib.org/doc/samtools-mpileup.html) module to generate the pileup files from the mapping process, which produces "pileup" textual format from an alignment.

**Output directory: `05-samtools`**

* `{sample_id}/{sample_id}.pileup`
  * Pileup file

## Variant calling and Annotation
### VarScan2
[VarScan 2](http://dkoboldt.github.io/varscan/) is a platform-independent software tool to detect variants in NGS data. In this pipeline, VarScan 2 is used in conjunction with SAMtools in order to call both high and low frequency variants.

**Output directory: `06-VarScan`**

* `{sample_id}/{sample_id}.vcf`
  * Variants file

### KGGSeq
[KGGSeq](http://pmglab.top/kggseq/index.htm) is a software platform constituted of Bioinformatics and statistical genetics functions making use of valuable biologic resources and knowledge for sequencing-based genetic mapping of variants/genes responsible for human diseases/traits.

**Output directory: `05-samtools`**

* `{sample_id}_all_annotated.tab`
  * Table with all the variants annotated.
* `{sample_id}_RB1_annotated.tab`
  * Table with RB1 gene the variants annotated.

## Stats
### Picard
[Picard](https://broadinstitute.github.io/picard/command-line-overview.html) is a set of command-line tools for manipulating high-throughput sequencing data. We use picard-tools in this pipeline to obtain mapping and coverage metrics.

**Output directory: `99-stats/picard`**

* `{sample_id}_hsMetrics.out`
  * Mapping metrics for each sample.
* `hsMetrics_all.out`
  * Metrics for all the samples.

### MultiQC
[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

**Output directory: `99-stats/`**
* `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
* `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
