# Output description for viral genome pipeline

## Pipeline overview

* [FastQC](#fastqc) v0.11.8 - read quality control.
* [Trimmomatic](#trimming) v.0.33 - adapter and low quality trimming.
* [BWA](#bwa) v0.7.12 - mapping against reference genome.
* [SAMtools](#samtools) v1.2 - Alignment result processing and unmapped reads selection.
* [Picard](#picard) v1.140 - Enrichment and alignment metrics.
* [Bcftools](#bcftools) v1.9 - Variant calling and consensus genome
* [SPADES](#spades) v3.7.1 - Viral genome assembly.
<!---** [QUAST](#quast) v4.1 - Assembly quality assessment.--->
* [PlasmidID](#plasmidid) v.1.4.1 - Visualization of the alignment.
<!---*** [BLAST](#blast) v2.10.0 - Local alignment.--->

Depending on the analysis, we will have some ANALYSIS_IDs. This ANALYSIS_IDs are going to be composed of the date of the analysis, the viral genome and the host genome (for example 20191112_GENOMEDENGUE_BUG and 20191112_GENOMEDENGUE_HUMAN) for the analysis of Dengue samples with Bug host and Human host. Inside these analysis you will find the folder corresponding to 04-mapping_host, 05_mapping_virus, 06-mapping_consensus, 07-assembly, 08-quast and 10-final_results, which are specific analysis for each host.

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
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used for removal of adapter contamination and trimming of low quality regions.
Parameters included for trimming are:
-  Nucleotides with phred quality < 10 in 3'end.
-  Mean phred quality < 20 in a 4 nucleotide window.
-  Read lenght < 50

**Results directory: `02-preprocessing`**
- Files:
   - `{sample_id}/{sample_id}_R[12]_filtered.fastq.gz`: contains high quality reads with both forward and reverse tags surviving.
   - `{sample_id}/{sample_id}_R[12]_unpaired.fastq.gz`: contains high quality reads with only forward or reverse tags surviving.

 **Note**:To see how your reads look after trimming, look at the FastQC reports in the ANALYSIS/{ANALYSIS_ID}/03-preprocQC directory

 **Note**:From now on, all the steps will be host specific.

## Mapping
### BWA
[BWA](http://bio-bwa.sourceforge.net/) or Burrows-Wheeler Aligner, is designed for mapping low-divergent sequence reads against reference genomes. The result alignment files are further processed with [SAMtools](http://samtools.sourceforge.net/), sam format is converted to bam, sorted and an index .bai is generated.

We mapped the fastq file agains both reference host genome and reference viral genome.

**Output directory: `0[4/5]-mapping_[host/virus]`**

* `{sample_id}_sorted.bam`
  * Sorted aligned bam file.
* `{sample_id}_sorted.bam.bai`
  * Index file for soreted aligned bam.
* `{sample_id}_flagstat.txt`
  * Mapping stats summary.

### Picard
[Picard](https://broadinstitute.github.io/picard/index.html) is a set of command line tools for manipulating high-throughput sequencing (HTS) data. In this case we used it to obtain mapping stats.

**Output directory: `0[4/5]-mapping_[host/virus]`**

* `{sample_id}.stats`
  * Picard metrics summary file for evaluating coverage and performance.

Picard documentation: [Picarddocs](https://broadinstitute.github.io/picard/command-line-overview.html)

## Variant calling and consensus genome
The first approach we used to generate the consensus viral genome was to call for variants between the mapped reads and the reference viral genome, and adding these variants to the reference viral genome.
## Bcftools
[Bcftools](http://samtools.github.io/bcftools/bcftools.html) mpileup command is used for generate a pileup for one the BAM files. In the pileup format each line represents a genomic position, consisting of chromosome name, 1-based coordinate, reference base, the number of reads covering the site, read bases, base qualities and alignment mapping qualities. Information on match, mismatch, indel, strand, mapping quality and start and end of a read are all encoded at the read base column.

The resulting variant calling vcf for haploid genomes is indexed and then the consensus genome is created adding the variants to the reference viral genome. This consensus genome was obtained using the predominant variants of the mapping file.

**Output directory: `06-mapping_consensus`**
* `{sample_id}_{reference_virus_name}/{sample_id}_{reference_virus_name}.vcf.gz`
  * Compressed format of the variants file.
* `{sample_id}_{reference_virus_name}/{sample_id}_{reference_virus_name}.vcf.gz.csi`
  * Compressed format of the indexed variants file.
* `{sample_id}_{reference_virus_name}/{sample_id}_{reference_virus_name}_consensus.fasta`
  * Consensus viral genome file generated from adding the variants called before to the viral reference genome. These variants are only the majoritarian variants, inlcuding only SNPs and small indels. This file is also contained in the 10-final_results folder as {sample_id}_{reference_virus_name}_consensus.fasta. This file is also contained in the ../RESULTS folder with the main results.

## Viral genome assembly
Other approach we used to generate the consensus viral genome was assembling the reads that didn't mapped to the host genome.
### SAMtools
In this section SAMtools was used to obtain the reads that didn't mapped with the host genome. This reads where sorted and converted to fastq files.

**Output directory: `07-assembly`**
* `{sample_id}_R[1,2]_unmapped.fastq`
  * fastq file with the reads that didn't mapped with the host genome.


### SPAdes
[SPAdes](https://kbase.us/applist/apps/kb_SPAdes/run_SPAdes/release?gclid=Cj0KCQiAt_PuBRDcARIsAMNlBdroQS7y2hPFuhagq1QPvQ39FcvGxbhtZwhn8YbxIB4LrGIHKjJ-iPwaAn_lEALw_wcB) is a de Bruijn graph-based assembler. We selected the reads that didn't mapped with the host genome and assembled them using SPAdes to create a viral genome assembly.

**Output directory: `07-assembly`**
* `{sample_id}/contigs.fasta`
  * Assembled contigs. This file is also contained in the 10-final_results folder as {sample_id}_consensus.fasta.
* `{sample_id}/scaffolds.fasta`
  * Assembled scaffolds.

<!---*### QUAST
[QUAST](http://bioinf.spbau.ru/quast) evaluates genome assemblies. We compared the reference genome with the contigs and scaffold assemblies. The html results can be opened with any browser (we recommend using Google Chrome).

**Output directory: `08-quast`**
* `quast_results/date/report.html`
  * Compressed format of the indexed variants file.
  * The meaning of the different metrics:
    * Contigs (≥ x bp): is total number of contigs of length ≥ x bp.
    * Total length (≥ x bp): is the total number of bases in contigs of length ≥ x bp.
    * Contigs: is the total number of contigs in the assembly.
    * Largest contig: is the length of the longest contig in the assembly.
    * Total length: is the total number of bases in the assembly.
    * Reference length: is the total number of bases in the reference genome.
    * GC (%): is the total number of G and C nucleotides in the assembly, divided by the total length of the assembly.
    * Reference GC (%): is the percentage of G and C nucleotides in the reference genome.
    * N50: is the length for which the collection of all contigs of that length or longer covers at least half an assembly.
    * NG50: is the length for which the collection of all contigs of that length or longer covers at least half the reference genome. This metric is computed only if the reference genome is provided.
    * N75 and NG75: are defined similarly to N50 but with 75 % instead of 50 %.
    * L50 (L75, LG50, LG75) is the number of contigs equal to or longer than N50 (N75, NG50, NG75). In other words, L50, for example, is the minimal number of contigs that cover half the assembly.--->

### PlasmidID
[PlasmidID](https://github.com/BU-ISCIII/plasmidID) was used to graphically represent the alignment of the reference genome with the assembly obtained with SPAdes. This helps to visualize the coverage of the reference genome in the assembly. To find more information about the output files go to: https://github.com/BU-ISCIII/plasmidID/wiki/Understanding-the-image:-track-by-track

**Output directory: `../RESULTS/`**
* `{sample_id}_{reference_viral_genome}.png`
  * PNG file with the visualization of the alignment between the assembled viral genome and the reference viral genome.
* `{sample_id}_consensus.fasta`
  * Fasta filte with the contigs that reconstruct the reference genome.

## Custom Analysis
### BLAST
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) finds regions of similarity between biological sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance. We performed a nucleotide BLAST (blastn) through the web server, against the Nucleotide Collection database. The results are not shown. The analysis was performed both using a viral genomes database (10-blast_virus) and a bacterial genomes database (11-blast_bacteria).
**Output directory: `10-blast_virus/ & 11-blast_bacteria`**
* `{sample_id}_blast_filt_header.txt`: txt file with the parsed BLAST results. It contains the following columns:
  * stitle: Subject title corresponds to the name of the genome that has a BLAST hit.
  * qaccver: Corresponds to the sequence on the query (sample) fasta file that has a hit with stitle.
  * saccver: Corresponds to the accession of stitle in the NCBI database.
  * pident: Percentage of identity between the subject (reference) and the query (sample).
  * length: length of the alignment.
  * mismatch: # of mismatches in the alignment.
  * gapopen: # of gaps that have been opened.
  * qstart: nucleotide where the alignment starts in the query.
  * qend: nucleotide where the alignment ends in the query.
  * sstart: nucleotide where the alignment starts in the subject.
  * send: nucleotide where the alignment ends in the subject.
  * evalue: is the number of expected hits of similar quality (score) that could be found just by chance.
  * bitscore: measures sequence similarity independent of query sequence length and database size and is normalized based on the rawpairwise alignment score.
  * slen: Subject nucleotide length.
  * qlen: Query nucleotide length.
  * qcovs: % of query coverage.
  * %cgAligned: % of contig aligned to the reference.
  * %refCovered: %  of the reference that is covered.
