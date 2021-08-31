<p>
<img src="../assets/BU_ISCIII_logo.png" alt="logo" width="200" align="left"/>
</p>
<br>
<br>

## Pipeline description
Pipeline for bacterial assembly and annotation using nextflow.

## Pipeline overview

* [FastQC](#fastqc) v0.11.8 - read quality control.
* [Trimmomatic](#trimming) v.0.33 - adapter and low quality trimming.
* [Unicycler](#unicycler) v.0.4.6 - prokaryote assembler.
* [Quast](#quast) v4.1 - assemblies quality control
* [Kmerfinder](#kmerfinder) v.3.1 - species and contamination determination.
* [chewBBACA](#chewbbca) v.2.5.5 - cg/wgMLST analysis.
* [GrapeTree](#grapetree) v.2.2 - Allelic profile visualization.
* [Hamming distance](#hamming-distance) - Calculate the hamming distance among the alleles.
* [Ariba](#ariba) v.2.14.4 - Search resistance genes.
* [SRST2](#srst2) v.0.1.8 - the presence of STs and/or reference genes.


**Note:**

Depending on the analysis, we could have some ANALYSIS_IDs. This ANALYSIS_IDs are going to be composed of the date of the analysis, and an analysis identification. You can find a README in de ANALYSIS folder with a brief description of the different analysis.

## Preprocessing
### FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output directory: `03-assembly/fastqc/`**

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

**Results directory: `01-preprocessing`**
- Files:
   - `{sample_id}/{sample_id}_R[12]_filtered.fastq.gz`: contains high quality reads with both forward and reverse tags surviving.
   - `{sample_id}/{sample_id}_R[12]_unpaired.fastq.gz`: contains high quality reads with only forward or reverse tags surviving.

 **Note**:To see how your reads look after trimming, look at the FastQC reports in the 03-preprocQC directory

### Kmerfinder
[Kmerfinder](https://cge.cbs.dtu.dk/services/KmerFinder) (2) is a software used for species indentification and the determination of possible contamination in the sample. We use this software using the bacterial database provided by the developers, and with the "winner takes it all" algorithm. You can check [here](https://cge.cbs.dtu.dk/services/KmerFinder/output.php) for a description of the columns in the output.

**Output directory: `02-kmerfinder/{sample_id}/`**

* `data.json`
  * results in json format.
* `results.spa`
  * results in spa format.
* `results.txt`
  * results in txt format. **This is the format you have to use if you are going to open it with excel.**

**NOTE**: You can also find in `99-stats` a summary of all samples results (`kmerfinder.csv`).

## Assembly
### Unicycler
[Unicycler](https://github.com/rrwick/Unicycler) (3) is an assembly pipeline for bacterial genomes. It can assemble Illumina-only read sets where it functions as a SPAdes-optimiser.

**Output directory:** `03-assembly/unicycler/{sample_id}/`

* `{sample_id}.fasta`: fasta file containing the assembled reads in form of contigs and scaffolds. This is the file we use for annotation and upstream analysis.
* `{sample_id}.gfa`: Graph files for the different assembly optimizations. This files can be used by advanced users with software like [Bandage](https://github.com/rrwick/Bandage)


### QUAST
[QUAST](http://bioinf.spbau.ru/quast) (4) evaluates genome assemblies. We compared the reference genome with the contigs and scaffold assemblies. The html results can be opened with any browser (we recommend using Google Chrome).

**Output directory:** `03-assembly/quast/quast_results/`
* `quast_results/results_{date}/report.html`
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
    * L50 (L75, LG50, LG75) is the number of contigs equal to or longer than N50 (N75, NG50, NG75). In other words, L50, for example, is the minimal number of contigs that cover half the assembly.

### chewBBACA
[chewBBACA](https://github.com/B-UMMI/chewBBACA) (5) is a comprehensive pipeline for the creation and validation of whole genome and core genome MultiLocus Sequence Typing (wg/cgMLST) schemas, providing an allele calling algorithm based on BLAST Score Ratio that can be run in multiprocessor settings and a set of functions to visualize and evaluate allele variation in the loci. chewBBACA performs the schema creation and allele calls on complete or draft genomes resulting from de novo assemblers.
**Output directory:** `04-chewbbaca/allele_calling/results_<date>/`
* `results_contgsInfo.txt`: This file contains the loci positions in the genomes analyzed. This information can be useful to study locus synteny in complete genomes, for example.
* `results_alleles.tsv`: This file provides the MLST allelic profile output. You can calculate a Minimum spanning tree using this file, for example using [phyloviz](http://www.phyloviz.net/)
* `results_statistics.txt`: statitics about the status of the locus found or not found in each sample. More information about the different columns can be found [here](https://github.com/B-UMMI/chewBBACA/wiki/2.-Allele-Calling). Important ones are EXC (allele found as an exact match with previous identified alleles), INF (inferred new alleles using prodigal) and LNF (loci not found in our sample - meaning no blast hits within the bSR threshold for allele asignment).

### GrapeTree
[GrapeTree](https://enterobase.readthedocs.io/en/latest/grapetree/grapetree-about.html) (6) is a fully interactive, tree visualization program within EnteroBase, which supports facile manipulations of both tree layout and metadata. It generates GrapeTree figures using the neighbor-joining (NJ) algorithm, the classical minimal spanning tree algorithm (MSTree) similar to PhyloViz, or an improved minimal spanning tree algorithm which we call MSTree V2.
**Output directory:** `05-grapetree/`
* `metadata.txt`: Table with the samples and the metadata information for colouring.
* `ms_tree.nwk`: Newick tree of the MSTree.
* `MSTree.svg`: SVG figure if the MSTree.

### Hamming distance
Hamming distance is a metric for comparing two binary data strings. While comparing two binary strings of equal length, Hamming distance is the number of bit positions in which the two bits are different.
**Output directory:** `06-hamming_dist/`
* `dist.txt`: Hamming distance matrix.

### ARIBA
[ARIBA](https://github.com/sanger-pathogens/ariba/wiki) (7) was used to search for resistance genes and plasmids in different databases, for which it uses local assemblies.

**Output directory:** `07-ariba/`
* `run/{sample_id}`:
  * `out.card1043-17run`: Results of the Card database.
  * `out.megares1043-17run`: Results of the Megares database.
  * `out.plasmidfinder1043-17run`: Results of the plasmidfinder database.
  * `out.vfdb_full1043-17run`: Results of the VFDB database.
  For each folder:
    * `{sample_id}report.tsv`: report with detected genes (field explanation in [ANEXI](#anexi))
    * `assembled_genes.fa.gz`: sequence for the assembled resistance/virulence/etc genes.
* `summary/`:
  * `out.summarycard.csv`: Summary table of the Card database.
  * `out.summarymegares.csv`: Summary table of the Megares database.
  * `out.summaryplasmidfinder.csv`: Summary table of the plasmidfinder database.
  * `out.summaryvfdb_full.csv`: Summary table of the VFDB database.

### SRST2
[SRST2](https://github.com/katholt/srst2) (8) is designed to take Illumina sequence data, a MLST database and/or a database of gene sequences (e.g. resistance genes, virulence genes, etc) and report the presence of STs and/or reference genes.
**Output directory:** `08-srst2/`
* `<sample_name>/<outputprefix>__mlst__<db>__results.txt`: Tab delimited format table with the STs. Each locus has a column in which the best scoring allele number is printed.
* `<sample_name>/<outputprefix>__<sample>.<db>.sorted.bam`: Sorted bowtie2 alignment of reads to each input database.
* `<sample_name>/<outputprefix>__<sample>.<db>.pileup`: samtools pileup of the alignment.
* `<sample_name>/<outputprefix>.log`: SRST2 log file.
* `spneumoniae_mlst_report.txt__compiledResults.txt`: Table summarizing all the database results. Is a combination of the MLST style table plus the tabulated gene summary.

shigella_mlst_report.txt__compiledResults.txt

**Note**: For mor information about the results files, visit [their resulsts website](https://github.com/katholt/srst2#mlst-results)

## Final results.
We provide a folder with the most relevant results, in most of the cases these are the results you are interested in, and the ones you need to check and interpret.
**Output directory: `RESULTS`**
  * `results_alleles.tsv`: cgMLST allele matrix from chewBBACA.
  * `spneumoniae_MSTree.svg`: Minimum spanning tree image created using GrapeTree. This image shows the MST created from the allele calling using the cgMLST created using Chewbbaca.
  * `spneumoniae_ms_tree.nwk`: Minimum spanning tree image created using GrapeTree. This newick tree shows the MST created from the allele calling using the cgMLST created using Chewbbaca.
  * `kmerfinder.csv`: Kmerfinder summary table.
  * `allele_dist_matrix.txt`: shows number of allele differences among the strains. This table can be opened in excel.
  * `out.summaryvfdb_full.csv`: virulence genes detected in isolates using vfdb.
  * `out.summaryplasmidfinder.csv`: plasmids detected using INC genes plasmidFinder database.
  * `out.summarycard.csv`: resistances detected in the isolates using Card database.
  * `spneumoniae_mlst_report.txt__compiledResults.txt`: Table summarizing all the database results. Is a combination of the MLST style table plus the tabulated gene summary.

## Bibliography
1. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
2. Rapid and precise alignment of raw reads against redundant databases with KMA Philip T.L.C. Clausen, Frank M. Aarestrup, Ole Lund.
3. Wick RR, Judd LM, Gorrie CL, Holt KE. Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 2017.
4. Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi and Glenn Tesler. QUAST: quality assessment tool for genome assemblies, Bioinformatics (2013) 29 (8): 1072-1075.
5. Silva M, Machado M, Silva D, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço J. 15/03/2018. M Gen 4(3): doi:10.1099/mgen.0.000166
6. Z Zhou, NF Alikhan, MJ Sergeant, N Luhmann, C Vaz, AP Francisco, JA Carrico, M Achtman (2018) “GrapeTree: Visualization of core genomic relationships among 100,000 bacterial pathogens”, Genome Res.
7. Martin Hunt, Alison E. Mather, Leonor Sánchez-Bus, Andrew J. Page, Julian Parkhill, Jacqueline A Keane, Simon R. Harris. ARIBA: rapid antimicrobial resistance genotyping directly from sequencing reads.  Microb Genom. 2017 Sep 4;3(10):e000131
8. Inouye, M., Dashnow, H., Raven, LA. et al. SRST2: Rapid genomic surveillance for public health and hospital microbiology labs. Genome Med 6, 90 (2014).


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
