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
* [PlasmidID](#plasmidid) v.1.4.1 - Assembly report and visualization.
* [ARIBA](#ariba) v.2.14.4 - Search resistance genes.
* [Chewbbaca](#chewbbaca) v.2.0.17.2 - cg/wgMLST construction and analysis software.

**Note:**

Depending on the analysis, we could have some ANALYSIS_IDs. This ANALYSIS_IDs are going to be composed of the date of the analysis, and an analysis identification. You can find a README in de ANALYSIS folder with a brief description of the different analysis.

## Preprocessing
### FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output directory: `03-assembly/fastqc`**

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

**Results directory: `03-assembly/trimming`**
- Files:
   - `{sample_id}/{sample_id}_R[12]_filtered.fastq.gz`: contains high quality reads with both forward and reverse tags surviving.
   - `{sample_id}/{sample_id}_R[12]_unpaired.fastq.gz`: contains high quality reads with only forward or reverse tags surviving.

 **Note**:To see how your reads look after trimming, look at the FastQC reports in the 03-preprocQC director

## Assembly
### Unicycler
[Unicycler](https://github.com/rrwick/Unicycler) (3) is an assembly pipeline for bacterial genomes. It can assemble Illumina-only read sets where it functions as a SPAdes-optimiser.

**Output directory:** `03-assembly/unicycler/{sample_id}`

* `{sample_id}.fasta`: fasta file containing the assembled reads in form of contigs and scaffolds. This is the file we use for annotation and upstream analysis.
* `{sample_id}.gfa`: Graph files for the different assembly optimizations. This files can be used by advanced users with software like [Bandage](https://github.com/rrwick/Bandage)


### QUAST
[QUAST](http://bioinf.spbau.ru/quast) (5) evaluates genome assemblies. We compared the reference genome with the contigs and scaffold assemblies. The html results can be opened with any browser (we recommend using Google Chrome).

**Output directory:** `03-assembly/quast/quast_results`
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

### PlasmidID

[PlasmidID](https://github.com/BU-ISCIII/plasmidID) was used to graphically represent the alignment of the reference genome relative to a given assembly. This helps to visualize the coverage of the reference genome in the assembly. To find more information about the output files refer to the [documentation](https://github.com/BU-ISCIII/plasmidID/wiki/Understanding-the-image:-track-by-track).

**Output directory:** `05-plasmidID/NO_GROUP/`

* `00_summary/FINAL_REPORT_NO_GROUP.tsv`: .tsv table with the  summary of the plasmids found in each sample and the times one plasmid appeared among all the samples.
* `{sample_id}/images/{sample_id}_{plasmid_id}.png`: PNG file with the visualization of the alignment between the viral assembly and the reference viral genome.
* `{sample_id}/data/`: Files used for drawing the circos images.
* `{sample_id}/database/`: Annotation files used for drawing the circos images.
* `{sample_id}/fasta_files`: Folder with fasta files that correspond to the selection of contigs/scaffolds required to reconstruct the reference genome generated in the `images/` folder.
* `{sample_id}/log/`: Log files.


### ARIBA

[ARIBA](https://github.com/sanger-pathogens/ariba/wiki) (6) was used to search for resistance genes and plasmids in different databases, for which it uses local assemblies.

**Output directory:** `06-ariba/`

* `run/{sample_id}`:
  * `out.card1043-17run`: Results of the Card database.
  * `out.megares1043-17run`: Results of the Megares database.
  * `out.plasmidfinder1043-17run`: Results of the plasmidfinder database.
  * `out.vfdb_full1043-17run`: Results of the VFDB database.
* `summary/`:
  * `out.summarycard.csv`: Summary table of the Card database.
  * `out.summarymegares.csv`: Summary table of the Megares database.
  * `out.summaryplasmidfinder.csv`: Summary table of the plasmidfinder database.
  * `out.summaryvfdb_full.csv`: Summary table of the VFDB database.

## Phylogenetic analysis
### Chewbbaca
[chewBBACA](https://github.com/B-UMMI/chewBBACA/wiki) (Silva et al, 2018) is a comprehensive pipeline for the creation and validation of whole genome and core genome MultiLocus Sequence Typing (wg/cgMLST) schemas, providing an allele calling algorithm based on BLAST Score Ratio that can be run in multiprocessor settings and a set of functions to visualize and evaluate allele variation in the loci. chewBBACA performs the schema creation and allele calls on complete or draft genomes resulting from de novo assemblers.

**Output directory:** `03-chewbbaca`

* `cgmlst_extraction`
    * `Presence_Abscence.tsv` - allele presence and absence matrix (1 or 0, respectively) for all the loci found in the `-i` file (including those loci and genomes that were flagged to be excluded).
    * `cgMLST.tsv` - matrix with the cgMLST allelic profile (already excluding the identified loci and genomes).
    * `cgMLSTschema.txt`- list of loci that constitute the redefined cgMLST schema. To use in the _Allele Calling_ engine, the full path to the loci FASTA must be added.
  * `allele_calling_cgmlst/results_*/`
    * `results_contgsInfo.txt`: This file contains the loci positions in the genomes analyzed. This information can be useful to study locus synteny in complete genomes, for example.
    * `results_alleles.tsv`: This file provides the MLST allelic profile output. You can calculate a Minimum spanning tree using this file, for example using [phyloviz](http://www.phyloviz.net/)
    * `results_statistics.txt`: statitics about the status of the locus found or not found in each sample. More information about the different columns can be found [here](https://github.com/B-UMMI/chewBBACA/wiki/2.-Allele-Calling). Important ones are EXC (allele found as an exact match with previous identified alleles), INF (inferred new alleles using prodigal) and LNF (loci not found in our sample - meaning no blast hits within the bSR threshold for allele asignment).

## Final results.
We provide a folder with the most relevant results, in most of the cases these are the results you are interested in, and the ones you need to check and interpret.
**Output directory: `RESULTS`**
  * `bfragilis_cgmlst_alleles.tsv`: cgMLST allele matrix.
  * `bfragilis_cgmlst_mst.svg`: Minimum spanning tree image created using GrapeTree. This image shows the MST created from the allele calling using the cgMLST created using Chewbbaca.
  * `bfragilis_cgmlst_mst.nwk`: Minimum spanning tree image created using GrapeTree. This newick tree shows the MST created from the allele calling using the cgMLST created using Chewbbaca.
  * `bfragilis_cgmlst_mst.json`: Minimum spanning tree image created using GrapeTree. This files shows the MST created from the allele calling using the cgMLST created using Chewbbaca, it can be used with GrapeTree software.
  * `bfragilis_hammingDistance.tsv`: shows number of allele differences among the strains.
  * `bfragilis_cardresistances.csv`: resistances detected in the isolates using Card database.
  * `bfragilis_vfdb.csv`: virulence genes detected in isolates using vfdb.
  * `bfragilis_plasmidfinder.csv`: plasmids detected using INC genes plasmidFinder database.
  * `{sample_id}_{plasmidref}_individual.png`: circos visualization of plasmid found in isolates.

## Bibliography
S1. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
2. Rapid and precise alignment of raw reads against redundant databases with KMA Philip T.L.C. Clausen, Frank M. Aarestrup, Ole Lund.
3. Wick RR, Judd LM, Gorrie CL, Holt KE. Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 2017.
4. Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics 2014 Jul 15;30(14):2068-9.
5. Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi and Glenn Tesler. QUAST: quality assessment tool for genome assemblies, Bioinformatics (2013) 29 (8): 1072-1075.
6. Martin Hunt, Alison E. Mather, Leonor Sánchez-Bus, Andrew J. Page, Julian Parkhill, Jacqueline A Keane, Simon R. Harris. ARIBA: rapid antimicrobial resistance genotyping directly from sequencing reads.  Microb Genom. 2017 Sep 4;3(10):e000131
