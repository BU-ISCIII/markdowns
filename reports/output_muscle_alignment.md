<p>
<img src="../assets/BU_ISCIII_logo.png" alt="logo" width="200" align="left"/>
</p>
<br>
<br>

## Pipeline description
Pipeline for multifasta alignment using MUSCLE.

## Pipeline overview

* [MUSCLE](#muscle) v3.8.1551 - Multifasta sequence alignment

**Note:**

Depending on the analysis, we could have some ANALYSIS_IDs. This ANALYSIS_IDs are going to be composed of the date of the analysis, and an analysis identification. You can find a README in de ANALYSIS folder with a brief description of the different analyses.

## Alignment
### MUSCLE
[MUSCLE](http://www.drive5.com/muscle/) MUSCLE is one of the best-performing multiple alignment programs according to published benchmark tests, with accuracy and speed that are consistently better than CLUSTALW. MUSCLE can align hundreds of sequences in seconds. Most users learn everything they need to know about MUSCLE in a few minutes—only a handful of command-line options are needed to perform common alignment tasks.

**Files:**

* `{file_name}_aligned.fasta`
  * Aligned FASTA file that can be visualized using any alignment visualizing program such as JalView or MegaX.

## References
1. Edgar R.C. MUSCLE: a multiple sequence alignment method with reduced time and space complexity. (2004) BMC Bioinformatics. 2004 Aug 19;5:113. doi: 10.1186/1471-2105-5-113.
