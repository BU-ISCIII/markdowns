# Output description for SRVIIER242 analysis.

## Pipeline overview

* [Samtools](#samtools) v1.9 - mpileup generation.
* [Custom parsing](#custom-parsing) - mpileup parsing.


**Note:**

Depending on the analysis, we could have some ANALYSIS_IDs. This ANALYSIS_IDs are going to be composed of the date of the analysis, and an analysis identification. You can find a README in de ANALYSIS folder with a brief description of the different analysis.

## Variant Calling
### Samtools mpileup
[Samtools mpileup](http://www.htslib.org/) is used for generating a pileup file for all he positions provided in the service.

**Output directory: `01-mpileup`**

* `{sample_id}.pileup`
  * file in pileup format with information for each position.

## Custom parsing

First we parse the pileup file into a tab format, then we merge all results in only one file and generate some statistics with a custom R script.

**Output directory: `02-mpileup2tab`**

* `{sample_id}_counts_var.tab`
  * pileup file parsed into tab format for each sample.

**Output directory: `03-pampuCaller`**

* `all_samples_variants_counts.tab`
  * calculations of frequency and stats for all samples.
  
  ## Final results
  Final result file can be found in `RESULTS` folder.
