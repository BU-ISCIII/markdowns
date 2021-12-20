# Output description for differ expression analysis

## Pipeline overview
* <a href="#deseq2">DESeq2</a> v1.18.1 - Differential expression analysis and plots

### DESeq2
Differential expression analysis with <a name="deseq2">DESeq2</a>. [DESeq2](https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf) <a href="#DESeq2_reference">[12]</a> is a Bioconductor package for R used for RNA-seq data analysis. The script included in the pipeline uses DESeq2 to normalize read counts and create a heatmap / dendrogram showing pairwise euclidean distance (sample similarity). It also creates other plots to evaluate the sample dispersion. It also provides PCA plots to evaluate sample grouping.

**MA plot**

<img src="/data/bi/pipelines/rnaseq-nf/docs/images/ma_plot.png" width="400" height="400" align="middle" alt="DESeq2 MA plot"/>

**Sample to sample heatmap**
<img src="/data/bi/pipelines/rnaseq-nf/docs/images/sample_to_sample.png" width="500" height="300" align="middle" alt="DESeq2 sample to sample heatmap"/>

**PCA plot**
<img src="/data/bi/pipelines/rnaseq-nf/docs/images/PCA_plot.png" width="500" height="300" align="middle" alt="DESeq2 PCA plot"/>

**Normalized Boxplot**
<img src="/data/bi/pipelines/rnaseq-nf/docs/images/boxplot.jpg" width="400" height="400" align="middle" alt="DESeq2 normalized boxplot"/>

**Cook Boxplot**
<img src="/data/bi/pipelines/rnaseq-nf/docs/images/cooks_boxplot.png" width="400" height="400" align="middle" alt="DESeq2 Cook boxplot"/>

**Dispersion Estimate**
<img src="/data/bi/pipelines/rnaseq-nf/docs/images/disp_calc.png" width="400" height="400" align="middle" alt="DESeq2 dispersion estimate"/>

**Pvalue test histogram**

<img src="/data/bi/pipelines/rnaseq-nf/docs/images/pvalue_hist.png" width="300" height="300" align="middle" alt="DESeq2 pvalue histogram"/>

**Top20 genes heatmap**
<img src="/data/bi/pipelines/rnaseq-nf/docs/images/heatmap_top20.png" width="500" height="500" align="middle" alt="DESeq2 top20 heatmap"/>

**Hierarchical clustering**
<img src="/data/bi/pipelines/rnaseq-nf/docs/images/hclust.png" width="700" height="500" align="middle" alt="DESeq2 hierarchical clustering"/>

**Standard desviation curve**
<img src="/data/bi/pipelines/rnaseq-nf/docs/images/standard_deviation.png" width="500" height="500" align="middle" alt="DESeq2 SD curve"/>


**Output files:**
* DE_matrix.txt
  * Comparative table with the differential expression of two conditions.
* maPlot.pdf
  * MA plot of the DESeq analysis results for all the samples
* heatmap_sample_to_sample.pdf
  * Heatmap with the euclidean distance between samples.
* plotPCA.pdf
  * PCA plot of the samples for the rlog and the vsd.
      * rlog refers to the regularized log transformation, which transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size.
      * vsd refers to variance stabilizing transformation (VST), which calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) and then transforms the count data (normalized by division by the size factors or normalization factors), yielding a matrix of values which are now approximately homoskedastic (having constant variance along the range of mean values). The transformation also normalizes with respect to library size.
* boxplot.pdf
  * PDF file with the box_plots
      * Box plot of the normalized Counts
      * Box plot of the counts cook distances to see if one sample is consistently higher than others.
* plotDispersions.pdf
  * PDF file with plots to analyze the dispersion of the samples
      * Dispersion calc is the per-gene dispersion estimate together with the fitted mean-dispersion relationship.
      * Histogram with the test of the differential expression pvalues
* hierarchical_clustering.pdf
  * PDF file with the hierarchical clustering of the samples. The input data comes from the normalization of the counts. For the normalization DESeq uses the normalization of the ratios where the counts are divided by sample-specific size factors determined by median ratio of gene counts relative ro geometric mean per gene.
* heatmapCount_top20.pdf
  * Heat map of the top 20 genes with the higher normalized mean count. The normalization is the same that the one of the hierarchical clustering.
* heatmapCount_all_genes.pdf
  * Heatmap of the normalized counts of all the genes.
* plotSD.pdf
  * Standard deviation of the transformed data, across samples, against the mean, using the shifted logarithm transformation.--->
