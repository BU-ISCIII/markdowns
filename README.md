<img src="./BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

# Markdowns
## Introduction
Repository with the list of markdowns developed by our group for the different analyses:

* [RNA-seq nextflow](https://github.com/BU-ISCIII/rnaseq-nf/blob/master/docs/output.md)
* [Bacterial Assembly nextflow](https://github.com/BU-ISCIII/bacterial_assembly-nf/blob/master/docs/output.md)
* [Low Frequency Variants nextflow](https://github.com/BU-ISCIII/panelLowFreq-nf/blob/master/doc/output.md)
* [Viral Genome Assembly](https://github.com/BU-ISCIII/markdowns/blob/master/reports/output_virus.md)
* [Bacterial assembly + snippy + canSNP](https://github.com/BU-ISCIII/markdowns/blob/master/reports/output_snippy_assembly_cansnp.md)
* [Exome pipeline](https://github.com/BU-ISCIII/exome_pipeline/blob/develop/doc/output.md)


## How to use
To convert the .md file into a html file you must run the following command:

```
Rscript markdown_to_html.r <imput_markdown.md> <output_report.html>
```

Once you have generated your html file, if you want to create a PDF file of the report you must install ['wkhtmltopdf'](https://wkhtmltopdf.org/) package. Then you should run the following command.

```
wkhtmltopdf --keep-relative-links <input_html.html> <output_pdf.pdf>
```
