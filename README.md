# Markdowns

![BU-ISCIII/markdowns](BU_ISCIII_logo.png =100x20)

## Introduction
Repository with the list of markdowns developed by our group for the different analyses:

* [RNA-seq nextflow](https://github.com/BU-ISCIII/rnaseq-nf/blob/master/docs/output.md)
* [Bacterial Assembly nextflow](https://github.com/BU-ISCIII/bacterial_assembly-nf/blob/master/docs/output.md)
* [Low Frequency Variants nextflow](https://github.com/BU-ISCIII/panelLowFreq-nf/blob/master/doc/output.md)
* [Viral Genome Assembly](https://github.com/BU-ISCIII/markdowns/output_virus.md)


## How to use
To convert the .md file into a html file you must run the following command:

```
Rscript markdowntohtml.r <imput_markdown.md> <output_report.html>
```

Once you have generated your html file, if you want to create a PDF file of the report you must install ['wkhtmltopdf'](https://wkhtmltopdf.org/) package. Then you should run the following command.

```
wkhtmltopdf --keep-relative-links <input_html.html> <output_pdf.pdf>
```
