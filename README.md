# Introduction
Repository with the markdowns of all of our analyses. 

## How to use
To convert the .md file into a html file you must run the following command:

```
Rscript markdowntohtml.r <imput_markdown.md> <output_report.html>
```

Once you have generated your html file, if you want to create a PDF file of the report you must install ['wkhtmltopdf'](https://wkhtmltopdf.org/) package. Then you should run the following command.

```
wkhtmltopdf --keep-relative-links <input_html.html> <output_pdf.pdf>
```
