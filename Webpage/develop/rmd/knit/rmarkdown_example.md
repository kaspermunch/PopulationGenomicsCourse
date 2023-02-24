---
title: Rmarkdown example
summary: This is an example of an Rmarkdown development
date: 2000-01-01
knit: (function(inputFile, encoding) { 
      rmarkdown::render(inputFile,
                        encoding=encoding,
                        output_format='all',
                        output_dir='./knit/')})
output:
  github_document: 
     preserve_yaml: TRUE
     html_preview: FALSE
     pandoc_args: [
      "--wrap", "none" # this is needed to not break admonitions
    ]
---

Rmarkdown example
================
2000-01-01

# Rmarkdown example

!!! note “Section Overview”

    &#128368; **Time Estimation:** X minutes  

    &#128172; **Learning Objectives:**    
        1. First item  
        2. Second item  

This is an example of how to make md pages from rmarkdown documents. The source code to generate this md document is in `develop/rmd/rmarkdown_example.Rmd`.

First we need to write the document as usual. Try to use the admonitions as in previous examples.

## Writing code

You can run code as usual in Rmarkdown. The chunk options can be very useful. By default, the code will not run unless explicity specified in the chunk (you can change this behaviour on the “knitr” chunk above).

``` r
library(tidyverse)
```

``` r
ggplot(data = txhousing, aes(x = sales, y = listings)) +
  geom_point()
```

<img src="./images/rmarkdown_example/ggplot-1.png" style="display: block; margin: auto;" />

## Adding images

You can add images using the `knitr::include_graphics()` option. The images should be in the same folder, otherwise they will break in the step below:

<img src="./images/favicon.png" width="853" style="display: block; margin: auto;" />

## Knitting the document and formatting item

To use this document in the webpage, you will have to knit it first as a github_document. This is already setup in the yaml menu above. The knitted files will be created in the `knit` folder. Unfortunately, knitting creates a couple of issues. First, we need to get rid of the knit options of the yaml file. Knitting will also create a header with the title, author and date on the newly created md document. Lastly, for some reason, knitting will change the normal quotes “” into cursive quotes “”. This will completely mess up the admonitions, and mkdocs will not display them properly.

In order to solve this issues, there is a small `modify_md.ipynb` file that will take all the created github_document (\*.md) files and fix these issues.

``` python
from os import listdir
import re

files = listdir("./knit/")

print(files)
```

``` python
for i in files:
    if i.endswith(".md"):
        with open("./knit/"+ i) as f:
            text = f.read()
            text = re.sub("\nknit((.|\n)*)\n[0-9]{4}-[0-9]{2}-[0-9]{2}\n", "\n---\n" ,text)
            text = re.sub("(”|“)","\"",text)
        with open("./knit/"+ i, "w") as f:
            f.write(text)
```

You can then move the knitted documents (and the images if needed) to the `develop` folder.
