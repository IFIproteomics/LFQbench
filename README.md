---
references:
- id: fenner2012a
  title: One-click science marketing
  author:
  - family: Fenner
    given: Martin
  container-title: Nature Materials
  volume: 11
  URL: 'http://dx.doi.org/10.1038/nmat3283'
  DOI: 10.1038/nmat3283
  issue: 4
  publisher: Nature Publishing Group
  page: 261-263
  type: article-journal
  issued:
    year: 2012
    month: 3
---


[<img src="https://raw.githubusercontent.com/IFIproteomics/LFQbench/master/logo.png">](https://github.com/IFIproteomics/LFQbench)

[![Build Status](https://travis-ci.org/IFIproteomics/LFQbench.svg)](https://travis-ci.org/IFIproteomics/LFQbench) [![DOI](https://zenodo.org/badge/15862/IFIproteomics/LFQbench.svg)](https://zenodo.org/badge/latestdoi/15862/IFIproteomics/LFQbench)
======

An [R package](https://github.com/IFIproteomics/LFQbench) to  standardize  and analyze the results of SWATH data and specially the output of all software tools such as: PeakView, Spectronaut, OpenSWATH, DIA-Umpire. LFQbench evaluates and represents graphically precision and accuracy of label free quantification experiments based on hybrid samples, providing software developers with an standardized set of reports that enable an in-depth evaluation of their software performance. This is just a test [fenner2012a]

At the moment most of the library can be use using a set of shell scripts that generate a final report with the information of the data. The final report contains the following sections:


### Installation  

First, we need to install `devtools`:  

    install.packages("devtools")
    library(devtools)
   
Then we just call  

    install_github("IFIproteomics/LFQbench")
    library(LFQbench)

##Examples
=================

```{r, engine='bash', count_lines}

# Download the library code using git
> git clone https://github.com/IFIproteomics/LFQbench

# Move to the scripts folder
> cd inst/scripts

# run one of the datasets examples
> sh exec.example.peakview.sh

# In the main folder of the dataset: ext/data/example_peakview
> ../../ext/data/peakview

#Different folder shows the plot and information: (i) The log folder contain the log information, (ii) the input the example information, (iii) the plot folder the resulted plot fomr the analysis

```

### How to cite

