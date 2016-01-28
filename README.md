[<img src="https://raw.githubusercontent.com/IFIproteomics/LFQbench/master/logo.png">](https://github.com/IFIproteomics/LFQbench)

[![Build Status](https://travis-ci.org/IFIproteomics/LFQbench.svg)](https://travis-ci.org/IFIproteomics/LFQbench) [![DOI](https://zenodo.org/badge/15862/IFIproteomics/LFQbench.svg)](https://zenodo.org/badge/latestdoi/15862/IFIproteomics/LFQbench)
======

### Description

LFQbench is an open source [R package](https://github.com/IFIproteomics/LFQbench) for the automated evaluation of label-free quantification performance. The evaluation bases on the interpretation of the quantitative analysis results of hybrid proteome samples prepared in known ratios<sup>[1]</sup>.

LFQbench calculates and represents graphically a set of qualitative and quantitative performance metrics like identification rates,  precision and accuracy of quantification, providing developers and end-users with a standardized set of reports to enable an in-depth performance evaluation of their software and analysis platforms.

### Installation  

First, we need to install `devtools`:  

    install.packages("devtools")
    library(devtools)
   
Then we just call  

    install_github("IFIproteomics/LFQbench")
    library(LFQbench)

### Examples

```{r, engine='bash', count_lines}
# Download the library code using git
> git clone https://github.com/IFIproteomics/LFQbench

# Move to the scripts folder
> cd inst/scripts

# run one of the datasets examples
> sh exec.example.peakview.sh

# In the main folder of the dataset: ext/data/example_peakview
> ../../ext/data/peakview

# Different folder shows the plot and information: 
# (i) the log folder contain the log information, 
# (ii) the input the example information, 
# (iii) the plot folder the resulted plot from the analysis
```

### References

[1] Kuharev, J., Navarro, P., Distler, U., Jahn, O. & Tenzer, S. In-depth evaluation of software tools for data-independent acquisition based label-free quantification. Proteomics 15, 3140–3151 (2015).

