# DNEA

*DNEA* is an R package to construct data-driven biological networks from -omics data. The package was specifically built for metabolomics and lipidomics data, but can be used on any approximately normal data sets. The package takes as input the raw peak intensity / concentrations for a set of features from two experimental conditions and uses the GLASSO algorithm to estimate metabolite-metabolite partial correlations.

## Installation

```r
if (!require("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("DNEA")
```

## Citation
The DNEA R package is accompanied by a peer-reviewed manuscript published in ***BMC Bioinformatics***. If using this software in your work, please cite:

Patsalis C, Iyer G, Brandenburg M, Karnovsky A, Michailidis G. DNEA: an R package for fast and versatile data-driven network analysis of metabolomics data. BMC Bioinformatics 25, 383 (2024). DOI:10.1186/s12859-024-05994-1.

## Motivation for DNEA

 Advancements in analytical methods, such as Liquid Chromatography-Mass Spectrometry (LC-MS), have enabled the high-throughput identification of hundreds to thousands of metabolites in biological samples, creating larger, more complex data sets that need to be analyzed. Pathway enrichment analysis is commonly used to identify the metabolic mechanism(s) underlying a disease state. However, there are several challenges present in analyzing metabolomics and lipidomics data sets using traditional bioinformatics tools. These conventional methods rely on well annotated pathway databases, which has proven to be a difficult task in metabolomics and lipidomics due to de novo identification of the metabolome. Further, metabolomics studies often identify exogenous compounds that cannot be mapped to a human pathway. 
  
## How it works

  To combat this, our group developed ***The Differential Network Expression Analysis (DNEA)*** algorithm outlined in [Ma et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748777/) and later implemented in the Filigree java-application presented in [Iyer et al 2021](https://pubmed.ncbi.nlm.nih.gov/33255384/). DNEA is a data-driven network analysis tool for biologial data that utilizes gaussian graphical models (GGM) to jointly estimate the biological networks of two conditions. Data-driven approaches to network analysis of metabolomics and lipidomics data has been difficult due to an abundance of features in comparison to the number of samples in the average -omics study, also known as the ***p >> n*** problem. R1 regularization, in combination with a novel method for stability selection, allow us to overcome this challenge by only keeping important metabolites in the model. This approach provides more robust results and enables the identification of true metabolite-metabolite interactions otherwise lost by reference-based methods. Networks constructed using DNEA can be clustered to identify metabolic modules, or sub networks, within the data and subsequently test them for enrichment across experimental conditions using the [netgsa R package](https://cran.rstudio.com/web/packages/netgsa/index.html). By using a data-driven approach to define the metabolic pathways, we improve the accuracy of pathway analysis.
  
## Improvements in this implementation

  The DNEA R package is the latest implementation of the algorithm and implements several enhancements to the workflow. GGM's require considerably more compute power and, as biological datasets grow, the available resources on the typical personal computer can be constraining. DNEA is designed to work on high-performance or cloud-computing machines to utilize all available computing cores through parallelization, while keeping the memory footprint at a minimum. The algorithm has also been modularized into several key steps. The user now has more control over model tuning and edge inclusion through filtering of the resulting networks based on the strength of the metabolite-metabolite associations. This enables tailoring of the algorithm to the input data.

## Tutorial

An walkthrough of the DNEA algorithm is stored in the package in the form of a vignette, and can be accessed by typing `browseVignettes("DNEA")` in the R console. Example data is also stored in the package and instructions on how to access it is available in the vignette.

## Notes

- More information about the previous implementation of DNEA, ***Filigree***, can be found [here](https://metscape.med.umich.edu/filigree.html). 
- Information about an associated tool we developed for the construction of a data-driven network for a single group, ***CorrelationCalculator***, can be found [here](https://metscape.med.umich.edu/calculator.html).
