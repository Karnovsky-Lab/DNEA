The example data provided in this package, as well as the data sets generated for analysis in the associated paper, were created using the R scripts provided in this folder. 

All data was downloaded from the Metabolomics workbench here: https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR000950

LMadjust.R -- This file contains LMadjust(), the function used to adjust the expression data for age and sex. Each metabolite is modeled using linear regression with the metabolite as the dependent variable and the covariates to adjust for as the dependent variables. The residuals of the model are then taken and transformed by adding the minimum value to every residual so that the new minimum is 0. 

TEDDY-metadata_process.R -- This file was used to process the metadata file (downloaded from metabolomics workbench).

TEDDY-T1D_process.R -- This file was used to create the TEDDY example data in this package as well as the Type 1 Diabetes (T1D) data sets used in the associated paper (downloaded from metabolomics workbench).

TEDDY-IA_process.R -- This file was used to create the Islet auto-antibody (IA) data sets used in the associated paper (downloaded from metabolomics workbench).