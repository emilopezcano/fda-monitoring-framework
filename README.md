<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](https://sandbox.zenodo.org/badge/1337052177.svg)](https://handle.test.datacite.org/10.5072/zenodo.588641)

# A Nonparametric Phase II Monitoring Framework for Functional Data via Mahalanobis Distances

The goal of this repository is to provide reproducible code for the
methods and case studies on the paper “A Nonparametric Phase II
Monitoring Framework for Functional Data via Mahalanobis Distances” by
Priscila Guayasamín, Miguel Flores, Emilio L. Cano, Javier M. Moguerza,
Javier Tarrío‐Saavedra, and Salvador Naya, submitted to the journal
*Quality and Reliability Engineering International*.

This readme file follows the structure of the paper, and it is divided
into sections that correspond to the sections of the paper. Each section
contains a brief description of the methods and case studies, as well as
links to the corresponding code and data.

## Environment setup

In order to reproduce the results of the paper, you can clone this
repository and open it with your preferred IDE, e.g., RStudio, Positron,
or even a terminal window and the R console.

This repository is configured to use the
[{renv}](https://rstudio.github.io/renv/) package to manage
dependencies. To install the required packages, you can run the
following code in your R console:

    renv::restore()

Please note that you need the {renv} package installed to use this
functionality.

Optionally, you can manually install the required packages, but complete
reproducibility is only guaranteed when using the same package versions
as those used in the paper. The required packages are the following
(note that two of them are not available on CRAN and must be installed
from github, for example using the {pak} package):

- {extrafont}
- {extrafont}
- {data.table}
- {fda.usc}
- {qcr}
- {imputeTS}
- {fdahotelling}, not in CRAN: <https://github.com/astamm/fdahotelling>
- {ggplot2}
- {readxl}
- {janitor}
- {dplyr}
- {tidyr}
- {parallel}
- {doSNOW}
- {doRNG}
- {patchwork}

The `sessionInfo()` output at the end of this file contains the versions
of the loaded packages used in the last update of this repository.

The script `R/000_setup.R` contains the code to setup the environment
and load the required packages. You can run it with the following code:

    source("R/000_setup.R")

In addition to loading the packages, the script also creates the
`params` object, a list that contains the parameters used throughout the
code. You can modify the parameters in this list to change the behavior
of the code. Furthermore, the helper functions in
`R/001_helper_functions.R` are loaded to provide additional
functionality. Some other expressions are run to configure input and
output options.

## Section 2: Methodology

Algorithm 1 in the paper describes the proposed methodology for Phase II
monitoring of functional data using Mahalanobis distances. The code is
implemeted …..

## Section 4: Results and Discussion

### Simulations

### Case study 1: Monitoring of *N**O*<sub>2</sub> concentration in Madrid (COVID-19 lockdown)

### Case study 2: Monitoring electricity consumption in a retail store (Panama City)

# Session Info

The following code chunk prints out the session information
corresponding to the last repository update, which is useful for
reproducibility.

    sessionInfo()
    #> R version 4.6.1 (2026-06-24)
    #> Platform: aarch64-apple-darwin23
    #> Running under: macOS Tahoe 26.6.1
    #> 
    #> Matrix products: default
    #> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
    #> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    #> 
    #> locale:
    #> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    #> 
    #> time zone: Europe/Madrid
    #> tzcode source: internal
    #> 
    #> attached base packages:
    #> [1] parallel  splines   stats     graphics 
    #> [5] grDevices utils     datasets  methods  
    #> [9] base     
    #> 
    #> other attached packages:
    #>  [1] patchwork_1.3.2        
    #>  [2] doRNG_1.8.6.3          
    #>  [3] rngtools_1.5.2         
    #>  [4] doSNOW_1.0.20          
    #>  [5] snow_0.4-4             
    #>  [6] iterators_1.0.14       
    #>  [7] foreach_1.5.2          
    #>  [8] tidyr_1.3.2            
    #>  [9] dplyr_1.2.1            
    #> [10] janitor_2.2.1          
    #> [11] readxl_1.5.0           
    #> [12] ggplot2_4.0.3          
    #> [13] fdahotelling_0.0.0.9000
    #> [14] imputeTS_3.4           
    #> [15] qcr_1.4                
    #> [16] mvtnorm_1.4-2          
    #> [17] qcc_2.7                
    #> [18] fda.usc_2.2.0          
    #> [19] knitr_1.51             
    #> [20] mgcv_1.9-4             
    #> [21] nlme_3.1-169           
    #> [22] fda_6.3.0              
    #> [23] deSolve_1.42           
    #> [24] fds_1.9                
    #> [25] RCurl_1.98-1.19        
    #> [26] rainbow_3.8            
    #> [27] pcaPP_2.0-5            
    #> [28] MASS_7.3-65            
    #> [29] data.table_1.18.4      
    #> [30] extrafont_0.20         
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] tidyselect_1.2.1   hdrcde_3.5.0      
    #>  [3] timeDate_4052.112  farver_2.1.2      
    #>  [5] S7_0.2.2           bitops_1.1-0      
    #>  [7] fastmap_1.2.0      pracma_2.4.6      
    #>  [9] digest_0.6.39      timechange_0.4.0  
    #> [11] lifecycle_1.0.5    cluster_2.1.8.2   
    #> [13] magrittr_2.0.5     compiler_4.6.1    
    #> [15] rlang_1.3.0        tools_4.6.1       
    #> [17] yaml_2.3.12        labeling_0.4.3    
    #> [19] mclust_6.1.3       xml2_1.6.0        
    #> [21] RColorBrewer_1.1-3 SuppDists_1.1-9.9 
    #> [23] KernSmooth_2.23-26 withr_3.0.3       
    #> [25] purrr_1.2.2        grid_4.6.1        
    #> [27] colorspace_2.1-3   extrafontdb_1.1   
    #> [29] scales_1.4.0       cli_3.6.6         
    #> [31] rmarkdown_2.31     generics_0.1.4    
    #> [33] stringr_1.6.0      forecast_9.0.2    
    #> [35] urca_1.3-4         cellranger_1.1.0  
    #> [37] vctrs_0.7.3        Matrix_1.7-5      
    #> [39] stinepack_1.5      pak_0.11.1        
    #> [41] glue_1.8.1         codetools_0.2-20  
    #> [43] ggtext_0.1.2       lubridate_1.9.5   
    #> [45] stringi_1.8.9      gtable_0.3.6      
    #> [47] tibble_3.3.1       pillar_1.11.1     
    #> [49] htmltools_0.5.9    R6_2.6.1          
    #> [51] ks_1.15.3          doParallel_1.0.17 
    #> [53] evaluate_1.0.5     lattice_0.22-9    
    #> [55] gridtext_0.1.6     snakecase_0.11.1  
    #> [57] renv_1.2.4         fracdiff_1.5-4    
    #> [59] Rcpp_1.1.2         Rttf2pt1_1.3.14   
    #> [61] xfun_0.60          kSamples_1.2-12   
    #> [63] zoo_1.9-0          pkgconfig_2.0.3

Run this code to update the README.md file from README.Rmd, if needed:

    rmarkdown::render("README.Rmd", output_format = "md_document")
