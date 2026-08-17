<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](https://sandbox.zenodo.org/badge/1337052177.svg)](https://handle.test.datacite.org/10.5072/zenodo.588641)

# A Nonparametric Phase II Monitoring Framework for Functional Data via Mahalanobis Distances

The goal of this repository is to provide reproducible code for the
methods and use cases on the paper “A Nonparametric Phase II Monitoring
Framework for Functional Data via Mahalanobis Distances” by Priscila
Guayasamín, Miguel Flores, Emilio L. Cano, Javier M. Moguerza, Javier
Tarrío‐Saavedra, and Salvador Naya, submitted to the journal *Quality
and Reliability Engineering International*.

This readme file follows the structure of the paper, and it is divided
into sections that correspond to the sections of the paper. Each section
contains a brief description of the methods and use cases, as well as
links to the corresponding code and data.

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
    #> [1] stats     graphics  grDevices utils     datasets 
    #> [6] methods   base     
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] jsonlite_2.0.0    compiler_4.6.1    renv_1.2.3       
    #>  [4] gitcreds_0.1.2    callr_3.8.0       credentials_2.0.3
    #>  [7] yaml_2.3.12       fastmap_1.2.0     R6_2.6.1         
    #> [10] pak_0.11.0        curl_7.1.0        httr2_1.3.0      
    #> [13] knitr_1.51        tibble_3.3.1      rprojroot_2.1.1  
    #> [16] openssl_2.4.2     pillar_1.11.1     rlang_1.3.0      
    #> [19] cachem_1.1.0      xfun_0.60         fs_2.1.0         
    #> [22] sys_3.4.3         pkgload_1.5.3     otel_0.2.0       
    #> [25] memoise_2.0.1     cli_3.6.6         withr_3.0.3      
    #> [28] magrittr_2.0.5    ps_1.9.3          processx_3.9.0   
    #> [31] digest_0.6.39     gh_1.6.1          rstudioapi_0.19.0
    #> [34] markdown_2.0      devtools_2.5.2    askpass_1.2.1    
    #> [37] gert_2.3.1        lifecycle_1.0.5   vctrs_0.7.3      
    #> [40] evaluate_1.0.5    glue_1.8.1        whisker_0.4.1    
    #> [43] sessioninfo_1.2.4 pkgbuild_1.4.8    rmarkdown_2.31   
    #> [46] purrr_1.2.2       tools_4.6.1       usethis_3.2.1    
    #> [49] pkgconfig_2.0.3   ellipsis_0.3.3    htmltools_0.5.9

Run this code to update the README.md file from README.Rmd, if needed:

    rmarkdown::render("README.Rmd", output_format = "md_document")
