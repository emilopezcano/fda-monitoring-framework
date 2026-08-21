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

Please note that the {renv} package must be installed to use this
functionality. It is automatically installed when starting R from the
repository folder for the first time. Furthermore, for Windows systems
it is needed to install the
[Rtools](https://cran.r-project.org/bin/windows/Rtools/), and for MacOS
systems [XCode Tools and gfortran](https://mac.r-project.org/tools/),
follow the links for instructions.

Optionally, you can manually install the required packages, but complete
reproducibility is only guaranteed when using the same package versions
as those used in the paper. The required packages are the following
(note that one of them is not available on CRAN and must be installed
from github, for example using the {pak} package):

- {extrafont}
- {data.table}
- {fda.usc}
- {qcr}
- {imputeTS}
- {fdahotelling}, not in CRAN: <https://github.com/astamm/fdahotelling>.
  Install using \`pak::pak(“astamm/fdahotelling”)
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
R and of the loaded packages used in the last update of this repository.
Even though everyting is likely to work fine with recent versions of R,
older ones may find issues. Thus, it is recommended to use the same as
the reported by `sessionInfo()` and recorded in the `renv.lock` file.

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

## Simulations

The simulations scripts are located at `R/simulations`. There is a
folder for each scenario, and within each folder there are scripts for
each simulation setting. The results of the simulations are saved in the
`results/simulations` folder. The script `R/simulations/000_setup_sim.R`
contains the code to setup the simulations, including the parameters and
the parallelization setup. The scripts for each simulation setting are
named according to the following pattern:
`<stat>_<method>_<n1>_<n2>_<scenario>.R`, where `<stat>` is the
statistic used, `<method>` is the method used, `<n1>` is the calibration
sample size, `<n2>` is the monitoring sample size, and `<scenario>` is
the scenario code.

Before running the individual simulation scripts, the
`R/simulations/000_setup_sim.R` script must be run to setup the
simulation common parameters, as well as the paralellization
configuration. You can run it with the following code:

    source("R/simulations/000_setup_sim.R")

Then, the simulation scripts can be run in order, or you can run them
individually. The results of the simulations are saved in the
`results/simulations` folder, and they can be loaded with the `load()`
function. The results are saved as RData files, and they contain the
objects `senal` and `potencia`, which correspond to the out-of-control
signal and the power of the test, respectively.

When the simulations are finished, it is advisable to stop the parallel
cluster to free up resources. You can do this with the following code:

    stopCluster(cl)

**NOTE**: Even though the simulations are paralelized, they may take a
long time to run, depending on the number of simulations and the
computational resources available. The parameters in the
`R/simulations/000_setup_sim.R` can be modified to change the number of
simulations, the number of bootstrap samples, and the number of
permutations.

The following expression can be used to run all the simulations in
order:

    source("R/simulations/000_setup_sim.R")
    for (f in list.files(
      "R/simulations/Scenario1A",
      pattern = "\\.R$",
      full.names = TRUE
    )) {
      source(f)
    }
    for (f in list.files(
      "R/simulations/Scenario1B",
      pattern = "\\.R$",
      full.names = TRUE
    )) {
      source(f)
    }
    for (f in list.files(
      "R/simulations/Scenario2A",
      pattern = "\\.R$",
      full.names = TRUE
    )) {
      source(f)
    }
    for (f in list.files(
      "R/simulations/Scenario2B",
      pattern = "\\.R$",
      full.names = TRUE
    )) {
      source(f)
    }

Nevertheless, most of the scripts take hours to run, so to reproduce
everything it is better to plan in advance the execution in appropriate
High Performance Computing (HPC) resources.

The results of each simulation can be visualized by loading the objects
and showing the results. For example, to visualize the results of the
simulation `l1std_bootstrap_150_100_1A.R`, you can run the following
code:

    load("results/simulations/l1std_boot_100_50_1A.RData")
    print(potencia_l1std_boot_100_50)
    print(mean(senal_l1std_boot_100_50[[1]]))

**NOTE**: The results may not be exactly the same as those reported in
the paper, since the simulations are based on random samples that can
vary across R versions and hardware architectures. However, they should
be very similar, and the conclusions should be the same.

## Case study 1: Monitoring of NO<sub>2</sub> concentration in Madrid (COVID-19 lockdown)

The original data for this case study is available in the
`data/case_study_2` folder. The code to reproduce the results of this
case study is available in the `R/case_study_2` folder. The script
`R/case_study_2/000_data_prep.R` prepares the data, whereas
`R/case_study_2/001_case_study_2_phaseI.R` and
`R/case_study_2/002_case_study_2_phaseII_monitoring.R` run the Phase I
and Phase II monitoring respectively. The results are saved in the
`results/case_study_2` folder. Thus, you can open and run the code of
the scripts line by line to check all the intermediate results, or run
the following code to run the scripts in order:

    source("R/case_study_1/000_data_prep.R")
    source("R/case_study_1/001_case_study_1_phaseI_paired.R")
    source("R/case_study_1/002_case_study_1_phaseII_monitoring_paired.R")

## Case study 2: Monitoring electricity consumption in a retail store (Panama City)

The original data for this case study is available in the
`data/case_study_1` folder. The code to reproduce the results of this
case study is available in the `R/case_study_1` folder. The script
`R/case_study_1/000_data_prep.R` prepares the data, whereas
`R/case_study_1/001_case_study_1_phaseI_paired.R` and
`R/case_study_1/002_case_study_1_phaseII_monitoring_paired.R` run the
Phase I and Phase II monitoring respectively. The results are saved in
the `results/case_study_1` folder. Thus, you can open and run the code
of the scripts line by line to check all the intermediate results, or
run the following code to run the scripts in order:

    source("R/case_study_2/000_data_prep.R")
    source("R/case_study_2/001_case_study_2_phaseI.R")
    source("R/case_study_2/002_case_study_2_phaseII_monitoring.R")

## ARL experiment

The code for the ARL experiment in the conclusions can be found in the
script `R/ARL/ARL_simulations.R`. It contains helper functions, data
generation, simulations, and ARL computations. The results are saved in
the `results/ARL_experiment` folder. The code can be run line by line
from the script to see the intermediate results, or all at once with the
following expression:

    source("R/ARL/ARL_simulations.R")

# Session Info

The following code chunk prints out the session information
corresponding to the last repository update, which is useful for
reproducibility.

    sessionInfo()
    #> R version 4.6.1 (2026-06-24)
    #> Platform: aarch64-apple-darwin23
    #> Running under: macOS Tahoe 26.6.2
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
    #>  [1] compiler_4.6.1  fastmap_1.2.0   cli_3.6.6      
    #>  [4] htmltools_0.5.9 tools_4.6.1     yaml_2.3.12    
    #>  [7] rmarkdown_2.31  knitr_1.51      xfun_0.60      
    #> [10] digest_0.6.39   rlang_1.3.0     renv_1.2.4     
    #> [13] evaluate_1.0.5

Run this code to update the README.md file from README.Rmd, if needed:

    rmarkdown::render("README.Rmd", output_format = "md_document")
