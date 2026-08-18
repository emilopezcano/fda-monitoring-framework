## System settings

### All warnings are ignored (change to 0 to show them)
options(warn = -1)

### Numbers are never shown in scientific notation
options(scipen = 99)

### Strings are not treated as factors when importing CSV files
options(stringsAsFactors = FALSE)


## This expression is not usually needed, uncomment in case you
## find issues with your locale
# Sys.setlocale("LC_TIME", "English")

## Load required packages

library(extrafont)
library(data.table)
library(fda.usc)
library(qcr)
library(imputeTS)
library(fdahotelling)
library(ggplot2)
library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(parallel)
library(doSNOW)
library(doRNG)
library(patchwork)

## Create parameters list

params <- list()

### Date of execution

params$date <- Sys.Date()

### Set random seeds

params$seed <- list(
    case_study_1 = list(
        hotelling = 202611,
        l1_std = 202612,
        l2_std = 202613,
        q = 202614
    ),
    case_study_2 = list(
        hotelling = 202621,
        l1_std = 202622,
        l2_std = 202623,
        q = 202624
    )
)

## Load helper functions

source("R/001_helper_functions.R")

## Load fonts for plots

loadfonts()
