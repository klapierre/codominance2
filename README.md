README
================
2025-02-26

## README

README file for codominance project

## File description

| file                                                              | note                                                                    |
|:------------------------------------------------------------------|:------------------------------------------------------------------------|
| `code/01_library.R`                                               | libraries to load                                                       |
| `code/02_functions.R`                                             | functions used throughout                                               |
| `code/03_generatingCodominaceData/03a_CoRRE_codominance.R`        | determining codominance for CoRRE database                              |
| `code/03_generatingCodominaceData/03b_GEx_codominance.R`          | determining codominance for GEx database                                |
| `code/03_generatingCodominaceData/03c_NutNet_codominance.R`       | determining codominance for NutNet database                             |
| `code/03_generatingCodominaceData/03d_combine_data.R`             | merging databases: generates allSppList.rds, codomSppList.rds, numCodomPlotYear.rds, and envData.rds  |
| `code/04_generatingTraitData/04a_NutNet_traits.R`                 | imputing continous and gathering categorical trait data for NutNet species not in CoRRE Trait Data      |
| `code/04_generatingTraitData/04b_allTraits_mergeFile.R`           |  merge trait data from CoRRE database and new NutNet trait data         |
| `code/05_Question01/05_Q1_site_modes.R`                           |  calculates num codominant modes per site in control plots only         |
| `code/06_Question02/06a_Q2_generate_codom_groups.R`               |  generates lists of codominating species that most commonly occur in each experiment         |
| `code/06_Question02/sample_gawdist`                               | sample code for calculating dendrogram trait distance `gawdist()`       |
