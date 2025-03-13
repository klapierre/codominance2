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
| `code/04_generatingTraitData/04_NutNet_traits.R`                  | imputing continous and gathering categorical trait data for NutNet species not in CoRRE Trait Data      |
| `code/05_Question01/05a_site_mode.R`                              |  calculates num codominant modes per site in control plots only         |
| `code/05_Question01/05b_map.R`                                    |  create dataset for map generation                                      |
| `code/05_Question01/05c_Q1_multinomial_model.R`                   |  multinomial model to determine drivers of codominance globally         |



| `code/analysis_q1.R`                                              | multinomial model relating \# of codominance to envrionmental variables |
| `code/format_data.R`                                              | entry point of data formatting; combining all codominance data.         |
| `code/map_mode.R`                                                 | generating source information for GIS                                   |
| `code/allTraits_mergeFile.R`                                      | getting all trait data together                                         |
| `code/sample_gawdist`                                             | sample code for calculating dendrogram trait distance `gawdist()`       |
