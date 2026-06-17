README
================
2026-06-17

## README

README file for codominance test dataset for peer review

## System requirements
R version 4.5.2 or later

## Installation guide
(1) Install R version 4.5.2 or later /n
(2) Run code in numeric order, pulling from the test data housed within the peer review repository

## Demo
(1)

## File description

| file                                                              | note                                                                    |
|:------------------------------------------------------------------|:------------------------------------------------------------------------|
| `code/01_library.R`                                               | libraries to load                                                       |
| `code/02_functions.R`                                             | functions used throughout                                               |
| `code/03_generatingCodominaceData/03a_CoRRE_codominance.R`        | determining codominance for CoRRE database                              |
| `code/03_generatingCodominaceData/03b_GEx_codominance.R`          | determining codominance for GEx database                                |
| `code/03_generatingCodominaceData/03c_NutNet_codominance.R`       | determining codominance for NutNet database                             |
| `code/03_generatingCodominaceData/03d_combine_data.R`             | merging databases: generates allSppList.rds, codomSppList.rds, numCodomPlotYear.rds, expInfo.rds, and envData.rds  |
| `code/04_generatingTraitData/04_allTraits_mergeFile.R`            | gathering trait data from CoRRE database: generates allTraits.rds       |
| `code/05_Question01/05a_Q1_format_mode.R`                         | calculates number of dominant species per site in control plots only    |
| `code/05_Question01/05b_Q1_multinomial_model.R`                   | multinomial model of number of dominant species per site                |
| `code/05_Question01/05c_Q1_figures.R`                             | generates map and multinomial model probability figures                 |
| `code/06_Question02/06a_Q2_generate_codom_groups.R`               | generates lists of codominating species that most commonly occur in each experiment         |
| `code/06_Question02/06b_Q2_trait_analysis.R`                      | calculates dendrogram of trait distances and determines RDFD for each codominant pair |
| `code/07_Question03/07a_Q3a_trt_modes.R`                          | calculates number of dominant species in treatment plots                |
| `code/07_Question03/07b_Q3b_sppOverlap.R`                         | determines species overlap in codominants with treatments               |
| `code/07_Question03/07c_Q3c_trait_analysis.R`                     | calculates dendrogram of trait distances and determines RDFD for each codominant pair in treatment plots and compares to control plots |
| `code/08_supplemental/08a_richnessRelationship.R`                 | determines RDFD and codominance relationships with species richness     |
| `code/08_supplemental/08b_plotSizeRelationship.R`                 | determines RDFD and codominance relationships with plot size            |
| `code/08_supplemental/08c_functionalGroups.R`                     | determines RDFD and codominance relationships with plant families       |

