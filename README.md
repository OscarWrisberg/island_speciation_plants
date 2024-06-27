# Island Speciation Plants

## Description

This repository contains data and scripts related to the study of plant speciation on islands. The project involves various analyses, including model fitting, plotting, and statistical testing, to understand the patterns and processes of plant speciation in island ecosystems.

## Project Structure

- **R Scripts and Markdown Files**: These files contain the core analysis scripts and documentation.
  - `All_species_mono_phylo_islands.Rmd`
  - `Cory_mono_phylo_islands.Rmd`
  - `Data_setup_and_functions.Rmd`
  - `Endemism_phylo_plot.Rmd`
  - `Max_height_calculation.R`
  - `Model_fitting.Rmd`
  - `Model_fitting_other_models.Rmd`
  - `Plotting_script.R`
  - `Scatter_plots.Rmd`
  - `Supplement_figures.Rmd`
  - `Total_sp_testing.R`
  - `Vif_all_models.R`
  - `islands_with_no_endemics_supp.R`
  - `volcanic_botanical_countries_models.Rmd`

- **Data Directory**: Contains datasets used in the analysis.
  - `Output_all_sp_test.rds`
  - `island_endmics2_BOR.rds`
  - `island_endmics_all_species.rds`
  - `output_for_test_all_sp`
  - `output_for_test_all_sp.rds`
  - `output_for_test_all_sp_subset_log`
  - `output_for_test_coryphoideae_subset_log`
  - `output_with_island_data.rds`

- **Figures Directory**: Contains figures generated from the analysis.
  - `Anagenesis_hump.pdf`
  - `Botanical_countries_map.pdf`
  - `scatter_plots_speciation_all_sp.pdf`
  - `scatter_plots_speciation_all_sp_facet.pdf`
  - And more...

- **LICENSE**: Licensing information.

## Installation

To run the analyses in this repository, you'll need to have R installed on your system. You can install R from [CRAN](https://cran.r-project.org/).

Additionally, you will need to install several R packages. You can install them by running the following command in your R console:

```R
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble"))
