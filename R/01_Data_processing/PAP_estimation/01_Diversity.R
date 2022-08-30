#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                   Estimate diversity
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# This script prepares data for diversity analyses and saves the diversity
#  results in the tibbles for later use in subsequent analyses.

#----------------------------------------------------------#
# 1. Set up  -----
#----------------------------------------------------------#

library(here)

# Load configuration
source(
    here::here("R/00_Config_file.R")
)

set.seed(1234)


#----------------------------------------------------------#
# 2. Load data -----
#----------------------------------------------------------#

# Input data is the full dataset of Asia with 209 sites.
data_pollen <-
    readr::read_rds(
        here::here(
            "Data/Input/Pollen_processed/Data_processed_2022-08-30.rds"
        )
    )


#----------------------------------------------------------#
# 3. Estimate Diversity -----
#----------------------------------------------------------#

data_work_diverity <-
    data_pollen %>%
    dplyr::mutate(
        PAP_diversity = purrr::map(
            .x = pollen_counts,
            .f = ~ REcopol::diversity_estimate(
                data_source = .x,
                sel_method = "taxonomic",
                rand = n_rand # [config_criteria]
            )
        )
    )


#----------------------------------------------------------#
# 4. Save -----
#----------------------------------------------------------#

data_diverity <-
    data_work_diverity %>%
    dplyr::select(dataset_id, PAP_diversity)

readr::write_rds(
    data_diverity,
    here::here(
        "Data/Processed/Diversity/Data_diverity_2022-07-30.rds"
    ),
    compress = "gz"
)
