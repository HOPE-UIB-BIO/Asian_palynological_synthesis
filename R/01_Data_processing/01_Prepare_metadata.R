#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                   Prepare metadata
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# This script prepares meta data information

#----------------------------------------------------------#
# 1. Set up  -----
#----------------------------------------------------------#

library(here)

# Load configuration
source(
    here::here("R/00_Config_file.R")
)

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
# 3. Prepare metadata -----
#----------------------------------------------------------#

data_meta <-
    data_pollen %>%
    # dataset beyond the date-line is removed
    dplyr::filter(long >= 0 & long <= 180) %>%
    # single dataset of "Tropical_Monsoon" zone is removed
    dplyr::filter(!ecozone_koppen_15 == "Tropical_Monsoon") %>% # 205 datasets
    # Modify the climate zones as suggested by John
    dplyr::mutate(
        Climate_zone = dplyr::case_when(
            ecozone_koppen_15 == "Arid_Desert" ~ "Arid",
            ecozone_koppen_15 == "Arid_Steppe" ~ "Arid",
            ecozone_koppen_15 == "Cold_Dry_Summer" ~ "Cold_Dry",
            ecozone_koppen_15 == "Cold_Dry_Winter" ~ "Cold_Dry",
            ecozone_koppen_15 == "Temperate_Dry_Summer" ~ "Temperate",
            ecozone_koppen_15 == "Temperate_Dry_Winter" ~ "Temperate",
            ecozone_koppen_15 == "Temperate_Without_dry_season" ~ "Temperate",
            ecozone_koppen_15 == "Polar_Tundra" ~ "Polar",
            TRUE ~ ecozone_koppen_15
        )
    ) %>%
    dplyr::select(
        dataset_id, handle, sitename,
        long, lat, altitude,
        depositionalenvironment,
        country,
        Climate_zone,
        source_of_data, data_publicity
    )

#----------------------------------------------------------#
# 4. Save -----
#----------------------------------------------------------#

readr::write_csv(
    x = data_meta,
    file = here::here("Outputs/Tables/Metadata.csv")
)

readr::write_rds(
     x = data_meta,
    file = here::here("Data/Processed/Metadata/Metadata-2022-09-15.rds")
)
