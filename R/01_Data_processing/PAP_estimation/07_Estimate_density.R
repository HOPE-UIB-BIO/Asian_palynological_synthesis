#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                   Density of change points
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# Create density of regression tree change points of the estimates of
#  Hill's diversity, DCCA1, and MRT temporally

#--------------------------------------------------------#
# 1. Setup  ----
#--------------------------------------------------------#

library(here)

# Load configuration
source(
    here::here("R/00_Config_file.R")
)

#--------------------------------------------------------#
# 2. Load the data of change point estimates ----
#--------------------------------------------------------#
data_change_points <-
    readr::read_rds(
        here::here(
            "Data/Processed/Partitions/PAP_change_points_2022-09-19.rds"
        )
    )

#--------------------------------------------------------#
# 3. Estimate density ----
#--------------------------------------------------------#

data_cp_density_merge <-
    get_density_for_all_vars(
        data_source = data_change_points,
        age_table = new_data_general # [config]
    )

#--------------------------------------------------------#
# 4. Save ----
#--------------------------------------------------------#

data_density <-
    data_cp_density_merge %>%
    dplyr::select(dataset_id, pap_density)

readr::write_rds(
    data_density,
    here::here(
        "Data/Processed/Partitions/PAP_density_2022-09-19.rds"
    ),
    compress = "gz"
)
