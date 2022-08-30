#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                     Estimate MVRT
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# Multivariate regression tree partition

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
# 3. Estimate MVRT -----
#----------------------------------------------------------#

data_work_mrt <-
  data_pollen %>%
  dplyr::mutate(
    PAP_mrt = purrr::map2(
      .x = pollen_percentages,
      .y = levels,
      .f = ~ REcopol::mv_regression_partition(
          data_source_counts = .x,
          data_source_levels = .y,
          rand = n_rand, # [config_criteria]
          transformation = transformation_coef # [config_criteria]
        )
    )
  )


#----------------------------------------------------------#
# 4. Save -----
#----------------------------------------------------------#

data_mrt <-
    data_work_mrt %>%
    dplyr::select(dataset_id, PAP_mrt)

readr::write_rds(
    data_mrt,
    here::here(
        "Data/Processed/MVRT/Data_mvrt_2022-07-30.rds"
    ),
    compress = "gz"
)
