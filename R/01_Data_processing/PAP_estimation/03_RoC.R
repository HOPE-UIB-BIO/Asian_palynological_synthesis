
#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                 Estimate Rate-of-Change
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# Estimation of vegetation rate-of-change

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
# 3. Estimate RoC -----
#----------------------------------------------------------#

data_work_roc <-
    data_pollen %>%
    dplyr::mutate(
        PAP_roc = purrr::pmap(
            .l = list(
                pollen_counts, # ..1
                levels, # ..2
                age_uncertainty_matrix, # ..3
                dataset_id # ..4
            ),
            .f = ~ {
                message(
                    msg = paste("dataset", ..4)
                )

                res <-
                    RRatepol::fc_estimate_RoC(
                        data_source_community = ..1,
                        data_source_age = ..2,
                        age_uncertainty = ..3,
                        smooth_method = smoothing_method, # [config_criteria]
                        smooth_N_points = min_points_smoothing, # [config_criteria]
                        smooth_N_max = max_points_smoothing, # [config_criteria]
                        smooth_age_range = age_range_smoothing, # [config_criteria]
                        Working_Units = working_units_selection, # [config_criteria]
                        bin_size = size_of_bin, # [config_criteria]
                        Number_of_shifts = n_mowing_windows, # [config_criteria]
                        bin_selection = which_level_select_in_bin, # [config_criteria]
                        rand = n_rand, # [config_criteria]
                        use_parallel = TRUE,
                        standardise = TRUE,
                        N_individuals = n_individuals_to_standardise, # [config_criteria]
                        tranform_to_proportions = TRUE,
                        DC = transformation_coef, # [config_criteria]
                        interest_threshold = age_max, # [config_criteria]
                        time_standardisation = size_of_bin # [config_criteria]
                    )

                return(res)
            }
        )
    )

#----------------------------------------------------------#
# 4. Save -----
#----------------------------------------------------------#

data_roc <-
    data_work_roc %>%
    dplyr::select(dataset_id, PAP_roc)

readr::write_rds(
    data_roc,
    here::here(
        "Data/Processed/RoC/Data_roc_2022-07-30.rds"
    ),
    compress = "gz"
)
