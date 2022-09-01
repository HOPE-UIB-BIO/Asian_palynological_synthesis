
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
                pollen_percentages, # ..2
                orig_pollen_percentage, # ..3
                levels, # ..4
                age_uncertainty_matrix, # ..5

                dataset_id # ..6
            ),
            .f = ~ {
                message(
                    msg = paste("dataset", ..6)
                )

                if (
                    ..3 == TRUE
                ) {
                    data_sel <- ..2 %>%
                        dplyr::mutate(
                            dplyr::across(
                                tidyselect:::where(is.numeric), ~ .x / 100
                            )
                        )
                } else {
                    data_sel <- ..1
                }

                try(
                    roc_res <-
                        RRatepol::fc_estimate_RoC(
                            data_source_community = data_sel,
                            data_source_age = ..4,
                            age_uncertainty = ..5,
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
                            standardise = !..3,
                            N_individuals = n_individuals_to_standardise, # [config_criteria]
                            tranform_to_proportions = !..3,
                            DC = transformation_coef, # [config_criteria]
                            interoc_rest_throc_reshold = age_max, # [config_criteria]
                            time_standardisation = size_of_bin # [config_criteria]
                        ),
                    silent = TRUE
                )

                if (!exists("roc_res")) {
                    roc_res <- NA
                }

                return(roc_res)
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
        "Data/Processed/RoC/Data_roc_2022-07-31.rds"
    ),
    compress = "gz"
)
