#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                    Vegetation turnover
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# DCCA

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

# Input data is the full dataset of Asia with 207 sites.
data_pollen <-
    readr::read_rds(
        here::here(
            "Data/Input/Data_processed_2022-08-30.rds"
        )
    )

#----------------------------------------------------------#
# 3. Estimate DCCA -----
#----------------------------------------------------------#

data_dcca <-
    data_pollen %>%
    dplyr::mutate(
        dcca = purrr::map2(
            .x = pollen_percentages,
            .y = levels,
            .f = ~ REcopol::fit_ordination(
                data_source_community = .x,
                data_source_predictors = .y,
                sel_method = "constrained",
                var_name_pred = "age",
                sel_complexity = "poly_2",
                transform_to_percentage = FALSE,
                tranformation = "none"
            )
        )
    )

data_dcca_proc <-
    data_dcca %>%
    dplyr::mutate(
        dcca_scores = purrr::map(
            .x = dcca,
            .f = ~ .x %>% purrr::pluck("case_r")
        ),
        dcca_grad_length = purrr::map_dbl(
            .x = dcca,
            .f = ~ .x %>% purrr::pluck("axis_1_grad_length")
        )
    )


#----------------------------------------------#
# 4. Save ----
#----------------------------------------------#

data_turnover <-
    data_dcca_proc %>%
    dplyr::select(
        dataset_id,
        dcca,
        dcca_scores,
        dcca_grad_length
    )

readr::write_rds(
    data_turnover,
    here::here(
        "Data/Processed/Turnover/Data_turnover_2022-09-29.rds"
    ),
    compress = "gz"
)
