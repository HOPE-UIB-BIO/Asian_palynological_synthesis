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
        dcca = purrr::pmap(
            .l = list(
                pollen_percentages,
                levels,
                dataset_id
            ),
            .f = run_dcca
        )
    )

# 15 warnings: "Failed to get eigenvalue(Ax1) greater than the 1st
#  unconstrained eigenvalue".

data_dcca_proc <-
    data_dcca %>%
    dplyr::mutate(
        dcca_scores = purrr::map(
            .x = dcca,
            .f = ~ .x$case.r %>%
                as.data.frame() %>%
                tibble::rownames_to_column("sample_id") %>%
                janitor::clean_names() %>%
                dplyr::as_tibble()
        ),
        dcca_grad_length = purrr::map_dbl(
            .x = dcca_scores,
            .f = ~ {
                dcca_range <-
                    .x %>%
                    purrr::pluck("axis_1") %>%
                    range()

                (dcca_range[2] - dcca_range[1]) %>%
                    abs() %>%
                    return()
            }
        )
    )


#----------------------------------------------------------#
# 4. Estimate DCCA -----
#----------------------------------------------------------#
data_dca <-
    data_pollen %>%
    dplyr::mutate(
        dca = purrr::map(
            .x = pollen_percentages,
            .f = ~ REcopol::fit_ordination(
                data_source_community = .x,
                data_source_predictors = NULL,
                sel_method = "dca"
            )
        )
    )

data_dca_proc <-
    data_dca %>%
    dplyr::mutate(
        dca_scores = purrr::map(
            .x = dca,
            .f = ~ vegan::scores(
                x = .x,
                display = "sites",
                scaling = "sites",
                origin = FALSE
            ) %>%
                as.data.frame() %>%
                tibble::rownames_to_column("sample_id") %>%
                janitor::clean_names() %>%
                tibble::as_tibble()
        ),
        dca_grad_length = purrr::map_dbl(
            .x = dca_scores,
            .f = ~ {
                dcca_range <-
                    .x %>%
                    purrr::pluck("dca1") %>%
                    range()

                (dcca_range[2] - dcca_range[1]) %>%
                    abs() %>%
                    return()
            }
        )
    )


#----------------------------------------------#
# 5. Save ----
#----------------------------------------------#

data_turnover <-
    dplyr::inner_join(
        data_dcca_proc %>%
            dplyr::select(dataset_id, dcca, dcca_scores, dcca_grad_length),
        data_dca_proc %>%
            dplyr::select(dataset_id, dca, dca_scores, dca_grad_length),
        by = "dataset_id"
    )

readr::write_rds(
    data_turnover,
    here::here(
        "Data/Processed/Turnover/Data_turnover_2022-09-14.rds"
    ),
    compress = "gz"
)
