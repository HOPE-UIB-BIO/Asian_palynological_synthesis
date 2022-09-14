#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                     Combine PAP
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# Combine PAPs ----

#----------------------------------------------------------#
# 1. Set up  -----
#----------------------------------------------------------#

library(here)

# Load configuration
source(
    here::here("R/00_Config_file.R")
)

#--------------------------------------------------------#
# 2. Load the data of PAP estimates ----
#--------------------------------------------------------#

data_diversity <-
    readr::read_rds(
        here::here(
            "Data/Processed/Diversity/Data_diverity_2022-07-30.rds"
        )
    )

data_mrt <-
    readr::read_rds(
        here::here(
            "Data/Processed/MVRT/Data_mvrt_2022-09-14.rds"
        )
    )

data_roc <-
    readr::read_rds(
        here::here(
            "Data/Processed/RoC/Data_roc_2022-07-31.rds"
        )
    )

data_turnover <-
    readr::read_rds(
        here::here(
            "Data/Processed/Turnover/Data_turnover_2022-09-14.rds"
        )
    )

data_levels <-
    readr::read_rds(
        here::here(
            "Data/Input/Data_processed_2022-08-30.rds"
        )
    ) %>%
    dplyr::select(
        dataset_id,
        levels
    )

#--------------------------------------------------------#
# 3. Combine and re-structture the PAPs ----
#--------------------------------------------------------#

subset_by_vector <-
    function(data_source, var_name, id_vec) {
        data_source %>%
            dplyr::mutate(
                !!var_name := purrr::map2(
                    .x = get(var_name),
                    .y = get(id_vec),
                    .f = ~ .x %>%
                        dplyr::filter(.data$sample_id %in% .y)
                )
            ) %>%
            return()
    }

data_combine_paps <-
    data_levels %>%
    dplyr::inner_join(
        data_diversity,
        by = "dataset_id"
    ) %>%
    dplyr::inner_join(
        data_mrt,
        by = "dataset_id"
    ) %>%
    dplyr::inner_join(
        data_roc,
        by = "dataset_id"
    ) %>%
    dplyr::inner_join(
        data_turnover,
        by = "dataset_id"
    ) %>%
    # in order to make sure we have same levels across all data
    dplyr::mutate(
        # get a list of intercept of all samples across data
        valid_sample_id = purrr::pmap(
            .l = list(
                PAP_diversity, # ..1
                mvrt_partitions, # ..2
                dcca_scores, # ..3
                dca_scores, # ..4
                levels # ..5
            ),
            .f = ~ dplyr::inner_join(
                ..1,
                ..2,
                by = "sample_id"
            ) %>%
                dplyr::inner_join(
                    ..3,
                    by = "sample_id"
                ) %>%
                dplyr::inner_join(..4,
                    by = "sample_id"
                ) %>%
                dplyr::inner_join(..5,
                    by = "sample_id"
                ) %>%
                purrr::pluck("sample_id")
        )
    ) %>%
    # subset all data.frames by the list of common sample_id
    subset_by_vector(
        var_name = "PAP_diversity",
        id_vec = "valid_sample_id"
    ) %>%
    subset_by_vector(
        var_name = "mvrt_partitions",
        id_vec = "valid_sample_id"
    ) %>%
    subset_by_vector(
        var_name = "dcca_scores",
        id_vec = "valid_sample_id"
    ) %>%
    subset_by_vector(
        var_name = "dca_scores",
        id_vec = "valid_sample_id"
    ) %>%
    subset_by_vector(
        var_name = "levels",
        id_vec = "valid_sample_id"
    )  %>% 
    dplyr::select(-valid_sample_id)

#--------------------------------------------------------#
# 4. Save data ----
#--------------------------------------------------------#
readr::write_rds(
    data_combine_paps,
    here::here(
        "Data/Processed/PAP_all/pap_all_2022-09-14.rds"
    ),
    compress = "gz"
)