#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                Prepare data for analyses
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# This script to prepare dataset for analyses

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
data_meta <-
    readr::read_rds(
        here::here(
            "Data/Processed/Metadata/Metadata-2022-09-15.rds"
        )
    )

data_combine_paps <-
    readr::read_rds(
        here::here(
            "Data/Processed/PAP_all/pap_all_2022-09-14.rds"
        )
    )

data_density <-
    readr::read_rds(
        here::here(
            "Data/Processed/Partitions/PAP_density_2022-09-15.rds"
        )
    )

#----------------------------------------------------------#
# 3. Prepare dataset -----
#----------------------------------------------------------#

data_join <-
    data_meta %>%
    dplyr::left_join(
        data_combine_paps,
        by = "dataset_id"
    ) %>%
    dplyr::left_join(
        data_density,
        by = "dataset_id"
    )

data_for_analyses <-
    data_join %>%
    dplyr::mutate(
        PAP_merge = purrr::pmap(
            .l = list(
                levels, # ..1
                PAP_diversity, # ..2
                dcca_scores, # ..3
                dca_scores # ..4
            ),
            .f = ~ ..1 %>%
                dplyr::select(
                    sample_id,
                    age
                ) %>%
                dplyr::inner_join(
                    ..2,
                    by = "sample_id"
                ) %>%
                dplyr::inner_join(
                    ..3 %>%
                        dplyr::select(
                            sample_id,
                            dcca_axis_1 = axis_1
                        ),
                    by = "sample_id"
                ) %>%
                dplyr::inner_join(
                    ..4 %>%
                        dplyr::select(
                            sample_id,
                            dca_axis_1 = dca1
                        ),
                    by = "sample_id"
                )
        )
    ) %>%
    dplyr::select(
        dataset_id,
        long, lat,
        Climate_zone,
        dcca_grad_length,
        dca_grad_length,
        mvrt_groups_n,
        PAP_merge,
        PAP_roc,
        pap_density
    )


#----------------------------------------------------------#
# 4. Save  -----
#----------------------------------------------------------#

readr::write_rds(
    data_for_analyses,
    file = here::here(
        "Data/Processed/Data_for_analyses/Data_for_analyses-2022-09-15.rds"
    ),
    compress = "gz"
)
