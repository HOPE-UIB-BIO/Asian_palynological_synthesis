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
            "Data/Processed/Partitions/PAP_change_points_2022-09-14.rds"
        )
    )

#--------------------------------------------------------#
# 3. Estimate density ----
#--------------------------------------------------------#

# helper function
get_density_rescaled <-
    function(data_source,
             age_table) {
        if (
            is.null(data_source)
        ) {
            res <-
                age_table %>%
                dplyr::mutate(
                    density = 0
                )
        } else {
            res <-
                REcopol::get_density(
                    data_source = data_source,
                    reflected = TRUE,
                    values_range = c(
                        min(age_table$age),
                        max(age_table$age)
                    ),
                    bw = 1000 / max(age_table$age),
                    n = max(age_table$age)
                ) %>%
                dplyr::mutate(
                    age = round(var),
                    density = scales::rescale(density, to = c(0, 1))
                ) %>%
                dplyr::filter(
                    age %in% age_table$age
                ) %>%
                dplyr::select(
                    age, density
                )
        }
        return(res)
    }

data_cp_density <-
    data_change_points %>%
    dplyr::mutate(
        mvrt_cp_density = purrr::map(
            .x = mvrt_cp,
            .f = ~ get_density_rescaled(
                data_source = .x,
                age_table = new_data_general # [config]
            )
        ),
        roc_cp_density = purrr::map(
            .x = roc_cp,
            .f = ~ get_density_rescaled(
                data_source = .x,
                age_table = new_data_general # [config]
            )
        ),
        dcca_cp_density = purrr::map(
            .x = dcca_cp,
            .f = ~ get_density_rescaled(
                data_source = .x,
                age_table = new_data_general # [config]
            )
        ),
        dca_cp_density = purrr::map(
            .x = dca_cp,
            .f = ~ get_density_rescaled(
                data_source = .x,
                age_table = new_data_general # [config]
            )
        ),
        n0_density = purrr::map(
            .x = diversity_cp,
            .f = ~ .x %>%
                dplyr::filter(var_name == "n0") %>%
                purrr::pluck("age") %>%
                get_density_rescaled(
                    data_source = .,
                    age_table = new_data_general # [config]
                )
        ),
        n1_density = purrr::map(
            .x = diversity_cp,
            .f = ~ .x %>%
                dplyr::filter(var_name == "n1") %>%
                purrr::pluck("age") %>%
                get_density_rescaled(
                    data_source = .,
                    age_table = new_data_general # [config]
                )
        ),
        n2_density = purrr::map(
            .x = diversity_cp,
            .f = ~ .x %>%
                dplyr::filter(var_name == "n2") %>%
                purrr::pluck("age") %>%
                get_density_rescaled(
                    data_source = .,
                    age_table = new_data_general # [config]
                )
        ),
        n1_divided_by_n0_density = purrr::map(
            .x = diversity_cp,
            .f = ~ .x %>%
                dplyr::filter(var_name == "n1_divided_by_n0") %>%
                purrr::pluck("age") %>%
                get_density_rescaled(
                    data_source = .,
                    age_table = new_data_general # [config]
                )
        ),
        n2_divided_by_n1_density = purrr::map(
            .x = diversity_cp,
            .f = ~ .x %>%
                dplyr::filter(var_name == "n2_divided_by_n1") %>%
                purrr::pluck("age") %>%
                get_density_rescaled(
                    data_source = .,
                    age_table = new_data_general # [config]
                )
        )
    )

data_cp_density_merge <-
    data_cp_density %>%
    dplyr::mutate(
        pap_density = purrr::pmap(
            .l = list(
                mvrt_cp_density, # ..1
                roc_cp_density, # ..2
                dcca_cp_density, # ..3
                dca_cp_density, # ..4
                n0_density, # ..5
                n1_density, # ..6
                n2_density, # ..7
                n1_divided_by_n0_density, # ..8
                n2_divided_by_n1_density # ..9
            ),
            .f = ~ new_data_general %>% # [config]
                dplyr::left_join(
                    ..1 %>%
                        dplyr::select(
                            age,
                            mvrt = density
                        ),
                    by = "age"
                ) %>%
                dplyr::left_join(
                    ..2 %>%
                        dplyr::select(
                            age,
                            roc = density
                        ),
                    by = "age"
                ) %>%
                dplyr::left_join(
                    ..3 %>%
                        dplyr::select(
                            age,
                            dcca = density
                        ),
                    by = "age"
                ) %>%
                dplyr::left_join(
                    ..4 %>%
                        dplyr::select(
                            age,
                            dca = density
                        ),
                    by = "age"
                ) %>%
                dplyr::left_join(
                    ..5 %>%
                        dplyr::select(
                            age,
                            n0 = density
                        ),
                    by = "age"
                ) %>%
                dplyr::left_join(
                    ..6 %>%
                        dplyr::select(
                            age,
                            n1 = density
                        ),
                    by = "age"
                ) %>%
                dplyr::left_join(
                    ..7 %>%
                        dplyr::select(
                            age,
                            n2 = density
                        ),
                    by = "age"
                ) %>%
                dplyr::left_join(
                    ..8 %>%
                        dplyr::select(
                            age,
                            n1_divided_by_n0 = density
                        ),
                    by = "age"
                ) %>%
                dplyr::left_join(
                    ..9 %>%
                        dplyr::select(
                            age,
                            n2_divided_by_n1 = density
                        ),
                    by = "age"
                )
        )
    )

#--------------------------------------------------------#
# 4. Save ----
#--------------------------------------------------------#

data_density <-
    data_cp_density_merge  %>% 
    dplyr::select(dataset_id, pap_density)

readr::write_rds(
    data_density,
    here::here(
        "Data/Processed/Partitions/PAP_density_2022-09-15.rds"
    ),
    compress = "gz"
)
