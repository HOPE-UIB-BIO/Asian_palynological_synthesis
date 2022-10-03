#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                 Temporal trends - estimation
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# Estimate of Hill's diversity, DCCA1, MRT and RoC temporally

#--------------------------------------------------------#
# 1. Setup  ----
#--------------------------------------------------------#

library(here)

# Load configuration
source(
    here::here("R/00_Config_file.R")
)

# change to TRUE if you need to run all the models (takes time!)
run_models <- FALSE

#--------------------------------------------------------#
# 2. Load the data of summary estimates ----
#--------------------------------------------------------#

data_all_estimates <-
    readr::read_rds(
        here::here(
            "Data/Processed/Data_for_analyses/Data_for_analyses-2022-09-29.rds"
        )
    )

data_change_points <-
    readr::read_rds(
        here::here(
            "Data/Processed/Partitions/PAP_change_points_2022-09-29.rds"
        )
    )

# diversity estimates
data_diversity <-
    data_all_estimates %>%
    dplyr::select(
        Climate_zone, dataset_id, PAP_merge
    ) %>%
    tidyr::unnest(PAP_merge) %>%
    add_age_bin(
        bin_size = bin_size # [config]
    ) %>%
    dplyr::filter(BIN >= 0) %>%
    dplyr::mutate(
        # adjust values so small that they are smaller than zero
        dcca_axis_1 = ifelse(dcca_axis_1 < 0, 0, dcca_axis_1)
    )

# RoC estimates
data_roc <-
    data_all_estimates %>%
    dplyr::select(
        Climate_zone, dataset_id, PAP_roc
    ) %>%
    tidyr::unnest(PAP_roc) %>%
    dplyr::rename(
        age = Age
    ) %>%
    add_age_bin(
        bin_size = bin_size # [config]
    ) %>%
    dplyr::filter(BIN >= 0)


#--------------------------------------------------------#
# 3. Variable - temporal trends - individual sequence ----
#--------------------------------------------------------#

# length of each sequences
data_sequence_length <-
    data_diversity %>%
    dplyr::group_by(dataset_id) %>%
    dplyr::summarise(
        age_min = min(age),
        age_max = max(age)
    ) %>%
    add_age_bin(
        .,
        age_var_name = "age_min",
        bin_var_name = "BIN_min",
        bin_size = time_step # [config]
    ) %>%
    add_age_bin(
        .,
        age_var_name = "age_max",
        bin_var_name = "BIN_max",
        bin_size = time_step, # [config]
        sel_method = "forward"
    )


#---------------------#
# - 3.1 fit the models -----
#---------------------#

if (
    run_models == TRUE
) {
    purrr::pwalk(
        .l = list(
            vars_table$var_name, # [config] # ..1
            vars_table$sel_error, # [config] # ..2
            vars_table$sel_data # [config]  ..3
        ),
        .f = ~ {
            var_sel <- ..1
            error_sel <- ..2
            data_sel <- ..3

            message(
                paste("fitting", var_sel)
            )

            # Fit GAM model
            data_mod <-
                REcopol:::fit_multiple_gams(
                    data_source = get(data_sel),
                    x_var = "age",
                    y_var = var_sel,
                    group_var = "dataset_id",
                    smooth_basis = "tp",
                    error_family = error_sel,
                    max_k = max_temporal_k # [config]
                )

            readr::write_rds(
                x = data_mod,
                file = paste0(
                    here::here(
                        "Data/Processed/Models/Per_sequence/"
                    ),
                    "/",
                    var_sel,
                    ".rds"
                ),
                compress = "gz"
            )
        }
    )
}

#---------------------#
# - 3.2 Predict the models -----
#---------------------#

data_per_seq_pred <-
    purrr::map_dfr(
        .x = vars_table$var_name, # [config]
        .f = ~ {
            var_sel <- .x

            message(var_sel)

            data_mod <-
                readr::read_rds(
                    file = paste0(
                        here::here(
                            "Data/Processed/Models/Per_sequence/"
                        ),
                        "/",
                        var_sel,
                        ".rds"
                    )
                )
            # Predict the model
            data_pred <-
                data_mod %>%
                dplyr::mutate(
                    pred_data = purrr::map(
                        .x = mod,
                        .f = ~ REcopol::predic_model(
                            model_source = .x,
                            data_source = new_data_general # [config]
                        ) %>%
                            dplyr::rename(
                                !!var_sel := fit
                            )
                    )
                )
            data_pred %>%
                dplyr::select(dataset_id, pred_data) %>%
                tidyr::unnest(pred_data) %>%
                return()
        }
    )


# restructure the data so all values are in one column and there is a column
#   with the name of the variable
# Filter out prediction! Keep only those prediction which is in the same time
#   bin as the data (do not predict more than 100 years away from data).
data_per_seq_pred_restruct <-
    data_per_seq_pred %>%
    dplyr::select(
        dplyr::any_of(
            c(
                "Climate_zone",
                "dataset_id",
                "age",
                "sd_error",
                "lwr",
                "upr",
                vars_table$var_name # [config]
            )
        )
    ) %>%
    tidyr::pivot_longer(
        cols = vars_table$var_name, # [config]
        names_to = "var_name",
        values_to = "var"
    ) %>%
    tidyr::drop_na(var) %>%
    tidyr::nest(
        data = c(age, lwr, upr, sd_error, var_name, var)
    ) %>%
    dplyr::left_join(
        data_sequence_length,
        by = "dataset_id"
    ) %>%
    dplyr::mutate(
        data_filtered = purrr::pmap(
            .l = list(data, BIN_min, BIN_max),
            .f = ~ ..1 %>%
                dplyr::filter(
                    age >= ..2 & age <= ..3
                )
        )
    ) %>%
    dplyr::select(
        dplyr::any_of(
            c(
                "dataset_id", "data_filtered"
            )
        )
    ) %>%
    tidyr::unnest(data_filtered) %>%
    # Set desired order of the facets
    rename_variable_names() %>%
    tidyr::nest(
        merge_data = -var_name
    )


#--------------------------------------------------------#
# 4. Variable - temporal trends - climate-zone ----
#--------------------------------------------------------#

#---------------------#
# - 4.1 fit the models -----
#---------------------#

if (
    run_models == TRUE
) {
    purrr::pwalk(
        .l = list(
            vars_table$var_name, # [config] # ..1
            vars_table$sel_error, # [config] # ..2
            vars_table$sel_data # [config] # ..3
        ),
        .f = ~ {
            var_sel <- ..1
            error_sel <- ..2
            data_sel <- ..3

            message(
                paste("fitting", var_sel)
            )

            purrr::walk(
                .x = climate_zone_vec, # [config]
                .f = ~ {
                    sel_ecozone <- .x

                    sel_data <-
                        get(data_sel) %>%
                        dplyr::filter(Climate_zone == sel_ecozone) %>%
                        dplyr::mutate(
                            dataset_id = as.factor(dataset_id)
                        )

                    message(sel_ecozone)

                    sel_data$dataset_id %>%
                        unique() %>%
                        length() %>%
                        message()

                    if (
                        nrow(sel_data) > 0
                    ) {
                        # Fit GAM model
                        data_mod <-
                            REcopol::fit_hgam(
                                x_var = "age",
                                y_var = var_sel,
                                group_var = "dataset_id",
                                smooth_basis = "tp",
                                data_source = sel_data,
                                error_family = error_sel,
                                sel_k = max_temporal_k, # [config]
                                common_trend = TRUE,
                                use_parallel = FALSE,
                                max_itiration = 200
                            )


                        readr::write_rds(
                            x = data_mod,
                            file = paste0(
                                here::here(
                                    "Data/Processed/Models/Per_ecozone"
                                ),
                                "/",
                                var_sel,
                                "_",
                                sel_ecozone,
                                ".rds"
                            ),
                            compress = "gz"
                        )

                        return(data_mod)
                    }
                }
            )
        }
    )
}


#---------------------#
# - 4.2 Predict the models -----
#---------------------#

data_per_ecozone_pred <-
    purrr::map2_dfr(
        .x = vars_table$var_name, # [config]
        .y = vars_table$sel_data, # [config]
        .f = ~ {
            var_sel <- .x
            data_sel <- .y

            message(var_sel)

            data_mod <-
                purrr::map(
                    .x = climate_zone_vec, # [config],
                    .f = ~ readr::read_rds(
                        file = paste0(
                            here::here(
                                "Data/Processed/Models/Per_ecozone"
                            ),
                            "/",
                            var_sel,
                            "_",
                            .x,
                            ".rds"
                        )
                    )
                )

            data_ecozone <-
                purrr::map_dfr(
                    .x = climate_zone_vec, # [config]
                    .id = "Climate_zone",
                    .f = ~ {
                        sel_ecozone <- .x

                        message(sel_ecozone)

                        sel_data <-
                            get(data_sel) %>%
                            dplyr::filter(Climate_zone == sel_ecozone) %>%
                            dplyr::mutate(
                                dataset_id = as.factor(dataset_id)
                            )

                        # Predict the models
                        data_pred <-
                            REcopol::predic_model(
                                model_source = data_mod[[sel_ecozone]],
                                data_source = new_data_general %>% # [config]
                                    dplyr::mutate(
                                        dataset_id = sel_data$dataset_id[1]
                                    ),
                                exclude = data_mod[[sel_ecozone]] %>%
                                    gratia::smooths() %>%
                                    stringr::str_subset(., "dataset_id")
                            ) %>%
                            dplyr::rename(
                                !!var_sel := fit
                            )

                        return(data_pred)
                    }
                )
            return(data_ecozone)
        }
    )

data_per_ecozone_pred_restruct <-
    data_per_ecozone_pred %>%
    dplyr::select(
        dplyr::any_of(
            c(
                "Climate_zone",
                "age",
                "lwr",
                "upr",
                vars_table$var_name # [config]
            )
        )
    ) %>%
    tidyr::pivot_longer(
        cols = vars_table$var_name, # [config]
        names_to = "var_name",
        values_to = "var"
    ) %>%
    tidyr::drop_na(var) %>%
    rename_variable_names() %>%
    tidyr::nest(
        merge_data = -var_name
    )

#--------------------------------------------------------#
# 5. Variable - temporal trends - continent ----
#--------------------------------------------------------#

#---------------------#
# - 5.1 fit the models -----
#---------------------#

if (
    run_models == TRUE
) {
    purrr::pwalk(
        .l = list(
            vars_table$var_name, # [config] # ..1
            vars_table$sel_error, # [config] # ..2
            vars_table$sel_data # [config] # ..3
        ),
        .f = ~ {
            var_sel <- ..1
            error_sel <- ..2
            data_sel <- ..3

            message(
                paste("fitting", var_sel)
            )

            # Fit GAM model
            data_mod <-
                REcopol::fit_hgam(
                    x_var = "age",
                    y_var = var_sel,
                    group_var = "dataset_id",
                    smooth_basis = "tp",
                    data_source = get(data_sel),
                    error_family = error_sel,
                    sel_k = max_temporal_k, # [config]
                    common_trend = TRUE,
                    use_parallel = FALSE,
                    max_itiration = 200
                )

            readr::write_rds(
                x = data_mod,
                file = paste0(
                    here::here(
                        "Data/Processed/Models/Per_continent"
                    ),
                    "/",
                    var_sel,
                    ".rds"
                ),
                compress = "gz"
            )
        }
    )
}

#---------------------#
# - 5.2 Predict the models -----
#---------------------#

data_per_continent_pred <-
    purrr::map2_dfr(
        .x = vars_table$var_name, # [config]
        .y = vars_table$sel_data, # [config]
        .f = ~ {
            var_sel <- .x
            data_sel <- .y

            message(var_sel)

            data_mod <-
                readr::read_rds(
                    file = paste0(
                        here::here(
                            "Data/Processed/Models/Per_continent/"
                        ),
                        "/",
                        var_sel,
                        ".rds"
                    )
                )

            sel_data <-
                get(data_sel)

            # Predict the model
            data_pred <-
                REcopol::predic_model(
                    model_source = data_mod,
                    data_source = new_data_general %>% # [config]
                        dplyr::mutate(
                            dataset_id = sel_data$dataset_id[1]
                        ),
                    exclude = data_mod %>%
                        gratia::smooths() %>%
                        stringr::str_subset(., "dataset_id")
                ) %>%
                dplyr::rename(
                    !!var_sel := fit
                )

            return(data_pred)
        }
    )

data_per_continent_pred_restruct <-
    data_per_continent_pred %>%
    dplyr::select(
        dplyr::any_of(
            c(
                "age",
                "lwr",
                "upr",
                vars_table$var_name # [config]
            )
        )
    ) %>%
    tidyr::pivot_longer(
        cols = vars_table$var_name, # [config]
        names_to = "var_name",
        values_to = "var"
    ) %>%
    tidyr::drop_na(var) %>%
    rename_variable_names() %>%
    tidyr::nest(
        merge_data = -var_name
    )


#--------------------------------------------------------#
# 6. Density - temporal trends - sequence ----
#--------------------------------------------------------#

data_cp_density_merge <-
    get_density_for_all_vars(
        data_source = data_change_points,
        age_table = new_data_general # [config]
    ) %>%
    dplyr::select(dataset_id, pap_density)

# restructure data
data_density_per_seq_restruct <-
    data_cp_density_merge %>%
    tidyr::unnest(pap_density) %>%
    add_age_bin(
        bin_size = bin_size # [config]
    ) %>%
    dplyr::filter(BIN >= 0) %>%
    tidyr::pivot_longer(
        cols = -c(dataset_id, age, BIN),
        names_to = "var_name",
        values_to = "var"
    ) %>%
    tidyr::drop_na(var) %>%
    dplyr::group_by(dataset_id, var_name) %>%
    dplyr::mutate(
        var = scales::rescale(var, to = c(0, 1))
    ) %>%
    dplyr::ungroup() %>%
    # Set desired order of the facets
    rename_variable_names() %>%
    tidyr::nest(
        merge_data = -var_name
    )


#--------------------------------------------------------#
# 7. Density - temporal trends - climate-zone ----
#--------------------------------------------------------#

data_change_points_per_ecozone <-
    data_all_estimates %>%
    dplyr::select(Climate_zone, dataset_id) %>%
    dplyr::inner_join(
        data_change_points,
        by = "dataset_id"
    ) %>%
    dplyr::group_by(Climate_zone) %>%
    dplyr::summarise(
        .groups = "drop",
        mvrt_cp = list(c(unlist(mvrt_cp))),
        diversity_cp = list(dplyr::bind_rows(diversity_cp)),
        roc_cp = list(c(unlist(roc_cp))),
        roc_pp = list(c(unlist(roc_pp))),
        dcca_cp = list(c(unlist(dcca_cp)))
    )

data_density_per_ecozone <-
    get_density_for_all_vars(
        data_source = data_change_points_per_ecozone,
        age_table = new_data_general # [config]
    )

data_density_per_ecozone_restruct <-
    data_density_per_ecozone %>%
    dplyr::select(Climate_zone, pap_density) %>%
    tidyr::unnest(pap_density) %>%
    tidyr::pivot_longer(
        cols = -c(Climate_zone, age),
        names_to = "var_name",
        values_to = "var"
    ) %>%
    tidyr::drop_na(var) %>%
    dplyr::group_by(Climate_zone, var_name) %>%
    dplyr::mutate(
        var = scales::rescale(var, to = c(0, 1), from = c(0, max(var)))
    ) %>%
    dplyr::ungroup() %>%
    # Set desired order of the facets
    rename_variable_names() %>%
    tidyr::nest(
        merge_data = -var_name
    )


#--------------------------------------------------------#
# 8. Density - temporal trends - continent ----
#--------------------------------------------------------#

data_change_points_per_continent <-
    data_change_points %>%
    dplyr::summarise(
        .groups = "drop",
        mvrt_cp = list(c(unlist(mvrt_cp))),
        diversity_cp = list(dplyr::bind_rows(diversity_cp)),
        roc_cp = list(c(unlist(roc_cp))),
        roc_pp = list(c(unlist(roc_pp))),
        dcca_cp = list(c(unlist(dcca_cp)))
    )

data_density_per_continent <-
    get_density_for_all_vars(
        data_source = data_change_points_per_continent,
        age_table = new_data_general # [config]
    )

data_density_per_continent_restruct <-
    data_density_per_continent %>%
    dplyr::select(pap_density) %>%
    tidyr::unnest(pap_density) %>%
    tidyr::pivot_longer(
        cols = -c(age),
        names_to = "var_name",
        values_to = "var"
    ) %>%
    tidyr::drop_na(var) %>%
    dplyr::group_by(var_name) %>%
    dplyr::mutate(
        var = scales::rescale(var, to = c(0, 1), from = c(0, max(var)))
    ) %>%
    dplyr::ungroup() %>%
    # Set desired order of the facets
    rename_variable_names() %>%
    tidyr::nest(
        merge_data = -var_name
    )

#--------------------------------------------------------#
# 9. Merge data for plotting together ----
#--------------------------------------------------------#

data_for_plotting <-
    dplyr::bind_rows(
        data_per_seq_pred_restruct %>%
            dplyr::mutate(
                grain = "sequence",
                var_type = "var"
            ),
        data_per_ecozone_pred_restruct %>%
            dplyr::mutate(
                grain = "climate-zone",
                var_type = "var"
            ),
        data_per_continent_pred_restruct %>%
            dplyr::mutate(
                grain = "continent",
                var_type = "var"
            ),
        data_density_per_seq_restruct %>%
            dplyr::mutate(
                grain = "sequence",
                var_type = "density"
            ),
        data_density_per_ecozone_restruct %>%
            dplyr::mutate(
                grain = "climate-zone",
                var_type = "density"
            ),
        data_density_per_continent_restruct %>%
            dplyr::mutate(
                grain = "continent",
                var_type = "density"
            ),
    )

readr::write_rds(
    data_for_plotting,
    file = here::here(
        "Data/Processed/Data_for_temporal_plotting/",
        "Data_for_temporal_plotting-2022-10-03.rds"
    ),
    compress = "gz"
)
