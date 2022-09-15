#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                     Temporal trends
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# Plot the estimates of Hill's diversity, DCCA1, MRT and RoC temporally

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
            "Data/Processed/Data_for_analyses/Data_for_analyses-2022-09-15.rds"
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
        dcca_axis_1 = ifelse(dcca_axis_1 < 0, 0, dcca_axis_1),
        dca_axis_1 = ifelse(dca_axis_1 < 0, 0, dca_axis_1)
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
# 3. Plot the temporal trends of individual sequence ----
#--------------------------------------------------------#

new_data_dataset <-
    expand.grid(
        age = age_vec,
        dataset_id = unique(data_all_estimates$dataset_id)
    ) %>%
    tibble::as_tibble()

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
        bin_var_name = "BIN_min"
    ) %>%
    add_age_bin(
        .,
        age_var_name = "age_max",
        bin_var_name = "BIN_max"
    )

# vector with names of climate zones
climate_zone_vec <-
    data_all_estimates %>%
    purrr::pluck("Climate_zone") %>%
    unique() %>%
    sort() %>%
    rlang::set_names()

# tables with names of variables, errors, and dataframes
vars_table <-
    tibble::tibble(
        var_name = c(
            "n0",
            "n1",
            "n2",
            "dcca_axis_1",
            "dca_axis_1",
            "n2_divided_by_n1",
            "n1_divided_by_n0",
            "ROC",
            "Peak"
        ),
        sel_error = c(
            rep("mgcv::Tweedie(p = 1.1)", 5),
            rep("mgcv::betar(link = 'logit')", 2),
            "mgcv::Tweedie(p = 1.1)",
            "stats::quasibinomial(link = 'logit')"
        ),
        sel_data = c(
            rep("data_diversity", 7),
            rep("data_roc", 2)
        )
    )


#---------------------#
# 3.1 fit the models -----
#---------------------#

# 3.1.A - using multiple GAMs -----
if (
    run_models == TRUE
) {
    purrr::pwalk(
        .l = list(
            vars_table$var_name, # ..1
            vars_table$sel_error, # ..2
            vars_table$sel_data # ..3
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
                    max_k = 10
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
                )
            )
        }
    )
}

data_per_seq_pred <-
    purrr::map_dfr(
        .x = vars_table$var_name,
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
                            data_source = new_data_general
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
#   bin as the data (do not predict more than 1000 years away from data).
#   In addtion, filer out all prediction with too high uncertaity (sd > 2)
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
                vars_table$var_name
            )
        )
    ) %>%
    tidyr::pivot_longer(
        cols = vars_table$var_name,
        names_to = "var_name",
        values_to = "var"
    ) %>%
    tidyr::drop_na(var) %>%
    # remove the high uncertainty
    dplyr::filter(sd_error < 2) %>%
    add_age_bin(.,
        bin_size = bin_size # [config]
    ) %>%
    tidyr::nest(
        data = c(age, BIN, lwr, upr, sd_error, var_name, var)
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
                    BIN >= ..2 & BIN <= ..3
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
    dplyr::mutate(
        var_name = dplyr::case_when(
            var_name == "n0" ~ "N0",
            var_name == "n1" ~ "N1",
            var_name == "n2" ~ "N2", 
            var_name == "n2_divided_by_n1" ~ "N2 divided by N1",
            var_name == "n1_divided_by_n0" ~ "N1 divided by N0",
            var_name == "ROC" ~ "RoC",
            var_name == "Peak" ~ "Peak-points",
            var_name == "dcca_axis_1" ~ "DCCA1",
            var_name == "dca_axis_1" ~ "DCA1",
            TRUE ~ var_name
        )
    )  %>% 
    dplyr::mutate(
        var_name = factor(var_name,
            levels = c(
                "DCCA1",
                "DCA1",
                "N0",
                "N1",
                "N2",
                "N2 divided by N1",
                "N1 divided by N0",
                "RoC",
                "Peak-points"
            )
        )
    )

#---------------------#
# 3.2 Plot results -----
#---------------------#

# Data for boxplot is binned at 1000 years
data_per_seq_pred_boxplot <-
    data_per_seq_pred_restruct %>%
    dplyr::group_by(dataset_id, var_name, BIN) %>%
    dplyr::summarise(
        .groups = "drop",
        var = mean(var)
    )

plot_sequence_boxplot <-
    data_per_seq_pred_boxplot %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = BIN,
            y = var
        )
    ) +
    ggplot2::geom_line(
        data = data_per_seq_pred_restruct,
        ggplot2::aes(
            x = age,
            group = dataset_id
        ),
        color = "#2CA388",
        alpha = 0.4
    ) +
    ggplot2::geom_violin(
        ggplot2::aes(group = BIN),
        col = NA,
        fill = "#993300",
        alpha = 0.5,
    ) +
    ggplot2::geom_boxplot(
        ggplot2::aes(group = BIN),
        col = "black",
        fill = "white",
        alpha = 1,
        width = 250,
        outlier.shape = NA
    ) +
    ggplot2::facet_wrap(
        ~var_name,
        scales = "free_y",
        ncol = 2
    ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(
        "Temporal diversity trends"
    ) +
    ggplot2::labs(
        x = "Time (yr BP)",
        y = "Estimate"
    ) +
    ggplot2::theme(
        strip.text.x = ggplot2::element_text(size = 16),
        plot.title = ggplot2::element_text(size = 16),
        axis.text.x = ggplot2::element_text(
            color = "black",
            size = 16,
            angle = 45,
            hjust = 1
        ),
        axis.text.y = ggplot2::element_text(color = "black", size = 16),
        axis.title = ggplot2::element_text(color = "black", size = 20),
        plot.margin = ggplot2::margin(0, 0.5, 0, 0, "cm") # t, r, b, l
    )

ggplot2::ggsave(
    plot_sequence_boxplot,
    filename = here::here(
        "Outputs/Figures/Temporal_squences.tiff"
    ),
    width = 18,
    height = 24,
    units = "cm",
    dpi = 100,
    compression = "lzw"
)

#--------------------------------------------------------#
# 4. Plot the temporal trends of climate-zone ----
#--------------------------------------------------------#

# Fit GAM models for each variable with a single common smoother plus group-level
#  smothers with diffring wiggliness (see Pedersen et al. 2019. Hierarchical
#  generalized additive models in ecology: an introduction with mgcv).

# 1. Explicitely include a random effect for the intercept (bs = "re" term).

# 2. Specify "m = 1" instead of "m = 2" for the group-level smoothers, which means
#  the marginal thin plate regression splines (TPRS) basis for this term will
#  penalize the squared  first derivative of the function, rather than the second
#  derivative. This, also, reduces colinearity between the global smoother and the
#  group-specific terms which occasionally leads to high uncertainty around the
#  global smoother. TPRS with m = 1 have a more restricted null space than m = 2
#  smoothers, so should not be as collinear with the global smoother.

# 3. Use restricted maximum likelihood (REML) [method = "REML"] to estimate model
#  coefficients and smoothing parameters, which gives a resonable fit to the data.

# 4. In the factorial random effect smoother (bs="re"), 'k' equals number of levels
#  in the the grouping variable. There are 5 climate zones, so k = 5.
#------------------------------------------------------------------------------#

#---------------------#
# 4.1 fit the models -----
#---------------------#

if (
    run_models == TRUE
) {
    purrr::pwalk(
        .l = list(
            vars_table_ecozone$var_name, # ..1
            vars_table_ecozone$sel_error, # ..2
            vars_table_ecozone$sel_data # ..3
        ),
        .f = ~ {
            var_sel <- ..1
            error_sel <- ..2
            data_sel <- ..3

            message(
                paste("fitting", var_sel)
            )

            data_ecozone <-
                purrr::map(
                    .x = climate_zone_vec,
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
                                    sel_k = 10,
                                    common_trend = TRUE,
                                    use_parallel = FALSE,
                                    max_itiration = 200
                                )

                            return(data_mod)
                        }
                    }
                )

            readr::write_rds(
                x = data_ecozone,
                file = paste0(
                    here::here(
                        "Data/Processed/Models/Per_ecozone"
                    ),
                    "/",
                    var_sel,
                    ".rds"
                )
            )
        }
    )
}

data_per_ecozone_pred <-
    purrr::map2_dfr(
        .x = vars_table_ecozone$var_name,
        .y = vars_table_ecozone$sel_data,
        .f = ~ {
            var_sel <- .x
            data_sel <- .y

            message(var_sel)

            data_mod <-
                readr::read_rds(
                    file = paste0(
                        here::here(
                            "Data/Processed/Models/Per_ecozone/"
                        ),
                        "/",
                        var_sel,
                        ".rds"
                    )
                )

            data_ecozone <-
                purrr::map_dfr(
                    .x = climate_zone_vec,
                    .id = "Climate_zone",
                    .f = ~ {
                        sel_ecozone <- .x

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
                                data_source = new_data_general %>%
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
                vars_table_ecozone$var_name
            )
        )
    ) %>%
    tidyr::pivot_longer(
        cols = vars_table_ecozone$var_name,
        names_to = "var_name",
        values_to = "var"
    ) %>%
    tidyr::drop_na(var)


#---------------------#
# 4.2 Plot the results -----
#---------------------#
# Assign unique colour to each climate zone

ecozone_pallete <-
    c(
        "#FFCC99",
        "#993300",
        "#FF6600",
        "#3399FF",
        "#00CCCC"
    ) %>%
    rlang::set_names(
        nm = c(
            "Arid",
            "Cold_Dry",
            "Cold_Without_dry_season",
            "Polar",
            "Temperate"
        )
    )

plot_estimates_per_ecozone <-
    data_per_ecozone_pred_restruct %>%
    # Set desired order of the facets
    dplyr::mutate(
        var_name = dplyr::case_when(
            var_name == "proportion_of_zones" ~ "MRT-zones proportion",
            var_name == "peak" ~ "Peak-points proportion",
            var_name == "N2_divided_by_N1" ~ "N2 divided by N1",
            var_name == "N1_divided_by_N0" ~ "N1 divided by N0",
            var_name == "roc" ~ "RoC",
            TRUE ~ var_name
        )
    ) %>%
    dplyr::mutate(
        var_name = factor(var_name,
            levels = c(
                "DCCA1",
                "N0",
                "N1",
                "N2",
                "N2 divided by N1",
                "N1 divided by N0",
                "RoC",
                "Peak-points proportion",
                "MRT-zones proportion"
            )
        )
    ) %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = age,
            y = var
        )
    ) +
    ggplot2::geom_ribbon(
        ggplot2::aes(
            ymin = lwr,
            ymax = upr,
            fill = Climate_zone,
        ),
        colour = NA,
        alpha = 0.1
    ) +
    ggplot2::geom_line(
        ggplot2::aes(
            group = Climate_zone,
            colour = Climate_zone
        ),
        size = 1
    ) +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(
        ~var_name,
        ncol = 2,
        scales = "free_y"
    ) +
    ggplot2::labs(
        x = "Time (yr BP)",
        y = "Estimate"
    ) +
    ggplot2::scale_fill_manual(values = ecozone_pallete) +
    ggplot2::scale_color_manual(values = ecozone_pallete) +
    ggplot2::theme(
        strip.text.x = ggplot2::element_text(size = 14),
        plot.title = ggplot2::element_text(
            color = "black",
            size = 18
        ),
        axis.text.x = ggplot2::element_text(
            color = "black",
            size = 16,
            angle = 45,
            hjust = 1
        ),
        axis.text.y = ggplot2::element_text(
            color = "black",
            size = 16
        ),
        axis.title = ggplot2::element_text(
            color = "black",
            size = 18
        ),
        legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"),
        legend.title = ggplot2::element_text(size = 14),
        axis.title.x = ggplot2::element_text(
            margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
        ),
        legend.text = ggplot2::element_text(size = 11)
    ) +
    ggplot2::guides(
        colour = ggplot2::guide_legend(
            title.position = "top",
            nrow = 1,
            byrow = TRUE
        )
    )

ggplot2::ggsave(
    plot_estimates_per_ecozone,
    filename =
        here::here(
            "Outputs/Figures/Temporal_ecozone.tiff"
        ),
    width = 18,
    height = 24,
    units = "cm",
    dpi = 400,
    compress = "lzw"
)
