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

# diversity estimates
data_diversity <-
    readr::read_rds(
        here::here(
            "Data/Processed/", "Full_estimates_210410.rds"
        )
    )

# RoC estimates
data_roc <-
    readr::read_rds(
        here::here(
            "Data/Processed/", "Full_roc_estimates_210410.rds"
        )
    ) %>%
    dplyr::rename(
        roc = ROC,
        age = Age
    )

#--------------------------------------------------------#
# 3. Plot the temporal trends of individual sequence ----
#--------------------------------------------------------#

# Values for predicting trends
age_vec <-
    seq(0, 12e3, by = 100)

# data for predicting
new_data_general <-
    tibble::tibble(
        age = age_vec
    )

new_data_dataset <-
    expand.grid(
        age = age_vec,
        unique(data_diversity$dataset_id)
    ) %>%
    tibble::as_tibble()

# tables with names of variables, errors, and dataframes
vars_table_indiv <-
    tibble::tibble(
        var_name = c(
            "N0",
            "N1",
            "N2",
            "DCCA1",
            "N2_divided_by_N1",
            "N1_divided_by_N0",
            "roc"
        ),
        sel_error = c(
            rep("mgcv::Tweedie(p = 1.1)", 4),
            rep("mgcv::betar(link = 'logit')", 2),
            "mgcv::Tweedie(p = 1.1)"
        ),
        sel_data = c(
            rep("data_diversity", 6),
            "data_roc"
        )
    )

vars_table_ecozone <-
    dplyr::bind_rows(
        vars_table_indiv,
        tibble::tibble(
            var_name = c(
                "proportion_of_zones",
                "peak"
            ),
            sel_error = c(
                "mgcv::betar(link = 'logit')",
                "stats::binomial(link = 'logit')"
            ),
            sel_data = c(
                "data_mrt_seq",
                "data_peak_seq"
            )
        )
    )

# MRT partition estimates
n_zones_per_seq <-
    data_diversity %>%
    dplyr::select(dataset_id, BIN, mrt_partiton) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(dataset_id) %>%
    dplyr::summarise(
        .groups = "drop",
        n_zones = unique(mrt_partiton) %>%
            length()
    )

data_mrt_seq <-
    data_diversity %>%
    dplyr::select(Climate_zone, dataset_id, BIN, mrt_partiton) %>%
    dplyr::group_by(Climate_zone, dataset_id, BIN) %>%
    summarise(
        .groups = "drop",
        zones_perbin = unique(mrt_partiton) %>%
            length(),
    ) %>%
    dplyr::left_join(
        n_zones_per_seq,
        by = "dataset_id"
    ) %>%
    dplyr::mutate(
        proportion_of_zones = (zones_perbin / n_zones),
        age = BIN
    )

readr::write_rds(
    data_mrt_seq,
    here::here("Data/Processed/mrt_per_seq.rds")
)

# Peak points
data_peak_seq <-
    data_roc %>%
    dplyr::mutate(
        BIN = floor(age / 1e3) * 1e3
    ) %>%
    dplyr::select(Climate_zone, dataset_id, BIN, Peak) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(Climate_zone, dataset_id, BIN) %>%
    dplyr::summarise(
        .groups = "drop",
        peak = max(Peak),
    ) %>%
    dplyr::rename(
        age = BIN
    )

readr::write_rds(
    data_peak_seq,
    here::here("Data/Processed/peak_points_per_seq.rds")
)



#---------------------#
# 3.1 fit the models -----
#---------------------#

# 3.1.A - using ecozone hGAM -----
# pros - error structure within ecozone -> fewer outliers (good prediction)
# cons - computationally more difucult. Inpose ecozone structure

climate_zone_vec <-
    data_diversity %>%
    purrr::pluck("Climate_zone") %>%
    unique() %>%
    sort() %>%
    rlang::set_names()

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
                    .id = "Climate_zone",
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
                                    common_trend = FALSE,
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
                        "Data/Processed/Models/Per_sequence/Individual_seq_hGAM"
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
    purrr::map2_dfr(
        .x = vars_table_indiv$var_name,
        .y = vars_table_indiv$sel_data,
        .f = ~ {
            var_sel <- .x
            data_sel <- .y

            data_mod <-
                readr::read_rds(
                    file = paste0(
                        here::here(
                            "Data/Processed/Models/Per_sequence/",
                            "Individual_seq_hGAM"
                        ),
                        "/",
                        var_sel,
                        ".rds"
                    )
                )

            message(var_sel)

            data_ecozone <-
                purrr::map_dfr(
                    .x = climate_zone_vec,
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

                        new_data_ecozone <-
                            expand.grid(
                                age = age_vec,
                                dataset_id = unique(sel_data$dataset_id)
                            ) %>%
                            tibble::as_tibble()

                        # Predict the model
                        data_pred <-
                            REcopol::predic_model(
                                model_source = data_mod[[sel_ecozone]],
                                data_source = new_data_ecozone
                            ) %>%
                            dplyr::rename(
                                !!var_sel := fit
                            )

                        return(data_pred)
                    }
                )
        }
    )

# alternative method to fit the models
# 3.1.B - using multiple GAMs -----
# pros - computationally easy, can use RRatepol
# cons - no general error structure -> high outliers
if (
    # it is turn off by default
    FALSE
) {
    if (
        run_models == TRUE
    ) {
        purrr::pmap_walk(
            .l = list(
                vars_table_indiv$var_name, # ..1
                vars_table_indiv$sel_error, # ..2
                vars_table_indiv$sel_data # ..3
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
                            "Data/Processed/Models/Per_sequence/",
                            "Individual_seq_mulitGAM"
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
            .x = vars_table_indiv$var_name,
            .f = ~ {
                var_sel <- .x

                data_mod <-
                    readr::read_rds(
                        file = paste0(
                            here::here(
                                "Data/Processed/Models/Per_sequence/",
                                "Individual_seq_mulitGAM"
                            ),
                            "/",
                            var_sel,
                            ".rds"
                        )
                    )
                # Predict the model
                data_pred <-
                    REcopol::predic_model(
                        model_source = data_mod,
                        data_source = new_data_dataset
                    ) %>%
                    dplyr::rename(
                        !!var_sel := fit
                    )
            }
        )
}

data_per_seq_pred_restruct <-
    data_per_seq_pred %>%
    dplyr::select(
        dplyr::any_of(
            c(
                "Climate_zone",
                "dataset_id",
                "age",
                "conf.low",
                "conf.high",
                vars_table_indiv$var_name
            )
        )
    ) %>%
    tidyr::pivot_longer(
        cols = vars_table_indiv$var_name,
        names_to = "var_name",
        values_to = "var"
    ) %>%
    tidyr::drop_na(var) %>%
    dplyr::filter(
        !(var_name == "roc" & var > 1.5)
    ) %>%
    dplyr::filter(
        !(var_name == "DCCA1" & var > 2)
    ) %>%
    dplyr::filter(
        !(var_name == "N1" & var > 20)
    ) %>%
    # Set desired order of the facets
    dplyr::mutate(
        var_name = dplyr::case_when(
            var_name == "N2_divided_by_N1" ~ "N2 divided by N1",
            var_name == "N1_divided_by_N0" ~ "N1 divided by N0",
            var_name == "roc" ~ "RoC",
            TRUE ~ var_name
        )
    )

#---------------------#
# 3.2 Plot results -----
#---------------------#

# Data for boxplot is binned at 1000 years
data_per_seq_pred_boxplot <-
    data_per_seq_pred_restruct %>%
    dplyr::filter(age %in% seq(0, 12e3, 1e3)) %>%
    dplyr::bind_rows(
        data_peak_seq %>%
            dplyr::select(
                Climate_zone, dataset_id, age,
                var = peak
            ) %>%
            dplyr::mutate(
                var_name = "Peak-points proportion"
            )
    ) %>%
    dplyr::bind_rows(
        data_mrt_seq %>%
            dplyr::select(
                Climate_zone, dataset_id, age,
                var = proportion_of_zones
            ) %>%
            dplyr::mutate(
                var_name = "MRT-zones proportion"
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
    )

data_per_seq_boxplot_median <-
    data_per_seq_pred_boxplot %>%
    dplyr::group_by(var_name, age) %>%
    dplyr::summarise(
        .groups = "drop",
        var_median = median(var)
    )


plot_sequence_boxplot <-
    data_per_seq_pred_boxplot %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = age,
            y = var
        )
    ) +
    ggplot2::geom_line(
        data = data_per_seq_pred_restruct %>%
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
            ),
        ggplot2::aes(
            group = dataset_id
        ),
        color = "#2CA388",
        alpha = 0.4
    ) +
    ggplot2::geom_violin(
        ggplot2::aes(group = age),
        col = NA,
        fill = "#993300",
        alpha = 0.5,
    ) +
    ggplot2::geom_boxplot(
        ggplot2::aes(group = age),
        col = "black",
        fill = "white",
        alpha = 1,
        width = 250,
        outlier.shape = NA
    ) +
    ggplot2::geom_line(
        data = data_per_seq_boxplot_median,
        ggplot2::aes(
            y = var_median
        ),
        size = 1,
        color = "black"
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
                "conf.low",
                "conf.high",
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
            ymin = conf.low,
            ymax = conf.high,
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
