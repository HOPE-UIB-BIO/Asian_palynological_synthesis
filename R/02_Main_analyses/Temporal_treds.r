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

number_of_cores <- parallel::detectCores()

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
# 3. Plot the temporal trends in N0, N1, N2 DCCA1, MRT and RoC ----
#--------------------------------------------------------#

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

mrt_seq_data <-
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

#---------------------#
# I Individual sequences
#---------------------#

# 3.1 Fit models -----

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
            rep("mgcv::Tweedie(p = 1.9)", 4),
            rep("mgcv::betar(link = 'logit')", 2),
            "mgcv::Tweedie(p = 1.9)"
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
            var_name = "proportion_of_zones",
            sel_error = "mgcv::betar(link = 'logit')",
            sel_data = "mrt_seq_data",
        )
    )


# IA - using multiple GAMs -----
# pros - computationally easy, can use RRatepol
# cons - no general error structure -> high outliers

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

gam_pred_all <-
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
                marginaleffects::predictions(
                    model = data_mod,
                    newdata = new_data_dataset
                ) %>%
                tibble::as_tibble() %>%
                dplyr::rename(
                    !!var_sel := predicted
                )
        }
    )


# IB - using ecozone hGAM -----
# pros - error structure within ecozone -> fewer outliers (good prediction)
# cons - computationally easier than continent (!)

climate_zone_vec <-
    data_diversity %>%
    purrr::pluck("Climate_zone") %>%
    unique() %>%
    sort() %>%
    rlang::set_names()

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
                    message(nrow(sel_data))

                    if (
                        nrow(sel_data) > 0
                    ) {
                        # Fit GAM model
                        data_mod <-
                            mgcv::bam(
                                get(var_sel) ~
                                    s(age,
                                        by = dataset_id,
                                        bs = "tp", m = 1
                                    ) +
                                    s(dataset_id,
                                        bs = "re",
                                        k = length(unique(sel_data$dataset_id))
                                    ),
                                family = error_sel,
                                method = "fREML",
                                data = sel_data,
                                control = mgcv::gam.control(trace = TRUE)
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

gam_pred_all <-
    purrr::map2_dfr(
        .x = vars_table_indiv$var_name,
        .y = vars_table$sel_data,
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

                        new_data_ecozone <-
                            expand.grid(
                                age = age_vec,
                                dataset_id = unique(sel_data$dataset_id)
                            ) %>%
                            tibble::as_tibble()

                        # Predict the model
                        data_pred <-
                            marginaleffects::predictions(
                                model = data_mod[[sel_ecozone]],
                                newdata = new_data_ecozone
                            ) %>%
                            tibble::as_tibble() %>%
                            dplyr::rename(
                                !!var_sel := predicted
                            )
                        return(data_pred)
                    }
                )
        }
    )

# 3.2 Plot results -----

# Data for boxplot is binned at 1000 years
boxplot_dat <-
    gam_pred_all %>%
    dplyr::filter(age %in% seq(0, 12e3, 1e3)) %>%
    group_by(dataset_id, age, sel_var)

plot_sequence_boxplot <-
    gam_pred_all %>%
    dplyr::mutate(
        # Set desired order of facets
        sel_var = factor(
            sel_var,
            levels = c(
                "DCCA1",
                "N0",
                "N1",
                "N2",
                "N2_divided_by_N1",
                "N1_divided_by_N0",
                "roc"
            )
        )
    ) %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = age,
            y = var
        )
    ) +
    ggplot2::geom_line(
        ggplot2::aes(
            group = dataset_id
        ),
        color = "#2CA388",
        alpha = 0.4
    ) +
    ggplot2::geom_boxplot(
        ggplot2::aes(group = age),
        data = boxplot_dat,
        col = "#993300", fill = "#993300",
        alpha = 0.5,
    ) +
    ggplot2::facet_wrap(
        ~sel_var,
        scales = "free_y",
        nrow = 4
    ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(
        "Temporal diversity trends (without binning the estimates)"
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
    filename =
        paste0(
            "07.Outputs_EDA_041021/modified_zones/figures_revised_210622/",
            "00.Diversity_estimates_composite_continent_no_bins_120822.tiff",
            sep = ""
        ),
    width = 18,
    height = 24,
    units = "cm",
    dpi = 100,
    compression = "lzw"
)


#---------------------#
# II. Plots at the climate-zone scale
#---------------------#


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

climate_zone_vec <-
    data_diversity %>%
    purrr::pluck("Climate_zone") %>%
    unique() %>%
    sort() %>%
    rlang::set_names()

# var_sel <- vars_table_ecozone$var_name[[8]]
# error_sel <- vars_table_ecozone$sel_error[[8]]
# data_sel <- vars_table_ecozone$sel_data[[8]]

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
                            mgcv::bam(
                                get(var_sel) ~
                                    s(age, bs = "tp") +
                                    s(age,
                                        by = dataset_id,
                                        bs = "tp", m = 1
                                    ) +
                                    s(dataset_id,
                                        bs = "re",
                                        k = length(unique(sel_data$dataset_id))
                                    ),
                                family = error_sel,
                                method = "fREML",
                                data = sel_data,
                                control = mgcv::gam.control(trace = TRUE)
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


gam_pred_ecozone <-
    purrr::map2_dfr(
        .x = vars_table_indiv$var_name,
        .y = vars_table$sel_data,
        .f = ~ {
            var_sel <- .x
            data_sel <- .y

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
                            marginaleffects::predictions(
                                model = data_mod[[sel_ecozone]],
                                newdata = new_data_general %>%
                                    dplyr::mutate(
                                        dataset_id = sel_data$dataset_id[1]
                                    ),
                                exclude = gratia::smooths(data_mod) %>%
                                    stringr::str_subset(., "dataset_id")
                            ) %>%
                            tibble::as_tibble() %>%
                            dplyr::rename(
                                !!var_sel := predicted
                            )
                        return(data_pred)
                    }
                )
        }
    )


# END OF REVIEW !


# Plot the estimates

# Assign unique colour to each climate zone

req_pallete <-
    c(
        "#FFCC99",
        "#993300",
        "#FF6600",
        "#3399FF",
        "#00CCCC"
    )

names(req_pallete) <-
    c(
        "Arid",
        "Cold_Dry",
        "Cold_Without_dry_season",
        "Polar",
        "Temperate"
    )

# Set desired order of the facets
combined_mod$Variable <-
    factor(combined_mod$Variable,
        levels = c(
            "DCCA1", "MRT_partition_percentage", "N0", "N1", "N2",
            "N2_divided_by_N1", "N1_divided_by_N0", "RoC"
        )
    )

main_plot <-
    combined_mod %>%
    ggplot(aes(x = age, y = var)) +
    geom_line(aes(
        group = Climate_zone,
        colour = Climate_zone
    ),
    size = 1
    ) +
    theme_classic() +
    facet_wrap(~Variable, ncol = 2, scales = "free_y") +
    labs(
        x = "Time (yr BP)",
        y = "Estimate/1000 yr"
    ) +
    scale_fill_manual(values = req_pallete) +
    scale_color_manual(values = req_pallete) +
    theme(
        strip.text.x = element_text(size = 14),
        plot.title = element_text(color = "black", size = 18),
        axis.text.x = element_text(
            color = "black",
            size = 16,
            angle = 45,
            hjust = 1
        ),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        legend.text = element_text(size = 11)
    ) +
    guides(colour = guide_legend(
        title.position = "top",
        nrow = 1,
        byrow = TRUE
    ))


ggsave(
    main_plot,
    filename =
        paste0(
            "07.Outputs_EDA_041021/modified_zones/figures_revised_210622/",
            "00.Diversity_estimates_composite_perzone_240422.tiff",
            sep = ""
        ),
    width = 18,
    height = 24,
    units = "cm",
    dpi = 400,
    compress = "lzw"
)



#--------------------------------------------------------#
# 4. Plot the temporal trends in RoC peaks ----
#--------------------------------------------------------#
# Detect the BINs with significant peak points for each dataset.
# Within each BIN, determine the number of datasets with at least one peak.
# Plot as histogram

# 4a. At continent level
roc_peak_continent <-
    dat_roc %>%
    dplyr::select(dataset_id,
        climate_zone = Climate_zone, BIN, age = Age,
        roc = ROC, peak = Peak
    ) %>%
    dplyr::mutate(peak = ifelse(peak == FALSE, 0, 1)) %>%
    arrange(dataset_id) %>%
    group_by(dataset_id, BIN) %>%
    summarise(.groups = "keep", peak_per_bin = sum(peak)) %>%
    ungroup() %>%
    dplyr::mutate(peak_presence = ifelse(peak_per_bin > 1, 1, peak_per_bin)) %>%
    group_by(BIN) %>%
    summarise(.groups = "keep", peak_counts = sum(peak_presence)) %>%
    dplyr::mutate(climate_zone = paste("Whole_continent"))



# 4b. At climate-zone level
roc_peak_climate_zone <-
    dat_roc %>%
    dplyr::select(dataset_id,
        climate_zone = Climate_zone, BIN, age = Age,
        roc = ROC, peak = Peak
    ) %>%
    dplyr::mutate(peak = ifelse(peak == FALSE, 0, 1)) %>%
    arrange(dataset_id) %>%
    group_by(dataset_id, BIN, climate_zone) %>%
    summarise(.groups = "keep", peak_per_bin = sum(peak)) %>%
    ungroup() %>%
    dplyr::mutate(peak_presence = ifelse(peak_per_bin > 1, 1, peak_per_bin)) %>%
    group_by(climate_zone, BIN) %>%
    summarise(.groups = "keep", peak_counts = sum(peak_presence))


req_pallete <-
    c(
        "#FFCC99",
        "#993300",
        "#FF6600",
        "#3399FF",
        "#00CCCC",
        "#339999"
    )

names(cbPalette1) <-
    c(
        "Arid",
        "Cold_Dry",
        "Cold_Without_dry_season",
        "Polar",
        "Temperate",
        "Whole_continent"
    )


dat_to_plot <-
    bind_rows(roc_peak_climate_zone, roc_peak_continent) %>%
    ggplot(aes(
        x = BIN, y = peak_counts,
        group = as_factor(climate_zone)
    )) +
    geom_bar(aes(fill = climate_zone),
        stat = "identity",
        position = "identity"
    ) +
    theme_classic() +
    scale_fill_manual(values = req_pallete) +
    facet_wrap(~climate_zone, ncol = 2, scales = "free_y") +
    labs(
        x = "Time (yr BP)",
        y = "Number of daasets with at least one peak in a time-bin",
        title = "No. of datasets with at least one peak in a BIN"
    ) +
    theme(
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(
            color = "black", size = 13,
            angle = 45, hjust = 1
        ),
        axis.text.y = element_text(color = "black", size = 13),
        axis.title = element_text(color = "black", size = 16),
        plot.title = element_text(color = "black", size = 16),
        legend.position = "none"
    )

ggsave(
    dat_to_plot,
    filename =
        paste(
            "07.Outputs_EDA_041021/modified_zones/figures_revised_210622/",
            "RoC_peaks.tiff",
            sep = ""
        ),
    width = 16,
    height = 24,
    units = "cm",
    dpi = 200,
    compress = "lzw"
)
