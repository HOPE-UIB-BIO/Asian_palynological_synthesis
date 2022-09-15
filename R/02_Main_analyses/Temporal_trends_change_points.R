#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#             Temporal changes of chnage points
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# Plot the temporal trends of densities of change points of the estimates of
#  Hill's diversity, DCCA1, and MRT.

#--------------------------------------------------------#
# 1. Setup  ----
#--------------------------------------------------------#

library(here)

# Load configuration
source(
    here::here("R/00_Config_file.R")
)

#--------------------------------------------------------#
# 2. Load the data and restrucutre ----
#--------------------------------------------------------#

data_change_points <-
    readr::read_rds(
        here::here(
            "Data/Processed/Partitions/PAP_change_points_2022-09-14.rds"
        )
    )

data_all_estimates <-
    readr::read_rds(
        here::here(
            "Data/Processed/Data_for_analyses/Data_for_analyses-2022-09-15.rds"
        )
    )

# restructure data
data_cp_per_seq <-
    data_all_estimates %>%
    dplyr::select(
        Climate_zone, dataset_id, pap_density
    ) %>%
    tidyr::unnest(pap_density) %>%
    add_age_bin(
        bin_size = bin_size # [config]
    ) %>%
    dplyr::filter(BIN >= 0) %>%
    tidyr::pivot_longer(
        cols = -c(Climate_zone, dataset_id, age, BIN),
        names_to = "var_name",
        values_to = "var"
    ) %>%
    # Set desired order of the facets
    dplyr::mutate(
        var_name = dplyr::case_when(
            var_name == "n0" ~ "N0",
            var_name == "n1" ~ "N1",
            var_name == "n2" ~ "N2",
            var_name == "n2_divided_by_n1" ~ "N2 divided by N1",
            var_name == "n1_divided_by_n0" ~ "N1 divided by N0",
            var_name == "roc" ~ "RoC",
            var_name == "mvrt" ~ "MVRT",
            var_name == "dcca" ~ "DCCA1",
            var_name == "dca" ~ "DCA1",
            TRUE ~ var_name
        )
    ) %>%
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
                "MVRT"
            )
        )
    )

#--------------------------------------------------------#
# 3. Plot per sequence ----
#--------------------------------------------------------#

# Data for boxplot is binned at 1000 years
data_cp_boxplot <-
    data_cp_per_seq %>%
    dplyr::group_by(Climate_zone, dataset_id, var_name, BIN) %>%
    dplyr::summarise(
        .groups = "drop",
        var = mean(var, na.rm = TRUE)
    )

plot_sequence_boxplot <-
    data_cp_boxplot %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = BIN,
            y = var
        )
    ) +
    ggplot2::geom_line(
        data = data_cp_per_seq,
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
        ncol = 2
    ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(
        "Temporal trends of change-points"
    ) +
    ggplot2::labs(
        x = "Time (yr BP)",
        y = "Density of chnage points"
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
        "Outputs/Figures/Temporal_cp_squences.tiff"
    ),
    width = 18,
    height = 24,
    units = "cm",
    dpi = 100,
    compression = "lzw"
)


#--------------------------------------------------------#
# 4. Estimate density per ecozone ----
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
        dcca_cp = list(c(unlist(dcca_cp))),
        dca_cp = list(c(unlist(dca_cp))),
    )


data_cp_density_per_ecozone <-
    get_density_for_all_vars(
        data_source = data_change_points_per_ecozone,
        age_table = new_data_general # [config]
    )

data_cp_density_per_ecozone_restructure <-
    data_cp_density_per_ecozone %>%
    dplyr::select(Climate_zone, pap_density) %>%
    tidyr::unnest(pap_density) %>%
    tidyr::pivot_longer(
        cols = -c(Climate_zone, age),
        names_to = "var_name",
        values_to = "var"
    ) %>%
    # Set desired order of the facets
    dplyr::mutate(
        var_name = dplyr::case_when(
            var_name == "n0" ~ "N0",
            var_name == "n1" ~ "N1",
            var_name == "n2" ~ "N2",
            var_name == "n2_divided_by_n1" ~ "N2 divided by N1",
            var_name == "n1_divided_by_n0" ~ "N1 divided by N0",
            var_name == "roc" ~ "RoC",
            var_name == "mvrt" ~ "MVRT",
            var_name == "dcca" ~ "DCCA1",
            var_name == "dca" ~ "DCA1",
            TRUE ~ var_name
        )
    ) %>%
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
                "MVRT"
            )
        )
    )


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

plot_ecozone <-
    data_cp_density_per_ecozone_restructure %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = age,
            y = var
        )
    ) +
    ggplot2::geom_ribbon(
        ggplot2::aes(
            ymin = 0,
            ymax = var,
            fill = Climate_zone
        )
    ) +
    ggplot2::facet_grid(
        Climate_zone ~ var_name,
        labeller = ggplot2::labeller(
            .cols = ggplot2::label_wrap_gen(10),
            .rows = ggplot2::label_wrap_gen(10)
        )
    ) +
    ggplot2::scale_color_manual(values = ecozone_pallete) +
     ggplot2::scale_fill_manual(values = ecozone_pallete) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(
        "Temporal trends of change-points"
    ) +
    ggplot2::labs(
        x = "Time (yr BP)",
        y = "Density of change points"
    ) +
    ggplot2::theme(
        legend.position = "none",
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
    plot_ecozone,
    filename = here::here(
        "Outputs/Figures/Temporal_cp_ecozone.tiff"
    ),
    width = 24,
    height = 18,
    units = "cm",
    dpi = 400,
    compression = "lzw"
)
