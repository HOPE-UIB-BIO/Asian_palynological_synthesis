#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#             Temporal changes of variable zones
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# Plot the relative proportion of the zones of the estimates of
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
# 2. Load the data of summary estimates ----
#--------------------------------------------------------#

# diversity estimates
data_diversity_zones <-
    readr::read_rds(
        here::here(
            "Data/Processed/PAP_all/pap_all_20220826"
        )
    ) %>%
    dplyr::select(
        Climate_zone, dataset_id, PAP_diversity
    ) %>%
    tidyr::unnest(PAP_diversity) %>%
    dplyr::select(-depth) %>%
    add_age_bin(
        bin_size = 1e3
    ) %>%
    dplyr::filter(BIN >= 0) %>%
    dplyr::select(
        Climate_zone,
        dataset_id,
        age,
        BIN,
        part_N0,
        part_N1,
        part_N2,
        part_N2_divided_by_N1,
        part_N1_divided_by_N0,
        rpart_dcca1
    )

data_mrt_seq <-
    readr::read_rds(
        here::here("Data/Processed/mrt_per_seq.rds")
    )

data_peak_seq <-
    readr::read_rds(
        here::here("Data/Processed/peak_points_per_seq.rds")
    )

#--------------------------------------------------------#
# 3. Estimate the presence of regression partitions per 1000-yr ----
#--------------------------------------------------------#

# PROCEDURE FOR ESTIMATING THE RELATIVE PERCENTAGE OF MRT PARTITIONS:

# Detect the partitions in each BIN for each dataset.
# Regard the presence of partition in each BIN as 1, and absence as 0.
# Sum the number of partitions in each BIN across all the datasets

# This, each BIN gives and estimate of number of datasets with atleast one
#  partition in that BIN.

data_diversity_zones_long <-
    data_diversity_zones %>%
    tidyr::pivot_longer(
        cols = -c(
            dataset_id,
            Climate_zone,
            age,
            BIN,
        ),
        names_to = "var_name",
        values_to = "zone"
    )

data_zones_per_sequence <-
    data_diversity_zones_long %>%
    split(.$var_name) %>%
    purrr::map(
        .f = ~ .x %>%
            dplyr::group_by(dataset_id) %>%
            dplyr::summarise(
                .groups = "drop",
                n_zones = unique(zone) %>%
                    length()
            )
    ) %>%
    purrr::map_dfr(
        .id = "var_name",
        .f = ~ dplyr::bind_rows(.x)
    )

data_zones_proportion <-
    data_diversity_zones_long %>%
    split(.$var_name) %>%
    purrr::map(
        .f = ~ .x %>%
            dplyr::group_by(var_name, Climate_zone, dataset_id, BIN) %>%
            dplyr::summarise(
                .groups = "drop",
                n_zones_per_bin = unique(zone) %>%
                    length()
            ) %>%
            dplyr::left_join(
                data_zones_per_sequence,
                by = c("dataset_id", "var_name")
            ) %>%
            dplyr::mutate(
                zones_proportion = n_zones_per_bin / n_zones
            )
    ) %>%
    purrr::map_dfr(
        .id = "var_name",
        .f = ~ dplyr::bind_rows(.x)
    ) %>%
    dplyr::bind_rows(
        data_mrt_seq %>%
            dplyr::mutate(
                var_name = "MRT-zones",
            ) %>%
            dplyr::select(
                Climate_zone,
                dataset_id,
                BIN,
                var_name,
                zones_proportion = proportion_of_zones
            )
    ) %>%
    dplyr::bind_rows(
        data_peak_seq %>%
            dplyr::mutate(
                var_name = "Peak-points"
            ) %>%
            dplyr::rename(
                zones_proportion = peak,
                BIN = age
            )
    ) %>%
    # Set desired sequence of the facets
    dplyr::mutate(
        var_name = dplyr::case_when(
            var_name == "part_N0" ~ "N0",
            var_name == "part_N1" ~ "N1",
            var_name == "part_N2" ~ "N2",
            var_name == "part_N2_divided_by_N1" ~ "N2 divided by N1",
            var_name == "part_N1_divided_by_N0" ~ "N1 divided by N0",
            var_name == "rpart_dcca1" ~ "DCCA1",
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
                "Peak-points",
                "MRT-zones"
            )
        )
    )

data_zone_median_per_bin <-
    data_zones_proportion %>%
    dplyr::group_by(var_name, Climate_zone, BIN) %>%
    dplyr::summarise(
        .groups = "drop",
        var_median = median(zones_proportion)
    )

# Assign unique colour to each climate zone
palette_ecozones <-
    c(
        "#FFCC99",
        "#993300",
        "#FF6600",
        "#3399FF",
        "#999999",
        "#00CCCC",
        "#99CC00",
        "#006600",
        "#996600"
    ) %>%
    rlang::set_names(
        nm = c(
            "Arid",
            "Cold_Dry",
            "Cold_Without_dry_season",
            "Polar",
            "Polar_Frost",
            "Temperate",
            "Tropical_Monsoon",
            "Tropical_Rainforest",
            "Tropical_Savannah"
        )
    )

plot_proportions_per_ecozone <-
    data_zones_proportion %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = BIN,
            y = zones_proportion
        )
    ) +
    ggplot2::geom_violin(
        ggplot2::aes(
            group = BIN,
            fill = Climate_zone
        ),
        col = NA,
        alpha = 0.75,
    ) +
    ggplot2::geom_boxplot(
        ggplot2::aes(
            group = BIN
        ),
        col = "black",
        fill = "white",
        outlier.shape = NA,
        alpha = 1,
        width = 250
    ) +
    ggplot2::geom_line(
        data = data_zone_median_per_bin,
        ggplot2::aes(
            y = var_median
        ),
        size = 1,
        color = "black"
    ) +
    ggplot2::theme_classic() +
    ggplot2::facet_grid(
        var_name ~ Climate_zone
    ) +
    ggplot2::scale_fill_manual(
        values = palette_ecozones
    ) +
    ggplot2::scale_y_continuous(
        breaks = c(0, 0.5, 1)
    ) +
    ggplot2::ggtitle("Proportion of zones per 1000-yrs") +
    ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(size = 14, colour = "black"),
        axis.text = ggplot2::element_text(size = 13, colour = "black"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.title = ggplot2::element_text(size = 16, colour = "black"),
        plot.title = ggplot2::element_text(size = 16, colour = "black")
    ) +
    ggplot2::labs(
        x = "Time (yr BP)",
        y = "Proportion"
    )

ggsave(
    plot_proportions_per_ecozone,
    filename = here::here(
        "Outputs/Figures/Temporal_proportions.tiff"
    ),
    width = 24,
    height = 20,
    units = "cm",
    dpi = 400,
    compress = "lzw"
)
