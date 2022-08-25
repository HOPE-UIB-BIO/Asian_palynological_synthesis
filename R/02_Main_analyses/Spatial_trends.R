#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                     Spatial trends
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# Plot the estimates of MRT zones, inertia and turnover as maps


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

summary_diversity <-
    readr::read_rds(
        here::here(
            "Data/Processed", "Summary_estimates_210410.rds"
        )
    ) %>%
    dplyr::filter(!long < 0) %>% # dataset beyond the date-line is removed
    dplyr::filter(long >= 0 & long <= 180)

#--------------------------------------------------------#
# 3. Plot the diagrams ----
#--------------------------------------------------------#

#--------------------------------------------------------#
# A. Base map with 5 modified Köppen-Geiger climate zones (in Beck et al. 2018) ----
#--------------------------------------------------------#

# Read the raster points from the geo-tiff file published in Beck et al. 2018

beck_raster_file <-
    raster::raster(
        here::here(
            "Data/Input/Biomes_spatial/Beck_KG_V1_present_0p083.tif"
        )
    )

# Read the raster value-climatic zone tranlation table
koppen_tranlation_table <-
    readr::read_csv(
        here::here("Data/Input/Biomes_spatial//koppen_link.csv")
    )

# Extract the required raster points

raster_df <-
    # Convert raster points into a dataframe
    as.data.frame(beck_raster_file, xy = TRUE) %>%
    tibble::as_tibble() %>%
    # Extract the rater points of only required area
    dplyr::filter(x > 20 & y > 0) %>% # for required lat, long
    dplyr::rename(raster_values = Beck_KG_V1_present_0p083) %>%
    dplyr::mutate(raster_values = round(raster_values, digits = 0)) %>%
    # Assign the names of climate zone to the raster values
    left_join(koppen_tranlation_table, by = c("raster_values")) %>%
    dplyr::filter(!raster_values == 0) %>%
    dplyr::rename(
        ecozone_koppen_30 = genzone,
        ecozone_koppen_15 = genzone_cluster,
        ecozone_koppen_5 = broadbiome
    ) %>%
    # Modify the climate zones as suggested by John
    dplyr::mutate(
        Climate_zone = dplyr::case_when(
            ecozone_koppen_15 == "Arid_Desert" ~ "Arid",
            ecozone_koppen_15 == "Arid_Steppe" ~ "Arid",
            ecozone_koppen_15 == "Cold_Dry_Summer" ~ "Cold_Dry",
            ecozone_koppen_15 == "Cold_Dry_Winter" ~ "Cold_Dry",
            ecozone_koppen_15 == "Temperate_Dry_Summer" ~ "Temperate",
            ecozone_koppen_15 == "Temperate_Dry_Winter" ~ "Temperate",
            ecozone_koppen_15 == "Temperate_Without_dry_season" ~ "Temperate",
            ecozone_koppen_15 == "Polar_Tundra" ~ "Polar",
            TRUE ~ ecozone_koppen_15
        )
    )

# Assign unique colour to each climate zone
cbPalette1 <-
    c(
        "#FFCC99",
        "#993300",
        # "#CC6600",
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

if (
    FALSE
) {
    # WIP use sf zoom to add country borders
    library(sf)
    library(rnaturalearth)
    worldmap <-
        rnaturalearth::ne_countries(
            scale = "medium",
            type = "map_units",
            returnclass = "sf"
        )

    map_asia <-
        sf::st_crop(
            worldmap,
            ymin = 5,
            ymax = 80,
            xmin = 30,
            xmax = 173
        )
}

# Base map (modified climate zones)
base_map <-
    summary_diversity %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = long,
            y = lat
        )
    ) +
    ggplot2::coord_fixed(
        # ylim = c(5.00, 80.00),
        # xlim = c(30.00, 173.00)
    ) +
    ggplot2::scale_x_continuous(
        limits = c(30, 173),
    ) +
    ggplot2::scale_y_continuous(
        limits = c(5, 80)
    ) +
    ggplot2::geom_tile(
        data = raster_df,
        aes(x = x, y = y, fill = Climate_zone),
        inherit.aes = FALSE, alpha = 0.75
    ) +
    ggplot2::scale_fill_manual(
        values = cbPalette1
    ) +
    ggplot2::labs(
        x = "Longitude",
        y = "Latitude",
        fill = "Climate zones"
    ) +
    ggplot2::theme_classic() +
    # ggplot2::borders(
    #    colour = "black",
    #     size = 0.2
    # ) +
    ggplot2::guides(
        fill = ggplot2::guide_legend(
            nrow = 3,
            byrow = TRUE,
            title.position = "top"
        ),
        size = ggplot2::guide_legend(
            nrow = 1,
            byrow = TRUE,
            title.position = "left"
        )
    ) +
    ggplot2::theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.direction = "horizontal",
        legend.key.size = unit(0.70, "cm"),
        legend.title = ggplot2::element_text(size = 14),
        legend.text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(
            color = "black",
            size = 16
        ),
        axis.text = ggplot2::element_text(
            colour = "black",
            size = 12
        ),
        plot.margin = ggplot2::margin(0, 0, 0.1, 0, "cm")
    ) # t, r, b, l

#--------------------------------------------------------#
# C.DCCA turnover ----
#--------------------------------------------------------#

plot_dcca <-
    plot_spatial_dist(
        data_source = summary_diversity,
        base_map = base_map,
        var_name = "DCCA1",
        lab_name = "DCCA1",
        plot_title = "Compositional turnover (DCCA axis 1)",
        error_family = "mgcv::Tweedie(p = 1.1, link = 'log')",
        side_scale_ratio = 0.25
    )

ggplot2::ggsave(
    filename = here::here("Outputs/Figures/Spatial_DCCA.tiff"),
    plot = plot_dcca,
    dpi = 400,
    width = 20,
    height = 16,
    units = "cm",
    compress = "lzw"
)

#--------------------------------------------------------#
# D. Raw MRT partitions (MRT partitions not divided by number of levels)
#--------------------------------------------------------#

plot_mrt <-
    plot_spatial_dist(
        data_source = summary_diversity,
        base_map = base_map,
        var_name = "mrt_groups_raw",
        lab_name = "MRT groups",
        plot_title = "MRT_partitions",
        error_family = "stats::poisson(link = 'log')",
        side_scale_ratio = 0.25
    )

ggplot2::ggsave(
    filename = here::here("Outputs/Figures/Spatial_MRT.tiff"),
    plot = plot_mrt,
    dpi = 400,
    width = 20,
    height = 16,
    units = "cm",
    compress = "lzw"
)

#--------------------------------------------------------#
# E.Sequences per ecozone ----
#--------------------------------------------------------#

pallete2 <-
    cbPalette1[
        c(
            "Arid",
            "Cold_Dry",
            "Cold_Without_dry_season",
            "Polar",
            "Temperate"
        )
    ]

plot_ecozone_counts <-
    summary_diversity %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = forcats::fct_infreq(Climate_zone),
            fill = Climate_zone
        )
    ) +
    ggplot2::geom_bar() +
    ggplot2::labs(
        x = "Climate zones",
        y = "Number of sequences"
    ) +
    ggplot2::geom_text(
        stat = "count",
        aes(label = ..count..),
        vjust = -0.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::scale_fill_manual(values = pallete2) +
    ggplot2::ggtitle("Sequences in each climate zone") +
    ggplot2::labs(fill = "Modified Köppen-Geiger climate zones", ) +
    ggplot2::theme(
        axis.title.y = ggplot2::element_text(
            color = "black",
            size = 16
        ),
        axis.title.x = ggplot2::element_text(
            color = "black",
            size = 16
        ),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(
            color = "black",
            size = 12,
        ),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = c(0.75, 0.75)
    )

ggplot2::ggsave(
    filename = here::here("Outputs/Figures/Ecozone_counts.tiff"),
    plot = plot_ecozone_counts,
    dpi = 400,
    width = 15,
    height = 15,
    units = "cm",
    compress = "lzw"
)
