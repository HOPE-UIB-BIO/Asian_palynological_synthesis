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

raster_file <-
    raster(paste("Spatial_hope_master_200421/Biomes/Beck_tif/",
        "Beck_KG_V1_present_0p0083.tif",
        sep = ""
    ))

# Read the raster value-climatic zone tranlation table
koppen_tranlation_table <-
    read_csv("Spatial_hope_master_200421/Biomes/Beck_tif/koppen_link.csv")

# Increase the storage capacity
# memory.limit(size = 56000) # no longer supported

# Extract the required raster points

# Reduce the raster dimension (by factor of 10), otherwise would be too big file
raster_file1 <- aggregate(raster_file, fact = 10)

raster_df <-
    # Convert raster points into a dataframe
    as.data.frame(raster_file1, xy = TRUE) %>%
    # Extract the rater points of only required area
    dplyr::filter(x > 20 & y > 0) %>% # for required lat, long
    dplyr::rename(raster_values = Beck_KG_V1_present_0p0083) %>%
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
        ecozone_koppen_15 =
            ifelse(ecozone_koppen_15 == "Arid_Desert", "Arid",
                ecozone_koppen_15
            ),
        ecozone_koppen_15 =
            ifelse(ecozone_koppen_15 == "Arid_Steppe", "Arid",
                ecozone_koppen_15
            ),
        ecozone_koppen_15 =
            ifelse(ecozone_koppen_15 == "Cold_Dry_Summer", "Cold_Dry",
                ecozone_koppen_15
            ),
        ecozone_koppen_15 =
            ifelse(ecozone_koppen_15 == "Cold_Dry_Winter", "Cold_Dry",
                ecozone_koppen_15
            ),
        ecozone_koppen_15 =
            ifelse(ecozone_koppen_15 == "Temperate_Dry_Summer", "Temperate",
                ecozone_koppen_15
            ),
        ecozone_koppen_15 =
            ifelse(ecozone_koppen_15 == "Temperate_Dry_Winter", "Temperate",
                ecozone_koppen_15
            ),
        ecozone_koppen_15 =
            ifelse(ecozone_koppen_15 == "Temperate_Without_dry_season",
                "Temperate", ecozone_koppen_15
            ),
        ecozone_koppen_15 =
            ifelse(ecozone_koppen_15 == "Polar_Tundra", "Polar",
                ecozone_koppen_15
            )
    ) %>%
    dplyr::rename(Climate_zone = ecozone_koppen_15)


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
    )

names(cbPalette1) <-
    c(
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
        ylim = c(5.00, 80.00),
        xlim = c(30.00, 173.00)
    ) +
    # geom_tile(data = raster_df,
    #          aes(x = x, y = y, fill = Climate_zone),
    #         inherit.aes = FALSE, alpha = 0.75) +
    ggplot2::scale_fill_manual(
        values = cbPalette1
    ) +
    ggplot2::labs(
        x = "Longitude",
        y = "Latitude",
        fill = "Climate zones"
    ) +
    ggplot2::theme_classic() +
    ggplot2::borders(
        colour = "black",
        size = 0.2
    ) +
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
# B.Distribution of sequences across ecozones ----
#--------------------------------------------------------#

xbp <-
    ggpubr::gghistogram(
        summary_diversity$long,
        fill = "#2CA388",
        color = "#2CA388",
        binwidth = 5,
        size = 0.1,
        alpha = 0.6
    ) +
    ggpubr::theme_transparent()

ybp <-
    ggpubr::gghistogram(
        summary_diversity$lat,
        fill = "#2CA388",
        color = "#2CA388",
        binwidth = 5,
        size = 0.1,
        alpha = 0.6
    ) +
    ggpubr::rotate() +
    ggpubr::theme_transparent()

xbp_grob <- ggplot2::ggplotGrob(xbp)
ybp_grob <- ggplot2::ggplotGrob(ybp)


#--------------------------------------------------------#
# C.DCCA turnover ----
#--------------------------------------------------------#

dcca_lat <-
    REcopol::fit_custom_gam(
        x_var = "lat",
        y_var = "DCCA1",
        error_family = "mgcv::Tweedie(p = 1.9)",
        smooth_basis = "tp",
        data_source = summary_diversity,
        sel_k = 10
    )

dcca_long <-
    REcopol::fit_custom_gam(
        x_var = "long",
        y_var = "DCCA1",
        error_family = "mgcv::Tweedie(p = 1.9)",
        smooth_basis = "tp",
        data_source = summary_diversity,
        sel_k = 10
    )

# Predict the models
pred_lat_dcca <-
    REcopol::predic_model(
        data_source = with(
            summary_diversity,
            expand.grid(lat = seq(min(lat), max(lat), by = 0.05))
        ),
        model_source = dcca_lat
    ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(
        DCCA1 = var
    )

pred_long_dcca <-
    REcopol::predic_model(
        data_source = with(
            summary_diversity,
            expand.grid(long = seq(min(long), max(long), by = 0.05))
        ),
        model_source = dcca_long
    ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(
        DCCA1 = var
    )

# Plotted individually

# Scatterplot
pl1 <-
    base_map +
    ggplot2::geom_point(aes(size = DCCA1),
        colour = "#000099",
        data = summary_diversity
    ) +
    ggplot2::scale_size_continuous(range = c(0, 5)) +
    ggplot2::ggtitle("Compositional turnover (DCCA axis 1)") +
    ggplot2::theme(
        legend.margin = ggplot2::margin(-0.2, 0.8, 0, 0, unit = "cm")
    )

# Latitudinal and longitudinal trends

lat_plot <-
    summary_diversity %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = lat,
            y = DCCA1
        )
    ) +
    ggplot2::geom_point(
        col = "#2CA388",
        alpha = 0.5
    ) +
    ggplot2::geom_line(
        size = 1,
        colour = "#0072B2",
        data = pred_lat_dcca
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.title = ggplot2::element_text(color = "black"),
        axis.text = ggplot2::element_text(color = "black")
    ) +
    ggplot2::labs(
        x = "Latitude", y = "DCCA1"
    )

long_plot <-
    summary_diversity %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = long,
            y = DCCA1
        )
    ) +
    ggplot2::geom_point(
        col = "#2CA388",
        alpha = 0.5
    ) +
    ggplot2::geom_line(
        size = 1,
        colour = "#0072B2",
        data = pred_long_dcca
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.title = ggplot2::element_text(color = "black"),
        axis.text = ggplot2::element_text(color = "black")
    ) +
    ggplot2::labs(
        x = "Longitude", y = "DCCA1"
    )


main_plot_dcca <-
    pl1 +
    inset(ggplotGrob(lat_plot), xmin = 28, xmax = 69, ymin = 5, ymax = 30) +
    inset(ggplotGrob(long_plot), xmin = 139, xmax = 180, ymin = 5, ymax = 30)


final_plot_dcca <-
    pl1 +
    # Add histograms of number of datasets
    annotation_custom(
        grob = xbp_grob,
        xmin = 18,
        xmax = 190,
        ymin = -2.4,
        ymax = 13
    ) +
    annotation_custom(
        grob = ybp_grob,
        xmin = 19.4,
        xmax = 34,
        ymin = 0.5,
        ymax = 85
    )


ggsave(final_plot_dcca,
    filename =
        paste(
            "07.Outputs_EDA_041021/modified_zones/figures_revised_210622/",
            "03.studyarea_Compositional_turnover_new.tiff",
            sep = ""
        ),
    dpi = 400,
    compress = "lzw"
)


#--------------------------------------------------------#
# D. Raw MRT partitions (MRT partitions not divided by number of levels)
#--------------------------------------------------------#

mrt_lat <-
    REcopol::fit_custom_gam(
        x_var = "lat",
        y_var = "mrt_groups_raw",
        error_family = "stats::poisson(link = 'log')",
        smooth_basis = "tp",
        data_source = summary_diversity,
        sel_k = 10
    )


mrt_long <-
    REcopol::fit_custom_gam(
        x_var = "long",
        y_var = "mrt_groups_raw",
        error_family = "stats::poisson(link = 'log')",
        smooth_basis = "tp",
        data_source = summary_diversity,
        sel_k = 10
    )

pred_lat_mrt <-
    REcopol::predic_model(
        data_source = with(
            summary_diversity,
            expand.grid(lat = seq(min(lat), max(lat), by = 0.05))
        ),
        model_source = mrt_lat
    ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(
        mrt_groups = var
    )

pred_long_mrt <-
    REcopol::predic_model(
        data_source = with(
            summary_diversity,
            expand.grid(long = seq(min(long), max(long), by = 0.05))
        ),
        model_source = mrt_long
    ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(
        mrt_groups = var
    )

# Scatterplot
pl2 <-
    base_map +
    ggplot2::geom_point(
        ggplot2::aes(size = mrt_groups_raw),
        colour = "#000099",
        data = summary_diversity
    ) +
    ggplot2::scale_size_continuous(range = c(1, 6)) +
    ggplot2::ggtitle("MRT_partitions") +
    ggplot2::labs(size = "MRT groups") +
    ggplot2::theme(
        legend.margin = ggplot2::margin(-0.2, 0.8, 0, 0, unit = "cm")
    )


# Latitudinal and longitudinal trends

lat_mrt <-
    summary_diversity %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = lat,
            y = mrt_groups_raw
        )
    ) +
    ggplot2::geom_point(col = "#2CA388", alpha = 0.5) +
    ggplot2::geom_line(
        ggplot2::aes(
            x = lat,
            y = mrt_groups
        ),
        size = 1,
        colour = "#0072B2",
        data = pred_lat_mrt
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.title = ggplot2::element_text(color = "black"),
        axis.text = ggplot2::element_text(color = "black")
    ) +
    ggplot2::labs(
        x = "Latitude", y = "MRT_groups"
    )

long_mrt <-
    summary_diversity %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = long, y = mrt_groups_raw
        )
    ) +
    ggplot2::geom_point(col = "#2CA388", alpha = 0.5) +
    ggplot2::geom_line(
        ggplot2::aes(x = long, y = mrt_groups),
        size = 1,
        colour = "#0072B2",
        data = pred_long_mrt
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.title = ggplot2::element_text(color = "black"),
        axis.text = ggplot2::element_text(color = "black")
    ) +
    ggplot2::labs(
        x = "Longitude",
        y = "MRT_groups"
    )

main_plot_mrt <-
    pl2 +
    inset(ggplotGrob(lat_mrt), xmin = 21.9, xmax = 71, ymin = 0, ymax = 30) +
    inset(ggplotGrob(long_mrt), xmin = 130, xmax = 180, ymin = 0, ymax = 30)

ggsave(main_plot_mrt,
    filename =
        paste(
            "07.Outputs_EDA_041021/modified_zones/figures_revised_210622/",
            "03.MRT_groups.tiff",
            sep = ""
        ),
    dpi = 400,
    compress = "lzw"
)


#--------------------------------------------------------#
# E.Sequences per ecozone ----
#--------------------------------------------------------#
# Sys.setlocale("LC_ALL", "German") # To paste special characters in ggplot properly
# options(encoding = "UTF-8")

pallete2 <-
    cbPalette1[c(
        "Arid",
        "Cold_Dry",
        "Cold_Without_dry_season",
        "Polar",
        "Temperate"
    )]

pl3 <-
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
    paste(
        "07.Outputs_EDA_041021/modified_zones/figures_revised_210622/",
        "02.Sequences_per_climate_zone.tiff",
        sep = ""
    ),
    plot = pl3,
    width = 15,
    height = 15,
    units = "cm",
    compress = "lzw"
)


#--------------------------------------------------------#
# All other analyses are similar to previous version!!!
#--------------------------------------------------------#
