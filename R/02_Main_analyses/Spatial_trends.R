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

data_spatial_trends <-
    readr::read_rds(
        here::here(
            "Data/Processed/Data_for_analyses/Data_for_analyses-2022-09-15.rds"
        )
    ) %>%
    dplyr::select(
        dataset_id,
        long, lat,
        Climate_zone,
        dcca_grad_length,
        dca_grad_length,
        mvrt_groups_n
    )

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


#--------------------------------------------------------#
# 3. Plot the diagrams ----
#--------------------------------------------------------#

#--------------------------------------------------------#
# A. Base map with 5 modified Köppen-Geiger climate zones (in Beck et al. 2018) ----
#--------------------------------------------------------#

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



# Base map (modified climate zones)
base_map <-
    data_spatial_trends %>%
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
    ggplot2::geom_tile(
        data = raster_df,
        aes(x = x, y = y, fill = Climate_zone),
        inherit.aes = FALSE, alpha = 0.75
    ) +
    ggplot2::scale_fill_manual(
        values = ecozone_pallete_full # [config]
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
        legend.key.size = unit(0.6, "cm"),
        legend.title = ggplot2::element_text(size = 12),
        legend.text = ggplot2::element_text(size = 11),
        axis.title = ggplot2::element_text(
            color = "black",
            size = 14
        ),
        axis.text = ggplot2::element_text(
            colour = "black",
            size = 12
        )
    )

base_map_seq_hist <-
  base_map +
   ggplot2::annotation_custom(
    grob =  ggplot2::ggplotGrob(
      ggpubr::gghistogram(
        data_spatial_trends$long,
        fill = "#2CA388",
        color = "#2CA388",
        binwidth = 5,
        size = 0.1,
        alpha = 0.7
      ) +
        ggpubr::theme_transparent()
    ),
    xmin = 17,
    xmax = 180,
    ymin = -2.8,
    ymax = 13
  ) +
   ggplot2::annotation_custom(
    grob =  ggplot2::ggplotGrob(
      ggpubr::gghistogram(
        data_spatial_trends$lat,
        fill = "#2CA388",
        color = "#2CA388",
        binwidth = 5,
        size = 0.1,
        alpha = 0.7
      ) +
        ggpubr::rotate() +
        ggpubr::theme_transparent()
    ),
    xmin = 18.5,
    xmax = 34,
    ymin = 0.3,
    ymax = 85
  )

#--------------------------------------------------------#
# C. turnover ----
#--------------------------------------------------------#

plot_dcca <-
    plot_spatial_dist(
        data_source = data_spatial_trends,
        base_map = base_map_seq_hist,
        var_name = "dcca_grad_length",
        lab_name = "DCCA gradient length",
        error_family = "mgcv::Tweedie(p = 1.1, link = 'log')"
    )

ggplot2::ggsave(
    filename = here::here("Outputs/Figures/Spatial_DCCA.tiff"),
    plot = plot_dcca,
    dpi = 400,
    width = 15,
    height = 15,
    units = "cm",
    compress = "lzw"
)

plot_dca <-
    plot_spatial_dist(
        data_source = data_spatial_trends,
        base_map = base_map_seq_hist,
        var_name = "dca_grad_length",
        lab_name = "DCA gradient length",
        error_family = "mgcv::Tweedie(p = 1.1, link = 'log')"
    )

ggplot2::ggsave(
    filename = here::here("Outputs/Figures/Spatial_DCA.tiff"),
    plot = plot_dca,
    dpi = 400,
    width = 15,
    height = 15,
    units = "cm",
    compress = "lzw"
)

#--------------------------------------------------------#
# D.  MVRT partitions 
#--------------------------------------------------------#

plot_mrt <-
    plot_spatial_dist(
        data_source = data_spatial_trends,
        base_map = base_map_seq_hist,
        var_name = "mvrt_groups_n",
        lab_name = "MRT groups",
        error_family = "stats::poisson(link = 'log')"
    )

ggplot2::ggsave(
    filename = here::here("Outputs/Figures/Spatial_MRT.tiff"),
    plot = plot_mrt,
    dpi = 400,
    width = 15,
    height = 15,
    units = "cm",
    compress = "lzw"
)

#--------------------------------------------------------#
# E.Sequences per ecozone ----
#--------------------------------------------------------#

plot_ecozone_counts <-
    data_spatial_trends %>%
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
    ggplot2::scale_fill_manual(
        values = ecozone_pallete # [config]
    ) +
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
