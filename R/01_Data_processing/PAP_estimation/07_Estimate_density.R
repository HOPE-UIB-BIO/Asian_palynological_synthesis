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
# Estimate density ----
#--------------------------------------------------------#
data_cp_density <-
    data_change_points %>%
    dplyr::mutate(
        mvrt_cp_density = purrr::map(
            .x = mvrt_cp,
            .f = ~ REcopol::get_density(
                data_source = .x,
                reflected = TRUE,
                values_range = c(
                    age_min, # [config]
                    age_max # [config]
                ),
                bw = 1000 / age_max, # [config]
                n = age_max
            ) %>%
                dplyr::mutate(
                    age = round(var)
                ) %>%
                dplyr::filter(
                    age %in% age_vec
                )
        )
    )


data_cp_density %>%
    dplyr::sample_n(size = 25) %>%
    dplyr::mutate(
        p_diversity = purrr::map2(
            .x = mvrt_cp_density,
            .y = mvrt_cp,
            .f = ~ .x %>%
                ggplot2::ggplot() +
                ggplot2::geom_line(
                    ggplot2::aes(
                        age, density
                    )
                ) +
                ggplot2::geom_vline(
                    xintercept = .y,
                    colour = "gray80"
                ) +
                ggplot2::coord_cartesian(
                    xlim = range(age_vec), # [config]
                    ylim = c(0, 5)
                )
        )
    ) %>%
    purrr::pluck("p_diversity") %>%
    ggpubr::ggarrange(
        plotlist = .
    )
