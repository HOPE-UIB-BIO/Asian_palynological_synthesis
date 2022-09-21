#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#              Temporal trends - plots
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

set_fig_width <- 6
set_fig_height <- 4

#--------------------------------------------------------#
# 2. Load the data of summary estimates ----
#--------------------------------------------------------#

data_for_plotting <-
    readr::read_rds(
        file = here::here(
            "Data/Processed/Data_for_temporal_plotting",
            "Data_for_temporal_plotting-2022-09-19.rds"
        )
    )

data_with_indiv_plots <-
    data_for_plotting %>%
    dplyr::mutate(
        ifigure = purrr::pmap(
            .l = list(
                var_name, # ..1
                merge_data, # ..2
                scope, # ..3
                what # ..4
            ),
            .f = ~ plot_temporal_trend(
                data_source = ..2,
                var_name = ..1,
                sel_scope = ..3,
                sel_type = ..4
            ) +
                ggpubr::rremove("xylab")
        )
    )

# create all possible combination
dummy_table <-
    data_for_plotting %>%
    tidyr::expand(var_name, scope, what) %>%
    dplyr::mutate(
        name_merge = paste(scope, what, sep = "-")
    ) %>%
    dplyr::mutate(
        order = dplyr::case_when(
            name_merge == "seq-var" ~ 1,
            name_merge == "ecozone-var" ~ 2,
            name_merge == "seq-density" ~ 3,
            name_merge == "ecozone-density" ~ 4,
            TRUE ~ 5
        ),
        var_name = factor(
            var_name,
            levels = c(
                "DCCA1",
                "DCA1",
                "MVRT",
                "RoC",
                "Peak-points",
                "N0",
                "N1",
                "N2",
                "N1 divided by N0",
                "N2 divided by N1"
            )
        )
    ) %>%
    dplyr::arrange(var_name, order)

# create an empty plot
dummy_plot <-
    ggplot2::ggplot() +
    ggplot2::geom_blank() +
    ggplot2::theme_void()

# use dummy table to have all posible figures,
#   replace the missing figures with empty
#   transform all plos to Grobs
data_with_indiv_plots_full <-
    dummy_table %>%
    dplyr::left_join(
        data_with_indiv_plots,
        by = c("var_name", "scope", "what")
    ) %>%
    dplyr::mutate(
        has_a_plot = purrr::map_lgl(
            .x = ifigure,
            .f = ~ !is.null(.x)
        ),
        ifigure_grob = purrr::map2(
            .x = ifigure,
            .y = has_a_plot,
            .f = ~ {
                if (
                    .y == TRUE
                ) {
                    ggplot2::ggplotGrob(.x) %>%
                        return()
                } else {
                    ggplot2::ggplotGrob(dummy_plot) %>%
                        return()
                }
            }
        )
    )

temporal_figure_full <-
    data_with_indiv_plots_full %>%
    dplyr::filter(
        !var_name %in% c("DCA1", "Peak-points")
    ) %>%
    purrr::pluck("ifigure_grob") %>%
    gridExtra::arrangeGrob(
        grobs = .,
        ncol = 4
    )

ggplot2::ggsave(
    plot = temporal_figure_full,
    filename =
        here::here(
            "Outputs/Figures/Temporal_full.tiff"
        ),
    width = ncol(temporal_figure_full) * set_fig_width,
    height = nrow(temporal_figure_full) * set_fig_height,
    units = "cm",
    dpi = 400,
    compress = "lzw"
)

temporal_figure_only_var <-
    data_with_indiv_plots_full %>%
    dplyr::filter(
        !var_name %in% c("DCA1", "Peak-points")
    ) %>%
    dplyr::filter(what == "var") %>%
    dplyr::filter(has_a_plot == TRUE) %>%
    purrr::pluck("ifigure_grob") %>%
    gridExtra::arrangeGrob(
        grobs = .,
        ncol = 2
    )

ggplot2::ggsave(
    plot = temporal_figure_only_var,
    filename =
        here::here(
            "Outputs/Figures/Temporal_only_var.tiff"
        ),
    width = ncol(temporal_figure_only_var) * set_fig_width,
    height = nrow(temporal_figure_only_var) * set_fig_height,
    units = "cm",
    dpi = 400,
    compress = "lzw"
)

temporal_figure_only_dens <-
    data_with_indiv_plots_full %>%
    dplyr::filter(
        !var_name %in% c("DCA1", "Peak-points")
    ) %>%
    dplyr::filter(what == "density") %>%
    dplyr::filter(has_a_plot == TRUE) %>%
    purrr::pluck("ifigure_grob") %>%
    gridExtra::arrangeGrob(
        grobs = .,
        ncol = 2
    )

ggplot2::ggsave(
    plot = temporal_figure_only_dens,
    filename =
        here::here(
            "Outputs/Figures/Temporal_only_dens.tiff"
        ),
    width = ncol(temporal_figure_only_dens) * set_fig_width,
    height = nrow(temporal_figure_only_dens) * set_fig_height,
    units = "cm",
    dpi = 400,
    compress = "lzw"
)

temporal_figure_only_turnover <-
    data_with_indiv_plots_full %>%
    dplyr::filter(
        var_name %in% c(
            "DCCA1",
            "RoC",
            "MVRT"
        )
    ) %>%
    purrr::pluck("ifigure_grob") %>%
    gridExtra::arrangeGrob(
        grobs = .,
        ncol = 4
    )

ggplot2::ggsave(
    plot = temporal_figure_only_turnover,
    filename =
        here::here(
            "Outputs/Figures/Temporal_only_turnover.tiff"
        ),
    width = ncol(temporal_figure_only_turnover) * set_fig_width,
    height = nrow(temporal_figure_only_turnover) * set_fig_height,
    units = "cm",
    dpi = 400,
    compress = "lzw"
)

temporal_figure_only_diversity <-
    data_with_indiv_plots_full %>%
    dplyr::filter(
        var_name %in% c(
            "N0",
            "N1",
            "N2",
            "N1 divided by N0",
            "N2 divided by N1"
        )
    ) %>%
    purrr::pluck("ifigure_grob") %>%
    gridExtra::arrangeGrob(
        grobs = .,
        ncol = 4
    )

ggplot2::ggsave(
    plot = temporal_figure_only_diversity,
    filename =
        here::here(
            "Outputs/Figures/Temporal_only_diversity.tiff"
        ),
    width = ncol(temporal_figure_only_diversity) * set_fig_width,
    height = nrow(temporal_figure_only_diversity) * set_fig_height,
    units = "cm",
    dpi = 400,
    compress = "lzw"
)
