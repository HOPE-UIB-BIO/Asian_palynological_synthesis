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
            "Data_for_temporal_plotting-2022-10-28.rds"
        )
    )

# table with y-axis limits
data_limits <-
    data_for_plotting %>%
    dplyr::distinct(var_name) %>%
    dplyr::mutate(
        var_type = "var"
    ) %>%
    dplyr::arrange() %>%
    dplyr::mutate(
        y_limits = list(
            c(0, 50), # N0
            c(0, 20), # N1
            c(0, 15), # N2
            c(0, 3), # DCCA1
            c(0, 1), # N2 divided by N1
            c(0, 1), # N1 divided by N0
            c(0, 1), # RoC,
            c(0, 1), # Peak-points,
            c(0, 1) # MVRT
        ),
        breaks = list(
            c(0, 5, 50), # N0
            c(0, 5, 15), # N1
            c(0, 3, 10), # N2
            c(0, 1.5, 3), # DCCA1
            c(0, 0.5, 1), # N2 divided by N1
            c(0, 0.5, 1), # N1 divided by N0
            c(0, 0.5, 1), # RoC,
            c(0, 0.5, 1), # Peak-points,
            c(0, 0.5, 1) # MVRT
        )
    )

data_trans <-
    data_for_plotting %>%
    dplyr::distinct(var_name) %>%
    dplyr::arrange() %>%
    dplyr::mutate(
        sel_trans = c(
            rep("pseudo_log", 3),
            rep("identity", 6)
        )
    )

data_for_plotting_with_limits <-
    data_for_plotting %>%
    dplyr::left_join(
        data_limits,
        by = c("var_name", "var_type")
    ) %>%
    dplyr::left_join(
        data_trans,
        by = "var_name"
    )

#--------------------------------------------------------#
# 3. Plot all variables in a grid ----
#--------------------------------------------------------#

sel_var_vec <-
    c(
        "N0",
        "N1",
        "N2",
        "N1 divided by N0",
        "N2 divided by N1",
        "DCCA1",
        "RoC",
        "MRT"
    )

# - 3.1 full grid with rescale Y-axis -----

fig_grid_all_rescale <-
    plot_temporal_grid(
        data_source = data_for_plotting_with_limits,
        sel_vars = sel_var_vec,
        sel_var_type = c("var", "density"),
        sel_grain = c("sequence", "climate-zone", "continent"),
        use_limits = TRUE,
        auto_y_breaks = TRUE,
        n_y_ticks = 3,
        use_trans = FALSE,
        plot_rmse = TRUE,
        plot_summary = FALSE,
        def_color = "#2CA388",
        def_of_color = "#A3882C",
        def_violin_color = "#882CA3",
        def_text_size = 16,
        heading_text_multiplier = 2
    )

# View results
grid::grid.newpage()
grid::grid.draw(fig_grid_all_rescale)

# save results
c(
    "tiff",
     "pdf"
) %>%
    purrr::walk(
        .f = ~ ggplot2::ggsave(
            plot = fig_grid_all_rescale,
            filename =
                paste0(
                    here::here(
                        "Outputs/Figures/Temporal_trends_full_rescale"
                    ),
                    ".",
                    .x
                ),
            # 2 types , 3 grains, +1 is for label
            width = ((2 * 3) + 1) * set_fig_width,
            # +1 is for label
            height = (length(sel_var_vec) + 1) * set_fig_height,
            units = "cm",
            dpi = 400,
            compress = "lzw"
        )
    )

# - 3.2 full grid with free Y-axis -----

fig_grid_all_free_y <-
    plot_temporal_grid(
        data_source = data_for_plotting_with_limits,
        sel_vars = sel_var_vec,
        sel_var_type = c("var", "density"),
        sel_grain = c("sequence", "climate-zone", "continent"),
        use_limits = FALSE,
        auto_y_breaks = TRUE,
        n_y_ticks = 3,
        use_trans = FALSE,
        plot_rmse = TRUE,
        plot_summary = FALSE,
        def_color = "#2CA388",
        def_of_color = "#A3882C",
        def_violin_color = "#882CA3",
        def_text_size = 16,
        heading_text_multiplier = 2
    )

# View results
grid::grid.newpage()
grid::grid.draw(fig_grid_all_free_y)

# save results
c(
    "tiff",
     "pdf"
) %>%
    purrr::walk(
        .f = ~ ggplot2::ggsave(
            plot = fig_grid_all_free_y,
            filename =
                paste0(
                    here::here(
                        "Outputs/Figures/Temporal_trends_full_free_y"
                    ),
                    ".",
                    .x
                ),
            # 2 types , 3 grains, +1 is for label
            width = ((2 * 3) + 1) * set_fig_width,
            # +1 is for label
            height = (length(sel_var_vec) + 1) * set_fig_height,
            units = "cm",
            dpi = 400,
            compress = "lzw"
        )
    )


# - 3.3 only var with rescale Y-axis -----

fig_grid_var_rescale <-
    plot_temporal_grid(
        data_source = data_for_plotting_with_limits,
        sel_vars = sel_var_vec,
        sel_var_type = "var",
        sel_grain = c("sequence", "climate-zone", "continent"),
        use_limits = TRUE,
        auto_y_breaks = TRUE,
        n_y_ticks = 3,
        use_trans = FALSE,
        plot_rmse = TRUE,
        plot_summary = FALSE,
        def_color = "#2CA388",
        def_of_color = "#A3882C",
        def_violin_color = "#882CA3",
        def_text_size = 16,
        heading_text_multiplier = 2
    )

# View results
grid::grid.newpage()
grid::grid.draw(fig_grid_var_rescale)

# save results
c(
    "tiff",
     "pdf"
) %>%
    purrr::walk(
        .f = ~ ggplot2::ggsave(
            plot = fig_grid_var_rescale,
            filename =
                paste0(
                    here::here(
                        "Outputs/Figures/Temporal_trends_var_rescale"
                    ),
                    ".",
                    .x
                ),
            # 3 grains, +1 is for label
            width = (3 + 1) * set_fig_width,
            # +1 is for label
            height = (length(sel_var_vec) + 1) * set_fig_height,
            units = "cm",
            dpi = 400,
            compress = "lzw"
        )
    )


# - 3.4 only var with free Y-axis -----

fig_grid_var_free_y <-
    plot_temporal_grid(
        data_source = data_for_plotting_with_limits,
        sel_vars = sel_var_vec,
        sel_var_type = "var",
        sel_grain = c("sequence", "climate-zone", "continent"),
        use_limits = FALSE,
        auto_y_breaks = TRUE,
        n_y_ticks = 3,
        use_trans = FALSE,
        plot_rmse = TRUE,
        plot_summary = FALSE,
        def_color = "#2CA388",
        def_of_color = "#A3882C",
        def_violin_color = "#882CA3",
        def_text_size = 16,
        heading_text_multiplier = 2
    )

# View results
grid::grid.newpage()
grid::grid.draw(fig_grid_var_free_y)

# save results
c(
    "tiff",
     "pdf"
) %>%
    purrr::walk(
        .f = ~ ggplot2::ggsave(
            plot = fig_grid_var_free_y,
            filename =
                paste0(
                    here::here(
                        "Outputs/Figures/Temporal_trends_var_free_y"
                    ),
                    ".",
                    .x
                ),
            # 3 grains, +1 is for label
            width = (3 + 1) * set_fig_width,
            # +1 is for label
            height = (length(sel_var_vec) + 1) * set_fig_height,
            units = "cm",
            dpi = 400,
            compress = "lzw"
        )
    )

# - 3.5 only density with rescale Y-axis -----

fig_grid_density <-
    plot_temporal_grid(
        data_source = data_for_plotting_with_limits,
        sel_vars = sel_var_vec,
        sel_var_type = "density",
        sel_grain = c("sequence", "climate-zone", "continent"),
        use_limits = FALSE,
        auto_y_breaks = TRUE,
        n_y_ticks = 3,
        use_trans = FALSE,
        plot_rmse = TRUE,
        plot_summary = FALSE,
        def_color = "#2CA388",
        def_of_color = "#A3882C",
        def_violin_color = "#882CA3",
        def_text_size = 16,
        heading_text_multiplier = 2
    )

# View results
grid::grid.newpage()
grid::grid.draw(fig_grid_density)

# save results
c(
    "tiff",
     "pdf"
) %>%
    purrr::walk(
        .f = ~ ggplot2::ggsave(
            plot = fig_grid_density,
            filename =
                paste0(
                    here::here(
                        "Outputs/Figures/Temporal_trends_density"
                    ),
                    ".",
                    .x
                ),
            # 3 grains, +1 is for label
            width = (3 + 1) * set_fig_width,
            # +1 is for label
            height = (length(sel_var_vec) + 1) * set_fig_height,
            units = "cm",
            dpi = 400,
            compress = "lzw"
        )
    )
