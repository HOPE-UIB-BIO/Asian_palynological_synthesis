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
            "Data_for_temporal_plotting-2022-10-03.rds"
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
        )
    )

data_for_plotting_with_limits <-
    data_for_plotting %>%
    dplyr::left_join(
        data_limits,
        by = c("var_name", "var_type")
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
        "MVRT"
    )

fig_all_grid <-
    plot_temporal_grid(
        data_source = data_for_plotting_with_limits,
        sel_vars = sel_var_vec,
        use_limits = TRUE,
        plot_rmse = TRUE,
        plot_summary = FALSE,
        def_color = "#2CA388",
        def_of_colo = "#A3882C",
        def_violin_color = "#882CA3",
        def_text_size = 16
    )


# View results
grid::grid.newpage()
grid::grid.draw(fig_all_grid)


#--------------------------------------------------------#
# 4. Save ----
#--------------------------------------------------------#

ggplot2::ggsave(
    plot = fig_all_grid,
    filename =
        here::here(
            "Outputs/Figures/Temporal_trends_wip.tiff"
        ),
    width = ((2 * 3) + 1) * set_fig_width, # 2 types , 3 grains, +1 is for label
    height = (length(sel_var_vec) + 1) * set_fig_height, # +1 is for label
    units = "cm",
    dpi = 400,
    compress = "lzw"
)
