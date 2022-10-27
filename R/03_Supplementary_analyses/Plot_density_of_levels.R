#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                 Estimate density of levels
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# This script prepares fifures with densities of levels

#----------------------------------------------------------#
# 1. Set up  -----
#----------------------------------------------------------#

library(here)

# Load configuration
source(
    here::here("R/00_Config_file.R")
)

#----------------------------------------------------------#
# 2. Load and prepare data -----
#----------------------------------------------------------#

# Input data is the full dataset of Asia with 209 sites.
data_meta <-
    readr::read_rds(
        here::here(
            "Data/Processed/Metadata/Metadata-2022-09-19.rds"
        )
    )

data_combine_paps <-
    readr::read_rds(
        here::here(
            "Data/Processed/PAP_all/pap_all_2022-09-29.rds"
        )
    )

data_levels <-
    dplyr::left_join(
        data_meta,
        data_combine_paps,
        by = "dataset_id"
    ) %>%
    dplyr::select(Climate_zone, dataset_id, levels) %>%
    tidyr::unnest(levels) %>%
    dplyr::filter(
        age >= age_min & age < age_max # [config_criteria]
    )

#----------------------------------------------------------#
# 3. Estimate densities -----
#----------------------------------------------------------#

# helper function
get_density <-
    function(data_source) {
        stats::density(
            x = data_source,
            bw = "nrd0",
            adjust = 1,
            kernel = "gaussian"
        ) %>%
            return()
    }

# - 3.1 sequences level -----
data_density_sequence <-
    data_levels %>%
    tidyr::nest(data_nest = -dataset_id) %>%
    dplyr::mutate(
        density = purrr::map(
            .x = data_nest,
            .f = ~ .x %>%
                purrr::pluck("age") %>%
                get_density()
        )
    ) %>%
    dplyr::mutate(
        data_nest = purrr::map(
            .x = density,
            .f = ~ tibble::tibble(
                age = .x$x,
                var = .x$y
            )
        )
    ) %>%
    dplyr::select(dataset_id, data_nest) %>%
    tidyr::unnest(data_nest) %>%
    dplyr::filter(
        age >= age_min & age < age_max # [config_criteria]
    )

# - 3.2 ecozone level -----
data_density_ecozone <-
    data_levels %>%
    tidyr::nest(data_nest = -Climate_zone) %>%
    dplyr::mutate(
        density = purrr::map(
            .x = data_nest,
            .f = ~ .x %>%
                purrr::pluck("age") %>%
                get_density()
        )
    ) %>%
    dplyr::mutate(
        data_nest = purrr::map(
            .x = density,
            .f = ~ tibble::tibble(
                age = .x$x,
                var = .x$y
            )
        )
    ) %>%
    dplyr::select(Climate_zone, data_nest) %>%
    tidyr::unnest(data_nest) %>%
    dplyr::filter(
        age >= age_min & age < age_max # [config_criteria]
    )


# - 3.3 cotinent level -----
density_continent <-
    data_levels$age %>%
    get_density()

data_density_continent <-
    tibble::tibble(
        age = density_continent$x,
        var = density_continent$y
    ) %>%
    dplyr::filter(
        age >= age_min & age < age_max # [config_criteria]
    )

#----------------------------------------------------------#
# 4. plot -----
#----------------------------------------------------------#

data_for_plotting <-
    tibble::tibble(
        var_name = "density of levels",
        grain = c("sequence", "climate-zone", "continent"),
        var_type = "density",
        y_limits = list(
            c(0, max(data_density_sequence$var)),
            c(0, max(data_density_ecozone$var)),
            c(0, max(data_density_continent$var))
        ),
        merge_data = list(
            data_density_sequence,
            data_density_ecozone,
            data_density_continent
        )
    )

fig_grid_density <-
    plot_temporal_grid(
        data_source = data_for_plotting,
        sel_vars = "density of levels",
        sel_var_type = "density",
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
                        "Outputs/Figures/Density_of_levels"
                    ),
                    ".",
                    .x
                ),
            # 3 grains, +1 is for label
            width = (3 + 1) * set_fig_width,
            # 1 + 1 is for label
            height = (1 + 1) * set_fig_height,
            units = "cm",
            dpi = 400,
            compress = "lzw"
        )
    )
