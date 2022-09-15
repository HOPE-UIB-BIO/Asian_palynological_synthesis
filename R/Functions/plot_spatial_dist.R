plot_spatial_dist <-
    function(data_source,
             base_map,
             var_name,
             lab_name,
             error_family) {

        # Scatterplot
        plot_main <-
            base_map +
            ggplot2::geom_point(
                data = data_source,
                ggplot2::aes(size = get(var_name)),
                colour = "#000099",
            ) +
            ggplot2::scale_size_continuous(range = c(0, 5)) +
            ggplot2::theme(
                legend.margin = ggplot2::margin(-0.3, 0.8, 0, 0, unit = "cm")
            ) +
            ggplot2::labs(
                size = lab_name
            )

        plot_main_legend <-
            ggpubr::get_legend(plot_main)

        plot_main_no_legend <-
            plot_main +
            ggpubr::rremove("legend") +
            ggplot2::theme(
                plot.margin = ggplot2::unit(
                    c(0, -0.1, -0.2, -0.5), # t, r, b, l
                    "cm"
                )
            )

        # Latitudinal and longitudinal trends
        data_pred_long <-
            make_figure_spatial_trend(
                data_source = data_source,
                side = "long",
                var_name = var_name,
                error_family = error_family
            )

        data_pred_lat <-
            make_figure_spatial_trend(
                data_source = data_source,
                side = "lat",
                var_name = var_name,
                error_family = error_family
            )

        lat_plot <-
            data_source %>%
            ggplot2::ggplot(
                ggplot2::aes(
                    x = lat,
                    y = get(var_name)
                )
            ) +
              ggplot2::geom_point(
                col = "#2CA388",
                size = 2,
                alpha = 0.8
            ) +
            ggplot2::geom_ribbon(
                data = data_pred_lat,
                ggplot2::aes(
                    ymax = upr,
                    ymin = lwr
                ),
                fill = "#0072B2",
                colour = "NA",
                alpha = 0.75
            ) +
            ggplot2::geom_line(
                data = data_pred_lat,
                size = 1.5,
                colour = "#0072B2",
            ) +
            ggplot2::theme_classic() +
            ggplot2::labs(x = "Latitude", y = lab_name) +
            ggplot2::theme(
                legend.position = "none",
                axis.title = ggplot2::element_text(
                    color = "black",
                    size = 10
                ),
                axis.text = ggplot2::element_text(
                    color = "black",
                    size = 9
                )
            )

        long_plot <-
            data_source %>%
            ggplot2::ggplot(
                ggplot2::aes(
                    x = long,
                    y = get(var_name)
                )
            ) +
              ggplot2::geom_point(
                col = "#2CA388",
                size = 2,
                alpha = 0.8
            ) +
            ggplot2::geom_ribbon(
                data = data_pred_long,
                ggplot2::aes(
                    ymax = upr,
                    ymin = lwr
                ),
                fill = "#0072B2",
                colour = "NA",
                alpha = 0.75
            ) +
            ggplot2::geom_line(
                data = data_pred_long,
                size = 1.5,
                colour = "#0072B2"
            ) +
            ggplot2::theme_classic() +
            ggplot2::labs(x = "Longitude", y = "") +
            ggplot2::theme(
                legend.position = "none",
                axis.title = ggplot2::element_text(
                    color = "black",
                    size = 10
                ),
                axis.text = ggplot2::element_text(
                    color = "black",
                    size = 9
                )
            )

        trends_merged <-
            ggpubr::ggarrange(
                lat_plot,
                long_plot,
                ncol = 2,
                align = "v"
            ) +
            ggplot2::theme(
                plot.margin = ggplot2::unit(
                    c(0, 0.2, -0.2, 0.2), # t, r, b, l
                    "cm"
                )
            )

        final_plot <-
            ggpubr::ggarrange(
                trends_merged,
                plot_main_no_legend,
                nrow = 2,
                ncol = 1,
                heights = c(0.25, 0.7)
            )

        final_plot_legend <-
            ggpubr::ggarrange(
                final_plot,
                plot_main_legend,
                nrow = 2,
                ncol = 1,
                heights = c(0.70, 0.3)
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                panel.border = ggplot2::element_blank()
            )

        return(final_plot_legend)
    }
