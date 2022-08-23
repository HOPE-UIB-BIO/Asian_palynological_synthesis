plot_spatial_dist <-
    function(data_source,
             base_map,
             var_name,
             lab_name,
             plot_title,
             error_family,
             side_scale_ratio = 0.3) {

        # Scatterplot
        plot_main <-
            base_map +
            ggplot2::geom_point(
                data = data_source,
                size = 6,
                colour = "white",
            ) +
            ggplot2::geom_point(
                data = data_source,
                ggplot2::aes(size = get(var_name)),
                colour = "#000099",
            ) +
            ggplot2::scale_size_continuous(range = c(0, 5)) +
            ggplot2::ggtitle(plot_title) +
            ggplot2::theme(
                legend.margin = ggplot2::margin(-0.2, 0.8, 0, 0, unit = "cm")
            ) +
            ggplot2::labs(
                size = lab_name
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

        plot_main_no_legend <-
            plot_main +
            ggpubr::rremove("legend")

        plot_main_legend <-
            ggpubr::get_legend(plot_main)

        plot_with_sides <-
            plot_main_no_legend +
            ggside::geom_xsidebar(
                data = data_pred_long$data_hist,
                ggplot2::aes(x = long, y = n_rescale),
                stat = "identity",
                fill = "#2CA388",
                color = "#2CA388",
                size = 0.1,
                alpha = 0.2
            ) +
            ggside::geom_xsidepoint(
                data = data_source,
                ggplot2::aes(x = long, y = get(var_name)),
                col = "#2CA388",
                alpha = 0.5
            ) +
            ggside::geom_xsidepath(
                data = data_pred_long$data_trend_pred,
                ggplot2::aes(x = long, y = get(var_name)),
                size = 1,
                colour = "#0072B2"
            ) +
            ggside::scale_xsidey_continuous(
                limits = c(
                    floor(min(data_pred_long$data_hist$n_rescale) * 0.9),
                    ceiling(max(data_pred_long$data_hist$n_rescale) * 1.1)
                )
            ) +
            ggside::geom_ysidebar(
                data = data_pred_lat$data_hist,
                ggplot2::aes(y = lat, x = n_rescale),
                stat = "identity",
                fill = "#2CA388",
                color = "#2CA388",
                size = 0.1,
                alpha = 0.2
            ) +
            ggside::geom_ysidepoint(
                data = data_source,
                ggplot2::aes(y = lat, x = get(var_name)),
                col = "#2CA388",
                alpha = 0.5
            ) +
            ggside::geom_ysidepath(
                data = data_pred_lat$data_trend_pred,
                ggplot2::aes(y = lat, x = get(var_name)),
                size = 1,
                colour = "#0072B2"
            ) +
            ggside::scale_ysidex_continuous(
                limits = c(
                    floor(min(data_pred_lat$data_hist$n_rescale) * 0.9),
                    ceiling(max(data_pred_lat$data_hist$n_rescale) * 1.1)
                )
            ) +
            ggplot2::theme(
                ggside.panel.scale = side_scale_ratio
            )

        final_plot <-
            ggpubr::ggarrange(
                plot_with_sides,
                plot_main_legend,
                nrow = 2,
                ncol = 1,
                heights = c(0.6, 0.4)
            )

        return(final_plot)
    }
