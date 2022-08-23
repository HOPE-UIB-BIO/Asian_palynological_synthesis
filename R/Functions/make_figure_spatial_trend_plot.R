make_figure_spatial_trend_plot <-
    function(data_source,
             side = c("long", "lat"),
             var_name,
             error_family) {
        model_trend <-
            REcopol::fit_custom_gam(
                x_var = side,
                y_var = var_name,
                error_family = error_family,
                smooth_basis = "tp",
                data_source = data_source,
                sel_k = 10
            )

        x_range <-
            data_source %>%
            purrr::pluck(side) %>%
            range()

        # Predict the models
        data_pred <-
            marginaleffects::predictions(
                model = model_trend,
                newdata = tibble::tibble(
                    !!side := seq(
                        from = min(x_range),
                        to = max(x_range)
                    )
                )
            ) %>%
            tibble::as_tibble() %>%
            dplyr::rename(
                !!var_name := predicted
            )

        y_range <-
            c(
                data_source %>%
                    purrr::pluck(var_name),
                data_pred %>%
                    purrr::pluck(var_name)
            ) %>%
            range()

        bin_size <- 5

        data_hist_raw <-
            data_source %>%
            dplyr::mutate(
                bin = floor(get(side) / bin_size) * bin_size
            ) %>%
            dplyr::group_by(bin) %>%
            dplyr::count() %>%
            dplyr::rename(
                !!side := bin
            )

        data_hist <-
            data_hist_raw %>%
            dplyr::mutate(
                n_rescale = scales::rescale(
                    n,
                    from = range(data_hist_raw$n),
                    to = y_range
                )
            )

        plot_res <-
            data_source %>%
            ggplot2::ggplot(
                ggplot2::aes(
                    x = get(side)
                )
            ) +
            ggplot2::geom_bar(
                data = data_hist,
                ggplot2::aes(y = n_rescale),
                stat = "identity",
                fill = "#2CA388",
                color = "#2CA388",
                size = 0.1,
                alpha = 0.2
            ) +
            ggplot2::geom_ribbon(
                data = data_pred,
                ggplot2::aes(
                    y = get(var_name),
                    ymin = conf.low,
                    ymax = conf.high
                ),
                fill = "#0072B2",
                alpha = 0.25,
                colour = NA
            ) +
            ggplot2::geom_point(
                ggplot2::aes(y = get(var_name)),
                col = "#2CA388",
                alpha = 0.5
            ) +
            ggplot2::geom_line(
                data = data_pred,
                ggplot2::aes(y = get(var_name)),
                size = 1,
                colour = "#0072B2"
            ) +
            ggplot2::theme_classic() +
            ggplot2::theme(
                legend.position = "none",
                axis.title = ggplot2::element_text(color = "black"),
                axis.text = ggplot2::element_text(color = "black")
            ) +
            ggplot2::labs(
                x = "", y = ""
            )

        if (
            side == "lat"
        ) {
            plot_res <-
                plot_res + 
                ggplot2::coord_flip()
        }

        return(plot_res)
    }