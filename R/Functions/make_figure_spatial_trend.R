make_figure_spatial_trend <-
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

        return(
            list(
                data_trend_pred = data_pred,
                data_hist = data_hist
            )
        )
    }
