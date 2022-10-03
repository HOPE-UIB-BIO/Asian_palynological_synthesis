plot_temporal_trend <-
    function(data_source,
             var_name,
             sel_grain = c("sequence", "climate-zone", "continent"),
             sel_type = c("var", "density"),
             y_limits,
             plot_rmse = TRUE,
             plot_summary = TRUE,
             def_color,
             def_of_color,
             def_violin_color,
             def_text_size = 16) {
        sel_grain <- match.arg(sel_grain)
        sel_type <- match.arg(sel_type)

        if (
            is.null(data_source)
        ) {
            return(NULL)
        }

        # add the name of the variable (to be plotted as facet)
        data_source <-
            data_source %>%
            dplyr::mutate(
                var_name = as.character(var_name)
            )

        if (
            missing(y_limits) ||
                is.null(y_limits) == TRUE ||
                is.na(y_limits) == TRUE
        ) {
            y_limits <-
                vector(mode = "numeric", length = 2)


            if (
                sel_type == "density"
            ) {
                y_limits[1] <- 0
                y_limits[2] <- 1
            } else {
                y_limits_raw <-
                    data_source %>%
                    dplyr::select(
                        dplyr::any_of(
                            c(
                                "var",
                                "upr",
                                "lwr"
                            )
                        )
                    ) %>%
                    unlist()

                y_limits[1] <-
                    min(y_limits_raw) * 0.9

                y_limits[2] <-
                    max(y_limits_raw) * 1.1
            }
        }

        # create a common visualisatin style for all figures
        p0 <-
            data_source %>%
            ggplot2::ggplot(
                ggplot2::aes(
                    x = age,
                    y = var
                )
            ) +
            ggplot2::theme_classic() +
            ggplot2::labs(
                x = "Time (yr BP)",
                y = "Estimate"
            ) +
            ggplot2::scale_x_continuous(
                breaks = seq(0, 12e3, 4e3),
                labels = seq(0, 12, 4)
            ) +
            ggplot2::scale_y_continuous(
                limits = y_limits,
                n.breaks = 3
            ) +
            ggplot2::theme(
                strip.text.x = ggplot2::element_text(
                    size = def_text_size
                ),
                axis.text.x = ggplot2::element_text(
                    color = "black",
                    size = def_text_size
                ),
                axis.text.y = ggplot2::element_text(
                    color = "black",
                    size = def_text_size
                ),
                axis.title = ggplot2::element_text(
                    color = "black",
                    size = def_text_size
                ),
                axis.title.x = ggplot2::element_text(
                    margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
                ),
                legend.position = "none"
            )


        if (
            sel_grain == "sequence"
        ) {
            data_seq <-
                data_source %>%
                dplyr::group_by(age) %>%
                dplyr::mutate(
                    var_diff = var - mean(var, na.rm = TRUE)
                ) %>%
                dplyr::group_by(dataset_id) %>%
                dplyr::mutate(
                    col_var = mean(var_diff^2) %>%
                        sqrt()
                ) %>%
                dplyr::ungroup()

            if (
                plot_rmse == TRUE
            ) {
                res <-
                    p0 +
                    ggplot2::geom_line(
                        data = data_seq,
                        ggplot2::aes(
                            x = age,
                            group = dataset_id,
                            color = col_var
                        ),
                        alpha = 0.2
                    ) +
                    ggplot2::scale_color_gradient(
                        high = def_of_color,
                        low = def_color
                    )
            } else {
                res <-
                    p0 +
                    ggplot2::geom_line(
                        data = data_seq,
                        ggplot2::aes(
                            x = age,
                            group = dataset_id
                        ),
                        color = def_color,
                        alpha = 0.2
                    )
            }

            if (
                plot_summary == TRUE
            ) {

                # Data for boxplot (mean of each time BIN for each dataset)
                data_boxplot <-
                    data_source %>%
                    add_age_bin() %>%
                    dplyr::group_by(dataset_id, BIN) %>%
                    dplyr::summarise(
                        .groups = "drop",
                        var = mean(var)
                    )

                res <-
                    res +
                    ggplot2::geom_violin(
                        data = data_boxplot,
                        ggplot2::aes(
                            x = BIN,
                            group = BIN
                        ),
                        col = NA,
                        fill = def_violin_color,
                        alpha = 0.25,
                    ) +
                    ggplot2::geom_boxplot(
                        data = data_boxplot,
                        ggplot2::aes(
                            x = BIN,
                            group = BIN
                        ),
                        col = "black",
                        fill = "white",
                        alpha = 0.5,
                        width = 250,
                        size = 0.1,
                        outlier.shape = NA
                    )
            }
        }

        if (
            sel_grain == "climate-zone"
        ) {
            p0 <-
                p0 +
                ggplot2::scale_fill_manual(
                    values = ecozone_pallete # [config]
                ) +
                ggplot2::scale_color_manual(
                    values = ecozone_pallete # [config]
                )

            if (
                sel_type == "var"
            ) {
                p1 <-
                    p0 +
                    ggplot2::geom_ribbon(
                        data = data_source,
                        ggplot2::aes(
                            ymin = lwr,
                            ymax = upr,
                            fill = Climate_zone,
                        ),
                        colour = NA,
                        alpha = 0.1
                    )
            } else {
                p1 <-
                    p0
            }

            res <-
                p1 +
                ggplot2::geom_line(
                    data = data_source,
                    ggplot2::aes(
                        group = Climate_zone,
                        colour = Climate_zone
                    ),
                    size = 1
                )
        }

        if (
            sel_grain == "continent"
        ) {
            if (
                sel_type == "var"
            ) {
                p1 <-
                    p0 +
                    ggplot2::geom_ribbon(
                        data = data_source,
                        ggplot2::aes(
                            ymin = lwr,
                            ymax = upr,
                        ),
                        fill = "gray50",
                        colour = NA,
                        alpha = 1
                    )
            } else {
                p1 <-
                    p0
            }

            res <-
                p1 +
                ggplot2::geom_line(
                    data = data_source,
                    size = 1
                )
        }

        return(res)
    }
