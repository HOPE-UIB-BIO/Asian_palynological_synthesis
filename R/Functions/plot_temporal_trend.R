plot_temporal_trend <-
    function(data_source,
             var_name,
             sel_scope = c("seq", "ecozone"),
             sel_type = c("var", "density")) {
        sel_scope <- match.arg(sel_scope)
        sel_type <- match.arg(sel_type)

        # add the name of the variable (to be plotted as facet)
        data_source <-
            data_source %>%
            dplyr::mutate(
                var_name = as.character(var_name)
            )

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
            ggplot2::scale_fill_manual(
                values = ecozone_pallete # [config]
            ) +
            ggplot2::scale_color_manual(
                values = ecozone_pallete # [config]
            ) +
            ggplot2::scale_x_continuous(
                breaks = seq(0, 12e3, 4e3),
                labels = seq(0, 12, 4)
            ) +
            ggplot2::theme(
                strip.text.x = ggplot2::element_text(
                    size = 14
                ),
                axis.text.x = ggplot2::element_text(
                    color = "black",
                    size = 16
                ),
                axis.text.y = ggplot2::element_text(
                    color = "black",
                    size = 16
                ),
                axis.title = ggplot2::element_text(
                    color = "black",
                    size = 18
                ),
                axis.title.x = ggplot2::element_text(
                    margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
                ),
                legend.position = "none"
            )

        if (
            sel_type == "density"
        ) {
            p0 <-
                p0 +
                ggplot2::scale_y_continuous(
                    breaks = c(0, 0.5, 1),
                    limits = c(0, 1)
                )
        }


        if (
            sel_scope == "seq"
        ) {

            # Data for boxplot (mean of each time BIN for each dataset)
            data_boxplot <-
                data_source %>%
                dplyr::group_by(dataset_id, BIN) %>%
                dplyr::summarise(
                    .groups = "drop",
                    var = mean(var)
                )

            res <-
                p0 +
                ggplot2::geom_line(
                    data = data_source,
                    ggplot2::aes(
                        x = age,
                        group = dataset_id
                    ),
                    color = "#2CA388",
                    alpha = 0.4
                ) +
                ggplot2::geom_violin(
                    data = data_boxplot,
                    ggplot2::aes(
                        x = BIN,
                        group = BIN
                    ),
                    col = NA,
                    fill = "#993300",
                    alpha = 0.5,
                ) +
                ggplot2::geom_boxplot(
                    data = data_boxplot,
                    ggplot2::aes(
                        x = BIN,
                        group = BIN
                    ),
                    col = "black",
                    fill = "white",
                    alpha = 1,
                    width = 250,
                    outlier.shape = NA
                )
        }

        if (
            sel_scope == "ecozone"
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
                            fill = Climate_zone,
                        ),
                        colour = NA,
                        alpha = 0.05
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

        return(res)
    }
