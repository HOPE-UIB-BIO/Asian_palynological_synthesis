plot_temporal_grid <-
    function(data_source,
             sel_vars,
             sel_var_type = c("var", "density"),
             sel_grain = c(
                 "sequence",
                 "climate-zone",
                 "continent"
             ),
             use_limits = TRUE,
             plot_rmse = TRUE,
             plot_summary = TRUE,
             def_color,
             def_of_colo,
             def_violin_color,
             def_text_size,
             heading_text_multiplier = 1.5) {

        # helper function to order grobs to grid
        draw_grid <-
            function(data_source,
                     n_cols,
                     n_rows) {
                draw_row <-
                    function(data_source,
                             fig_index) {
                        Reduce(
                            f = cbind,
                            x = data_source[fig_index]
                        ) %>%
                            return()
                    }

                dummy_matrix <-
                    1:(n_cols * n_rows) %>%
                    matrix(.,
                        nrow = n_cols,
                        ncol = n_rows
                    ) %>%
                    t()

                dummy_list <-
                    vector("list", length = n_rows)

                for (i in 1:n_rows) {
                    dummy_list[[i]] <-
                        draw_row(
                            data_source,
                            fig_index = as.numeric(dummy_matrix[i, ])
                        )
                }

                Reduce(
                    f = rbind,
                    x = dummy_list
                ) %>%
                    return()
            }

        # helper function to plot vector as gglot
        plot_name <-
            function(var_name,
                     sel_size = 16) {
                tibble::tibble(
                    x = 1,
                    y = 1,
                    label = var_name
                ) %>%
                    ggplot2::ggplot(
                        ggplot2::aes(
                            x, y,
                            label = label
                        ),
                        size = sel_size
                    ) +
                    ggplot2::geom_label() +
                    ggplot2::theme_void() %>%
                    return()
            }

        if (
            use_limits == FALSE
        ) {
            data_source <-
                data_source %>%
                dplyr::mutate(
                    y_limits = NA
                )
        }

        rows_names <-
            sel_vars

        n_rows <-
            length(rows_names)

        cols_names <-
            sel_grain

        n_cols <-
            length(cols_names) * length(sel_var_type)

        data_sub <-
            data_source %>%
            dplyr::filter(
                var_name %in% sel_vars
            ) %>%
            dplyr::filter(
                var_type %in% sel_var_type
            ) %>%
            dplyr::filter(
                grain %in% sel_grain
            )

        data_with_indiv_plots <-
            data_sub %>%
            dplyr::mutate(
                ifigure = purrr::pmap(
                    .l = list(
                        var_name, # ..1
                        merge_data, # ..2
                        grain, # ..3
                        var_type, # ..4
                        y_limits # ..5
                    ),
                    .f = ~ plot_temporal_trend(
                        data_source = ..2,
                        var_name = ..1,
                        sel_grain = ..3,
                        sel_type = ..4,
                        y_limits = ..5,
                        plot_rmse = plot_rmse,
                        plot_summary = plot_summary,
                        def_color = def_color,
                        def_of_colo = def_of_colo,
                        def_violin_color = def_violin_color,
                        def_text_size = def_text_size
                    )
                )
            )

        # create all possible combination
        dummy_table <-
            data_with_indiv_plots %>%
            dplyr::select(
                var_name, grain, var_type
            ) %>%
            dplyr::mutate(
                var_name = factor(
                    var_name,
                    levels = sel_vars
                ),
                grain = factor(
                    grain,
                    levels = sel_grain
                ),
                var_type = factor(
                    var_type,
                    levels = sel_var_type
                )
            ) %>%
            tidyr::expand(
                var_name, var_type, grain
            ) %>%
            dplyr::arrange(var_name, var_type, grain)

        # create an empty plot
        dummy_plot <-
            ggplot2::ggplot() +
            ggplot2::geom_blank() +
            ggplot2::theme_void()

        dummy_plot_grob <-
            ggplot2::ggplotGrob(dummy_plot)

        # use dummy table to have all posible figures,
        #   replace the missing figures with empty
        #   transform all plos to Grobs
        suppressWarnings(
            data_with_indiv_plots_full <-
                dummy_table %>%
                dplyr::left_join(
                    data_with_indiv_plots,
                    by = c("var_name", "var_type", "grain")
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
                                no_axis <-
                                    .x +
                                    ggpubr::rremove("xylab")

                                ggplot2::ggplotGrob(no_axis) %>%
                                    return()
                            } else {
                                return(dummy_plot_grob)
                            }
                        }
                    )
                )
        )

        data_fig_all <-
            data_with_indiv_plots_full %>%
            purrr::pluck("ifigure_grob") %>%
            draw_grid(
                data_source = .,
                n_cols = n_cols,
                n_rows = n_rows
            )

        cols_names_grob_list <-
            cols_names %>%
            purrr::map(
                .f = ~ {
                    sel_plot <-
                        plot_name(
                            var_name = .x,
                            sel_size = def_text_size * heading_text_multiplier
                        )
                    plot_grob <-
                        ggplot2::ggplotGrob(sel_plot) %>%
                        return()
                }
            )

        cols_names_fig_list_full <-
            c(
                list(dummy_plot_grob),
                rep(cols_names_grob_list, length(sel_var_type))
            ) %>%
            Reduce(cbind, .)

        rows_names_grob_list <-
            rows_names %>%
            purrr::map(
                .f = ~ {
                    sel_plot <-
                        plot_name(
                            var_name = .x,
                            sel_size = def_text_size * heading_text_multiplier
                        )
                    plot_grob <-
                        ggplot2::ggplotGrob(sel_plot) %>%
                        return()
                }
            ) %>%
            Reduce(rbind, .)

        fin <-
            rbind(
                cols_names_fig_list_full,
                cbind(
                    rows_names_grob_list,
                    data_fig_all
                )
            )
        
        return(fin)
    }
