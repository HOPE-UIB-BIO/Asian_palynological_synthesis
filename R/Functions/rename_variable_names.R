rename_variable_names <-
    function(data_source) {
        data_source %>%
            dplyr::mutate(
                var_name = dplyr::case_when(
                    var_name == "n0" ~ "N0",
                    var_name == "n1" ~ "N1",
                    var_name == "n2" ~ "N2",
                    var_name == "n2_divided_by_n1" ~ "N2 divided by N1",
                    var_name == "n1_divided_by_n0" ~ "N1 divided by N0",
                    var_name == "ROC" ~ "RoC",
                    var_name == "roc" ~ "RoC",
                    var_name == "peakpoints" ~ "Peak-points",
                    var_name == "Peak" ~ "Peak-points",
                    var_name == "dcca_axis_1" ~ "DCCA1",
                    var_name == "dcca" ~ "DCCA1",
                    var_name == "dca_axis_1" ~ "DCA1",
                    var_name == "dca" ~ "DCA1",
                    var_name == "mvrt" ~ "MVRT",
                    TRUE ~ var_name
                )
            ) %>%
            dplyr::mutate(
                var_name = factor(var_name,
                    levels = c(
                        "DCCA1",
                        "DCA1",
                        "MVRT",
                        "N0",
                        "N1",
                        "N2",
                        "N2 divided by N1",
                        "N1 divided by N0",
                        "RoC",
                        "Peak-points"
                    )
                )
            ) %>%
            return()
    }
