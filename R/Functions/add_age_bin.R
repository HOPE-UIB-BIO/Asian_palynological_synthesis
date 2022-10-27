add_age_bin <-
    function(data_source,
             age_var_name = "age",
             bin_var_name = "BIN",
             bin_size = 1000,
             sel_method = c("backward", "forward")) {
        sel_method <- match.arg(sel_method)

        if (
            sel_method == "backward"
        ) {
            data_source %>%
                dplyr::mutate(
                    !!bin_var_name := floor(get(age_var_name) / bin_size) * bin_size
                ) %>%
                return()
        } else if (
            sel_method == "forward"
        ) {
            data_source %>%
                dplyr::mutate(
                    !!bin_var_name := ceiling(get(age_var_name) / bin_size) * bin_size
                ) %>%
                return()
        }
    }
