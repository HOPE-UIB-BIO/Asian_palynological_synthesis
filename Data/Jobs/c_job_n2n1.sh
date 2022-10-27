#!/bin/tcsh
#SBATCH -J c_n2n1 # (job name)
#SBATCH -o c_n2n1.%j.out # (name of the job output file, %j expands to the job name)
#SBATCH -N 10 # (Number of requested nodes)
#SBATCH --ntasks-per-node 32 # (Number of requested cores per node)
#SBATCH -t 24:00:00 # (Requested wall time)

module load gnu10 R

R --no-save <<EOF
Sys.setenv(LANG = "en")

# define packages
package_list <-
    c(
        "dplyr",
        "here",
        "mgcv",
        "purrr",
        "readr",
        "REcopol",
        "rlang",
        "tibble",
        "tidyr"
    )

# load all packages
sapply(package_list, library, character.only = TRUE)

message("load function")

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

# Values for predicting trends

time_step <- 100

age_min <- 0

age_max <- 12e3

age_vec <-
    seq(age_min, age_max, by = time_step)

# data for predicting
new_data_general <-
    tibble::tibble(
        age = age_vec
    )

bin_size <- 1000

# tables with names of variables, errors, and dataframes
vars_table <-
    tibble::tibble(
        var_name = c(
            "n0",
            "n1",
            "n2",
            "dcca_axis_1",
            "n2_divided_by_n1",
            "n1_divided_by_n0",
            "ROC",
            "Peak"
        ),
        sel_error = c(
            rep("mgcv::Tweedie(p = 1.1)", 4),
            rep("mgcv::betar(link = 'logit')", 2),
            "mgcv::Tweedie(p = 1.1)",
            "stats::quasibinomial(link = 'logit')"
        ),
        sel_data = c(
            rep("data_diversity", 6),
            rep("data_roc", 2)
        )
    )


max_temporal_k <- 24

message("load data")

# Load the data
data_all_estimates <-
    readr::read_rds(
        here::here(
            "Data/Input/Data_for_analyses-2022-09-29.rds"
        )
    )

print(dim(data_all_estimates))

# diversity estimates
data_diversity <-
    data_all_estimates %>%
    dplyr::select(
        Climate_zone, dataset_id, PAP_merge
    ) %>%
    tidyr::unnest(PAP_merge) %>%
    add_age_bin(
        bin_size = bin_size # [config]
    ) %>%
    dplyr::filter(BIN >= 0) %>%
    dplyr::mutate(
        # adjust values so small that they are smaller than zero
        dcca_axis_1 = ifelse(dcca_axis_1 < 0, 0, dcca_axis_1)
    )

print(dim(data_diversity))

# RoC estimates
data_roc <-
    data_all_estimates %>%
    dplyr::select(
        Climate_zone, dataset_id, PAP_roc
    ) %>%
    tidyr::unnest(PAP_roc) %>%
    dplyr::rename(
        age = Age
    ) %>%
    add_age_bin(
        bin_size = bin_size # [config]
    ) %>%
    dplyr::filter(BIN >= 0)

print(dim(data_roc))

message("fitting model")

var_n <- 5

var_sel <- vars_table[var_n, "var_name", drop = TRUE]
error_sel <- vars_table[var_n, "sel_error", drop = TRUE]
data_sel <- vars_table[var_n, "sel_data", drop = TRUE]

# Fit GAM model
data_mod <-
    REcopol::fit_hgam(
        x_var = "age",
        y_var = var_sel,
        group_var = "dataset_id",
        smooth_basis = "tp",
        data_source = get(data_sel),
        error_family = error_sel,
        sel_k = max_temporal_k, # [config]
        common_trend = TRUE,
        use_parallel = TRUE,
        max_itiration = 200
    )

message("saving model")

readr::write_rds(
    x = data_mod,
    file = paste0(
        here::here(
            "Data/Output"
        ),
        "/",
        var_sel,
        ".rds"
    ),
    compress = "gz"
)

q()
EOF
