#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#                     Project setup
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#
# Configuration script with the variables that should be consistent throughout
#   the whole repo. It loads packages, defines important variables,
#   authorises the user, and saves config file.


#----------------------------------------------------------#
# 1. Load packages -----
#----------------------------------------------------------#

if (!exists("update_repo_packages")) {
  update_repo_packages <- TRUE
}

if (update_repo_packages == TRUE) {

  # install RRatepol from github
  if (!exists("already_installed_rratepol")) {
    already_installed_rratepol <- FALSE
  }

  if (already_installed_rratepol == FALSE) {
    devtools::install_github("HOPE-UIB-BIO/R-Ratepol-package",
      quiet = FALSE,
      upgrade = FALSE
    )
    already_installed_rratepol <- TRUE
  }

  # install REcopol from GitHub
  if (!exists("already_installed_recopol")) {
    already_installed_recopol <- FALSE
  }

  if (already_installed_recopol == FALSE) {
    devtools::install_github("HOPE-UIB-BIO/R-Ecopol-package",
      quiet = FALSE,
      upgrade = FALSE
    )
    already_installed_recopol <- TRUE
  }

  if (!exists("already_synch")) {
    already_synch <- FALSE
  }

  if (already_synch == FALSE) {
    library(here)
    # synchronise the package versions
    renv::restore(lockfile = here::here("renv/library_list.lock"))
    already_synch <- TRUE

    # save snapshot of package versions
    # renv::snapshot(lockfile =  "renv/library_list.lock")  # do only for update
  }
}

# define packages
package_list <-
  c(
    "assertthat",
    "devtools",
    "dplyr",
    "forcats",
    "ggplot2",
    "ggpubr",
    "ggside",
    "here",
    "mgcv",
    "purrr",
    "readr",
    "raster",
    "REcopol",
    "RRatepol",
    "renv",
    "rgdal",
    "roxygen2",
    "stringr",
    "tibble",
    "tidyr",
    "usethis"
  )

# load all packages
sapply(package_list, library, character.only = TRUE)


#----------------------------------------------------------#
# 2. Define space -----
#----------------------------------------------------------#

current_date <- Sys.Date()

# project directory is set up by 'here' package, Adjust if needed
current_dir <- here::here()


#----------------------------------------------------------#
# 3. Load functions -----
#----------------------------------------------------------#

# get vector of general functions
fun_list <-
  list.files(
    path = "R/Functions/",
    pattern = ".R",
    recursive = TRUE
  )

# source them
if (
  length(fun_list > 0)
) {
  sapply(
    here::here(
      paste0("R/Functions/", fun_list, sep = "")
    ),
    source
  )
}

#----------------------------------------------------------#
# 4. Authorise the user -----
#----------------------------------------------------------#

# if applicable

#----------------------------------------------------------#
# 5. Define variables -----
#----------------------------------------------------------#

# 5.1. Analystical criteria -----

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

# 5.2. PAP estimation specification -----

n_rand <- 1000

transformation_coef <- "chisq"

smoothing_method <- "age.w"
min_points_smoothing <- 5
max_points_smoothing <- 9
age_range_smoothing <- 500
working_units_selection <- "MW"
size_of_bin <- 500
n_mowing_windows <- 5
which_level_select_in_bin <- "random"
standardise_data <- TRUE
n_individuals_to_standardise <- 150
tranform_to_proportions <- TRUE

# 5.3. GAM Modelling 
# tables with names of variables, errors, and dataframes
vars_table <-
    tibble::tibble(
        var_name = c(
            "n0",
            "n1",
            "n2",
            "dcca_axis_1",
            "dca_axis_1",
            "n2_divided_by_n1",
            "n1_divided_by_n0",
            "ROC",
            "Peak"
        ),
        sel_error = c(
            rep("mgcv::Tweedie(p = 1.1)", 5),
            rep("mgcv::betar(link = 'logit')", 2),
            "mgcv::Tweedie(p = 1.1)",
            "stats::quasibinomial(link = 'logit')"
        ),
        sel_data = c(
            rep("data_diversity", 7),
            rep("data_roc", 2)
        )
    )

#----------------------------------------------------------#
# 6. Graphical options -----
#----------------------------------------------------------#

## examples
# set ggplot output
ggplot2::theme_set(
  ggplot2::theme_classic()
)

# define general
text_size <- 10
line_size <- 0.1

# define output sizes
image_width <- 16
image_height <- 12
image_units <- "cm"

climate_zone_vec <-
  c(
    "Arid",
    "Cold - Dry",
    "Cold - Without dry season",
    "Polar",
    "Polar - Frost",
    "Temperate",
    "Tropical - Monsoon",
    "Tropical - Rainforest",
    "Tropical - Savannah"
  )


# define pallets
ecozone_pallete_full <-
  c(
    "#FFCC99",
    "#993300",
    "#FF6600",
    "#3399FF",
    "#999999",
    "#00CCCC",
    "#99CC00",
    "#006600",
    "#996600"
  ) %>%
  rlang::set_names(
    nm = climate_zone_vec
  )

ecozone_pallete <-
  ecozone_pallete_full[
    c(
      "Arid",
      "Cold - Dry",
      "Cold - Without dry season",
      "Polar",
      "Temperate"
    )
  ]

# define common color
