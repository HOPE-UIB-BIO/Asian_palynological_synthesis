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
    "here",
    "mgcv",
    "purrr",
    "readr",
    "REcopol",
    "RRatepol",
    "renv",
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

# define pallets

# define common color
