#----------------------------------------------------------#
#
#
#              Asian palynological synthesis
#
#             Temporal changes of chnage points
#
#
#                 K. Bhatta, O. Mottl
#                         2022
#
#----------------------------------------------------------#

# Plot the temporal trends of densities of change points of the estimates of
#  Hill's diversity, DCCA1, and MRT.

#--------------------------------------------------------#
# 1. Setup  ----
#--------------------------------------------------------#

library(here)

# Load configuration
source(
    here::here("R/00_Config_file.R")
)

#--------------------------------------------------------#
# 2. Load the data ----
#--------------------------------------------------------#

data_change_points <-
    readr::read_rds(
        here::here(
            "Data/Processed/Partitions/PAP_change_points_2022-09-14.rds"
        )
    )
data_all_estimates <-
    readr::read_rds(
        here::here(
            "Data/Processed/Data_for_analyses/Data_for_analyses-2022-09-15.rds"
        )
    )

#--------------------------------------------------------#
# 3. Restructure data ----
#--------------------------------------------------------#
# restructure data
data_cp_restructured <-
  data_all_estimates %>%
  dplyr::select(Climate_zone, dataset_id) %>%
  dplyr::inner_join(data_change_points, 
                    by = "dataset_id") %>%
  dplyr::mutate(
    mvrt_cp = purrr::map(mvrt_cp, ~ .x %>% 
                          unlist() %>% 
                          tibble::enframe(name = NULL,
                                          value = "age") %>%
                          dplyr::mutate(var_name = paste("MRT")
                                        )
                        ),
    roc_cp = purrr::map(roc_cp, ~ .x %>% 
                          unlist() %>% 
                          tibble::enframe(name = NULL,
                                          value = "age") %>%
                          dplyr::mutate(var_name = paste("RoC")
                          )
                        ),
    dcca_cp = purrr::map(dcca_cp, ~ .x %>% 
                          unlist() %>% 
                          tibble::enframe(name = NULL,
                                          value = "age") %>%
                          dplyr::mutate(var_name = paste("DCCA1")
                          )
                         ),
    dca_cp = purrr::map(dca_cp, ~ .x %>% 
                              unlist() %>% 
                              tibble::enframe(name = NULL,
                                              value = "age") %>%
                              dplyr::mutate(var_name = paste("DCA1")
                              )
                         )
    ) %>% 
  gather(key = "var",
         value = "change_points",
         -c(dataset_id, Climate_zone)
         ) %>% 
  unnest(change_points)

                  
# Select only required variables  
req_var <-
  c("DCCA1",
    "MRT",
    "n0",
    "n1",
    "n2",
    "n2_divided_by_n1",
    "n1_divided_by_n0",
    "RoC")

data_cp <- 
  data_cp_restructured %>% 
  dplyr::filter(var_name %in% req_var) %>% 
  dplyr::select(-var) %>% 
  # Set desired name and sequence of the facets
  dplyr::mutate(
    var_name = dplyr::case_when(
      var_name == "n0" ~ "N0",
      var_name == "n1" ~ "N1",
      var_name == "n2" ~ "N2",
      var_name == "n2_divided_by_n1" ~ "N2 divided by N1",
      var_name == "n1_divided_by_n0" ~ "N1 divided by N0",
      TRUE ~ var_name
    )
  ) %>% 
  dplyr::mutate(
    var_name = factor(var_name,
                      levels = c(
                        "DCCA1",
                        "MRT",
                        "N0",
                        "N1",
                        "N2",
                        "N2 divided by N1",
                        "N1 divided by N0",
                        "RoC"
                      )
    )
  ) %>% 
  add_age_bin(
    bin_size = bin_size # [config]
    ) %>%
    dplyr::filter(BIN >= 0) %>%
    tidyr::drop_na(var_name) 

#--------------------------------------------------------#
# 4. Estimate change point density for whole continent ----
#--------------------------------------------------------#
data_density_continet <-
  data_cp %>%
  dplyr::mutate(var = paste(var_name)) %>% 
  group_by(var) %>% 
  nest() %>% 
  dplyr::mutate(
    density_data = 
      purrr::map(data, 
                 .f = ~ REcopol::get_density(
                   data_source = .x$age,
                   reflected = TRUE,
                   values_range = c(
                     min(new_data_general$age), # new_data_general [Config]
                     max(new_data_general$age)
                   ), 
                   bw = 1000 / max(new_data_general$age),
                   n = max(new_data_general$age)
                 )
      )
  ) %>% 
  arrange(
    factor(var, 
           levels = c(
             "DCCA1",
             "MRT",
             "N0",
             "N1",
             "N2",
             "N2 divided by N1",
             "N1 divided by N0",
             "RoC"
           )
    )
  )

#--------------------------------------------------------#
# 5. Plot the data ----
#--------------------------------------------------------#
whole_continent <-
  data_density_continet %>%
  dplyr::mutate(
    plot =
      purrr::map2(
        .x = density_data,
        .y = data,
        .f = ~ .x %>%
          ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(x = var, y = density)) +
          ggplot2::coord_cartesian(xlim = range(.y$BIN)) +
          labs(x = "",
               y = "") +
          ggplot2::ggtitle(
            paste0(
              unique(.y$var_name)
            )
          ) + 
          theme(
            axis.text = element_text(color = "black", size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_text(color = "black", size = 11),
            plot.title = element_text(color = "black", size = 10),
            plot.margin = unit(c(0, 0.1, 0.26, 0.26), "cm") # t, r, b, l
          )
      )
  ) %>%
  purrr::pluck("plot") %>%
  ggpubr::ggarrange(plotlist = .,
                    ncol = 4,
                    nrow = 2) +
  theme(
    plot.margin = unit(c(0, 0, -0.5, -0.3), "cm"
    )
  )

library(grid)
plot_continent <-
  annotate_figure(
    whole_continent,
    left = textGrob(
      "Density of regression tree change points",
      rot = 90,
      vjust = 1,
      gp = gpar(cex = 1.2)
    ),
    bottom = textGrob("Time (yr BP)")
  )

ggsave(
  filename =
    here::here(
      "Outputs/Figures/Change_point_density_continent_170922.tiff"
    ),
  plot = plot_continent,
  dpi = 200,
  width = 18,
  height = 10,
  units = "cm",
  compress = "lzw"
)

#--------------------------------------------------------#
# 6. Estimate change point density for climate zones ----
#--------------------------------------------------------#
data_density_ecozone <-
  data_cp %>%
  dplyr::mutate(variable = paste(var_name)) %>% 
  group_by(Climate_zone) %>% 
  nest() %>% 
  ungroup() %>% 
  dplyr::mutate(
    per_zone_dat =
      purrr::map(
        data,
        ~ .x %>%
          group_by(variable) %>%
          nest() %>% 
          ungroup()
      ),
    density_data_nested =
      purrr::map(per_zone_dat, ~ .x %>%
                   dplyr::mutate(
                     density_data = 
                       purrr::map(data, 
                                  .f = ~ REcopol::get_density(
                                    data_source = .x$age,
                                    reflected = TRUE,
                                    values_range = c(
                                      min(new_data_general$age), # new_data_general [Config]
                                      max(new_data_general$age)
                                    ), 
                                    bw = 1000 / max(new_data_general$age),
                                    n = max(new_data_general$age)
                                  )
                              )
                           )
                        )
                    ) 

#--------------------------------------------------------#
# Plot the data ----
#--------------------------------------------------------#
full_dat_density <- 
  data_density_ecozone %>% 
  dplyr::select(Climate_zone, density_data_nested) %>% 
  unnest(density_data_nested) %>% 
  dplyr::select(-data) %>% 
  unnest(density_data) 

full_dat_density$variable <- 
  factor(full_dat_density$variable, 
         levels = c(
           "DCCA1",
           "MRT",
           "N0",
           "N1",
           "N2",
           "N2 divided by N1",
           "N1 divided by N0",
           "RoC")
  )

# Function to wrap the facet labels ('label_wrap_gen()' does not work for 
#  multiple facets together for an unknown reason

my_label <- 
  function(x) {
    # Wrap var names
    names(x)[[1]] <- stringr::str_wrap(gsub("_", " ", names(x)[[1]]), 11)
    # Wrap value labels
    x[[1]] <- stringr::str_wrap(gsub("_", " ", x[[1]]), 12)
    # Call label both with sep "\n"
    label_value(x)
  }

# library "grid" should be installed for facet_grid()
# library "lemon" should be installed for modifying axes of each panel
per_climate_zone <-
  full_dat_density %>%
  ggplot2::ggplot(
    aes(group = Climate_zone)
  ) +
  
  ggplot2::geom_line(
    aes(x = var,
        y = density,
        colour = Climate_zone)) +
  
  ggplot2::scale_x_continuous(
    limits = c(0, 12500),
    breaks = seq(0, 12500, 3000)
  ) +
  
  facet_rep_grid(variable ~ Climate_zone,
                 scale = "free_y",
                 labeller = my_label) +
  
  coord_capped_cart(bottom = 'both', left = 'both') +
  ggplot2::scale_colour_manual(values = ecozone_pallete) + # [config]
  
  labs(
    x = "Time (yr BP)", 
    y = "Density of regression tree change points"
  )  +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 15, angle = 0, hjust = 0),
    strip.text.y = element_text(size = 15, angle = -90, hjust = 0),
    axis.text = element_text(color = "black", size = 17),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(color = "black", size = 24),
    panel.border = element_blank(), 
    axis.line = element_line(),
    panel.spacing = unit(-2.3, "lines")
  )


ggsave(
  filename =
    here::here(
      "Outputs/Figures/Change_point_density_climate_zone_test_170922.tiff"
    ),
  plot = per_climate_zone,
  dpi = 200,
  width = 20,
  height = 30,
  units = "cm",
  compress = "lzw"
)
  
  
  
