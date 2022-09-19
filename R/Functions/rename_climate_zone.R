rename_climate_zone <-
    function(data_source) {
        data_source %>%
            # Modify the climate zones as suggested by John
            dplyr::mutate(
                Climate_zone = dplyr::case_when(
                    ecozone_koppen_15 == "Arid_Desert" ~ "Arid",
                    ecozone_koppen_15 == "Arid_Steppe" ~ "Arid",
                    ecozone_koppen_15 == "Cold_Dry_Summer" ~ "Cold - Dry",
                    ecozone_koppen_15 == "Cold_Dry_Winter" ~ "Cold - Dry",
                    ecozone_koppen_15 == "Cold_Without_dry_season" ~ "Cold - Without dry season",
                    ecozone_koppen_15 == "Polar_Frost" ~ "Polar - Frost",
                    ecozone_koppen_15 == "Polar_Tundra" ~ "Polar",
                    ecozone_koppen_15 == "Temperate_Dry_Summer" ~ "Temperate",
                    ecozone_koppen_15 == "Temperate_Dry_Winter" ~ "Temperate",
                    ecozone_koppen_15 == "Temperate_Without_dry_season" ~ "Temperate",
                    ecozone_koppen_15 == "Tropical_Monsoon" ~ "Tropical - Monsoon",
                    ecozone_koppen_15 == "Tropical_Rainforest" ~ "Tropical - Rainforest",
                    ecozone_koppen_15 == "Tropical_Savannah" ~ "Tropical - Savannah",
                    TRUE ~ ecozone_koppen_15
                )
            ) %>%
            return()
    }
