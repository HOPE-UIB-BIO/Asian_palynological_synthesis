run_dcca <-
  function(data_pollen, data_levels, dataset_id) {
    
    cat(paste0(dataset_id), "\n")
    
    sp <-
      data_pollen %>%
      as.data.frame() %>%
      column_to_rownames("sample_id") %>% 
      round(., digits = 2)
    
    env <- data_levels$age
    
    # This should be the location where 'canoco.exe' file is located   
    output <-
      select.DCCA.can(
        here::here("R/Functions/DCCA"),
        sp,
        env)
    
    return(output)
  }
