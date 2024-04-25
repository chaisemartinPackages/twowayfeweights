#' Internal workhorse function for creating the return object of a
#' `twowayfeweights()` call.
#' 
#' @param dat A data frame, as per the return object from
#'   `twowayfeweights_calculate()`.
#' @param beta Coefficient value of the treatment variable ("D"), again as per
#'   the return object of `twowayfeweights_calculate()`.
#' @param controls A vector indicating the column names of random weights.
#' @param treatments A vector indicating the column names of other treatments.
#' @returns A list.
#' @details This function is normally run directly after
#'   `twowayfeweights_calculate()`.
#' @importFrom magrittr %>%
#' @noRd
twowayfeweights_result = function(dat, beta, random_weights, treatments = NULL) {
  suppressWarnings({ 
  .data = NULL
  
  # Two distinct cases/workflows:
  #  1) No other treatments,
  #  2) With other treatments
  
  if (is.null(treatments)) {
    
    # Avoid overcounting of positive and negative weights close to 0
    limit_sensitivity = 10^(-10)
    dat$weight_result = ifelse(dat$weight_result < limit_sensitivity & dat$weight_result > -limit_sensitivity, 0, dat$weight_result)
    ret = twowayfeweights_summarize_weights(dat, "weight_result")
    
    W_mean = weighted.mean(dat$W, dat$nat_weight)
    # Modif. Diego: DoF adjustment to the sd of w_gt
    M = sum(dat$nat_weight != 0)
    W_sd = sqrt(sum(dat$nat_weight * (dat$W - W_mean)^2, na.rm = TRUE)) * sqrt(M/(M - 1))
    sensibility = abs(beta) / W_sd
    
    dat_result = dat %>% 
      dplyr::select(.data$T, .data$G, .data$weight_result) %>% 
      dplyr::rename(weight = .data$weight_result)
    
    ret$dat_result = dat_result
    ret$beta = beta
    ret$sensibility = sensibility
    if (length(random_weights) > 0) {
      ret$mat = twowayfeweights_test_random_weights(dat, random_weights)
    }
    
    if (ret$sum_minus < 0) {
      dat_sens = dat %>%
        dplyr::filter(.data$weight_result != 0) %>%
        dplyr::arrange(dplyr::desc(.data$W)) %>% 
        dplyr::mutate(P_k = 0, S_k = 0, T_k = 0)
      
      # Modif. Diego: Replaced the previous two loops with build-in routines
      N = nrow(dat_sens)
      dat_sens = dat_sens[order(dat_sens$W, -dat_sens$G, -dat_sens$T),]
      dat_sens$Wsq = dat_sens$nat_weight * dat_sens$W^2
      dat_sens$P_k = cumsum(dat_sens$nat_weight) 
      dat_sens$S_k = cumsum(dat_sens$weight_result) 
      dat_sens$T_k = cumsum(dat_sens$Wsq) 
      dat_sens = dat_sens[order(-dat_sens$W, dat_sens$G, dat_sens$T),]
      dat_sens = dat_sens %>% 
        dplyr::mutate(sens_measure2 = abs(beta) / sqrt(.data$T_k + .data$S_k^2 / (1 - .data$P_k))) %>%
        dplyr::mutate(indicator = as.numeric(.data$W < - .data$S_k / (1 - .data$P_k)))
      dat_sens$indicator[1] = 0
      dat_sens = dat_sens %>%
        dplyr::mutate(indicator_l = dplyr::lag(.data$indicator, default = -1))
      dat_sens = dat_sens %>%
        dplyr::rowwise() %>%
        dplyr::mutate(indicator=max(.data$indicator, .data$indicator_l))
      total_indicator = sum(dat_sens$indicator)
      sensibility2 = dat_sens$sens_measure2[N - total_indicator + 1]
      ret$sensibility2 = sensibility2

      # Since, with one treatment, we could have either D or D0 as the main treatment, 
      # the row below computes the number of cells such that their treatment is different than 0
      ret$tot_cells = sum(as.numeric(dat$nat_weight != 0), na.rm = TRUE)
    }
    
  } else {
    limit_sensitivity = 10^(-10)
    for (v in c("result", treatments)) {
      dat[[paste0("weight_",v)]] = ifelse(dat[[paste0("weight_",v)]] < limit_sensitivity & dat[[paste0("weight_",v)]] > -limit_sensitivity, 0, dat[[paste0("weight_",v)]])
    }
    
    columns = c("T", "G", "weight_result")
    ret = twowayfeweights_summarize_weights(dat, "weight_result")
    ret$tot_cells = sum(as.numeric(dat$nat_weight != 0), na.rm = TRUE)
    
    if (length(random_weights) > 0) {
      ret$mat = twowayfeweights_test_random_weights(dat, random_weights)
    }
    
    for (treatment in treatments) {
      varname = fn_treatment_weight_rename(treatment)
      columns = c(columns, varname)
      ret2 = twowayfeweights_summarize_weights(dat, varname)
      ret[[treatment]] <- ret2
      ret[[treatment]]$tot_cells = sum(as.numeric(dat[[treatment]] != 0), na.rm = TRUE)
    }
    dat_result = dat %>% 
      dplyr::select_at(dplyr::vars(columns)) %>% 
      dplyr::rename(weight = .data$weight_result)
    
    ret$beta = beta
    ret[["dat_result"]] <- dat_result
    
  }
  })
  
  
  return(ret)
  
}