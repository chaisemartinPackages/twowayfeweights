#' Internal workhorse function for calculating the twoway FE weights
#' 
#' @param dat A data frame.
#' @param type Desired type of calculation.
#' @param controls A vector of controls.
#' @param treatments Additional treatments.
#' @returns A list.
#' @importFrom fixest feols
#' @importFrom magrittr %>%
#' @importFrom stats as.formula resid weighted.mean weights
#' @noRd
twowayfeweights_calculate = function(
    dat, 
    type = c("feTR", "fdTR", "feS", "fdS"),
    controls = NULL, 
    treatments = NULL
) {
  
  suppressWarnings({
  .data = NULL
  
  type = match.arg(type)
  
  # GM: Already have this check higher up, but useful for debugging
  if (!is.null(treatments) && type != "feTR") {
    stop("When the `other_treatments` argument is specified, you need to specify `type = 'feTR'` too.")
  }
  
  type_TR = type %in% c("feTR", "fdTR")
  type_fe = type %in% c("feTR", "feS")
  
  if (type_TR) {
    DVAR = if (type=="feTR") "D" else if (type=="fdTR") "D0"
    mean_D = weighted.mean(dat[[DVAR]], dat$weights, na.rm = TRUE)
  }

  obs = sum(dat$weights)
  gdat = dat %>%
    dplyr::group_by(.data$G, .data$T) %>%
    dplyr::summarise(P_gt = sum(.data$weights)) %>% dplyr::ungroup()
  dat = dat %>% 
    dplyr::left_join(gdat, by=c("T", "G")) %>% 
    dplyr::mutate(P_gt = .data$P_gt / obs)
  if (type_TR) {
    dat = dplyr::mutate(dat, nat_weight = .data$P_gt * .data[[DVAR]] / mean_D)
  }

  # GM: Conveniently set formula using fixest macros (see: ?dsb)  
  if (is.null(controls)) controls = 1
  fes = "Tfactor"
  if (type_fe) fes = c("G", fes)
  
  # Add non-NULL treatment vars
  xvars = c(controls, treatments)
  
  # GM: could we make this work through if(!type_TR) ?
  if (type=="fdS") {
    denom.lm = feols(D ~ .[xvars] | .[fes], data = subset(dat, weights!=0), weights = dat$weights)
  } else {
    denom.lm = feols(D ~ .[xvars] | .[fes], data = dat, weights = dat$weights)
  }
  
  EPS_VAR = if (type_fe) "eps_1" else "eps_2"
  
  # GM: could we make this if(!type_TR), combined with weights !=0 above?
  if (type_fe || type=="fdS") {
    dat[[EPS_VAR]] = resid(denom.lm)
  } else if (type == "fdTR") {
    dat[[EPS_VAR]] = resid(denom.lm, na.rm = FALSE)
    dat[[EPS_VAR]] = ifelse(is.na(dat[[EPS_VAR]]), 0, dat[[EPS_VAR]])
  }
  
  
  # Beta reg ----
  
  ## Adjust data first
  if (type=="feTR") {

    dat[["eps_1_E_D_gt"]] = dat[[EPS_VAR]] * dat[[DVAR]]
    if (is.null(treatments)) {
      denom_W = weighted.mean(dat[["eps_1_E_D_gt"]], dat[["weights"]], na.rm = TRUE)
    } else {
      
      denom_W = mean(dat[["eps_1_E_D_gt"]], na.rm = TRUE)
    }
    
    dat = dat %>% 
      dplyr::mutate(W = .data[[EPS_VAR]] * mean_D / denom_W) %>% 
      dplyr::mutate(weight_result = .data[["W"]] * .data[["nat_weight"]])
    
    if (!is.null(treatments)) {
      for (treatment in treatments) {
        varname = fn_treatment_weight_rename(treatment)
        dat[[varname]] = dat[["W"]] * dat[["P_gt"]] * dat[[treatment]] / mean_D
      }
    }
    
    dat = dat %>%
      dplyr::select(-.data[[EPS_VAR]], -.data[["P_gt"]])
    
  } else if (type=="feS") {
    
    dat = dat %>% 
      dplyr::mutate(eps_1_weight = .data[[EPS_VAR]] * .data$weights) %>%
      dplyr::arrange(.data$G, .data$Tfactor) %>%
      dplyr::group_by(.data$G) %>%
      dplyr::mutate(E_eps_1_g_ge_aux = rev(cumsum(rev(.data$eps_1_weight)))) %>%
      dplyr::mutate(weights_aux = rev(cumsum(rev(.data$weights)))) %>%
      dplyr::mutate(E_eps_1_g_ge = .data$E_eps_1_g_ge_aux / .data$weights_aux) %>% dplyr::ungroup()
    
  } else if (type=="fdTR") {
    
    dat = dat %>% 
      dplyr::mutate(eps_2 = ifelse(is.na(.data[[EPS_VAR]]), 0, .data[[EPS_VAR]])) # dup with l. 55?
    
  }
  
  ## New regression
  xvars = c("D", xvars)
  if (type=="fdS") {
    beta.lm = feols(Y ~ .[xvars] | .[fes], data = subset(dat, weights != 0), weights = dat$weights, only.coef = TRUE)
  } else {
    beta.lm = feols(Y ~ .[xvars] | .[fes], data = dat, weights = dat$weights, only.coef = TRUE)
  }
  beta = beta.lm[["D"]]
  
  if (type == "feTR") {
    
    # * Keeping only one observation in each group * period cell
    # This should be done after this function
    # bys `group' `time': gen group_period_unit=(_n==1)	
    # 	drop if group_period_unit==0
    # 	drop group_period_unit
    dat = dat %>%
      dplyr::group_by(.data$G, .data$Tfactor) %>%
      dplyr::filter(dplyr::row_number(.data$D) == 1) %>% dplyr::ungroup()
    
  } else if (type == "fdTR") {
    
    dat = dat %>% 
      dplyr::arrange(.data$G, .data$TFactorNum) %>%
      dplyr::group_by(.data$G) %>% 
      dplyr::mutate(w_tilde_2 = ifelse(.data$TFactorNum + 1 == dplyr::lead(.data$TFactorNum), .data$eps_2 - dplyr::lead(.data$eps_2) * (dplyr::lead(.data$P_gt) / .data$P_gt), NA)) %>%
      dplyr::mutate(w_tilde_2 = ifelse(is.na(.data$w_tilde_2) | is.infinite(.data$w_tilde_2), .data$eps_2, .data$w_tilde_2)) %>%
      dplyr::mutate(w_tilde_2_E_D_gt = .data$w_tilde_2 * .data$D0) %>% dplyr::ungroup()
    
    denom_W = weighted.mean(dat$w_tilde_2_E_D_gt, dat$P_gt, na.rm = TRUE)
    dat = dat %>% 
      # dplyr::mutate(W = .data$w_tilde_2 * mean_D0 / denom_W) %>% 
      dplyr::mutate(W = .data$w_tilde_2 * mean_D / denom_W) %>% 
      dplyr::mutate(weight_result = .data$W * .data$nat_weight)
    dat = dat %>%
      dplyr::select(-.data$eps_2, -.data$P_gt, -.data$w_tilde_2, -.data$w_tilde_2_E_D_gt)
    
  } else if (type=="feS") {
    
    dat = dat %>% 
      dplyr::arrange(.data$G, .data$Tfactor) %>%
      dplyr::group_by(.data$G) %>% 
      dplyr::mutate(delta_D = ifelse(.data$TFactorNum - 1 == dplyr::lag(.data$TFactorNum), .data$D - dplyr::lag(.data$D), NA)) %>%
      dplyr::filter(!is.na(.data$delta_D)) %>%
      dplyr::mutate(abs_delta_D = abs(.data$delta_D)) %>%
      dplyr::mutate(s_gt = dplyr::case_when(.data$delta_D > 0 ~ 1,
                                            .data$delta_D < 0 ~ -1,
                                            TRUE ~ 0)) %>%
      dplyr::mutate(nat_weight = .data$P_gt * .data$abs_delta_D) %>% dplyr::ungroup()
    
    P_S = sum(dat$nat_weight, na.rm = TRUE)
    dat = dat %>% 
      dplyr::mutate(nat_weight = .data$nat_weight / P_S) %>%
      dplyr::mutate(om_tilde_1 = .data$s_gt * .data$E_eps_1_g_ge / .data$P_gt)
    
    denom_W = weighted.mean(dat$om_tilde_1, dat$nat_weight, na.rm = TRUE)
    dat = dat %>%
      dplyr::mutate(W = .data$om_tilde_1 / denom_W) %>%
      dplyr::mutate(weight_result = .data$W * .data$nat_weight) %>%
      dplyr::select(-.data$eps_1, -.data$P_gt, -.data$om_tilde_1, -.data$E_eps_1_g_ge,
                    -.data$E_eps_1_g_ge_aux, -.data$weights_aux, -.data$abs_delta_D, -.data$delta_D)
    
    
  } else if (type=="fdS") {
    dat = dat %>%
      dplyr::mutate(s_gt = dplyr::case_when(.data$D > 0 ~ 1,
                                            .data$D < 0 ~ -1,
                                            TRUE ~ 0)) %>%
      dplyr::mutate(abs_delta_D = abs(.data$D)) %>%
      dplyr::mutate(nat_weight = .data$P_gt * .data$abs_delta_D)
    
    P_S = sum(dat$nat_weight)
    dat = dat %>% 
      dplyr::mutate(nat_weight = .data$nat_weight / P_S) %>%
      dplyr::mutate(W = .data$s_gt * .data$eps_2)
    denom_W = weighted.mean(dat$W, dat$nat_weight, na.rm = TRUE)
    dat = dat %>% 
      dplyr::mutate(W = .data$W / denom_W) %>% 
      dplyr::mutate(weight_result = .data$W * .data$nat_weight) %>%
      dplyr::select(-.data$eps_2, -.data$P_gt, -.data$abs_delta_D)
  }
  
  return(list(dat = dat, beta = beta))
  })
  
}
  