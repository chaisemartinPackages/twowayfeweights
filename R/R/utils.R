#' @importFrom rlang sym

printf <- function(...) print(sprintf(...))

fn_ctrl_rename <- function(x) paste("ctrl", x, sep="_")
get_controls_rename <- function(controls) unlist(lapply(controls, fn_ctrl_rename))
fn_treatment_rename <- function(x) paste("OT", x, sep="_")

get_treatments_rename <- function(treatments) {
  unlist(lapply(treatments, fn_treatment_rename))
}

fn_treatment_weight_rename <- function(x) paste("weight_", x, sep = "")
fn_random_weight_rename <- function(x) paste("RW", x, sep="_")
get_random_weight_rename <- function(ws) unlist(lapply(ws, fn_random_weight_rename))


##
# twowayfeweights_rename_var
twowayfeweights_rename_var <- function(df, Y, G, T, D, D0, controls, treatments, random_weights) {
  
  controls_rename <- get_controls_rename(controls)
  treatments_rename <- get_treatments_rename(treatments)
  
  if (length(random_weights) > 0) {
    random_weight_rename <- get_random_weight_rename(random_weights)
    random_weight_df <- df[, random_weights, drop = FALSE]
    # random_weight_df <- df %>% dplyr::select(all_off(random_weights))
    colnames(random_weight_df) <- random_weight_rename
  }
  
  original_names = c(Y, G, T, D, controls, treatments)
  new_names = c("Y", "G", "T", "D", controls_rename, treatments_rename)
  
  if (!is.null(D0)) {
    original_names = c(original_names, D0)
    new_names = c(new_names, "D0")
  }
  
  df <- data.frame(df) %>% dplyr::select_at(dplyr::vars(original_names))
  colnames(df) <- new_names
  
  if (length(random_weights) > 0) {
    df <- cbind(df, random_weight_df)
  }
  
  return(df)
}

# ##
# # twowayfeweights_normalize_var
# twowayfeweights_normalize_var <- function(df, varname) {
#   
#   .data = NULL
#   
#   var = rlang::sym(varname)
#   sdf <- df %>%
#     dplyr::group_by(.data$G, .data$T) %>%
#     dplyr::summarise(tmp_mean_gt = mean(!!var), tmp_sd_gt = stats::sd(!!var))
#   
#   tmp_sd_gt_sum = sum(sdf$tmp_sd_gt, na.rm=TRUE)
#   if (tmp_sd_gt_sum > 0) {
#     df <- df %>% 
#       dplyr::left_join(sdf, by=c("T", "G")) %>%
#       dplyr::mutate(!!var := .data$tmp_mean_gt) %>%
#       dplyr::select(-.data$tmp_mean_gt) %>%
#       dplyr::select(-.data$tmp_sd_gt)
#   }
#   
#   return(list(retcode = (tmp_sd_gt_sum > 0), df = df))
# }


##
# twowayfeweights_transform
twowayfeweights_transform <- function(df, controls, weights, treatments) {
  
  .data = NULL
  
  ret = twowayfeweights_normalize_var(df, "D")
  if (ret$retcode) {
    df <- ret$df
    printf("The treatment variable in the regression varies within some group * period cells.")
    printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
    printf("with a group * period level treatment.")
    printf("The command will replace the treatment by its average value in each group * period.")
    printf("The results below apply to the two-way fixed effects regression with that treatment variable.")
  }
  
  for (control in controls) {
    ret = twowayfeweights_normalize_var(df, control)
    if (ret$retcode) {
      df <- ret$df
      printf("The control variable %s in the regression varies within some group * period cells.", control)
      printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
      printf("with controls apply to group * period level controls.")
      printf("The command will replace replace control variable %s by its average value in each group * period.", control)
      printf("The results below apply to the regression with control variable %s averaged at the group * period level.", control)
    }
  }
  
  for (treatment in treatments) {
    ret = twowayfeweights_normalize_var(df, treatment)
    if (ret$retcode) {
      df <- ret$df
      printf("The other treatment variable %s in the regression varies within some group * period cells.", treatment)
      printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
      printf("with several treatments apply to group * period level controls.")
      printf("The command will replace replace other treatment variable %s by its average value in each group * period.", treatment)
      printf("The results below apply to the regression with other treatment variable %s averaged at the group * period level.", treatment)
    }
  }
  
  if (is.null(weights)) {
    df$weights <- 1
  } else {
    df$weights <- weights
  }
  
  df$Tfactor <- factor(df$T)
  TfactorLevels <- length(levels(df$Tfactor))
  df <- df %>% dplyr::mutate(TFactorNum = as.numeric(factor(.data$Tfactor, labels = seq(1:TfactorLevels))))
  
  return(df)
}


##
# twowayfeweights_filter
twowayfeweights_filter <- function(df, Y, G, T, D, D0, cmd_type, controls, treatments) {
  .data = NULL
  # Remove rows with NA values
  if (cmd_type != "fdTR") {
    df <- df %>%
      dplyr::mutate(tag = rowSums(dplyr::across(.cols = c(Y, G, T, D, controls, treatments), .fns = is.na))) %>%
      dplyr::filter(.data$tag == 0) %>%
      dplyr::select(-.data$tag)
  } else {
    df <- df %>%
      dplyr::mutate(tag1 = rowSums(dplyr::across(.cols = c(D, T, Y), .fns = is.na))) %>%
      dplyr::mutate(tag2 = rowSums(dplyr::across(.cols = c(D0), .fns = is.na))) %>%
      dplyr::filter(.data$tag1 == 0 | .data$tag2 == 0)
    
    if (length(controls) > 0) {
      df <- df %>%
        dplyr::mutate(tag3 = rowSums(dplyr::across(.cols = controls, .fns = is.na))) %>%
        dplyr::filter(.data$tag1 == 1 | .data$tag3 == 0) %>%
        dplyr::select(-.data$tag3)
    }
    df <- df %>% dplyr::select(-.data$tag1, -.data$tag2)
  }
  return(df)
}



##
# twowayfeweights_summarize_weights
twowayfeweights_summarize_weights <- function(df, var_weight) {
  
  weight_plus <- df[[var_weight]][df[[var_weight]] > 0 & !is.na(df[[var_weight]])]
  nr_plus <- length(weight_plus)
  sum_plus <- sum(weight_plus, na.rm = TRUE)
  
  weight_minus <- df[[var_weight]][df[[var_weight]] < 0 & !is.na(df[[var_weight]])]
  nr_minus <- length(weight_minus)
  sum_minus <- sum(weight_minus, na.rm = TRUE)
  
  nr_weights <- nr_plus + nr_minus
  
  return(
    list(
      nr_plus    = nr_plus,
      nr_minus   = nr_minus,
      nr_weights = nr_weights,
      sum_plus   = sum_plus,
      sum_minus  = sum_minus
    )
  )

}

##
# twowayfeweights_test_random_weights
twowayfeweights_test_random_weights <- function(df, random_weights) {
  
  .data = NULL
  
  mat <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(mat) <- c("Coef", "SE", "t-stat", "Correlation")
  df_filtered <- df %>% dplyr::filter(is.finite(.data$W))
  df_filtered_sub <- subset(df_filtered, df_filtered$nat_weight != 0) #Modif. Diego: added extra line to solve note in R CMD Check
 
  for (v in random_weights) {
    formula <- sprintf("%s ~ W", v)
    # rw.lm = estimatr::lm_robust(formula = as.formula(formula), data = df_filtered_sub, weights = df_filtered_sub$nat_weight, clusters = df_filtered_sub$G, se_type = "stata")
    # beta <- rw.lm$coefficients[["W"]]
    # se <- rw.lm$std.error[["W"]]
    # r2 <- rw.lm$r.squared
    rw_lm = fixest::feols(fml = as.formula(formula), data = df_filtered_sub, weights = ~nat_weight, vcov = ~G)
    beta = stats::coef(rw_lm)[["W"]]
    se = sqrt(diag(stats::vcov(rw_lm)))[["W"]]
    r2 = fixest::r2(rw_lm)["r2"]
    
    mat[v, ] <- c(beta, se, beta/se, if (beta > 0) { sqrt(r2) } else { -sqrt(r2) })
  }
  
  return(mat)
}

twowayfeweights_result <- function(df, beta, random_weights) {
  
  .data = NULL

  # Avoid overcounting of positive and negative weights close to 0
  limit_sensitivity <- 10^(-10)
  df$weight_result <- ifelse(df$weight_result < limit_sensitivity & df$weight_result > -limit_sensitivity, 0, df$weight_result)
  ret <- twowayfeweights_summarize_weights(df, "weight_result")

  W_mean <- weighted.mean(df$W, df$nat_weight)
  # Modif. Diego: DoF adjustment to the sd of w_gt
  M <- sum(df$nat_weight != 0)
  W_sd <- sqrt(sum(df$nat_weight * (df$W - W_mean)^2, na.rm = TRUE)) * sqrt(M/(M - 1))
  sensibility <- abs(beta) / W_sd
  
  df_result <- df %>% 
    dplyr::select(.data$T, .data$G, .data$weight_result) %>% 
    dplyr::rename(weight = .data$weight_result)
  
  ret$df_result = df_result
  ret$beta = beta
  ret$sensibility = sensibility
  if (length(random_weights) > 0) {
    ret$mat = twowayfeweights_test_random_weights(df, random_weights)
  }
  
  if (ret$sum_minus < 0) {
    df_sens <- df %>%
      dplyr::filter(.data$weight_result != 0) %>%
      dplyr::arrange(dplyr::desc(.data$W)) %>% 
      dplyr::mutate(P_k = 0, S_k = 0, T_k = 0)
    
    # Modif. Diego: Replaced the previous two loops with build-in routines
    N = nrow(df_sens)
    df_sens <- df_sens[order(df_sens$W, -df_sens$G, -df_sens$T),]
    df_sens$Wsq <- df_sens$nat_weight * df_sens$W^2
    df_sens$P_k <- cumsum(df_sens$nat_weight) 
    df_sens$S_k <- cumsum(df_sens$weight_result) 
    df_sens$T_k <- cumsum(df_sens$Wsq) 
    df_sens <- df_sens[order(-df_sens$W, df_sens$G, df_sens$T),]
    df_sens <- df_sens %>% 
      dplyr::mutate(sens_measure2 = abs(beta) / sqrt(.data$T_k + .data$S_k^2 / (1 - .data$P_k))) %>%
      dplyr::mutate(indicator = as.numeric(.data$W < - .data$S_k / (1 - .data$P_k)))
    df_sens$indicator[1] = 0
    df_sens <- df_sens %>%
      dplyr::mutate(indicator_l = dplyr::lag(.data$indicator, default = -1))
    df_sens <- df_sens %>%
      dplyr::rowwise() %>%
      dplyr::mutate(indicator=max(.data$indicator, .data$indicator_l))
    total_indicator <- sum(df_sens$indicator)
    sensibility2 <- df_sens$sens_measure2[N - total_indicator + 1]
    ret$sensibility2 = sensibility2
  }
  return(ret)
  
}

# ##
# # twowayfeweights_print_results
# twowayfeweights_print_results <- function(cmd_type, r, d, summary_measures, twfe, random_weights) {
#   treat = dplyr::case_when(cmd_type == "feTR" || cmd_type == "fdTR" ~ "ATT",
#                     cmd_type == "feS" || cmd_type == "fdS" ~ "LATE", 
#                     TRUE ~ "BLANK")
#   assumption = dplyr::case_when(cmd_type == "feTR" || cmd_type == "fdTR" ~ "Under the common trends assumption",
#                          cmd_type == "feS" || cmd_type == "fdS" ~ "Under the common trends, treatment monotonicity, and if groups' treatment effect does not change over time", 
#                          TRUE ~ "BLANK")
#   
#   tot_weights <- r$nr_plus + r$nr_minus
#   tot_sums <- round(r$sum_plus + r$sum_minus, 4)
#   
#   cat("\n")
#   print(paste(rep("-",68), collapse= ""))
#   printf("%s, ", assumption)
#   printf("beta estimates a weighted sum of %d %ss.", r$nr_weights, treat)
#   printf("%d %s receive a positive weight, and %d receive a negative weight.", r$nr_plus, treat, r$nr_minus)
#   print(paste(rep("-",68), collapse= ""))
#   
#   ##Doulo
#   print_out <- huxtable::hux(c(paste0("Treat. var: ", d),"Positive weights", "Negative weights", "Total")
#                 ,c(paste0(treat, "s"), round(r$nr_plus, 2), round(r$nr_minus, 2), tot_weights),
#                 c(paste0("\U03A3", " weights"), round(r$sum_plus, 4), round(r$sum_minus, 4), tot_sums),
#                 add_colnames = FALSE
#   )
#   
#   environment_res <- list(r$sum_minus, r$sensibility, r$beta)
#   environment_names <- c("sum_neg_w", "lb_se_te", "beta")
#   
#   print_out = print_out %>% 
#     huxtable::set_all_padding(4) %>% 
#     huxtable::set_outer_padding(0) %>% 
#     #set_number_format(2) %>% 
#     huxtable::set_bottom_border(row = 1, col = huxtable::everywhere) %>% 
#     huxtable::set_top_border(row = 1, col = huxtable::everywhere) %>% 
#     huxtable::set_bottom_border(row = 3, col = huxtable::everywhere) %>% 
#     huxtable::set_width(2) %>% 
#     huxtable::set_align(huxtable::everywhere, 2:3, "center") %>% 
#     huxtable::set_bold(1, huxtable::everywhere) %>% 
#     huxtable::set_bottom_border(row = 4, col = huxtable::everywhere)
#   
#   colnames(print_out) = NULL
#   huxtable::print_screen(print_out) 
#   
#   ###If we want to add the option export_excel
#   # export_path  = "export_path.xlsx"
#   # if(length(export_path)>0){
#   #   quick_xlsx(print_out %>%set_width(1.5), file = export_path)
#   # }
#   
#   if (!is.null(summary_measures)) {
#       subscr <- substr(cmd_type, 1, 2)
#       print("Summary Measures:")
#       printf("TWFE Coefficient (\U03B2_%s): %.4f", subscr, twfe)
#       printf("min \U03C3(\U0394) compatible with \U03B2_%s and \U0394_TR = 0: %.4f", subscr, r$sensibility)
#     if (r$sum_minus < 0) {
#       printf("min \U03C3(\U0394) compatible with \U03B2_%s and \U0394_TR of a different sign: %.4f", subscr, r$sensibility2)
#       environment_res <- append(environment_res, r$sensibility2)
#       environment_names <- c(environment_names, "lb_se_te2") 
#     } 
#     print("Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)")
#   }
#   
#   if (length(random_weights) > 0) {
#     cat("\n")
#     print("Test random weights")
#     print(r$mat)
#     environment_res <- append(environment_res, r$mat)
#     environment_names <- c(environment_names, "randomweightstest")
#   }
#   
#   cat("\n")
#   print("The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N. 101043899).")
#   
#   names(environment_res) <- environment_names
#   list2env(environment_res, envir = .GlobalEnv) ## BIG OOF
#   r$df_result
# }

# ##
# # twowayfeweights_calculate_fetr_other_treatment
# twowayfeweights_calculate_fetr_other_treatment <- function(df, controls, treatments) {
#   mean_D <- weighted.mean(df$D, df$weights, na.rm = TRUE)
#   obs <- sum(df$weights)
#   gdf <- df %>% dplyr::group_by(.data$G, .data$T) %>% dplyr::summarise(P_gt = sum(.data$weights))
#   df <- df %>% 
#     dplyr::left_join(gdf, by=c("T", "G")) %>% 
#     dplyr::mutate(P_gt = .data$P_gt / obs) %>% 
#     dplyr::mutate(nat_weight = .data$P_gt * .data$D / mean_D)
#   
#   vars = c(controls, treatments)
#   formula = paste(vars, collapse = " + ")
#   formula = paste("D ~ ", formula, sep = "")
#   formula = paste(formula, " | G + Tfactor", sep = "")
#   denom.lm <- fixest::feols(as.formula(formula), data = df, weights = df$weights)
#   df$eps_1 <- df$D - predict(denom.lm, df)
#   df$eps_1_E_D_gt <- df$eps_1 * df$D
#   denom_W <- mean(df$eps_1_E_D_gt, na.rm = TRUE)
#   
#   df <- df %>% 
#     dplyr::mutate(W = .data$eps_1 * mean_D / denom_W) %>% 
#     dplyr::mutate(weight_result = .data$W * .data$nat_weight)
#   
#   for (treatment in treatments) {
#     varname = fn_treatment_weight_rename(treatment)
#     df <- df %>% dplyr::mutate(!!rlang::sym(varname) := .data$W * .data$P_gt * !!rlang::sym(treatment) / mean_D)
#   }
#   
#   df <- df %>% dplyr::select(-.data$eps_1, -.data$P_gt)
#   
#   formula = paste(vars, collapse = " + ")
#   formula = paste("Y ~ D + ", formula, sep = "")
#   formula = paste(formula, " | G + Tfactor", sep = "")
#   beta.lm <- fixest::feols(as.formula(formula), data = df, weights = df$weights)
#   beta <- as.numeric(coef(beta.lm)["D"])
#   
#   return(list(df = df, beta = beta))
# }

##
#twowayfeweights_result_other_treatment
twowayfeweights_result_other_treatment <- function(df, treatments, beta, random_weights) {
  
  .data = NULL
  
  columns <- c("T", "G", "weight_result")
  ret <- twowayfeweights_summarize_weights(df, "weight_result")
  
  if (length(random_weights) > 0) {
    ret$mat = twowayfeweights_test_random_weights(df, random_weights)
  }
  
  for (treatment in treatments) {
    varname = fn_treatment_weight_rename(treatment)
    columns <- c(columns, varname)
    ret2 <- twowayfeweights_summarize_weights(df, varname)
    ret[[treatment]] <- ret2
  }
  df_result <- df %>% 
    dplyr::select_at(dplyr::vars(columns)) %>% 
    dplyr::rename(weight = .data$weight_result)
  
  ret$beta = beta
  ret[["df_result"]] <- df_result
  return(ret)
}


##
# twowayfeweights_print_result_other_treatment
twowayfeweights_print_result_other_treatment <- function(r, treatments, d, twfe, random_weights) {

  tot_weights <- r$nr_plus + r$nr_minus
  tot_sums <- round(r$sum_plus + r$sum_minus, 4)
  

  cat("\n")
  print(paste(rep("-",75), collapse= ""))
  printf("Under the common trends assumption, beta estimates the sum of several terms.")
  printf("The first term is a weighted sum of %d ATTs of the treatment.", r$nr_weights)
  printf("%d ATTs receive a positive weight, and %d receive a negative weight.", r$nr_plus, r$nr_minus)
  print(paste(rep("-",75), collapse= ""))
  
  ##Doulo
  print_out <- huxtable::hux(c(paste0("Treat. var: ", d),"Positive weights", "Negative weights", "Total")
                ,c(paste0("ATT", "s"), round(r$nr_plus, 2), round(r$nr_minus, 2), tot_weights),
                c(paste0("\U03A3", " weights"), round(r$sum_plus, 4), round(r$sum_minus, 4), tot_sums),
                add_colnames = FALSE
  )
  environment_res <- list(r$sum_minus, twfe)
  environment_names <- c("sum_neg_w", "beta")
  
  print_out = print_out %>% 
    huxtable::set_all_padding(4) %>% 
    huxtable::set_outer_padding(0) %>% 
    #set_number_format(2) %>% 
    huxtable::set_bottom_border(row = 1, col = huxtable::everywhere) %>% 
    huxtable::set_top_border(row = 1, col = huxtable::everywhere) %>% 
    huxtable::set_bottom_border(row = 3, col = huxtable::everywhere) %>% 
    huxtable::set_width(2) %>% 
    huxtable::set_align(huxtable::everywhere, 2:3, "center") %>% 
    huxtable::set_bold(1, huxtable::everywhere) %>% 
    huxtable::set_bottom_border(row = 4, col = huxtable::everywhere)
  
  colnames(print_out) = NULL
  huxtable::print_screen(print_out) 
  
  ##
  
  j <- 1
  for (treatment in treatments) {
    r2 = r[[treatment]]
    ot_treat <- sub("OT_", "", treatment) 
    tot_weights <- r2$nr_plus + r2$nr_minus
    tot_sums <- round(r2$sum_plus + r2$sum_minus, 4)
    
    cat("\n")
    print(paste(rep("-",108), collapse= ""))
    printf("The next term is a weighted sum of %d ATTs of treatment %s included in the other_treatments option.", r2$nr_weights, treatment)
    printf("%d ATTs receive a positive weight, and %d receive a negative weight.", r2$nr_plus, r2$nr_minus)
    print(paste(rep("-",108), collapse= ""))
    
    ##Doulo
    print_out <- huxtable::hux(c(paste0("Other treat.: ", ot_treat),"Positive weights", "Negative weights", "Total")
                  ,c(paste0("ATT", "s"), round(r2$nr_plus, 2), round(r2$nr_minus, 2), tot_weights),
                  c(paste0("\U03A3", " weights"), round(r2$sum_plus, 4), round(r2$sum_minus, 4), tot_sums),
                  add_colnames = FALSE
    )
    
    oth_title <- sprintf("sum_neg_w_othertreatment%.0f", j)
    environment_res <- append(environment_res, r2$sum_minus)
    environment_names <- c(environment_names, oth_title)

    print_out = print_out %>% 
      huxtable::set_all_padding(4) %>% 
      huxtable::set_outer_padding(0) %>% 
      #set_number_format(2) %>% 
      huxtable::set_bottom_border(row = 1, col = huxtable::everywhere) %>% 
      huxtable::set_top_border(row = 1, col = huxtable::everywhere) %>% 
      huxtable::set_bottom_border(row = 3, col = huxtable::everywhere) %>% 
      huxtable::set_width(2) %>% 
      huxtable::set_align(huxtable::everywhere, 2:3, "center") %>% 
      huxtable::set_bold(1, huxtable::everywhere) %>% 
      huxtable::set_bottom_border(row = 4, col = huxtable::everywhere)
    
    colnames(print_out) = NULL
    huxtable::print_screen(print_out) 
  }
  
  if (length(random_weights) > 0) {
    cat("\n")
    print("Test random weights")
    print(r$mat)
    environment_res <- append(environment_res, r$mat)
    environment_names <- c(environment_names, "randomweightstest")
  }
  cat("\n")
  print("The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N. 101043899).")
  
  names(environment_res) <- environment_names
  list2env(environment_res, envir = .GlobalEnv) ## BIG OOF

  return(r$df_result)
}
