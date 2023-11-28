#Install and Import all the packages
packages <- c("data.table", "dplyr", "logr", "estimatr", "huxtable", "magrittr", "fixest")

ds_libraries <- function(packages){
  for(package in packages){
    if(!require(package, character.only = TRUE)){
      install.packages(package, dependencies = TRUE)
    }
    #Load package
    library(package, character.only = TRUE)
  }
}
ds_libraries(packages)

options(dplyr.summarise.inform = FALSE)
suppressWarnings({
setFixest_notes(FALSE)

printf <- function(...)print(sprintf(...))
fn_ctrl_rename <- function(x) { paste("ctrl", x, sep="_") }
get_controls_rename <- function(controls) { unlist(lapply(controls, fn_ctrl_rename)) }
fn_treatment_rename <- function(x) { paste("OT", x, sep="_")}
get_treatments_rename <- function(treatments) { unlist(lapply(treatments, fn_treatment_rename)) }
fn_treatment_weight_rename <- function(x) { paste("weight_", x, sep = "") }
fn_random_weight_rename <- function(x) { paste("RW", x, sep="_")}
get_random_weight_rename <- function(ws) { unlist(lapply(ws, fn_random_weight_rename)) }

twowayfeweights_rename_var <- function(df, Y, G, T, D, D0, controls, treatments, random_weights) {
  controls_rename <- get_controls_rename(controls)
  treatments_rename <- get_treatments_rename(treatments)
  
  if (length(random_weights) > 0) {
    random_weight_rename <- get_random_weight_rename(random_weights)
    #random_weight_df <- df[random_weights]
    random_weight_df <- df %>% select(all_of(random_weights))
    colnames(random_weight_df) <- random_weight_rename
    
  }
  
  
  original_names = c(Y, G, T, D, controls, treatments)
  new_names = c("Y", "G", "T", "D", controls_rename, treatments_rename)
  
  if (!is.null(D0)) {
    original_names = c(original_names, D0)
    new_names = c(new_names, "D0")
  }
  
  df <- data.frame(df) %>% select_at(vars(original_names))
  colnames(df) <- new_names
  
  if (length(random_weights) > 0) {
    df <- cbind(df, random_weight_df)
  }
  
  df
}

twowayfeweights_normalize_var <- function(df, varname){
  var = sym(varname)
  sdf <- df %>%
    group_by(G, T) %>%
    summarise(tmp_mean_gt = mean(!!var), tmp_sd_gt = sd(!!var))
  
  tmp_sd_gt_sum = sum(sdf$tmp_sd_gt, na.rm=TRUE)
  if (tmp_sd_gt_sum > 0) {
    df <- df %>% 
      left_join(sdf, by=c("T", "G")) %>%
      mutate(!!var := tmp_mean_gt) %>%
      select(-tmp_mean_gt) %>%
      select(-tmp_sd_gt)
  }
  
  list(retcode = (tmp_sd_gt_sum > 0), df = df)
}

twowayfeweights_transform <- function(df, controls, weights, treatments) {
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
  df <- df %>% mutate(TFactorNum = as.numeric(factor(Tfactor, labels = seq(1:TfactorLevels))))
  
  df
}

twowayfeweights_filter <- function(df, Y, G, T, D, D0, cmd_type, controls, treatments) {
  # Remove rows with NA values
  if (cmd_type != "fdTR") {
    df <- df %>%
      mutate(tag = rowSums(across(.cols = c(Y, G, T, D, controls, treatments), .fns = is.na))) %>%
      filter(tag == 0) %>%
      select(-tag)
  } else {
    df <- df %>%
      mutate(tag1 = rowSums(across(.cols = c(D, T, Y), .fns = is.na))) %>%
      mutate(tag2 = rowSums(across(.cols = c(D0), .fns = is.na))) %>%
      filter(tag1 == 0 | tag2 == 0)
    
    if (length(controls) > 0) {
      df <- df %>%
        mutate(tag3 = rowSums(across(.cols = controls, .fns = is.na))) %>%
        filter(tag1 == 1 | tag3 == 0) %>%
        select(-tag3)
    }
    df <- df %>% select(-tag1, -tag2)
  }
  df
}

twowayfeweights_calculate_fetr <- function(df, controls) {
  mean_D <- weighted.mean(df$D, df$weights, na.rm = TRUE)
  obs <- sum(df$weights)
  gdf <- df %>% group_by(G, T) %>% summarise(P_gt = sum(weights)) #removed 
  df <- df %>% 
    left_join(gdf, by=c("T", "G")) %>% 
    mutate(P_gt = P_gt / obs) %>% 
    mutate(nat_weight = P_gt * D / mean_D)
  
  if (length(controls) == 0) {
    formula = "D ~ 1"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("D ~ ", formula, sep = "")
  }
  formula = paste(formula, " | G + Tfactor", sep = "")
  denom.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  df$eps_1 <- residuals(denom.lm)
  df$eps_1_E_D_gt <- df$eps_1 * df$D
  denom_W <- weighted.mean(df$eps_1_E_D_gt, df$weights, na.rm = TRUE)
  
  df <- df %>% 
    mutate(W = eps_1 * mean_D / denom_W) %>% 
    mutate(weight_result = W * nat_weight) %>%
    select(-eps_1, -P_gt)
  
  if (length(controls) == 0) {
    formula = "Y ~ D"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("Y ~ D + ", formula, sep = "")
  }
  formula = paste(formula, " | G + Tfactor", sep = "")
  beta.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  beta <- as.numeric(coef(beta.lm)["D"])
  
  # * Keeping only one observation in each group * period cell
  # This should be done after this function
  # bys `group' `time': gen group_period_unit=(_n==1)	
  # 	drop if group_period_unit==0
  # 	drop group_period_unit
  df <- df %>%
    group_by(G, Tfactor) %>%
    filter(row_number(D) == 1)
  
  list(df = df, beta = beta)
}

twowayfeweights_calculate_fdtr <- function(df, controls) {
  mean_D0 <- weighted.mean(df$D0, df$weights, na.rm = TRUE)
  obs <- sum(df$weights)
  gdf <- df %>% group_by(G, T) %>% summarise(P_gt = sum(weights))
  df <- df %>% 
    left_join(gdf, by=c("G", "T")) %>% 
    mutate(P_gt = P_gt / obs) %>% 
    mutate(nat_weight = P_gt * D0 / mean_D0)
  
  if (length(controls) == 0) {
    formula = "D ~ 1"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("D ~ ", formula, sep = "")
  }
  formula = paste(formula, " | Tfactor", sep = "")
  denom.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  df$eps_2 <- residuals(denom.lm, na.rm = FALSE)
  # df$eps_2 <- df$D - predict(denom.lm, df)
  
  df <- df %>% mutate(eps_2 = ifelse(is.na(eps_2), 0, eps_2))
  

  if (length(controls) == 0) {
    formula = "Y ~ D"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("Y ~ D + ", formula, sep = "")
  }
  formula = paste(formula, " | Tfactor", sep = "")
  beta.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  beta <- as.numeric(coef(beta.lm)["D"])
  
  df <- df %>% 
    arrange(G, TFactorNum) %>%
    group_by(G) %>% 
    mutate(w_tilde_2 = ifelse(TFactorNum + 1 == lead(TFactorNum), eps_2 - lead(eps_2) * (lead(P_gt) / P_gt), NA)) %>%
    mutate(w_tilde_2 = ifelse(is.na(w_tilde_2) | is.infinite(w_tilde_2), eps_2, w_tilde_2)) %>%
    mutate(w_tilde_2_E_D_gt = w_tilde_2 * D0)
  
  denom_W <- weighted.mean(df$w_tilde_2_E_D_gt, df$P_gt, na.rm = TRUE)
  df <- df %>% 
    mutate(W = w_tilde_2 * mean_D0 / denom_W) %>% 
    mutate(weight_result = W * nat_weight)
  df <- df %>%
    select(-eps_2, -P_gt, -w_tilde_2, -w_tilde_2_E_D_gt)
  list(df = df, beta = beta)
}

twowayfeweights_calculate_fes <- function(df, controls) {
  obs <- sum(df$weights)
  gdf <- df %>% group_by(G, T) %>% summarise(P_gt = sum(weights))
  df <- df %>% 
    left_join(gdf, by=c("T", "G")) %>% 
    mutate(P_gt = P_gt / obs)
  
  if (length(controls) == 0) {
    formula = "D ~ 1"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("D ~ ", formula, sep = "")
  }
  formula = paste(formula, " | G + Tfactor", sep = "")
  denom.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  df$eps_1 <- residuals(denom.lm)
  # df$eps_1 <- df$D - predict(denom.lm, df)
  
  df <- df %>% 
    mutate(eps_1_weight = eps_1 * weights) %>%
    arrange(G, Tfactor) %>%
    group_by(G) %>%
    mutate(E_eps_1_g_ge_aux = rev(cumsum(rev(eps_1_weight)))) %>%
    mutate(weights_aux = rev(cumsum(rev(weights)))) %>%
    mutate(E_eps_1_g_ge = E_eps_1_g_ge_aux / weights_aux)
  
  if (length(controls) == 0) {
    formula = "Y ~ D"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("Y ~ D + ", formula, sep = "")
  }
  formula = paste(formula, " | G + Tfactor", sep = "")
  beta.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  beta <- as.numeric(coef(beta.lm)["D"])
  
  # * Keeping only one observation in each group * period cell
  #   bys `group' `time': gen group_period_unit=(_n==1)	
  # 	drop if group_period_unit==0
  # 	drop group_period_unit
  
  # df <- df %>% 
  #   group_by(G, Tfactor) %>%
  #   summarize(P_gt, nat_weight)
  
  df <- df %>% 
    arrange(G, Tfactor) %>%
    group_by(G) %>% 
    mutate(delta_D = ifelse(TFactorNum - 1 == lag(TFactorNum), D - lag(D), NA)) %>%
    filter(!is.na(delta_D)) %>%
    mutate(abs_delta_D = abs(delta_D)) %>%
    mutate(s_gt = case_when(delta_D > 0 ~ 1,
                            delta_D < 0 ~ -1,
                            TRUE ~ 0)) %>%
    mutate(nat_weight = P_gt * abs_delta_D)
  
  P_S = sum(df$nat_weight, na.rm = TRUE)
  df <- df %>% 
    mutate(nat_weight = nat_weight / P_S) %>%
    mutate(om_tilde_1 = s_gt * E_eps_1_g_ge / P_gt)
  
  denom_W = weighted.mean(df$om_tilde_1, df$nat_weight, na.rm = TRUE)
  df <- df %>%
    mutate(W = om_tilde_1 / denom_W) %>%
    mutate(weight_result = W * nat_weight) %>%
    select(-eps_1, -P_gt, -om_tilde_1, -E_eps_1_g_ge,
           -E_eps_1_g_ge_aux, -weights_aux, -abs_delta_D, -delta_D)
  list(df = df, beta = beta)
}

twowayfeweights_calculate_fds <- function(df, controls) {
  obs <- sum(df$weights)
  gdf <- df %>% group_by(G, T) %>% summarise(P_gt = sum(weights))
  df <- df %>% 
    left_join(gdf, by=c("T", "G")) %>% 
    mutate(P_gt = P_gt / obs)
  
  if (length(controls) == 0) {
    formula = "D ~ 1"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("D ~ ", formula, sep = "")
  }
  formula = paste(formula, " | Tfactor", sep = "")
  denom.lm <- feols(as.formula(formula), data = subset(df, weights != 0), weights = df$weights)
  df$eps_2 <- residuals(denom.lm)
  # df$eps_2 <- df$D - predict(denom.lm, df)
  
  if (length(controls) == 0) {
    formula = "Y ~ D"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("Y ~ D + ", formula, sep = "")
  }
  formula = paste(formula, " | Tfactor", sep = "")
  beta.lm <- feols(as.formula(formula), data = subset(df, weights != 0), weights = df$weights)
  beta <- as.numeric(coef(beta.lm)["D"])
  
  # * Keeping only one observation in each group * period cell
  #   bys `group' `time': gen group_period_unit=(_n==1)	
  # 	drop if group_period_unit==0
  # 	drop group_period_unit
  
  # df <- df %>% 
  #   group_by(G, Tfactor) %>%
  #   summarize(P_gt, nat_weight)
  
  df <- df %>%
    mutate(s_gt = case_when(D > 0 ~ 1,
                            D < 0 ~ -1,
                            TRUE ~ 0)) %>%
    mutate(abs_delta_D = abs(D)) %>%
    mutate(nat_weight = P_gt * abs_delta_D)
  
  P_S = sum(df$nat_weight)
  df <- df %>% 
    mutate(nat_weight = nat_weight / P_S) %>%
    mutate(W = s_gt * eps_2)
  denom_W = weighted.mean(df$W, df$nat_weight, na.rm = TRUE)
  df <- df %>% 
    mutate(W = W / denom_W) %>% 
    mutate(weight_result = W * nat_weight) %>%
    select(-eps_2, -P_gt, -abs_delta_D)
  
  list(df = df, beta = beta)
}

twowayfeweights_summarize_weights <- function(df, var_weight) {
  
  weight_plus <- df[[var_weight]][df[[var_weight]] > 0 & !is.na(df[[var_weight]])]
  nr_plus <- length(weight_plus)
  sum_plus <- sum(weight_plus, na.rm = TRUE)
  
  weight_minus <- df[[var_weight]][df[[var_weight]] < 0 & !is.na(df[[var_weight]])]
  nr_minus <- length(weight_minus)
  sum_minus <- sum(weight_minus, na.rm = TRUE)
  
  nr_weights <- nr_plus + nr_minus
  
  list(nr_plus = nr_plus, nr_minus = nr_minus, nr_weights = nr_weights, 
       sum_plus = sum_plus, sum_minus = sum_minus)
}

twowayfeweights_test_random_weights <- function(df, random_weights) {
  mat <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(mat) <- c("Coef", "SE", "t-stat", "Correlation")
  df_filtered <- df %>% filter(is.finite(W))
 
  for (v in random_weights) {
    formula <- sprintf("%s ~ W", v)
    rw.lm = lm_robust(formula =as.formula(formula), data = subset(df_filtered, nat_weight != 0), weights = nat_weight, clusters = G, se_type = "stata")
    beta <- rw.lm$coefficients[["W"]]
    se <- rw.lm$std.error[["W"]]
    r2 <- rw.lm$r.squared
    
    mat[v, ] <- c(beta, se, beta/se, if (beta > 0) { sqrt(r2) } else { -sqrt(r2) })
  }
  
  return(mat)
}

twowayfeweights_result <- function(df, beta, random_weights) {

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
    select(T, G, weight_result) %>% 
    rename(weight = weight_result)
  
  ret$df_result = df_result
  ret$beta = beta
  ret$sensibility = sensibility
  if (length(random_weights) > 0) {
    ret$mat = twowayfeweights_test_random_weights(df, random_weights)
  }
  
  if (ret$sum_minus < 0) {
    df_sens <- df %>%
      filter(weight_result != 0) %>%
      arrange(desc(W)) %>% 
      mutate(P_k = 0, S_k = 0, T_k = 0)
    
    # Modif. Diego: Replaced the previous two loops with build-in routines
    N = nrow(df_sens)
    df_sens <- df_sens[order(df_sens$W, -df_sens$G, -df_sens$T),]
    df_sens$Wsq <- df_sens$nat_weight * df_sens$W^2
    df_sens$P_k <- cumsum(df_sens$nat_weight) 
    df_sens$S_k <- cumsum(df_sens$weight_result) 
    df_sens$T_k <- cumsum(df_sens$Wsq) 
    df_sens <- df_sens[order(-df_sens$W, df_sens$G, df_sens$T),]
    df_sens <- df_sens %>% 
      mutate(sens_measure2 = abs(beta) / sqrt(T_k + S_k^2 / (1 - P_k))) %>%
      mutate(indicator = as.numeric(W < - S_k / (1 - P_k)))
    df_sens$indicator[1] = 0
    df_sens <- df_sens %>%
      mutate(indicator_l = lag(indicator, default = -1))
    df_sens <- df_sens %>%
      rowwise() %>%
      mutate(indicator=max(indicator, indicator_l))
    total_indicator <- sum(df_sens$indicator)
    sensibility2 <- df_sens$sens_measure2[N - total_indicator + 1]
    ret$sensibility2 = sensibility2
  }
  ret
  
}

twowayfeweights_print_results <- function(cmd_type, r, d, summary_measures, twfe, random_weights) {
  treat = case_when(cmd_type == "feTR" || cmd_type == "fdTR" ~ "ATT",
                    cmd_type == "feS" || cmd_type == "fdS" ~ "LATE", 
                    TRUE ~ "BLANK")
  assumption = case_when(cmd_type == "feTR" || cmd_type == "fdTR" ~ "Under the common trends assumption",
                         cmd_type == "feS" || cmd_type == "fdS" ~ "Under the common trends, treatment monotonicity, and if groups' treatment effect does not change over time", 
                         TRUE ~ "BLANK")
  
  tot_weights <- r$nr_plus + r$nr_minus
  tot_sums <- round(r$sum_plus + r$sum_minus, 4)
  
  cat("\n")
  print(paste(rep("-",68), collapse= ""))
  printf("%s, ", assumption)
  printf("beta estimates a weighted sum of %d %ss.", r$nr_weights, treat)
  printf("%d %s receive a positive weight, and %d receive a negative weight.", r$nr_plus, treat, r$nr_minus)
  print(paste(rep("-",68), collapse= ""))
  
  ##Doulo
  print_out <- hux(c(paste0("Treat. var: ", d),"Positive weights", "Negative weights", "Total")
                ,c(paste0(treat, "s"), round(r$nr_plus, 2), round(r$nr_minus, 2), tot_weights),
                c(paste0("\U03A3", " weights"), round(r$sum_plus, 4), round(r$sum_minus, 4), tot_sums),
                add_colnames = FALSE
  )
  
  environment_res <- list(r$sum_minus, r$sensibility, r$beta)
  environment_names <- c("sum_neg_w", "lb_se_te", "beta")
  
  print_out = print_out %>% 
    set_all_padding(4) %>% 
    set_outer_padding(0) %>% 
    #set_number_format(2) %>% 
    set_bottom_border(row = 1, col = everywhere) %>% 
    set_top_border(row = 1, col = everywhere) %>% 
    set_bottom_border(row = 3, col = everywhere) %>% 
    set_width(2) %>% 
    set_align(everywhere, 2:3, "center") %>% 
    set_bold(1, everywhere) %>% 
    set_bottom_border(row = 4, col = everywhere)
  
  colnames(print_out) = NULL
  print_screen(print_out) 
  
  ###If we want to add the option export_excel
  # export_path  = "export_path.xlsx"
  # if(length(export_path)>0){
  #   quick_xlsx(print_out %>%set_width(1.5), file = export_path)
  # }
  
  if (!is.null(summary_measures)) {
      subscr <- substr(cmd_type, 1, 2)
      print("Summary Measures:")
      printf("TWFE Coefficient (\U03B2_%s): %.4f", subscr, twfe)
      printf("min \U03C3(\U0394) compatible with \U03B2_%s and \U0394_TR = 0: %.4f", subscr, r$sensibility)
    if (r$sum_minus < 0) {
      printf("min \U03C3(\U0394) compatible with \U03B2_%s and \U0394_TR of a different sign: %.4f", subscr, r$sensibility2)
      environment_res <- append(environment_res, r$sensibility2)
      environment_names <- c(environment_names, "lb_se_te2") 
    } 
    print("Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)")
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
  list2env(environment_res, envir = .GlobalEnv)
  r$df_result
}


twowayfeweights_calculate_fetr_other_treatment <- function(df, controls, treatments) {
  mean_D <- weighted.mean(df$D, df$weights, na.rm = TRUE)
  obs <- sum(df$weights)
  gdf <- df %>% group_by(G, T) %>% summarise(P_gt = sum(weights))
  df <- df %>% 
    left_join(gdf, by=c("T", "G")) %>% 
    mutate(P_gt = P_gt / obs) %>% 
    mutate(nat_weight = P_gt * D / mean_D)
  
  vars = c(controls, treatments)
  formula = paste(vars, collapse = " + ")
  formula = paste("D ~ ", formula, sep = "")
  formula = paste(formula, " | G + Tfactor", sep = "")
  denom.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  df$eps_1 <- df$D - predict(denom.lm, df)
  df$eps_1_E_D_gt <- df$eps_1 * df$D
  denom_W <- mean(df$eps_1_E_D_gt, na.rm = TRUE)
  
  df <- df %>% 
    mutate(W = eps_1 * mean_D / denom_W) %>% 
    mutate(weight_result = W * nat_weight)
  
  for (treatment in treatments) {
    varname = fn_treatment_weight_rename(treatment)
    df <- df %>% mutate(!!sym(varname) := W * P_gt * !!sym(treatment) / mean_D)
  }
  
  df <- df %>% select(-eps_1, -P_gt)
  
  formula = paste(vars, collapse = " + ")
  formula = paste("Y ~ D + ", formula, sep = "")
  formula = paste(formula, " | G + Tfactor", sep = "")
  beta.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  beta <- as.numeric(coef(beta.lm)["D"])
  
  list(df = df, beta = beta)
}

twowayfeweights_result_other_treatment <- function(df, treatments, beta, random_weights) {
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
    select_at(vars(columns)) %>% 
    rename(weight = weight_result)
  
  ret$beta = beta
  ret[["df_result"]] <- df_result
  ret
}

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
  print_out <- hux(c(paste0("Treat. var: ", d),"Positive weights", "Negative weights", "Total")
                ,c(paste0("ATT", "s"), round(r$nr_plus, 2), round(r$nr_minus, 2), tot_weights),
                c(paste0("\U03A3", " weights"), round(r$sum_plus, 4), round(r$sum_minus, 4), tot_sums),
                add_colnames = FALSE
  )
  environment_res <- list(r$sum_minus, twfe)
  environment_names <- c("sum_neg_w", "beta")
  
  print_out = print_out %>% 
    set_all_padding(4) %>% 
    set_outer_padding(0) %>% 
    #set_number_format(2) %>% 
    set_bottom_border(row = 1, col = everywhere) %>% 
    set_top_border(row = 1, col = everywhere) %>% 
    set_bottom_border(row = 3, col = everywhere) %>% 
    set_width(2) %>% 
    set_align(everywhere, 2:3, "center") %>% 
    set_bold(1, everywhere) %>% 
    set_bottom_border(row = 4, col = everywhere)
  
  colnames(print_out) = NULL
  print_screen(print_out) 
  
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
    print_out <- hux(c(paste0("Other treat.: ", ot_treat),"Positive weights", "Negative weights", "Total")
                  ,c(paste0("ATT", "s"), round(r2$nr_plus, 2), round(r2$nr_minus, 2), tot_weights),
                  c(paste0("\U03A3", " weights"), round(r2$sum_plus, 4), round(r2$sum_minus, 4), tot_sums),
                  add_colnames = FALSE
    )
    
    oth_title <- sprintf("sum_neg_w_othertreatment%.0f", j)
    environment_res <- append(environment_res, r2$sum_minus)
    environment_names <- c(environment_names, oth_title)

    print_out = print_out %>% 
      set_all_padding(4) %>% 
      set_outer_padding(0) %>% 
      #set_number_format(2) %>% 
      set_bottom_border(row = 1, col = everywhere) %>% 
      set_top_border(row = 1, col = everywhere) %>% 
      set_bottom_border(row = 3, col = everywhere) %>% 
      set_width(2) %>% 
      set_align(everywhere, 2:3, "center") %>% 
      set_bold(1, everywhere) %>% 
      set_bottom_border(row = 4, col = everywhere)
    
    colnames(print_out) = NULL
    print_screen(print_out) 
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
  list2env(environment_res, envir = .GlobalEnv)

  r$df_result
}


#' twowayfeweights
#'
#' @param df the data frame for input
#' @param Y the name of Y variable
#' @param G the name of group variable
#' @param T the name of time variable
#' @param D the name of treatment variable
#' @param cmd_type the type of command, including fetr, fdtr, fes, fds
#' @param D0 the name of the mean of the treatment in group g and at period t
#' @param controls the list of names of control variables, empty if not specified
#' @param weights a column of data that replaces the default weight
#' @param other_treatments the list of other treatment variables to include in the regression other than D
#' @param test_random_weights weights when this option is specified, the command estimates the correlation
#'                            between each variable in varlist and the weights
#' @covariance
#' @average_effect
#' @param parallel parallelly perform bootstrap
#'
#' @return twowayfeweights returns data frame that contains the following columns
#'         T: time variable
#'         G: group variable
#'         weight: the result of the weight
#'
#' @export

twowayfeweights <- function(df, Y, G, T, D, type, D0 = NULL, summary_measures = NULL, controls = c(), weights = NULL, other_treatments = c(), test_random_weights = c(), path = NULL) {
  if (!is.null(D0) && type != "fdTR") {
    printf("Type fdTR requires D0 defined")
    return(c())
  }
  
  if (length(other_treatments) > 0 && type != "feTR") {
    printf("When the other_treatments option is specified, you need to specify the type(feTR) option.")
    return(c())
  }
  
  controls_rename <- get_controls_rename(controls)
  treatments_rename <- get_treatments_rename(other_treatments)
  random_weight_rename <- get_random_weight_rename(test_random_weights)
  
  df_renamed <- twowayfeweights_rename_var(df, Y, G, T, D, D0, controls, other_treatments, test_random_weights)
  df_transformed <- twowayfeweights_transform(df_renamed, controls_rename, weights, treatments_rename)
  df_filtered <- twowayfeweights_filter(df_transformed, Y, G, T, D, D0, type, controls_rename, treatments_rename)

  if (length(other_treatments) == 0) {
    res <- if (type == "feTR") {
      twowayfeweights_calculate_fetr(df_filtered, controls_rename)
    } else if (type == "fdTR") {
      twowayfeweights_calculate_fdtr(df_filtered, controls_rename)
    } else if (type == "feS") {
      twowayfeweights_calculate_fes(df_filtered, controls_rename)
    } else if (type == "fdS") {
      twowayfeweights_calculate_fds(df_filtered, controls_rename)
    }
    
    res <- twowayfeweights_result(res$df, res$beta, random_weight_rename)
    df_result <- twowayfeweights_print_results(type, res, D, summary_measures, res$beta, random_weight_rename)
  } else {
    res <- twowayfeweights_calculate_fetr_other_treatment(df_filtered, controls_rename, treatments_rename)
    res <- twowayfeweights_result_other_treatment(res$df, treatments_rename, res$beta, random_weight_rename)
    df_result <- twowayfeweights_print_result_other_treatment(res, treatments_rename, D, res$beta, random_weight_rename)
  }
  
  if (length(path) != 0) {
    write.csv(df_result, path)
  }
  
}

})