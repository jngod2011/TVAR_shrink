#library(vars); data("Canada"); setwd("~/Documents/Uni/Fall_2017/TVAR_shrink/testing_R/");
#y <- diff(Canada); p = 1:3; type = "const"; exogen <- NULL; info_crit <- T
#trim <- 0.1; thVar <- 1; thDelay <- 1; nth <- 1; trim <- .1; thValue <- NULL
#ff <- est_var(diff(Canada), p = 2, type = "const", nth = 0)
# fff <- TVAR(diff(Canada), lag = 2, include = "const", nthresh = 1)
#p <- 1:3; thDelay <- 2; y <- diff(Canada); info_crit <- T; nth = 0
#thValue <- 0
est_var <- function(y, p = 1, type = c("const", "trend", "both", "none"), exogen = NULL, 
                    nth = 0, thDelay = 1, thVar = 1, trim = 0.1, thValue = NULL, 
                    info_crit = TRUE) {
  #y...data
  #p...lags (if info_crit == FALSE); possible lag lengths if (info_criteria == TRUE)
  #exogen...exogenous variables
  #nth...number thresholds: nth = 0 -> linear VAR; nth = 1 -> two regimes,... (ignored with given thValue)
  #thDelay...threshold delay
  #thVar...position or name of threshold variable
  #thValue...manually defined threshold value (grid search by default -> NULL by default); nth will be ignored
  y_orig <- as.matrix(y)
  nobs <- nrow(y_orig) #of complete sample!
  K <- ncol(y)
  if (K < 2) stop("Data is not multivariate!\n")
  colnames(y) <- make.names(colnames(y))
  type <- match.arg(type)
  if (info_crit & length(p) == 1) {
    lag_lengths <- 1:p 
  } else if (!(info_crit) & length(p) > 1) {
    lag_lengths <- p[1]
    warning("Info Criteria is equal to FALSE! Only the first lag length is evaluated!")
  } else {
    lag_lengths <- p
  }
  
  return_value <- list() #for all lag lengths
  #p <- lag_lengths[1]
  for (p in lag_lengths) { #all possible lag lengths
    return_value_lag <- list() #for particular lag length
    y <- y_orig[-c(1:p), ]
    sample_nobs <- nobs - p  #to be returned as sample size
    
    lag_y <- embed(y_orig, dimension = p + 1)[, -c(1:K)]
    names_rhs <- paste0("l", rep(1:p, each = K), "_", colnames(y)) 
    if (type == "none") {
      rhs <- lag_y
      colnames(rhs) <- names_rhs
    } else if (type == "const") {
      rhs <- cbind(1, lag_y)
      colnames(rhs) <- c("const", names_rhs)
    } else if (type == "trend") {
      rhs <- cbind((p+1):nobs, lag_y)
      colnames(rhs) <- c("trend", names_rhs)
    } else if (type == "both") {
      rhs <- cbind(1, (p+1):nobs, lag_y)
      colnames(rhs) <- c("const", "trend", names_rhs)
    }
    if (!(is.null(exogen))) {
      exogen <- exogen[-c(1:p), ]
      names_exogen <- paste("exog", colnames(exogen), sep = "_")
      colnames(exogen) <- names_exogen
      rhs <- cbind(rhs, exogen)
    }
    
    thresh_info <- list(nth = nth, thDelay = thDelay, thVar = thVar, trim = trim, thValue = thValue)
    var_info <- list(p = p, K = K, nobs = sample_nobs, type = type, rhs = rhs, y_orig = y_orig, 
                     info_crit = info_crit)
    
    if (nth != 0) {
      var_res <- est_threshold(y, thresh_info, var_info)
      var_res_add_info <- var_res
      var_res_add_info$reg_res <- NULL
      var_res <- var_res$reg_res
      return_value_lag <- list(var_res = var_res, var_res_add_info = var_res_add_info,
                           var_info = var_info, thresh_info = thresh_info, y = y)
    } else { #linear VAR
      var_res_add_info <- list()
      var_res <- lm.fit(x = rhs, y = y)
      cov_matrix <- t(var_res$residuals) %*% var_res$residuals / var_res$df.residual 
      coeff_by_regime <- t(var_res$coefficients) 
      stab_var <- stability_var(coeff_by_regime)
      var_res_add_info$coeff_array <- coeff_by_regime
      var_res_add_info$cov_matrix <- cov_matrix
      var_res_add_info$stability_results <- stab_var
      var_res_add_info$aic <- log(det(cov_matrix)) + 2*p*K^2/nrow(y) #Luetkephol 147ff
      var_res_add_info$bic <- log(det(cov_matrix)) + log(nrow(y))/nrow(y)*p*K^2
      return_value_lag <- list(var_res = var_res, var_res_add_info = var_res_add_info,
                           var_info = var_info, thresh_info = thresh_info, y = y)
    }
    return_value[[paste0("lag_", p)]] <- return_value_lag
  }
  
  if (!(is.null(info_crit))) {
    aic <- bic <- numeric()
    for (i in lag_lengths) {
      aic <- c(aic, return_value[[paste0("lag_", i)]]$var_res_add_info$aic)
      bic <- c(bic, return_value[[paste0("lag_", i)]]$var_res_add_info$bic)
    }
    return_value$info_choice <- matrix(c(which.min(aic), lag_lengths[which.min(aic)], 
                                         which.min(bic), lag_lengths[which.min(bic)]), nrow = 2, 
                                       dimnames = list(c("Model", "Lag_Length"), c("AIC", "BIC")))
  }
 
  return(return_value)
}

#lhs <- y; thresh_info <- thresh_info; var_info <- var_info
est_threshold <- function(lhs, thresh_info, var_info) {
  thVar <- thresh_info$thVar
  thDelay <- thresh_info$thDelay
  rhs <- var_info$rhs
  y_orig <- var_info$y_orig
  nth <- thresh_info$nth
  K <- var_info$K
  info_crit <- var_info$info_crit
  p <- var_info$p
  
  if(is.character(thVar)) { 
    thVar <- which(colnames(lhs) == thVar)
    if(length(thVar) == 0) stop("No endogenous variable called ", thresh_info$thVar, " found!\n")
  }
  if(thVar > ncol(y_orig)) 
    stop("Threshold Vaiable has to be smaller or equal to ", ncol(y_orig), ".")
  
  if (thDelay > p) {
    warning("thDelay Value is larger than number of lags of the model! Sample size decreases to ",
            var_info$nobs - thresh_info$thDelay + var_info$p, " observations")
    lhs <- lhs[-c(1:(thDelay - var_info$p)), ]
    rhs <- rhs[-c(1:(thDelay - var_info$p)), ]
  }
  
  if (is.null(thresh_info$thValue)) {#check next line again
    possible_thVals <- y_orig[(nrow(y_orig) - nrow(lhs) - thDelay + 1) : (nrow(y_orig) - thDelay), thVar]
    #ordered! -> left column always smaller than right in next step
    #n choose k; order doesn't matter; no repition
    useable_thVals <- t(combn(sort(unique(possible_thVals)) , nth)) 
    
    all_ssr <- numeric()
    estimation_thresh <- numeric() #saving all possible (estimateable) thresholds
    for (step_th in 1:nrow(useable_thVals)) { #loop over all possible thresholds from combn
      big_rhs <- matrix(0, nrow = nrow(lhs), ncol = (nth + 1) * ncol(rhs))
      names_big_rhs <- paste0(colnames(rhs), "_reg_", rep(1 : (nth+1), each = ncol(rhs)))
      colnames(big_rhs) <- names_big_rhs
      thresh_combination <- useable_thVals[step_th, ]
      
      #recording obs in each regime
      obs_in_regime <- matrix(0, nrow = nrow(lhs), ncol = nth+1)
      colnames(obs_in_regime) <- paste0("regime_", 1:(nth+1))
      for (k in 1:nrow(rhs)) { #finding suitable regime
        curr_thresh <- possible_thVals[k] #possible threshold contains all the thresh variables in right order
        which_regime <- sum(curr_thresh > thresh_combination) + 1
        start_pos <- ncol(rhs) * (which_regime-1) + 1 #start pos in big_rhs
        end_pos <- ncol(rhs) * (which_regime-1) + ncol(rhs)
        big_rhs[k, start_pos:end_pos] <- rhs[k, ]
        obs_in_regime[k, which_regime] <- 1
      }
      perc_in_regime <- colSums(obs_in_regime) / sum(obs_in_regime)
      trim_criteria <- all(perc_in_regime > thresh_info$trim) #has to be true to fulfill trim value requirement
      
      if (trim_criteria) {
        estimation_thresh <- rbind(estimation_thresh, thresh_combination) #saving all "estimateable thresh vals"
        reg_res <- lm.fit(x = big_rhs, y = lhs)
        ssr <- t(as.numeric(reg_res$residuals)) %*% as.numeric(reg_res$residuals)
        all_ssr <- c(all_ssr, ssr)
        #info criteria
        cov_all_resids <- t(reg_res$residuals) %*% reg_res$residuals / reg_res$df.residual #seems to work, but test again
        aic <- log(det(cov_all_resids)) + 2*p*K^2/nrow(lhs) #Luetkephol 147ff
        bic <- log(det(cov_all_resids)) + log(nrow(lhs))/nrow(lhs)*p*K^2
        
        if (!(exists("best_ssr"))) {#first attempted value?
          best_ssr <- ssr
          regime_mask <- obs_in_regime
          est_res <- list(reg_res = reg_res, best_thresh = thresh_combination, big_rhs = big_rhs,
                          perc_in_regime = perc_in_regime, ssr = best_ssr, aic = aic, bic = bic)
        } else if (ssr < best_ssr) {
          best_ssr <- ssr
          regime_mask <- obs_in_regime
          est_res <- list(reg_res = reg_res, best_thresh = thresh_combination, big_rhs = big_rhs,
                          perc_in_regime = perc_in_regime, ssr = best_ssr, aic = aic, bic = bic)
        }
      }
    }
  } else { #threshold value manually given
    nth <- length(thresh_info$thValue)  #ignore nth if thValue is given 
    big_rhs <- matrix(0, nrow = nrow(lhs), ncol = (nth + 1) * ncol(rhs))
    names_big_rhs <- paste0(colnames(rhs), "_reg_", rep(1 : (nth+1), each = ncol(rhs)))
    colnames(big_rhs) <- names_big_rhs
    
    #all thresh values in right order
    possible_thVals <- y_orig[(nrow(y_orig) - nrow(lhs) - thDelay + 1) : (nrow(y_orig) - thDelay), thVar]
    obs_in_regime <- matrix(0, nrow = nrow(lhs), ncol = nth + 1) #counting the obs in each regime
    for (k in 1:length(possible_thVals)) {
      which_regime <- sum(possible_thVals[k] > thresh_info$thVal) + 1
      start_pos <- ncol(rhs) * (which_regime-1) + 1
      end_pos <- ncol(rhs) * (which_regime-1) + ncol(rhs)
      big_rhs[k, start_pos:end_pos] <- rhs[k, ]
      obs_in_regime[k, which_regime] <- 1
    }
    perc_in_regime <- colSums(obs_in_regime) / sum(obs_in_regime)
    regime_mask <- obs_in_regime #for cov matrix and stability
    
    reg_res <- lm.fit(x = big_rhs, y = lhs)
    ssr <- t(as.numeric(reg_res$residuals)) %*% as.numeric(reg_res$residuals)
    cov_all_resids <- t(reg_res$residuals) %*% reg_res$residuals / (nrow(rhs) - ncol(big_rhs))
    aic <- log(det(cov_all_resids)) + 2*p*K^2/nrow(lhs) #Luetkephol 147ff
    bic <- log(det(cov_all_resids)) + log(nrow(lhs))/nrow(lhs)*p*K^2
    
    est_res <- list(reg_res = reg_res, best_thresh = NULL, big_rhs = big_rhs,
                    perc_in_regime = perc_in_regime, ssr = ssr, aic = aic, bic = bic)
  }
  
  #Cov matrix and stability
  resids <- est_res$reg_res$residuals
  cov_m_names <- list(colnames(lhs), colnames(lhs), paste0("regime_", 1:(nth+1)))
  cov_matrix <- array(0, dim = c(K, K, nth+1), dimnames = cov_m_names)
  coeff_by_regime <- coeff_2_regime_array(est_res) #coefficients in array by regime
  stability_results <- list()
  for (i in 1:(nth+1)) {
    obs_for_cov_matrix <- which(regime_mask[,i] == 1) #which obs belong to regime
    df_resid <- sum(regime_mask[,i]) - ncol(rhs) 
    curr_resid <- resids[obs_for_cov_matrix, ]
    cov_matrix[,,i] <- t(curr_resid) %*% curr_resid / df_resid 
    stability_results[[paste0("regime_", i)]] <- stability_var(coeff_by_regime[,,i])
  }
  
  if(is.null(thresh_info$thVal)) 
    add_info <- list(all_thresh = estimation_thresh, all_ssr = all_ssr)
  else 
    add_info <- NULL
  
  est_res$cov_matrix <- cov_matrix
  est_res$stability_results <- stability_results
  est_res$add_info <- add_info
  est_res$coeff_array <- coeff_by_regime
  est_res$regime_mask <- regime_mask
  return(est_res)
}

