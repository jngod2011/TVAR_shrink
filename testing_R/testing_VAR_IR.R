library(vars)
library(tsDyn)
library(ggplot2)
library(ggthemes)
library(dplyr)

rm(list = ls())

theme_abc <- theme_hc() +
  theme(axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) + 
  theme(plot.title = element_text(size=8, face="bold")) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))+
  theme(legend.position="bottom", legend.text=element_text(size=8)) 

setwd("/home/alex/Documents/Uni/Fall_2017/TVAR_shrink/frf_TVAR/")

source("../testing_R/est_var.R")
source("../testing_R/aux_functions.R")
source("../testing_R/est_girf_v2.R")
data("Canada")
can <- est_var(diff(Canada), p = 2, type = "const", nth = 1, trim = 0.2)
can.girf <- est_girf(can)
plot(can.girf$girf$regime_1[,4], type = "l")


can1 <- est_var(diff(Canada), p = 1, type = "const", nth = 1, thValue = 0)
can2 <- est_var(diff(Canada), p = 1, type = "const", nth = 1, thVar = "prod")
imp <- est_girf(can)

dim(imp$raw_girf)
reg1_ir <- imp$raw_girf[,4,1,,]
dim(reg1_ir)
df_all <- data.frame(horizon = numeric(), obs = numeric())
for (i in 1:dim(reg1_ir)[1]) {
  curr <- as.numeric(reg1_ir[i,,])
  df <- data.frame(horizon = i-1, obs = curr)
  df_all <- rbind(df_all, df)
}

df_all$horizon <- as.factor(df_all$horizon)
ggplot(df_all, aes(horizon, obs)) + 
  geom_boxplot(fill = alpha("blue", .1), outlier.colour = "orange") +
  theme_abc



resids_by_reg <- get_regime_residuals(ff)
resids_by_reg$regime_1

P1 <- t(chol(ff$var_res_add_info$cov_matrix[,,1]))
P2 <- t(chol(ff$var_res_add_info$cov_matrix[,,2]))
P1_inv <- solve(P1)
P2_inv <- solve(P2)

struc_err1 <- t(P1_inv %*% t(resids_by_reg$regime_1))
struc_err2 <- t(P2_inv %*% t(resids_by_reg$regime_2))

cov1 <- t(struc_err1) %*% struc_err1 / (sum(ff$var_res_add_info$regime_mask[,1]) - ncol(ff$var_info$rhs))
cov2 <- t(struc_err2) %*% struc_err2 / (sum(ff$var_res_add_info$regime_mask[,2]) - ncol(ff$var_info$rhs))
eye <- diag(4)
colnames(eye) <- colnames(cov1); rownames(eye) <- rownames(cov2)
all.equal(cov1, eye) #works!
all.equal(cov2, eye) #works!

##############################################################
###############################################################

#frf Data loaded
frf_data <- read.table("./frf-transformed.txt", sep = "\t", header = TRUE)
frf_est_data <- frf_data[117:224, c("delta_g", "delta_gdp", "debt_gdp", "d_inflation", "delta_i", "SPREAD")]

est_frf <- est_var(y = frf_est_data, p = 3, type = "const", nth = 1, thDelay = 2, thVar = 6, info_crit = T)
source("../testing_R/est_girf_v3.R")
imp <- est_girf(est_res = est_frf, horizon = 19, shocked_variable = "delta_g", sample_hist = T, model_choice = 1,
                seed = 1)
plot(imp$girf[1,,2], type = "l"); abline(h = 0, col = "red")
#imp_ortho <- est_irf(est_frf)
#plot(imp_ortho$regime_1$delta_g_TO_delta_g, type= "l")
est_can <- est_var(diff(Canada), p = 3, nth = 1)
imp <- est_girf(est_can, model_choice = 2)
plot(imp$girf[1,,1], type = "l")
abline(h = 0, col = "red")

plot(imp$girf$regime_2[,6], type = "l")
abline(h = 0, col = "red")

frf_est_data2 <- frf_data[117:224, c("delta_g", "delta_gdp", "debt_gdp", "d_inflation", "delta_i", "MIX")]
est_frf <- est_var(y = frf_est_data, p = 1, type = "const", nth = 1, thDelay = 2, thVar = 6)
est_frf$var_res_add_info$best_thresh
imp <- est_girf(est_res = est_frf, horizon = 19, shocked_variable = "delta_g", shock_size = 1, 
                n.history = 200, n.boot = 500, seed = 1)

plot(imp$girf$regime_2[,1], type = "l")
abline(h = 0, col = "red")

dim(imp$raw_girf)
reg1_ir <- imp$raw_girf[,2,1,,]
dim(reg1_ir)
df_all <- data.frame(horizon = numeric(), obs = numeric())
for (i in 1:dim(reg1_ir)[1]) {
  curr <- as.numeric(reg1_ir[i,,])
  df <- data.frame(horizon = i-1, obs = curr)
  df_all <- rbind(df_all, df)
}

df_all$horizon <- as.factor(df_all$horizon)
ggplot(df_all, aes(horizon, obs)) + 
  geom_boxplot(fill = alpha("blue", .1), outlier.colour = "orange") +
  theme_abc

te <- VAR(frf_est_data, ic = "AIC")
cov_m <- summary(te)$covres
P <- t(chol(cov_m))
P_inv <- solve(P)
struc_err <- t(P_inv %*% t(residuals(te)))
t(struc_err) %*% struc_err / (te$obs - 7) #works out; diagonal = 1



#VAR example for ordering: (because FRF seems to imply the opposite)
#ordered 1st: immediate effect on 1,2,3,4
#ordered 2nd: immediate effect on 2,3,4 ...
data("Canada")
Canada <- Canada[, c("prod", "e", "U", "rw")]
p1ct <- VAR(Canada, p = 2, type = "const")
ff <- est_var(diff(Canada), p = 2, type = "const", nth = 1)

regime_obs <- which(ff$var_res$regime_mask[,2] == 1)
regime_resid <- ff$var_res$reg_res$residuals[regime_obs, ]
P <- t(chol(ff$var_res$cov_matrix[,,2]))
P_inv <- solve(P)
struc_err <- t(P_inv %*% t(regime_resid))
t(struc_err) %*% struc_err / (sum(ff$var_res$regime_mask[,2]) - ncol(ff$var_res$coeff_array[,,2]))


summ_p1ct <- summary(p1ct)
p1ct.irf <- irf(p1ct, response = "e", n.ahead = 48, boot = TRUE)
#plot(p1ct.irf)
#View(vars:::irf.varest) #impulse response code

tvar_res <- TVAR(diff(Canada), include = "const", model = "TAR", lag = 2)
tvar_res$coefficients


#EXOG VARIABLES HAVE TO START WITH EXOG!!!!

tvar_res <- TVAR(diff(Canada), include = "const", model = "TAR", lag = 2, trim = .1)
#Own VAR
#y <- diff(Canada); p = 2; type = "const"; exogen <- NULL
#trim <- 0.1; thVal <- NULL
#y = frf_est_data; p = 1; type = "const"; nth = 1; thDelay = 2; thVar = 6; exogen <- NULL
est_var <- function(y, p = 1, type = c("const", "trend", "both", "none"), exogen = NULL, 
                    nth = 0, thDelay = 1, thVar = 1, trim = 0.1, thVal = NULL) {
  #y...data
  #p...lags
  #exogen...exogenous variables
  #nth...number thresholds: nth = 0 -> linear VAR; nth = 1 -> two regimes,...
  #thDelay...threshold delay
  #thVar...position or name of threshold variable
  #thVal...manually defined threshold value (grid search by default -> NULL by default)
  
  y_orig <- as.matrix(y)
  y <- y_orig[-c(1:p), ]
  nobs <- nrow(y_orig)
  K <- ncol(y)
  sample_nobs <- nobs - p  #to be returned as sample size
  if (K < 2) stop("Data is not multivariate!\n")
  colnames(y) <- make.names(colnames(y))
  type <- match.arg(type)
  
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
  
  thresh_info <- list(nth = nth, thDelay = thDelay, thVar = thVar, trim = trim, thVal = thVal)
  var_info <- list(p = p, K = K, nobs = sample_nobs, type = type, rhs = rhs, y_orig = y_orig)
  ### all the threshold BS here
  if (nth != 0) {
    thresh_res <- est_threshold(y, thresh_info, var_info)
    return_value <- list(var_res = thresh_res, y = y, var_info = var_info, thresh_info = thresh_info)
  } else { #linear VAR
    var_res <- lm.fit(x = rhs, y = y)
    cov_matrix <- t(var_res$residuals) %*% var_res$residuals / var_res$df.residual 
    stab_var <- stability_var(var_res)
    return_value <- list(var_res = var_res, cov_matrix = cov_matrix, roots_A = stab_var, y = y,
                         var_info = var_info, thresh_info = thresh_info)
  }
  return(return_value)
}

#lhs <- y; rhs <- rhs; thresh_info <- thresh_info
est_threshold <- function(lhs, thresh_info, var_info) {
  thVar <- thresh_info$thVar
  thDelay <- thresh_info$thDelay
  rhs <- var_info$rhs
  y_orig <- var_info$y_orig
  if(is.character(thVar)) {
    thVar <- which(colnames(y) == thVar)
    if(length(thVar) == 0) stop("No endogenous variable called ", thresh_info$thVar, " found!\n")
  }
  
  if (thresh_info$thDelay > var_info$p) {
    warning("thDelay Value is larger than number of lags of the model! Sample size decreases to ",
            var_info$nobs - thresh_info$thDelay + var_info$p, " observations")
    lhs <- lhs[-c(1:(thresh_info$thDelay - var_info$p)), ]
    rhs <- rhs[-c(1:(thresh_info$thDelay - var_info$p)), ]
  }
  
  if (is.null(thresh_info$thVal)) {#check next line again
    possible_thVals <- var_info$y_orig[(nrow(var_info$y_orig) - nrow(lhs) - thDelay + 1) : 
                                         (nrow(y_orig) - thDelay), thVar]
    trimming_values <- quantile(possible_thVals, c(thresh_info$trim, 1 - thresh_info$trim))
    useable_thVals <- possible_thVals[possible_thVals >= trimming_values[1] & 
                                        possible_thVals <= trimming_values[2]]
    
    best_ssr <- all_ssr <- numeric()
    for (curr_thresh in useable_thVals) {
      observ_reg <- which(possible_thVals <= curr_thresh) # <= as in Galvao (2003)
      low_lhs <- lhs[observ_reg, ]
      low_rhs <- rhs[observ_reg, ]
      reg_low <- lm.fit(x = low_rhs, y = low_lhs)
      ssr_low <- t(as.numeric(reg_low$residuals)) %*% as.numeric(reg_low$residuals)
      
      observ_reg <- which(possible_thVals > curr_thresh) 
      up_lhs <- lhs[observ_reg, ]
      up_rhs <- rhs[observ_reg, ]
      reg_up <- lm.fit(x = up_rhs, y = up_lhs)
      ssr_up <- t(as.numeric(reg_up$residuals)) %*% as.numeric(reg_up$residuals)
      
      all_ssr <- c(all_ssr, ssr_up + ssr_low)
      
      if (length(best_ssr) == 0) {#first attempted value?
        best_ssr <- ssr_low + ssr_up
        est_res <- list(reg_low = reg_low, reg_up = reg_up, best_thresh = curr_thresh,
                        low_lhs = low_lhs, low_rhs = low_rhs, up_lhs = up_lhs, up_rhs = up_rhs)
      } else if ((ssr_low + ssr_up) < best_ssr) {
        best_ssr <- ssr_low + ssr_up
        est_res <- list(reg_low = reg_low, reg_up = reg_up, best_thresh = curr_thresh,
                        low_lhs = low_lhs, low_rhs = low_rhs, up_lhs = up_lhs, up_rhs = up_rhs)
      }
    }
  } else { #threshold value manually given
    observ_reg <- which(possible_thVals <= thresh_info$thVal) 
    low_lhs <- lhs[observ_reg, ]
    low_rhs <- rhs[observ_reg, ]
    reg_low <- lm.fit(x = low_rhs, y = low_lhs)
    
    observ_reg <- which(possible_thVals > curr_thresh) 
    up_lhs <- lhs[observ_reg, ]
    up_rhs <- rhs[observ_reg, ]
    reg_up <- lm.fit(x = up_rhs, y = up_lhs)
    
    est_res <- list(reg_low = reg_low, reg_up = reg_up, best_thresh = NULL,
                    low_lhs = low_lhs, low_rhs = low_rhs, up_lhs = up_lhs, up_rhs = up_rhs)
  }
  
  cov_matrix_low <- t(est_res$reg_low$residuals) %*% est_res$reg_low$residuals / est_res$reg_low$df.residual
  cov_matrix_up <- t(est_res$reg_up$residuals) %*% est_res$reg_up$residuals / est_res$reg_up$df.residual
  stab_var_low <- stability_var(est_res$reg_low)
  stab_var_up <- stability_var(est_res$reg_up)
  
  est_res[["cov_matrix_low"]] <- cov_matrix_low
  est_res[["cov_matrix_up"]] <- cov_matrix_up
  est_res[["stab_var_low"]] <- stab_var_low
  est_res[["stab_var_up"]] <- stab_var_up
  if(is.null(thresh_info$thVal)) 
    add_info <- list(possible_thVals = possible_thVals, useable_thVals = useable_thVals, all_ssr = all_ssr)
  else 
    add_info <- NULL
  
  est_res[["add_info"]] <- add_info
  
  return(est_res)
  
}

#PLOTS
plt_data <- data.frame(Threshold_Number = 1:length(useable_thVals), Threshold_Value = useable_thVals,
                       SSR = all_ssr)
min_ssr <- min(all_ssr)
plt_data$min_ssr <- "0"
plt_data$min_ssr[plt_data$SSR == min_ssr] <- "1"
ggplot(plt_data, aes(Threshold_Number, Threshold_Value)) + geom_line() + 
  xlab("") + ylab("Threshold Value") + geom_hline(yintercept = est_res$best_thresh, col = "red") +
  theme_abc

ggplot(plt_data, aes(Threshold_Value, SSR, col = min_ssr)) + geom_point() + xlab("") + 
  theme_abc + theme(legend.position="none") + scale_colour_manual(values = c("darkblue", "orange")) 

plt_data1 <- data.frame(New_order = 1:length(possible_thVals), Threshold_Value = sort(possible_thVals))
ggplot(plt_data1, aes(New_order, Threshold_Value)) + geom_line() + 
  xlab("") + ylab("Threshold Value") + geom_hline(yintercept = est_res$best_thresh, col = "red") + 
  geom_hline(yintercept = trimming_values, col = "green") + 
  theme_abc

# plt_series <- y
# ind_low <- which(plt_series %in% est_res$low_lhs) 
# ind_low <- ind_low[ind_low <= nrow(plt_series)]
# plt_series <- as.data.frame(plt_series)
# plt_series$regime <- "up"
# plt_series$regime[ind_low] <- "low"
# plt_series$date <- 1:nrow(plt_series)
# 
# ggplot(plt_series, aes(date, plt_series[,1])) + geom_line() + 
#   geom_vline(xintercept = regime)



#est_res <- var_res
stability_var <- function(est_res) { #est_res...from lm.fit
  #stability analysis
  coeff_matrix <- est_res$coefficients
  res_lag <- rownames(coeff_matrix)[(substr(rownames(coeff_matrix),1,1) == "l")]
  pos_lag <- which((substr(rownames(coeff_matrix),1,1) == "l"))
  min_pos <- min(pos_lag)
  max_lag <- max(as.numeric(substr(res_lag, 2, 2)))
  K <- ncol(coeff_matrix)
  p <- length(pos_lag) / K
  
  A <- numeric() #constructing the big VAR(1) matrix
  for (i in 1:p) {
    curr_matrix <- coeff_matrix[((i-1)*K + min_pos) : ((i-1)*K + min_pos + (K-1)), ] 
    A <- cbind(A, curr_matrix)
  }
  if(p > 1) { #in case of VAR(1) we are already done
    for (i in 1:(p-1)) {
      curr_matrix <- matrix(0, nrow = K, ncol = ncol(A))
      for (j in 1:K)  
        curr_matrix[j, (i-1)*K + j] <- 1
      A <- rbind(A, curr_matrix)
    }
  }

  root_vals <- Mod(eigen(A)$values)
  return(root_vals)
}



Canada_diff <- diff(Canada)

p1ct <- VAR(Canada_diff, p = 3, type = "const")
summary(p1ct)$roots
irf_p1ct <- irf(p1ct, impulse = "prod", response = "prod", n.ahead = 20)
plot(irf_p1ct)
est_res <- est_var(diff(Canada), p = 2, type = "const", nth = 1)
est_res$roots_A

#dec <- t(chol(cov_matrix))
#dec %*% t(dec) - cov_matrix

#assuming intercept only for the moment
horizon <- 20; shocked_variable = 1; n.history = 100; n.boot = 100; shock_size <- 1; sim_interval = c(.05, .95)
est_girf <- function(est_res, horizon = 20, shocked_variable = 1, shock_size = 1, n.history = 100, n.boot = 100,
                     seed.hist = NULL, seed.boot = NULL, sim_interval = c(.05, .95)) {
  cl <- match.call()
  p <- est_res$var_info$p
  K <- est_res$var_info$K
  var_res <- est_res$var_res
  redform_err <- var_res$residuals
  y_orig <- est_res$var_info$y_orig
  
  pos_lagged_vars <- which(substr(rownames(var_res$coefficients), 1,1) == "l")
  pos_fixed_vars <- which(rownames(var_res$coefficients) %in% c("const", "trend"))
  pos_exog_vars <- which(substr(rownames(var_res$coefficients), 1, 4) == "exog")
  
  coeff_fixed <- t(var_res$coefficients[pos_fixed_vars, , drop = F])
  fixed_cols <- matrix(rep(1, ncol(coeff_fixed)))
    
  coeff_matrix <- array(0, dim = c(K, K, p))
  for (i in 1:p) {
    row_pos <- pos_lagged_vars[((i-1)*K + 1) : ((i-1)*K + K)]
    coeff_matrix[,,i] <- t(var_res$coefficients[row_pos, 1:K]) #we have to transpose here bc of lm.fit format!
  }
  
  P <- t(chol(est_res$cov_matrix))
  P_inv <- solve(P)
  str_err <- t(P_inv %*% t(redform_err)) 
  
  girf <- matrix(0, nrow = horizon + 1, ncol = K)
  girf_i <- array(0, dim = c(horizon + 1, K, n.history))
  for (i in 1:n.history) {
    hist_pos <- sample(p:nrow(y_orig), size = 1)
    start_hist <- y_orig[hist_pos:(hist_pos - p + 1), ] #has to be done differently in TVAR case
    girf_ij <- array(0, dim = c(horizon + 1, K, n.boot))
    for (j in 1:n.boot) {
      pos_boot_err <- sample(1:nrow(str_err), size = horizon + 1, replace = TRUE)
      non_shock_err <- str_err[pos_boot_err, ]
      shock_err <- non_shock_err
      shock_err[1, shocked_variable] <- shock_err[1, shocked_variable] + shock_size
      redform_non_shock <- t(P %*% t(non_shock_err))
      redform_shock <- t(P %*% t(shock_err))
      
      sys_non_shock <- sys_shock <- matrix(0, nrow = horizon + 1, ncol = K)
      curr_hist_non_shock <- curr_hist_shock <- matrix(start_hist, ncol = K) #some problems with one lag otherwise
      for (curr_step in 1:(horizon+1)) {
        for (k in 1:p) {
          sys_non_shock[curr_step, ] <- sys_non_shock[curr_step, ] + 
            coeff_matrix[,,k] %*% t(curr_hist_non_shock[k,, drop = F])
          sys_shock[curr_step, ] <- sys_shock[curr_step, ] +
            coeff_matrix[,,k] %*% t(curr_hist_shock[k,, drop = F])
        }
        sys_non_shock[curr_step, ] <- sys_non_shock[curr_step, ] + coeff_fixed %*% fixed_cols +
           t(redform_non_shock[curr_step,, drop = F])
        sys_shock[curr_step, ] <- sys_shock[curr_step] + coeff_fixed %*% fixed_cols + 
           t(redform_shock[curr_step,, drop = F])
  
        curr_hist_non_shock <- rbind(sys_non_shock[curr_step, ], curr_hist_non_shock[-p, ])
        curr_hist_shock <- rbind(sys_shock[curr_step, ], curr_hist_shock[-p, ])
        
      }
      girf_ij[,,j] <- sys_shock - sys_non_shock
    }
    girf_i[,,i] <- apply(girf_ij, MARGIN = c(1,2), mean)
  }
  
  girf <- apply(girf_i, MARGIN = c(1,2), mean)
  low_girf <- apply(girf_i, MARGIN = c(1,2), quantile, sim_interval[1])
  high_girf <- apply(girf_i, MARGIN = c(1,2), quantile, sim_interval[2])
  ret_list <- list(girf = girf, low_girf = low_girf, high_girf = high_girf, Call = cl)
  
  return(ret_list)
}

p1ct <- VAR(diff(Canada), p = 4, type = "const")
serial.test(p1ct)
#sum_p1ct <- summary(p1ct)

est_res <- est_var(diff(Canada), p = 2, type = "const")
# p_inv <- solve(t(chol(est_res$cov_matrix)))
# struc <- t(p_inv %*% t(est_res$var_res$residuals))
# cov(struc)
girf_res <- est_girf(est_res, shocked_variable = 1, n.history = 100, n.boot = 100)
irf_res <- est_irf(est_res)
#girf_res$Call
plot(girf_res$girf[,4], type = "l")
lines(irf_res$e_TO_rw, col = "green")
abline(h = 0, col ="red")
plot(irf_res$prod_TO_prod, type="l")

sum_p1ct$covres  - est_res$cov_matrix
pp <- t(chol(sum_p1ct$covres))
all.equal(pp %*% t(pp), sum_p1ct$covres)
structural_error <- t(solve(pp) %*% t(residuals(p1ct)))
cov(structural_error)

plot(imp_resp[,1], type = "l")




p1ct <- VAR(diff(Canada), p = 2, type = "const")
summ_p1ct <- summary(p1ct)
p1ct.irf <- irf(p1ct, impulse = "prod", response = "e", n.ahead = 20, boot = TRUE, ortho = T)
#a <- p1ct.irf$irf$prod[2]
plot(p1ct.irf)
a <- p1ct.irf$irf$prod 

#est_res <- est_var(diff(Canada), p = 2, type = "const")
#impulse <- "prod"; response <- "e"; horizon <- 20; ortho <- T; cumulative <- FALSE
#boot <- TRUE; ci <- .95; runs <- 100; seed <- NULL
est_irf <- function(est_res, impulse = NULL, response = NULL, horizon = 20, ortho = TRUE, cumulative = FALSE,
                    boot = TRUE, ci = .95, runs = 100, seed = NULL) {
  K <- est_res$var_info$K
  p <- est_res$var_info$p
  variable_names <- colnames(est_res$var_res$coefficients)
  pos_lagged_vars <- which(substr(rownames(est_res$var_res$coefficients), 1,1) == "l")
  
  #impulse vectors
  impulse_vecs <- list()
  if (is.character(impulse)) { 
    count <- 1
    for (i in 1:length(impulse)) {
      pos_imp <- which(variable_names == impulse[i])
      if (length(pos_imp) == 0) stop("Variable: ", impulse[i], " does not exist!")
      curr_imp <- vector("numeric", length = K)
      curr_imp[pos_imp] <- 1
      impulse_vecs[[count]] <- curr_imp
      count <- count + 1
    }
  } else if (is.numeric(impulse)) {
    count <- 1
    for (i in impulse) { 
      if (i > K) stop("There are only: ", K, " variables in the VAR!")
      curr_imp <- vector("numeric", length = K)
      curr_imp[i] <- 1
      impulse_vecs[[count]] <- curr_imp
      count <- count + 1
    }
  } else if (is.list(impulse)) {
    count <- 1
    for (i in impulse) { 
      i <- unlist(i)
      if (is.numeric(i)) {
        if(length(i) != K) stop("Vector length has to be equal to ", K, " !")
        impulse_vecs[[count]] <- i
      } else { #character case
        pos_imp <- which(variable_names %in% i)
        if (length(pos_imp) != length(i)) stop("Variables: ", i, " do not exist!")
        curr_imp <- vector("numeric", length = K)
        curr_imp[pos_imp] <- 1
        impulse_vecs[[count]] <- curr_imp
      }
      count <- count + 1
    }
  } else if (is.null(impulse)) {#all impulses
    impulse <- diag(K)
    for (i in 1:K)
      impulse_vecs[[i]] <- impulse[,i] 
  }
  
  #response vectors
  response_vecs <- list()
  if (is.numeric(response)) {
    count <- 1
    for (i in response) {
      if (i > K) stop("There are only: ", K, " variables in the VAR!")
      curr_imp <- vector("numeric", length = K)
      curr_imp[i] <- 1
      response_vecs[[count]] <- curr_imp
      count <- count + 1
    }
  } else if (is.character(response)) {
    count <- 1
    for (i in response) {
      pos_resp <- which(variable_names == i)
      if (length(pos_resp) == 0) stop("Variable: ", impulse[i], " does not exist!")
      curr_resp <- vector("numeric", length = K)
      curr_resp[pos_resp] <- 1
      response_vecs[[count]] <- curr_resp
      count <- count + 1
    }
  } else if (is.null(response)) {
    responses <- diag(K)
    for (i in 1:K)
      response_vecs[[i]] <- responses[,i] 
  }
  
  #put A matrices in array
  coeff_matrix <- array(0, dim = c(K, K, p))
  for (i in 1:p) { 
    row_pos <- pos_lagged_vars[((i-1)*K+1) : ((i-1)*K+K)] 
    coeff_matrix[,,i] <- t(est_res$var_res$coefficients[row_pos, 1:K]) #we have to transpose here bc of lm.fit format!
  }
  
  #MA terms
  phi <- array(0, dim = c(K, K, horizon+1))
  phi[,,1] <- diag(4)
  for (i in 2:(horizon+1)) { #calculating the MA terms (phi; Luetkepohl p. 23)
    which_A <- 1
    for (j in (i-1):(i-p)) {
      if (j > 0) 
        phi[,,i] <- phi[,,i] + phi[,,j] %*% coeff_matrix[,,which_A]
      which_A <- which_A + 1
    }
  }
  
  if (ortho) { #orthogonalized impulses (theta; Luetkepohl p. 59)
    P <- t(chol(est_res$cov_matrix))
    theta <- array(0, dim = dim(phi))
    for (i in 1:(horizon+1)) 
      theta[,,i] <- phi[,,i] %*% P
  }
  
  #responses finally...
  return_value <- list()
  for (imp in impulse_vecs) { 
    for (rsp in response_vecs) {
      imp <- unlist(imp)
      rsp <- unlist(rsp)
      name_list <- paste(variable_names[which(imp == 1)],
                          variable_names[which(rsp == 1)],
                          sep = "_TO_")
      curr <- vector("numeric" , length = horizon+1)
      if (ortho) {
        for (i in 1:(horizon+1)) 
          curr[i] <- rsp %*% theta[,,i] %*% imp
        return_value[[name_list]] <- curr
      } else {
        for (i in 1:(horizon+1)) 
          curr[i] <- rsp %*% phi[,,i] %*% imp
        return_value[[name_list]] <- curr
      }
    }
  }
  return(return_value)
}



