#setwd("/home/alex/Documents/Uni/Fall_2017/TVAR_shrink/frf_TVAR/")
#source("../testing_R/est_var.R"); library(vars)
#data("Canada"); est_res <- est_var(diff(Canada), p = 2, type = "const", nth = 1)
#impulse <- response <- NULL; horizon <- 20; ortho <- TRUE; comulative <- FALSE
#boot <- TRUE; ci <- .95; runs <- 100; seed <- NULL;
ortho_imp <- est_irf(est_res)
#est_res <- est_var(diff(Canada), p = 2, type = "const", nth = 0)
est_irf <- function(est_res, impulse = NULL, response = NULL, horizon = 20, ortho = TRUE, cumulative = FALSE,
                    boot = TRUE, ci = .95, runs = 100, seed = NULL) {
  K <- est_res$var_info$K
  p <- est_res$var_info$p
  coeff_array <- est_res$var_res_add_info$coeff_array
  cov_matrix <- est_res$var_res_add_info$cov_matrix
  variable_names <- colnames(est_res$var_res$coefficients)
  if(length(dim(coeff_array)) == 2) { #linear VAR
    coeff_array_new <- array(0, dim = c(dim(coeff_array), 1), dimnames = list(rownames(coeff_array),
                                                                              colnames(coeff_array), 1))
    cov_matrix_new <- array(0, dim = c(dim(cov_matrix), 1), dimnames = list(rownames(cov_matrix),
                                                                           colnames(cov_matrix), 1))
    coeff_array_new[,, 1] <- coeff_array
    coeff_array <- coeff_array_new
    cov_matrix_new[,, 1] <- cov_matrix
    cov_matrix <- cov_matrix_new
  } 
  no_regime <- dim(coeff_array)[3]
  pos_lagged_vars <- which(substr(colnames(coeff_array[,,1]), 1,1) == "l")
  
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
  
  return_value <- list()
  for (reg in 1:no_regime) { #responses for each regime
    #put A matrices in array
    name_for_regime <- paste0("regime_", reg)
    curr_coeff <- coeff_array[,,reg]
    coeff_matrix <- array(0, dim = c(K, K, p))
    for (i in 1:p) { 
      row_pos <- pos_lagged_vars[((i-1)*K+1) : ((i-1)*K+K)] 
      coeff_matrix[,,i] <- est_res$var_res$coefficients[row_pos, 1:K]
    }
    
    #MA terms
    phi <- array(0, dim = c(K, K, horizon+1))
    phi[,,1] <- diag(K)
    for (i in 2:(horizon+1)) { #calculating the MA terms (phi; Luetkepohl p. 23)
      which_A <- 1
      for (j in (i-1):(i-p)) {
        if (j > 0) 
          phi[,,i] <- phi[,,i] + phi[,,j] %*% coeff_matrix[,,which_A]
        which_A <- which_A + 1
      }
    }
    
    if (ortho) { #orthogonalized impulses (theta; Luetkepohl p. 59)
      P <- t(chol(cov_matrix[,, reg]))
      theta <- array(0, dim = dim(phi))
      for (i in 1:(horizon+1)) 
        theta[,,i] <- phi[,,i] %*% P
    }
    
    #responses finally...
    return_reg <- list()
    #imp <- impulse_vecs[1]; rsp <- response_vecs[1]
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
          return_reg[[name_list]] <- curr
        } else {
          for (i in 1:(horizon+1)) 
            curr[i] <- rsp %*% phi[,,i] %*% imp
          return_reg[[name_list]] <- curr
        }
      }
    }
    if (no_regime > 1) {
      return_value[[name_for_regime]] <- return_reg 
    } else{
      return_value <- return_reg
    }
  }
  return(return_value)
}
