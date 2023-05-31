
##### SVM ####

# Load required libraries
library(TMB)
setwd("/Users/thoma/Documents/skole/Master/Masteroppgave")

# Compile C++ functions
compile("SVM_ARp.cpp")
dyn.load(dynlib("SVM_ARp"))

compile("SVM_ARq.cpp")
dyn.load(dynlib("SVM_ARq"))

compile("SVM_ARMA.cpp")
dyn.load(dynlib("SVM_ARMA"))
var(exp())
# The makeSVM function creates an AD model object using the specified ARMA order and initial parameter values
makeSVM <- function(Y, arma_order = c(p = 1, q = 1), Mu = -6, Sigma = 0.1, avg = 0.00, burn_in = 100) {
  X <- rep(0.0, length(Y) + burn_in)
  df2 <- data.frame(Y = as.numeric(Y))
  
  # Check the specified ARMA order and create the appropriate model object
  if (arma_order["q"] == 0) {
    params <- list(X = X, PHI0 = rep(0.0, arma_order["p"]), Mu = Mu, logSigma = log(Sigma), avg = 0.00)
    model <- MakeADFun(df2, params, random = c("X"), DLL = "SVM_ARp", silent = F, ADreport = F)
  } else if (arma_order["p"] == 0) {
    params <- list(X = X, THETA0 = rep(0.0, arma_order["q"]), Mu = Mu, logSigma = log(Sigma), avg = 0.00)
    model <- MakeADFun(df2, params, random = c("X"), DLL = "SVM_ARq", silent = F, ADreport = F)
  } else {
    params <- list(X = X, PHI0 = rep(0.0, arma_order["p"]), THETA0 = rep(0.0, arma_order["q"]), Mu = Mu, logSigma = log(Sigma), avg = 0.00)
    model <- MakeADFun(df2, params, random = c("X"), DLL = "SVM_ARMA", silent = F, ADreport = F)
  }
  
  return(model)
}

# The mysvm function fits a stochastic volatility model to the observed time series using the specified ARMA order and initial parameter values
mysvm <- function(observed, arma_order = c(p = 1, q = 1), Mu = -6, Sigma = 0.1, avg = 0.00, burn_in = 100) {
  # Create the AD model object
  model <- makeSVM(observed, arma_order = arma_order, Mu = Mu, Sigma = Sigma, avg = avg, burn_in = burn_in)
  
  # Optimize the model parameters using the provided objective function and gradient
  fit <- nlminb(model$par, model$fn, model$gr)
  
  return(list(model = model, fit = fit))
}


# The adaptive_burn_in function computes the appropriate burn-in period based on the provided ARMA order and coefficients
adaptive_burn_in <- function(PHI, q) {
  # Calculate the minimum modulus of the roots of the AR polynomial
  minroots <- min(Mod(polyroot(c(1, -PHI))))
  
  # Check if the AR part of the model is stationary
  if (minroots <= 1) {
    stop("'ar' part of model is not stationary")
  }
  
  # Compute the adaptive burn-in period based on the minimum modulus of the roots and the MA order
  burn_in = length(PHI) + q + ifelse(length(PHI) > 0, ceiling(6 / log(minroots)), 0)
  return(burn_in)
}

# The mysvm_adaptive_burnin function fits a stochastic volatility model using an adaptive burn-in period
mysvm_adaptive_burnin <- function(y, arma_order = c(p = 1, q = 1), Mu = -6, Sigma = 0.1, avg = 0.00, initial_burn_in = 30) {
  # Set the initial burn-in period
  burn_in = initial_burn_in
  
  # Fit the SVM model with the initial burn-in period
  ret <- mysvm(y, arma_order = arma_order, Mu = Mu, Sigma = Sigma, avg = avg, burn_in = burn_in)
  
  # If the AR order is greater than 0, update the burn-in period adaptively
  if (arma_order["p"] > 0) {
    PHI <- get_parameters(ret$fit$par[1:p])
    burn_in <- adaptive_burn_in(PHI, arma_order["q"])
    if (burn_in > 1000) {
      burn_in = 1000
    }
    print(burn_in)
    # Refit the SVM model with the updated adaptive burn-in period
    ret <- mysvm(y, arma_order = arma_order, Mu = Mu, Sigma = Sigma, avg = avg, burn_in = burn_in)
  }
  
  # Store the final burn-in period in the results
  ret$burn_in = burn_in
  return(ret)
}

mysvm_p0_q0 <- function(y,par){
  lambda = par[1]
  sigma = exp(par[2]/2)
  n = length(y)
  epsilon <- y[-1] - y[-n] - lambda
  loglik <- sum(dnorm(epsilon, 0, sigma, log = TRUE))
  return(-loglik)
}

fit_mysvm_p0_q0 <- function(y, model = mysvm_p0_q0) {
  # Initial values for the optimization
  par_init <- c(0, -6)
  
  # Optimization using the optim function
  opt_result <- optim(
    par = par_init,
    fn = model,
    y = y,
    control = list(maxit = 10000),
    hessian = TRUE
  )
  
  # Return the optimization result
  return(opt_result)
}


# mapping from unconstrained to constrained parameters
pacf2param <- function(x) {
  n = length(x)
  r = x
  #print(r)
  y = matrix(0, nrow = n,ncol = n)
  for (k in 1:n) {
    i = 1
    while (i <= (k-1)) {
      y[k,i] <- y[k-1, i] - r[k]*y[k-1,k-i]
      i = i + 1
    }
    y[k,k] = r[k]
  }
  return(y)
}
get_parameters<- function(param_vector) {
  p = param_vector/sqrt(1 + param_vector^2)
  param = pacf2param(p)
  return(param[length(p),])
}



##### DATA #### 

#install.packages("quantmod")
library(quantmod)

symbols <- c("TSLA", "SPY", "AAPL", "MSFT", "AMZN", "GOOGL", "META", "BRK-B", "JPM", "JNJ")

end_date <- as.Date("2022-12-30") # Assuming the last trading day of 2022 is 2022-12-30
start_date <- end_date - 2000 # Download more data, as we will later filter down to the last 1000 trading days

prices_list <- lapply(symbols, function(symbol) {
  tryCatch({
    getSymbols(symbol, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
  }, error = function(e) {
    message(paste("Error downloading", symbol))
    return(NULL)
  })
})

# Extract the adjusted close prices
adj_prices <- lapply(prices_list, function(x) {
  if (!is.null(x)) {
    Ad(x)
  } else {
    NULL
  }
})

# Keep only the last 1000 trading days
adj_prices_last_1000 <- lapply(adj_prices, function(x) {
  if (!is.null(x)) {
    tail(x, 500)
  } else {
    NULL
  }
})

# Calculate the natural logarithm of the closing prices
log_adj_prices <- lapply(adj_prices_last_1000, function(x) {
  if (!is.null(x)) {
    log(x)
  } else {
    NULL
  }
})

log_adj_prices
combined_log_prices <- do.call(cbind, log_adj_prices)
colnames(combined_log_prices) <- symbols
combined_log_prices <- data.frame(date = index(log_adj_prices[[1]]), combined_log_prices)
combined_log_prices$date <- NULL
length(combined_log_prices)

extract_parameters_and_loglikelihood_1d <- function(OneDim_model_results) {
  model <- OneDim_model_results
  parameters <- model$fit$par
  # Extract phi and theta values if they exist, otherwise set to NA
  if ("PHI0" %in% names(parameters)) {
    phi <- get_parameters(parameters[names(parameters)=="PHI0"])
  } else {
    phi <- NA
  }
  if ("THETA0" %in% names(parameters)) {
    theta <- get_parameters(parameters[names(parameters)=="THETA0"])
  } else {
    theta <- NA
  }
  
  # Create phi and theta columns with NA values
  phi_columns <- setNames(as.list(rep(NA, length(phi))), paste0("phi", 1:length(phi)))
  theta_columns <- setNames(as.list(rep(NA, length(theta))), paste0("theta", 1:length(theta)))

  # Fill in the phi and theta values
  for (i in seq_along(phi)) phi_columns[[paste0("phi", i)]] <- phi[i]
  for (i in seq_along(theta)) theta_columns[[paste0("theta", i)]] <- theta[i]
  
  # Extract other model parameters and loglikelihood
  Mu <- parameters["Mu"]
  sigma <-exp(parameters["logSigma"])
  lambda <- parameters["avg"]
  loglikelihood <- -model$fit$objective
  
  # Create a data frame for the current model
  model_df <- data.frame(mu = Mu,
                         phi_columns,
                         theta_columns,
                         sigma = sigma,
                         lambda = lambda,
                         loglikelihood = loglikelihood)
  
  # Add the current model data frame to the result data frame

  return(model_df)
}

extract_parameters_and_loglikelihood <- function(model_results) {
  # Initialize an empty data frame to store the results
  result_df <- data.frame()
  
  # Loop through all models in model_results
  for (model_name in names(model_results)) {
    model <- model_results[[model_name]]
    parameters <- model$fit$par
    
    # Extract phi and theta values if they exist, otherwise set to NA
    phi <- if ("PHI0" %in% names(parameters)) get_parameters(parameters["PHI0"]) else NA
    theta <- if ("THETA0" %in% names(parameters)) get_parameters(parameters["THETA0"]) else NA
    
    # Create phi and theta columns with NA values
    phi_columns <- setNames(as.list(rep(NA, length(phi))), paste0("phi", 1:length(phi)))
    theta_columns <- setNames(as.list(rep(NA, length(theta))), paste0("theta", 1:length(theta)))
    
    # Fill in the phi and theta values
    for (i in seq_along(phi)) phi_columns[[paste0("phi", i)]] <- phi[i]
    for (i in seq_along(theta)) theta_columns[[paste0("theta", i)]] <- theta[i]
    
    # Extract other model parameters and loglikelihood
    Mu <- parameters["Mu"]
    sigma <- exp(parameters["logSigma"])
    lambda <- parameters["avg"]
    loglikelihood <- model$fit$objective
    
    # Create a data frame for the current model
    model_df <- data.frame(symbol = model_name,
                           phi_columns,
                           theta_columns,
                           mu = Mu,
                           sigma = sigma,
                           lambda = lambda,
                           loglikelihood = loglikelihood)
    
    # Add the current model data frame to the result data frame
    result_df <- rbind(result_df, model_df)
  }
  
  # Set the "symbol" column as the index and remove it
  row.names(result_df) <- result_df$symbol
  result_df$symbol <- NULL
  
  return(result_df)
}


svm_result_df <- extract_parameters_and_loglikelihood(model_results)
svm_result_df = parameter_list



##### GARCH #####


# Unconstrained log-likelihood function for the GARCH model
garch_loglik_unconst <- function(par, y, p = 1, q = 1) {
  # Transform unconstrained parameters to constrained parameters
  val2 <- garch_uncsonts_to_const(par, p = p, q = q)
  
  # Extract individual parameters
  alpha0 <- val2[1]
  alpha <-  val2[2:(1+p)]
  beta <- val2[(2+p):(1+p+q)]
  lambda <- val2[length(val2)]
  
  m <- max(p, q)
  n <- length(y)
  
  # Initialize sigma2 and epsilon
  sigma2 <- rep(alpha0, n)
  epsilon <- numeric(n)
  
  # Set initial values for epsilon
  epsilon[1] <- sqrt(alpha0) # / (1 - sum(alpha) - sum(beta)))
  epsilon[-1] <- y[-1] - y[-n] - lambda
  
  # Calculate sigma2 and the log-likelihood
  for (i in (m+1):n) {
    # Handle the case when p = 0 or q = 0
    alpha_term <- if (p > 0) sum(alpha * epsilon[(i-1):(i-p)]^2) else 0
    beta_term <- if (q > 0) sum(beta * sigma2[(i-1):(i-q)]) else 0
    
    sigma2[i] <- alpha0 + alpha_term + beta_term
  }
  
  # Compute the log-likelihood
  loglik <- sum(dnorm(epsilon[-(1:m)], 0, sqrt(sigma2[-(1:m)]), log = TRUE))
  
  return(-loglik)
}

simple_model <- function(par, y) {
  n <- length(y)
  mu <- par[1]
  lambda <- par[2]
  epsilon <- y[-1] - y[-n] - lambda
  loglik <- sum(dnorm(epsilon, exp(mu/2), log = TRUE))
  return(-loglik)
}


# Function to fit GARCH(p, q) model to a single time series
fit_garch <- function(y, p = 1, q = 1, model = garch_loglik_unconst) {
  # Initial values for the optimization
  par_init <- rep(0, p + q + 2)
  
  # Optimization using the optim function
  opt_result <- optim(
    par = par_init,
    fn = model,
    y = y,
    p = p,
    q = q,
    control = list(maxit = 10000),
    hessian = TRUE
  )
  
  # Return the optimization result
  return(opt_result)
}

garch_results <- lapply(combined_log_prices, fit_garch, p = 1, q = 2, model = garch_loglik_unconst)


extract_garch_parameters_and_loglikelihood_1d <- function(garch_results_1d, p=1, q=1) {
  model <- garch_results_1d
  parameters <- model$par
  
  # Map the unconstrained parameters to constrained parameters
  constrained_params <- garch_uncsonts_to_const(parameters, p = p, q = q)
  alpha0 <- constrained_params[1]
  
  if (p > 0) {
    alpha <- constrained_params[2:(1+p)]
  } else {
    alpha <- NA
  }
  
  if (q > 0) {
    beta <- constrained_params[(2+p):(1+p+q)]
  } else {
    beta <- NA
  }
  
  lambda <- constrained_params[length(constrained_params)]
  
  loglikelihood <- -model$value
  
  # Create alpha and beta columns with NA values
  alpha_columns <- setNames(as.list(rep(NA, length(alpha))), paste0("alpha", 1:length(alpha)))
  beta_columns <- setNames(as.list(rep(NA, length(beta))), paste0("beta", 1:length(beta)))
  
  # Fill in the alpha and beta values
  if (!is.null(alpha)) {
    for (i in seq_along(alpha)) {
      alpha_columns[[paste0("alpha", i)]] <- alpha[i]
    }
  }
  
  if (!is.null(beta)) {
    for (i in seq_along(beta)) {
      beta_columns[[paste0("beta", i)]] <- beta[i]
    }
  }
  
  # Create a data frame for the current model
  model_df <- data.frame(
                         alpha0 = alpha0,
                         alpha_columns,
                         beta_columns,
                         lambda = lambda,
                         loglikelihood = loglikelihood)
  
  # Add the current model data frame to the result data frame
  return(model_df)
}



# Fit GARCH(p, q) model to all time series in combined_log_prices
extract_garch_parameters_and_loglikelihood <- function(garch_results, p = 1, q = 1) {
  # Initialize an empty data frame to store the results
  result_df <- data.frame()
  
  # Loop through all models in garch_results
  for (model_name in names(garch_results)) {
    model <- garch_results[[model_name]]
    parameters <- model$par
    
    # Map the unconstrained parameters to constrained parameters
    constrained_params <- garch_uncsonts_to_const(parameters, p = p, q = q)
    alpha0 <- constrained_params[1]
    
    if (p > 0) {
      alpha <- constrained_params[2:(1+p)]
    } else {
      alpha <- NA
    }
    
    if (q > 0) {
      beta <- constrained_params[(2+p):(1+p+q)]
    } else {
      beta <- NA
    }
    
    lambda <- constrained_params[length(constrained_params)]
    
    loglikelihood <- -model$value
    
    # Create alpha and beta columns with NA values
    alpha_columns <- setNames(as.list(rep(NA, length(alpha))), paste0("alpha", 1:length(alpha)))
    beta_columns <- setNames(as.list(rep(NA, length(beta))), paste0("beta", 1:length(beta)))
    
    # Fill in the alpha and beta values
    if (!is.null(alpha)) {
      for (i in seq_along(alpha)) {
        alpha_columns[[paste0("alpha", i)]] <- alpha[i]
      }
    }
    
    if (!is.null(beta)) {
      for (i in seq_along(beta)) {
        beta_columns[[paste0("beta", i)]] <- beta[i]
      }
    }
    
    # Create a data frame for the current model
    model_df <- data.frame(symbol = model_name,
                           alpha0 = alpha0,
                           alpha_columns,
                           beta_columns,
                           lambda = lambda,
                           loglikelihood = loglikelihood)
    
    # Add the current model data frame to the result data frame
    result_df <- rbind(result_df, model_df)
  }
  # Set the symbol as the index
  rownames(result_df) = result_df$symbol
  result_df$symbol <- NULL
  return(result_df)
}



# Extract GARCH parameters and loglikelihood
garch_result_df <- extract_garch_parameters_and_loglikelihood(garch_results, p = 1, q = 2)
#garch_result_df2 <- extract_garch_parameters_and_loglikelihood(garch_results2)



# Initialize lists to store results for SVM and GARCH models
svm_results_list <- list()
garch_results_list <- list()

# Remove the date column from the combined_log_prices data frame
combined_log_prices$date <- NULL

# Create a progress bar

library(utils)
total_iterations <- (4*4) - 1
progress_bar <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration_count <- 0
svm_result = lapply(combined_log_prices, mysvm_adaptive_burnin, arma_order = c(p = 1, q = 1))
svm_result_df = extract_parameters_and_loglikelihood(svm_result)
garch_result <- lapply(combined_log_prices, fit_garch, p = 2, q = 1)
garch_result_df <- extract_garch_parameters_and_loglikelihood(garch_result, p = 2, q = 1)
cbind(-svm_result_df$loglikelihood, garch_result_df$loglikelihood)
dir.create("svm_results")
dir.create("garch_results")
##### grid search #####

# Loop through all combinations of p and q, except for p=0 and q=0
for (q in 0:3) {
  for (p in 0:3) {
    if (p == 0 && q == 0) next
    
    # Fit SVM models with adaptive burn-in
    start_time <- Sys.time()
    svm_result <- lapply(combined_log_prices, mysvm_adaptive_burnin, arma_order = c(p = p, q = q), Mu = 0, Sigma = 2, avg = 0.00)
    end_time <- Sys.time()
    svm_result_df <- extract_parameters_and_loglikelihood(svm_result)
    svm_result_df$time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    svm_results_list[[paste0("SVM_p_", p, "_q_", q)]] <- svm_result_df
    
    write.csv(svm_result_df, file = paste0("svm_results/SVM_p_", p, "_q_", q, ".csv"), row.names = TRUE)
    # Fit GARCH models
    start_time <- Sys.time()
    garch_result <- lapply(combined_log_prices, fit_garch, p = p, q = q)
    end_time <- Sys.time()
    garch_result_df <- extract_garch_parameters_and_loglikelihood(garch_result, p = p, q = q)
    garch_result_df$time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    garch_results_list[[paste0("GARCH_p_", p, "_q_", q)]] <- garch_result_df
    
    write.csv(garch_result_df, file = paste0("garch_results/GARCH_p_", p, "_q_", q, ".csv"), row.names = TRUE)
    # Update the progress bar
    iteration_count <- iteration_count + 1
    setTxtProgressBar(progress_bar, iteration_count)
  }
}

random_model <- lapply(combined_log_prices, fit_mysvm_p0_q0)
random_model$TSLA$value
AIC_00_model <- list(AIC = list(), lambda = list(), mu = list())
for(n in names(combined_log_prices)){
  AIC_00_model$AIC[[n]] <- 2*random_model[[n]]$value + 4
  AIC_00_model$lambda[[n]] <- random_model[[n]]$par[1]
  AIC_00_model$mu[[n]] <- random_model[[n]]$par[2]
}

a <- t(data.frame(AIC_00_model$AIC, row.names = "AIC"))
l <- t(data.frame(AIC_00_model$lambda, row.names = "lambda"))
m <- t(data.frame(AIC_00_model$mu, row.names = "mu"))
cbind(a,l,m)

# Close the progress bar
close(progress_bar)

# Combine the GARCH and SVM results lists into one list
combined_results_list <- c(garch_results_list, svm_results_list)

# Get the number of rows in the combined_log_prices dataframe
nrows <- ncol(combined_log_prices)

# Initialize an empty dataframe to store the loglikelihood values with the same number of rows as combined_log_prices
loglikelihood_df <- data.frame(matrix(ncol = 0, nrow = nrows))

# Loop through the combined results list
for (model_name in names(combined_results_list)) {
  # Extract the loglikelihood values for the current model
  loglikelihood_values <- combined_results_list[[model_name]]$loglikelihood
  
  # If the model is an SVM model, make the loglikelihood values positive
  if (grepl("SVM", model_name)) {
    loglikelihood_values <- -1 * loglikelihood_values
  }
  
  # Add the loglikelihood values to the dataframe with the appropriate model name
  loglikelihood_df[[model_name]] <- loglikelihood_values
}


# Initialize an empty dataframe to store the AIC values
svm_aic_df <- data.frame(matrix(ncol = 0, nrow = nrows))
garch_aic_df <- data.frame(matrix(ncol = 0, nrow = nrows))



# Counter for the row index in num_params_df
row_idx <- 1
p = 1
q = 1

# Loop through all combinations of p and q, except for p=0 and q=0
for (q in 0:3) {
  for (p in 0:3) {
    if (p == 0 && q == 0) next
    
    # Get the number of parameters for the current SVM and GARCH models
    svm_num_params <- 3 + p + q
    garch_num_params <- 2 + p + q
    
    # Get the loglikelihood values for the current SVM and GARCH models
    svm_loglikelihood <- loglikelihood_df[[paste0("SVM_p_", p, "_q_", q)]]
    garch_loglikelihood <- loglikelihood_df[[paste0("GARCH_p_", p, "_q_", q)]]
    
    # Calculate the AIC values for the current SVM and GARCH models
    svm_aic <- 2 * svm_num_params - 2 * svm_loglikelihood
    garch_aic <- 2 * garch_num_params - 2 * garch_loglikelihood
    
    # Add the AIC values to the aic_df dataframe with the appropriate model names
    svm_aic_df[[paste0("SVM_p_", p, "_q_", q, "_AIC")]] <- svm_aic
    garch_aic_df[[paste0("GARCH_p_", p, "_q_", q, "_AIC")]] <- garch_aic
    
    # Increment the row index for num_params_df
    row_idx <- row_idx + 1
  }
}

# Set the row names of the dataframe to match the original data
rownames(svm_aic_df) <- colnames(combined_log_prices)
rownames(garch_aic_df) <- colnames(combined_log_prices)
# Print the AIC dataframe
svm_aic_df <- t(svm_aic_df)
garch_aic_df <- t(garch_aic_df)

write.csv(svm_aic_df, "svm_aic.csv", row.names = TRUE)
write.csv(garch_aic_df, "garch_aic.csv", row.names = TRUE)

find_best_models_svm <- function(aic_values) {
  # Find the index of the minimum AIC value for each stock
  min_aic_indices <- apply(aic_values, 2, which.min)
  
  # Find the minimum AIC value for each stock
  min_aic_values <- apply(aic_values, 2, min)
  
  # Extract the best model names
  best_models <- rownames(aic_values)[min_aic_indices]
  
  # Extract p and q values from the model names
  p_values <- as.numeric(gsub("SVM_p_([0-9]+)_q_[0-9]+_AIC", "\\1", best_models))
  q_values <- as.numeric(gsub("SVM_p_[0-9]+_q_([0-9]+)_AIC", "\\1", best_models))
  
  # Combine the best models and their corresponding AIC scores into a data frame
  best_models_aic <- data.frame(p = p_values, q = q_values, AIC = min_aic_values)
  
  return(best_models_aic)
}
find_best_models_garch <- function(aic_values) {
  # Find the index of the minimum AIC value for each stock
  min_aic_indices <- apply(aic_values, 2, which.min)
  
  # Find the minimum AIC value for each stock
  min_aic_values <- apply(aic_values, 2, min)
  
  # Extract the best model names
  best_models <- rownames(aic_values)[min_aic_indices]
  
  # Extract p and q values from the model names
  p_values <- as.numeric(gsub("GARCH_p_([0-9]+)_q_[0-9]+_AIC", "\\1", best_models))
  q_values <- as.numeric(gsub("GARCH_p_[0-9]+_q_([0-9]+)_AIC", "\\1", best_models))
  
  # Combine the best models and their corresponding AIC scores into a data frame
  best_models_aic <- data.frame(p = p_values, q = q_values, AIC = min_aic_values)
  
  return(best_models_aic)
}

svm_best_models <- find_best_models_svm(svm_aic_df)
garch_best_models <- find_best_models_garch(garch_aic_df)
# Create a list of lists containing p and q values
model_list <- lapply(seq_along(p_values), function(i) list(p = p_values[i], q = q_values[i]))


# Assign the list of lists to the Model column in the data frame
best_models_aic$Model <- model_list





rolling_validation_splits <- function(combined_log_prices, train_window_size, test_window_size = 1) {
  train_splits <- list()
  test_splits <- list()
  
  # Loop through all the symbols in combined_log_prices
  for (symbol in names(combined_log_prices)) {
    data <- combined_log_prices[[symbol]]
    n <- length(data)
    total_iterations <- (n - train_window_size - test_window_size + 1)
    
    # Create lists to store the train and test splits for the current symbol
    symbol_train_splits <- list()
    symbol_test_splits <- list()
    
    # Perform the rolling cross-validation split for the current symbol
    for (i in seq(1, total_iterations, by = test_window_size)) {
      train_data <- data[1:(i + train_window_size - 1)]
      test_data <- data[(i + train_window_size):(i + train_window_size + test_window_size - 1)]
      
      symbol_train_splits[[length(symbol_train_splits) + 1]] <- train_data
      symbol_test_splits[[length(symbol_test_splits) + 1]] <- test_data
    }
    
    # Add the train and test splits for the current symbol to the final output lists
    train_splits[[symbol]] <- symbol_train_splits
    test_splits[[symbol]] <- symbol_test_splits
  }
  
  return(list(train = train_splits, val = test_splits))
}

cross_validation_sets <- rolling_validation_splits(combined_log_prices, train_window_size = 350, test_window_size = 25)



evaluate_models <- function(cross_validation_sets, svm_p_list, svm_q_list, garch_p_list, garch_q_list) {
  train_sets <- cross_validation_sets$train
  test_sets <- cross_validation_sets$val
  
  cond_ll_results <- list(svm = list(), garch = list())
  
  for (symbol in names(train_sets)) {
    train_symbol_splits <- train_sets[[symbol]]
    test_symbol_splits <- test_sets[[symbol]]
    svm_p <- svm_p_list[symbol]
    svm_q <- svm_q_list[symbol]
    garch_p <- garch_p_list[symbol]
    garch_q <- garch_q_list[symbol]
    names(svm_p) <- NULL
    names(svm_q) <- NULL
    names(garch_p) <- NULL
    names(garch_q) <- NULL
    svm_cond_ll_list <- c()
    garch_cond_ll_list <- c()
    
    for (i in seq_along(train_symbol_splits)) {
      train_data <- train_symbol_splits[[i]]
      test_data <- test_symbol_splits[[i]]
      print(length(train_data))
      # Fit the SVM model
      svm_result <- mysvm_adaptive_burnin(train_data, arma_order = c(p = svm_p, q = svm_q), Mu = 0, Sigma = 2, avg = 0.00)
      ll_train <- svm_result$fit$objective
      print(svm_result$burn_in)
      svm_test_train_model <- makeSVM(c(train_data, test_data),arma_order = c(p = svm_p, q = svm_q), burn_in = svm_result$burn_in)
      ll_test_train <- svm_test_train_model$fn(svm_result$fit$par)
      svm_cond_ll <- -ll_test_train + ll_train
      
      # Fit the GARCH model
      garch_result <- fit_garch(train_data, p = garch_p, q = garch_q)
      ll_train <- -garch_result$value
      ll_test_train <- -garch_loglik_unconst(garch_result$par, c(train_data, test_data), p = garch_p, q = garch_q)
      garch_cond_ll <- ll_test_train - ll_train
      
      svm_cond_ll_list <- c(svm_cond_ll_list, svm_cond_ll)
      garch_cond_ll_list <- c(garch_cond_ll_list, garch_cond_ll)
      print(svm_cond_ll)
      print(garch_cond_ll)
    }
    
    cond_ll_results$svm[[symbol]] <- svm_cond_ll_list
    cond_ll_results$garch[[symbol]] <- garch_cond_ll_list
  }
  
  return(cond_ll_results)
}

svm_p <- svm_best_models$p
svm_q <- svm_best_models$q
garch_p <- garch_best_models$p
garch_q <- garch_best_models$q
test <- svm_p[names(cross_validation_sets$train)[1]]
names(test) <- NULL

names(svm_p) <- names(combined_log_prices)
names(svm_q) <- names(combined_log_prices)
names(garch_p) <- names(combined_log_prices)
names(garch_q) <- names(combined_log_prices)

ret <- evaluate_models(cross_validation_sets, svm_p, svm_q, garch_p, garch_q)

ret

##### Fit Best Models ####

Final_models <- list(svm = list(), garch = list())
Filal_models_ret <- list(svm = list(), garch = list())
combined_log_prices$date <- NULL
names(combined_log_prices)[1]
#i = "TSLA"
for (i in names(combined_log_prices)) {
  
  
  y <- combined_log_prices[[i]]
  pe <- svm_p[i]
  qe <- svm_q[i]
  names(pe) <- NULL
  names(qe) <- NULL
  
  start_time <- Sys.time()
  print(pe)
  svm_result <- mysvm_adaptive_burnin(y, arma_order = c(p = pe, q = qe), Mu = -7, Sigma = 0.05, avg = 0.00)
  print(svm_result$fit$par)
  Filal_models_ret$svm[[i]] <- svm_result
  end_time <- Sys.time()
  result_df <- extract_parameters_and_loglikelihood_1d(svm_result)
  result_df$time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  print(result_df)
  Final_models$svm[[i]] <- result_df
  
  p <- garch_p[i]
  q <- garch_q[i]
  names(p) <- NULL
  names(q) <- NULL
  
  start_time <- Sys.time()
  garch_result <- fit_garch(y, p, q)
  end_time <- Sys.time()
  result_df <- extract_garch_parameters_and_loglikelihood_1d(garch_result, p, q)
  result_df$time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  Final_models$garch[[i]] <- result_df
  Filal_models_ret$garch[[i]] <-garch_result
}

test <- data.frame()
for (i in (names(Final_models$svm))) {
  test <- rbind(test, Final_models$svm[[i]])
}


svm_result$fit$objective
Final_models$garch

deltavar
svm_result$model
##### find STD of MLE #####


### STD GARCH ###
library(emdbook)
y
p <-
optimised_garch <- fit_garch(y, )

p <- garch_p[i]
q <- garch_q[i]
names(p) <- NULL
names(q) <- NULL

garch_result <- fit_garch(y, p, q)


garch_hess <- garch_result$hessian
is 
ML_estimates <- garch_uncsonts_to_const_parameters(garch_result$par, p, q)
fisher_info <- solve(garch_hess)
var_delta <- deltavar(garch_uncsonts_to_const_parameters,
                      meanval = garch_result$par,
                      Sigma = fisher_info)



# Function to map unconstrained parameters to constrained parameters
garch_uncsonts_to_const_parameters()


garch_uncsonts_to_const_parameters <- function(par, p = 1, q = 1) {
  if (q == 0){
    alpha0 <- exp(par[1]) 
    alpha <- exp(par[2:(1+p)])/(1 + sum(exp(par[2:(1 + p + q)])))
    lambda <- par[length(par)]
    return(c(alpha0, alpha, lambda))
  }
  else if (p == 0) {
    alpha0 <- exp(par[1]) 
    beta <- exp(par[(2):(1+q)])/(1 + sum(exp(par[2:(1 + p + q)])))
    lambda <- par[length(par)]
    return(c(alpha0, beta, lambda))
  }
  else{
    alpha0 <- exp(par[1]) 
    alpha <- exp(par[2:(1+p)])/(1 + sum(exp(par[2:(1 + p + q)])))
    beta <- exp(par[(2+p):(1+p+q)])/(1 + sum(exp(par[2:(1 + p + q)])))
    lambda <- par[length(par)]
    return(c(alpha0, alpha, beta,lambda))
  }
}

garch_uncsonts_to_const_parameters(garch_result$par, p, q)



#install.packages("numDeriv")
library(numDeriv)

garch_transform_wrapper <- function(par) {
  p <- p # Define p
  q <- q # Define q
  garch_uncsonts_to_const_parameters(par, p, q)
}

#install.packages("pracma")
library(pracma)

jacobian_matrix <- function(f, x) {
  n <- length(x)
  J <- jacobian(f, x)
  return(J)
}



params <- garch_result$par

unconst_params <- garch_uncsonts_to_const_parameters(params, p = 3, q = 2)

jac_mat <- jacobian_matrix(garch_transform_wrapper, params)
fisher_mat <- solve(garch_hess)

asym_cov_matrix <- jac_mat %*% solve(fisher_mat) %*% t(jac_mat)


constrained_sd <- sqrt(diag(asym_cov_matrix))

result <- data.frame(Estimate = unconst_params, Std.Error = constrained_sd)


###### GARCH results_list_sd #####

# rename garch rows
rename_garch_rows <- function(df,p, q) {
  n <- rownames(df)
  if(p==0) {
    n[1] <- paste0("alpha", 0)
    n[(2+p):(1+p+q)] <- paste0("beta", 1:q)
    n[2+p+q] <- "lambda"
  }else if(q == 0) {
    n[1] <- paste0("alpha", 0)
    n[(2):(1+p)] <- paste0("alpha", 1:p)
    n[2+p+q] <- "lambda"
  }else{
    n[1] <- paste0("alpha", 0)
    n[(2):(1+p)] <- paste0("alpha", 1:p)
    n[(2+p):(1+p+q)] <- paste0("beta", 1:q)
    n[2+p+q] <- "lambda"
  }
  rownames(df) <- n
  return(df)
}

garch_results_list_sd <- list()

# loop through and find the standard deviation of the parameter estimates for all garch models
for (i in names(combined_log_prices)) {
  
  y <- combined_log_prices[[i]]
  
  p <- garch_p[i]
  q <- garch_q[i]
  names(p) <- NULL
  names(q) <- NULL
  garch_result <- Filal_models_ret$garch[[i]]
  garch_hess <- garch_result$hessian
  params <- garch_result$par
  
  unconst_params <- garch_uncsonts_to_const_parameters(params, p = p, q = q)

  
  jac_mat <- jacobian_matrix(garch_transform_wrapper, params)
  fisher_mat <- solve(garch_hess)
  print(jac_mat)
  print(fisher_mat)
  asym_cov_matrix <- jac_mat %*% fisher_mat %*% t(jac_mat)
  
  
  constrained_sd <- sqrt(diag(asym_cov_matrix))
  
  result <- data.frame(Std.Error = constrained_sd)
  result <- rename_garch_rows(result, p, q)
  garch_results_list_sd[[i]] <- result
}

garch_results_list_sd
###### SVM results_list_sd #####

# function that maps unconstrained to constrained parameters
svm_unconstrained_to_constrained <- function(par, p, q) {
  if(p == 0) {
    theta <- get_parameters(par[1:q])
    Mu <- par[(p+q+1)]
    sigma <- exp(par[p+q+2])
    lambda <- par[p+q+3]
    return(c(theta, Mu, sigma, lambda))
  } else if(q == 0) {
    phi <- get_parameters(par[(1):(p)])
    Mu <- par[(p+q+1)]
    sigma <- exp(par[p+q+2])
    lambda <- par[p+q+3]
    return(c(phi, Mu, sigma, lambda))
  }
  else{
    phi <- get_parameters(par[(1):(p)])
    theta <- get_parameters(par[(1+p):(q+p)])
    Mu <- par[(p+q+1)]
    sigma <- exp(par[p+q+2])
    lambda <- par[p+q+3]
    return(c(phi, theta ,Mu, sigma, lambda))
  }
}

# Transform wrapper
svm_transform_wrapper <- function(par) {
  p <- p # Define p
  q <- q # Define q
  svm_unconstrained_to_constrained(par, p, q)
}

# Renames all rows in the dataframe for a SVM
rename_svm_rows <- function(df,p, q) {
  n <- rownames(df)
  if(p==0) {
    n[1:q] <- paste0("theta", 1:q)
    n[p+q+1] <- "mu"
    n[p+q+2] <- "sigma"
    n[p+q+3] <- "lambda"
  }else if(q == 0) {
    n[1:p] <- paste0("phi", 1:p)
    n[p+q+1] <- "mu"
    n[p+q+2] <- "sigma"
    n[p+q+3] <- "lambda"
  }else{
    n[1:p] <- paste0("phi", 1:p)
    n[(1+p):(p+q)] <- paste0("theta", 1:q)
    n[p+q+1] <- "mu"
    n[p+q+2] <- "sigma"
    n[p+q+3] <- "lambda"
  }
  rownames(df) <- n
  return(df)
}



svm_results_list_sd <- list()

# loop through and find the standard deviation of the parameter estimates for all garch models
for(i in names(combined_log_prices)) {
  
  y <- combined_log_prices[[i]]
  
  p <- svm_p[i]
  q <- svm_q[i]
  names(p) <- NULL
  names(q) <- NULL
  
  svm_result <- Filal_models_ret$svm[[i]]
  r <- sdreport(svm_result$model)
  #summary(r, "report")
  covariance_mat <- r$cov.fixed
  params <- svm_result$fit$par
  
  const_params <- svm_unconstrained_to_constrained(params, p = p, q = q)

  jac_mat <- jacobian_matrix(svm_transform_wrapper, params)
  
  asym_cov_matrix <- jac_mat %*% covariance_mat %*% t(jac_mat)
  
  
  constrained_sd <- sqrt(diag(asym_cov_matrix))
  
  result <- data.frame(Std.div = constrained_sd)
  result <- rename_svm_rows(result, p, q)
  svm_results_list_sd[[i]] <- result
}








#### SVM Theoretical vs Empirical ACF####

#install.packages("ltsa")
library(ltsa)


empirical_acf <- function(lambda, y) {
  
  epsilon <- y[-1] - y[-length(y)] - lambda
  
  epsilon2 <- epsilon^2
  
  emp_acvf <- acf(epsilon2, type = "covariance", plot = F)
  emp_acf <- acf(epsilon2, plot = F)
  return(list(acvf = emp_acvf, acf = emp_acf))
}


Filal_models_ret$svm[""]
test <- Filal_models_ret$svm[["TSLA"]]
test
param <- test$fit$par
lambda <- param[length(param)]# for svm
y <- combined_log_prices[["TSLA"]]
#acf(diff(combined_log_prices[["AMZN"]])^2)
empirical_and_theoretical_ACF_ACVF <- function(symbol, plot_ACF = F, plot_ACVF = F) {
  y <- combined_log_prices[[symbol]]
  pe <- svm_p[symbol]
  qe <- svm_q[symbol]
  model <- Filal_models_ret$svm[[symbol]]
  param <- model$fit$par
  param <- svm_unconstrained_to_constrained(param, pe, qe)
  lambda <- param[pe + qe + 3]
  emp_acf <- empirical_acf(lambda, y)
  
  if(pe == 0 ) {
    theta <- param[(1):(pe+qe)]
    sigma2 <- param[(pe+qe+2)]^2
    emp_acvf_arma <- tacvfARMA(theta  =theta, sigma2 = sigma2, maxLag = 15)
  } else if (qe == 0) {
    phi <- param[1:pe]
    sigma2 <- param[(pe+qe+2)]^2
    emp_acvf_arma <- tacvfARMA(phi = phi, sigma2 = sigma2, maxLag = 15)
  } else{
    phi <- param[1:pe]
    theta <- param[(pe+1):(pe+qe)]
    sigma2 <- param[(pe+qe+2)]^2
    emp_acvf_arma <- tacvfARMA(phi = phi, theta = theta, sigma2 = sigma2, maxLag = 15)
  }
  mu <- param[(pe+qe+1)]
  svm_emp_acvf <-exp(2*mu + 2*emp_acvf_arma[1])*(exp(emp_acvf_arma) - 1)
  svm_emp_acf <- svm_emp_acvf/svm_emp_acvf[1]
  
  lags <- 0:(length(svm_emp_acvf) - 1)
  significance_level_acf <- qnorm((1 + 0.95)/2)/sqrt(length(y)-1)
  significance_level_acvf <- significance_level_acf*emp_acf$acvf$acf[1]
  if (plot_ACF == T) {
    # Create a vector of lags
    tit <- paste0("SVM(", pe, ",", qe, "), ACF Plot, ", symbol)
    y_lim <- c(min(c(emp_acf$acvf$acf[1:length(svm_emp_acf)], svm_emp_acvf)),1)
    # Plot the first set of ACVF values
    plot(lags, svm_emp_acf, type = "h", xlab = "Lag", ylab = "ACF", main = tit,ylim = y_lim, cex.main = 0.7)
    abline(h = significance_level_acf,lty = 2)
    abline(h = 0)
    # Add the second set of ACVF values
    lines(lags, emp_acf$acf$acf[1:length(svm_emp_acf)], col = "red", type = "p")
    
    # Add a legend
    legend("topright", legend = c("Theoretic ACF", "empirical ACF"), col = c("black", "red"), lty = 1, cex = 0.5)
  }
  if (plot_ACVF == T) {
    y_lim <- c(min(c(emp_acf$acvf$acf[1:length(svm_emp_acf)], svm_emp_acvf)),max(c(emp_acf$acvf$acf[1:length(svm_emp_acf)]), svm_emp_acvf))
    # Plot the first set of ACVF values
    tit <- paste0("SVM(", pe, ",", qe, "), ACVF Plot, ", symbol)
    plot(lags, svm_emp_acvf, type = "h", xlab = "Lag", ylab = "ACVF", main = tit, ylim=y_lim, cex.main = 0.7)
    
    abline(h = significance_level_acvf, lty = 2)
    abline(h = 0)
    # Add the second set of ACVF values
    lines(lags, emp_acf$acvf$acf[1:length(svm_emp_acf)], col = "red", type = "p")
    
    # Add a legend
    legend("topright", legend = c("Theoretic ACVF", "empirical ACVF"), col = c("black", "red"), lty = 1, cex = 0.5)
  }
  theo = list(ACF = svm_emp_acf, ACVF = svm_emp_acvf, ARMA_ACVF = emp_acvf_arma)
  emp = list(ACF = emp_acf$acf$acf[1:length(svm_emp_acf)], ACVF = emp_acf$acf$acf[1:length(svm_emp_acf)])
  ret = list(theoretical = theo, empirical = emp)
  return(ret)
}


emp_theo_list = list()
par(mfrow = c(2, 1))
for (i in names(combined_log_prices)) {
  emp_theo_list[[i]] = empirical_and_theoretical_ACF_ACVF(i, plot_ACF = T, plot_ACVF = T)
}
empirical_and_theoretical_ACF_ACVF("SPY", plot = F)

Filal_models_ret$svm[["META"]]$fit$par

#### GARCH Theoretical vs Empirical ACF####
qnorm(0.975)

garch_empirical_and_theoretical_ACF_ACVF <- function(symbol, plot_ACF = F, plot_ACVF = F) {
  y <- combined_log_prices[[symbol]]
  pe <- garch_p[symbol]
  qe <- garch_q[symbol]
  model <- Filal_models_ret$garch[[symbol]]
  param <- model$par
  param <- garch_uncsonts_to_const_parameters(param, pe, qe)
  lambda <- param[pe + qe + 2]
  print(lambda)
  emp_acf <- empirical_acf(lambda, y)
  print(emp_acf)
  
  n[1] <- paste0("alpha", 0)
  n[(2):(1+p)] <- paste0("alpha", 1:p)
  n[(2+p):(1+p+q)] <- paste0("beta", 1:q)
  n[2+p+q] <- "lambda"
  
  
  
  if(pe == 0) {
    alpha0 <- param[1]
    beta <- param[(2):(1+qe)]
    theta = -beta
    phi = beta
    emp_acvf_arma <- ARMAacf(ar = phi,ma = theta, lag.max = 15)
  } else if (qe == 0) {
    alpha0 <- param[1]
    alpha <- param[(2):(1+pe)]
    phi = alpha
    emp_acvf_arma <- ARMAacf(ar = phi, lag.max = 15)
  } else{
    alpha0 <- param[1]
    alpha <- param[(2):(1+pe)]
    beta <- param[(2+pe):(1+pe+qe)]
    m = max(pe, qe)
    alph = numeric(m)
    alph[1:pe] <- alpha
    bet = numeric(m)
    bet[1:qe] = beta
    phi = alph + bet
    theta = -beta
    emp_acvf_arma <- ARMAacf(ar = phi, ma = theta, lag.max = 15) #HVA ER SIMGA2?
  }
  
  lags <- 0:(length(emp_acvf_arma) - 1)
  significance_level_acf <- qnorm((1 + 0.95)/2)/sqrt(length(y)-1)
  #significance_level_acvf <- significance_level_acf*emp_acf$acvf$acf[1]
  if (plot_ACF == T) {
    # Create a vector of lags
    tit <- paste0("GARCH(", pe, ",", qe, "), ACF Plot, ", symbol)
    y_lim <- c(min(c(emp_acf$acvf$acf[1:length(emp_acvf_arma)], emp_acvf_arma)),1)
    # Plot the first set of ACVF values
    plot(lags, emp_acvf_arma, type = "h", xlab = "Lag", ylab = "ACF", main = tit,ylim = y_lim, cex.main = 0.7)
    abline(h = significance_level_acf,lty = 2)
    abline(h = 0)
    # Add the second set of ACVF values
    lines(lags, emp_acf$acf$acf[1:length(emp_acvf_arma)], col = "red", type = "p")
    
    # Add a legend
    legend("topright", legend = c("Theoretic ACF", "empirical ACF"), col = c("black", "red"), lty = 1, cex = 0.5)
  }
  theo = list(ACF = emp_acvf_arma) #, ACVF = svm_emp_acvf)
  emp = list(ACF = emp_acf$acf$acf[1:length(emp_acvf_arma)] )#, ACVF = emp_acf$acf$acf[1:length(svm_emp_acf)])
  ret = list(theoretical = theo, empirical = emp)
  return(ret)
}


emp_theo_list = list(svm = list(), garch = list())
par(mfrow = c(2, 1))
for (i in names(combined_log_prices)) {
  emp_theo_list$svm[[i]] = empirical_and_theoretical_ACF_ACVF(i, plot_ACF = T)
  emp_theo_list$garch[[i]] = garch_empirical_and_theoretical_ACF_ACVF(i, plot_ACF = T)
}
gamma_0 <- list()
for (i in names(combined_log_prices)) {
  #gamma_0[[i]] <- emp_theo_list$svm[[i]]$theoretical$ACVF[1]
  #emp_theo_list$garch[[i]] = garch_empirical_and_theoretical_ACF_ACVF(i, plot_ACF = T)
  gamma_0[[i]] <- emp_theo_list[[i]]$theoretical$ARMA_ACVF[1]
}

marginal_variance_SVM <- list()
marginal_variance_GARCH <- list()
empirical_varaince <- list()
for (i in names(combined_log_prices)){
  marginal_variance_SVM[[i]] <- exp(Final_models$svm[[i]]$mu + gamma_0[[i]]/2)
  
  a0 <- Final_models$garch[[i]]$alpha0
  l <- length(Final_models$garch[[i]])
  params <- Final_models$garch[[i]][-c((l-2):l)]
  params <- params[-1]
  params <- params[!is.na(params)]
  marginal_variance_GARCH[[i]] <- a0/(1-sum(params))
  v <- diff(combined_log_prices[[i]]) - mean(diff(combined_log_prices[[i]]))
  empirical_varaince[[i]] <- var(v)
}

var(combined_log_prices[[i]])
df = rbind(data.frame(marginal_variance_SVM), data.frame(marginal_variance_GARCH), data.frame(empirical_varaince))
rownames(df) <- c("SVM", "GARCH", "Empirical")
df <- t(df)


