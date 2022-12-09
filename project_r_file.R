# Project
# install.packages("latex2exp")
library(latex2exp)
# install.packages("ggiplot")
library(styler)
library(formatR)
library(TMB)
library(ggplot2)
require(gridExtra)
library(TTR)
library(quantmod)

##### Generate SV-process #####
stochastic_volatility_series_AR1 <- function(n = 1000, theta = 0, phi = 0.5, mu = 0, sigma = 1, y0 = 0) {
  m <- 100 # fade in period
  v <- rnorm(m + n, mean = 0, sd = 1)
  # print(v)
  u <- rnorm(m + n, mean = 0, sd = sigma)
  X <- numeric(m + n)
  Y <- numeric(m + n)
  X[1] <- mu # assume initial log volatility is mu
  Y[1] <- y0 # assume initial value of the series is y0
  for (i in 2:(m + n)) {
    X[i] <- mu + phi * (X[i - 1] - mu) + u[i]
    Y[i] <- Y[i - 1] + theta + exp(X[i] / 2) * v[i]
  }
  # print(length(x))
  return(list(x = X[(m + 1):(m + n)], y = Y[(m + 1):(m + n)]))
}



# generate and plot a instance of the series
set.seed(1)
simulated_1 <- stochastic_volatility_series_AR1(1000, y0 = 100, theta = 0.1, mu = 1, sigma = 0.1, phi = 0.95)

df <- data.frame(simulated_1)

# function for plotting observations
plot_obs <- function(y) {
  ggplot(data = data.frame(y = y), aes(x = 1:length(y), y = y)) +
    geom_line() +
    labs(title = "Observations", x = "time step, t", y = "Observed value, y")
}

# function for plotting latent variables (AR(1))
plot_latent <- function(x) {
  ggplot(data = data.frame(x = x), aes(x = 1:length(x), y = x)) +
    geom_line() +
    labs(title = "Latent variable", x = "time step, t", y = "latent variable, x")
}

# Plotting observed values and AR(1) process in the same figure
plot_both <- function(y, x) {
  p1_obs <- plot_obs(y)
  p1_latent <- plot_latent(x)
  grid.arrange(p1_obs, p1_latent, ncol = 1)
}

# Figure 4
plot_both(simulated_1$y, simulated_1$x)

set.seed(2)
# Big sigma
simulated_2 <- stochastic_volatility_series_AR1(n = 1000, theta = 0.1, phi = 0.95, mu = 1, sigma = 1)
p2 <- plot_obs(simulated_2$y)

set.seed(2)
# Small sigma
simulated_3 <- stochastic_volatility_series_AR1(n = 1000, theta = 0.1, phi = 0.95, mu = 1, sigma = .5)
p3 <- plot_obs(simulated_3$y)

# Figure 5
grid.arrange(p2, p3, ncol = 1)
# make a new model with other model parameters



# Fit the parretameters using TMB

# Compile C++ file and load the model
compile("SVM_lastlast.cpp")
dyn.load(dynlib("SVM_lastlast"))

# Makes the SV model for observations
makeSVM <- function(observed, initial = list(Mu = 0, Theta = 0, Phi = 0.0, Sigma = 1)) {
  X <- rep(0.0, length(observed))
  df2 <- data.frame(Y = observed)
  params <- list(X = X, Mu = initial$Mu, Theta = initial$Theta, logitPhi = 2 / (1 + exp(-initial$Phi)) - 1, logSigma = log(initial$Sigma))
  model <- MakeADFun(df2, params, random = c("X"), DLL = "SVM_lastlast", silent = T, ADreport = FALSE)
  return(model)
}

# Estimates the MLE´s for the SVM
mysvm <- function(observed, initial = list(Mu = 0, Theta = 0, Phi = 0.0, Sigma = 1)) {
  model <- makeSVM(observed, initial)

  # fit the model
  fit <- nlminb(model$par, model$fn, model$gr)
  # summary(fit, "report")
  return(list(model = model, fit = fit))
}



# Testing sequences of different length
test_n <- function(B = 100, N = c(15, 100, 1000, 3000, 5000, 10000)) {
  # Parameter values
  theta <- 0.1
  phi <- -0.95
  mu <- 1
  sigma <- 0.1

  # loop through the parameters
  r_mat <- array(dim = c(4, 2, length(N), B))
  # different choices of n
  for (n in 1:length(N)) {
    # parametric bootstrap
    for (b in 1:B) {
      ret <- stochastic_volatility_series_AR1(n = N[n], theta = theta, phi = phi, mu = mu, sigma = sigma)
      y <- ret$y
      ret2 <- mysvm(y)
      sd <- sdreport(ret2$model)
      sum2 <- summary(sd, "report")
      sum1 <- summary(sd)
      # sum1[c(1,4),]
      r_mat[, , n, b] <- rbind(sum1[c(1, 4), ], sum2)
    }
  }
  return(r_mat)
}

# run model
varying_n <- test_n()
# takes in a matrix and returns a data frame used for plotting

# function format input values for plotting
fix_input <- function(mat) {
  new_m <- apply(mat, 1, quantile, probs = c(0.05, 0.5, 0.95))
  df <- data.frame(t(new_m))
  df$index <- seq(1, nrow(df))
  return(df)
}

# plot
plot_vars <- function(data, true_value, x, xlab = "index", ylab = "mle") {
  ggplot(data, aes(y = X50., x = x)) +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = X5., ymax = X95.), color = "blue", alpha = 0.5) +
    geom_hline(aes(yintercept = true_value), color = "red") +
    xlab(xlab) +
    ylab(ylab)
}

# format input values for plotting in function plot_param_estimation_difflength_obs
fix_input_plot_all <- function(mat) {
  df_mu <- fix_input(mat[1, 1, , ])
  df_theta <- fix_input(mat[2, 1, , ])
  df_sigma <- fix_input(mat[3, 1, , ])
  df_phi <- fix_input(mat[4, 1, , ])
  return(list(df_mu = df_mu, df_theta = df_theta, df_sigma = df_sigma, df_phi = df_phi))
}


# plotting parameter estimates for varying N
plot_param_estimation_difflength_obs <- function(mat) {
  l <- fix_input_plot_all(mat)
  N <- c(15, 100, 1000, 3000, 5000, 10000)
  p4 <- plot_vars(l$df_mu, 1, log(N), xlab = TeX("$\\log(n)$"), ylab = TeX("$\\hat{\\mu}$"))
  p5 <- plot_vars(l$df_theta, .1, log(N), xlab = TeX("$\\log(n)$"), ylab = TeX("$\\hat{\\theta}$"))
  p6 <- plot_vars(l$df_sigma, 0.1, log(N), xlab = TeX("$\\log(n)$"), ylab = TeX("$\\hat{\\sigma}$"))
  p7 <- plot_vars(l$df_phi, 0.95, log(N), xlab = TeX("$\\log(n)$"), ylab = TeX("$\\hat{\\phi}$"))
  grid.arrange(p4, p5, p6, p7, ncol = 2)
}

# Figure 7
plot_param_estimation_difflength_obs(varying_n)


# simulate and fit parameter to different instances of theta
test_theta <- function(B = 100, theta = c(-0.1, 0, 0.1)) {
  # Parameter choices
  N <- 1000 # length of ts
  phi <- 0.95
  mu <- 1
  sigma <- 0.1
  # do a grid search and save parameters
  r_mat <- array(dim = c(4, 2, length(theta), B))
  # loop through the parameters
  for (n in 1:length(theta)) {
    # parametric bootstrap
    for (b in 1:B) {
      ret <- stochastic_volatility_series_AR1(n = N, theta = theta[n], phi = phi, mu = mu, sigma = sigma)
      y <- ret$y
      ret2 <- mysvm(y)
      sd <- sdreport(ret2$model)
      sum2 <- summary(sd, "report")
      sum1 <- summary(sd)
      # sum1[c(1,4),]
      r_mat[, , n, b] <- rbind(sum1[c(1, 4), ], sum2)
    }
  }
  return(r_mat)
}
# returning matrix of parameter estimates
varying_theta <- test_theta()

# Plotting different theta values with confint
plot_theta_estimates <- function(mat_theta) {
  df_theta <- fix_input(mat_theta[2, 1, , ])
  colors <- c("True" = "red", "MLE" = "blue")
  df_theta$true <- c(-0.1, 0, 0.1)
  ggplot(df_theta, aes(y = X50., x = true)) +
    geom_point(aes(color = "MLE")) +
    geom_errorbar(aes(ymin = X5., ymax = X95., color = "MLE"), alpha = 0.5) +
    geom_point(aes(y = true, color = "True")) +
    xlab(TeX("$\\theta$")) +
    ylab(TeX("$\\hat{\\theta}$")) +
    labs(title = TeX("Estimated vs True $\\theta$"), colour = "Legend") +
    scale_color_manual(values = colors)
}

# Figure 8
plot_theta_estimates(varying_theta)

# test for varying eta values
test_eta <- function(B) {
  N <- 1000 # length of ts
  omega2 <- c(0.1, 0.3, 0.8, 2)
  theta <- 0
  phi <- -0.5
  mu <- 1
  sigma <- sqrt(omega2 * (1 - phi^2))
  # do a grid search and save parameters

  # loop through the parameters
  r_mat <- array(dim = c(4, 2, length(sigma), B))
  for (n in 1:length(sigma)) {
    for (b in 1:B) {
      ret <- stochastic_volatility_series_AR1(n = N, theta = theta, phi = phi, mu = mu, sigma = sigma[n])
      y <- ret$y
      ret2 <- mysvm(y)
      sd <- sdreport(ret2$model)
      sum2 <- summary(sd, "report")
      sum1 <- summary(sd)
      # sum1[c(1,4),]
      r_mat[, , n, b] <- rbind(sum1[c(1, 4), ], sum2)
    }
  }
  return(r_mat)
}
# matrix of parameter estimates for eta
varying_eta <- test_eta(B = 100)

# Figure 12
par(mfrow = c(2, 2))
hist(varying_eta[4, 1, 1, ], n = 50, xlim = c(-1, 1), main = TeX("$\\eta^2$ = 0.1, $\\phi = -0.5$"), xlab = TeX("$\\hat{\\phi}$"), freq = FALSE)
abline(v = -0.5, col = "red", )
hist(varying_eta[4, 1, 2, ], n = 50, xlim = c(-1, 1), main = TeX("$\\eta^2$ = 0.3, $\\phi = -0.5$"), xlab = TeX("$\\hat{\\phi}$"), freq = FALSE)
abline(v = -0.5, col = "red")
hist(varying_eta[4, 1, 3, ], n = 50, xlim = c(-1, 1), main = TeX("$\\eta^2$ = 0.8, $\\phi = -0.5$"), xlab = TeX("$\\hat{\\phi}$"), freq = FALSE)
abline(v = -0.5, col = "red")
hist(varying_eta[4, 1, 4, ], n = 50, xlim = c(-1, 1), main = TeX("$\\eta^2$ = 2, $\\phi = -0.5$"), xlab = TeX("$\\hat{\\phi}$"), freq = FALSE)
abline(v = -0.5, col = "red")


# simulate and fit parameter to different instances of phi
test_phi <- function(B = 100, N = 1000) {
  # length of ts
  theta <- 0.1
  # 0.99, 0.95, 0.5, 0, -0.5,
  phi <- c(0.99, 0.95, 0.5, 0, -0.5, -0.95, -0.99)
  mu <- 1
  sigma <- 0.1
  # do a grid search and save parameters

  # loop through the parameters
  r_mat <- array(dim = c(4, 2, length(phi), B))
  for (n in 1:length(phi)) {
    for (b in 1:B) {
      ret <- stochastic_volatility_series_AR1(n = N, theta = theta, phi = phi[n], mu = mu, sigma = sigma)
      y <- ret$y
      ret2 <- mysvm(y)
      sd <- sdreport(ret2$model)
      sum2 <- summary(sd, "report")
      sum1 <- summary(sd)
      # sum1[c(1,4),]
      r_mat[, , n, b] <- rbind(sum1[c(1, 4), ], sum2)
    }
  }
  return(r_mat)
}

# matrix of parameter estimates for phi
varying_phi <- test_phi()

# Plotting varying phi with confint
plot_phi_estimates <- function(mat_phi) {
  # plot values for different phi
  df_phi <- fix_input(mat_phi[4, 1, , ])
  colors <- c("True" = "red", "MLE" = "blue")

  df_phi$true <- c(0.99, 0.95, 0.5, 0, -0.5, -0.95, -0.99)
  ggplot(df_phi, aes(y = X50., x = true)) +
    geom_point(aes(color = "MLE")) +
    geom_errorbar(aes(ymin = X5., ymax = X95., color = "MLE"), alpha = 0.5) +
    geom_point(aes(y = true, color = "True")) +
    xlab(TeX("$\\phi$")) +
    ylab(TeX("$\\hat{\\phi}$")) +
    labs(title = TeX("Estimated vs True $\\phi$"), colour = "Legend") +
    scale_color_manual(values = colors)
}

# Figure 11
plot_phi_estimates(varying_phi)

# simulate and fit parameter to different instances of sigma
test_sigma <- function(B = 100) {
  N <- 1000 # length of ts
  theta <- 0.1
  phi <- 0.95
  mu <- 1
  sigma <- c(2, 1, 0.5, 0.1, 0.01)
  # do a grid search and save parameters

  # loop through the parameters
  r_mat <- array(dim = c(4, 2, length(sigma), B))
  for (n in 1:length(sigma)) {
    for (b in 1:B) {
      ret <- stochastic_volatility_series_AR1(n = N, theta = theta, phi = phi, mu = mu, sigma = sigma[n])
      y <- ret$y
      ret2 <- mysvm(y)
      sd <- sdreport(ret2$model)
      sum2 <- summary(sd, "report")
      sum1 <- summary(sd)
      # sum1[c(1,4),]
      r_mat[, , n, b] <- rbind(sum1[c(1, 4), ], sum2)
    }
  }
  return(r_mat)
}

# matrix of parameter estimates for sigma
varying_sigma <- test_sigma()

# Plotting varying sigma with confint
plot_sigma_estimates <- function(mat_sigma) {
  df_sigma <- fix_input(mat_sigma[3, 1, , ])
  colors <- c("True" = "red", "MLE" = "blue")
  df_sigma$true <- c(2, 1, 0.5, 0.1, 0.01)
  ggplot(df_sigma, aes(y = X50., x = true)) +
    geom_point(aes(color = "MLE")) +
    geom_errorbar(aes(ymin = X5., ymax = X95., color = "MLE"), alpha = 0.5) +
    geom_point(aes(y = true, color = "True")) +
    xlab(TeX("$\\sigma$")) +
    ylab(TeX("$\\hat{\\sigma}$")) +
    labs(title = TeX("Estimated vs True $\\sigma$"), colour = "Legend") +
    scale_color_manual(values = colors)
}

# Figure 10
plot_sigma_estimates(varying_sigma)

# simulate and fit parameter to different instances of mu
test_mu <- function(B = 100) {
  N <- 1000 # length of ts
  theta <- 0.1
  phi <- 0.95
  mu <- c(10, 1, 0, -1, -10)
  sigma <- 0.1
  # do a grid search and save parameters

  # loop through the parameters
  r_mat <- array(dim = c(4, 2, length(mu), B))
  for (n in 1:length(mu)) {
    for (b in 1:B) {
      ret <- stochastic_volatility_series_AR1(n = N, theta = theta, phi = phi, mu = mu[n], sigma = sigma)
      y <- ret$y
      ret2 <- mysvm(y)
      sd <- sdreport(ret2$model)
      sum2 <- summary(sd, "report")
      sum1 <- summary(sd)
      # sum1[c(1,4),]
      r_mat[, , n, b] <- rbind(sum1[c(1, 4), ], sum2)
    }
  }
  return(r_mat)
}

# matrix of parameter estimates for mu
varying_mu <- test_mu()

# Plotting varying mu with confint
plot_mu_estimates <- function(mat_mu) {
  df_mu <- fix_input(mat_mu[1, 1, , ])
  colors <- c("True" = "red", "MLE" = "blue")

  df_mu$true <- c(10, 1, 0, -1, -10)
  ggplot(df_mu, aes(y = X50., x = true)) +
    geom_point(aes(color = "MLE")) +
    geom_errorbar(aes(ymin = X5., ymax = X95., color = "MLE"), alpha = 0.5) +
    geom_point(aes(y = true, color = "True")) +
    xlab(TeX("$\\mu$")) +
    ylab(TeX("$\\hat{\\mu}$")) +
    labs(title = TeX("Estimated vs True $\\mu$"), colour = "Legend") +
    scale_color_manual(values = colors)
}

# Figure 9
plot_mu_estimates(varying_mu)


#### Validation function definitions ####

# function that splits into rolling CV sets
validation_splits <- function(mat, len_val = 100, init_length = 500) {
  train <- c()
  test <- c()
  n <- length(mat)
  length_val <- n - init_length
  s <- seq(from = init_length, to = n, by = len_val)
  for (i in 1:(length(s) - 1)) {
    tr <- mat[1:s[i]]
    length(tr) <- n
    val <- mat[(s[i] + 1):s[i + 1]]
    length(val) <- len_val
    train <- cbind(train, tr)
    test <- cbind(test, val)
  }

  if (s[i + 1] != n) {
    tr <- mat[1:s[i + 1]]
    length(tr) <- n
    val <- mat[(s[i + 1] + 1):n]
    length(val) <- len_val
    train <- cbind(train, tr)
    test <- cbind(test, val)
  }


  return(list(train = train, val = test))
}

# logit(x)
glogitinv <- function(x) (exp(x) - 1) / (1 + exp(x))

# Calculating conditional log-lik for SVM
validate_SVM <- function(train, train_val) {
  m1 <- mysvm(train)
  thetaHat <- m1$fit$par

  ll1 <- m1$model$fn(thetaHat)
  ll1
  m2 <- makeSVM(train_val)
  ll2 <- m2$fn(thetaHat)
  # likelihood is negative by construction
  cond_ll <- -ll2 + ll1
  cond_ll
  parameters <- m1$fit$par
  parameters[2] <- glogitinv(parameters[2])
  parameters[3] <- exp(parameters[3])
  names(parameters) <- c("mu", "phi", "sigma", "theta")
  return(list(ll = cond_ll, parameters = parameters))
}

# log-lik GARCH
garch11ll <- function(theta, y) {
  omega <- theta[1]
  alpha1 <- theta[2]
  beta1 <- theta[3]
  thetaHat <- theta[4]
  n <- length(y)

  ll <- 0
  sigma2 <- numeric(n)
  epsilon <- numeric(n)
  epsilon[1] <- sqrt(omega)
  sigma2[1] <- omega
  epsilon[2:n] <- y[2:n] - y[1:(n - 1)] - thetaHat
  for (t in 2:n) {
    sigma2[t] <- omega + alpha1 * epsilon[t - 1]^2 + beta1 * sigma2[t - 1] # conditional variance
    ll <- ll + dnorm(epsilon[t], 0, sqrt(sigma2[t]), log = TRUE) # f(epsilon_t | history up to time t-1)
  }
  return(ll)
}

# log-lik GARCH function for unconstrained parameters
garch11ll_unconstrained <- function(theta, y) {
  omega <- exp(theta[1])
  alpha1 <- exp(theta[2]) / (1 + exp(theta[3]) + exp(theta[2]))
  beta1 <- exp(theta[3]) / (1 + exp(theta[3]) + exp(theta[2]))
  thetaHat <- theta[4]
  n <- length(y)

  ll <- 0
  sigma2 <- numeric(n)
  epsilon <- numeric(n)
  epsilon[1] <- sqrt(omega)
  sigma2[1] <- omega
  epsilon[2:n] <- y[2:n] - y[1:(n - 1)] - thetaHat
  for (t in 2:n) {
    sigma2[t] <- omega + alpha1 * epsilon[t - 1]^2 + beta1 * sigma2[t - 1] # conditional variance
    ll <- ll + dnorm(epsilon[t], 0, sqrt(sigma2[t]), log = TRUE) # f(epsilon_t | history up to time t-1)
  }
  return(ll)
}


# Convert GARCH unconstrain values to true values
theta_convert <- function(theta_raw) {
  omega <- exp(theta_raw[1])
  alpha <- exp(theta_raw[2]) / (1 + exp(theta_raw[3]) + exp(theta_raw[2]))
  beta <- exp(theta_raw[3]) / (1 + exp(theta_raw[3]) + exp(theta_raw[2]))
  thet <- theta_raw[4]
  return(c(omega, alpha, beta, thet))
}

myGARCH11 <- function(y, init = c(0, 0, 0, 0)) {
  r <- optim(init, fn = garch11ll_unconstrained, y = y, control = list(fnscale = -1), hessian = TRUE)
  ML_estimates <- theta_convert(r$par)
  fisher_info <- solve(-r$hessian)
  var_delta <- deltavar(c(exp(omega), exp(alpha) / (1 + exp(alpha) + exp(beta)), exp(beta) / (1 + exp(alpha) + exp(beta)), theta),
    meanval = c(omega = r$par[1], alpha = r$par[2], beta = r$par[3], theta = r$par[4]),
    Sigma = fisher_info
  )
  return(list(optim_ret = r, MLE = ML_estimates, SD = sqrt(var_delta)))
}


# Calculating conditional log-lik for GARCH
validate_GARCH <- function(train, train_val) {
  # RUGARCH til å finne parametere
  # hva skal vi gjøre med RW?

  # finne theta først
  ret <- myGARCH11(train)
  ll1 <- ret$optim_ret$value
  theta <- ret$MLE
  names(theta) <- c("omega", "alpha", "beta", "theta")
  ll2 <- garch11ll(theta, train_val)
  ll <- ll2 - ll1
  return(list(ll = ll, parameters = theta))
}




# Calculating M(n,n+h)
rolling_validate <- function(val_list) {
  cond_ll_svm <- c()
  cond_ll_garch <- c()
  param_svm <- c()
  param_garch <- c()
  for (t in 1:ncol(val_list$train)) {
    t1 <- val_list$train[, t]
    val <- val_list$val[, t]
    length(t1) <- length(t1) - sum(is.na(t1))
    length(val) <- length(val) - sum(is.na(val))
    t2 <- c(t1, val)
    svm <- validate_SVM(t1, t2)
    garch <- validate_GARCH(t1, t2)
    cond_ll_svm <- c(cond_ll_svm, svm$ll)
    cond_ll_garch <- c(cond_ll_garch, garch$ll)
    param_svm <- rbind(param_svm, svm$parameters)
    param_garch <- rbind(param_garch, garch$parameters)
  }

  return(list(svm_ll = cond_ll_svm, garch_ll = cond_ll_garch, svm_param = param_svm, garch_param = param_garch))
}

##### Tesla stock #####

# retrieve daily tesla data
tesla <- getSymbols("TSLA", scr = "yahoo", auto.assign = FALSE)
tesla <- tesla[1:(nrow(tesla) - 16), ]
tail(tesla)

# log-scale
log_tsla <- log(tesla$TSLA.Close)[(nrow(tesla) - 999):nrow(tesla)]

plot_ts <- function(df) {
  colnames(df) <- "daily_return"
  df$date <- rownames(df)
  df$date <- as.Date(df$date)
  ggplot(data = df, aes(y = daily_return, x = date)) +
    geom_line()
}

# plot log-daily closing prices of Tesla stock
p_tesla <- plot_ts(data.frame(log_tsla))
p_tesla <- p_tesla + labs(title = "log Tesla Stock Price", x = "Date", y = "log(Price)")


# make matrix
mat_tesla <- matrix(log_tsla)

# fit SVM and return MLE
estimate_parameters <- function(mat) {
  ret <- mysvm(mat)
  r <- sdreport(ret$model)
  sum1 <- summary(r)[c(1, 4), ]
  sum2 <- summary(r, "report")
  estimates <- rbind(sum1, sum2)
  return(estimates)
}
# MLE
estimates_tesla <- estimate_parameters(mat_tesla)

tesla_params_SVM <- estimates_tesla[, 1]




# Split into training av validation set Tesla
splits_tesla10 <- validation_splits(mat_tesla, len_val = 10, init_length = 500)

# Perform validation
tsla_lls10 <- rolling_validate(splits_tesla10)

# plot diff log-lik
df <- data.frame(svm = tsla_lls10$svm_ll, garch = tsla_lls10$garch_ll)


# plot parameter convergence
df_svm <- data.frame(tsla_lls10$svm_param)
df_garch <- data.frame(tsla_lls10$svm_param)

# df = data.frame(para = tsla_lls10$svm_ll, garch = tsla_lls10$garch_ll)



# geom_line(aes(y = cumsum(garch), colour = "garch"))


# final parameter values


final_sgarch_tesla <- myGARCH11(mat_tesla) # In table 4
parameters_final <- estimate_parameters(mat_tesla) # In table 4


#### SPY #####

SPY <- getSymbols("SPY", scr = "yahoo", auto.assign = FALSE)
SPY <- SPY[1:(nrow(SPY) - 16), ]
tail(SPY)
# plot the returns
log_spy <- log(SPY$SPY.Close[(nrow(SPY) - 999):nrow(SPY)])

# plot log-daily closing prices of SPY ETF
p_spy <- plot_ts(data.frame(log_spy))
p_spy <- p_spy + labs(title = "log SPY ETF Price", x = "Date", y = "log(Price)")

# Figure 6
grid.arrange(p_tesla, p_spy, ncol = 1)

# make matrix
mat_spy <- matrix(log_spy)


# Split into training av validation set SPY
splits_SPY10 <- validation_splits(mat_spy, len_val = 10, init_length = 500)
SPY_lls10 <- rolling_validate(splits_SPY10)


svm_par_spy <- estimate_parameters(mat_spy) # In table 4
final_sgarch_spy <- myGARCH11(mat_spy, init = c(-1, -1, -1, 0)) # In table 4



# plot diff log-lik

# Tesla
tesla_svm_vs_garch <- data.frame(svm = tsla_lls10$svm_ll, garch = tsla_lls10$garch_ll)
p1 <- ggplot(tesla_svm_vs_garch, aes(x = seq(500, 990, 10))) +
  geom_line(aes(y = svm - garch), ) +
  geom_hline(yintercept = 0, color = "red") +
  labs(x = "n", y = "M", title = TeX("Tesla: $M(n, n+10)$"))

# SPY
SPY_svm_vs_garch <- data.frame(svm = SPY_lls10$svm_ll, garch = SPY_lls10$garch_ll)
p2 <- ggplot(SPY_svm_vs_garch, aes(x = seq(500, 990, 10))) +
  geom_line(aes(y = svm - garch), ) +
  geom_hline(yintercept = 0, color = "red") +
  labs(x = "n", y = "M", title = TeX("SPY: $M(n, n+10)$"))


L_tesla <- apply(tesla_svm_vs_garch, 2, sum)
L_spy <- apply(SPY_svm_vs_garch, 2, sum)

SCORE_Table <- data.frame(rbind(L_tesla, L_spy))
SCORE_Table[, 3] <- SCORE_Table[, 1] - SCORE_Table[, 2]
colnames(SCORE_Table) <- c("SVM", "GARCH", "SCORE")
rownames(SCORE_Table) <- c("Tesla", "SPY")

# Table 3
SCORE_Table

# plot parameter convergence
tsla_svm <- data.frame(tsla_lls10$svm_param)
tsla_garch <- data.frame(tsla_lls10$garch_param)

spy_svm <- data.frame(SPY_lls10$svm_param)
spy_garch <- data.frame(SPY_lls10$garch_param)


svm_plot_valid_mle <- function(df_svm) {
  p_mu <- ggplot(df_svm, aes(x = seq(500, 990, 10))) +
    geom_line(aes(y = mu)) +
    labs(x = "observations, n", y = TeX("$\\hat{\\mu}$"))
  p_phi <- ggplot(df_svm, aes(x = seq(500, 990, 10))) +
    geom_line(aes(y = phi)) +
    labs(x = "observations, n", y = TeX("$\\hat{\\phi}$"))

  p_sigma <- ggplot(df_svm, aes(x = seq(500, 990, 10))) +
    geom_line(aes(y = sigma)) +
    labs(x = "observations, n", y = TeX("$\\hat{\\sigma}$"))

  p_theta <- ggplot(df_svm, aes(x = seq(500, 990, 10))) +
    geom_line(aes(y = theta)) +
    labs(x = "observations, n", y = TeX("$\\hat{\\theta}$"))
  grid.arrange(p_mu, p_phi, p_sigma, p_theta, ncol = 2)
}
# Figure 14
svm_plot_valid_mle(tsla_svm)
# Figure 16
svm_plot_valid_mle(spy_svm)

garch_plot_valid_mle <- function(df_garch) {
  p_omega <- ggplot(df_garch, aes(x = seq(500, 990, 10))) +
    geom_line(aes(y = omega)) +
    labs(x = "observations, n", y = TeX("$\\hat{\\omega}$"))
  p_alpha <- ggplot(df_garch, aes(x = seq(500, 990, 10))) +
    geom_line(aes(y = alpha)) +
    labs(x = "observations, n", y = TeX("$\\hat{\\alpha}$"))

  p_beta <- ggplot(df_garch, aes(x = seq(500, 990, 10))) +
    geom_line(aes(y = beta)) +
    labs(x = "observations, n", y = TeX("$\\hat{\\beta}$"))

  p_theta <- ggplot(df_garch, aes(x = seq(500, 990, 10))) +
    geom_line(aes(y = theta)) +
    labs(x = "observations, n", y = TeX("$\\hat{\\theta}$"))
  grid.arrange(p_omega, p_alpha, p_beta, p_theta, ncol = 2)
}
# Figure 15
garch_plot_valid_mle(tsla_garch)
# Figure 16
garch_plot_valid_mle(spy_garch)
