## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
library(ccostr)
library(ggplot2)
library(knitr)
library(parallel)
library(msm)

## ------------------------------------------------------------------------
set.seed(123)
sim <- simCostData(n = 1000, dist = "unif", censor = "light", cdist = "exp", L = 10)
est <- ccmean(sim$censoredCostHistory)
est

## ----fig.height=3, fig.width=7-------------------------------------------
plot(est) + geom_hline(yintercept = 40000, linetype = "dotted", size = 1)

## ---- eval = FALSE-------------------------------------------------------
#  
#  nSim   <- 1000
#  nYears <- 10
#  indv   <- 1000 # increating individuals increases computing time exponential
#  ## true mean for unif is 40000 and exp is 35956
#  unif_light <- lapply(1:nSim, function(x) simCostData(n = indv, dist = "unif", censor = "light", cdist = "exp", L = nYears))
#  unif_heavy <- lapply(1:nSim, function(x) simCostData(n = indv, dist = "unif", censor = "heavy", cdist = "exp", L = nYears))
#  exp_light  <- lapply(1:nSim, function(x) simCostData(n = indv, dist = "exp",  censor = "light", cdist = "exp", L = nYears))
#  exp_heavy  <- lapply(1:nSim, function(x) simCostData(n = indv, dist = "exp",  censor = "heavy", cdist = "exp", L = nYears))

## ---- eval = FALSE-------------------------------------------------------
#  nCores <- parallel::detectCores() - 1
#  cl <- makeCluster(nCores)
#  clusterExport(cl = cl, c("unif_light", "unif_heavy", "exp_light", "exp_heavy"))
#  invisible(clusterEvalQ(cl = cl, {library(dplyr)
#                                   library(ccostr)
#                                   library(data.table)
#                                   library(survival)}))
#  est_unif_light <- parLapply(cl, unif_light, function(x) ccmean(x$censoredCostHistory, L = 10))
#  est_unif_heavy <- parLapply(cl, unif_heavy, function(x) ccmean(x$censoredCostHistory, L = 10))
#  est_exp_light  <- parLapply(cl, exp_light,  function(x) ccmean(x$censoredCostHistory, L = 10))
#  est_exp_heavy  <- parLapply(cl, exp_heavy,  function(x) ccmean(x$censoredCostHistory, L = 10))
#  stopCluster(cl)

## ---- eval = FALSE-------------------------------------------------------
#  results_unif_light <- do.call(rbind, lapply(est_unif_light, function(x) x[[3]]))
#  results_unif_heavy <- do.call(rbind, lapply(est_unif_heavy, function(x) x[[3]]))
#  results_exp_light  <- do.call(rbind, lapply(est_exp_light,  function(x) x[[3]]))
#  results_exp_heavy  <- do.call(rbind, lapply(est_exp_heavy,  function(x) x[[3]]))
#  results_true <- data.frame("unif_light" = 40000,
#                             "unif_heavy" = 40000,
#                             "exp_light"  = 35956,
#                             "exp_heavy"  = 35956)
#  results_bias <- data.frame("unif_light" = (colMeans(results_unif_light)),
#                             "unif_heavy" = (colMeans(results_unif_heavy)),
#                             "exp_light"  = (colMeans(results_exp_light)),
#                             "exp_heavy"  = (colMeans(results_exp_heavy)))
#  results <- rbind(results_true, results_bias)
#  row.names(results) <- c("true_mean", colnames(results_unif_light))
#  
#  results_bias <- cbind(round(results[,c(1,2)] - 40000,2),
#                        round(results[,c(3,4)] - 35956,2))
#  
#  cov_unif_light <- do.call(rbind, lapply(est_unif_light, function(x) ifelse(x[[4]][5,] >= 40000 & x[[4]][4,] <= 40000, 1, 0)))
#  cov_unif_heavy <- do.call(rbind, lapply(est_unif_heavy, function(x) ifelse(x[[4]][5,] >= 40000 & x[[4]][4,] <= 40000, 1, 0)))
#  cov_exp_light  <- do.call(rbind, lapply(est_exp_light,  function(x) ifelse(x[[4]][5,] >= 35956 & x[[4]][4,] <= 35956, 1, 0)))
#  cov_exp_heavy  <- do.call(rbind, lapply(est_exp_heavy,  function(x) ifelse(x[[4]][5,] >= 35956 & x[[4]][4,] <= 35956, 1, 0)))
#  results_coverage <- data.frame("unif_light" = (colMeans(cov_unif_light, na.rm = T)),
#                                 "unif_heavy" = (colMeans(cov_unif_heavy, na.rm = T)),
#                                 "exp_light"  = (colMeans(cov_exp_light,  na.rm = T)),
#                                 "exp_heavy"  = (colMeans(cov_exp_heavy,  na.rm = T)))

## ------------------------------------------------------------------------
load(system.file("extdata", "results.Rdata", package = "ccostr"))
kable(results)

## ------------------------------------------------------------------------
kable(results_bias)

## ------------------------------------------------------------------------
kable(results_coverage, digits = 3)

