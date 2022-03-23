mixed.solve=function (y, Z = NULL, K = NULL, X = NULL, method = "REML", bounds = c(1e-09, 
                                                                                   1e+09), SE = FALSE, return.Hinv = FALSE) 
{
  pi <- 3.14159
  n <- length(y)
  y <- matrix(y, n, 1)
  not.NA <- which(!is.na(y))
  if (is.null(X)) {
    p <- 1
    X <- matrix(rep(1, n), n, 1)
  }
  p <- ncol(X)
  if (is.null(p)) {
    p <- 1
    X <- matrix(X, length(X), 1)
  }
  if (is.null(Z)) {
    Z <- diag(n)
  }
  m <- ncol(Z)
  if (is.null(m)) {
    m <- 1
    Z <- matrix(Z, length(Z), 1)
  }
  stopifnot(nrow(Z) == n)
  stopifnot(nrow(X) == n)
  if (!is.null(K)) {
    stopifnot(nrow(K) == m)
    stopifnot(ncol(K) == m)
  }
  Z <- as.matrix(Z[not.NA, ])
  X <- as.matrix(X[not.NA, ])
  n <- length(not.NA)
  y <- matrix(y[not.NA], n, 1)
  XtX <- crossprod(X, X)
  rank.X <- qr(XtX)$rank
  if (rank.X < p) {
    stop("X not full rank")
  }
  XtXinv <- solve(XtX)
  S <- diag(n) - tcrossprod(X %*% XtXinv, X)
  if (n <= m + p) {
    spectral.method <- "eigen"
  }
  else {
    spectral.method <- "cholesky"
    if (!is.null(K)) {
      diag(K) <- diag(K) + 1e-06
      B <- try(chol(K), silent = TRUE)
      if (inherits(B, what = "try-error")) {
        stop("K not positive semi-definite.")
      }
    }
  }
  if (spectral.method == "cholesky") {
    if (is.null(K)) {
      ZBt <- Z
    }
    else {
      ZBt <- tcrossprod(Z, B)
    }
    svd.ZBt <- svd(ZBt, nu = n)
    U <- svd.ZBt$u
    phi <- c(svd.ZBt$d^2, rep(0, n - m))
    SZBt <- S %*% ZBt
    svd.SZBt <- try(svd(SZBt), silent = TRUE)
    if (inherits(svd.SZBt, what = "try-error")) {
      svd.SZBt <- svd(SZBt + matrix(1e-10, nrow = nrow(SZBt), 
                                    ncol = ncol(SZBt)))
    }
    QR <- qr(cbind(X, svd.SZBt$u))
    Q <- qr.Q(QR, complete = TRUE)[, (p + 1):n]
    R <- qr.R(QR)[p + 1:m, p + 1:m]
    ans <- try(solve(t(R^2), svd.SZBt$d^2), silent = TRUE)
    if (inherits(ans, what = "try-error")) {
      spectral.method <- "eigen"
    }
    else {
      theta <- c(ans, rep(0, n - p - m))
    }
  }
  if (spectral.method == "eigen") {
    offset <- sqrt(n)
    if (is.null(K)) {
      Hb <- tcrossprod(Z, Z) + offset * diag(n)
    }
    else {
      Hb <- tcrossprod(Z %*% K, Z) + offset * diag(n)
    }
    Hb.system <- eigen(Hb, symmetric = TRUE)
    phi <- Hb.system$values - offset
    if (min(phi) < -1e-06) {
      stop("K not positive semi-definite.")
    }
    U <- Hb.system$vectors
    SHbS <- S %*% Hb %*% S
    SHbS.system <- eigen(SHbS, symmetric = TRUE)
    theta <- SHbS.system$values[1:(n - p)] - offset
    Q <- SHbS.system$vectors[, 1:(n - p)]
  }
  omega <- crossprod(Q, y)
  omega.sq <- omega^2
  if (method == "ML") {
    f.ML <- function(lambda, n, theta, omega.sq, phi) {
      n * log(sum(omega.sq/(theta + lambda))) + sum(log(phi + 
                                                          lambda))
    }
    soln <- optimize(f.ML, interval = bounds, n, theta, omega.sq, 
                     phi)
    lambda.opt <- soln$minimum
    df <- n
  }
  else {
    f.REML <- function(lambda, n.p, theta, omega.sq) {
      n.p * log(sum(omega.sq/(theta + lambda))) + sum(log(theta + 
                                                            lambda))
    }
    soln <- optimize(f.REML, interval = bounds, n - p, theta, 
                     omega.sq)
    lambda.opt <- soln$minimum
    df <- n - p
  }
  Vu.opt <- sum(omega.sq/(theta + lambda.opt))/df
  Ve.opt <- lambda.opt * Vu.opt
  Hinv <- U %*% (t(U)/(phi + lambda.opt))
  W <- crossprod(X, Hinv %*% X)
  beta <- array(solve(W, crossprod(X, Hinv %*% y)))
  rownames(beta) <- colnames(X)
  if (is.null(K)) {
    KZt <- t(Z)
  }
  else {
    KZt <- tcrossprod(K, Z)
  }
  KZt.Hinv <- KZt %*% Hinv
  u <- array(KZt.Hinv %*% (y - X %*% beta))
  if (is.null(K)) {
    rownames(u) <- colnames(Z)
  }
  else {
    rownames(u) <- rownames(K)
  }
  LL = -0.5 * (soln$objective + df + df * log(2 * pi/df))
  if (!SE) {
    if (return.Hinv) {
      return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, 
                  u = u, LL = LL, Hinv = Hinv))
    }
    else {
      return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, 
                  u = u, LL = LL))
    }
  }
  else {
    Winv <- solve(W)
    beta.SE <- array(sqrt(Vu.opt * diag(Winv)))
    rownames(beta.SE) <- rownames(beta)
    WW <- tcrossprod(KZt.Hinv, KZt)
    WWW <- KZt.Hinv %*% X
    if (is.null(K)) {
      u.SE <- array(sqrt(Vu.opt * (rep(1, m) - diag(WW) + 
                                     diag(tcrossprod(WWW %*% Winv, WWW)))))
    }
    else {
      u.SE <- array(sqrt(Vu.opt * (diag(K) - diag(WW) + 
                                     diag(tcrossprod(WWW %*% Winv, WWW)))))
    }
    rownames(u.SE) <- rownames(u)
    if (return.Hinv) {
      return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, 
                  beta.SE = beta.SE, u = u, u.SE = u.SE, LL = LL, 
                  Hinv = Hinv))
    }
    else {
      return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, 
                  beta.SE = beta.SE, u = u, u.SE = u.SE, LL = LL))
    }
  }
}


#library(tidyverse)
#library(rrBLUP)
library(data.table)

t = commandArgs(trailingOnly = T)
t = as.numeric(t)
set.seed(t)

avg_merge=read.csv("data/2020merge_avg_narm.csv")

avg_gp = avg_merge[,-c(13,17,19:24)]
avg_gp = avg_gp[order(avg_gp[,1]),]

################ only for 20220222 ###############
avg_gp <- avg_gp[,c("GenotypeID", "TotalGrainMassGrams_NEmean","Yield_MImean")]

# read marker
setDTthreads(threads = 0, restore_after_fork = NULL, throttle = NULL)
snp_markers  = fread("WiDiv_maf005_725geno.num.txt")
snp_markers = snp_markers[order(snp_markers[,1]),]
snp_markers2 = snp_markers[1:nrow(snp_markers),2:ncol(snp_markers)]


rrblup_full = function(phe, i, snp_markers2){
  # phe is the dataframe contains phenotypic data, i is the column of phe, snp_markers2 is a matrix
  pheno2 = as.matrix(phe[,i]) # has to be matrix, even though one column, otherwise throws errors when splitting data
  
  pheno_training_data = as.matrix(pheno2)
  snp_training_data = snp_markers2 # Converting to matrix here causes a vector as result which creates problems in model
  
  
  snp_training_data = sapply(snp_training_data, as.numeric) # Matrix must be numeric
  snp_training_data = as.matrix(snp_training_data)
  
  trained_model = mixed.solve(y=pheno_training_data,Z = snp_training_data)
  
  marker_effects = as.matrix(trained_model$u) # marker effects center around zero
  BLUE = as.vector(trained_model$beta) # BLUE is baseline of which marker effects center
  
  predicted_train = as.matrix(snp_training_data) %*% marker_effects # rrBLUPS
  return(as.vector((predicted_train[,1])+BLUE))
}

#do full rrblup
full_train_result=avg_gp
full_train_result$predicted_yield = rrblup_full(avg_gp, 2, snp_markers2)

write.table(full_train_result,"data/rrblupresult/rrblup_full.csv",quote = F,row.names = F,sep = ",")
# 
# avg_yield <- avg_gp[,c("GenotypeID", "TotalGrainMassGrams_NEmean","Yield_MImean")]
# full_yield_result <- avg_yield
# full_yield_result[2] <- rrblup_full(avg_yield, 2, snp_markers2)

# bootstrap 80% for 10 times
g_number = nrow(avg_gp)
full_train_result_bs=list()
corr_ob_full = c()
corr_gp_full = c()
for (t in 1:20) {
  sp <- sample(1:g_number,as.integer(0.8*g_number))
  snp_sp <- snp_markers2[sp,]
  avg_gp_bs = avg_gp[sp,]
  full_train_result_bs[[t]]=avg_gp_bs
  for (i in 2:ncol(avg_gp)) {
    full_train_result_bs[[t]][c(1:length(sp)),i] = rrblup_full(avg_gp_bs, i, snp_sp)
  }
  #write.table(full_train_result_bs[[t]],paste0("data/rrblupresult/NE_rrblupfull_bs.csv",t),quote = F,row.names = F,sep = ",")
  # correlations
  MIyield = avg_gp_bs$Yield_MImean
  NEyield = avg_gp_bs$TotalGrainMassGrams_NEmean
  GPyield = full_train_result_bs[[t]]$TotalGrainMassGrams_NEmean
  yield_comparisons = data.frame(MIyield=MIyield,NEyield=NEyield,GPyield=GPyield)
  write.table(yield_comparisons,paste0("data/rrblupresult/NE_rrblupfullbs.csv",t),quote=F,col.names = T,row.names = F, sep=",")
  #corr_ob_full[t] = cor(NEyield, MIyield)
  #corr_gp_full[t] = cor(GPyield, MIyield)
}
#corrs = data.frame(observedvalue = corr_ob_full, gpvalue = corr_gp_full)
#write.table(corrs, "data/rrblupresult/rrblup_full_corrs.csv",quote=F,row.names = F,sep=",")

##############################################################################
# cross validation 5-fold 20 times

rrblup_cv = function(phe, i, snpmarkers2, k) {
  g_number = nrow(phe)
  random.order = sample(1:g_number,g_number)
  validate = matrix(random.order, ncol = k)
  predicted_cv_result=phe
  for (k in 1:k) {
    pheno2 = as.matrix(phe[,i])

    validation_entries = validate[,k] # select the the validation entries
    training_entries = as.vector(validate[,-k])

    pheno_training_data = pheno2[training_entries,]
    snp_training_data = snp_markers2[training_entries,] # Converting to matrix here causes a vector as result which creates problems in model

    pheno_validate_data = pheno2[validation_entries,]
    snp_validate_data = snp_markers2[validation_entries,]

    snp_training_data = sapply(snp_training_data, as.numeric) # Matrix must be numeric
    snp_training_data = as.matrix(snp_training_data)

    snp_validate_data = sapply(snp_validate_data, as.numeric) # Matrix must be numeric
    snp_validate_data = as.matrix(snp_validate_data,K=NULL)

    trained_model = mixed.solve(y=pheno_training_data,Z = snp_training_data)
    summary(trained_model)

    marker_effects = as.matrix(trained_model$u) # marker effects center around zero
    BLUE = as.vector(trained_model$beta) # BLUE is baseline of which marker effects center

    predicted_train = as.matrix(snp_training_data) %*% marker_effects # rrBLUPS
    predicted_validate = as.matrix(snp_validate_data) %*% marker_effects   # rrBLUPS

    predicted_train_result = as.vector((predicted_train[,1])+BLUE)
    predicted_cv_result[validation_entries,i] = as.vector((predicted_validate[,1])+BLUE)
  }
  return(predicted_cv_result[,i])
}


#### for 0224 only
avg_gp <- avg_yield

predicted_cv_list=list()
corr_gp_cv = c()
corr_gp_cvNE = c()
for (t in 1:10) {
  set.seed(t)
  predicted_cv_list[[t]]=avg_gp
  for (i in 2:(ncol(avg_gp)-1)) {
    predicted_cv_list[[t]]$Yield_NEpredicted = rrblup_cv(avg_gp, i, snp_markers2, 5)
  }
  write.table(predicted_cv_list[[t]],paste0("data/rrblupresult/NE_rrblupcv.csv",t),quote=F,col.names = T,row.names = F,sep = ",")
  # correlations
  # MIyield = avg_gp$Yield_MImean
  # NEyield = avg_gp$TotalGrainMassGrams_NEmean
  # GPyield = predicted_cv_list[[t]]$TotalGrainMassGrams_NEmean
  # yield_comparisons = data.frame(MIyield=MIyield,NEyield=NEyield,GPyield=GPyield)
  # write.table(yield_comparisons,paste0("data/rrblupresult/2020NEMI_rrblupcv_comparisons.csv",t),quote=F,col.names = T,row.names = F, sep=",")
  # corr_gp_cv[t] = cor(GPyield, MIyield)
  # corr_gp_cvNE[t] = cor(GPyield, NEyield)
}
#corrs = data.frame(gpcvMI = corr_gp_cv, gpcvNE = corr_gp_cvNE)
#write.table(corrs, "data/rrblupresult/rrblup_cv_corrs.csv",quote=F,row.names = F,sep=",")
