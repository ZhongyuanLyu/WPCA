########## Missing blocks ##############
# train_ratio <- 0.8
# block_row <- round(sqrt(1-train_ratio) * N)
# block_col <- round(sqrt(1-train_ratio) * Times)
# start_row <- sample(1:(N - block_row + 1), 1)
# start_col <- sample(1:(Times - block_col + 1), 1)
# train_mat <- matrix(1, nrow = N, ncol = Times)
# train_mat[start_row:(start_row + block_row - 1), 
#            start_col:(start_col + block_col - 1)] <- 0
# X_train <- X*train_mat
# X_test <- X*(1-train_mat)



# data <- fredqd(file = "FRED-QD-2025-02.csv", date_start = NULL, date_end = NULL, transform = FALSE)

# col_na_prop <- apply(is.na(data), 2, mean)
# data_select <- data[, (col_na_prop < 0.05)]
# data_bal <- na.omit(data_select)
# rownames(data_bal) <- data_bal[,1]
# X_bal <- data_bal[,2:ncol(data_bal)]
# X <- t(as.matrix(X_bal))
# N <- nrow(X)
# Times <- ncol(X)
# X <- scale(X)
# 
# ########## Missing at random ##############
# train_ratio <- 0.7
# train_mat <- matrix(runif(N * Times) < train_ratio, nrow = N, ncol = Times)
# X_train <- X*train_mat
# X_test <- X*(1-train_mat)
# 
# ########### APCA ##############
# r <- 3
# svd_res <- svd(X_train, nu = r, nv = r)
# if (r>1){
#   X_PCA <- svd_res$u %*%  diag(svd_res$d[1:r]) %*% t(svd_res$v)/train_ratio
# } else {
#   X_PCA <- svd_res$u %*% t(svd_res$v) * svd_res$d[r] /train_ratio
# }
# 
# U_heteroPCA <- heteroPCA(X_train, r, T0 = 10)
# 
# res_APCA <- CV_APCA(X_train, r, p_star = 0.8, J = 10)
# # res_APCA$gamma
# # res_APCA$err
# 
# sum((X_test - (res_APCA$U %*% t(res_APCA$U) %*% X_train)/train_ratio * (1-train_mat))^2)/sum(X_test^2)
# sum((X_test - X_PCA * (1-train_mat))^2)/sum(X_test^2)
# sum((X_test - (U_heteroPCA %*%  t(U_heteroPCA) %*% X_train)/train_ratio * (1-train_mat))^2)/sum(X_test^2)



library(readr)
source("functions.R")
set.seed(1234)
# data <- fredmd(file = "FRED-MD-2025-02.csv", date_start = NULL, date_end = NULL, transform = FALSE)
data <- fredqd(file = "FRED-QD-2025-02.csv", date_start = NULL, date_end = NULL, transform = FALSE)
col_na_prop <- apply(is.na(data), 2, mean)
data_select <- data[, (col_na_prop < 0.05)]
data_bal <- na.omit(data_select)
rownames(data_bal) <- data_bal[,1]
X_bal <- data_bal[,2:ncol(data_bal)]
X <- t(as.matrix(X_bal))
N <- nrow(X)
Times <- ncol(X)
# est_r <- ratio_based_r(X, 100)

one_sim <- function(X, train_ratio, r){
  N <- nrow(X)
  Times <- ncol(X)
  err <- rep(NA, 3)
  ########## Missing at random ##############
  train_mat <- matrix(runif(N * Times) < train_ratio, nrow = N, ncol = Times)
  X_train <- X*train_mat
  X_test <- X*(1-train_mat)
  ########### PCA ##############
  svd_res <- svd(X_train, nu = r, nv = r)
  if (r>1){
    X_PCA <- svd_res$u %*%  diag(svd_res$d[1:r]) %*% t(svd_res$v)/train_ratio
  } else {
    X_PCA <- svd_res$u %*% t(svd_res$v) * svd_res$d[r] /train_ratio
  }
  ########### HeteroPCA ##############
  U_heteroPCA <- heteroPCA(X_train, r, T0 = 10)
  ########### PCA ##############
  res_APCA <- CV_APCA(X_train, r, p_star = 0.8, J = 10)
  err[1] <- sum((X_test - (res_APCA$U %*% t(res_APCA$U) %*% X_train)/train_ratio * (1-train_mat))^2)/sum(X_test^2)
  err[2]<- sum((X_test - X_PCA * (1-train_mat))^2)/sum(X_test^2)
  err[3] <- sum((X_test - (U_heteroPCA %*%  t(U_heteroPCA) %*% X_train)/train_ratio * (1-train_mat))^2)/sum(X_test^2)
  return(err)
}

n_rep <- 100
r <- 3
train_ratio <- 0.9
errs <- replicate(n_rep, one_sim(X, train_ratio, r))
save(errs, file = "error_qd_r3_tr_0_9.RData")
rowMeans(errs)
apply(errs, 1, sd)










