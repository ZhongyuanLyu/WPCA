library(ggplot2)
library(dplyr)
library(tidyr)
source("functions.R")
set.seed(1234)
sim_Times <- function(N = 100, Times,  rho, N_rep = 1){
  r <- 3;  
  err_res <- matrix(NA, nrow = N_rep, ncol = 6)
  for (ii in 1:N_rep) {
    ############ Generate F #############
    B0=matrix(0,nrow=r)
    diagPhivec <- seq(0.9,0.9, length.out=r)
    # Phi <- diag(diagPhivec)
    U <- svd(matrix(rnorm(r*r),nrow=r))$u
    V <- svd(matrix(rnorm(r*r),nrow=r))$u
    Phi <- U%*%diag(diagPhivec)%*%t(V)
    sigvare=diag(1,r,r)-Phi%*%t(Phi)
    Fmat=VAR1.sim(B0,Phi,Times,sigvare)
    Fmat=ts(Fmat)
    ############ Generate L #############
    svd_L <- svd(matrix(rnorm(N*r),nrow=N))
    L <- svd_L$u %*% diag(svd_L$d)
    ############ Generate M #############
    M <- L%*%t(Fmat) 
    U_true <- svd(M,nu=r)$u
    V_true <- svd(M,nv=r)$v
    
    ############ Isotropic Covariance #############
    E<- t(mvrnorm(Times, rep(0,N), 10*diag(N)))
    
    ############ Diagonal Covarinace #############
    # E<- t(mvrnorm(Times, rep(0,N), diag(runif(N, 1,20))))
    
    
    ############ General Covarinace #############
    # sigeC <- diag(runif(N, 1,20)) + rho*(rep(1,N)%*%t(rep(1,N))-diag(N))
    # E<- t(mvrnorm(Times, rep(0,N), sigeC))
    ############ Generate Data ############ 
    X <- M+E
    
    U_PCA <- svd(X, nu = r, nv = r)$u
    U_HPCA <- heteroPCA(X, r, 10)
    U_APCA <- CV_APCA(X, r, p_star = 0.8, J = 10)$U
    # U_APCA <- APCA(X, 0, r)
    
    err_res[ii,1] <- sqrt(sum((U_PCA%*%t(U_PCA)-U_true%*%t(U_true))^2)/sqrt(r))
    err_res[ii,2] <- sqrt(sum((U_HPCA%*%t(U_HPCA)-U_true%*%t(U_true))^2)/sqrt(r))
    err_res[ii,3] <- sqrt(sum((U_APCA%*%t(U_APCA)-U_true%*%t(U_true))^2)/sqrt(r))
    
    V_PCA <- svd(U_PCA %*% t(U_PCA) %*% X, nu = r, nv = r)$v
    V_HPCA <- svd(U_HPCA %*% t(U_HPCA) %*% X, nu = r, nv = r)$v
    V_APCA <- svd(U_APCA %*% t(U_APCA) %*% X, nu = r, nv = r)$v
    
    
    err_res[ii,4] <- sqrt(sum((V_PCA%*%t(V_PCA)-V_true%*%t(V_true))^2)/sqrt(r))
    err_res[ii,5] <- sqrt(sum((V_HPCA%*%t(V_HPCA)-V_true%*%t(V_true))^2)/sqrt(r))
    err_res[ii,6] <- sqrt(sum((V_APCA%*%t(V_APCA)-V_true%*%t(V_true))^2)/sqrt(r))
    
  }
  return(err_res)
}


################ Plot with T #####################
# num_par <- 2
# Times_list <- seq(500,1000, length.out = num_par)
# res_list <- vector("list", length = num_par)
# for (i in 1:num_par) {
#   res_list[[i]] <- sim_Times(Times_list[[i]], rho = 0.6, N_rep = 5)
# }
# save(res_list, file = "sim_T_rho0_6.Rdata")

# means_list <- lapply(res_list, colMeans)
# means_df <- do.call(rbind, means_list) %>%      
#   as.data.frame() %>%
#   mutate(Times = Times_list)         
# colnames(means_df)[1:3] <- c("PCA", "HeteroPCA", "APCA")
# means_long <- means_df %>%
#   pivot_longer(
#     cols      = -Times,
#     names_to  = "variable",
#     values_to = "mean"
#   )
# 
# ggplot(means_long, aes(x = Times, y = mean, color = variable)) +
#   geom_line(size = 1) +
#   geom_point(size = 5) +
#   labs(
#     x     = "T",
#     y     = "Subspace Estimation Error",
#     color = "Method"
#   ) +
#   theme_minimal() +
#     theme(
#       panel.border = element_rect(color = "black", fill = NA, size = 2),
#       legend.title = element_blank(),
#       legend.text = element_text(size=40),
#       legend.position = "right",
#       legend.background = element_rect(fill='transparent'),
#       axis.text.x = element_text(size = 35, colour = "black"),
#       axis.text.y = element_text(size = 40),
#       axis.title.x = element_blank(),
#       axis.title.y = element_text(size = 40)
#     )



################ Plot with Methods #####################
# N <- 100
# Times <- 500
# err_res <- sim_Times(Times, rho = 0.6, N_rep = 100)
# cat("(N,T)=(",N,",",Times,"):",colMeans(err_res))
# apply(err_res, 2, sd)


# 
# methods <- c("PCA", "HeteroPCA", "APCA")
# err_df <- data.frame(
#   Method = rep(methods, each = N_rep),  # Repeat each K value for all repetitions
#   Error = as.vector(err_res)                  # Flatten the running times
# )
# 
# ggplot(err_df, aes(x = factor(Method), y = Error)) +
#   geom_boxplot(fill = "skyblue", color = "darkblue") +
#   ylab("Estimation error") +
#   theme_minimal() +
#   theme(
#     panel.border = element_rect(color = "black", fill = NA, size = 2),
#     # plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
#     axis.text.x = element_text(size = 35, colour = "black"),
#     axis.text.y = element_text(size = 40),
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(size = 40)
#   )


N = 100
num_par <- 4
Times_list <- c(250, 500, 750, 1000)
res_list <- vector("list", length = num_par)
for (i in 1:num_par) {
  res_list[[i]] <- sim_Times(N, Times_list[[i]], rho = 0.6, N_rep = 100)
  cat("(N,T)=(",N,",",Times_list[[i]],") mean:",colMeans(res_list[[i]]),"\n")
  cat("(N,T)=(",N,",",Times_list[[i]],") std:",apply(res_list[[i]], 2, sd),"\n")
}
save(res_list, file = "sim_est_random_N100_isotropic.Rdata")


for (i in 1:num_par) {
  cat("(N,T)=(",N,",",Times_list[[i]],") mean:",round(colMeans(res_list[[i]]),3),"\n")
  cat("(N,T)=(",N,",",Times_list[[i]],") std:",round(apply(res_list[[i]], 2, sd),3),"\n")
}
