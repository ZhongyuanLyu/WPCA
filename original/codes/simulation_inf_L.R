library(ggplot2)
library(dplyr)
library(tidyr)
library(expm)
source("functions.R")
set.seed(123)
N = 100
Times <- 800  
rho = 0.6
r <- 3;  
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
############ General Covarinace #############
sigeC <- diag(runif(N, 1,20)) + rho*(rep(1,N)%*%t(rep(1,N))-diag(N))

########### Compute Covaraince ###########
gamma <- 0
Q <- diag(gamma, nrow = Times)
diag(Q[-nrow(Q),-1]) = (1-gamma)
Q[lower.tri(Q)]  <- t(Q)[lower.tri(Q)]

MQMt <- M %*% Q %*% t(M)
ind <- 0.5 * N

############ Simulation ############        
N_rep <- 500
ind <-  0.5*N
LResone <- matrix(NA, nrow = r, ncol = N_rep)
LResnorm <- matrix(NA, nrow = r, ncol = N_rep)
for (ii in 1:N_rep) {
  E<- t(mvrnorm(Times, rep(0,N), sigeC))
  X <- M+E
  U_WPCA <- APCA(X, 0, r, symmetric = TRUE)
  svd_nUUtX <- svd(U_WPCA %*% t(U_WPCA) %*% X / sqrt(Times), nu = r, nv = r)
  rotation_cov <- compute_rotation_and_covariance_L(svd_nUUtX$u, svd_nUUtX$v, Fmat, L, M/sqrt(Times), 
                                                    MQMt, sigeC, diag(Times), Q, ind)
  LRes <- svd_nUUtX$u %*% diag(svd_nUUtX$d[1:r]) - L %*% rotation_cov$R_LQ
  LResnorm[, ii] <-  LRes[ind,] %*% solve(sqrtm(rotation_cov$Sigma_L))
}


jj <- 3
ggplot(data.frame(x = LResnorm[jj,]), aes(sample = LResnorm[jj,])) +
  stat_qq(size = 2) +
  stat_qq_line(color = "steelblue", size = 1) +
  labs(title = "Normal Q–Q Plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal(base_size = 30) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),  # panel border
    panel.background = element_blank()
  )



df <- data.frame(EstimationError = LResnorm[jj, ])
ggplot(df, aes(x = EstimationError)) +
  geom_histogram(aes(y = ..density..),
                 bins  = 12,
                 fill  = "royalblue3",
                 color = "white") +
  stat_function(fun = dnorm,
                args = list(mean = 0, sd = 1),
                color = "black",
                size  = 1) +
  labs(x     = "Estimation Error",
       y     = "Frequency",
       title = "Histogram of Loading Errors") +
  theme_minimal(base_size = 30)  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),  # panel border
    panel.background = element_blank()
  )


