library(matrixNormal)
library(MASS)


fredqd <- function(file = "", date_start = NULL, date_end = NULL, transform = TRUE) {
  # Error checking
  if (!is.logical(transform))
    stop("'transform' must be logical.")
  if ((class(date_start) != "Date") && (!is.null(date_start)))
    stop("'date_start' must be Date or NULL.")
  if ((class(date_end) != "Date") && (!is.null(date_end)))
    stop("'date_end' must be Date or NULL.")
  
  if (class(date_start) == "Date") {
    if (as.numeric(format(date_start, "%d")) != 1)
      stop("'date_start' must be Date whose day is 1.")
    if (!as.numeric(format(date_start, "%m")) %in% c(3,6,9,12))
      stop("'date_start' must be Date whose month is March, June,
           September, or December.")
    if (date_start < as.Date("1959-03-01"))
      stop("'date_start' must be later than 1959-03-01.")
  }
  
  if (class(date_end) == "Date") {
    if (as.numeric(format(date_end, "%d")) != 1)
      stop("'date_end' must be Date whose day is 1.")
    if (!as.numeric(format(date_end, "%m")) %in% c(3,6,9,12))
      stop("'date_end' must be Date whose month is March, June,
           September, or December.")
  }
  
  
  
  # Prepare raw data
  rawdata <- readr::read_csv(file, col_names = FALSE, col_types = cols(X1 = col_date(format = "%m/%d/%Y")),
                             skip = 3)
  
  rawdata <- as.data.frame(rawdata)
  row_to_remove = c()
  for (row in (nrow(rawdata)-20):nrow(rawdata)){
    if(!any(is.finite(unlist(rawdata[row, ])))){
      row_to_remove = c(row_to_remove,row)# remove NA rows
    }
  }
  if(length(row_to_remove)>0){
    rawdata = rawdata[-row_to_remove,]
  }
  
  attrdata <- utils::read.csv(file, header = FALSE, nrows = 3)
  header <- c("date", unlist(attrdata[1,2:ncol(attrdata)]))
  colnames(rawdata) <- header
  
  
  # Import tcode tcodes is an internal data of the R package
  tcode <- unlist(attrdata[3,2:ncol(attrdata)])
  
  
  # Subfunction transxf: data transformation based on tcodes
  transxf <- function(x, tcode) {
    # Number of observations (including missing values)
    n <- length(x)
    
    # Value close to zero
    small <- 1e-06
    
    # Allocate output variable
    y <- rep(NA, n)
    y1 <- rep(NA, n)
    
    # TRANSFORMATION: Determine case 1-7 by transformation code
    if (tcode == 1) {
      # Case 1 Level (i.e. no transformation): x(t)
      y <- x
      
    } else if (tcode == 2) {
      # Case 2 First difference: x(t)-x(t-1)
      y[2:n] <- x[2:n] - x[1:(n - 1)]
      
    } else if (tcode == 3) {
      # case 3 Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
      y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
      
    } else if (tcode == 4) {
      # case 4 Natural log: ln(x)
      if (min(x, na.rm = TRUE) > small)
        y <- log(x)
      
    } else if (tcode == 5) {
      # case 5 First difference of natural log: ln(x)-ln(x-1)
      if (min(x, na.rm = TRUE) > small) {
        x <- log(x)
        y[2:n] <- x[2:n] - x[1:(n - 1)]
      }
      
    } else if (tcode == 6) {
      # case 6 Second difference of natural log:
      # (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
      if (min(x, na.rm = TRUE) > small) {
        x <- log(x)
        y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
      }
      
    } else if (tcode == 7) {
      # case 7 First difference of percent change:
      # (x(t)/x(t-1)-1)-(x(t-1)/x(t-2)-1)
      y1[2:n] <- (x[2:n] - x[1:(n - 1)])/x[1:(n - 1)]
      y[3:n] <- y1[3:n] - y1[2:(n - 1)]
    }
    
    return(y)
  }
  
  
  # Transform data
  if (transform) {
    # Apply transformations
    N <- ncol(rawdata)
    data <- rawdata
    data[, 2:N] <- NA
    
    # Perform transformation using subfunction transxf (see below for
    # details)
    for (i in 2:N) {
      temp <- transxf(rawdata[, i], tcode[i - 1])
      data[, i] <- temp
    }
    
  } else {
    data <- rawdata
  }
  
  
  # Null case of date_start and date_end
  if (is.null(date_start))
    date_start <- as.Date("1959-03-01")
  if (is.null(date_end))
    date_end <- data[, 1][nrow(data)]
  
  
  # Subset data
  index_start <- which.max(data[, 1] == date_start)
  index_end <- which.max(data[, 1] == date_end)
  
  outdata <- data[index_start:index_end, ]
  class(outdata) <- c("data.frame", "fredqd")
  return(outdata)
  
}



procrustes_align <- function(A, B) {
  sv <- svd(t(A) %*% B)
  R  <- sv$v %*% t(sv$u)
}


compute_rotation_and_covariance_L <- function(U_hat, V_hat, F_true, L_true, Mnorm, MQMt, Sigma_C, Sigma_T, Q, ind = NULL) {
  N   <- nrow(L_true)
  r   <- ncol(L_true)
  TT   <- ncol(Sigma_T)
  
  svd_MQMt  <- svd(MQMt, nu=r, nv=r)
  U_bar <- svd_MQMt$u
  
  svd_Mnorm <- svd(Mnorm, nu=r, nv=r)
  U_true <- svd_Mnorm$u            
  V_true <- svd_Mnorm$v            
  
  RV <- procrustes_align(V_hat, V_true)  
  O_bar <- procrustes_align(U_true, U_bar)
  
  O_F <- t(O_bar) %*% diag(1/svd_MQMt$d[1:r]) %*% O_bar %*% diag(svd_Mnorm$d[1:r]) %*% t(RV)
    
  # Form the "alignment" matrix for U
  RU <- procrustes_align(U_hat, U_true)  
  
  # Rotation matrix for loadings:
  B <- sqrt(TT) * solve(crossprod(F_true)) %*%  t(F_true) %*% V_true
  R_LQ <- t(solve(B %*% t(RV)))        
  # R_LQ <- t(solve(B)) %*% t(RV) 
  
  out_term <- t(O_F) %*% diag(svd_Mnorm$d[1:r]) %*% t(V_true) %*% Q
  
  if (!is.null(ind)) {
    Sigma_L <- TT * Sigma_C[ind,ind] * out_term %*% Sigma_T %*% t(out_term)
  } else {
    Sigma_L   <- vector("list", N)
    for (i in seq_len(N)) {
      Sigma_L[[i]] <- TT * Sigma_C[i,i] * out_term %*% Sigma_T %*% t(out_term)
    }
  }
  

  
  list(R_LQ = R_LQ, Sigma_L = Sigma_L)
}


compute_rotation_and_covariance_F <- function(V_hat, F_true, M, Sigma_C, Sigma_T, ind = NULL) {
  TT   <- nrow(F_true)
  r   <- ncol(F_true)
  N   <- nrow(M)
  
  svd_M  <- svd(M, nu=r, nv=r)
  U_true <- svd_M$u            # N×r
  Sigma_true <- diag(svd_M$d[1:r])  # r×r
  V_true <- svd_M$v            # T×r
  
  # 4. Form the "alignment" matrix for V
  RV <- procrustes_align(V_hat, V_true)  #
  
  # 5. Compute the linking matrix B from V_true and F_true:
  B <- sqrt(TT) * solve(crossprod(F_true)) %*%  t(F_true) %*% V_true
  
  # 6. Rotation matrix for factors:
  R_FQ <- B %*% t(RV)         # r×r  
  
  # 7. The theoretical covariance for each time t, from Theorem 9:
  #    Σ_{F,t} = [Σ_T]_{t,t} · RV %*% Σ_true^{-1} U_trueᵀ Σ_C U_true Σ_true^{-1} %*% t(RV)
  Sigma_inv   <- diag(1/svd_M$d[1:r])
  middle_term <- Sigma_inv %*% t(U_true) %*% Sigma_C %*% U_true %*% Sigma_inv
  
  if (!is.null(ind)) {
    Sigma_F <- Sigma_T[ind,ind] * (RV %*% middle_term %*% t(RV))
  } else {
    Sigma_F   <- vector("list", TT)
    for (t in seq_len(TT)) {
      Sigma_F[[t]] <- Sigma_T[t,t] * (RV %*% middle_term %*% t(RV))
    }
  }
  
  
  # covariance matrices Σ_{F,t} for t=1,…,T 
  
  list(R_FQ = R_FQ, Sigma_F = Sigma_F)
}



generateKFolds<- function(data, k) {
  n <- nrow(data)
  d <- ncol(data)
  total <- n * d
  assignments <- sample(rep(1:k, length.out = total))
  folds <- vector("list", k)
  for (i in 1:k) {
    indicator <- as.integer(assignments == i)
    folds[[i]] <- matrix(indicator, nrow = n, ncol = d)
  }
  return(folds)
}

H_mat <- function(X) {
  diag(X) <- 0
  return(X)
}

heteroPCA <- function(R, K, T0 = 10) {
  M <- H_mat(R %*% t(R))
  M_no_diag <- M
  if (T0 > 0) {
    for (t in 1:T0) {
      svd_res <- svd(M, nu = K, nv = K)
      if (K == 1) {
        M_bar <- svd_res$u %*% t(svd_res$v) * svd_res$d
      } else {
        M_bar <- svd_res$u %*% diag(svd_res$d[1:K]) %*% t(svd_res$v)
      }
      M <- M_no_diag + diag(diag(M_bar))
    }
  }
  U_hat <- svd(M, nu = K, nv = K)$u
  return(U_hat)
}



VAR1.sim <- function(B0,Phi,n,sigvare)
{
  nob <- 50+n
  k = nrow(B0)
  x <- matrix(0,nrow=nob,ncol=k)
  tem <- chol(sigvare)
  for (i in 2:nob)
    x[i,] <-  B0 + Phi %*% matrix(x[i-1,]) + crossprod(tem, matrix(rnorm(k),nrow=k))
  return(x[51:nob,])
}




APCA <- function(X, gamma, r, symmetric = TRUE){
  Times <- ncol(X)
  Q <- diag(gamma, nrow = Times)
  if (symmetric == TRUE){
    diag(Q[-nrow(Q),-1]) = (1-gamma)
    Q[lower.tri(Q)]  <- t(Q)[lower.tri(Q)]
    Uhat <- svd(X%*%Q%*%t(X),nu=r)$u
  } else{
    diag(Q[-nrow(Q),-1]) = (1-gamma)
    Uhat <- svd(X%*%Q%*%t(X),nu=r)$u
  }
  return(Uhat)
}



create_folds <- function(data, k) {
  nC <- ncol(data)
  folds <- split(seq_len(nC), rep(1:k, length.out = nC))
  return(folds)
}


###### CV #####
CV_APCA <- function(X, r, p_star = 0.8, grid_len  = 10,  J = 5){
  N <- nrow(X)
  Times <- ncol(X)
  gamma_list = seq(0, 1, length.out = grid_len)
  err_mat = matrix(0, nrow = J, ncol = grid_len)
  for (jj in 1:J){
    mask <- matrix(runif(N * Times) < p_star, nrow = N, ncol = Times)
    Y_train <- X*mask
    Y_test <- X*(1-mask)
    for (j in 1:length(gamma_list)){
      Ugammahat_j <- APCA(Y_train, gamma_list[j], r)
      err_mat[jj, j] <- sum((Y_test - (Ugammahat_j %*% t(Ugammahat_j) %*% Y_train / p_star) * (1-mask))^2)/sum(1-mask)
    }
  }
  err_vec <- colMeans(err_mat)
  best_gamma <-  gamma_list[which.min(err_vec)]
  U_best <- APCA(X, best_gamma, r)
  return(list(U = U_best, gamma =  best_gamma, err = err_vec))
}

ratio_based_r <- function(X, R, alpha = 0){
  Times <- ncol(X)
  Q <- diag(0, nrow = Times)
  diag(Q[-nrow(Q),-1]) = 1
  svs <- svd(X%*%t(X) + alpha * X%*%Q%*%t(X))$d
  gap_list <- rep(NA, R-1)
  for (j in 1:(R-1)){
    gap_list[j] <- svs[j+1]/svs[j]
  }
  return(sort(gap_list, index.return = TRUE)[[2]][1:3])
}
