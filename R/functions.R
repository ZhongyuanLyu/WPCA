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
  rawdata <- readr::read_csv(file, col_names = FALSE, col_types = readr::cols(X1 = readr::col_date(format = "%m/%d/%Y")),
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
  R  <- sv$u %*% t(sv$v)
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
  O_bar <- procrustes_align(U_bar, U_true)
  
  Sigma_bar_inv <- diag(
    1 / svd_MQMt$d[seq_len(r)], nrow = r, ncol = r
  )
  Sigma_true <- diag(
    svd_Mnorm$d[seq_len(r)], nrow = r, ncol = r
  )
  O_F <- t(O_bar) %*% Sigma_bar_inv %*% O_bar %*% Sigma_true %*% t(RV)
    
  # Form the "alignment" matrix for U
  RU <- procrustes_align(U_hat, U_true)  
  
  # Rotation matrix for loadings:
  B <- sqrt(TT) * solve(crossprod(F_true)) %*%  t(F_true) %*% V_true
  R_LQ <- t(solve(B %*% t(RV)))        
  # R_LQ <- t(solve(B)) %*% t(RV) 
  
  C_L_base <- TT * Sigma_true %*% t(V_true) %*% Q %*% Sigma_T %*% Q %*% V_true %*% Sigma_true
  
  if (!is.null(ind)) {
    C_L <- Sigma_C[ind,ind] * C_L_base
    Sigma_L <- t(O_F) %*% C_L %*% O_F
    W_L <- symmetric_matrix_invsqrt(C_L) %*% solve(t(O_F))
  } else {
    Sigma_L   <- vector("list", N)
    W_L <- vector("list", N)
    for (i in seq_len(N)) {
      C_L <- Sigma_C[i,i] * C_L_base
      Sigma_L[[i]] <- t(O_F) %*% C_L %*% O_F
      W_L[[i]] <- symmetric_matrix_invsqrt(C_L) %*% solve(t(O_F))
    }
  }
  

  
  list(R_LQ = R_LQ, Sigma_L = Sigma_L, W_L = W_L)
}


compute_rotation_and_covariance_F <- function(V_hat, F_true, M, Sigma_C, Sigma_T, ind = NULL) {
  TT   <- nrow(F_true)
  r   <- ncol(F_true)
  N   <- nrow(M)
  
  svd_M  <- svd(M, nu=r, nv=r)
  U_true <- svd_M$u            # N×r
  Sigma_true <- diag(svd_M$d[seq_len(r)], nrow = r, ncol = r)  # r×r
  V_true <- svd_M$v            # T×r
  
  # 4. Form the "alignment" matrix for V
  RV <- procrustes_align(V_hat, V_true)  #
  
  # 5. Compute the linking matrix B from V_true and F_true:
  B <- sqrt(TT) * solve(crossprod(F_true)) %*%  t(F_true) %*% V_true
  
  # 6. Rotation matrix for factors:
  R_FQ <- B %*% t(RV)         # r×r  
  
  # 7. The theoretical covariance for each time t, from Theorem 9:
  #    Σ_{F,t} = [Σ_T]_{t,t} · RV %*% Σ_true^{-1} U_trueᵀ Σ_C U_true Σ_true^{-1} %*% t(RV)
  Sigma_inv <- diag(1 / svd_M$d[seq_len(r)], nrow = r, ncol = r)
  C_F_base <- t(U_true) %*% Sigma_C %*% U_true
  D_F <- RV %*% Sigma_inv
  
  if (!is.null(ind)) {
    C_F <- Sigma_T[ind,ind] * C_F_base
    Sigma_F <- D_F %*% C_F %*% t(D_F)
    W_F <- symmetric_matrix_invsqrt(C_F) %*% solve(D_F)
  } else {
    Sigma_F   <- vector("list", TT)
    W_F <- vector("list", TT)
    for (t in seq_len(TT)) {
      C_F <- Sigma_T[t,t] * C_F_base
      Sigma_F[[t]] <- D_F %*% C_F %*% t(D_F)
      W_F[[t]] <- symmetric_matrix_invsqrt(C_F) %*% solve(D_F)
    }
  }
  
  
  # covariance matrices Σ_{F,t} for t=1,…,T 
  
  list(R_FQ = R_FQ, Sigma_F = Sigma_F, W_F = W_F)
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
        M_bar <- svd_res$u %*% t(svd_res$v) * svd_res$d[1L]
      } else {
        M_bar <- svd_res$u %*%
          diag(svd_res$d[seq_len(K)], nrow = K, ncol = K) %*%
          t(svd_res$v)
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
  return(x[51:nob, , drop = FALSE])
}




APCA <- function(X, gamma, r, symmetric = TRUE){
  Times <- ncol(X)
  Q <- diag(gamma, nrow = Times)
  if (symmetric == TRUE){
    if (Times > 1L) Q[cbind(seq_len(Times - 1L), seq.int(2L, Times))] = (1-gamma)
    Q[lower.tri(Q)]  <- t(Q)[lower.tri(Q)]
    XQXt <- X%*%Q%*%t(X)
    Uhat <- eigen((XQXt+t(XQXt))/2,symmetric=TRUE)$vectors[,seq_len(r),drop=FALSE]
  } else{
    if (Times > 1L) Q[cbind(seq_len(Times - 1L), seq.int(2L, Times))] = (1-gamma)
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
# base_observation_probability is the retention probability that produced X
# before the current CV split (one for a complete input matrix).
CV_APCA <- function(X, r, p_star = 0.8, grid_len = 10, J = 5,
                    observed_mask = NULL,
                    base_observation_probability = 1) {
  X <- as.matrix(X)
  N <- nrow(X)
  Times <- ncol(X)

  if (N < 1L || Times < 1L) {
    stop("X must have at least one row and one column.")
  }
  if (!is.numeric(X)) {
    stop("X must be numeric.")
  }
  if (length(r) != 1L || !is.numeric(r) || !is.finite(r) ||
      r < 1L || r > min(N, Times) || r != floor(r)) {
    stop("r must be a positive integer no larger than min(dim(X)).")
  }
  if (length(p_star) != 1L || !is.numeric(p_star) || !is.finite(p_star) ||
      p_star <= 0 || p_star >= 1) {
    stop("p_star must be a finite scalar strictly between 0 and 1.")
  }
  if (length(base_observation_probability) != 1L ||
      !is.numeric(base_observation_probability) ||
      !is.finite(base_observation_probability) ||
      base_observation_probability <= 0 || base_observation_probability > 1) {
    stop("base_observation_probability must be a finite scalar in (0, 1].")
  }
  if (length(grid_len) != 1L || !is.numeric(grid_len) || !is.finite(grid_len) ||
      grid_len < 1L || grid_len > .Machine$integer.max ||
      grid_len != floor(grid_len)) {
    stop("grid_len must be a positive integer.")
  }
  if (length(J) != 1L || !is.numeric(J) || !is.finite(J) ||
      J < 1L || J > .Machine$integer.max || J != floor(J)) {
    stop("J must be a positive integer.")
  }

  if (is.null(observed_mask)) {
    observed_mask <- matrix(TRUE, nrow = N, ncol = Times)
  } else {
    if (!is.matrix(observed_mask) || !identical(dim(observed_mask), dim(X))) {
      stop("observed_mask must be a matrix with the same dimensions as X.")
    }
    if (is.logical(observed_mask)) {
      if (anyNA(observed_mask)) {
        stop("observed_mask must not contain missing values.")
      }
    } else if (is.numeric(observed_mask)) {
      if (any(!is.finite(observed_mask)) ||
          any(!(observed_mask %in% c(0, 1)))) {
        stop("A numeric observed_mask must contain only finite zeros and ones.")
      }
      observed_mask <- observed_mask != 0
    } else {
      stop("observed_mask must be a logical matrix or a numeric zero-one matrix.")
    }
  }

  n_observed <- sum(observed_mask)
  if (n_observed < 2L) {
    stop("observed_mask must contain at least two observed entries for cross-validation.")
  }
  if (any(!is.finite(X[observed_mask]))) {
    stop("X must be finite at every entry selected by observed_mask.")
  }

  # Values outside observed_mask are unavailable to both model fitting and scoring.
  X_observed <- X
  X_observed[!observed_mask] <- 0

  r <- as.integer(r)
  grid_len <- as.integer(grid_len)
  J <- as.integer(J)
  gamma_list <- seq(0, 1, length.out = grid_len)
  err_mat <- matrix(NA_real_, nrow = J, ncol = grid_len)
  for (jj in seq_len(J)) {
    attempts <- 0L
    repeat {
      attempts <- attempts + 1L
      split_mask <- matrix(runif(N * Times) < p_star,
                           nrow = N, ncol = Times)
      inner_train_mask <- observed_mask & split_mask
      inner_valid_mask <- observed_mask & !split_mask
      if (any(inner_train_mask) && any(inner_valid_mask)) {
        break
      }
      if (attempts >= 1000L) {
        stop("Could not construct non-empty inner training and validation sets.")
      }
    }

    Y_train <- X_observed * inner_train_mask
    n_valid <- sum(inner_valid_mask)
    for (j in seq_along(gamma_list)) {
      Ugammahat_j <- APCA(Y_train, gamma_list[j], r)
      prediction <- Ugammahat_j %*% t(Ugammahat_j) %*% Y_train /
        (base_observation_probability * p_star)
      err_mat[jj, j] <- sum((X_observed[inner_valid_mask] -
                              prediction[inner_valid_mask])^2) / n_valid
    }
  }
  err_vec <- colMeans(err_mat)
  best_gamma <- gamma_list[which.min(err_vec)]
  U_best <- APCA(X_observed, best_gamma, r)
  return(list(U = U_best, gamma = best_gamma, err = err_vec))
}

ratio_based_r <- function(X, R) {
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X) < 2L || ncol(X) < 1L ||
      any(!is.finite(X))) {
    stop("X must be a finite numeric matrix with at least two rows and one column.")
  }
  if (length(R) != 1L || !is.numeric(R) || !is.finite(R) ||
      R < 1L || R > .Machine$integer.max || R != floor(R)) {
    stop("R must be a positive integer.")
  }

  R <- as.integer(R)
  if (R >= nrow(X)) {
    stop("R must be smaller than nrow(X) so that sigma_{R+1} is available.")
  }

  singular_values <- svd(tcrossprod(X), nu = 0L, nv = 0L)$d
  denominators <- singular_values[seq_len(R)]
  if (any(!is.finite(denominators)) || any(denominators <= 0)) {
    stop("The first R singular values of X X' must be finite and positive.")
  }

  ratios <- singular_values[seq.int(2L, R + 1L)] / denominators
  as.integer(which.min(ratios))
}
