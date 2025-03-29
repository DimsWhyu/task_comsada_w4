# TUGAS PRAKTIKUM W4 - Komputasi Sains Data
# Nama: Dimas Wahyu Saputra
# NRP : 5052231022

# 3. Menggabungkan metode Newton-Rhapson dan IRLS
logistic_regression <- function(X, y, method = "NR", tol = 1e-5, max.iter = 1000){
  beta <- rep(0, ncol(X))

  for(i in 1:max.iter){
    p <- 1 / (1 + exp(-X %*% beta))
    W <- diag(as.vector(p * (1 - p)))

    if (method == "NR") {
      grad <- t(X) %*% (y - p)
      H <- -t(X) %*% W %*% X
      beta_new <- beta - solve(H) %*% grad
    } else if (method == "IRLS") {
      z <- X %*% beta + solve(W) %*% (y - p)
      xtw <- t(X) %*% W
      xtwx_inv <- solve(t(X) %*% W %*% X)
      beta_new <- xtwx_inv %*% (xtw %*% z)
    } else {
      stop("Method must be 'NR' or 'IRLS'")
    }

    if(sqrt(sum((beta_new - beta)^2)) < tol){
      cat('Converged in', i, 'iterations using', method, '\n')
      return(list(method = method, beta = beta_new, fit = 1 / (1 + exp(-X %*% beta_new))))
    }
    beta <- beta_new
  }

  cat('Reached Max Iteration without full convergence using', method, '\n')
  return(list(method = method, beta = beta, fit = 1 / (1 + exp(-X %*% beta))))
}

