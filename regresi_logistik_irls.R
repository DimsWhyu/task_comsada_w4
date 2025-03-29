# TUGAS PRAKTIKUM W4 - Komputasi Sains Data
# Nama: Dimas Wahyu Saputra
# NRP : 5052231022

# 2. Regresi Logistik dengan metode IRLS
# IRLS Algorithm
irls <- function(X, y, tol=1e-5, max.iter=1000){
  beta <- rep(0, ncol(X))
  for(i in 1:max.iter){
    p <- 1/(1+exp(-X %*% beta))
    W <- diag(as.vector(p * (1 - p)))
    z <- X %*% beta + solve(W) %*% (y - p)
    xtw <- t(X) %*% W
    xtwx_inv <- solve(t(X) %*% W %*% X)
    beta_new <- xtwx_inv %*% (xtw %*% z)

    if(sqrt(sum((beta_new - beta)^2)) < tol){
      cat('Converged in', i, 'iterations\n')
      return(list(method = "IRLS", beta = beta_new, fit = 1 / (1 + exp(-X %*% beta_new))))
    }
    beta <- beta_new
  }
  cat('Reached Max Iteration without full convergence\n')
  return(list(method = "IRLS", beta = beta, fit = 1 / (1 + exp(-X %*% beta))))
}
