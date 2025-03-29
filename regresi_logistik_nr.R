# TUGAS PRAKTIKUM W4 - Komputasi Sains Data
# Nama: Dimas Wahyu Saputra
# NRP : 5052231022

# 1. Regresi Logistik dengan metode Newton-Rhapson
log_likelihood <- function(beta, X, y){
  p <- 1/(1+exp(-X%*%beta))
  ll <- sum(y*log(p)+(1-y)*log(1-p))
  return(ll)
}

gradient <- function(beta,X,y){
  p <- 1/(1+exp(-X%*%beta))
  grad <- t(X)%*%(y-p)
  return(grad)
}

hessian <- function(beta, X){
  p <- 1/(1+exp(-X%*%beta))
  W <- diag(as.vector(p*(1-p)))
  H <- -t(X)%*%W%*%X
  return(H)
}

# NR Algorithm
nr <- function(X, y, tol = 1e-5, max.iter=1000){
  beta <- rep(0, ncol(X))
  for(i in 1:max.iter){
    grad <- gradient(beta, X, y)
    H <- hessian(beta, X)
    beta_new <- beta - solve(H) %*% grad
    if(sqrt(sum((beta_new - beta)^2)) < tol){
      cat('Convergence reached in', i, 'Iterations\n')
      return(list(method = "Newton-Raphson", beta = beta_new, fit = 1 / (1 + exp(-X %*% beta_new))))
    }
    beta <- beta_new
  }
  warning('Maximum Iteration Reached Without Convergence')
  return(list(method = "Newton-Raphson", beta = beta, fit = 1 / (1 + exp(-X %*% beta))))
}
