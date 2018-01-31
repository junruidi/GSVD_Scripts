# # all functions for 4th order tensor OD
# 
# Xvv = function(X,v){
#   dimX = dim(X)[1]
#   
#   z = NULL
#   for(i in 1:dimX){
#     aa = 0
#     for(j in 1:dimX){
#       for(l in 1:dimX){
#         for(k in 1:dimX){
#           aa = aa + X[j,l,k,i] * v[j] * v[l] * v[k]
#         }
#       }
#     }
#     z = c(z,aa)
#   }
#   return(z)
# }
# 
# Xvvv = function(X,v){
#   lm  = dim(X)[1]
#   zz = 0
#   for(j4 in 1:lm){
#    for(j3 in 1:lm){
#      for(j2 in 1:lm){
#        for(j1 in 1:lm){
#          zz = zz +  X[j1, j2, j3,j4] * v[j1] * v[j2] * v[j3] * v[j4]
#         }
#       }
#     }
#   }
#   return(zz)
# }
library(rTensor)
Tensor_vector = function(A,v,m){
  tnsr = as.tensor(A)
  V = matrix(v,ncol = 1)
  #m = length(dim(A))
  c.t = tnsr
  for(i in 1:m){
    c.t = ttm(c.t,t(V),m=i)
  }
  result = c(c.t@data)
}

power_itr = function(theta, X,N){
  for(i in 1:N){
    #next_itr = Xvv(X, v = theta)
    next_itr = Tensor_vector(A = X, v = theta, m = 3)
    theta = next_itr/sqrt(sum(next_itr^2))
  }
  return(theta)
}

#select the best theta, lambda, and deflate the original tensor
Est_PI = function(X, L=10, N=10){
  p = dim(X)[1]
  
  theta_list = list()
  for ( t in 1:L){
    vc = rnorm(p)
    theta_0 = vc/sqrt(sum(vc^2))
    theta_list[[t]] = power_itr(theta = theta_0, X=X, N=N)
  }
  
  lambda_list = NULL
  for(t in 1:L){
    #lambda_t = Xvvv(X = X, v = theta_list[[t]])
    lambda_t = Tensor_vector(A = X,v = theta_list[[t]], m = 4)
    lambda_list = c(lambda_list,lambda_t)
  }
  
  ind = which.max(abs(lambda_list))
  theta_tau = theta_list[[ind]]
  
  theta_hat = power_itr(theta = theta_tau, X=X, N=N)
  #lambda_hat = Xvvv(X = X, v = theta_hat)
  lambda_hat = Tensor_vector(A = X,v = theta_hat, m = 4)
  
  
  def_X = X - lambda_hat * theta_hat %o% theta_hat %o% theta_hat %o% theta_hat
  
  result = list("theta_hat" = theta_hat, "lambda_hat" = lambda_hat, "def_X" = def_X)
  return(result)
} 

library(far)
Rob_TSM = function(X,L = 10, N = 10, rank = 10){
  p = dim(X)[1]
  result.1 = Est_PI(X = X, L = L, N = N)
  eigenv = matrix(result.1$theta_hat,nrow = length(result.1$theta_hat))
  eigenl = result.1$lambda_hat
  X = result.1$def_X
  
  for( m in 2:rank){
    result.m = Est_PI(X = X, L = L, N = N)
    eigenv = cbind(eigenv, result.m$theta_hat)
    eigenl = c(eigenl, result.m$lambda_hat)
    #eigenv = orthonormalization(eigenv,basis = F,norm = T)
    X = result.m$def_X
  }
  eigen_v = eigenv
  eigen_l = eigenl
  
  
  eigen_v = eigen_v[,order(abs(eigen_l),decreasing = T)]
  eigen_l = eigen_l[order(abs(eigen_l),decreasing = T)]
  
  
  result = list("eigenv" = eigen_v, "eigenl" = eigen_l)
}
