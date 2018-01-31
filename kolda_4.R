#kolda 4
library(rTensor)
Tensor_vector = function(A,v){
  tnsr = as.tensor(A)
  V = matrix(v,ncol = 1)
  m = length(dim(A))
  c.t = tnsr
  for(i in 1:m){
    c.t = ttm(c.t,t(V),m=i)
  }
  result = c.t@data
}


kolda_4 = function(A){
  p = dim(A)[1]
  
  beta = runif(p)
  B = matrix(0,ncol = p, nrow = p)
  for(i in 1:p){
    for(j in 1:p){
      B = B + beta[i] * A[,,i,j]
    }
  }
  v = eigen(B)$vectors
  
  lambda = NULL
  for(j in 1:ncol(v)){
    lambda = c(lambda,Tensor_vector(A,v[,j]))
  }
  
  eigen_v = v
  eigen_l = lambda
  
  
  eigen_v = eigen_v[,order(abs(eigen_l),decreasing = T)]
  eigen_l = eigen_l[order(abs(eigen_l),decreasing = T)]
  
  return(list("eigen_v" = eigen_v, "eigen_l" = eigen_l))
}