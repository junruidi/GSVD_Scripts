#power tensor iteration on vectors

#power iteration for the first eigen vecotr
power_itr1 = function(theta,Y,N,order){
  n = nrow(Y)
  d = 1/n
  
  for(k in 1:N){
    w = (Y %*% theta)^(order-1)
    next_itr = as.vector(t(w) %*% Y)
    next_itr = d * next_itr
    
    thetai = next_itr/sqrt(sum(next_itr^2))
    err = sum((thetai - theta)^2)
    theta = thetai
    if( err <= 1e-6){break}
  }
  return(theta)
}

#power iteration for all following eigen vectors
power_itrk = function(theta,eigenv, eigenl,Y,N,order){
  n = nrow(Y)
  d = 1/n
  
  for(k in 1:N){
    w1 = (Y %*% theta)^(order-1)
    next_itr1 = as.vector(t(w1) %*% Y)
    
    w2 = (t(eigenv) %*% theta)^(order - 1)
    if(length(eigenl) == 1){next_itr2 = as.vector(eigenv %*% eigenl %*% w2)}
    if(length(eigenl) >1){next_itr2 = eigenv %*% diag(eigenl) %*% w2}
    
    next_itr = d * next_itr1 - next_itr2
    thetai = next_itr/sqrt(sum(next_itr^2))
    err = sum((thetai - theta)^2)
    theta = thetai
    if( err <= 1e-6){break}
  }
  return(theta)
}

#get first lambda
Est_PI1 = function(Y, L=10, N=10,order){
  p = ncol(Y)
  theta_list = list()
  for ( t in 1:L){
    vc = rnorm(p)
    theta_0 = vc/sqrt(sum(vc^2))
    theta_list[[t]] = power_itr1(theta = theta_0, Y = Y, N=N, order = order)
  }
  
  n = nrow(Y)
  d = 1/n
  lambda_list = NULL
  
  for(t in 1:L){
    a = (Y %*% theta_list[[t]])^order
    lambda_t = d * sum(a)
    lambda_list = c(lambda_list,lambda_t)
  }
  
  ind = which.max(abs(lambda_list))
  theta_tau = theta_list[[ind]]
  
  theta_hat = power_itr1(theta = theta_tau,Y = Y, N=N, order = order)
  
  a = (Y %*% theta_hat)^order
  lambda_hat = d * sum(a)
  
  result = list("theta_hat" = theta_hat, "lambda_hat" = lambda_hat)
  return(result)
} 


#get all following lambda

Est_PIk = function(Y,eigenv, eigenl,L=10, N=10,order){
  p = ncol(Y)
  theta_list = list()
  for ( t in 1:L){
    vc = rnorm(p)
    theta_0 = vc/sqrt(sum(vc^2))
    theta_list[[t]] = power_itrk(theta = theta_0, eigenv, eigenl, Y = Y, N=N, order = order)
  }
  
  n = nrow(Y)
  d = 1/n
  lambda_list = NULL
  
  for(t in 1:L){
    a =  sum((Y %*% theta_list[[t]])^order)
    if(length(eigenl) == 1){ b =  sum(eigenl * ((t(eigenv) %*% theta_list[[t]])^order))}
    if(length(eigenl) >1){b =  sum(diag(eigenl) %*% ((t(eigenv) %*% theta_list[[t]])^order))}
    lambda_t = d * a - b
    lambda_list = c(lambda_list,lambda_t)
  }
  
  ind = which.max(abs(lambda_list))
  theta_tau = theta_list[[ind]]
  theta_hat = power_itrk(theta = theta_tau,eigenv, eigenl,Y = Y, N=N, order = order)
  
  a =  sum((Y %*% theta_hat)^order)
  if(length(eigenl) == 1){ b =  sum(eigenl * ((t(eigenv) %*% theta_hat)^order))}
  if(length(eigenl) >1){b =  sum(diag(eigenl) %*% ((t(eigenv) %*% theta_hat)^order))}
  
  lambda_hat = d * a - b
  
  result = list("theta_hat" = theta_hat, "lambda_hat" = lambda_hat)
  return(result)
} 



Rob_TSM = function(Y,L = 10, N = 10,order, p = ncol(Y)){

  result.1 = Est_PI1(Y = Y, L = L, N = N, order = order)
  eigenv = matrix(result.1$theta_hat,nrow = length(result.1$theta_hat))
  eigenl = result.1$lambda_hat
  
  for( m in 2:p){
    result.m = Est_PIk(Y = Y,eigenv = eigenv, eigenl = eigenl, L = L, N = N, order = order)
    eigenv = cbind(eigenv, result.m$theta_hat)
    eigenl = c(eigenl, result.m$lambda_hat)
    print(m)
  }
  
  eigenv = eigenv[,order(abs(eigenl),decreasing = T)]
  eigenl = eigenl[order(abs(eigenl),decreasing = T)]
  
  
  result = list("eigenv" = eigenv, "eigenl" = eigenl)
}