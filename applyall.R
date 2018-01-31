library(rTensor)
#all functions to be used regarding our analysis 
#1. moment/gram tensors
MGT3<- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  mm <- 1 : m
  nn <- 1 : n
  c3 <- array (0, c (m, m, m))
  for (i in nn) {
    c3 <- c3 + outer (outer (x [i, ], x [i, ]), x [i,])
  }
  return (c3)
}


MGT4 <- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  mm <- 1 : m
  nn <- 1 : n
  c4 <- array (0, c (m, m, m, m))
  for (i in nn) {
    c4 <- c4 + outer (outer (x [i, ], x [i, ]), outer (x [i, ], x [i, ]))
  }
  return (c4)
}

#2. cumulant tensors
third_cumulant_p<- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  mm <- 1 : m
  nn <- 1 : n
  c3 <- array (0, c (m, m, m))
  for (i in nn) {
    c3 <- c3 + outer (outer (x [i, ], x [i, ]), x [i,])
  }
  return (c3)
}


four_cumulants_p <- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  mm <- 1 : m
  nn <- 1 : n
  r2 <- crossprod (x)
  c4 <- array (0, c (m, m, m, m))
  for (i in nn) {
    c4 <- c4 + outer (outer (x [i, ], x [i, ]), outer (x [i, ], x [i, ]))
  }
  for (i in mm) for (j in mm) for (k in mm) for (l in mm)
  {
    s4 <- r4 [i, j, k, l]
    s22 <- r2 [i, j] * r2 [k, l] + r2 [i, k] * r2 [j, l] + r2 [j, k] * r2 [i, l]
    
    c4 [i, j, k, l] <- s4 - s22
  }
  return (c4)
}

#3. HOEVD
hoevd = function(x,rank){
  tnsr = as.tensor(x)
  unfold = k_unfold(tnsr,m = 1)@data
  u = svd(unfold, nu = rank)$u
  return(u)
}


#4. HOSVD for gram tensors
Gram3_hosvd = function(y){
  n = dim(y)[1]
  p = dim(y)[2]
  
  svd_y = svd(y)
  u2 = svd_y$u
  v2 = svd_y$v
  s2 = diag(svd_y$d)
  
  vvv = MGT3(v2)
  middle = s2 %*% k_unfold(as.tensor(vvv),m=1)@data %*% (s2 %x% s2)
  middle_u = svd(middle)$u
  
  u3 = u2 %*% middle_u
  return(u3)
}

Gram4_hosvd = function(y){
  n = dim(y)[1]
  p = dim(y)[2]
  
  svd_y = svd(y)
  u2 = svd_y$u
  v2 = svd_y$v
  s2 = diag(svd_y$d)
  
  vvvv = MGT4(v2)
  middle = s2 %*% k_unfold(as.tensor(vvvv),m=1)@data %*% (s2 %x% s2 %x% s2)
  middle_u = svd(middle)$u
  
  u4 = u2 %*% middle_u
  return(u4)
}