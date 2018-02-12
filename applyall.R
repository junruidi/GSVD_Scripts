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


four_cumulants_direct <- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  mm <- 1 : m
  nn <- 1 : n
  r1 <- colSums (x) / n
  r2 <- crossprod (x) / n
  r3 <- array (0, c (m, m, m))
  r4 <- array (0, c (m, m, m, m))
  for (i in nn) {
    r3 <- r3 + outer (outer (x [i, ], x [i, ]), x [i,])
    r4 <- r4 + outer (outer (x [i, ], x [i, ]), outer (x [i, ], x [i, ]))
  }
  r3<-r3 / n
  r4<-r4 / n
  c2<-r2 - outer (r1, r1)
  c3<-r3
  c4<-r4
  for (i in mm) for (j in mm) for (k in mm) for (l in mm)
  {
    s4 <- r4 [i, j, k, l]
    s31 <- r3 [i, j, k] * r1 [l] + r3 [i, j, l] * r1 [k] + r3 [i, k, l] * r1 [j] + r3 [j, k, l] * r1 [i]
    s22 <- r2 [i, j] * r2 [k, l] + r2 [i, k] * r2 [j, l] + r2 [j, k] * r2 [i, l]
    s211 <- r2 [i, j] * r1 [k] * r1 [l] + r2 [i, k] *
      r1 [j] * r1 [l] + r2 [i, l] * r1 [k] * r1 [j] +
      r2 [j, k] * r1 [i] * r1 [l] + r2 [j, l] * r1 [i] * 
      r1 [k] + r2 [k, l] * r1 [i] * r1 [j]
    s1111 <- r1 [i] * r1 [j] * r1 [k] * r1 [l]
    
    c4 [i, j, k, l] <- s4 - s31 - s22 + 2 * s211 - 6 * s1111
  }
  return (c4)
}

#3. HOEVD
hoevd = function(x,rank){
  tnsr = as.tensor(x)
  unfold = k_unfold(tnsr,m = 1)@data
  u = svd(unfold, nu = rank)$u
  ulist = list(u,u,u)
  z = ttl(tnsr, lapply(ulist, t), ms = 1:tnsr@num_modes)@data
  return(list( u = u, z = z))
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