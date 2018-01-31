require("partitions")
source("~/Dropbox/Junrui Di/tensor analysis/multi-cumulant/script/cumulant/apl.R")


raw_moments_upto_p <- function (x, p = 4) {
   n <- nrow (x)
   m <- ncol (x)
   if (p == 1) {
     return (c (1, apply (x, 2, mean)))
     }
   y <- array (0, rep (m + 1, p))
   for (i in 1 : n) {
     xi <- c (1, x[i, ])
     z <- xi
     for (s in 2:p) {
       z <- outer (z, xi)
       }
     y <- y + z
     }
   return (y / n)
   }

cumulants_from_raw_moments <- function (raw) {
  dimr <- dim (raw)
  nvar <- dimr[1]
  cumu <- array (0, dimr)
  nele <- prod (dimr)
  ldim <- length (dimr)
  spp <<- rep (list (0), ldim)
  qpp <<- rep (0, ldim)
  rpp <<- rep (list (0), ldim)
  for (i in 1 : ldim) {
    spp[[i]] <<- setparts (i)
    qpp[[i]] <<- factorial (i)
    if (i %% 2) {
      qpp[[i]] <<- -qpp[[i]]}
    rpp[[i]] <<- aplSelect (raw, c (rep (list (1 : nvar), i), rep (list (1 : 1), ldim - i)))
    }
  qpp <<- c(1, qpp)
  for (i in 2 : nele) {
    ind <- aplEncode (i, dimr)
    cumu[i] <- one_cumulant_from_raw_moments (ind, raw)
    }
  return (cumu)
  }

cumulants_upto_p <- function (x, p = 4) {
  return (cumulants_from_raw_moments (raw_moments_upto_p
                                         (x, p)))}

first_four_cumulants <- function (x) {
 cumu <- cumulants_upto_p (x)
  nsel <- 2 : dim (cumu)[1]
  return (list (c1 = cumu[1, 1, 1, nsel],
                   c2 = cumu[1, 1, nsel, nsel],
                   c3 = cumu[1, nsel, nsel, nsel],
                   c4 = cumu[nsel, nsel, nsel, nsel]))
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
  for (i in mm) for (j in mm) for (k in mm) {
    s3 <- r3 [i, j, k]
    s21 <- r2 [i, j] * r1 [k] + r2 [i, k] * r1 [j] + r2[j, k] * r1 [i]
    s111 <- r1 [i] * r1 [j] * r1 [k]
    c3 [i, j, k] <- s3 - s21 + 2 * s111
    }
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
   return (list(c1 = r1, c2 = c2, c3 = c3, c4 = c4))
}