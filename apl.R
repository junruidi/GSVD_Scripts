aplSelect <- function(a, x, drop = FALSE) {
  sa <- aplShape(a)
  ra <- aplRank(a)
  sz <- sapply(x, length)
  z <- array(0, sz)
  nz <- prod(sz)
  for (i in 1:nz) {
    ivec <- aplEncode(i, sz)
    jvec <- vector()
    for (j in 1:ra)
      jvec <- c(jvec, x[[j]][ivec[j]])
    z[i] <- a[aplDecode(jvec, sa)]
  }
  if (drop)
    return(drop(z))
  else
    return(z)
}

#  drop

aplDrop <- function(a, x, drop = FALSE) {
  sa <- aplShape(a)
  ra <- aplRank(a)
  y <- as.list(rep(0, ra))
  for (i in 1:ra) {
    ss <- sa[i]
    xx <- x[i]
    sx <- ss + xx
    if (xx >= 0)
      y[[i]] <- (xx + 1):ss
    if (xx < 0)
      y[[i]] <- 1:sx
  }
  return(aplSelect(a, y, drop))
}

#  take

aplTake <- function(a, x, drop = FALSE) {
  sa <- aplShape(a)
  ra <- aplRank(a)
  y <- as.list(rep(0, ra))
  for (i in 1:ra) {
    ss <- sa[i]
    xx <- x[i]
    sx <- ss + xx
    if (xx > 0)
      y[[i]] <- 1:xx
    if (xx < 0)
      y[[i]] <- (sx + 1):ss
  }
  return(aplSelect(a, y, drop))
}

# reduce vector

aplRDV <- function(x, f = "+") {
  if (length(x) == 0)
    return(x)
  s <- x[1]
  if (length(x) == 1)
    return(s)
  for (i in 2:length(x))
    s <- match.fun(f)(s, x[i])
  return(s)
}

# scan vector

aplSCV <- function(x, f = "+") {
  if (length(x) <= 1)
    return(x)
  return(sapply(1:length(x), function(i)
    aplRDV(x[1:i], f)))
}

# inner product vector

aplIPV <- function(x, y, f = "*", g = "+") {
  if (length(x) != length(y))
    stop("Incorrect vector length")
  if (length(x) == 0)
    return(x)
  z <- match.fun(f)(x, y)
  return(aplRDV(z, g))
}

#  expand vector

aplEXV <- function(x, y) {
  z <- rep(0, length(y))
  m <- which(y == TRUE)
  if (length(m) != length(x))
    stop("Incorrect vector length")
  z[m] <- x
  return(z)
}

#  expand

aplExpand <- function(x, y, axis = 1) {
  if (is.vector(x))
    return(aplEXV(x, y))
  d <- dim(x)
  m <- which(y == TRUE)
  n <- length (y)
  e <- d
  e[axis] <- n
  if (length(m) != d[axis])
    stop("Incorrect dimension length")
  z <- array(0, e)
  for (i in 1:prod(d)) {
    k <- aplEncode (i, d)
    k[axis] <- m[k[axis]]
    z[aplDecode (k, e)] <- x[i]
  }
  return (z)
}

#  compress/replicate vector

aplCRV <- function(x, y) {
  n <- aplShape(x)
  m <- aplShape(y)
  if (m == 1)
    y <- rep(y, n)
  if (length(y) != n)
    stop("Length Error")
  z <- vector()
  for (i in 1:n)
    z <- c(z, rep(x[i], y[i]))
  return(z)
}

#  compress/replicate

aplReplicate <- function(x, y, k = aplRank (y)) {
  if (is.vector(x))
    return(aplCRV(x, y))
  sx <- aplShape(x)
  sy <- aplShape(y)
  sk <- sx[k]
  if (max(sy) == 1)
    y <- rep(y, sk)
  if (length(y) != sk)
    stop("Length Error")
  sz <- sx
  sz[k] <- sum(y)
  nz <- prod(sz)
  gg <- aplCRV(1:sk, y)
  z <- array(0, sz)
  for (i in 1:nz) {
    jvec <- aplEncode(i, sz)
    jvec[k] <- gg[jvec[k]]
    z[i] <- x[aplDecode(jvec, sx)]
  }
  return(z)
}

#  rotate vector

aplRTV <- function(a, k) {
  n <- aplShape(a)
  if (k > 0)
    return(c(a[-(1:k)], a[1:k]))
  if (k < 0)
    return(c(a[(n + k + 1):n], a[1:(n + k)]))
  return(a)
}

#  rotate

aplRotate <- function(a, b, axis = aplRank (a)) {
  if (is.vector(a))
    return(aplRTV(a, b))
  sa <- aplShape(a)
  sx <- aplShape(b)
  if (max(sx) == 1)
    b <- array(b, sa[-axis])
  if (!identical(sa[-axis], aplShape(b)))
    stop("Index Error")
  z <- array(0, sa)
  sz <- sa
  nz <- prod(sz)
  sk <- sz[axis]
  for (i in 1:nz) {
    ivec <- aplEncode(i, sz)
    xx <- b[aplDecode(ivec[-axis], sx)]
    ak <- rep(0, sk)
    for (j in 1:sk) {
      jvec <- ivec
      jvec[axis] <- j
      ak[j] <- a[aplDecode(jvec, sz)]
    }
    bk <- aplRTV(ak, xx)
    for (j in 1:sk) {
      jvec <- ivec
      jvec[axis] <- j
      z[aplDecode(jvec, sz)] <- bk[j]
    }
  }
  return(z)
}

# transpose -- will be overwritten by the C version

aplTranspose <- function(a, x = rev(1:aplRank(a))) {
  sa <- aplShape(a)
  ra <- aplRank(a)
  if (length(x) != ra)
    stop("Length Error")
  rz <- max(x)
  sz <- rep(0, rz)
  for (i in 1:rz)
    sz[i] <- min(sa[which(x == i)])
  nz <- prod(sz)
  z <- array(0, sz)
  for (i in 1:nz)
    z[i] <- a[aplDecode(aplEncode(i, sz)[x], sa)]
  return(z)
}

#  representation -- will be overwritten by the C version

aplEncode <- function(rrr, base) {
  b <- c(1, butLast(cumprod(base)))
  r <- rep(0, length(b))
  s <- rrr - 1
  for (j in length(base):1) {
    r[j] <- s %/% b[j]
    s <- s - r[j] * b[j]
  }
  return(1 + r)
}

#  base value -- will be overwritten by the C version

aplDecode <- function(ind, base) {
  b <- c(1, butLast(cumprod(base)))
  return(1 + sum(b * (ind - 1)))
}

# get

aplGet <- function(a, cell) {
  dims <- dim(a)
  n <- length(dims)
  b <- 0
  if (any(cell > dims) || any(cell < 1))
    stop("No such cell")
  return(a[aplDecode(cell, dims)])
}

# set

aplSet <- function(a, b, cell) {
  dims <- dim(a)
  n <- length(dims)
  if (any(cell > dims) || any(cell < 1))
    stop("No such cell")
  a[aplDecode(cell, dims)] <- b
  return(a)
}

#  join

aplJoin <- function(a, b, axis = 1) {
  if (is.vector(a) && is.vector(b))
    return(c(a, b))
  sa <- aplShape(a)
  sb <- aplShape(b)
  ra <- aplRank(a)
  rb <- aplRank(b)
  if (ra != rb)
    stop("Rank error in aplJoin")
  if (!identical(sa[-axis], sb[-axis]))
    stop("Shape error")
  sz <- sa
  sz[axis] <- sz[axis] + sb[axis]
  nz <- prod(sz)
  u <- unit(axis, ra)
  z <- array(0, sz)
  for (i in 1:nz) {
    ivec <- aplEncode(i, sz)
    if (ivec[axis] <= sa[axis])
      z[i] <- a[aplDecode(ivec, sa)]
    else
      z[i] <- b[aplDecode(ivec - sa[axis] * u, sb)]
  }
  return(z)
}

# ravel

aplRavel <- function(a) {
  as.vector(a)
}

# outer product

aplOuterProduct <- function(x, y, f = "*") {
  return(outer(x, y, f))
}

# shape

aplShape <- function(a) {
  if (is.vector(a))
    return(length(a))
  return(dim(a))
}

#  rank

aplRank <- function(a) {
  aplShape(aplShape(a))
}

# reshape

aplReshape <- function(a, d) {
  return(array(a, d))
}



# inner product -- will be overwritten by the C version

aplInnerProduct <- function(a, b, f = "*", g = "+") {
  sa <- aplShape(a)
  sb <- aplShape(b)
  ra <- aplRank(a)
  rb <- aplRank(b)
  ia <- 1:(ra - 1)
  ib <- (ra - 1) + (1:(rb - 1))
  ff <- match.fun(f)
  gg <- match.fun(g)
  ns <- last(sa)
  nt <- first(sb)
  if (ns != nt)
    stop("Incompatible array dimensions")
  sz <- c(butLast(sa), butFirst(sb))
  nz <- prod(sz)
  z <- array(0, sz)
  for (i in 1:nz) {
    ivec <- aplEncode(i, sz)
    for (j in 1:ns) {
      aa <- a[aplDecode(c(ivec[ia], j), sa)]
      bb <- b[aplDecode(c(j, ivec[ib]), sb)]
      tt <- ff(aa, bb)
      if (j == 1)
        z[i] <- tt
      else
        z[i] <- gg(z[i], tt)
    }
  }
  return(z)
}

# member of

aplMemberOf <- function(a, b) {
  sa <- aplShape(a)
  sb <- aplShape(b)
  na <- prod(sa)
  nb <- prod(sb)
  z <- array (0, sa)
  for (i in 1:na) {
    aa <- a[i]
    for (j in 1:nb)
      if (aa == b[j])
        z[i] <- 1
  }
  return(z)
}

# utilities below

first <- function(x) {
  return(x[1])
}

butFirst <- function(x) {
  return(x[-1])
}

last <- function(x) {
  return(x[length(x)])
}

butLast <- function(x) {
  return(x[-length(x)])
}

unit <- function(i, n) {
  ifelse(i == (1:n), 1, 0)
}
