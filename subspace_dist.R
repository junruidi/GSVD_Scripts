subspace_dist = function(A,B){
  ATB = t(A) %*% B
  sigma = svd(ATB)$d
  theta = acos(sigma)
  
  grassmann = sqrt(sum(theta^2))
  asimov = theta[length(theta)]
  BC = sqrt(1 - prod(sigma^2))
  chordal = sqrt(sum((sin(theta))^2))
  FS = acos(prod(sigma))
  martin = sqrt(log(prod(1/sigma^2)))
  proc = 2*sqrt(sum((sin(theta/2))^2))
  proj = sin(theta[length(theta)])
  spectral = sin(theta[length(theta)]/2)
  
  dist = c(grassmann,asimov,BC,chordal,FS,martin,proc,proj,spectral)
  return(dist)
}