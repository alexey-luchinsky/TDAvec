computePS=function(D, homDim, p, scaleSeq){
  # D - N by 3 matrix (columns contain dimension, birth and death values respectively)
  # homDim - homological dimension (0 for H0, 1 for H1, etc.)
  # p - power of the weights for the silhouette function
  # scaleSeq - sequence of scale values for vectorization
  L <- length(scaleSeq)
  phi <- numeric(length = L)
  D <- D[D[, 1] == homDim,2:3,drop = FALSE]
  if (nrow(D)==0) return(phi)
  pp <- (D[,2]-D[,1])^p
  w <- pp/sum(pp)
  for (i in 1:L)
    phi[i] <- sum(w*pmax(pmin(scaleSeq[i] - D[,1], D[,2] - scaleSeq[i]),0))
  # returned object
  phi
}

