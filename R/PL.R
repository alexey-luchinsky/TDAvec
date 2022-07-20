computePL=function(D, homDim, k, scaleSeq){
  # D - N by 3 matrix (columns contain dimension, birth and death values respectively)
  # homDim - homological dimension (0 for H0, 1 for H1, etc.)
  # k - order of landscape function
  # scaleSeq - sequence of scale values for vectorization
  L <- length(scaleSeq)
  lambda <- numeric(length = L)
  D <- D[D[, 1] == homDim,2:3,drop = FALSE]
  if (nrow(D)==0) return(lambda)
  for (i in 1:L)
    lambda[i] <- sort(pmax(pmin(scaleSeq[i] - D[,1], D[,2] - scaleSeq[i]),0),decreasing = TRUE)[k]
  # returned object
  lambda
}

