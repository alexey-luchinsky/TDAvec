computePI <- function(D,homDim,res,sigma,minB,maxB,minP,maxP){
  # D - N by 3 matrix (columns contain dimension, birth and persistence values respectively)
  # PS for H0
  PSurfaceH0 = function(point) {
    y = point[2]
    out2 = pnorm(y_upper, mean = y, sd = sigma) - pnorm(y_lower,mean = y, sd = sigma)
    wgt = y/maxP * (y < maxP) + 1 * (y >= maxP)
    return(out2 * wgt)
  }
  # PS for Hk, k>0
  PSurfaceHk = function(point) {
    x = point[1]
    y = point[2]
    out1 = pnorm(x_upper, mean = x, sd = sigma) - pnorm(x_lower,mean = x, sd = sigma)
    out2 = pnorm(y_upper, mean = y, sd = sigma) - pnorm(y_lower,mean = y, sd = sigma)
    wgt = y/maxP * (y < maxP) + 1 * (y >= maxP)
    return(out1 %o% out2 * wgt)
  }
  # Body of computePI()
  D <- D[D[,1]==homDim,2:3,drop=F]
  if (nrow(D)==0) return(out=numeric(length = res^2))
  dy = (maxP-minP)/res
  y_lower = seq(minP,maxP-dy,by=dy)
  y_upper = y_lower + dy
  if (sum(D[,1])==0) {
    Psurf_mat = apply(D, 1, PSurfaceH0)
  } else{
    dx = (maxB-minB)/res
    x_lower = seq(minB,maxB-dx,by=dx)
    x_upper = x_lower + dx
    Psurf_mat = apply(D,1,PSurfaceHk)
  }
  out = rowSums(Psurf_mat)
  return(out)
}