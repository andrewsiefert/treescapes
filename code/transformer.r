# transform and backtransform functions -----------------------------------

transform <- function(x, trans, log = F) {
  if (isTRUE(log)) x <- log(x)
  x_s <- (x - attr(trans, "scaled:center"))/attr(trans, "scaled:scale")
  return(x_s)
}

backtransform <- function(x_s, trans, log = F) {
  x <- x_s * attr(trans, "scaled:scale") + attr(trans, "scaled:center")
  if(isTRUE(log)) x <- exp(x)
  return(x)
}