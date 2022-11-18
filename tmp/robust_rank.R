
robust_rank <- function(x) {
  robust_ranks <- matrix(x+rnorm(length(x)*100, 0, sd(x)), ncol=100, nrow=length(x))
  return(rank(round(apply(robust_ranks,1,mean))))
}
