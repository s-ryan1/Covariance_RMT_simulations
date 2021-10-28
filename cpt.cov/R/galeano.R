#' @import purrr
#' @import rlang 

#' @export
galeano.cusum.stat = function(X, minseglen){
    n = nrow(X)
    p = ncol(X)
    Delta_Werk = p*(p+1)/2 
    sigma = cov(X) 
    X = purrr::map(seq(1,n), ~X[.x,])
    X = map_dbl(X, ~.x%*%solve(sigma, .x))
    A = cumsum(X)
    A = seq(1,n)*(A/seq(1,n) - A[n]/n)/(sqrt(2*p*n))
    A[1:minseglen] = NA
    A[(n-minseglen):n] = NA
    return(abs(A))
    #cusum = imap(A, ~(.y/sqrt(2*p*n))*(A[[.y]]/.y - A[[n]]/n))
    #return(c(rep(NA,p+minseglen),(aue.cusum - Delta_Werk/4)/sqrt(Delta_Werk/8), rep(NA,p+minseglen)))
}

###' @export
#sample.brownian.bridge.int = function(p){
#    B = purrr::map(seq(1, p*(p+1)/2), ~sde::BBridge()^2)
#    return(max(purrr::reduce(B, `+`)))
#}
