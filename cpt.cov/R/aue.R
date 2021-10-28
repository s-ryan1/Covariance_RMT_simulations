#' @import purrr
#' @import rlang 

aue.cusum.calculater = function(X){
    A = map(X, ~.x%*%t(.x))
    A = accumulate(A, `+`)
    A = map(A, vech)
    n = length(X)
    function(tau){
        return((1/(sqrt(n)))*(A[[tau]] - (tau/(n))*(A[[n]])))
    }
}

calculate.aue.cov = function(X){
    n = length(X)
    sigma = purrr::map(X, ~.x%*%t(.x)) %>% purrr::reduce(`+`)
    X = map(X, ~vech(.x%*%t(.x)))
    a = map(X, ~.x%*%t(.x)) %>% purrr::reduce(`+`)
    return(a/n - (1/n)*vech(sigma)%*%t((1/n)*vech(sigma)))
    #lvar = map(X, ~vech(.x%*%t(.x))) %>% reduce(rbind) %>% 
#        sandwich::lrvar(type="Newey-West", prewhite=T, adjust=T) #%>% 
#    return(n*lvar)
}

vech = function(Y){
    return(as.vector(Y[upper.tri(Y, diag=TRUE)]))
}

#' @export
aue.stat = function(X, minseglen){
    n = nrow(X)
    p = ncol(X)
    Delta_Werk = p*(p+1)/2 
    X = purrr::map(seq(1,n), ~X[.x,])
    calculate.aue.cusum = aue.cusum.calculater(X)
    aue.cov= calculate.aue.cov(X)
    aue.cusum = purrr::map(seq(p+minseglen+1,n-p-minseglen), calculate.aue.cusum)
    aue.cusum = purrr::map_dbl(aue.cusum, ~.x%*%solve(aue.cov, .x))
    return(c(rep(NA,p+minseglen),(aue.cusum - Delta_Werk/4)/sqrt(Delta_Werk/8), rep(NA,p+minseglen)))
}

#' @export
sample.brownian.bridge.int = function(p){
    B = purrr::map(seq(1, p*(p+1)/2), ~sde::BBridge()^2)
    return(max(purrr::reduce(B, `+`)))
}
