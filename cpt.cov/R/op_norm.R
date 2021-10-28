wang.cusum.calculater = function(X){
    A = purrr::map(X, ~.x%*%t(.x))
    A = purrr::accumulate(A, `+`)
    n = length(X)
    function(tau){
        return(sqrt((n-tau)/(n*tau))*A[[tau]] - sqrt((tau)/(n*(n-tau)))*(A[[n]] - A[[tau]]))
    }
}

#' @export
wang.stat = function(X, minseglen){
    X = as.data.frame(t(X))
    n = length(X)
    p = length(X[[1]])
    calculate.wang.cusum = wang.cusum.calculater(X)
    aue.cusum = map(seq(p+minseglen+1,n-p-minseglen), ~calculate.wang.cusum(.x))
    #print(aue.cusum)
    return(c(rep(NA,p+minseglen),map_dbl(aue.cusum, ~max(abs(eigen(.x)[[1]]))),rep(NA,p+minseglen)))
}

#' @export
wang.threshold = function(data){
    n=nrow(data); p=ncol(data)
    B_sq = data %>% cov %>% eigen %>% pluck("values") %>% max
    return((B_sq/2)*sqrt(p*log(n)))
}
