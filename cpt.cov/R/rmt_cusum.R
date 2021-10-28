## @export
#cusum.estimator = function(data, epsilon, f){
#    n = nrow(data)
#    products = purrr::map(as.data.frame(t(data)), ~.x%*%t(.x))
#    forward.cumsum = purrr::accumulate(products, ~.x+.y)
#    backward.cumsum = purrr::map(forward.cumsum, ~ -.x + forward.cumsum[[nrow(data)]])
#
#    function(tau){
#        A = forward.cumsum[[tau]] - (tau/n)*backward.cumsum[[tau]]
#        return((1/sqrt(n))*sum(diag(A)))
#    }
#}


## @export
#matrix.dist.test.stat = function(data,epsilon){
#    p = ncol(data)
#    n = nrow(data)
#    t = seq(p+1,n-p-1)
#    #t = keep(t, ~p/.x < .5)
#    #t = keep(t, ~p/(n-.x) < .5)
#    estimate.cusum = cusum.estimator(data,epsilon,mean)
#    test.stat = purrr::map_dbl(t, estimate.cusum)
#    dimension.over.length = purrr::map(t, ~c(p/.x, p/(n-.x)))
#    trace = map(dimension.over.length, ~exec(calculate.expected.trace,!!!.x) )
#    values = map_lgl(trace, ~length(.x[[1]])==0)
#    bias = map_dbl(dimension.over.length, ~exec(asymptotic.bias,!!!.x))
#    variance = map_dbl(dimension.over.length, ~exec(asymptotic.variance,!!!.x) )
#    trace = map_if(trace, values, ~NA) %>% map_dbl(~.x[[1]])
#    bias = map_if(bias, values, ~NA)  
#    variance= map_if(variance, values, ~NA) 
#    test.stat = pmap_dbl(list(test.stat, trace, bias, variance), ~(p*(..1 - ..2) - ..3)/sqrt(..4))
#        return(c(rep(NA,p),test.stat, rep(NA,p)))
#}
