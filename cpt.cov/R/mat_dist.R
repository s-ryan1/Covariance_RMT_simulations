#' @export
covariance.distance.estimator = function(data, epsilon, f){
    products = purrr::map(as.data.frame(t(data)), ~.x%*%t(.x))
    forward.cumsum = purrr::accumulate(products, ~.x+.y)
    backward.cumsum = purrr::map(forward.cumsum, ~ -.x + forward.cumsum[[nrow(data)]])
    function(tau){
        sigma1 =  (1/tau)*forward.cumsum[[tau]] + epsilon*diag(ncol(data))
        sigma2 =  (1/(nrow(data)-tau))*backward.cumsum[[tau]] + epsilon*diag(ncol(data))
        output = tryCatch(cov.dist(sigma1, sigma2,f),
    error = function(error_message){return(NA)})
        return(output)
    }
}

#' @export
cov.dist= function(sigma1, sigma2, f){
    A = geigen::geigen(sigma2, sigma1, symmetric=TRUE)$values 
    return(rlang::exec(f,(A-1)^2) + rlang::exec(f,(1/A-1)^2))
}

#' @export
biased.matrix.dist.test.stat = function(data, epsilon){
    p = ncol(data)
    n = nrow(data)
    t = seq(p+1,n-p-1)
    estimate.covariance.distance = covariance.distance.estimator(data,epsilon,mean)
    test.stat = purrr::map_dbl(t, estimate.covariance.distance)
    #dimension.over.length = purrr::map(t, ~c(p/.x, p/(n-.x)))
    #trace = map(dimension.over.length, ~exec(calculate.expected.trace,!!!.x) )
    #values = map_lgl(trace, ~length(.x[[1]])==0)
    #bias = map_dbl(dimension.over.length, ~exec(asymptotic.bias,!!!.x))
    #variance = map_dbl(dimension.over.length, ~exec(asymptotic.variance,!!!.x) )
    #trace = map_if(trace, values, ~NA) %>% map_dbl(~.x[[1]])
    #bias = map_if(bias, values, ~NA)  
    #variance= map_if(variance, values, ~NA) 
    #test.stat = pmap_dbl(list(test.stat, trace, bias, variance), ~(p*(..1 - ..2) - ..3)/sqrt(..4))
        return(c(rep(NA,p),p*test.stat, rep(NA,p)))
}

#' @export
matrix.dist.test.stat = function(data, minseglen){
    p = ncol(data)
    n = nrow(data)
    t = seq(p+minseglen+1,n-p-minseglen)
    #t = keep(t, ~p/.x < .5)
    #t = keep(t, ~p/(n-.x) < .5)
    estimate.covariance.distance = covariance.distance.estimator(data,0,mean)
    test.stat = purrr::map_dbl(t, estimate.covariance.distance)
    dimension.over.length = purrr::map(t, ~c(p/.x, p/(n-.x)))
    trace = map(dimension.over.length, ~exec(calculate.expected.trace,!!!.x) )
    values = map_lgl(trace, ~length(.x[[1]])==0)
    bias = map_dbl(dimension.over.length, ~exec(asymptotic.bias,!!!.x))
    variance = map_dbl(dimension.over.length, ~exec(asymptotic.variance,!!!.x) )
    trace = map_if(trace, values, ~NA) %>% map_dbl(~.x[[1]])
    bias = map_if(bias, values, ~NA)  
    variance= map_if(variance, values, ~NA) 
    test.stat = pmap_dbl(list(test.stat, trace, bias, variance), ~(p*(..1 - ..2) - ..3)/sqrt(..4))
        return(c(rep(NA,p+minseglen),abs(test.stat), rep(NA,p+minseglen)))
}


function.prod = function(f1,f2){
    return(function(x){ return(f1(x)*f2(x))})
}  

#' @export 
calculate.expected.trace = function(y1, y2){
    asymptotic.pdf = fisher.esd(y1, y2)
    integrand = function.prod(function(x){(1-x)^2 + (1-1/x)^2}, asymptotic.pdf)
    asymptotic.supports = construct.fisher.support(y1, y2)
    safe_integral = safely(integrate)
    integral = exec(safe_integral, integrand, `!!!`(asymptotic.supports))[[1]]
    return(integral)
}

#' @export
asymptotic.bias = function(y1,y2){
    h = sqrt(y1 + y2 - y1*y2)
    K_1 = 2*h*(1+h^2)/(1-y2)^4 - 2*h/(1-y2)^2
    J_1 = 2*h*(1+h^2)/(1-y1)^4 - 2*h/(1-y1)^2

    return(2*(h^2 - y2^2)/(1-y2)^4 + 2*K_1*y2/h + 2*(h^2 - y1^2)/(1-y1)^4+ 2*J_1*y1/h)

    #K_2 = 2*(1-h^2)/(1-y2)^4 - 2*h/(1-y2)^2
    #J_1 = -2*(1-y2)^2
    #J_2 = (1-y2)^4

    #return(2*(h^2 - y2^2) + 2*K_2*y2 + J_2*(2*(3*h^2+1)/(h^2-1)^4 + (y2-1)^2/(h^2-1)^2) +
    #    2*J_1*(h^2-y2)/((1-y2)*(h^2-1)^2))
}

#' @export
asymptotic.variance = function(y1,y2){
    h = sqrt(y1 + y2 - y1*y2)

    K_21 = 2*h*(1+h^2)/(1-y2)^4 - 2*h/(1-y2)^2
    K_22 = 2*h*(1+h^2)/(1-y1)^4 - 2*h/(1-y1)^2
    K_31 = h^2/(1-y2)^4
    K_32 = h^2/(1-y1)^4
    J_1 = -2*(1-y2)^2
    J_2 = (1-y2)^4

    var_x = K_21^2 + 2*K_31^2
    var_y = K_22^2 + 2*K_32^2


    cov_xy = J_1*K_21/h + 
        J_1*K_21/(h*(h^2-1)) + 
        (-J_1*K_31*(h^2+1)/h^2) + 
        (-J_1*K_31/(h^2*(h^2-1))) +
        J_2*K_21*2*h/(h^2-1)^3 +
        J_2*K_31/h^2 +
        J_2*K_31*((1-3*h^2)/(h^2*(h^2-1)^3))
         #J_1*h*(K_21 - K_3*h)/(h^2-1) + J_2*h*(2*K_21 + K_3*h*(h^2-3))/(h^2-1)^3 
     
    return(3*(var_x + var_y + 2*cov_xy))
    #J_1 = -2*(1-y2)^2
    #J_2 = (1-y2)^4
    #L_1 = J_1*h/a
    #L_2 = 2*J_2*h/a^3
    #L_3 = 2*J_2*h^4/a^3
    #
    #print(c(h,a, K_2, K_3, J_1, J_2, L_1, L_2, L_3))

    #print(c(
    #2*J_2*K_2*h/(a^3),  
    #+2*J_2*K_3*h^2*(h^2-3)/a^3,  
    #+2*J_2*L_1*(h^4-1)/a^6,
    #+J_2*L_2*h*(-3*h^4+h^2 +3)/a^7, 
    #+J_2*L_3*h^2*(2*h^4 + h^2 -3)/a^7,  
    #+J_1*K_2*h/a , 
    #-2*J_1*K_3*h^2/a, 
    #+J_1*L_1*h/a^3 , 
    #-J_1*L_2*h/a^4 ,
    #+J_1*L_3*h^2/a^4,  
    #+K_2*(K_2 + L_1 + L_2), 
    #+K_3*(L_3 + 2*K_3 - 2*h*L_1 - 3*h*L_2)))
    #
    #limiting.var= ( 
    #    2*J_2*K_2*h/a^3  
    #    +2*J_2*K_3*h^2*(h^2-3)/a^3  
    #    +2*J_2*L_1*(h^4-1)/a^6 
    #    +J_2*L_2*h*(-3*h^4+h^2+3)/a^7 
    #    +J_2*L_3*h^2*(2*h^4 + h^2 -3)/a^7  
    #    +J_1*K_2*h/a  
    #    -2*J_1*K_3*h^2/a 
    #    +J_1*L_1*h/a^3  
    #    -J_1*L_2*h/a^4 
    #    +J_1*L_3*h^2/a^4  
    #    +K_2*(K_2 + L_1 + L_2) 
    #    +K_3*(L_3 + 2*K_3 - 2*h*L_1 - 3*h*L_2)
    #    )
    return(limiting.var)
}

#' @export
limiting.expectation = function(n, p){
    t = seq(p+1,n-p)
    dimension.over.length = purrr::map(t, ~c(p/.x, p/(n-.x)))
    trace = map(dimension.over.length, ~exec(calculate.expected.trace,!!!.x) )
    values = map_lgl(trace, ~length(.x[[1]])==0)
    trace = map_if(trace, values, ~NA) %>% map_dbl(~.x[[1]])
        return(trace)
}

#' @export
limiting.bias = function(n,p){
    t = seq(p+1,n-p)
    dimension.over.length = purrr::map(t, ~c(p/.x, p/(n-.x)))
    bias = map_dbl(dimension.over.length, ~exec(asymptotic.bias,!!!.x))
    return(bias)
}

#' @export
limiting.variance = function(n,p){
    t = seq(p+1,n-p)
    dimension.over.length = purrr::map(t, ~c(p/.x, p/(n-.x)))
    variance = map_dbl(dimension.over.length, ~exec(asymptotic.variance,!!!.x) )
    return(variance)
}

#' @export
calculate.minseglen = function(n, p, func, alpha=2, constraint=1/n){
    #alpha 
    if(func == "linear"){
        return(alpha*p)
    }  else if( func == "log-linear"){
        return(log(n)*p)
    } else if( func == "constrained-gradient"){
        grad = limiting.variance(n,p) %>% diff %>% abs 
        return(p + min(which(grad < constraint)))
    }
}
