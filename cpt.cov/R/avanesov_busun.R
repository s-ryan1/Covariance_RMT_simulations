calculate.T.mat = function(sigma, lambda){
    prec.mat = sigma %>% glasso::glasso(lambda) %>% `$`(wi)
    return(prec.mat + t(prec.mat) - t(prec.mat)%*%sigma%*%prec.mat)
}

av.bun.stat.calculater = function(X, lambda){
    n = nrow(X)
    p = ncol(X)

    cumsum = matrix.cumsum(X)
    forward.cumsum = cumsum[[1]]
    backward.cumsum = cumsum[[2]]

    #Precision matrix for entire data 
    prec = glasso::glasso(forward.cumsum[[n]], lambda) %>% pluck("wi")
    S = map(seq(1,p), function(i){map_dbl(seq(1,p), ~prec[.x,.x]*prec[i,i] + prec[i,.x]^2)}) %>% reduce(c) 

    function(window){
        sigma1 =  map(seq(window+1, n-window),~(1/window)*(forward.cumsum[[.x-1]] - forward.cumsum[[.x-window]])) %>% map(calculate.T.mat, lambda) %>% map(as.vector)
        sigma2 =  map(seq(window+1, n-window),~(1/window)*(forward.cumsum[[.x+window]] - forward.cumsum[[.x-1]])) %>% map(calculate.T.mat, lambda) %>% map(as.vector)
        diff.mat = map2(sigma1, sigma2, ~ abs(.x - .y)) 
        return(c(rep(NA,window),map_dbl(diff.mat, ~max(sqrt(window/2)*sqrt(diag(1/S))%*%.x)), rep(NA, window)))
    }
}

#' @export
av.bun.stat = function(X, window, lambda){
    n = nrow(X)
    p = ncol(X)
    calculate.av.bun.stat = av.bun.stat.calculater(X, lambda)
    av.bun.stat = map(window, ~calculate.av.bun.stat(.x))
    av.bun.stat = reduce(av.bun.stat, cbind) 
    return(av.bun.stat)
    
    #test.value = av.bun.stat %>% apply(2,max, na.rm=TRUE)
    #cpt.est = av.bun.stat %>% apply(2,which.max)
    #return(cbind(window,cpt.est, test.value))
    #return(c(rep(NA,p+buffer),av.bun.stat,rep(NA,p+buffer)))
}

matrix.cumsum = function(X){
    n = nrow(X)
    products = purrr::map(as.data.frame(t(X)), ~.x%*%t(.x))
    forward.cumsum = purrr::accumulate(products, ~.x+.y)
    backward.cumsum = purrr::map(forward.cumsum, ~ -.x + forward.cumsum[[n]])
    return(list(forward.cumsum = forward.cumsum, backward.cumsum=backward.cumsum))
}
