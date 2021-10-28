#' @export
whiten.test.stat = function(X){
    n = nrow(X)
    X = whiten(X)
    X = map(seq(1,ncol(X)), ~c(0,X[,.x]) %>%  `^`(2) %>% cumsum)  
    Y = map(X, ~.x[n+1] - .x)
    test.stat = map2(X,Y, ~map_dbl(seq(2, n-2), function(i){n*log(.x[n+1]/n) - i*log(.x[i+1]/i) - (n-i)*log(.y[i+1]/(n-i))})) 
    test.stat = test.stat %>% reduce(`+`)
    return(test.stat)
}
