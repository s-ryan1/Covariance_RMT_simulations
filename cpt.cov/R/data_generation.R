#' @export
generate.diff.mat = function(p){
    #Old version
    A = prcomp(rWishart(1, p, diag(p))[,,1])[[2]]
    eigs = rgamma(p, shape=5, scale=1/5)#runif(p,1/gamma, gamma)
    #eigs = map_dbl(eigs, random.invert)
    return(t(A)%*%diag(eigs)%*%A)

    #Using Wang Asymptotics
    #A = prcomp(rWishart(1, p, diag(p))[,,1])[[2]]
    #eigenvalues = c(1, runif(p-1,0,1))*sample(c(1,-1), p, replace=TRUE)
    #return(A%*%diag(eigenvalues)%*%t(A))
}


#' @export
gen.data = function(n, p, m, minseglen=30, cpts=NA, sigma=NA, error.dist="normal"){
    if(error.dist=="normal"){noise = rnorm(n*p, 0,1) %>% matrix(ncol=p)}
    if(error.dist=="poisson"){noise = rpois(n*p,1) %>% matrix(ncol=p)}
    if(error.dist=="exponential"){noise = rexp(n*p,1) %>% matrix(ncol=p)}
    if(error.dist=="bernoulli"){
        noise = runif(n*p) %>% matrix(ncol=p); 
        noise[noise<.5] = 0; noise[noise>=.5] = 1
    }
    if(error.dist=="uniform"){noise = runif(n*p) %>% matrix(ncol=p)}
    if(error.dist=="discrete"){
        noise = sample(seq(1,10), n*p, replace=T) %>% matrix(ncol=p)}
    if(error.dist=="chi_square"){
        noise = rnorm(n*p, 0,1) %>% `^`(2) %>% matrix(ncol=p)
    }
    if(error.dist=="student"){noise =rt(n*p, 5) %>% matrix(ncol=p)}
    if(error.dist=="cauchy"){noise=rcauchy(n*p) %>% matrix(ncol=p)}
	if(length(sigma)-1 != m & !(is.na(sigma[1]))){stop("Need 1 covariance matrix per segment")}
    if(length(cpts) != m & !(is.na(cpts[1]))){stop("Length of changepoint vector should equal number of changes")}
    if(is.na(cpts)){
        cpts = gen.cpt.locs(n,m,minseglen)
    }
    if(is.na(sigma[1])){
        sigma1 = (1/p)*rWishart(1, p, diag(p))[,,1]
        if(m>0){
            diff.mats = purrr::map(seq(1, m), ~generate.diff.mat(p))
            diff.mats = append(list(sigma1), diff.mats)
            sigma = purrr::accumulate(diff.mats, function(X,Y){B = expm::sqrtm(X); return(B%*%Y%*%B)})
        }else {
            sigma=list(sigma1)
        }
    }
    if(m>0){
            lengths = diff(c(0, cpts, n))
    } else {
            lengths = c(n)
    }

    if(m>0){
            locs= c(0, cpts, n)
    } else {
            locs = c(0,n)
    }

    data = map(seq(0,m), ~noise[(locs[.x+1]+1):locs[.x+2],])
    data = map2(data, sigma, ~.x%*%expm::sqrtm(.y))

    #purrr:map(lengths, 
    #data = purrr::map2(lengths, sigma, ~MASS::mvrnorm(.x,rep(0,p),.y))
    return(list(data=data, sigma=sigma, cpts=cpts))
}

#' @export
gen.cpt.locs <- function(n, m, minseglen = 30){
	#Generates m changepoint  locations for a multivariate time series
	
	#Catch errors if n,m or minseglen are non integers
	if(n%%1 != 0){ stop("n must be integer valued")}
	if(m%%1 != 0){ stop("m must be integer valued")}
	if(minseglen%%1 != 0){ stop("minseglen must be integer valued")}
	if(n<m){stop("More changepoints then data")} 
	if(n/m<minseglen){stop("Minseglen is greater than average segment length")} 

	if(m<1){return(NA)}
    X = seq(minseglen, n-minseglen)
    remove.segment = segment.remover(X, minseglen)
    cptlocs = purrr::map_dbl(seq(1,m), ~remove.segment()) 
	
    return(sort(cptlocs))
}

#' @export
uniform.distance = function(gamma){
    return(1  + (1/3)*gamma^2 + (1/gamma)*log(1-gamma) - (1/gamma)*log(1+gamma) + 1/(1-gamma^2))
}

segment.remover = function(X, minseglen){
    function(){
        i=sample(X,1)
        X <<- X[X<i-minseglen | X>i+minseglen]
        return(i)
    }
}

symmetric.random.matrix = function(p){
    c = MASS::mvrnorm(p*(p+1)/2 - p, 0, 1)     
    b = matrix(0, p,p)
    b[upper.tri(b, diag=FALSE)] = c
    return(b+t(b))
}

random.invert = function(x){
    x^(sample(c(-1,1),1))
}
