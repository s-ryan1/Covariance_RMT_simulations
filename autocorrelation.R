library(cpt.cov)
library(expm)
library(parallel)
library(purrr)
library(rlang)
library(tidyverse)

run.methods = function(data, n, p, m, i){
    result = list()
    minseglen=floor(p*log(n))
    if(p < 15){
           mat.dist=matrix.dist.test.stat(data,minseglen)
           wang= wang.stat(data,minseglen)
           aue = aue.stat(data,minseglen)
    result[[1]] = 
    c(i,n,p, m,"Fisher", max(mat.dist, na.rm=T),
        which.max(mat.dist))
    result[[2]] = 
    c(i,n,p, m,"Wang", max(wang,na.rm=T),
        which.max(wang))
    result[[3]] = 
    c(i,n,p, m,"Aue", 
        max(aue,na.rm=T),
        which.max(aue))
    } else {
           mat.dist=matrix.dist.test.stat(data,minseglen)
           wang= wang.stat(data,minseglen)
    result[[1]] = 
    c(i,n,p,m,"Fisher", 
        max(mat.dist,na.rm=T),
        which.max(mat.dist))
    result[[2]] = 
    c(i,n,p,m,"Wang", 
        max(wang,na.rm=T),
        which.max(wang))
    }

    result = result %>% reduce(rbind)
    colnames(result) = c("iter", "n", "p", "m", "method", "value", "cpt.est")
    result = result %>% as_tibble %>% 
        mutate(iter = as.integer(iter),  n = as.integer(n), 
               p = as.integer(p), m = as.integer(m), 
               value= as.double(value), cpt.est=as.integer(cpt.est))
    return(result)	
}

run.sim = function(i, n, p, delta, phi){

        #set.seed(1234)
        #wish = rnorm(10*p*p,0, 1) %>% matrix(ncol=p) %>% cov
        #wish = list(wish, delta*wish)


	set.seed(353*n + 541*p + 7*i)
        noise = rnorm(n*p,0, 1) %>% matrix(ncol=p)
        noise = list(noise[1:(n/2),], delta*noise[-c(1:n/2),])
        #noise = map2(noise, wish, ~.x%*%sqrtm(.y))

        data = map(noise, noise_to_var, phi) %>% reduce(rbind)
        y=matrix.dist.test.stat(data, max(30,4*p))
        return(c(n=n, p=p, delta=delta, phi=phi,test.value = max(y,na.rm=T), cpt.est=which.max(y)))
}

noise_to_var = function(Z, phi){
    p=ncol(Z);n=nrow(Z);
    Y = matrix(0, nrow=n+1, ncol=p)
    Y[1,] = rnorm(p,0,1)
    for(i in 2:(n+1)){
        Y[i,] = phi*Y[i-1,] + Z[i-1,]
    }
    return(Y[-1,])
}



n = c(2000)
p = c(50)
delta = c(1,1.1,1.2)
phi = c(0,.1,.3,.6,.9)
run = seq(0,9)*100

parameters = cross(list(n,p,delta,phi,run))
#colnames(params) = c("n", "p", "delta", "phi")


args = commandArgs(trailingOnly=TRUE)
iter = as.numeric(args[1])
n= parameters[[iter]][[1]]
p= parameters[[iter]][[2]]
delta=parameters[[iter]][[3]]
phi=parameters[[iter]][[4]]
run= parameters[[iter]][[5]]

result=map(seq(1,100), ~run.sim(run+.x, n, p,delta, phi))

string = sprintf("~/Documents/covariance/autocorrelation/Results/run%d.Rds",iter)
saveRDS(result, string)
#

#,

