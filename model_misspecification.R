library(cpt.cov)
library(expm)
library(purrr)
library(rlang)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
iter = as.numeric(args[1])

#run.methods = function(data, n, p, m, i, distribution){
#    result = list()
#    minseglen=floor(p*log(n))
#    if(p < 15){
#           mat.dist=matrix.dist.test.stat(data,minseglen)
#           wang= wang.stat(data,minseglen)
#		safe_aue = safely(aue.stat)
#		aue = safe_aue(data,minseglen)
#if(length(aue$result)==0){return(NULL)}else{aue = aue$result} 
#    result[[1]] = 
#    c(i,n,p, m,"Fisher", max(mat.dist, na.rm=T),
#        which.max(mat.dist), distribution)
#    result[[2]] = 
#    c(i,n,p, m,"Wang", max(wang,na.rm=T),
#        which.max(wang), distribution)
#    result[[3]] = 
#    c(i,n,p, m,"Aue", 
#        max(aue,na.rm=T),
#        which.max(aue), distribution)
#    } else {
#           mat.dist=matrix.dist.test.stat(data,minseglen)
#           wang= wang.stat(data,minseglen)
#    result[[1]] = 
#    c(i,n,p,m,"Fisher", 
#        max(mat.dist,na.rm=T),
#        which.max(mat.dist), distribution)
#    result[[2]] = 
#    c(i,n,p,m,"Wang", 
#        max(wang,na.rm=T),
#        which.max(wang), distribution)
#    }
#
#    result = result %>% reduce(rbind)
#    colnames(result) = c("iter", "n", "p", "m", "method", "value", "cpt.est", "error.dist")
#    result = result %>% as_tibble %>% 
#        mutate(iter = as.integer(iter),  n = as.integer(n), 
#               p = as.integer(p), m = as.integer(m), 
#               value= as.double(value), cpt.est=as.integer(cpt.est))
#    return(result)	
#}
#
#sim.runner = function(n, p){
#	function(i){
#	print(i)
#        distributions = c("normal", "exponential", "uniform", 
#                          "student")
#        result = list() 
#        for(k in 1:4){
#            set.seed(2903*n + 89*p + 5*i)
#            result[[k]] = gen.data(n,p,0,error.dist=distributions[k])$data %>% reduce(rbind) %>% run.methods(n,p,0,i,distributions[k]) 
#            set.seed(2903*n + 89*p + 5*i)
#            result[[k + 4]] = gen.data(n,p,1,error.dist=distributions[k], cpts=n/2)$data %>% reduce(rbind) %>% run.methods(n,p,1,i,distributions[k]) 
#        }
#print(result)
#	result = result %>% discard(~length(.x)==0)
#        result %>% reduce(rbind) %>% 
#            return
#    }
#}

run.sim = function(i, n, p, delta){

    result=list()

	set.seed(353*n + 541*p + 7*i)
    gaussian.noise = rnorm(n*p,0, 1) %>% matrix(ncol=p)
    gaussian.noise = list(gaussian.noise[1:(n/2),], delta*gaussian.noise[-c(1:n/2),])
    #noise = map2(noise, wish, ~.x%*%sqrtm(.y))
    data = gaussian.noise %>% reduce(rbind)
    out = matrix.dist.test.stat(data, 4*p)
    result[[1]] = list(dist="gaussian", value=max(out,na.rm=T), 
                  cpt.est=which.max(out))


	set.seed(353*n + 541*p + 7*i)
    exp.noise = rexp(n*p) %>% matrix(ncol=p)
    exp.noise = list(exp.noise[1:(n/2),]-1, delta*exp.noise[-c(1:n/2),]-1)
    data = exp.noise %>% reduce(rbind)
    out = matrix.dist.test.stat(data, 4*p)
    result[[2]] = list(dist="exponential", value=max(out,na.rm=T), 
                  cpt.est=which.max(out))

    set.seed(353*n + 541*p + 7*i)
    unif.noise = runif(n*p) %>% matrix(ncol=p)
    unif.noise = list(unif.noise[1:(n/2),]-1/2, delta*unif.noise[-c(1:n/2),]-1/2)
    data = unif.noise %>% reduce(rbind)
    out = matrix.dist.test.stat(data, 4*p)
    result[[3]] = list(dist="uniform", value=max(out,na.rm=T), 
                  cpt.est=which.max(out))


    set.seed(353*n + 541*p + 7*i)
    student.noise = rt(n*p,5) %>% matrix(ncol=p)
    student.noise = list(student.noise[1:(n/2),]-1/2, delta*student.noise[-c(1:n/2),]-1/2)
    data = student.noise %>% reduce(rbind)
    out = matrix.dist.test.stat(data, 4*p)
    result[[4]] = list(dist="student", value=max(out,na.rm=T), 
                  cpt.est=which.max(out))


    return(result)
}

#### Old Run
#n = c(1000,2000,5000)
#p = c(3,25,100,10)

n = c(200,500,1000,2000,5000)
p= c(10,50,100)
run = seq(0,9)*100
delta=c(1,1.2)



#p = p[iter%%4 + 1]
#n = n[floor(iter/4)%%3 +1]
parameters = cross(list(n,p,delta,run)) 
n= parameters[[iter]][[1]]
p= parameters[[iter]][[2]]
delta=parameters[[iter]][[3]]
run= parameters[[iter]][[4]]
result=map(seq(1,100), ~run.sim(run+.x, n, p,delta))

string = sprintf("~/Documents/covariance/model_misspecification/Results/model_mis%d.Rds",  iter) 
saveRDS(result,string)
