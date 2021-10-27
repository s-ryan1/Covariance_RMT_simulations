library(furrr)
library(cpt.cov)
library(purrr)
library(tibble)
library(tidyr)
library(rlang)
#library(tidyverse)
#source("~/Documents/covariance/proj_functions.R")

#plan(multiprocess, gc=TRUE)


args = commandArgs(trailingOnly=TRUE)
iter = as.numeric(args[1])

sim.runner = function(n, p){
	#cpts = c(floor(n/3))
	#set.seed(353*n + 541*p)
	#gamma = p^(3/4)#optimize(function(x){ abs(uniform.distance(x)-p^(-1/2)) }, c(0,1))[[1]]
        #diff.mats = purrr::map(seq(1, 1), ~generate.diff.mat(p))
delta=1.15
        sigma1 = (1/p)*rWishart(1, p, diag(p))[,,1]
	sigma1 = expm::sqrtm(sigma1)
        #diff.mats = append(list(sigma1), diff.mats)
        #sigma = purrr::accumulate(diff.mats, function(X,Y){ B = expm::sqrtm(X); return(B%*%Y%*%B)})

    #noise = map2(noise, wish, ~.x%*%sqrtm(.y))



	function(i){
		#set.seed(152342*i)
		#data = gen.data(n=n, p=p, m=1, minseglen=1, cpts=cpts, sigma=sigma) %>% pluck("data") %>% reduce(rbind)

	set.seed(353*n + 541*p + 7*i)
    noise = rnorm(n*p,0, 1) %>% matrix(ncol=p)
    noise = list(noise[1:floor((n/3)),], delta*noise[-c(1:floor(n/3)),])
    data = noise %>% reduce(rbind)
	data = data%*%sigma1

		
        if(p < 20){
            result = tibble(t = seq(1,n), 
                            n=n,p=p,
                    ratio=matrix.dist.test.stat(data, 0),
                    wang = wang.stat(data,0),
                    galeano = cpt.cov:::galeano.cusum.stat(data, 0),
                    aue = aue.stat(data, 0))
        }else{
            result = tibble(t = seq(1,n), 
                            n=n,p=p,
                    ratio=matrix.dist.test.stat(data, 0),
                    wang = wang.stat(data,0),
                    galeano = cpt.cov:::galeano.cusum.stat(data, 0))
        }
        return(result)
					
	}
}

params = list(c(500,15), c(2500,100))
run = seq(0,9)*100
#
##p = 100
#
n=params[[iter%%2+1]][1]
p=params[[iter%%2+1]][2]
run.sim = sim.runner(n,p)
run=run[iter%%10+1]
result=map(seq(1,10), ~run.sim(run+.x )) %>% tibble %>% unnest


#result=map(run.sim, ~future_map(seq(1,100), .x))
string = sprintf("~/Documents/covariance/single_change/Results/run%d.Rds", iter) 
saveRDS(result,string)
