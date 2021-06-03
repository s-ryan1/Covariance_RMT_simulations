library(cpt.cov)
library(tidyr)
library(tibble)
library(purrr)

args = commandArgs(trailingOnly=TRUE)
iter = as.numeric(args[1])

run.sim = function(i, n, p){
    delta = 1.15
	set.seed(353*n + 541*p + 7*i)
    noise = rnorm(n*p,0,1) %>% matrix(ncol=p) 		
    noise1 = list(noise[1:(n/2),], delta*noise[-c(1:n/2),])
    noise1 = noise1 %>% reduce(rbind) 
    data = list(noise, noise - apply(noise,2,mean), 
                noise1, noise1 - apply(noise1,2,mean)) 
    out = map(data, matrix.dist.test.stat, 4*p) 
    
    values = out %>% map_dbl(max,na.rm=T)
    cpts = out %>% map_int(which.max)
    settings=c("Null-Mean Known", "Null-Mean Unknown", 
         "Alt-Mean Known", "Alt-Mean Unknown")
    out = tibble(n=n,p=p,settings=settings, values=values, cpt.est=cpts) 
    return(out)
}


params = list(c(500,15), c(2500,100))
run = seq(0,9)*100
#
##p = 100
#
n=params[[iter%%2+1]][1]
p=params[[iter%%2+1]][2]
run=run[iter%%10+1]
result=map(seq(1,10), ~run.sim(run+.x, n, p)) %>% tibble %>% unnest
#
#
string = sprintf("~/Documents/covariance/mean_correction/Results/mean_correction%d.Rds", iter) 
saveRDS(result, string)
