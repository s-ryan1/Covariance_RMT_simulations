library(cpt.cov)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
iter = as.numeric(args[1])

run.sim = function(i, n, p){

	if(n < 5*p*log(n)){return(NA) }
    sigma=list()
    
    set.seed(1234)
    B = rnorm(10*p,0,1) %>% matrix(ncol=p) %>% cov %>% eigen %>% pluck("vectors")	
    sigma[[1]] = runif(p,.1,10)  
    
    diff=1
    set.seed(541*p +7*i)
    for(k in 2:5){
        sigma[[k]] = map_dbl(seq(1,p), 
                         ~runif(1, max(c(sigma[[k-1]][.x]-diff,0.1)),
                                sigma[[k-1]][.x]+diff))
        b = sample(seq(1,p),1)
        sigma[[k]][b] = sigma[[k-1]][b] + diff 
    }
    sigma = map(sigma, diag)
    sigma = map(sigma, ~B%*%((.x))%*%t(B))

	#print(c(n, floor(p*log(n))))

    cpts =  gen.cpt.locs(n, 4, max(30,floor(p*log(n))))
    lengths= diff(c(0,cpts,n))
    noise = map(lengths, ~rnorm(.x*p, 0, 1)) %>% map(matrix, ncol=p)
    data = map2(noise, sigma, ~(.x)%*%.y) %>% reduce(rbind) 


    result.fisher = bin.seg(data, c(0,n), matrix.dist.test.stat, threshold=qnorm(1-.05/(500^2)), minseglen=40, c())
    fisher = result.fisher %>% bin.seg.to.cpt(qnorm(1-.05/(500^2)))
    wang.thresh = wang.threshold(data)
    result.wang = bin.seg(data, c(0,n), wang.stat, threshold=wang.thresh, minseglen=p*log(n), c())
    wang = result.wang %>% bin.seg.to.cpt(wang.thresh)
    result.galeano = bin.seg(data, c(0,n), cpt.cov:::galeano.cusum.stat, threshold=qnorm(1-.05/(500^2)), minseglen=40, c())
    galeano = result.galeano %>% bin.seg.to.cpt(log(n))
    #print(fisher$cpts)



    fisher.cpt.error = detection.rates(fisher$cpts, cpts, 20)
    wang.cpt.error = detection.rates(wang$cpts, cpts, 20)
    galeano.cpt.error = detection.rates(galeano$cpts, cpts, 20)
	#set.seed(353*n + 541*p + 7*i)
    #return(max(matrix.dist.test.stat(data, 4*p), na.rm=T))

    if(p<=20 & p*(p+1) < 2*n){
        safe.bin.seg = safely(bin.seg)
        result.aue = safe.bin.seg(data, c(0,n), aue.stat, threshold=qnorm(.95), minseglen=p*log(n), c())
        if(length(result.aue[[1]])>0){
            aue = result.aue[[1]] %>% bin.seg.to.cpt(qnorm(.95))
            aue.cpt.error = detection.rates(aue$cpts, cpts, 20)
        }
        models = list(fisher$cpts, wang$cpts, aue$cpts, galeano$cpts)
        rates = map(models, detection.rates, cpts, 20)
        m_hat = map_dbl(rates, ~.x$m)
        TDR = map_dbl(rates, ~.x$TDR)
        FDR = map_dbl(rates, ~.x$FDR)
        mae = map_dbl(models, MAE,  cpts, data) 
        smae = map_dbl(models, spectral.error,  cpts, data) 
        result = tibble(
                        method=c("Ratio", "Wang", "Aue", "Galeano"),
                        TDR=TDR, FDR=FDR, m_hat=m_hat, 
                        MAE=mae, SMAE=smae, n=n, p=p, power=diff)
        result = result %>% gather("metric", "value",-method,-n,-p,-power)
        return(result)
    }

    models = list(fisher$cpts, wang$cpts, galeano$cpts)
print(cpts)
print(models)
    rates = map(models, detection.rates, cpts, 20)
    m_hat = map_dbl(rates, ~.x$m)
    TDR = map_dbl(rates, ~.x$TDR)
    FDR = map_dbl(rates, ~.x$FDR)
    mae = map_dbl(models, MAE,  cpts, data) 
    smae = map_dbl(models, spectral.error,  cpts, data) 
    result = tibble(
                    method=c("Ratio", "Wang", "Galeano"),
                    TDR=TDR, FDR=FDR, m_hat=m_hat, 
                    MAE=mae, SMAE=smae, n=n, p=p, power=diff)
    result = result %>% gather("metric", "value",-method,-n,-p,-power)
    return(result)
}


n = c(500,1000,2000,5000)
p= c(3,10,30,100)
run = seq(0,9)*100

parameters = cross3(n,p,run) 
n= parameters[[iter]][[1]]
p= parameters[[iter]][[2]]
run= parameters[[iter]][[3]]

safe_run = safely(run.sim)
result=map(seq(1,100), ~safe_run(run+.x, n, p))
#
#
string = sprintf("~/Documents/covariance/wang_multi/Results/run%d.Rds", iter) 
saveRDS(result, string)

#result=result %>% map(~.x$result) %>% reduce(rbind) %>% 
#    group_by(method,metric) %>% 
    #summarise(value=mean(value))
