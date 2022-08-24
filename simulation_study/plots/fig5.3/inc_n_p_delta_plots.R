library(tidyverse)

bonferoni = function(n,alpha){ qnorm(1-.05/n)}


n = c(200,500,1000,2000,5000)
p= c(10,50,100)
delta=c(1.05,1.1,1.15,1.2)
run = seq(0,9)*100
parameters = cross(list(n,p,delta,run)) 
parameters=parameters[-1]
results = readRDS("~/Documents/covariance/alt/full_results.Rds")
results = results[-1]
results = map2(parameters, results, ~list(n=.x[[1]], p=.x[[2]], delta=.x[[3]],result=.y))
results = results %>% discard(~length(.x[[4]][[2]])>0) 
results = results %>% map(function(X){X[[4]]=X$result$result$result;X})
results = results %>%  map(as_tibble) %>% 
    map(~.x %>% mutate(value = map_dbl(result, ~.x[1]), 
                       cpt.est= map_dbl(result, ~.x[2])) %>% 
    select(-result)) 
results = results %>% reduce(rbind)


results %>% filter(n==2000 & p <100) %>% 
    ggplot(aes(value)) + geom_histogram(bins=50) + 
    facet_grid(p~.) + 
    xlab("Test Statistic Value") #+
    theme(axis.text = element_text(size=15))
ggsave("~/Documents/covariance/resubmission/figures/alt_value.pdf", width=3.5, height=2)

#summary = results %>% mutate(threshold=bonferoni(n,.05)) %>% 
#    group_by(n,p,delta) %>% 
#    summarise(TPR = mean(value>threshold),
#              cpt.error = mean(abs(cpt.est-n/2))) %>%ungroup 
#    
#    %>% 
#   ggplot(aes(n,TPR, colour=factor(delta))) + geom_point(alpha=.8) + 
#    facet_grid(.~p) 
#
#FPR %>% ggplot(aes(n,FPR, colour=factor(p))) + geom_point()
#
#results %>% mutate(cpt.error = abs(cpt.est-n/2)) %>% 
#    group_by(n,p,delta) %>% 
#    summarise(error = quantile(cpt.error,.5)) %>% ungroup %>% 
#    ggplot(aes(p,error,colour=factor(delta))) + geom_point() + 
#    facet_grid(n~.,scales="free_y")
#
#summary %>% filter(delta==1.1 & n > 200) %>% 
#    ggplot(aes(factor(n), TPR, colour=factor(p), group=factor(p))) +
#    geom_line()
#
#summary %>% filter(delta==1.1 & n > 200) %>% 
#    ggplot(aes(factor(n), cpt.error, colour=factor(p), group=factor(p))) +
#    geom_line()


results %>% 
    mutate(threshold=bonferoni(n,.05)) %>% 
    filter(!(n==500 & p==50)) %>% 
    filter(!(n==1000 & p==100)) %>% 
    filter(delta==1.1 & n > 200) %>% 
    filter(value > threshold) %>% 
    filter(n==5000) %>% 
    mutate(cpt.error=abs(cpt.est-n/2),
           p=factor(p)) %>% 
    filter(cpt.error < 100) %>% 
    ggplot(aes(p, cpt.error)) +
    ylab("Changepoint Error") + 
    geom_boxplot() + 
    theme_bw() + 
    theme(axis.title = element_text(size=8)) + 
ggsave("~/Documents/covariance/resubmission/figures/alt_p.pdf", 
    width=3.5,height=1.5)

results %>% 
    mutate(threshold=bonferoni(n,.05)) %>% 
    filter(!(n==500 & p==50)) %>% 
    filter(!(n==1000 & p==100)) %>% 
    filter(delta==1.15 & n > 200) %>% 
    filter(value > threshold) %>% 
    filter(p==50) %>% 
    mutate(cpt.error=abs(cpt.est-n/2),
           n=factor(n)) %>% 
    ggplot(aes(n, cpt.error)) +
    ylab("Changepoint Error") + 
    geom_boxplot() + 
    theme_bw() + 
    theme(axis.title = element_text(size=8))
ggsave("~/Documents/covariance/resubmission/figures/alt_n.pdf", 
    width=3.5,height=1.5)


results %>% 
    mutate(threshold=bonferoni(n,.05)) %>% 
    filter(!(n==500 & p==50)) %>% 
    filter(!(n==1000 & p==100)) %>% 
    filter(p==50 & n > 200) %>% 
    filter(value > threshold) %>% 
    filter(n==2000) %>% 
    mutate(cpt.error=abs(cpt.est-n/2),
           delta=factor(delta)) %>% 
    ggplot(aes(delta, cpt.error)) +
    ylab("Changepoint Error") + 
    geom_boxplot() + 
    theme_bw() +
    theme(axis.title = element_text(size=8))
ggsave("~/Documents/covariance/resubmission/figures/alt_delta.pdf", 
    width=3.5,height=1.5)

#summary %>% filter(delta==1.1 & p == 50) %>%  
#    ggplot(aes(n, TPR)) + geom_point()
