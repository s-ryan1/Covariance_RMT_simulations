library(tidyverse)

bonferoni = function(n,alpha){ qnorm(1-.05/n)}


n = c(2000)
p= c(50)
run = seq(0,9)*100
delta = c(1,1.1,1.2)
parameters = cross(list(n,p,delta,run)) 

results = readRDS("~/Documents/covariance/model_misspecification/full_results.Rds")

#results = map2(parameters, results, 
#               ~list(n=.x[[1]], p=.x[[2]], delta=.x[[3]],result=.y))
#results = results %>% discard(~length(.x[[4]][[2]])>0) 
#results = results %>% map(function(X){X[[4]]=X[[4]][[1]];X})
results = map(results, ~.x$result) %>%
    map( map, transpose) %>% 
    modify_depth(3,unlist) %>% 
    map(map,as_tibble) %>% 
    map(reduce,rbind)
results = map2(parameters, results, 
     ~.y %>% mutate(n=.x[[1]], p=.x[[2]], delta=.x[[3]])) %>% 
    reduce(rbind)

null = results %>% filter(delta==1) %>% 
    mutate(threshold=bonferoni(2000,.05))

null %>% group_by(dist) %>% 
    summarise(FPR = mean(value>threshold))

alt = results %>% filter(delta > 1) %>%
    mutate(threshold=bonferoni(2000,.05))

results %>% filter(delta < 1.2) %>% 
    filter(value < 50) %>% 
    ggplot(aes(dist, value,fill=factor(delta))) + geom_boxplot()

alt %>% group_by(dist, delta) %>% 
    summarise(TPR = mean(value>threshold))

alt %>% group_by(dist, delta) %>% 
    filter(value > threshold) %>% 
    summarise(
              cpt_error = mean(abs(cpt.est -1000)))

alt %>% mutate(cpt_error = abs(cpt.est-1000)) %>% 
    mutate(dist = factor(dist, 
                         levels=c("gaussian","uniform", 
                                 "exponential", "student"))) %>% 
    filter(delta==1.1) %>% 
    filter(cpt_error < 300) %>% 
    ggplot(aes(dist, cpt_error)) + 
    geom_boxplot() + 
    scale_colour_grey() + 
    ylab("Changepoint Error") + 
    theme(axis.title.x = element_blank(),
            axis.title = element_text(size=9)) +
    theme_bw()
ggsave("~/Documents/covariance/resubmission/figures/cpt_error_model.pdf", 
width = 3.5, height =1.6)

results %>% mutate(cpt_error = abs(cpt.est-1000)) %>% 
    mutate(dist = factor(dist, 
                         levels=c("gaussian","uniform", 
                                 "exponential", "student"))) %>% 
    filter(delta<1.2) %>% 
#    filter(cpt_error < 300) %>% 
    filter(value < 25) %>%  
    mutate(delta=factor(delta)) %>%
    mutate(delta=fct_recode(delta, null="1", alt="1.1")) %>%
    ggplot(aes(dist, value, fill=delta)) + 
    geom_boxplot() + 
    scale_fill_grey() + 
    ylab("Test Value") +
    geom_hline(yintercept=bonferoni(2000,.95), linetype="dashed") + 
    theme(axis.title.x = element_blank(),
            axis.title = element_text(size=9)) +
    theme_bw()
ggsave("~/Documents/covariance/resubmission/figures/sep_model.pdf", 
width = 4.15, height =1.6)



results %>% filter(n==2000 & p <100) %>% 
    ggplot(aes(result)) + geom_histogram(bins=50) + 
    facet_grid(p~.) + 
    xlab("Test Statistic Value") #+
    theme(axis.text = element_text(size=15))
ggsave("~/Documents/covariance/resubmission/figures/null_value.pdf", width=5, height=4)

FPR=results %>% mutate(threshold=bonferoni(n,.05)) %>% 
    group_by(n,p) %>% 
    summarise(FPR = mean(result>threshold)) %>%ungroup

FPR %>% ggplot(aes(n,FPR, colour=factor(p))) + geom_point()

