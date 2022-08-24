library(tidyverse)
library(cowplot)


bonferoni = function(n,alpha){ qnorm(1-.05/n)}

results =readRDS("~/Documents/covariance/autocorrelation/full_results.Rds")

results = results %>% map(~.x$result) %>% map(reduce, rbind) %>% 
    reduce(rbind) %>% as_tibble

results = results %>% mutate(threshold=bonferoni(n,.05))

results %>% filter(delta==1) %>% 
    mutate(detect = test.value > bonferoni(2000,.05)) %>% 
    group_by(phi) %>% 
    summarise(FPR = mean(detect))  %>%
    ungroup #%>%
    ggplot(aes(factor(phi),FPR)) + geom_point() 


results %>% filter(delta==1 & phi < .9) %>% 
    ggplot(aes(factor(phi), test.value)) + geom_boxplot() #+ 
    scale_y_log10()


results %>% filter(delta>1) %>% 
    mutate(cpt_error = abs(cpt.est - n/2)) %>%
    ggplot(aes(factor(phi), cpt_error, colour=factor(delta))) + geom_boxplot() #+ 
    scale_y_log10()

results = map(results, ~.x$result) %>%
    map( map, transpose) %>% 
    modify_depth(3,unlist) %>% 
    map(map,as_tibble) %>% 
    map(reduce,rbind)
results = map2(parameters, results, 
     ~.y %>% mutate(n=.x[[1]], p=.x[[2]], delta=.x[[3]])) %>% 
    reduce(rbind)



df = df %>% group_by(n,p,delta,phi) %>% 
    summarise(test.value = quantile(test.value, c(.95)),
               cpt.error = quantile(abs(cpt.est-n/2), c(.5))) 


df %>% 
    filter(delta==1 & p%in%c(10,25,100) & phi!=.1 & n<5000 & n>100) %>% 
    ggplot(aes(n,test.value)) + geom_point() + 
    facet_grid(phi~p,scales="free_y") + scale_y_log10()

df %>% 
    filter(delta == 2 & p %in% c(10,25,100) & phi!=.1 & n < 5000) %>% 
    ggplot(aes(n,cpt.error)) + geom_point() + 
    facet_grid(phi~p,scales="free_y") #+ scale_y_log10()

df %>% 
    filter(delta == 2 & p %in% c(10,25,100) & phi!=.1 & n < 5000) %>% 
    ggplot(aes(n,test.value)) + geom_point() + 
    facet_grid(phi~p,scales="free_y") #+ scale_y_log10()


df %>% filter(delta == 2 & p %in% c(10,25,100) & phi!=.1 & n < 5000) %>% 
    ggplot(aes(n,cpt.error)) + geom_point() + 
    facet_grid(phi~p,scales="free_y") #+ scale_y_log10()

X = df %>% 
    filter(delta==1 & p %in% c(10,25,100)  & n==1000) %>%
    ungroup %>%
   select(p,phi,test.value) 
X=X$test.value %>% matrix(ncol=3)
colnames(X) = c("p=10","p=25","p=100") 
X=cbind(phi,X)

X = df %>% 
    filter(delta==2 & p %in% c(10,25,100)  & n==1000) %>%
    ungroup %>%
   select(p,phi,cpt.error) 
X=X$cpt.error %>% matrix(ncol=3)
colnames(X) = c("p=10","p=25","p=100") 
X=cbind(phi,X)

results %>% 
    filter(delta < 1.2 & phi < .9) %>% 
    ggplot(aes(factor(phi), test.value, colour=factor(delta))) + geom_boxplot() #+ scale_y_log10()

results %>% 
    filter(delta==1.1 &  phi<.9) %>% 
    filter(abs(cpt.est - n/2)<600) %>% 
    mutate(phi=factor(phi), cpt_error = abs(cpt.est-n/2), 
           delta=factor(delta)) %>%
    ggplot(aes(phi, cpt_error)) + 
    ylab("Changepoint Error") + 
        geom_boxplot() +# scale_y_log10()
    theme(legend.title= element_text(size=9),
            axis.title = element_text(size=9)) + 
    theme_bw()
ggsave("~/Documents/covariance/resubmission/figures/ar_cpt_error.pdf", 
width = 3.5, height =1.8)

simulated_threshold = results %>% filter(delta==1) %>% 
    group_by(phi) %>% 
    summarise(threshold = quantile(test.value, .95)) %>% 
    pull(threshold)

results %>% group_by(phi) %>% 
    nest %>% ungroup %>%  
    mutate(threshold=simulated_threshold) %>% 
    unnest %>% filter(delta>1) %>% 
#    filter(test.value>threshold) %>% 
    mutate(cpt_error = cpt.est - threshold) %>%
    ggplot(aes(factor(phi), cpt.est, colour=delta)) + geom_boxplot() 

results %>% group_by(phi) %>% 
    filter(delta == 1) %>% 
    summarise(bias=median(cpt.est-n/2))

results %>% filter(delta=0)

results %>% filter(delta > 1 & phi==.6) %>% 
    filter(test.value > 35.7) 


%>% 
    ggplot(aes(phi, abs(cpt.est-n/2), fill=delta)) + geom_boxplot()
    #group_by(delta) %>% 
    #summarise() 


results %>% 
    filter(delta<=1.1 &  phi<.9) %>% 
    mutate(phi=factor(phi), cpt_error = abs(cpt.est-n/2), 
           delta=factor(delta)) %>%
    mutate(delta=fct_recode(delta, null="1", alt="1.1")) %>%
    ggplot(aes(phi, test.value, fill=delta)) + 
        geom_boxplot() +
       scale_fill_grey() +
       ylab("Test Value") +
    theme(legend.title= element_text(size=9),
            axis.title = element_text(size=9)) + 
    geom_hline(yintercept=bonferoni(2000,.95), linetype="dashed") +
    theme_bw()
ggsave("~/Documents/covariance/resubmission/figures/ar_sep.pdf", 
width = 4.15, height =1.8)

results %>% filter(delta==1) %>% 
    group_by(phi) %>% 
    summarise(FPR = mean(test.value > threshold)) %>%
    xtable()


results %>% filter(delta==1.1) %>% 
    group_by(phi) %>% 
    summarise(TPR = mean(test.value > threshold)) %>%
    xtable()


results %>% filter(delta==1.2) %>% 
    group_by(phi) %>% 
    summarise(cpt_error = mean(abs(cpt.est-n/2))) %>%
    xtable()
