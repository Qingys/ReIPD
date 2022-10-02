library(tidyverse)
library(IPDfromKM)
library(survival)
library(survminer)


##################################################################
###import point
extract_ipd <- function(file0, file1,
                        timeSeq, 
                        numRisk0, numRisk1,
                        events0=NULL, events1=NULL,
                        totalnum0=NULL, totalnum1=NULL,
                        SurvOrCum='Surv',
                        maxy = 100){
   ###import data
   if (SurvOrCum=='Cum') {
 points_0 <- read_delim(file0, col_names = c('time', 'risk'), delim = '\t') %>% 
   mutate(time=as.numeric(time), risk=as.numeric(risk)) %>% 
   mutate(risk=100-risk)
 points_1 <- read_delim(file1, col_names = c('time', 'risk'), delim = '\t') %>% 
   mutate(time=as.numeric(time), risk=as.numeric(risk)) %>% 
   mutate(risk=100-risk) } else {
      points_0 <- read_delim(file0, col_names = c('time', 'surv'), delim = '\t') %>% 
         mutate(time=as.numeric(time), surv=as.numeric(surv))
      points_1 <- read_delim(file1, col_names = c('time', 'surv'), delim = '\t') %>% 
         mutate(time=as.numeric(time), surv=as.numeric(surv))
   }
 ###preprocess
 preprocess_0 <- preprocess(dat = points_0, nrisk = numRisk0, trisk = timeSeq,
                          maxy = maxy, totalpts = totalnum0)
 preprocess_1 <- preprocess(dat = points_1, nrisk = numRisk1, trisk = timeSeq,
                            maxy = maxy, totalpts = totalnum1)
 ###get data
 ipd_0 <- getIPD(preprocess_0, armID = 0, tot.events = events0)
 ipd_1 <- getIPD(preprocess_1, armID = 1, tot.events = events1)
 ipd <- bind_rows(ipd_0$IPD, ipd_1$IPD)
 print('------------------------Group==0----------------------------------------')
 summary(ipd_0)
 print('------------------------Group==1----------------------------------------')
 summary(ipd_1)
 print('------------------------validation by cox model-------------------------')
 ###Cox model
 survreport(ipd1 = ipd_0$IPD, ipd2 = ipd_1$IPD, arms = 2)
 coxph(Surv(time, status) ~ treat, data = ipd) %>% summary() %>% print()
 ###
 km <- survfit(Surv(time, status) ~ treat, data = ipd)
 if (is.numeric(timeSeq)) {
    ### have time seq
    if (SurvOrCum=='Cum') {
       p <- ggsurvplot(km, 
                       fun = 'event',
                       palette = 'lancet',
                       data = ipd,
                       risk.table = T,
                       censor=F,
                       surv.scale='percent',
                       xlim=c(0, timeSeq[length(timeSeq)]),
                       ylim=c(0, max(c(100-points_0$risk, 100-points_1$risk))/100),
                       break.x.by = timeSeq[2])
    } else {
       p <- ggsurvplot(km, 
                       fun = 'pct',
                       palette = 'lancet',
                       data = ipd,
                       risk.table = T,
                       censor=F,
                       surv.scale='percent',
                       xlim=c(0, timeSeq[length(timeSeq)]),
                       ylim=c(min(c(points_0$surv, points_1$surv)), 100),
                       break.x.by = timeSeq[2])
    }
 } else {
    ###not have time seq
 if (SurvOrCum=='Cum') {
 p <- ggsurvplot(km, 
                 fun = 'event',
                 palette = 'lancet',
                 surv.scale='percent',
                 data = ipd,
                 risk.table = T,
                 censor=F)
 } else {
    p <- ggsurvplot(km, 
                    fun = 'pct',
                    palette = 'lancet',
                    surv.scale='percent',
                    data = ipd,
                    risk.table = T,
                    censor=F)
 }
}
 tiff('validated-KM.tiff', res=300, width = 300*8, height = 300*6)
 print(p)
 dev.off()
 return(ipd)
}



######test
ipd <- extract_ipd('file1.txt', 'file2.txt',
                   timeSeq = c(0,12,24,36,48,60),
                   numRisk0 = c(1456,1190,1041,916,757,294),
                   numRisk1 = c(1454,1186,1050,917,750,302),
                   totalnum0 =1456, totalnum1 =1454,
                   SurvOrCum = 'Cum',
                   events0 =92, events1 =59)
write.csv(ipd, 'ipd.csv')






################################covariates################################
library(tidyverse)
library(survival)
ipd <- read_csv('ipd.csv')
covariates <- read_csv('covariate.csv')

###search function
reconstruct_Covs <- function(data, covariates){
   ###treat==0, status==1
   #sampling proportion
   sample_prop_1 <- covariates %>%
      filter(treat==0, status==1) %>% 
      summarise(proportions(n)) %>% pull(`proportions(n)`)
   index_1 <- data %>% 
      filter(treat==0, status==1) %>% 
      sample_frac(sample_prop_1[1], replace = F) %>% 
      pull(X1)
   ###treat==0, status==0
   #sampling proportion
   sample_prop_2 <- covariates %>%
      filter(treat==0, status==0) %>% 
      summarise(proportions(n)) %>% pull(`proportions(n)`)
   index_2 <- data %>% 
      filter(treat==0, status==0) %>% 
      sample_frac(sample_prop_1[1], replace = F) %>% 
      pull(X1)
   ###treat==1, status==1
   #sampling proportion
   sample_prop_3 <- covariates %>%
      filter(treat==1, status==1) %>% 
      summarise(proportions(n)) %>% pull(`proportions(n)`)
   index_3 <- data %>% 
      filter(treat==1, status==1) %>% 
      sample_frac(sample_prop_1[1], replace = F) %>% 
      pull(X1)
   ###treat==1, status==0
   #sampling proportion
   sample_prop_4 <- covariates %>%
      filter(treat==1, status==0) %>% 
      summarise(proportions(n)) %>% pull(`proportions(n)`)
   index_4 <- data %>% 
      filter(treat==1, status==0) %>% 
      sample_frac(sample_prop_1[1], replace = F) %>% 
      pull(X1)
   ###Cov naming
   data <- data %>% mutate(Cov1=case_when(
      X1 %in% index_1 ~ 1,
      X1 %in% index_2 ~ 1,
      X1 %in% index_3 ~ 1,
      X1 %in% index_4 ~ 1,
      TRUE ~ 0
   ))
   ###run cox model
   coxm_1 <- coxph(Surv(time, status) ~ treat, subset =(Cov1==1), data = data)
   coxm_0 <- coxph(Surv(time, status) ~ treat, subset =(Cov1==0), data = data)
   ###derive log HR 
   logHR_1 <- coxm_1$coefficients   
   logHR_0 <- coxm_0$coefficients
   ###import published log HR
   value_names <- table(covariates$value) %>% names()
   logHR_Pub_1 <- covariates %>% filter(value==value_names[2]) %>% pull(logest)
   logHR_Pub_0 <- covariates %>% filter(value==value_names[1]) %>% pull(logest)
   ###calculate
   MAE <- abs(logHR_1-logHR_Pub_1[1])+abs(logHR_0-logHR_Pub_0[1])
   names(MAE) <- 'mean absolute error'
   return(list(re_ipd=data,
               MAE=MAE))
}

###test
reconstruct_Covs(ipd, covariates)

###searching
i <- 0
MAE <- 1
###search condition: MAE < 0.01
while(MAE>0.01){
   m1 <- reconstruct_Covs(ipd, covariates)
   MAE <- m1$MAE
   print(MAE)
   i <- i+1
   print(i)
   if (i==5000) break
}

###validation
m1$re_ipd
coxph(Surv(time, status) ~ treat, subset = (Cov1==1), data = m1$re_ipd) %>% summary()
coxph(Surv(time, status) ~ treat, subset = (Cov1==0), data = m1$re_ipd) %>% summary()


