library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(broom)
library(scales)
library(reghelper)
library(modelr)
library(quantreg)

#################### 1) PEDS stats ####################
mollisols <- read_rds("2.Output/optimsocconcdata.rds") %>% 
  filter(rsq >= 0.8, beta < 1000) %>% 
  mutate(k = log(beta)) %>% 
  unnest(data) %>% 
  nest(-id, -EP, -use, -beta, -beta_se, -k, -rsq)

## Full model for SOC concentrations
r1 <- rq(k ~ use*EP, mollisols %>% 
           transform(use = factor(use, levels = c("A", "N"))), tau = 0.5)
summary.rq(r1, iid = T, se = "iid")


mollisolsbd <- read_rds("2.Output/optimsocstockdata.rds") %>% 
  filter(rsq >= 0.8, beta < 1000) %>% 
  mutate(k = log(beta)) %>% 
  unnest(data) %>% 
  nest(-id, -EP, -use, -beta, -beta_se, -k, -rsq)

## Full model for SOC stocks
r1 <- rq(k ~ use*EP, mollisolsbd %>% 
           transform(use = factor(use, levels = c("A", "N"))), tau = 0.5)
summary.rq(r1, iid = T, se = "iid")


eval <- left_join(mollisolsbd %>% dplyr::select(-data),
                  mollisols %>% dplyr::select(-data), by = c("id", "EP", "use"),
                  suffix = c("stock", "conc")) %>% filter(rsqconc >= 0.8)

### Model evaluation
## comparing betas from SOC concentrations and SOC stocks
cor.test(eval$kconc, eval$kstock, alternative = "two.sided",  method = "spearman", exact = F)


### T-test for the differences in Zsoc
zconc <- read_rds("2.Output/socconc_param.rds") %>% dplyr::select(-data)

source("ResampleMeanDiff.R")

restconc <- ResampMeanDiff(zconc$SOCz[zconc$use == "N" & zconc$rsq >= 0.8], 
                           zconc$SOCz[zconc$use == "A" & zconc$rsq >= 0.8], n = 10000)


zbd <- read_rds("2.Output/socstocks_param.rds") %>% dplyr::select(-data)


restbd <- ResampMeanDiff(zbd$SOCz[zbd$use == "N" & zbd$rsq >= 0.8], 
                         zbd$SOCz[zbd$use == "A" & zbd$rsq >= 0.8], n = 10000)

resamres <- data.frame(bind_rows(c(class = "conc", use = "N", avg = restconc$mnx),
                                 c(class = "conc", use = "A", avg = restconc$mny),
                                 c(class = "stock", use = "N", avg = restbd$mnx),
                                 c(class = "stock", use = "A", avg = restbd$mny)))

write_csv(resamres, "2.Output/resamplingttest_peds.csv")
