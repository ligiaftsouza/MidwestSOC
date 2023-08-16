library(tidyverse)
library(rgdal)
library(raster)
library(broom)
library(prism)
library(modelr)
library(sp)
library(scales)
library(viridis)
library(patchwork)

#################### 1) Fitting the PEDS profiles to the Nelder-Mead optimizer ####################
### Loading the optimization function
source("NelderMead_optim.R")

### Function to calculate R-squared
rsq <- function(par,d)
{
  soci <- exp(par[1])
  socf <- soci*(atan(par[2])+pi/2)/pi
  beta <- exp(par[3])
  soc <- d$soc
  depth <- d$depth
  
  resids <- soc-socf-(soci-socf)*exp(-beta*depth)
  total <- sum((soc - mean(soc))^2)
  ### calculation of R2
  rsq <- 1 - (sum(resids^2)/total)

  results <- list("rsq" = rsq)
  return(results)
}

### Loading the data
soc <- read_rds("1.Data/soils_midwest.rds")
stock <- read_rds("1.Data/soilsbd_midwest.rds")

### Loop to fit the Nelder-Mead optimizer
for(dat in c("soc", "stock")){
  d <- get(dat)
  allres_const2<-list()
  allconvergence_const2<-c()
  
  if(dat == "stock"){
    d <- d %>% unnest(data) %>% dplyr::select(-soc) %>% rename(soc = cstock) %>% 
      nest(-id, -lat, -long, -EP, -use)
  }
  
  for (counter in 1:(dim(d)[1]))
  {
    print(paste("Working on constrained-2 optimization",counter,"of",dim(d)[1]))
    allres_const2[[counter]]<-nd_optim(d = d$data[[counter]][,c("soc", "depth")], beta0 = 1, mkplot = F)
    
    if (class(allres_const2[[counter]])!="try-error")
    { allconvergence_const2<-c(allconvergence_const2,allres_const2[[counter]]$convergence)
    }else
    {allconvergence_const2<-c(allconvergence_const2,999)
    }
  }
  allconvergence_const2 
  sum(allconvergence_const2!=0) #so they all converged
  
  #***organize the results
  allrestab<-data.frame(soci=NA*numeric(dim(d)[1]),socf=NA,beta=NA,start=NA,end=NA)
  for (counter in 1:dim(d)[1])
  {#tabulate unconstrained parameters results
    par<-allres_const2[[counter]]$par
    
    soci<-exp(par[1])
    socf<-soci*(atan(par[2])+pi/2)/pi
    beta<-exp(par[3])
    
    allrestab[counter,1]<-soci
    allrestab[counter,2]<-socf
    allrestab[counter,3]<-beta
    allrestab[counter,4]<-allres_const2[[counter]]$start_value
    allrestab[counter,5]<-allres_const2[[counter]]$value
  }
  
  ### calculating r-squared
  data <- bind_cols(d, allrestab)
  for(p in c(1:length(allres_const2))){
    par <- allres_const2[[p]]$par
    res <- rsq(par, d$data[p] %>% data.frame)
    data[p, "rsq"] <- res[1]
  }
  
  ### saving the new dataset
  write_rds(data, paste0("2.Output/data", dat, "optimizer_rsq.rds"))
}

rm(list = ls())

############################### 2) Predicting SOC depth distribution ##############################
soc <- read_rds("2.Output/optimsocconcdata.rds")
soc <- soc %>% mutate(mxsoc = map_dbl(data, ~pull(., soc)[1]),
                      socdeep = map_dbl(data, ~pull(., soc)[length(pull(., soc))]),
                      exc = ifelse(mxsoc >= socdeep, "ok", "remove")) %>% 
  filter(exc == "ok")
depth <- seq(1, 500, 1)

soc_preds <- soc %>% 
  unnest(data) %>% nest(-id, -lat, -long, -EP, -use, -soci, -socf, -beta, -beta_se, -rsq) %>% 
  mutate(pred_depth = list(depth)) %>% unnest(pred_depth) %>% 
  mutate(pred_soc = socf + (soci - socf)*exp(-beta*pred_depth),
         pred_depth = ifelse(round(pred_soc, 5) <= round(socf, 5), NA, pred_depth)) %>%
  filter(!is.na(pred_depth)) %>% 
  nest(-id, -EP, -use, -beta, -beta_se, -rsq) %>% 
  mutate(SOCz = map_dbl(data, ~max(pull(., pred_depth), na.rm = T)))

write_rds(soc_preds, "2.Output/socconc_param.rds")


soc_preds <- read_rds("2.Output/socconc_param.rds")
range01 <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}

preds <- soc_preds %>%
  unnest(data) %>% 
  mutate(pdepth = ifelse(pred_depth <= SOCz, pred_depth, NA)) %>% 
  filter(!is.na(pdepth)) %>% 
  nest(-id, -EP, -use) %>% 
  mutate(soc_scale = map(data, ~range01(pull(., pred_soc)))) %>% unnest(soc_scale, data) %>% 
  nest(-id, -EP, -use, -beta, -rsq, -pred_depth, -soc_scale)


stocks <- read_rds("2.Output/optimsocstockdata.rds")
stocks <- stocks %>% mutate(mxsoc = map_dbl(data, ~pull(., soc)[1]),
                            socdeep = map_dbl(data, ~pull(., soc)[length(pull(., soc))]),
                            exc = ifelse(mxsoc >= socdeep, "ok", "remove")) %>% 
  filter(exc == "ok")
depth <- seq(1, 500, 1)

stocks_preds <- stocks %>% 
  unnest(data) %>% nest(-id, -lat, -long, -EP, -use, -soci, -socf, -beta, -beta_se, -rsq) %>% 
  mutate(pred_depth = list(depth)) %>% unnest(pred_depth) %>% 
  mutate(pred_soc = socf + (soci - socf)*exp(-beta*pred_depth),
         pred_depth = ifelse(round(pred_soc, 5) <= round(socf, 5), NA, pred_depth)) %>%
  filter(!is.na(pred_depth)) %>% 
  nest(-id, -EP, -use, -beta, -beta_se, -rsq) %>% 
  mutate(SOCz = map_dbl(data, ~max(pull(., pred_depth), na.rm = T)))

write_rds(stocks_preds, "2.Output/socstocks_param.rds")

range01 <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}

stocks_preds <- read_rds("2.Output/socstocks_param.rds")

preds_st <- stocks_preds %>% unnest(data) %>% 
  mutate(pdepth = ifelse(pred_depth <= SOCz, pred_depth, NA)) %>% 
  filter(!is.na(pdepth)) %>% 
  nest(-id, -EP, -use) %>% 
  mutate(soc_scale = map(data, ~range01(pull(., pred_soc)))) %>% unnest(soc_scale, data) %>% 
  nest(-id, -EP, -use, -beta, -rsq, -pred_depth, -soc_scale)

rm(list = ls())


###################### 3) Fitting the Kansas data to the Nelder Mead optimizer ####################
soc <- rows_update(read_csv("1.Data/soc.csv") %>% mutate(horizon_top = 1),
                   read_csv("1.Data/soc.csv") %>% filter(is.na(Horizon)) %>%  
                     mutate(horizon_top = case_when(MidDepth == 2.5 ~ 0,
                                                    MidDepth == 10 ~ 5,
                                                    MidDepth == 22.5 ~ 15,
                                                    MidDepth == 52.5 ~ 30,
                                                    MidDepth == 112.5 | MidDepth == 97.5 ~ 75)), by = "SampleID") %>% 
  filter(!str_detect(SampleID, "KNZ-A0"))
kspits <- read_csv("1.Data/nrcs.csv") %>% mutate(cstock = soc*BD) %>% rename(Horizon = horizon)

soc <- left_join(soc, kspits %>% dplyr::select(Site, Use, Horizon, BD, cstock, horizon_top),
                 by = c("Site", "Use", "Horizon")) %>% 
  mutate(horizon_top = ifelse(horizon_top.x == 1, horizon_top.y, horizon_top.x)) %>% 
  dplyr::select(-horizon_top.x, -horizon_top.y) %>% 
  rename(depth = horizon_top, soc = SOC) %>% 
  filter(!str_detect(SampleID, "LGN_N_S075|TRG_A_S120")) %>% 
  nest(-Site, -Use, -ep) 

### Loading the optimization function
source("NelderMead_optim.R")

d <- soc
allres_const2<-list()
allconvergence_const2<-c()
for (counter in 1:(dim(d)[1]))
{
  print(paste("Working on constrained-2 optimization",counter,"of",dim(d)[1]))
  allres_const2[[counter]]<-nd_optim(d = d$data[[counter]][,c(13, 4)], beta0 = 1, 
                                     mkplot = F)
  
  if (class(allres_const2[[counter]])!="try-error")
  { allconvergence_const2<-c(allconvergence_const2,allres_const2[[counter]]$convergence)
  }else
  {allconvergence_const2<-c(allconvergence_const2,999)
  }
}
allconvergence_const2 
sum(allconvergence_const2!=0) #so they all converged

#***organize the results
allrestab<-data.frame(soci=NA*numeric(dim(d)[1]),socf=NA,beta=NA,start=NA,end=NA)
for (counter in 1:dim(d)[1])
{#tabulate unconstrained parameters results
  par<-allres_const2[[counter]]$par
  
  soci<-exp(par[1])
  socf<-soci*(atan(par[2])+pi/2)/pi
  beta<-exp(par[3])
  
  allrestab[counter,1]<-soci
  allrestab[counter,2]<-socf
  allrestab[counter,3]<-beta
  allrestab[counter,4]<-allres_const2[[counter]]$start_value
  allrestab[counter,5]<-allres_const2[[counter]]$value
}

rsq <- function(par,d)
{
  soci <- exp(par[1])
  socf <- soci*(atan(par[2])+pi/2)/pi
  beta <- exp(par[3])
  soc <- d$soc
  depth <- d$depth
  
  resids <- soc-socf-(soci-socf)*exp(-beta*depth)
  total <- sum((soc - mean(soc))^2)
  ### calculation of R2
  rsq <- 1 - (sum(resids^2)/total)
  
  results <- list("rsq" = rsq)
  
  return(results)
}

### calculating r-squared
data <- bind_cols(d, allrestab)
for(p in c(1:length(allres_const2))){
  par <- allres_const2[[p]]$par
  res <- rsq(par, d$data[p] %>% data.frame)
  data[p, "rsq"] <- res[1]
  data[p, "beta_se"] <- res[2]
  data[p, "beta_sd"] <- res[3]
}

depths <- seq(0, 500, 1)

socks_preds <- data %>% mutate(pred_depths = list(depths)) %>% unnest(pred_depths) %>% 
  mutate(soc_pred =  socf + (soci - socf)*exp(-beta*pred_depths),
         pred_depth = ifelse(round(soc_pred, 5) <= round(socf, 5), NA, pred_depths)) %>%
  filter(!is.na(pred_depth)) %>% 
  nest(-Site, -Use, -ep, -beta, -beta_se, -rsq) %>% 
  mutate(SOCz = map_dbl(data, ~max(pull(., pred_depth), na.rm = T)))

write_rds(socks_preds, "2.Output/optimsoc_ks.rds")
rm(list = ls())

