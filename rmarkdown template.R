---
  title: "Function-Valued Trait Approach to Evaluating Photosynthesis"
author: "Rebekah Davis and Eric Goolsby"
date: "2023-12-17"
output: html_document
---
  ## Set R Library, Load Data, and Set Environment
  Set R library, load data, and results from long analysis 
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(xlsx)
library(readxl)
library(tidyverse)
library(dplyr) 
source("R/load_data.R")
source("R/mv.R")
load("results/res.RData")
```

## Equations to Model Photosynthesis
Store equations 1-6, 8, 9, 11 from Lobo et al. 2013
```{r}
#Defining light levels based on measurements (PARi) and based off a continuous range (PARi_fine).
PARi <- c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500)
PARi_fine <- seq(0,6000,by = .1)

#Equation 1- rectangular hyperbola Michaelis-Menten based model
eq1 <- function(pars = 
                  c(Pgmax = 19.5,phi_I0 = .0899,Rd = 1.8),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1]
)
{
  if(!missing(dat))
  {
    PARi <- dat$PARi
    if(any(PARi==0)) if(any(!is.na(dat$A[PARi==0]))) pars[["Rd"]] <- abs(dat$A[which(!is.na(dat$A[PARi==0]))][1]) 
    pars[["Pgmax"]] <- max(dat$A,na.rm = TRUE) + pars[["Rd"]]
    pars <- optim(par = pars,function(pars) sum(((dat$A) - (((pars[["phi_I0"]]*dat$PARi*pars[["Pgmax"]])/(pars[["phi_I0"]]*dat$PARi+pars[["Pgmax"]]))-pars[["Rd"]]))^2),control=list(fnscale=1))$par
    try({ pars <- optim(par = pars,function(pars) sum(((dat$A) - (((pars[["phi_I0"]]*dat$PARi*pars[["Pgmax"]])/(pars[["phi_I0"]]*dat$PARi+pars[["Pgmax"]]))-pars[["Rd"]]))^2),control=list(fnscale=1),method = "BFGS")$par },silent = TRUE)
    
    try({ 
      mod <- nls(control = nls.control(maxiter = 1000),A~((phi_I0*PARi*Pgmax)/(phi_I0*PARi+Pgmax))-Rd,data = dat,start = pars)
      pars <- coef(mod)
    },silent = TRUE)
  }
  # parameters to optimize
  Pgmax <- pars[["Pgmax"]]
  phi_I0 <- pars[["phi_I0"]]
  Rd <- pars[["Rd"]]
  
  # predicted values
  pred <- ((phi_I0*PARi*Pgmax)/(phi_I0*PARi+Pgmax))-Rd
  if(return == "predict") return(pred)
  
  # calculated quantities
  Icomp <- Rd*Pgmax/((phi_I0*(Pgmax-Rd)))
  x <- c(.05,.1,.15,.25,.5,.75,.85,.9,.95)
  PARi_fine <- seq(0,6000,by = .1)
  pred_fine <- ((phi_I0*PARi_fine*Pgmax)/(phi_I0*PARi_fine+Pgmax))-Rd
  
  Isat_x <- (Rd*Pgmax*(x-1)-x*(Pgmax^2))/(phi_I0*(Pgmax*(x-1)+Rd*(1-x)))
  
  lo_inds <- which(PARi_fine <= (max(PARi_fine)-50))
  hi_inds <- lo_inds + (length(PARi_fine) - max(lo_inds))
  diffs <- (pred_fine[hi_inds] - pred_fine[lo_inds])
  if(any(diffs <= .067)) Imax <- PARi_fine[which(diffs <= .067)[1]] else Imax <- max(PARi_fine)
  
  Pmax <- max(pred_fine)
  P_Imax <- pred_fine[which.min(abs(Imax - PARi_fine))[1]]
  phi_Icomp <- (phi_I0*(Pgmax^2))/((phi_I0*Icomp+Pgmax)^2)
  I0_Icomp <- if(Icomp <= 0) I0_Icomp <- 1:2 else which(PARi_fine <= Icomp)   
  Icomp_I200 <- which(PARi_fine >= Icomp & PARi_fine <= 200)
  
  phi_I0_Icomp <- as.numeric(coef(lm(pred_fine[I0_Icomp] ~ PARi_fine[I0_Icomp]))[2])
  phi_Icomp_I200 <- as.numeric(coef(lm(pred_fine[Icomp_I200] ~ PARi_fine[Icomp_I200]))[2])
  
  
  P2500 <- pred_fine[PARi_fine == 2500]
  
  Ix <- sapply(x,function(x){
    ind_x <- which((pred_fine - (x*P2500))>0 & PARi_fine < PARi_fine[which.max(pred_fine[PARi_fine <= 2500])])[1]
    if(is.na(ind_x)) return(NA)
    ind_x <- c((ind_x-1):(ind_x+1))
    Ix_mod <- lm(y~x,data.frame(x=PARi_fine[ind_x],y=pred_fine[ind_x]))
    (((x*P2500)-coef(Ix_mod)[[1]])/coef(Ix_mod)[[2]])
  })
  
  # Store calculated parameters  
  calc <- c(Pgmax = Pgmax,
            phi_I0 = phi_I0,
            Rd = Rd,
            Icomp = Icomp,
            Isat_25 = Isat_x[x==.25],Isat_50 = Isat_x[x==.5],
            Isat_75 = Isat_x[x==.75],Isat_85 = Isat_x[x==.85],
            Isat_90 = Isat_x[x==.9],
            Isat_95 = Isat_x[x==.95],
            Psat_25 = .25*(Pgmax-Rd),Psat_50 = .5*(Pgmax-Rd),
            Psat_75 = .75*(Pgmax-Rd),Psat_85 = .85*(Pgmax-Rd),
            Psat_90 = .90*(Pgmax-Rd),
            Psat_95 = .95*(Pgmax-Rd),
            I5 = Ix[x==.05],I10 = Ix[x==.1],I15 = Ix[x==.15],
            I25 = Ix[x==.25],I50 = Ix[x==.5],
            I75 = Ix[x==.75],I85 = Ix[x==.85],
            I90 = Ix[x==.9],
            I95 = Ix[x==.95],
            P25 = .25*Pgmax-Rd,P50 = .5*Pgmax-Rd,
            P75 = .75*Pgmax-Rd,P85 = .85*Pgmax-Rd,
            P90 = .90*Pgmax-Rd,
            P95 = .95*Pgmax-Rd,
            Imax = Imax,
            Pmax = Pmax,
            phi_Icomp = phi_Icomp,
            phi_I0_Icomp = phi_I0_Icomp,
            phi_Icomp_I200 = phi_Icomp_I200,
            P_Imax = P_Imax,
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500])  
  if(return == "calc") return(calc)
  
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}

# Equation 2- rectangular hyperbola Michaelis-Menten based model
eq2 <- function(pars = 
                  c(Pgmax = 19.5,I50 = 216.4,Rd = 1.8),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1]
)
{
  if(!missing(dat))
  {
    PARi <- dat$PARi
    if(any(PARi==0)) if(any(!is.na(dat$A[PARi==0]))) pars[["Rd"]] <- abs(dat$A[which(!is.na(dat$A[PARi==0]))][1]) 
    pars[["Pgmax"]] <- max(dat$A,na.rm = TRUE) + pars[["Rd"]]
    
    pars <- optim(par = pars,function(pars) sum(((dat$A) - (((pars[["Pgmax"]]*dat$PARi)/(dat$PARi+pars[["I50"]]))-pars[["Rd"]]))^2),control=list(fnscale=1))$par
    try({
      pars <- optim(par = pars,function(pars) sum(((dat$A) - (((pars[["Pgmax"]]*dat$PARi)/(dat$PARi+pars[["I50"]]))-pars[["Rd"]]))^2),control=list(fnscale=1))$par
    },silent = TRUE)
    try({
      mod <- nls(control = nls.control(maxiter = 1000),A~((Pgmax*PARi)/(PARi+I50))-Rd,data = dat,start = pars)
      pars <- coef(mod)
    },silent = TRUE)
  }
  
  # parameters to optimize
  Pgmax <- pars[["Pgmax"]]
  I50 <- pars[["I50"]]
  Rd <- pars[["Rd"]]
  
  # predicted values
  pred <- ((Pgmax*PARi)/(PARi+I50))-Rd
  if(return == "predict") return(pred)
  
  phi_I0 <- (I50*Pgmax)/((0+I50)^2)
  
  # calculated quantities
  Icomp <- (I50*Rd)/(Pgmax-Rd)
  x <- c(.05,.1,.15,.25,.5,.75,.85,.9,.95)
  PARi_fine <- seq(0,6000,by = .1)
  pred_fine <- ((Pgmax*PARi_fine)/(PARi_fine+I50))-Rd
  Isat_x <- (I50*(x*Rd-x*Pgmax-Rd))/(Pgmax*(x-1)+Rd*(1-x))
  lo_inds <- which(PARi_fine <= (max(PARi_fine)-50))
  hi_inds <- lo_inds + (length(PARi_fine) - max(lo_inds))
  
  diffs <- (pred_fine[hi_inds] - pred_fine[lo_inds])
  if(any(diffs <= .067)) Imax <- PARi_fine[which(diffs <= .067)[1]] else Imax <- max(PARi_fine)
  
  Pmax <- max(pred_fine)
  P_Imax <- pred_fine[which.min(abs(Imax -PARi_fine))[1]]
  phi_Icomp <-(I50*Pgmax)/((Icomp+I50)^2)
  I0_Icomp <- if(Icomp <= 0) I0_Icomp <- 1:2 else which(PARi_fine <= Icomp)   
  Icomp_I200 <- which(PARi_fine >= Icomp & PARi_fine <= 200)
  
  phi_I0_Icomp <- as.numeric(coef(lm(pred_fine[I0_Icomp] ~ PARi_fine[I0_Icomp]))[2])
  phi_Icomp_I200 <- as.numeric(coef(lm(pred_fine[Icomp_I200] ~ PARi_fine[Icomp_I200]))[2])
  
  P2500 <- pred_fine[PARi_fine == 2500]
  
  Ix <- sapply(x,function(x){
    ind_x <- which((pred_fine - (x*P2500))>0 & PARi_fine < PARi_fine[which.max(pred_fine[PARi_fine <= 2500])])[1]
    if(is.na(ind_x)) return(NA)
    ind_x <- c((ind_x-1):(ind_x+1))
    Ix_mod <- lm(y~x,data.frame(x=PARi_fine[ind_x],y=pred_fine[ind_x]))
    (((x*P2500)-coef(Ix_mod)[[1]])/coef(Ix_mod)[[2]])
  })
  
  calc <- c(Pgmax = Pgmax,
            phi_I0 = phi_I0,
            Rd = Rd,
            Icomp = Icomp,
            Isat_25 = Isat_x[x==.25],Isat_50 = Isat_x[x==.5],
            Isat_75 = Isat_x[x==.75],Isat_85 = Isat_x[x==.85],
            Isat_90 = Isat_x[x==.9],
            Isat_95 = Isat_x[x==.95],
            Psat_25 = .25*(Pgmax-Rd),Psat_50 = .5*(Pgmax-Rd),
            Psat_75 = .75*(Pgmax-Rd),Psat_85 = .85*(Pgmax-Rd),
            Psat_90 = .90*(Pgmax-Rd),
            Psat_95 = .95*(Pgmax-Rd),
            I5 = Ix[x==.05],I10 = Ix[x==.1],I15 = Ix[x==.15],
            I25 = Ix[x==.25],I50 = Ix[x==.5],
            I75 = Ix[x==.75],I85 = Ix[x==.85],
            I90 = Ix[x==.9],
            I95 = Ix[x==.95],
            P25 = .25*Pgmax-Rd,P50 = .5*Pgmax-Rd,
            P75 = .75*Pgmax-Rd,P85 = .85*Pgmax-Rd,
            P90 = .90*Pgmax-Rd,
            P95 = .95*Pgmax-Rd,
            Imax = Imax,
            Pmax = Pmax,
            phi_Icomp = phi_Icomp,
            phi_I0_Icomp = phi_I0_Icomp,
            phi_Icomp_I200 = phi_Icomp_I200,
            P_Imax = P_Imax,
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500])
  
  if(return == "calc") return(calc)
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}

# Equation 3- rectangular hyperbola Michaelis-Menten based model
eq3 <- function(pars = 
                  c(Pgmax = 19.5,phi_I0 = .0493,Rd = 1.8),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1]
)
{
  if(!missing(dat))
  {
    PARi <- dat$PARi
    if(any(PARi==0)) if(any(!is.na(dat$A[PARi==0]))) pars[["Rd"]] <- abs(dat$A[which(!is.na(dat$A[PARi==0]))][1]) 
    pars[["Pgmax"]] <- max(dat$A,na.rm = TRUE) + pars[["Rd"]]
    
    pars <- optim(par = pars,function(pars) sum(((dat$A) - ((pars[["phi_I0"]]*dat$PARi*pars[["Pgmax"]])/(((pars[["Pgmax"]]^2)+(pars[["phi_I0"]]^2)*(dat$PARi^2))^0.5)-pars[["Rd"]]))^2),control=list(fnscale=1))$par
    try({
      pars <- optim(par = pars,function(pars) sum(((dat$A) - ((pars[["phi_I0"]]*dat$PARi*pars[["Pgmax"]])/(((pars[["Pgmax"]]^2)+(pars[["phi_I0"]]^2)*(dat$PARi^2))^0.5)-pars[["Rd"]]))^2),control=list(fnscale=1))$par
    },silent = TRUE)
    
    try({
      mod <- nls(control = nls.control(maxiter = 1000),A~(phi_I0*PARi*Pgmax)/(((Pgmax^2)+(phi_I0^2)*(PARi^2))^0.5)-Rd,data = dat,start = pars)
      pars <- coef(mod)
    },silent = TRUE)
  }
  
  # parameters to optimize
  Pgmax <- pars[["Pgmax"]]
  phi_I0 <- pars[["phi_I0"]]
  Rd <- pars[["Rd"]]
  
  # predicted values
  pred <- (phi_I0*PARi*Pgmax)/(((Pgmax^2)+(phi_I0^2)*(PARi^2))^0.5)-Rd
  if(return == "predict") return(pred)
  
  # calculated quantities
  Icomp <- (Rd*Pgmax/phi_I0)*((1/((Pgmax^2)-(Rd^2)))^0.5)
  x <- c(.05,.1,.15,.25,.5,.75,.85,.9,.95)
  PARi_fine <- seq(0,6000,by = .1)
  pred_fine <- (phi_I0*PARi_fine*Pgmax)/(((Pgmax^2)+(phi_I0^2)*(PARi_fine^2))^0.5)-Rd
  Isat_x <- ((((Pgmax-Rd)*x)+Rd)*Pgmax)/((((Pgmax^2)-((((Pgmax-Rd)*x)+Rd)^2))^0.5)*phi_I0)
  phi_Icomp <- (phi_I0*(Pgmax^3))/(((phi_I0^2)*(Icomp^2)+(Pgmax^2))^(3/2))
  lo_inds <- which(PARi_fine <= (max(PARi_fine)-50))
  hi_inds <- lo_inds + (length(PARi_fine) - max(lo_inds))
  
  diffs <- (pred_fine[hi_inds] - pred_fine[lo_inds])
  if(any(diffs <= .067)) Imax <- PARi_fine[which(diffs <= .067)[1]] else Imax <- max(PARi_fine)
  
  Pmax <- max(pred_fine)
  P_Imax <- pred_fine[which.min(abs(Imax -PARi_fine))[1]]
  I0_Icomp <- if(Icomp <= 0) I0_Icomp <- 1:2 else which(PARi_fine <= Icomp)   
  Icomp_I200 <- which(PARi_fine >= Icomp & PARi_fine <= 200)
  
  phi_I0_Icomp <- as.numeric(coef(lm(pred_fine[I0_Icomp] ~ PARi_fine[I0_Icomp]))[2])
  phi_Icomp_I200 <- as.numeric(coef(lm(pred_fine[Icomp_I200] ~ PARi_fine[Icomp_I200]))[2])
  
  P2500 <- pred_fine[PARi_fine == 2500]
  
  Ix <- sapply(x,function(x){
    ind_x <- which((pred_fine - (x*P2500))>0 & PARi_fine < PARi_fine[which.max(pred_fine[PARi_fine <= 2500])])[1]
    if(is.na(ind_x)) return(NA)
    ind_x <- c((ind_x-1):(ind_x+1))
    Ix_mod <- lm(y~x,data.frame(x=PARi_fine[ind_x],y=pred_fine[ind_x]))
    (((x*P2500)-coef(Ix_mod)[[1]])/coef(Ix_mod)[[2]])
  })
  
  calc <- c(Pgmax = Pgmax,
            phi_I0 = phi_I0,
            Rd = Rd,
            Icomp = Icomp,
            Isat_25 = Isat_x[x==.25],Isat_50 = Isat_x[x==.5],
            Isat_75 = Isat_x[x==.75],Isat_85 = Isat_x[x==.85],
            Isat_90 = Isat_x[x==.9],
            Isat_95 = Isat_x[x==.95],
            Psat_25 = .25*(Pgmax-Rd),Psat_50 = .5*(Pgmax-Rd),
            Psat_75 = .75*(Pgmax-Rd),Psat_85 = .85*(Pgmax-Rd),
            Psat_90 = .90*(Pgmax-Rd),
            Psat_95 = .95*(Pgmax-Rd),
            I5 = Ix[x==.05],I10 = Ix[x==.1],I15 = Ix[x==.15],
            I25 = Ix[x==.25],I50 = Ix[x==.5],
            I75 = Ix[x==.75],I85 = Ix[x==.85],
            I90 = Ix[x==.9],
            I95 = Ix[x==.95],
            P25 = .25*Pgmax-Rd,P50 = .5*Pgmax-Rd,
            P75 = .75*Pgmax-Rd,P85 = .85*Pgmax-Rd,
            P90 = .90*Pgmax-Rd,
            P95 = .95*Pgmax-Rd,
            Imax = Imax,
            Pmax = Pmax,
            phi_Icomp = phi_Icomp,
            phi_I0_Icomp = phi_I0_Icomp,
            phi_Icomp_I200 = phi_Icomp_I200,
            P_Imax = P_Imax,
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500])
  
  if(return == "calc") return(calc)
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}

# Equation 4- hyperbolic tangent based model
eq4 <- function(pars = 
                  c(Pgmax = 19.5,phi_I0 = .0493,Rd = 1.8),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1]
)
{
  if(!missing(dat))
  {
    PARi <- dat$PARi
    if(any(PARi==0)) if(any(!is.na(dat$A[PARi==0]))) pars[["Rd"]] <- abs(dat$A[which(!is.na(dat$A[PARi==0]))][1]) 
    pars[["Pgmax"]] <- max(dat$A,na.rm = TRUE) + pars[["Rd"]]
    pars <- optim(par = pars,function(pars) sum(((dat$A) - ((pars[["Pgmax"]]*tanh(pars[["phi_I0"]]*dat$PARi/pars[["Pgmax"]])-pars[["Rd"]])))^2),control=list(fnscale=1))$par
    try({
      pars <- optim(par = pars,function(pars) sum(((dat$A) - ((pars[["Pgmax"]]*tanh(pars[["phi_I0"]]*dat$PARi/pars[["Pgmax"]])-pars[["Rd"]])))^2),control=list(fnscale=1))$par
    },silent = TRUE)
    try({
      mod <- nls(control = nls.control(maxiter = 1000),A~Pgmax*tanh(phi_I0*PARi/Pgmax)-Rd,data = dat,start = pars)
      pars <- coef(mod)
    },silent = TRUE)
  }
  
  # parameters to optimize
  Pgmax <- pars[["Pgmax"]]
  phi_I0 <- pars[["phi_I0"]]
  Rd <- pars[["Rd"]]
  
  pred <- Pgmax*tanh(phi_I0*PARi/Pgmax)-Rd
  Icomp <- atanh(Rd/Pgmax)*Pgmax/phi_I0
  if(Icomp < 0) Icomp <- 0
  # predicted values
  if(return == "predict") return(pred)
  
  # calculated quantities
  x <- c(.05,.1,.15,.25,.5,.75,.85,.9,.95)
  PARi_fine <- seq(0,6000,by = .1)
  pred_fine <- Pgmax*tanh(phi_I0*PARi_fine/Pgmax)-Rd
  Isat_x <- atanh((((Pgmax-Rd)*x)+Rd)/Pgmax)*Pgmax/phi_I0
  phi_Icomp <- phi_I0*(1/cosh(phi_I0*Icomp/Pgmax))^2
  lo_inds <- which(PARi_fine <= (max(PARi_fine)-50))
  hi_inds <- lo_inds + (length(PARi_fine) - max(lo_inds))
  
  diffs <- (pred_fine[hi_inds] - pred_fine[lo_inds])
  if(any(diffs <= .067)) Imax <- PARi_fine[which(diffs <= .067)[1]] else Imax <- max(PARi_fine)
  
  Pmax <- max(pred_fine)
  P_Imax <- pred_fine[which.min(abs(Imax -PARi_fine))[1]]
  I0_Icomp <- if(Icomp <= 0) I0_Icomp <- 1:2 else which(PARi_fine <= Icomp)   
  Icomp_I200 <- which(PARi_fine >= Icomp & PARi_fine <= 200)
  
  
  phi_I0_Icomp <- as.numeric(coef(lm(pred_fine[I0_Icomp] ~ PARi_fine[I0_Icomp]))[2])
  phi_Icomp_I200 <- as.numeric(coef(lm(pred_fine[Icomp_I200] ~ PARi_fine[Icomp_I200]))[2])
  
  P2500 <- pred_fine[PARi_fine == 2500]
  
  Ix <- sapply(x,function(x){
    ind_x <- which((pred_fine - (x*P2500))>0 & PARi_fine < PARi_fine[which.max(pred_fine[PARi_fine <= 2500])])[1]
    if(is.na(ind_x)) return(NA)
    ind_x <- c((ind_x-1):(ind_x+1))
    Ix_mod <- lm(y~x,data.frame(x=PARi_fine[ind_x],y=pred_fine[ind_x]))
    (((x*P2500)-coef(Ix_mod)[[1]])/coef(Ix_mod)[[2]])
  })
  
  calc <- c(Pgmax = Pgmax,
            phi_I0 = phi_I0,
            Rd = Rd,
            Icomp = Icomp,
            Isat_25 = Isat_x[x==.25],Isat_50 = Isat_x[x==.5],
            Isat_75 = Isat_x[x==.75],Isat_85 = Isat_x[x==.85],
            Isat_90 = Isat_x[x==.9],
            Isat_95 = Isat_x[x==.95],
            Psat_25 = .25*(Pgmax-Rd),Psat_50 = .5*(Pgmax-Rd),
            Psat_75 = .75*(Pgmax-Rd),Psat_85 = .85*(Pgmax-Rd),
            Psat_90 = .90*(Pgmax-Rd),
            Psat_95 = .95*(Pgmax-Rd),
            I5 = Ix[x==.05],I10 = Ix[x==.1],I15 = Ix[x==.15],
            I25 = Ix[x==.25],I50 = Ix[x==.5],
            I75 = Ix[x==.75],I85 = Ix[x==.85],
            I90 = Ix[x==.9],
            I95 = Ix[x==.95],
            P25 = .25*Pgmax-Rd,P50 = .5*Pgmax-Rd,
            P75 = .75*Pgmax-Rd,P85 = .85*Pgmax-Rd,
            P90 = .90*Pgmax-Rd,
            P95 = .95*Pgmax-Rd,
            Imax = Imax,
            Pmax = Pmax,
            phi_Icomp = phi_Icomp,
            phi_I0_Icomp = phi_I0_Icomp,
            phi_Icomp_I200 = phi_Icomp_I200,
            P_Imax = P_Imax,
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500])
  
  if(return == "calc") return(calc)
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}

# Equation 5- hyperbolic tangent based model
eq5 <- function(pars = 
                  c(Pgmax = 15.5,Isat = 359.2,Rd = .9),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1]
)
{
  if(!missing(dat))
  {
    PARi <- dat$PARi
    if(any(PARi==0)) if(any(!is.na(dat$A[PARi==0]))) pars[["Rd"]] <- abs(dat$A[which(!is.na(dat$A[PARi==0]))][1]) 
    pars[["Pgmax"]] <- max(dat$A,na.rm = TRUE) + pars[["Rd"]]
    pars <- optim(par = pars,function(pars) sum(((dat$A) - (pars[["Pgmax"]]*tanh(dat$PARi/pars[["Isat"]])-pars[["Rd"]]))^2),control=list(fnscale=1))$par
    try(pars <- optim(par = pars,function(pars) sum(((dat$A) - (pars[["Pgmax"]]*tanh(dat$PARi/pars[["Isat"]])-pars[["Rd"]]))^2),control=list(fnscale=1),method="BFGS")$par,silent=TRUE)
    try({
      mod <- nls(control = nls.control(maxiter = 1000),A~Pgmax*tanh(PARi/Isat)-Rd,data = dat,start = pars)
      pars <- coef(mod)
    },silent=TRUE)
  }
  
  # parameters to optimize
  Pgmax <- pars[["Pgmax"]]
  Isat <- pars[["Isat"]]
  Rd <- pars[["Rd"]]
  
  pred <- Pgmax*tanh(PARi/Isat)-Rd
  Icomp <- Isat*atanh(Rd/Pgmax)
  
  # predicted values
  if(return == "predict") return(pred)
  
  # calculated quantities
  x <- c(.05,.1,.15,.25,.5,.75,.85,.9,.95)
  PARi_fine <- seq(0,6000,by = .1)
  pred_fine <- Pgmax*tanh(PARi_fine/Isat)-Rd
  Isat_x <- atanh((((Pgmax-Rd)*x)+Rd)/Pgmax)*Isat
  phi_I0 <- (Pgmax/Isat)*(1/cosh(0/Isat))^2
  phi_Icomp <- (Pgmax/Isat)*(1/cosh(Icomp/Isat))^2
  lo_inds <- which(PARi_fine <= (max(PARi_fine)-50))
  hi_inds <- lo_inds + (length(PARi_fine) - max(lo_inds))
  
  diffs <- (pred_fine[hi_inds] - pred_fine[lo_inds])
  if(any(diffs <= .067)) Imax <- PARi_fine[which(diffs <= .067)[1]] else Imax <- max(PARi_fine)
  
  Pmax <- max(pred_fine)
  P_Imax <- pred_fine[which.min(abs(Imax-PARi_fine))[1]]
  I0_Icomp <- if(Icomp <= 0) I0_Icomp <- 1:2 else which(PARi_fine <= Icomp)   
  Icomp_I200 <- which(PARi_fine >= Icomp & PARi_fine <= 200)
  
  phi_I0_Icomp <- as.numeric(coef(lm(pred_fine[I0_Icomp] ~ PARi_fine[I0_Icomp]))[2])
  phi_Icomp_I200 <- as.numeric(coef(lm(pred_fine[Icomp_I200] ~ PARi_fine[Icomp_I200]))[2])
  
  P2500 <- pred_fine[PARi_fine == 2500]
  
  Ix <- sapply(x,function(x){
    ind_x <- which((pred_fine - (x*P2500))>0 & PARi_fine < PARi_fine[which.max(pred_fine[PARi_fine <= 2500])])[1]
    if(is.na(ind_x)) return(NA)
    ind_x <- c((ind_x-1):(ind_x+1))
    Ix_mod <- lm(y~x,data.frame(x=PARi_fine[ind_x],y=pred_fine[ind_x]))
    (((x*P2500)-coef(Ix_mod)[[1]])/coef(Ix_mod)[[2]])
  })
  
  calc <- c(Pgmax = Pgmax,
            phi_I0 = phi_I0,
            Rd = Rd,
            Icomp = Icomp,
            Isat_25 = Isat_x[x==.25],Isat_50 = Isat_x[x==.5],
            Isat_75 = Isat_x[x==.75],Isat_85 = Isat_x[x==.85],
            Isat_90 = Isat_x[x==.9],
            Isat_95 = Isat_x[x==.95],
            Psat_25 = .25*(Pgmax-Rd),Psat_50 = .5*(Pgmax-Rd),
            Psat_75 = .75*(Pgmax-Rd),Psat_85 = .85*(Pgmax-Rd),
            Psat_90 = .90*(Pgmax-Rd),
            Psat_95 = .95*(Pgmax-Rd),
            I5 = Ix[x==.05],I10 = Ix[x==.1],I15 = Ix[x==.15],
            I25 = Ix[x==.25],I50 = Ix[x==.5],
            I75 = Ix[x==.75],I85 = Ix[x==.85],
            I90 = Ix[x==.9],
            I95 = Ix[x==.95],
            P25 = .25*Pgmax-Rd,P50 = .5*Pgmax-Rd,
            P75 = .75*Pgmax-Rd,P85 = .85*Pgmax-Rd,
            P90 = .90*Pgmax-Rd,
            P95 = .95*Pgmax-Rd,
            Imax = Imax,
            Pmax = Pmax,
            phi_Icomp = phi_Icomp,
            phi_I0_Icomp = phi_I0_Icomp,
            phi_Icomp_I200 = phi_Icomp_I200,
            P_Imax = P_Imax,
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500])
  
  if(return == "calc") return(calc)
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}

#Equation 6- non-rectangular hyperbola based model
eq6 <- function(pars = 
                  c(Pgmax = 15.5,phi_I0 = .0493,theta = .433,Rd = .9),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1]
)
{
  if(!missing(dat))
  {
    # if(any(dat$PARi==0)) if(!is.na(dat$A[which(dat$PARi==0)[1]])) pars[["Rd"]] <- log(-dat$A[which(dat$PARi==0)[1]])
    # pars[["Pgmax"]] <- max(dat$A,na.rm=TRUE)
    PARi <- dat$PARi
    if(any(PARi==0)) if(any(!is.na(dat$A[PARi==0]))) pars[["Rd"]] <- abs(dat$A[which(!is.na(dat$A[PARi==0]))][1]) 
    pars[["Pgmax"]] <- max(dat$A,na.rm = TRUE) + pars[["Rd"]]
    
    pars[["Rd"]] <- log(pars[["Rd"]])
    
    pars <- optim(par = pars,function(pars) sum(((dat$A) - ((((dat$PARi*pars[["phi_I0"]]+pars[["Pgmax"]])-((((pars[["phi_I0"]]*dat$PARi+pars[["Pgmax"]])^2)-(4*pars[["phi_I0"]]*pars[["Pgmax"]]*pars[["theta"]]*dat$PARi))^0.5))/(2*pars[["theta"]]))-exp(pars[["Rd"]])))^2),control=list(fnscale=1))$par
    try({ pars <- optim(par = pars,function(pars) sum(((dat$A) - ((((dat$PARi*pars[["phi_I0"]]+pars[["Pgmax"]])-((((pars[["phi_I0"]]*dat$PARi+pars[["Pgmax"]])^2)-(4*pars[["phi_I0"]]*pars[["Pgmax"]]*pars[["theta"]]*dat$PARi))^0.5))/(2*pars[["theta"]]))-exp(pars[["Rd"]])))^2),control=list(fnscale=1),method="BFGS")$par},silent = TRUE)
    try({
      mod <- nls(control = nls.control(maxiter = 1000),A~(((PARi*phi_I0+Pgmax)-((((phi_I0*PARi+Pgmax)^2)-(4*phi_I0*Pgmax*theta*PARi))^0.5))/(2*theta))-exp(Rd),data = dat,start = pars)
      pars <- coef(mod)
    },silent = TRUE)
    pars[["Rd"]] <- exp(pars[["Rd"]])
  }
  
  # parameters to optimize
  Pgmax <- pars[["Pgmax"]]
  phi_I0 <- pars[["phi_I0"]]
  theta <- pars[["theta"]]
  Rd <- pars[["Rd"]]
  
  pred <- (((PARi*phi_I0+Pgmax)-((((phi_I0*PARi+Pgmax)^2)-(4*phi_I0*Pgmax*theta*PARi))^0.5))/(2*theta))-Rd
  Icomp <- (Rd*(theta*Rd-Pgmax))/(phi_I0*(Rd-Pgmax))
  
  # predicted values
  if(return == "predict") return(pred)
  
  # calculated quantities
  x <- c(.05,.1,.15,.25,.5,.75,.85,.9,.95)
  PARi_fine <- seq(0,6000,by = .1)
  pred_fine <- (((PARi_fine*phi_I0+Pgmax)-((((phi_I0*PARi_fine+Pgmax)^2)-(4*phi_I0*Pgmax*theta*PARi_fine))^0.5))/(2*theta))-Rd
  Isat_x <- (Pgmax*(x*Pgmax+Rd*(1-x))-theta*(x*(Pgmax-Rd)+Rd)^2)/(phi_I0*(Pgmax*(1-x)+Rd*(x-1)))
  phi_Icomp <- (phi_I0/(2*theta))*(1-((phi_I0*Icomp+Pgmax-2*theta*Pgmax)/((((phi_I0*Icomp+Pgmax)^2)-4*theta*phi_I0*Icomp*Pgmax)^0.5)))
  lo_inds <- which(PARi_fine <= (max(PARi_fine)-50))
  hi_inds <- lo_inds + (length(PARi_fine) - max(lo_inds))
  
  diffs <- (pred_fine[hi_inds] - pred_fine[lo_inds])
  if(any(diffs <= .067)) Imax <- PARi_fine[which(diffs <= .067)[1]] else Imax <- max(PARi_fine)
  
  Pmax <- max(pred_fine)
  P_Imax <- pred_fine[which.min(abs(Imax -PARi_fine))[1]]
  I0_Icomp <- if(Icomp <= 0) I0_Icomp <- 1:2 else which(PARi_fine <= Icomp)   
  Icomp_I200 <- which(PARi_fine >= Icomp & PARi_fine <= 200)
  
  phi_I0_Icomp <- as.numeric(coef(lm(pred_fine[I0_Icomp] ~ PARi_fine[I0_Icomp]))[2])
  phi_Icomp_I200 <- as.numeric(coef(lm(pred_fine[Icomp_I200] ~ PARi_fine[Icomp_I200]))[2])
  
  P2500 <- pred_fine[PARi_fine == 2500]
  
  Ix <- sapply(x,function(x){
    ind_x <- which((pred_fine - (x*P2500))>0 & PARi_fine < PARi_fine[which.max(pred_fine[PARi_fine <= 2500])])[1]
    if(is.na(ind_x)) return(NA)
    ind_x <- c((ind_x-1):(ind_x+1))
    Ix_mod <- lm(y~x,data.frame(x=PARi_fine[ind_x],y=pred_fine[ind_x]))
    (((x*P2500)-coef(Ix_mod)[[1]])/coef(Ix_mod)[[2]])
  })
  
  calc <- c(Pgmax = Pgmax,
            phi_I0 = phi_I0,
            Rd = Rd,
            Icomp = Icomp,
            Isat_25 = Isat_x[x==.25],Isat_50 = Isat_x[x==.5],
            Isat_75 = Isat_x[x==.75],Isat_85 = Isat_x[x==.85],
            Isat_90 = Isat_x[x==.9],
            Isat_95 = Isat_x[x==.95],
            Psat_25 = .25*(Pgmax-Rd),Psat_50 = .5*(Pgmax-Rd),
            Psat_75 = .75*(Pgmax-Rd),Psat_85 = .85*(Pgmax-Rd),
            Psat_90 = .90*(Pgmax-Rd),
            Psat_95 = .95*(Pgmax-Rd),
            I5 = Ix[x==.05],I10 = Ix[x==.1],I15 = Ix[x==.15],
            I25 = Ix[x==.25],I50 = Ix[x==.5],
            I75 = Ix[x==.75],I85 = Ix[x==.85],
            I90 = Ix[x==.9],
            I95 = Ix[x==.95],
            P25 = .25*Pgmax-Rd,P50 = .5*Pgmax-Rd,
            P75 = .75*Pgmax-Rd,P85 = .85*Pgmax-Rd,
            P90 = .90*Pgmax-Rd,
            P95 = .95*Pgmax-Rd,
            Imax = Imax,
            Pmax = Pmax,
            phi_Icomp = phi_Icomp,
            phi_I0_Icomp = phi_I0_Icomp,
            phi_Icomp_I200 = phi_Icomp_I200,
            P_Imax = P_Imax,
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500],
            theta = theta)
  
  if(return == "calc") return(calc)
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}

# Equation 8- Exponential Based Model
eq8 <- function(pars = 
                  c(Pgmax = 16.2,phi_I0 = .0597,Rd = 1.3),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1]
)
{
  if(!missing(dat))
  {
    PARi <- dat$PARi
    if(any(PARi==0)) if(any(!is.na(dat$A[PARi==0]))) pars[["Rd"]] <- abs(dat$A[which(!is.na(dat$A[PARi==0]))][1]) 
    pars[["Pgmax"]] <- max(dat$A,na.rm = TRUE) + pars[["Rd"]]
    
    pars <- optim(par = pars,function(pars) sum(((dat$A) - ((pars[["Pgmax"]]*(1-exp(-pars[["phi_I0"]]*dat$PARi/pars[["Pgmax"]])))-pars[["Rd"]]))^2),control=list(fnscale=1))$par
    
    try({
      mod <- nls(control = nls.control(maxiter = 1000),A~(Pgmax*(1-exp(-phi_I0*PARi/Pgmax)))-Rd,data = dat,start = pars)
      pars <- coef(mod)
    },silent = TRUE)
  }
  
  # parameters to optimize
  Pgmax <- pars[["Pgmax"]]
  phi_I0 <- pars[["phi_I0"]]
  Rd <- pars[["Rd"]]
  
  pred <- (Pgmax*(1-exp(-phi_I0*PARi/Pgmax)))-Rd
  
  # predicted values
  if(return == "predict") return(pred)
  
  Icomp <- (Pgmax/phi_I0)*(-log(1-Rd/Pgmax))
  
  # calculated quantities
  x <- c(.05,.1,.15,.25,.5,.75,.85,.9,.95)
  PARi_fine <- seq(0,6000,by = .1)
  pred_fine <- (Pgmax*(1-exp(-phi_I0*PARi_fine/Pgmax)))-Rd
  Isat_x <- (Pgmax/phi_I0)*(-log(1-(x*(Pgmax-Rd)+Rd)/Pgmax))
  phi_Icomp <- phi_I0*exp(-phi_I0*Icomp/Pgmax)
  lo_inds <- which(PARi_fine <= (max(PARi_fine)-50))
  hi_inds <- lo_inds + (length(PARi_fine) - max(lo_inds))
  
  diffs <- (pred_fine[hi_inds] - pred_fine[lo_inds])
  if(any(diffs <= .067)) Imax <- PARi_fine[which(diffs <= .067)[1]] else Imax <- max(PARi_fine)
  
  Pmax <- max(pred_fine)
  P_Imax <- pred_fine[which.min(abs(Imax -PARi_fine))[1]]
  I0_Icomp <- if(Icomp <= 0) I0_Icomp <- 1:2 else which(PARi_fine <= Icomp)   
  Icomp_I200 <- which(PARi_fine >= Icomp & PARi_fine <= 200)
  
  phi_I0_Icomp <- as.numeric(coef(lm(pred_fine[I0_Icomp] ~ PARi_fine[I0_Icomp]))[2])
  phi_Icomp_I200 <- as.numeric(coef(lm(pred_fine[Icomp_I200] ~ PARi_fine[Icomp_I200]))[2])
  
  P2500 <- pred_fine[PARi_fine == 2500]
  
  Ix <- sapply(x,function(x){
    ind_x <- which((pred_fine - (x*P2500))>0 & PARi_fine < PARi_fine[which.max(pred_fine[PARi_fine <= 2500])])[1]
    if(is.na(ind_x)) return(NA)
    ind_x <- c((ind_x-1):(ind_x+1))
    Ix_mod <- lm(y~x,data.frame(x=PARi_fine[ind_x],y=pred_fine[ind_x]))
    (((x*P2500)-coef(Ix_mod)[[1]])/coef(Ix_mod)[[2]])
  })
  
  calc <- c(Pgmax = Pgmax,
            phi_I0 = phi_I0,
            Rd = Rd,
            Icomp = Icomp,
            Isat_25 = Isat_x[x==.25],Isat_50 = Isat_x[x==.5],
            Isat_75 = Isat_x[x==.75],Isat_85 = Isat_x[x==.85],
            Isat_90 = Isat_x[x==.9],
            Isat_95 = Isat_x[x==.95],
            Psat_25 = .25*(Pgmax-Rd),Psat_50 = .5*(Pgmax-Rd),
            Psat_75 = .75*(Pgmax-Rd),Psat_85 = .85*(Pgmax-Rd),
            Psat_90 = .90*(Pgmax-Rd),
            Psat_95 = .95*(Pgmax-Rd),
            I5 = Ix[x==.05],I10 = Ix[x==.1],I15 = Ix[x==.15],
            I25 = Ix[x==.25],I50 = Ix[x==.5],
            I75 = Ix[x==.75],I85 = Ix[x==.85],
            I90 = Ix[x==.9],
            I95 = Ix[x==.95],
            P25 = .25*Pgmax-Rd,P50 = .5*Pgmax-Rd,
            P75 = .75*Pgmax-Rd,P85 = .85*Pgmax-Rd,
            P90 = .90*Pgmax-Rd,
            P95 = .95*Pgmax-Rd,
            Imax = Imax,
            Pmax = Pmax,
            phi_Icomp = phi_Icomp,
            phi_I0_Icomp = phi_I0_Icomp,
            phi_Icomp_I200 = phi_Icomp_I200,
            P_Imax = P_Imax,
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500])
  
  if(return == "calc") return(calc)
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}

#Equation 9- Exponential based model
eq9 <- function(pars = 
                  c(Pgmax = 22,Icomp = 10,k = .0015,Rd = .1),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1]
)
{
  if(!missing(dat))
  {
    PARi <- dat$PARi
    if(any(PARi==0)) if(any(!is.na(dat$A[PARi==0]))) pars[["Rd"]] <- abs(dat$A[which(!is.na(dat$A[PARi==0]))][1]) 
    pars[["Pgmax"]] <- max(dat$A,na.rm = TRUE) + pars[["Rd"]]
    pars <- optim(par = pars,function(pars) sum(((dat$A) - (pars[["Pgmax"]]*(1-exp(-pars[["k"]]*(dat$PARi-pars[["Icomp"]])))-(pars[["Rd"]])))^2),control=list(fnscale=1))$par
    try({ pars <- optim(par = pars,function(pars) sum(((dat$A) - (pars[["Pgmax"]]*(1-exp(-pars[["k"]]*(dat$PARi-pars[["Icomp"]])))-(pars[["Rd"]])))^2),control=list(fnscale=1),method = "BFGS")$par },silent = TRUE)
    pars <- optim(par = pars,function(pars) sum(((dat$A) - (pars[["Pgmax"]]*(1-exp(-pars[["k"]]*(dat$PARi-pars[["Icomp"]])))-(pars[["Rd"]])))^2),control=list(fnscale=1))$par
    pars <- optim(par = pars,function(pars) sum(((dat$A) - (pars[["Pgmax"]]*(1-exp(-pars[["k"]]*(dat$PARi-pars[["Icomp"]])))-(pars[["Rd"]])))^2),control=list(fnscale=1))$par
    
    try({
      mod <- nls(control = nls.control(maxiter = 1000),formula = A~Pgmax*(1-exp(-k*(PARi-Icomp)))-Rd,
                 data = dat,
                 start = pars,
                 lower = c(Pgmax=5,Icomp=5,k=.0001,Rd=0),
                 upper = c(Pgmax=60,Icomp=150,k=.02,Rd=5),
                 algorithm = "port")
      pars <- coef(mod)
      
    },silent = TRUE)
  }
  
  # parameters to optimize
  Pgmax <- pars[["Pgmax"]]
  Icomp <- pars[["Icomp"]]
  Rd <- pars[["Rd"]]
  k <- pars[["k"]]
  
  pred <- Pgmax*(1-exp(-k*(PARi-Icomp)))-Rd
  
  # predicted values
  if(return == "predict") return(pred)
  
  # calculated quantities
  x <- c(.05,.1,.15,.25,.5,.75,.85,.9,.95)
  PARi_fine <- seq(0,6000,by = .1)
  pred_fine <- Pgmax*(1-exp(-k*(PARi_fine-Icomp)))-Rd
  Isat_x <- Icomp-(log(1-((x*(Pgmax-Rd)+Rd)/Pgmax)))/k
  phi_I0 <- Pgmax*k*exp(-k*(Icomp))
  phi_Icomp <- Pgmax*k*exp(-k*(Icomp-Icomp))
  lo_inds <- which(PARi_fine <= (max(PARi_fine)-50))
  hi_inds <- lo_inds + (length(PARi_fine) - max(lo_inds))
  
  diffs <- (pred_fine[hi_inds] - pred_fine[lo_inds])
  if(any(diffs <= .067)) Imax <- PARi_fine[which(diffs <= .067)[1]] else Imax <- max(PARi_fine)
  
  Pmax <- max(pred_fine)
  P_Imax <- pred_fine[which.min(abs(Imax -PARi_fine))[1]]
  I0_Icomp <- if(Icomp <= 0) I0_Icomp <- 1:2 else which(PARi_fine <= Icomp)   
  Icomp_I200 <- which(PARi_fine >= Icomp & PARi_fine <= 200)
  
  phi_I0_Icomp <- as.numeric(coef(lm(pred_fine[I0_Icomp] ~ PARi_fine[I0_Icomp]))[2])
  phi_Icomp_I200 <- as.numeric(coef(lm(pred_fine[Icomp_I200] ~ PARi_fine[Icomp_I200]))[2])
  
  P2500 <- pred_fine[PARi_fine == 2500]
  
  Ix <- sapply(x,function(x){
    ind_x <- which((pred_fine - (x*P2500))>0 & PARi_fine < PARi_fine[which.max(pred_fine[PARi_fine <= 2500])])[1]
    if(is.na(ind_x)) return(NA)
    ind_x <- c((ind_x-1):(ind_x+1))
    Ix_mod <- lm(y~x,data.frame(x=PARi_fine[ind_x],y=pred_fine[ind_x]))
    (((x*P2500)-coef(Ix_mod)[[1]])/coef(Ix_mod)[[2]])
  })
  
  calc <- c(Pgmax = Pgmax,
            phi_I0 = phi_I0,
            Rd = Rd,
            Icomp = Icomp,
            Isat_25 = Isat_x[x==.25],Isat_50 = Isat_x[x==.5],
            Isat_75 = Isat_x[x==.75],Isat_85 = Isat_x[x==.85],
            Isat_90 = Isat_x[x==.9],
            Isat_95 = Isat_x[x==.95],
            Psat_25 = .25*(Pgmax-Rd),Psat_50 = .5*(Pgmax-Rd),
            Psat_75 = .75*(Pgmax-Rd),Psat_85 = .85*(Pgmax-Rd),
            Psat_90 = .90*(Pgmax-Rd),
            Psat_95 = .95*(Pgmax-Rd),
            I5 = Ix[x==.05],I10 = Ix[x==.1],I15 = Ix[x==.15],
            I25 = Ix[x==.25],I50 = Ix[x==.5],
            I75 = Ix[x==.75],I85 = Ix[x==.85],
            I90 = Ix[x==.9],
            I95 = Ix[x==.95],
            P25 = .25*Pgmax-Rd,P50 = .5*Pgmax-Rd,
            P75 = .75*Pgmax-Rd,P85 = .85*Pgmax-Rd,
            P90 = .90*Pgmax-Rd,
            P95 = .95*Pgmax-Rd,
            Imax = Imax,
            Pmax = Pmax,
            phi_Icomp = phi_Icomp,
            phi_I0_Icomp = phi_I0_Icomp,
            phi_Icomp_I200 = phi_Icomp_I200,
            P_Imax = P_Imax,
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500],
            k = k)
  
  if(return == "calc") return(calc)
  
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}

# Equation 11- Ye model (2007), a modified rectangular hyperbola model able to fit the photoinhibition stage and estimate Pgmax and Isat.
eq11 <- function(pars = 
                   c(phi_I0_Icomp = .0756,beta = .0000432,gamma = .0039,Icomp = 22.6),
                 dat,
                 PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                 return = c("predict","calc","all")[1]
)
{
  if(!missing(dat))
  {
    try({
      mod <- nls(control = nls.control(maxiter = 1000),A~phi_I0_Icomp*((1-beta*PARi)/(1+gamma*PARi))*(PARi-Icomp),data = dat,start = pars,lower = c(phi_I0_Icomp = .01,beta = 1e-5,gamma = 1e-5,Icomp = 5),upper = c(phi_I0_Icomp = .2,beta = 1,gamma = 1,Icomp = 150),algorithm = "port")
      pars <- coef(mod)
    },silent = TRUE)
    
    pars[["beta"]] <- log(pars[["beta"]])
    pars[["gamma"]] <- log(pars[["gamma"]])
    
    try({
      mod <- nls(control = nls.control(maxiter = 1000),A~phi_I0_Icomp*((1-exp(beta)*PARi)/(1+exp(gamma)*PARi))*(PARi-Icomp),data = dat,start = pars)
      pars <- coef(mod)
    },silent = TRUE)
    
    PARi <- dat$PARi
    if(!any(dat$PARi==0))
    {
      I0 <- coef(lm(A~PARi,data = dat[order(dat$PARi)[1:2],]))[1]
      if(I0<0) pars[["Icomp"]] <- approx(x = c(I0,dat$A),y = c(0,dat$PARi),xout = 0)$y
    } else
    {
      I0 <- coef(lm(A~PARi,data = dat[order(dat$PARi)[1:2],]))[1]
      pars[["Icomp"]] <- approx(x = dat$A,y = dat$PARi,xout = 0)$y
    }
    
    pars[["phi_I0_Icomp"]] <- abs(I0) / pars[["Icomp"]]
    pars <- optim(par = pars,function(pars) sum(((dat$A) - (pars[["phi_I0_Icomp"]]*((1-exp(pars[["beta"]])*dat$PARi)/(1+exp(pars[["gamma"]])*dat$PARi))*(dat$PARi-pars[["Icomp"]])))^2),control=list(fnscale=1))$par
    try({ pars <- optim(par = pars,function(pars) sum(((dat$A) - (pars[["phi_I0_Icomp"]]*((1-exp(pars[["beta"]])*dat$PARi)/(1+exp(pars[["gamma"]])*dat$PARi))*(dat$PARi-pars[["Icomp"]])))^2),control=list(fnscale=1),method = "BFGS")$par },silent = TRUE)
    
    try({
      mod <- nls(control = nls.control(maxiter = 1000),A~phi_I0_Icomp*((1-exp(beta)*PARi)/(1+exp(gamma)*PARi))*(PARi-Icomp),data = dat,start = pars)
      pars <- coef(mod)
    },silent = TRUE)
    
    pars[["beta"]] <- exp(pars[["beta"]])
    pars[["gamma"]] <- exp(pars[["gamma"]])
    
    try({
      mod <- nls(control = nls.control(maxiter = 1000),A~phi_I0_Icomp*((1-beta*PARi)/(1+gamma*PARi))*(PARi-Icomp),data = dat,start = pars,lower = c(phi_I0_Icomp = .01,beta = 1e-32,gamma = 1e-32,Icomp = 5),upper = c(phi_I0_Icomp = .2,beta = 1,gamma = 1,Icomp = 150),algorithm = "port")
      pars <- coef(mod)
    },silent = TRUE)
  }
  
  # parameters to optimize
  phi_I0_Icomp <- pars[["phi_I0_Icomp"]]
  Icomp <- pars[["Icomp"]]
  beta <- pars[["beta"]]
  gamma <- pars[["gamma"]]
  
  pred <- phi_I0_Icomp*((1-beta*PARi)/(1+gamma*PARi))*(PARi-Icomp)
  
  # predicted values
  if(return == "predict") return(pred)
  
  # calculated quantities
  x <- c(.05,.1,.15,.25,.5,.75,.85,.9,.95)
  PARi_fine <- seq(0,6000,by = .1)
  pred_fine <- phi_I0_Icomp*((1-beta*PARi_fine)/(1+gamma*PARi_fine))*(PARi_fine-Icomp)
  Isat <- ((((beta+gamma)*(1+gamma*Icomp)/beta)^0.5)-1)/gamma
  Rd <- phi_I0_Icomp*Icomp
  Pgmax <- phi_I0_Icomp*((1-beta*Isat)/(1+gamma*Isat))*(Isat-Icomp)+Rd
  Isat_x <- (((phi_I0_Icomp*beta*Icomp)-(x*(Pgmax-Rd)*gamma)+phi_I0_Icomp)-(((x*(Pgmax-Rd)*gamma)-(phi_I0_Icomp*beta*Icomp)-phi_I0_Icomp)^2-(4*phi_I0_Icomp*beta*(phi_I0_Icomp*Icomp+x*(Pgmax-Rd))))^0.5)/(2*phi_I0_Icomp*beta)
  phi_Icomp <- phi_I0_Icomp*((1-2*beta*Icomp-beta*gamma*(Icomp^2)+(gamma+beta)*Icomp)/(1+gamma*Icomp)^2)
  phi_I0 <- phi_I0_Icomp*((1-2*beta*0-beta*gamma*(0^2)+(gamma+beta)*0)/(1+gamma*0)^2)
  lo_inds <- which(PARi_fine <= (max(PARi_fine)-50))
  hi_inds <- lo_inds + (length(PARi_fine) - max(lo_inds))
  
  diffs <- (pred_fine[hi_inds] - pred_fine[lo_inds])
  if(any(diffs <= .067)) Imax <- PARi_fine[which(diffs <= .067)[1]] else Imax <- max(PARi_fine)
  
  Pmax <- max(pred_fine)
  P_Imax <- pred_fine[which.min(abs(Imax -PARi_fine))[1]]
  I0_Icomp <- if(Icomp <= 0) I0_Icomp <- 1:2 else which(PARi_fine <= Icomp)   
  Icomp_I200 <- which(PARi_fine >= Icomp & PARi_fine <= 200)
  
  phi_I0_Icomp <- as.numeric(coef(lm(pred_fine[I0_Icomp] ~ PARi_fine[I0_Icomp]))[2])
  phi_Icomp_I200 <- as.numeric(coef(lm(pred_fine[Icomp_I200] ~ PARi_fine[Icomp_I200]))[2])
  
  P2500 <- pred_fine[PARi_fine == 2500]
  
  Ix <- sapply(x,function(x){
    ind_x <- which((pred_fine - (x*P2500))>0 & PARi_fine < PARi_fine[which.max(pred_fine[PARi_fine <= 2500])])[1]
    if(is.na(ind_x)) return(NA)
    ind_x <- c((ind_x-1):(ind_x+1))
    Ix_mod <- lm(y~x,data.frame(x=PARi_fine[ind_x],y=pred_fine[ind_x]))
    (((x*P2500)-coef(Ix_mod)[[1]])/coef(Ix_mod)[[2]])
  })
  
  calc <- c(Pgmax = Pgmax,
            phi_I0 = phi_I0,
            Rd = Rd,
            Icomp = Icomp,
            Isat_25 = Isat_x[x==.25],Isat_50 = Isat_x[x==.5],
            Isat_75 = Isat_x[x==.75],Isat_85 = Isat_x[x==.85],
            Isat_90 = Isat_x[x==.9],
            Isat_95 = Isat_x[x==.95],
            Psat_25 = .25*(Pgmax-Rd),Psat_50 = .5*(Pgmax-Rd),
            Psat_75 = .75*(Pgmax-Rd),Psat_85 = .85*(Pgmax-Rd),
            Psat_90 = .90*(Pgmax-Rd),
            Psat_95 = .95*(Pgmax-Rd),
            I5 = Ix[x==.05],I10 = Ix[x==.1],I15 = Ix[x==.15],
            I25 = Ix[x==.25],I50 = Ix[x==.5],
            I75 = Ix[x==.75],I85 = Ix[x==.85],
            I90 = Ix[x==.9],
            I95 = Ix[x==.95],
            P25 = .25*Pgmax-Rd,P50 = .5*Pgmax-Rd,
            P75 = .75*Pgmax-Rd,P85 = .85*Pgmax-Rd,
            P90 = .90*Pgmax-Rd,
            P95 = .95*Pgmax-Rd,
            Imax = Imax,
            Pmax = Pmax,
            phi_Icomp = phi_Icomp,
            phi_I0_Icomp = phi_I0_Icomp,
            phi_Icomp_I200 = phi_Icomp_I200,
            P_Imax = P_Imax,
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500],
            beta = beta,
            gamma = gamma)
  
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  
  if(return == "calc") return(calc)
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}
```

## Results for all equations 
```{r}
eqX <- paste("eq",c(1,2,3,4,5,6,8,9,11),sep = "")
par_names <- c(names(eq9(return = "calc")),"beta","gamma")
# par_names <- par_names[-grep("Isat",par_names)]
res <- array(dim = c(9,length(par_names)+length(PARi),length(inds)))
dimnames(res) <- list(paste("eq",c(1:6,8:9,11),sep = ""),c(par_names,paste("PAR_",PARi,sep="")),inds)
mse <- r2 <- matrix(NA,nrow = length(inds),ncol = 9,dimnames = list(inds,eqX))

for(e in 1:9)
{
  for(i in 1:length(inds))
  {
    try({
      out <- get(eqX[e])(dat = dat[dat$SpeciesRepDate==inds[i],],return = "all")
      res[e,,i] <- c(out$calc[par_names],out$pred_fine[which(PARi_fine %in% PARi)])
      r2[i,e] <- out$fit[["r2"]]
      mse[i,e] <- out$fit[["mse"]]
    },silent = TRUE)
  }
}
```

## Ancestral State Reconstruction with Calcualted Values
```{r}
# Storing data to graph
i <- 1
f <- function(eqX,i)
{
  fit <- eqX(dat=dat[dat$SpeciesRepDate==inds[i],],return = "all")
  calc <- fit$calc
  x <- c(
    rep(0,10),
    calc[["Icomp"]]+seq(-1e-2,1e-2,by = 1e-2),
    rep(calc[["I25"]],3),
    rep(calc[["I75"]],3),
    rep(2500,3)
  )
  
  y <- c(
    rep(-calc[["Rd"]],10),
    calc[["phi_Icomp"]]*seq(-1e-2,1e-2,by = 1e-2),
    rep(.25*calc[["Pmax_obs"]],3),
    rep(.75*calc[["Pmax_obs"]],3),
    rep(calc[["Pmax_obs"]],3)
  )
  
  x <- c(
    rep(calc[["I15"]],2),
    rep(calc[["I25"]],2),
    rep(calc[["I85"]],2),
    rep(calc[["I95"]],3),
    rep(2500,4)
  )
  
  y <- c(
    rep(.15*calc[["Pmax_obs"]],2),
    rep(.25*calc[["Pmax_obs"]],2),
    rep(.85*calc[["Pmax_obs"]],2),
    rep(.95*calc[["Pmax_obs"]],3),
    rep(calc[["Pmax_obs"]],4)
  )
  
  plot(dat[dat$SpeciesRepDate==inds[i],2:3],main = paste(i,inds[i]))
  points(x,y,pch=19)
  if(any(is.na(x) | is.na(y)))
  {
    exclude <- which(is.na(x) | is.na(y))
    x <- x[-exclude]
    y <- y[-exclude]
  }
  
  recon <- eqX(dat=data.frame(PARi=x,A=y),return = "all")
  
  points(PARi_fine,recon$pred_fine,type='l',col='red',lwd=2)
  points(PARi_fine,fit$pred_fine,type='l',col='blue',lty=2)
  
  recon
}
# dev.off()
# i=sample(length(inds),1)

#
if(FALSE) for(i in 1:length(inds))
{
  print(i)
  jpeg(filename = paste("results/",i,".jpeg",sep=""),width = 800,height = 1000)
  par(mfrow=c(3,3))
  f(eq1,i)
  f(eq2,i)
  f(eq3,i)
  f(eq4,i)
  f(eq5,i)
  f(eq6,i)
  f(eq8,i)
  f(eq9,i)
  f(eq11,i)
  dev.off()
}

# Evolutionary Models

p <- p_star <- p_star2 <- setNames(vector(mode = "list",length = length(unique(ind_species))),unique(ind_species))
p <- p_star <- p_star2 <- setNames(rep(list(p),9)[[1]],eqX)
# save(p_i,file = "results/multivariate_BM.RData")
# save(p_i_star,file = "results/multivariate_star.RData")
# save(p_i_star2,file = "results/multivariate_star2.RData")
load(file = "results/multivariate_BM.RData")
load(file = "results/multivariate_star.RData")
load(file = "results/multivariate_star2.RData")


res_i_BM <- res_BM <- res_i_star <- res_star <- res_i_star2 <- res_star2 <- array(dim = c(9,length(par_names)+length(PARi),length(inds)))
dimnames(res_BM) <- dimnames(res_star) <- dimnames(res_star2) <- list(paste("eq",c(1:6,8:9,11),sep = ""),c(par_names,paste("PAR_",PARi,sep = "")),inds)
dimnames(res_i_BM) <- dimnames(res_i_star) <- dimnames(res_i_star2) <- dimnames(res_BM)


for(e in 1:9)
{
  for(i in 1:length(unique(ind_species)))
  {
    sp <- ind_species[which(ind_species==unique(ind_species)[i])][1]
    cat(e,i,unique(ind_species)[i],"\n")
    dat_long_i <- dat_long
    cols <- 2:9 #c(2:8,10)
    dat_long_i[which(dat_long$species==unique(ind_species)[i]),colnames(dat_long)[cols]] <- NA
    # p_i[[sp]] <- phylopars(dat_long_i,tree,model = 'BM')
    # p_i_star[[sp]] <- phylopars(data.frame(species = inds,dat_long_i[,-1]),tree = phytools::starTree(inds,rep(1,length(inds))),model = 'star')
    # p_i_star2[[sp]] <- phylopars(dat_long_i,phytools::starTree(unique(ind_species),rep(1,length(unique(ind_species)))),model = 'star')
    
    par_dat <- data.frame(species=ind_species,t(res[eqX[e],c("I15","I25","I85","I95","Pmax_obs"),]),row.names = NULL)
    
    par_dat_i <- par_dat
    cols <- 2:(ncol(par_dat)-1)
    par_dat_i[which(par_dat$species==unique(ind_species)[i]),colnames(par_dat)[cols]] <- NA
    
    p[[eqX[e]]][[sp]] <- phylopars(par_dat_i,tree,model='BM')
    p_star[[eqX[e]]][[sp]] <- phylopars(data.frame(species = inds,par_dat_i[,-1]),tree = phytools::starTree(inds,rep(1,length(inds))),model = 'star')
    p_star2[[eqX[e]]][[sp]] <- phylopars(par_dat_i,phytools::starTree(unique(ind_species),rep(1,length(unique(ind_species)))),model = 'star')
    
    if(i>1)
    {
      temp <- phylopars(par_dat_i,tree,model='BM',phylocov_start = p[[eqX[e]]][[unique(ind_species)[i-1]]]$pars[[1]],phenocov_start = p[[eqX[e]]][[unique(ind_species)[i-1]]]$pars[[2]])
      if(logLik(temp)[[1]] > logLik(p[[eqX[e]]][[sp]])[[1]]) p[[eqX[e]]][[sp]] <- temp
    }
    
    # save.image("results/curve_recons.RData")
    # dev.off()
    for(j in 1:length(which(unique(ind_species)[i]==ind_species)))
      
    {
      ind <- inds[which(ind_species==unique(ind_species)[i])][j]
      
      f(eqX = get(eqX[e]),i = which(ind_species==unique(ind_species)[i])[j])
      
      points(PARi,dat_long[which(dat_long$species==unique(ind_species)[i]),2:(ncol(dat_long))][j,],pch = 19,type='b',lwd=2,col = "darkorange")
      
      points(PARi,unlist(p_i[[sp]]$ind_recon[which(dat_long$species==unique(ind_species)[i]),][j,-1]),pch=19,col = 'blue',type='l')
      out <- get(eqX[e])(dat = data.frame(PARi=PARi,A=unlist(p_i[[sp]]$ind_recon[which(dat_long$species==unique(ind_species)[i]),][j,-1])),return = "all")
      res_i_BM[eqX[e],par_names,ind] <- out$calc[par_names]
      res_i_BM[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]
      
      points(PARi,unlist(p_i_star[[sp]]$anc_recon[which(dat_long$species==unique(ind_species)[i]),,drop = FALSE][j,]),pch=19,col = 'red',type='l')
      out <- get(eqX[e])(dat = data.frame(PARi=PARi,A=unlist(p_i_star[[sp]]$anc_recon[which(dat_long$species==unique(ind_species)[i]),,drop = FALSE][j,])),return = "all")
      res_i_star[eqX[e],par_names,ind] <- out$calc[par_names]
      res_i_star[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]
      
      points(PARi,unlist(p_i_star2[[sp]]$ind_recon[which(dat_long$species==unique(ind_species)[i]),1:(ncol(dat_long))][j,-1]),pch=19,col = 'purple',type='l')
      out <- get(eqX[e])(dat = data.frame(PARi=PARi,A=unlist(p_i_star2[[sp]]$ind_recon[which(dat_long$species==unique(ind_species)[i]),1:(ncol(dat_long))][j,-1])),return = "all")
      res_i_star2[eqX[e],par_names,ind] <- out$calc[par_names]
      res_i_star2[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]
      
      x <- c(unlist(p[[eqX[e]]][[sp]]$ind_recon[which(dat_long$species==unique(ind_species)[i]),][j,-c(1,ncol(p[[eqX[e]]][[sp]]$ind_recon))]),2500)
      x <- c(rep(x[1:3],2),rep(x[4],3),rep(x[5],4))
      y <- c(.15,.25,.85,.95,1)*unlist(p[[eqX[e]]][[sp]]$ind_recon[which(dat_long$species==unique(ind_species)[i]),][j,ncol(p[[eqX[e]]][[sp]]$ind_recon)])
      y <- c(rep(y[1:3],2),rep(y[4],3),rep(y[5],4))
      
      out <- get(eqX[e])(dat = data.frame(PARi=x,A=y),return = "all")
      res_BM[eqX[e],par_names,ind] <- out$calc[par_names]
      res_BM[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]
      
      points(PARi_fine,out$pred_fine,lty=3,type='l',lwd=3)
      
      x <- c(unlist(p_star[[eqX[e]]][[sp]]$anc_recon[ind,,drop = FALSE][1,-ncol(p_star[[eqX[e]]][[sp]]$anc_recon)]),2500)
      x <- c(rep(x[1:3],2),rep(x[4],3),rep(x[5],4))
      y <- c(.15,.25,.85,.95,1)*unlist(p_star[[eqX[e]]][[sp]]$anc_recon[ind,,drop = FALSE][1,ncol(p_star[[eqX[e]]][[sp]]$anc_recon)])
      y <- c(rep(y[1:3],2),rep(y[4],3),rep(y[5],4))
      
      out <- get(eqX[e])(dat = data.frame(PARi=x,A=y),return = "all")
      res_star[eqX[e],par_names,ind] <- out$calc[par_names]
      res_star[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]
      
      points(PARi_fine,out$pred_fine,lty=3,type='l',lwd=3,col = "darkgreen")
      
      x <- c(unlist(p_star2[[eqX[e]]][[sp]]$ind_recon[which(p_star2[[eqX[e]]][[sp]]$ind_recon[,1]==sp),,drop = FALSE][j,-c(1,ncol(p_star2[[eqX[e]]][[sp]]$ind_recon))]),2500)
      x <- c(rep(x[1:3],2),rep(x[4],3),rep(x[5],4))
      y <- c(.15,.25,.85,.95,1)*unlist(p_star2[[eqX[e]]][[sp]]$ind_recon[which(p_star2[[eqX[e]]][[sp]]$ind_recon[,1]==sp),,drop = FALSE][j,ncol(p_star2[[eqX[e]]][[sp]]$ind_recon)])
      y <- c(rep(y[1:3],2),rep(y[4],3),rep(y[5],4))
      
      out <- get(eqX[e])(dat = data.frame(PARi=x,A=y),return = "all")
      res_star2[eqX[e],par_names,ind] <- out$calc[par_names]
      res_star2[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]
      
      points(PARi_fine,out$pred_fine,lty=3,type='l',lwd=3,col = "darkblue")
      
    }
  }
}

res_vs_res_BM <- t(sapply(1:9,function(i) sapply(1:ncol(res),function(j) cor(res[i,j,],res_i_BM[i,j,],use = "pair"))))
colnames(res_vs_res_BM) <- colnames(res)
rownames(res_vs_res_BM) <- rownames(res)

round(res_vs_res_BM^2,2)[,c("Rd","Pgmax",paste("PAR_",PARi,sep=""))]

# save.image("results/res.RData")
```

## Evaluating Models
```{r}

pars <- c("Pmax_obs","Pgmax","phi_I0","phi_Icomp","phi_I0_Icomp","phi_Icomp_I200","Rd","Icomp",
          "I15","I25","I85","I95",
          "PAR_0","PAR_50","PAR_100","PAR_250","PAR_500","PAR_1000","PAR_1500","PAR_2000","PAR_2500")

r2_impute <- suppressWarnings(sapply(1:length(eqX),function(e) {
  setNames(round(sapply(1:ncol(res),function(i) cor(t(
    res[eqX[e],,])[,i],
    t(res_BM[eqX[e],,])[,i]
  ))^2,5),colnames(res))}))[pars,]
colnames(r2_impute) <- eqX
r2_impute

mse_impute <- suppressWarnings(sapply(1:length(eqX),function(e) {
  setNames(round(sapply(1:ncol(res),function(i) 
  {
    mean((t(res[eqX[e],,])[,i] - t(res_BM[eqX[e],,])[,i])^2)
  }),5),colnames(res))[pars]
}))
colnames(mse_impute) <- eqX
mse_impute

# R2 between parameters on observed vs imputed using parameters
suppressWarnings(cbind(
  setNames(round(sapply(1:ncol(res),function(i) cor(t(
    res["eq6",,])[,i],
    t(res_BM["eq6",,])[,i]
  ))^2,5),colnames(res))
)) 

# Examine MSE on phylogeny
# (log sqrt transformed so take this with grain of salt)
contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = log(mse[,"eq6"]^.5)))[,1])

# visuaize distribution of parameter of interest from an equation
hist(res["eq6","I95",])


phylo_fit <- phylopars(data.frame(species = ind_species,trait = res["eq6","phi_Icomp_I200",]),tree = tree,model = 'lambda')


range(phylo_fit$anc_recon[,1])
range(res["eq6","phi_I0",])

plot(tree)

nodelabels(frame = "none",col = "red")
```

## Plotting Results
```{r}
res[,,"Agrestis_1_29/10/19"]

plot(dat_long[,c("PARi_100")],sapply(1:length(inds),function(i) eq11(dat = dat[dat$SpeciesRepDate==inds[i],],return = "all")$pred_fine[PARi_fine==100]))
abline(a=0,b=1)


t(res["eq11",c("Rd","phi_I0_Icomp","I85"),])

range(r2)
mse

res_i_BM[1,,]
res_i_star

plot(PARi_fine,eq1(dat = dat[dat$SpeciesRepDate == inds[1],],return = "all")$pred_fine,type='l',xlim = c(0,2500))
```

## Quesions to Address from the Analysis
```{r}
# 1. Which model type best fits photosynthetic light response curves in sunflowers?
# fit 9 different model types; table comparing average mse and r2; representative figure
mse[,"eq6"]
r2[,"eq6"]

# 2. Does photosynthetic rate at different light levels exhibit phylogenetic signal?
# test phylogenetic signal at each level
data.frame(PARi = PARi,phy_signal = 
             sapply(2:ncol(dat_long),function(i) phylopars(dat_long[,c(1,i)],tree,model = 'lambda')$model$lambda))

dat_long=as.data.frame(dat_long)
vals <- sapply(1:length(inds),function(i) eq6(dat = dat[dat$SpeciesRepDate==inds[i],],return = "all")$pred_fine[PARi_fine%in%seq(0,2500,by = 50)])
a=sapply(1:nrow(vals),function(i) {
  phy_fit <- phylopars(data.frame(species=ind_species,trait=vals[i,]),tree,model='lambda')
  # phy_fit$model$lambda
  # phy_fit$pars[[1]]/(phy_fit$pars[[1]]+phy_fit$pars[[2]])
})


plot(seq(0,2500,by = 50),a,type='l',ylim = c(0,1))
points(PARi,sapply(2:ncol(dat_long),function(i) phylopars(dat_long[,c(1,i)],tree,model = 'lambda')$model$lambda),pch=19)


# 3. Do parameters extracted from photosynthetic light response curves exhibit phylogenetic signal?
# rows = equations; columns = parameters
setNames(sapply(1:length(pars),function(i){
  phylopars(data.frame(species = ind_species,trait = t(res["eq6",,])[,pars[i]],row.names = NULL),tree = tree,model = 'lambda')$model$lambda
}),pars)


# 4. Can we use phylogenetic comparative methods to predict light response curves in new species given only survey measurements?
# R2 and mse for LOO-CV

suppressWarnings(cbind(
  # R2 between parameters on observed vs imputed using parameters
  setNames(round(sapply(1:ncol(res),function(i) cor(t(
    res["eq6",,])[,i],
    t(res_BM["eq6",,])[,i]
  ))^2,5),colnames(res)),
  
  # R2 between parameters on observed vs imputed using values
  setNames(round(sapply(1:ncol(res),function(i) cor(t(
    res["eq6",,])[,i],
    t(res_i_BM["eq6",,])[,i]
  ))^2,5),colnames(res))))


suppressWarnings(cbind(
  # MSE between parameters on observed vs imputed using parameters
  setNames(round(sapply(1:ncol(res),function(i) 
    mean((t(res["eq6",,])[,i] - t(res_BM["eq6",,])[,i])^2
    )),3),colnames(res)),
  # MSE between parameters on observed vs imputed using values
  setNames(round(sapply(1:ncol(res),function(i)
    mean((t(res["eq6",,])[,i] - t(res_i_BM["eq6",,])[,i])^2
    )),3),colnames(res))))

```

## Defining Important Objects and Paramters
```{r}
# Pgmax = Pmax + Rd

# Pmax_obs <- calc[["Pmax_obs"]]
# x=c(0,Icomp,200,I50,2500)
# y=c(-Rd,0,phi_Icomp_I200*(200-Icomp),.5*(Pgmax-Rd),Pmax_obs)


# Rd gives: (0,-Rd)
# phi_I0 gives: (.01,-Rd+phi_I0*.01) and (-.01,-Rd-phi_I0*.01)
# Icomp gives: (Icomp,0)
# phi_Icomp gives: (Icomp-.01,phi_Icomp*(-.01)) and (Icomp+.01,phi_Icomp*.01)
# I50 gives: (I50,.5*Pgmax-Rd)
# I85 gives: (I85,.85*Pgmax-Rd)
# I90 gives: (I90,.90*Pgmax-Rd)
# I95 gives: (I95,.95*Pgmax-Rd)
# Pmax_obs gives: (2500,Pmax_obs)

# res - parameters from actual curves

# res_BM - parameters imputed using evolutionary model
# res_star - parameters imputed using means
# res_star2 - parameters imputed using random intercept

# res_i_BM - values imputed using evolutionary model
# res_i_star - values imputed using means
# res_i_star2 - values imputed using random intercept
```