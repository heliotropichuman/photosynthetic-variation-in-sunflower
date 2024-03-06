library(Rphylopars)
source("R/load_data.R")
source("R/mv.R")

# If resuming the analysis, start with the following code, and skip to line 1445
load(file = "saved_state_1_25_2024.RData")

# Irradiance should be considered at observed levels and over a continuous range
PARi <- c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500)
PARi_fine <- seq(0,6000,by = .1)
pars <- c("Pmax_obs","Pgmax","phi_I0","phi_Icomp","phi_I0_Icomp","phi_Icomp_I200","Rd","Icomp",
          "I15","I25","I85","I95",
          "PAR_0","PAR_50","PAR_100","PAR_250","PAR_500","PAR_1000","PAR_1500","PAR_2000","PAR_2500")

# NLS models of photosynthetic response, see Lobo et al. 2013
eq1 <- function(pars = 
                  c(Pgmax = 19.5,phi_I0 = .0899,Rd = 1.8),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1])
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
            # Imax_obs = PARi_fine[which.max(pred_fine[PARi_fine <= 2500])],Pmax_obs = max(pred_fine[PARi_fine <= 2500]))
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500])  
  if(return == "calc") return(calc)
  
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}

eq2 <- function(pars = 
                  c(Pgmax = 19.5,I50 = 216.4,Rd = 1.8),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1])
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

eq3 <- function(pars = 
                  c(Pgmax = 19.5,phi_I0 = .0493,Rd = 1.8),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1])
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

eq4 <- function(pars = 
                  c(Pgmax = 19.5,phi_I0 = .0493,Rd = 1.8),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1])
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

eq5 <- function(pars = 
                  c(Pgmax = 15.5,Isat = 359.2,Rd = .9),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1])
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

eq6 <- function(pars = 
                  c(Pgmax = 15.5,phi_I0 = .0493,theta = .433,Rd = .9),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1])
{
  if(!missing(dat))
  {
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
            # Imax_obs = PARi_fine[which.max(pred_fine[PARi_fine <= 2500])],Pmax_obs = max(pred_fine[PARi_fine <= 2500]),
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500],
            theta = theta)
  
  if(return == "calc") return(calc)
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}

eq8 <- function(pars = 
                  c(Pgmax = 16.2,phi_I0 = .0597,Rd = 1.3),
                dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1])
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
    # pargrid <- expand.grid(Pgmax=seq(5,60,length=10),Icomp=seq(5,150,length=10),k=exp(seq(log(.0001),log(.02),length=10)),Rd=seq(0,5,length=10))
    # out <- apply(pargrid,1,function(pars) sum(((dat$A) - (pars[["Pgmax"]]*(1-exp(-pars[["k"]]*(dat$PARi-pars[["Icomp"]])))-(pars[["Rd"]])))^2))
    
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
            # Imax_obs = PARi_fine[which.max(pred_fine[PARi_fine <= 2500])],Pmax_obs = max(pred_fine[PARi_fine <= 2500]),
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500],
            k = k)
  
  if(return == "calc") return(calc)
  
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}

eq11 <- function(pars = 
                  c(phi_I0_Icomp = .0756,beta = .0000432,gamma = .0039,Icomp = 22.6),
                 dat,
                PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                return = c("predict","calc","all")[1])
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
            # Imax_obs = PARi_fine[which.max(pred_fine[PARi_fine <= 2500])],Pmax_obs = max(pred_fine[PARi_fine <= 2500]),
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500],
            beta = beta,
            gamma = gamma)
  
  fit <- c(r2 = cor(dat$A,pred,use = "pair")^2,
           mse = mean((dat$A-pred)^2,na.rm = TRUE))
  
  
  if(return == "calc") return(calc)
  
  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}


##### Defining Parameters and Analysis Objects ####

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



##### Assessing Photosynthetic Equations #####

# Compare results from all equations at once 
eqX <- paste("eq",c(1,2,3,4,5,6,8,9,11),sep = "") 
par_names <- c(names(eq9(return = "calc")),"beta","gamma")
##  To narrow the list of parameters returned, use: 
# par_names <- par_names[-grep("Isat",par_names)]

# Create an array of results 
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

res
r2
mse

##### Plot and Print NLS Models for all Samples#####

# eqX <- eq1
# dev.off()
f <- function(eqX,i,title = "") 
{
  # Rd, Pgmax, Icomp, phi_Icomp,
  fit <- eqX(dat=dat[dat$SpeciesRepDate==inds[i],],return = "all")
  calc <- fit$calc
  
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
  
 
  plot(dat[dat$SpeciesRepDate==inds[i],2:3],main = paste(i,inds[i],title),xlim = c(0,2500),ylim = c(-4,max(dat[dat$SpeciesRepDate==inds[i],3],na.rm=TRUE) + 2))
  points(x,y,pch=19)
  if(any(is.na(x) | is.na(y)))
  {
    exclude <- which(is.na(x) | is.na(y))
    x <- x[-exclude]
    y <- y[-exclude]
  }
  
  eqX <- get("eq6")
  recon <- eqX(dat=data.frame(PARi=x,A=y),return = "all")
  
  points(PARi_fine,recon$pred_fine,type='l',col='red',lwd=2)
  points(PARi_fine,fit$pred_fine,type='l',col='blue',lty=2)
  # cat(fit$calc[["Rd"]],recon$calc[["Rd"]],"\n")
  recon
}
# dev.off()
# i=sample(length(inds),1)
for(i in 1:length(inds))
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

##### Evolutionary Models #####
# The following code takes about 2.25 hours to run #
# If it's been run before, load saved file. Otherwise, run save() functions 

# p_i <- p_i_star <- p_i_star2 <- 
p <- p_star <- p_star2 <- setNames(vector(mode = "list",length = length(unique(ind_species))),unique(ind_species))
p <- p_star <- p_star2 <- setNames(rep(list(p),9)[[1]],eqX)
# save(p_i,file = "results/multivariate_BM.RData")
# save(p_i_star,file = "results/multivariate_star.RData")
# save(p_i_star2,file = "results/multivariate_star2.RData")
load(file = "results/multivariate_BM.RData")
load(file = "results/multivariate_star.RData") # all individuals are unique
load(file = "results/multivariate_star2.RData") # all individuals retain their ID but the phylogeny isn't used

res_i_BM <- res_BM <- res_i_star <- res_star <- res_i_star2 <- res_star2 <- array(dim = c(9,length(par_names)+length(PARi),length(inds)))
dimnames(res_BM) <- dimnames(res_star) <- dimnames(res_star2) <- list(paste("eq",c(1:6,8:9,11),sep = ""),c(par_names,paste("PAR_",PARi,sep = "")),inds)
dimnames(res_i_BM) <- dimnames(res_i_star) <- dimnames(res_i_star2) <- dimnames(res_BM)
          
# big computational step vvv 
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
    par_dat_i[which(par_dat$species==unique(ind_species)[i]),colnames(par_dat)[ncol(par_dat)]] <- dat_long$PARi_2500[which(par_dat$species==unique(ind_species)[i])]
    if(any(is.na(dat_long$PARi_2500[which(par_dat$species==unique(ind_species)[i])])))
    {
      na_ind <- which(is.na(dat_long$PARi_2500[which(par_dat$species==unique(ind_species)[i])]))
      for(z in 1:length(na_ind))
      {
        par_dat_i[which(par_dat$species==unique(ind_species)[i])[na_ind],colnames(par_dat)[ncol(par_dat)]] <- max(dat_long[which(par_dat$species==unique(ind_species)[i])[na_ind],-1],na.rm = TRUE)
      }
    }
    
    try(p[[eqX[e]]][[sp]] <- phylopars(par_dat_i,tree,model='BM'),silent = TRUE)
    try(p_star[[eqX[e]]][[sp]] <- phylopars(data.frame(species = inds,par_dat_i[,-1]),tree = phytools::starTree(inds,rep(1,length(inds))),model = 'star'),silent = TRUE)
    try(p_star2[[eqX[e]]][[sp]] <- phylopars(par_dat_i,phytools::starTree(unique(ind_species),rep(1,length(unique(ind_species)))),model = 'star'),silent = TRUE)
    
    if(e>1)
    {
      try({
        temp <- phylopars(par_dat_i,tree,model='BM',phylocov_start = p[[eqX[e-1]]][[sp]]$pars[[1]],phenocov_start = p[[eqX[e-1]]][[sp]]$pars[[2]])
        if(is.null(p[[eqX[e]]][[sp]])) p[[eqX[e]]][[sp]] <- temp
        if(logLik(temp)[[1]] > logLik(p[[eqX[e]]][[sp]])[[1]]) p[[eqX[e]]][[sp]] <- temp
      },silent = TRUE)
      
      try({
        temp <- phylopars(data.frame(species = inds,par_dat_i[,-1]),tree = phytools::starTree(inds,rep(1,length(inds))),model = 'star',phylocov_start = p_star[[eqX[e-1]]][[sp]]$pars[[1]])
        if(is.null(p_star[[eqX[e]]][[sp]])) p_star[[eqX[e]]][[sp]] <- temp
        if(logLik(temp)[[1]] > logLik(p_star[[eqX[e]]][[sp]])) p_star[[eqX[e]]][[sp]] <- temp
      },silent = TRUE)
      
      try({
        temp <- phylopars(par_dat_i,phytools::starTree(unique(ind_species),rep(1,length(unique(ind_species)))),model = 'star',phylocov_start = p_star2[[eqX[e-1]]][[sp]]$pars[[1]],phenocov_start = p_star2[[eqX[e-1]]][[sp]]$pars[[2]])
        if(is.null(p_star2[[eqX[e]]][[sp]])) p_star2[[eqX[e]]][[sp]] <- temp
        if(logLik(temp)[[1]] > logLik(p_star2[[eqX[e]]][[sp]])) p_star2[[eqX[e]]][[sp]] <- temp
      },silent = TRUE)
      
    }
    
    if(i>1)
    {
      try({
        temp <- phylopars(par_dat_i,tree,model='BM',phylocov_start = p[[eqX[e]]][[unique(ind_species)[i-1]]]$pars[[1]],phenocov_start = p[[eqX[e]]][[unique(ind_species)[i-1]]]$pars[[2]])
        if(is.null(p[[eqX[e]]][[sp]])) p[[eqX[e]]][[sp]] <- temp
        if(logLik(temp)[[1]] > logLik(p[[eqX[e]]][[sp]])[[1]]) p[[eqX[e]]][[sp]] <- temp
      },silent = TRUE)
      
      try({
        temp <- phylopars(data.frame(species = inds,par_dat_i[,-1]),tree = phytools::starTree(inds,rep(1,length(inds))),model = 'star',phylocov_start = p_star[[eqX[e]]][[unique(ind_species)[i-1]]]$pars[[1]])
        if(is.null(p_star[[eqX[e]]][[sp]])) p_star[[eqX[e]]][[sp]] <- temp
        if(logLik(temp)[[1]] > logLik(p_star[[eqX[e]]][[sp]])) p_star[[eqX[e]]][[sp]] <- temp
      },silent = TRUE)
      
      try({
        temp <- phylopars(par_dat_i,phytools::starTree(unique(ind_species),rep(1,length(unique(ind_species)))),model = 'star',phylocov_start = p_star2[[eqX[e]]][[unique(ind_species)[i-1]]]$pars[[1]],phenocov_start = p_star2[[eqX[e]]][[unique(ind_species)[i-1]]]$pars[[2]])
        if(is.null(p_star2[[eqX[e]]][[sp]])) p_star2[[eqX[e]]][[sp]] <- temp
        if(logLik(temp)[[1]] > logLik(p_star2[[eqX[e]]][[sp]])) p_star2[[eqX[e]]][[sp]] <- temp
      },silent = TRUE)

    }
    
    # save.image("results/curve_recons.RData")
    # dev.off()
    for(j in 1:length(which(unique(ind_species)[i]==ind_species)))
      
    {
      ind <- inds[which(ind_species==unique(ind_species)[i])][j]
      jpeg(filename = paste("results/impute_",e,"_",i,"_",j,".jpeg",sep=""),width = 800,height = 1000)
      
      par(mfrow=c(3,3))

      # imputation from values assuming BM
      f(eqX = get(eqX[e]),i = which(ind_species==unique(ind_species)[i])[j],title = paste("\neq",e," values"," BM",sep = ""))
       
      out <- get(eqX[e])(dat = data.frame(PARi=PARi,A=unlist(p_i[[sp]]$ind_recon[which(dat_long$species==unique(ind_species)[i]),][j,-1])),return = "all")
      points(PARi_fine,out$pred_fine,type='l',lwd = 2,col = 'blue')
      
      res_i_BM[eqX[e],par_names,ind] <- out$calc[par_names]
      res_i_BM[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]
      
      # imputation from values assuming all plants independent
      f(eqX = get(eqX[e]),i = which(ind_species==unique(ind_species)[i])[j],title = paste("\neq",e," values"," ind_plants",sep = ""))
      out <- get(eqX[e])(dat = data.frame(PARi=PARi,A=unlist(p_i_star[[sp]]$anc_recon[which(dat_long$species==unique(ind_species)[i]),,drop = FALSE][j,])),return = "all")
      points(PARi_fine,out$pred_fine,type='l',lwd = 2,col = 'red')
      
      res_i_star[eqX[e],par_names,ind] <- out$calc[par_names]
      res_i_star[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]

      # imputation from values assuming all species independent      
      f(eqX = get(eqX[e]),i = which(ind_species==unique(ind_species)[i])[j],title = paste("\neq",e," values"," ind_species",sep = ""))
      out <- get(eqX[e])(dat = data.frame(PARi=PARi,A=unlist(p_i_star2[[sp]]$ind_recon[which(dat_long$species==unique(ind_species)[i]),1:(ncol(dat_long))][j,-1])),return = "all")
      points(PARi_fine,out$pred_fine,type='l',lwd = 2,col = 'purple')
      
      res_i_star2[eqX[e],par_names,ind] <- out$calc[par_names]
      res_i_star2[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]
      
      x <- c(unlist(p[[eqX[e]]][[sp]]$ind_recon[which(dat_long$species==unique(ind_species)[i]),][j,-c(1,ncol(p[[eqX[e]]][[sp]]$ind_recon))]),2500)
      x <- c(rep(x[1:3],2),rep(x[4],3),rep(x[5],4))
      y <- c(.15,.25,.85,.95,1)*unlist(p[[eqX[e]]][[sp]]$ind_recon[which(dat_long$species==unique(ind_species)[i]),][j,ncol(p[[eqX[e]]][[sp]]$ind_recon)])
      y <- c(rep(y[1:3],2),rep(y[4],3),rep(y[5],4))
      
      out <- get(eqX[e])(dat = data.frame(PARi=x,A=y),return = "all")
      res_BM[eqX[e],par_names,ind] <- out$calc[par_names]
      res_BM[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]
      
      # imputation from parameters assuming BM
      f(eqX = get(eqX[e]),i = which(ind_species==unique(ind_species)[i])[j],title = paste("\neq",e," pars"," BM",sep = ""))
      points(PARi_fine,out$pred_fine,type='l',lwd = 2,col = "blue")
      
      x <- c(unlist(p_star[[eqX[e]]][[sp]]$anc_recon[ind,,drop = FALSE][1,-ncol(p_star[[eqX[e]]][[sp]]$anc_recon)]),2500)
      x <- c(rep(x[1:3],2),rep(x[4],3),rep(x[5],4))
      y <- c(.15,.25,.85,.95,1)*unlist(p_star[[eqX[e]]][[sp]]$anc_recon[ind,,drop = FALSE][1,ncol(p_star[[eqX[e]]][[sp]]$anc_recon)])
      y <- c(rep(y[1:3],2),rep(y[4],3),rep(y[5],4))
      
      out <- get(eqX[e])(dat = data.frame(PARi=x,A=y),return = "all")
      res_star[eqX[e],par_names,ind] <- out$calc[par_names]
      res_star[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]
    
      # imputation from parameters assuming all plants independent
      f(eqX = get(eqX[e]),i = which(ind_species==unique(ind_species)[i])[j],title = paste("\neq",e," pars"," ind_plants",sep = ""))
      points(PARi_fine,out$pred_fine,type='l',lwd = 2,col = "red")
      
      x <- c(unlist(p_star2[[eqX[e]]][[sp]]$ind_recon[which(p_star2[[eqX[e]]][[sp]]$ind_recon[,1]==sp),,drop = FALSE][j,-c(1,ncol(p_star2[[eqX[e]]][[sp]]$ind_recon))]),2500)
      x <- c(rep(x[1:3],2),rep(x[4],3),rep(x[5],4))
      y <- c(.15,.25,.85,.95,1)*unlist(p_star2[[eqX[e]]][[sp]]$ind_recon[which(p_star2[[eqX[e]]][[sp]]$ind_recon[,1]==sp),,drop = FALSE][j,ncol(p_star2[[eqX[e]]][[sp]]$ind_recon)])
      y <- c(rep(y[1:3],2),rep(y[4],3),rep(y[5],4))
      
      # imputation from parameters assuming all species independent
      f(eqX = get(eqX[e]),i = which(ind_species==unique(ind_species)[i])[j],title = paste("\neq",e," pars"," ind_species",sep = ""))
      out <- get(eqX[e])(dat = data.frame(PARi=x,A=y),return = "all")
      res_star2[eqX[e],par_names,ind] <- out$calc[par_names]
      res_star2[eqX[e],paste("PAR_",PARi,sep=""),ind] <- out$pred_fine[which(PARi_fine %in% PARi)]
      
      points(PARi_fine,out$pred_fine,type='l',lwd = 2,col = 'purple')
      dev.off()
    }
  }
}
# After completing evolutionary modeling, save all results in the environemnt with the following: 
save.image(file = "saved_state_1_25_2024.RData") 



##### Assessing Evolutionary Models #####

# How well did the BM model perform? 
res_vs_res_BM <- t(sapply(1:9,function(i) sapply(1:ncol(res),function(j) cor(res[i,j,],res_i_BM[i,j,],use = "pair"))))
colnames(res_vs_res_BM) <- colnames(res)
rownames(res_vs_res_BM) <- rownames(res)

pars <- c("Pmax_obs","Pgmax","phi_I0","phi_Icomp","phi_I0_Icomp","phi_Icomp_I200","Rd","Icomp",
          "I15","I25","I85","I95",
          "PAR_0","PAR_50","PAR_100","PAR_250","PAR_500","PAR_1000","PAR_1500","PAR_2000","PAR_2500")


##### Evaluating FVT and PGNLS #####

# Which NLS model best fits photosynthetic light response curves in sunflowers?
  # fit 9 different model types; table comparing average mse and r2; representative figure
eq <- "eq6"
mse[,eq]
r2[,eq]

# visualize distribution of parameters traits 
hist(res[eq,"Pmax_obs",])
hist(res[eq,"Pgmax",])
hist(res[eq,"phi_I0",])
hist(res[eq,"phi_Icomp",])
hist(res[eq,"phi_I0_Icomp",])
hist(res[eq,"phi_Icomp_I200",]) #slightly skewed
hist(res[eq,"Rd",])
hist(res[eq,"Icomp",]) #slightly skewed
hist(res[eq,"I15",]) #slightly skewed
hist(res[eq,"I25",])
hist(res[eq,"I85",]) #slightly skewed
hist(res[eq,"I95",])
hist(res[eq,"PAR_0",])
hist(res[eq,"PAR_50",]) #slightly skewed
hist(res[eq,"PAR_100",])
hist(res[eq,"PAR_250",])
hist(res[eq,"PAR_500",]) #slightly skewed
hist(res[eq,"PAR_1000",])
hist(res[eq,"PAR_2000",])
hist(res[eq,"PAR_2500",])

# Does photosynthetic rate at different light levels exhibit phylogenetic signal?
  # test phylogenetic signal at each level
data.frame(PARi = PARi,phy_signal = 
sapply(2:ncol(dat_long),function(i) phylopars(dat_long[,c(1,i)],tree,model = 'lambda')$model$lambda))

dat_long=as.data.frame(dat_long)
vals <- sapply(1:length(inds),function(i) get(eq)(dat = dat[dat$SpeciesRepDate==inds[i],],return = "all")$pred_fine[PARi_fine%in%seq(0,2500,by = 50)])
a=sapply(1:nrow(vals),function(i) {
  phy_fit <- phylopars(data.frame(species=ind_species,trait=vals[i,]),tree,model='lambda')
  phy_fit$model$lambda
  # phy_fit$pars[[1]]/(phy_fit$pars[[1]]+phy_fit$pars[[2]])
})

# vals ???

# plots phyogenetic signal (Pagel's lambda) across the smooth function of PARi levels, as well as the signal for the observed A at each PARi
plot(seq(0,2500,by = 50),a,type='l',ylim = c(0,1))
points(PARi,sapply(2:ncol(dat_long),function(i) phylopars(dat_long[,c(1,i)],tree,model = 'lambda')$model$lambda),pch=19)


# Do parameters extracted from photosynthetic light response curves exhibit phylogenetic signal?
setNames(sapply(1:length(pars),function(i){
  phylopars(data.frame(species = ind_species,trait = t(res[eq,,])[,pars[i]],row.names = NULL),tree = tree,model = 'lambda')$model$lambda
}),pars)
# rows = equations; columns = parameters


# Can we use phylogenetic comparative methods to predict light response curves in new species given only survey measurements?
  # R2 and mse for LOO-CV

suppressWarnings(cbind(
  setNames(round(sapply(1:ncol(res),function(i) cor(t(
    res[eq,,])[,i],
    t(res_BM[eq,,])[,i]
    ))^2,5),colnames(res)),
  # left column:  prediction R2 between parameters on observed vs imputed using parameters
  # right column: prediction R2 between parameters on observed vs imputed using values
  
setNames(round(sapply(1:ncol(res),function(i) cor(t(
  res[eq,,])[,i],
  t(res_i_BM[eq,,])[,i]
  ))^2,5),colnames(res))))
# left column:  prediction MSE between parameters on observed vs imputed using parameters
# right column: prediction MSE between parameters on observed vs imputed using values

suppressWarnings(cbind(
  setNames(round(sapply(1:ncol(res),function(i) 
    mean((t(res[eq,,])[,i] - t(res_BM[eq,,])[,i])^2
         )),3),colnames(res)),
  setNames(round(sapply(1:ncol(res),function(i)
    mean((t(res[eq,,])[,i] - t(res_i_BM[eq,,])[,i])^2
  )),3),colnames(res))))

# prediction R2 for each parameter for each equation (using parameters)
r2_impute <- suppressWarnings(sapply(1:length(eqX),function(e) {
  setNames(round(sapply(1:ncol(res),function(i) cor(t(
  res[eqX[e],,])[,i],
  t(res_BM[eqX[e],,])[,i]
))^2,5),colnames(res))}))[pars,]
colnames(r2_impute) <- eqX
r2_impute

# prediction MSE for each parameter for each equation (using parameters)
mse_impute <- suppressWarnings(sapply(1:length(eqX),function(e) {
  setNames(round(sapply(1:ncol(res),function(i) 
  {
    mean((t(res[eqX[e],,])[,i] - t(res_BM[eqX[e],,])[,i])^2)
  }),5),colnames(res))[pars]
  }))
colnames(mse_impute) <- eqX
mse_impute

# prediction R2 for each parameter for each equation (using values)
r2_impute <- suppressWarnings(sapply(1:length(eqX),function(e) {
  setNames(round(sapply(1:ncol(res),function(i) cor(t(
    res[eqX[e],,])[,i],
    t(res_i_BM[eqX[e],,])[,i]
  ))^2,5),colnames(res))}))[pars,]
colnames(r2_impute) <- eqX
r2_impute

# prediction MSE for each parameter for each equation (using values)
mse_impute <- suppressWarnings(sapply(1:length(eqX),function(e) {
  setNames(round(sapply(1:ncol(res),function(i) 
  {
    mean((t(res[eqX[e],,])[,i] - t(res_i_BM[eqX[e],,])[,i])^2)
  }),5),colnames(res))[pars]
}))
colnames(mse_impute) <- eqX
mse_impute



###### Examine results on phylogeny #####
# How well can you reconstruct the tree with each parameter?

# distribution of mse 
phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = log(mse[,eq]^.5)))[,1])
# (log sqrt transformed so take this with grain of salt)

phylo_fit <- phylopars(data.frame(species = ind_species,trait = res[eq,"phi_Icomp_I200",]),tree = tree,model = 'lambda')
range(phylo_fit$anc_recon[,1])
range(res[eq,"phi_I0",])

# plot anc recon 
# note: numbers represent ancestral states
plot(tree)
nodelabels(frame = "none",col = "red")

# Pgmax
phylo_fit <- phylopars(data.frame(species = ind_species,trait = res[eq,"Pgmax",]),tree = tree,model = 'lambda')
range(phylo_fit$anc_recon[,1])
range(res[eq,"Pgmax",])

# Pmax_obs
phylo_fit <- phylopars(data.frame(species = ind_species,trait = res[eq,"Pmax_obs",]),tree = tree,model = 'lambda')
range(phylo_fit$anc_recon[,1])
range(res[eq,"Pmax_obs",])

# PAR_2500	
phylo_fit <- phylopars(data.frame(species = ind_species,trait = res[eq,"PAR_2500",]),tree = tree,model = 'lambda')
range(phylo_fit$anc_recon[,1])
range(res[eq,"PAR_2500",])



#### Traits on Trees ####

# Imputed trait distributions 
BM_parameters <- res_BM[6,,]
BM_parameters <- as.data.frame(t(BM_parameters))

phi1 <- BM_parameters$phi_I0
phi2 <- BM_parameters$phi_I0_Icomp
phi3 <- BM_parameters$phi_Icomp
phi4 <- BM_parameters$phi_Icomp_I200
Pgmax <- BM_parameters$Pgmax
Icomp <- BM_parameters$Icomp
Rd <- BM_parameters$Rd
PAR_2500 <- BM_parameters$PAR_2500
Imax <- BM_parameters$Imax

phi_I0_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = phi1))[,1], title = 'phiI0')
phi_I0_Icomp_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = phi2))[,1],)
phi_Icomp_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = phi3))[,1])
phi_Icomp_I200_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = phi4))[,1])
Pgmax_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = Pgmax))[,1])
Icomp_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = Icomp))[,1])
Rd_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = Rd))[,1])
PAR_2500_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = PAR_2500))[,1])
Imax_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = Imax))[,1])

plot(phi_I0_ContMap,fsize = c(.5,.5),leg.txt="phi_I0")
plot(phi_I0_Icomp_ContMap,fsize = c(.5,.5),leg.txt="phi_I0_Icomp")
plot(phi_Icomp_ContMap,fsize = c(.5,.5),leg.txt="phi_Icomp")
plot(phi_Icomp_I200_ContMap,fsize = c(.5,.5),leg.txt="phi_Icomp_I200")
plot(Pgmax_ContMap,fsize = c(.5,.5),leg.txt="Pgmax")
plot(Icomp_ContMap,fsize = c(.5,.5),leg.txt="Icomp")
plot(Rd_ContMap,fsize = c(.5,.5),leg.txt="Rd")
plot(PAR_2500_ContMap,fsize = c(.5,.5),leg.txt="PAR_2500")
plot(Imax_ContMap,fsize = c(.5,.5),leg.txt="Imax")

# measured trait distributions 
measured_parameters <- res[6,,]
measured_parameters <- as.data.frame(t(measured_parameters))

res_phi1 <- measured_parameters$phi_I0
res_phi2 <- measured_parameters$phi_I0_Icomp
res_phi3 <- measured_parameters$phi_Icomp
res_phi4 <- measured_parameters$phi_Icomp_I200
res_Pgmax <- measured_parameters$Pgmax
res_Icomp <- measured_parameters$Icomp
res_Rd <- measured_parameters$Rd
res_PAR_2500 <- measured_parameters$PAR_2500
res_Imax <- measured_parameters$Imax

res_phi_I0_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = res_phi1))[,1], title = 'phiI0')
res_phi_I0_Icomp_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = res_phi2))[,1],)
res_phi_Icomp_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = res_phi3))[,1])
res_phi_Icomp_I200_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = res_phi4))[,1])
res_Pgmax_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = res_Pgmax))[,1])
res_Icomp_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = res_Icomp))[,1])
res_Rd_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = res_Rd))[,1])
res_PAR_2500_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = res_PAR_2500))[,1])
res_Imax_ContMap <- phytools::contMap(tree,Rphylopars:::convert_to_means(data.frame(species = ind_species,rank = res_Imax))[,1])

plot(res_phi_I0_ContMap,fsize = c(.5,.5),leg.txt="phi_I0")
plot(res_phi_I0_Icomp_ContMap,fsize = c(.5,.5),leg.txt="phi_I0_Icomp")
plot(res_phi_Icomp_ContMap,fsize = c(.5,.5),leg.txt="phi_Icomp")
plot(res_phi_Icomp_I200_ContMap,fsize = c(.5,.5),leg.txt="phi_Icomp_I200")
plot(res_Pgmax_ContMap,fsize = c(.5,.5),leg.txt="Pgmax")
plot(res_Icomp_ContMap,fsize = c(.5,.5),leg.txt="Icomp")
plot(res_Rd_ContMap,fsize = c(.5,.5),leg.txt="Rd")
plot(res_PAR_2500_ContMap,fsize = c(.5,.5),leg.txt="PAR_2500")
plot(res_Imax_ContMap,fsize = c(.5,.5),leg.txt="Imax")


#### Comparing observed and predicted FVT ####

# compare evolutionary rates (annuals and per)
# heatmap
# sensitivity
#
# correlations with genome size
#
# seasonality

