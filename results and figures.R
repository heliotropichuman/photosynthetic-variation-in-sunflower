library(Rphylopars)
library(xlsx)
load(file = "saved_state_1_25_2024.RData")
source("R/load_data.R")
source("R/mv.R")


## Equation of best fit ####

# Which model type best fits photosynthetic light response curves in sunflowers?
# fit 9 different model types; table comparing average mse and r2; representative figure

# The goal: Make a heat map matrix of R2 and mse for all equations 
r2
mse
inds

# attempt 0 
res <- array(dim = c(9,length(par_names)+length(PARi),length(inds)))
dimnames(res) <- list(paste("eq",c(1:6,8:9,11),sep = ""),c(par_names,paste("PAR_",PARi,sep="")),inds)
ErrorMatrix <- array(ErrorMatrix, dim = length(mse)+length(r2), dimnames = list(inds,eqX))


ErrorMatrix <- array(dim = c(9,length(mse)+length(r2),length(inds)))
dimnames(ErrorMatrix) <- list(paste("eq",c(1:6,8:9,11),sep = ""),c(mse, r2,inds))

#Attempt 1
image(1:ncol(mse), 1:nrow(r2), t(mse), col = terrain.colors(60), axes = FALSE)
axis(1, 1:ncol(mse), colnames(mse))
axis(2, 1:nrow(r2), rownames(r2))
for (x in 1:ncol(mse))
  for (y in 1:nrow(r2))
    text(x, y, m[y,x])

# attampt 2
library(ggplot2)
heatmap(mse,r2)
heatmap(r2)
heatmap(mse, Rowv= NA, Colv=NA)
heatmap(r2, Rowv= NA, Colv=NA)
heatmap(mse, Rowv= NA, Colv=NA, main = "Mean Squared Error", xlab = "Model Equation", ylab = "Sample")
heatmap(r2, Rowv= NA, Colv=NA)

# attempt 3
install.packages("ggcorrplot")
library(ggcorrplot)
ggcorrplot(cor(mse, r2))


# attempt 4
ErrorMatrix<- rbind(mse, r2)

eq <- "eq6"
mse[,eq]
r2[,eq]

results <- res[eq,,]
average_parameters <- apply(results, length(dat_long), mean)



## Phylogenetic Signal in the Curve #####
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

# plots phyogenetic signal (Pagel's lambda) across the smooth function of PARi levels, as well as the signal for the observed A at each PARi
plot(seq(0,2500,by = 50),a,type='l',ylim = c(0,1))
points(PARi,sapply(2:ncol(dat_long),function(i) phylopars(dat_long[,c(1,i)],tree,model = 'lambda')$model$lambda),pch=19)



## Phylogenetic Signal in the Parameters #####
# Do parameters extracted from photosynthetic light response curves exhibit phylogenetic signal?
# rows = equations; columns = parameters
setNames(sapply(1:length(pars),function(i){
  phylopars(data.frame(species = ind_species,trait = t(res[eq,,])[,pars[i]],row.names = NULL),tree = tree,model = 'lambda')$model$lambda
}),pars)




## PGLS with Survey Measurements #####
# Can we use phylogenetic comparative methods to predict light response curves in new species given only survey measurements?
# R2 and mse for LOO-CV

# left column:  prediction R2 between parameters on observed vs imputed using parameters
# right column: prediction R2 between parameters on observed vs imputed using values
suppressWarnings(cbind(
  setNames(round(sapply(1:ncol(res),function(i) cor(t(
    res[eq,,])[,i],
    t(res_BM[eq,,])[,i]
  ))^2,5),colnames(res)),
  
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



## Traits on Trees ####

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

# font size was adjusted from fsize = c(.75,.5) to c(.5,.5)
# When copying plot, preview window was set to 450 wide and 400 tall

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


## LRC Comparisons ####

# treat A at each PARi level as a "trait"
PARi <- c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500)
PARi_fine <- seq(0,2500,length=500)
species_and_nodes <- anc_recon$species

lightcurve <- function(pars,PARi)
{
  pars <- as.numeric(pars)
  phi <- pars[1]
  beta <- pars[2]
  gamma <- pars[3]
  Icomp <- pars[4]
  phi*((1-beta*PARi)/(1+gamma*PARi))*(PARi-Icomp)
}

LH <- data.frame(species = c("H.agrestis", "H.annuus", "H.argophyllus", "H.atrorubens", 
                             "H.cusickii", "H.debilis", "H.divaricatus", "H.floridanus", "H.giganteus", 
                             "H.gracilentus", "H.grosseserratus", "H.heterophyllus", "H.longifolius", 
                             "H.maximiliani", "H.microcephalus", "H.mollis", "H.neglectus", 
                             "H.nuttalliissp.nuttallii", "H.occidentalis", "H.petiolaris", 
                             "H.praecox", "H.radula", "H.silphioides", "H.winteri"),
                 clade = c("NA", "A", "A", "SE", "LP", "A", "LP", "SE", "LP", "P", "LP", 
                           "SE", "SE", "LP", "LP", "SE", "A", "LP", "P", "A", "A", "SE", 
                           "SE", "A"),
                 lifehist_2cats = c("A", "A", "A", "P", "P", "A", "P", "P", "P", "P", "P", "P", 
                                    "P", "P", "P", "P", "A", "P", "P", "A", "A", "P", "P", "A"),
                 lifehist_3cats = c("A", "A", "A", "SE", "LP", "A", "LP", "SE", "LP", "P", "LP", 
                                    "SE", "SE", "LP", "LP", "SE", "A", "LP", "P", "A", "A", "SE", 
                                    "SE", "A"), #replace ) with , to include lifehist_4cats
                 lifehist_4cats = c("Annual", "Annual", "Annual", "Erect Perennial", "Erect Perennial", "Annual", "Erect Perennial", "Erect Perennial", "Erect Perennial", "Erect Perennial", "Erect Perennial", "Basal Rosette", 
                                    "Basal Rosette", "Erect Perennial", "Erect Perennial", "Erect Perennial", "Annual", "Erect Perennial", "Basal Rosette", "Annual", "Annual", "Basal Rosette", "Erect Perennial", "Annual"))


### Plot NLS #### 
plot_curve <- function(current_node,values = TRUE,col = "black",lwd = 3,type = 'l',lty = 1,add = FALSE,anc_recon = BM_mv_anc_recon)
{
  if(values)
  {
    current_A <- as.numeric(anc_recon[anc_recon$species == current_node,-1])
    current_dat <- data.frame(PARi = PARi,A = current_A)
    current_mod <- nls(A~phi*((1-beta*PARi)/(1+gamma*PARi))*(PARi-Icomp),data = current_dat,start = c(phi=.0756,beta=0,gamma=.0039,Icomp=22.5))
    
    fitted_A <- predict(current_mod,newdata = data.frame(PARi = PARi_fine))
  } else
  {
    current_mod <- unlist(nls_BM_recon[nls_BM_recon$species == current_node,-1])
    fitted_A <- with(data = as.list(current_mod),expr = phi*((1-beta*PARi_fine)/(1+gamma*PARi_fine))*(PARi_fine-Icomp))
  }
  plot_fn <- if(add) points else plot
  plot_fn(x = PARi_fine,
          y = fitted_A,
          type = type, # type = 'l' means to plot a line instead of points
          xlab = "PARi",
          ylab = "A",
          main = current_node, # title
          ylim = range(anc_recon[,-1]), # make the y plotting span the range of all observed values so curves are comparable
          lwd = lwd, # line width (default = 1, 3 is thicker, etc)
          lty = lty, # line type: options include solid (lty = 1) dashed(lty = 2) dotted(lty = 3)
          col = col) # set color of line
}

# example: plot H. annuus mean curve based on values
plot_curve("H.annuus",values = TRUE,col = "blue",lwd = 3,lty = 2)
# now overlay the H. annuus curve based on mean nls coefficients (overlay by setting add = TRUE)
plot_curve("H.annuus",values = FALSE,col = "red",lwd = 3,lty = 2,add = TRUE)

root <- as.character(length(tree$tip.label)+1)

# do the same thing for the root node of the tree: 28
plot_curve(root,values = TRUE,col = "blue",lwd = 3,lty = 2)
plot_curve(root,values = FALSE,col = "green",lwd = 3,lty = 2,add = TRUE)

# Example for loop to plot every curve individually (including ancestral nodes)
for(i in 1:length(species_and_nodes))
{
  current_node <- species_and_nodes[i]
  plot_curve(current_node)
}

### Plot curves with LH #### 
# (treating A at a given PARi as a trait)

# raw species average
merge(dat %>% filter(complete.cases(.)),LH,by="species") %>% 
  filter(species %in% tree$tip.label) %>% 
  group_by(species) %>% 
  summarize(
    A = predict(
      nls(formula = A~phi*((1-beta*PARi)/(1+gamma*PARi))*(PARi-Icomp),
          start = c(phi=.0756,beta=0,gamma=.0039,Icomp=22.5)),
      newdata = data.frame(PARi=PARi_fine)),
    LH = lifehist_2cats[1]) %>% 
  mutate(PARi=PARi_fine) %>% 
  ungroup() %>%
  ggplot(mapping=aes(x = PARi,y = A,color=LH,shape=species)) + 
  stat_smooth(se = FALSE) + 
  theme_classic(base_size = 15) + 
  labs(color = "Life History")

# "star phylogeny" (shrinkage) reconstruction based on treating A at each PARi as separate traits
star_mv_anc_recon_long |> 
  merge(y = LH,by="species") |> 
  filter(species %in% tree$tip.label) |> 
  group_by(species) |> 
  summarize(
    A = predict(
      nls(formula = A~phi*((1-beta*PARi)/(1+gamma*PARi))*(PARi-Icomp),
          start = c(phi=.0756,beta=0,gamma=.0039,Icomp=22.5)),
      newdata = data.frame(PARi=PARi_fine)),
    LH = lifehist_2cats[1]) |> 
  mutate(PARi=PARi_fine) |> 
  ungroup() |>
  ggplot(mapping=aes(x = PARi,y = A,color=LH,shape=species)) + 
  stat_smooth(se = FALSE) + 
  theme_classic(base_size = 15) + 
  labs(color = "Life History")

# BM reconstruction based on treating A at each PARi as separate traits
BM_mv_anc_recon_long |> 
  merge(y = LH,by="species") |> 
  filter(species %in% tree$tip.label) |> 
  group_by(species) |> 
  summarize(
    A = predict(
      nls(formula = A~phi*((1-beta*PARi)/(1+gamma*PARi))*(PARi-Icomp),
          start = c(phi=.0756,beta=0,gamma=.0039,Icomp=22.5)),
      newdata = data.frame(PARi=PARi_fine)),
    LH = lifehist_2cats[1]) |> 
  mutate(PARi=PARi_fine) |> 
  ungroup() |>
  ggplot(mapping=aes(x = PARi,y = A,color=LH,shape=species)) + 
  stat_smooth(se = FALSE) + 
  theme_classic(base_size = 15) + 
  labs(color = "Life History")

### Plot curves with LH predicted from parameters #### 
# (treating the nls regression parameters as "traits")

# raw species averages
raw_curves_fine <- data.frame(species=tree$tip.label,A=t(sapply(1:length(tree$tip.label),function(X) rowMeans(apply(nls_coefs[nls_coefs$species==tree$tip.label[X],-1,drop=FALSE],1,function(X) lightcurve(X,PARi_fine)))))) |> pivot_longer(cols=-1) |> mutate(PARi = PARi_fine[as.numeric(gsub("A.","",name))]) |> rename(A=value) |> select(-name) |> select(species,PARi,A)

merge(raw_curves_fine,LH,by="species") |> 
  ggplot(mapping = aes(x = PARi,y = A,color = lifehist_2cats,shape=species)) + 
  stat_smooth(se = FALSE) +
  theme_classic(base_size = 15) + 
  labs(color = "Life History")

# phylo signal
# heritability

# multivariate trait (reconstruct light curves from raw data at each PARi)
# fit curves to reconstructed A's at each PARi (traits = A at each PARi, Amax, Rd, phi, Icomp)

# nls curves (reconstruct parameters from fitted parameters)


# "star" (shirnkage) reconstruction of regression parameteres
star_curves_fine <- data.frame(species=tree$tip.label,A=t(sapply(1:length(tree$tip.label),function(X) rowMeans(apply(nls_star_recon[nls_star_recon$species==tree$tip.label[X],-1,drop=FALSE],1,function(X) lightcurve(X,PARi_fine)))))) |> pivot_longer(cols=-1) |> mutate(PARi = PARi_fine[as.numeric(gsub("A.","",name))]) |> rename(A=value) |> select(-name) |> select(species,PARi,A)

merge(star_curves_fine,LH,by="species") |> 
  ggplot(mapping = aes(x = PARi,y = A,color = lifehist_2cats,shape=species)) + 
  stat_smooth(se = FALSE) +
  theme_classic(base_size = 15) + 
  labs(color = "Life History")

# BM reconstruction of regression parameteres
BM_curves_fine <- data.frame(species=tree$tip.label,A=t(sapply(1:length(tree$tip.label),function(X) rowMeans(apply(nls_BM_recon[nls_BM_recon$species==tree$tip.label[X],-1,drop=FALSE],1,function(X) lightcurve(X,PARi_fine)))))) |> pivot_longer(cols=-1) |> mutate(PARi = PARi_fine[as.numeric(gsub("A.","",name))]) |> rename(A=value) |> select(-name) |> select(species,PARi,A)

merge(BM_curves_fine,LH,by="species") |> 
  ggplot(mapping = aes(x = PARi,y = A,color = lifehist_2cats,shape=species)) + 
  stat_smooth(se = FALSE) +
  theme_classic(base_size = 15) + 
  labs(color = "Life History")


## Making Correlations ####
m2016a <- read.xlsx(file = "processed_data/mason_data.xlsx",sheetName = "m2016a")
m2016a_nozero <- read.xlsx(file = "processed_data/mason_data.xlsx",sheetName = "m2016a_nozero")

get_cors <- function(dataset,use_minR2 = TRUE)
{
  photodat <- data.frame(species = tree$tip.label,Rd_from_coefs = Rd_from_coefs[tree$tip.label],Rd_from_values = Rd_from_values[tree$tip.label],Amax_from_coefs = Amax_from_coefs[tree$tip.label],Amax_from_values = Amax_from_values[tree$tip.label],phi_from_coefs = phi_from_coefs[tree$tip.label],Icomp_from_coefs = Icomp_from_coefs[tree$tip.label])
  dataset <- merge(photodat,dataset,by="species")
  rownames(dataset) <- dataset$species
  
  cor_results <- data.frame(photo_trait = character(),predictor = character(),photo_col = numeric(),predictor_col = numeric(),R2 = numeric(),minR2 = numeric(),pval = numeric(),slope = numeric(),intercept = numeric(),lambda = numeric())
  
  library(phylolm)
  for(i in 2:ncol(photodat))
  {
    for(j in (ncol(photodat)+1):(ncol(dataset)))
    {
      if(colnames(dataset)[j]=="growth_form") break
      mod <- phylolm(as.formula(paste(colnames(dataset)[i],"~",colnames(dataset)[j])),phy=tree,data=dataset,model='lambda')
      if(mod$r.squared>.25)
      {
        minR2 <- round(min(sapply(dataset[!is.na(dataset[,j]),"species"],function(tip) phylolm(as.formula(paste(colnames(dataset)[i],"~",colnames(dataset)[j])),phy=drop.tip(tree,tip),data=dataset,model='lambda')$r.squared)),3)
        if(use_minR2) if(minR2 < .25) next
        if(!use_minR2) if(mod$r.squared < .25) next
        cor_results <- rbind(cor_results,
                             data.frame(photo_trait = colnames(dataset)[i],predictor = colnames(dataset)[j],photo_col = i,predictor_col = j,R2 = round(mod$r.squared,3),minR2 = minR2,pval = round(summary(mod)$coef[2,4],5),slope = round(mod$coefficients[2],7),intercept = round(mod$coefficients[1],4),lambda = round(mod$optpar,2)))
        plot(dataset[,j],dataset[,i],xlab = colnames(dataset)[j],ylab = colnames(dataset)[i],pch=19,main=paste("R2 =",round(mod$r.squared,3),"- minR2 =",minR2))
        abline(mod)
      }
    }
  }
  cor_results
}
# get_cors(m2016a)
# get_cors(m2016a_nozero,use_minR2 = TRUE)
get_cors(m2016a_nozero,use_minR2 = TRUE)

########################################################

# access fitted curve parameters for a given plant
res[,,"Agrestis_1_29/10/19"]


