library(gofNetwork)
library(mixer)
library(doMC); 

rm(list=ls())
load("ZacharyTernary.Rdata")

 Y = ZacharyTernarySBM$Net+0
 X = ZacharyTernarySBM$EdgeCovar

id_alone = which(rowSums(Y) + colSums(Y) == 0) # find vertices not connected to the network

if(length(id_alone)) { # not connected vertices are removed to avoid numerical issues in gofNetwork
  Y_ = Y[-id_alone, -id_alone]
  X_ = X[-id_alone, -id_alone, ]
} else{
  Y_ = Y
  X_ = X
}

qmin<-1
qmax<-10

# gof2 : with covariates
res<-gofNetwork(Y_, X_, dim(X_)[3], qmin, qmax, maxit=100, epsconv=1e-6, nbrepeat=1, ncores=1)
plot(res)

pr<-computeProbBayes(res)
cat("gof2 (with covariates) : estimation of p(H_0 | Y) : ", pr$pH0_Y, "\n")

m<-computegraphon(res, L=100, nbsim=1000, prob=TRUE)
plot.graphon(m, zlim=c(0,1))

# mixer : no covariates
res2<-mixer(Y_, qmin=qmin, qmax=qmax, method="bayesian")
plot(res2, frame=c(1,2))

pr2<-computeProbBayes(res2)
cat("mixer (no covariates) : estimation of p(H_0 | Y) : ", pr2$pH0_Y, "\n")

m2<-computegraphon(res2, L=100, nbsim=1000, prob=TRUE)
plot.graphon(m2)
