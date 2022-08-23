
#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestión Moderna de Portafolio
# Capitulo 5: CVaR 
# Copyright 2022
#-----------------------------------------------------------
#-----------------------------------------------------------
#
# Notas:
# 1. Requiere cargar primero la funcion para importar precios  
#    de Yahoo Finance: "func_precios"
#-----------------------------------------------------------
#-----------------------------------------------------------

# Ejemplo 5.2
# Medida CVaR

rm(list=ls())

library(zoo)
library(xts)
library(quantmod)
library(quadprog)
library(ROI)

## Acciones seleccionadas
activos <- c("AAPL","AMZN","GOOG","MSFT")
fechai <- '2009-12-01'
fechaf <- '2021-12-31'
periodicidad <- "monthly"   

precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,]

mu <- colMeans(retornos)
cov <- cov(retornos)


# PMVG de Markowitz
n <- ncol(retornos) 
ones <- rep(1,n)

cov_inv_1 <- solve(cov, ones) 
(wpmvg <- (1/as.numeric(ones %*% cov_inv_1)) * cov_inv_1)
(rpmvg <- mu%*%wpmvg)
(sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg))


##--------------------------------------------------------------
## Modelo CVaR
##--------------------------------------------------------------

## Usando ROI y glpk

library("ROI")
cvarOpt = function(retornos, alpha=0.05, rmin=0, wmin=0, wmax=1, w.sum=1)
{
    require(Rglpk)
    s = nrow(retornos)
    Amat = rbind(cbind(rbind(1,mu),matrix(data=0,nrow=2,ncol=s+1)),
                 cbind(retornos,diag(s),1))
    objL = c(rep(0,n), rep(-1/(alpha*s), s), -1)
    bvec = c(w.sum,rmin,rep(0,s))
    # Direction vector
    dir.vec = c("==",">=",rep(">=",s))
    # Box: bounds on weights
    bounds = list(lower = list(ind = 1:n, val = rep(wmin,n)),
                  upper = list(ind = 1:n, val = rep(wmax,n)))
    res = Rglpk_solve_LP(obj=objL, mat=Amat, dir=dir.vec, rhs=bvec,
                         types=rep("C",length(objL)), max=T, bounds=bounds)
    w = as.numeric(res$solution[1:n])
    return(list(w=w,status=res$status))
}

solcvar <- cvarOpt(coredata(retornos), alpha=0.05, rmin=0, wmin=0, wmax=1, w.sum=1)
(wpcvar <- solcvar[["w"]])

# --------------------------------------------------------------------------------------------------

## Fuente: Vana, Schwendinger y Hochreiter (2016)
## ROML Package:

if (!require("ROI")) install.packages("ROI"); library("ROI")
if (!require("ROML")) install.packages("ROML"); library("ROML")
if (!require("ROML.portfolio")) install.packages("ROML.portfolio"); library("ROML.portfolio")
if (!require("ROI.plugin.glpk")) install.packages("ROI.plugin.glpk"); library("ROI.plugin.glpk")
if (!require("ROI.plugin.quadprog")) install.packages("ROI.plugin.quadprog"); library("ROI.plugin.quadprog")

## CVaR
m <- model()
m$variable(portfolio, lb = 0) # the portfolio choice vector; 
m$minimize( cvar(portfolio, 0.95) )
m$subject_to( budget_norm(portfolio) )
opt <- optimize(m, solver="glpk", 
                data=list(returns = coredata(retornos))) 
opt[["solution"]][1:4]

## --------------------------------------------------------------
## Comparacion CVaR y PMVG

pesos <- cbind(wpcvar,wpmvg)

windows()
barplot(t(pesos),beside=TRUE,ylim=c(0,0.8),cex.names = 0.85,ylab="Part. (%)")
legend("topright",legend = c("CVaR","Markowitz"),
       fill=c("black","gray"),cex=0.8)

## ----------------------------
## Portafolio CVaR
rpcvar <- mu%*%wpcvar
sigmapcvar <- sqrt(t(wpcvar)%*%cov%*%wpcvar)

# Plano Riesgo-Retorno
windows()
plot(sigma,mu, xlim=c(min(sigma)*0.8,max(sigma)*1.1), 
     ylim=c(min(mu)*0.8,max(mu)*1.1),
     cex=0.6,xlab="Riesgo",ylab="Retorno esperado",bty="L")
points(sigmapo,rpo,type='l')
text(sigma,mu,labels=activos,pos = 4, cex=0.7)
points(sigmapmvg,rpmvg, col="black",pch=20)
text(sigmapmvg,rpmvg,labels="PMVG",pos = 2, cex=0.7)
text(sigmapcvar,rpcvar,labels="CVaR",pos = 4, cex=0.7)
points(sigmapcvar,rpcvar, col="black",pch=19)

legend("bottomright",legend = c("FE con cortos","FE sin cortos"),
       lty=c(2,1),cex=0.9)


