
#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestion Moderna de Portafolio
# Capitulo 5: Medidas alternativas de riesgo - Semivarianza
# Copyright 2022
#-----------------------------------------------------------
#-----------------------------------------------------------
#
# Notas:
# 1. Requiere cargar primero la funcion para importar precios  
#    de Yahoo Finance: "func_precios"
#-----------------------------------------------------------
#-----------------------------------------------------------

# Ejemplo 5.1
# Medida Semivarianza 

rm(list=ls())
library(zoo)
library(xts)
library(quantmod)
library(quadprog)

## Acciones seleccionadas
activos <- c("AAPL","AMZN","GOOG","MSFT")
fechai <- '2009-12-01'
fechaf <- '2021-12-31'
periodicidad <- "monthly"   

precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,]

mu <- colMeans(retornos)*12
cov <- cov(retornos)*12
var <- diag(cov)
sigma <- sqrt(var)
corr <- cor(retornos)

## Semi-covarianza
h <- 0
semiret <- pmin(retornos,h)
semicov <- cov(semiret)*12
semivar <- diag(semicov)
semisigma <- sqrt(semivar)

## Diferencias MV y Media-Semivarianza

windows()
plot(sigma,mu, xlim=c(0,max(sigma)*1.2), ylim=c(0,max(mu)*1.2),cex=0.6,
     xlab="Riesgo",ylab="Retorno",font = 6)
points(semisigma,mu,pch=19)

## ---------------------------------------------------------
## FE de Markowitz
n <- ncol(retornos) 
ones <- rep(1,n)
x <- t(mu)%*%solve(cov)%*%mu
y <- t(mu)%*%solve(cov)%*%ones
z <- t(ones)%*%solve(cov)%*%ones
d <- x*z - y*y
g <- (solve(cov,ones)%*%x-solve(cov,mu)%*%y)%*%solve(d)
h <- (solve(cov,mu)%*%z-solve(cov,ones)%*%y)%*%solve(d)

# Ejempo: Rp: 0.12
rp <- 0.25
woptim <- g + h*rp

# Para toda la FE
rpmin <- min(mu)
rpmax <- max(mu) 
nport <- 100

j <- seq(rpmin,rpmax, length=nport) 
wpo <- matrix(c(0), ncol=n, nrow=nport) 
rpo <- matrix(c(0), nrow=nport)
sigmapo <- matrix(c(0), nrow=nport)
wj <- 0

for(i in 1:nport){
    wj <- g + h*j[i] 
    wpo[i,] <- t(wj)
    rpo[i,] <- t(wj)%*%mu
    sigmapo[i,] <- sqrt(t(wj)%*%cov%*%wj)
}

# ------------------------------------------------------
# Solucion Maximo Sharpe
rf <- 0.0
r <- mu -rf 
Z <- solve(cov,r)
sumaZ <- sum(Z)
wpt <- Z/sumaZ

rpt <- t(wpt)%*%mu
sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)

windows()
plot(sigma,mu, xlim=c(0,max(sigma)*1.2), ylim=c(0.1,max(mu)*1.2),cex=0.6,
     xlab="Riesgo",ylab="Retorno",font = 6)
points(sigmapo,rpo,type='l',lty=2)
points(sigmapt,rpt,pch=19,cex=1)
text(sigmapt,rpt,labels="T",pos = 2, cex=1)

## --------------------------------------------------------
## --------------------------------------------------------

## Portafolio Optimo Semivarianza
n <- ncol(retornos) 
ones <- rep(1,n)
xs <- t(mu)%*%solve(semicov)%*%mu
ys <- t(mu)%*%solve(semicov)%*%ones
zs <- t(ones)%*%solve(semicov)%*%ones
ds <- xs*zs - ys*ys
gs <- (solve(semicov,ones)%*%xs-solve(semicov,mu)%*%ys)%*%solve(ds)
hs <- (solve(semicov,mu)%*%zs-solve(semicov,ones)%*%ys)%*%solve(ds)

# Ejempo: Rp: 0.12
rp <- 0.25
woptims <- gs + hs*rp

rpmin <- min(mu)
rpmax <- max(mu)
nport <- 100

j <- seq(rpmin,rpmax, length=nport) 
wpos <- matrix(c(0), ncol=n, nrow=nport) 
rpos <- matrix(c(0), nrow=nport)
sigmapos <- matrix(c(0), nrow=nport)
wj <- 0

for(i in 1:nport){
    wj <- gs + hs*j[i] 
    wpos[i,] <- t(wj)
    rpos[i,] <- t(wj)%*%mu
    sigmapos[i,] <- sqrt(t(wj)%*%semicov%*%wj)
}

# Modelo de Sortino
rf <- rf
r <- mu -rf
Zs <- solve(semicov,r)
sumaZs <- sum(Zs)
wpts <- Zs/sumaZs

rpts <- t(wpts)%*%mu
sigmapts <- sqrt(t(wpts)%*%semicov%*%wpts)

## --------------------------------------------------------
##  Figuras

# Frontera eficiente: FE

windows()
plot(sigma,mu, xlim=c(0,max(sigma)*1.2), ylim=c(0.1,max(mu)*1.2),cex=0.6,
     xlab="Riesgo",ylab="Retorno",font = 6)
points(semisigma,mu,pch=19,cex=0.7)
points(sigmapo,rpo,type='l',lty=2)
points(sigmapos,rpos,type='l')
points(sigmapt,rpt,pch=19,cex=1)
text(sigmapt,rpt,labels="T",pos = 2, cex=1)
points(sigmapts,rpts,pch=19,cex=1)
text(sigmapts,rpts,labels="S",pos = 2, cex=1)
legend("bottomright",legend = c("FE Markowitz","FE Semivarianza"),
       lty=c(2,1),cex=0.8)

# Pesos optimos del portafolio y el tangente 

woptimo <- cbind(woptim,woptims)
rownames(woptimo) <- activos

windows()
barplot(t(woptimo), beside=TRUE,ylim=c(-0.2,0.9),cex.names = 0.85,
        ylab="Part. (%)")
legend("topright",legend = c("Markowitz","Semivarianza"),
       fill=c("black","gray"),cex=0.8)

# comparacion Sharpe y Sortino
wtangente <- cbind(wpt,wpts)
windows()
barplot(t(wtangente ), beside=TRUE,ylim=c(-0.2,0.8),cex.names = 0.85,
        ylab="Part. (%)")
legend("topright",legend = c("Sharpe","Sortino"),
       fill=c("black","gray"),cex=0.8)

#-------------------------------------------------------
# Programa QP - Restricciones en los cortos
library(quadprog)

# Portafolio Minima-Varianza de Markowitz
nport <- 100
Dmat <- cov*2
Amat <- cbind(mu,rep(1,n),diag(1,nrow=n))
dvec <- rep(0,n)

j <- seq(min(mu)*1.001,max(mu)*0.999,length=nport)

wpo3 <- matrix(0,nrow=nport,ncol=n)
sigmapo3 <- matrix(0,nrow=nport)
rpo3 <- matrix(0,nrow=nport)

for(i in 1:nport){
    bvec <- c(j[i],1,rep(0,n)) # restricciones de igualdad del sistema
    result <- solve.QP(Dmat, dvec, Amat, bvec, meq=2) #meq=restriccones de igualdad
    wj <- result[["solution"]]
    wpo3[i,] <- wj
    rpo3[i] <- mu%*%wj
    sigmapo3[i] <- sqrt(t(wj)%*%cov%*%wj)
}

#--------------------------------------------------------
# Portafolio Tangente 

sharpe_port <- (rpo3-rf)/sigmapo3
maxsharpe <- max(sharpe_port)
dif <- diff(sharpe_port)
length(dif[dif > 0])

windows()
plot(sharpe_port, type='l',xlab="Portafolio",ylab="Coef. Sharpe",
     ylim=c(min(sharpe_port),max(sharpe_port)*1.1))
points(length(dif[dif > 0]),maxsharpe,pch=19)
text(length(dif[dif > 0]),maxsharpe,labels="Máx. Sharpe",pos = 3, cex=1)


tabla <- cbind(sharpe_port,rpo3,sigmapo3)
colnames(tabla) <- c("Sharpe","Retorno","Riesgo")
colnames(wpo) <- activos

tabla2 <- cbind(tabla,wpo3)

sort_tabla2 <- tabla2[order(-tabla2[,1]),]
port_maxsharpe <- round(sort_tabla2[1,],4)

wpt
wpt2 <- port_maxsharpe[4:length(port_maxsharpe)]
rpt2 <- port_maxsharpe[2]
sigmapt2 <- port_maxsharpe[3]

# Plano Riesgo-Retorno
windows()
plot(sigma,mu, xlim=c(min(sigma)*0.6,max(sigma)*1.1), 
     ylim=c(min(mu)*0.8,max(mu)*1.2),
     cex=0.6,xlab="Riesgo",ylab="Retorno esperado")
points(sigmapo,rpo,type='l',lty=2,lwd=2)
points(sigmapt,rpt,pch=19)
points(sigmapt2,rpt2,pch=19)
text(sigma,mu,labels=activos,pos = 4, cex=0.7)
text(sigmapt,rpt,labels="T",pos = 2, cex=1)
text(sigmapt2,rpt2,labels="T*",pos = 4, cex=1)
points(sigmapo3,rpo3,col="black",lwd=2,type="l")
legend("bottomright",legend = c("FE con cortos","FE sin cortos"),
       lty=c(2,1),cex=0.9)

names(wpt2) <- activos 

windows()
barplot(t(wpt2), beside = TRUE, cex.names = 0.75,
        ylab="Part. (%)")

# ----------------------------------------------------------------
# ----------------------------------------------------------------
