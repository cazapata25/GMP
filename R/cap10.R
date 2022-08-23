#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestión Moderna de Portafolio
# Capitulo 10: Portafolio Internacional
# Copyright 2022
#-----------------------------------------------------------
#-----------------------------------------------------------

# Ejemplo 10.1
## Retorno sistematico: BR + X
## Varianza sistematica: V(BR+X)

rm(list=ls())
library(zoo)
library(xts)
library(quantmod)

# Informacion para los mercados locales
fechai <- '2015-12-01'
fechaf <- '2020-12-31'
periodicidad <- "monthly"   
activos <- c('^DJI','^N100','^N225','USDEUR=X','USDJPY=X') 
#activos <- c('^FTSE','^GSPC','^GDAXI','GBPUSD=X','GBPEUR=X')

precios.mk <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.mk))[-1,]

(rx <- colMeans(retornos)*12)
(cov <- cov(retornos)*12)
var <- diag(cov)
sigma <- sqrt(var)
(corr <- cor(retornos))
    
## Importar la info descargada desde la web-page de MSCI
# library(readr)
# msci <- read_csv(".../msci.csv")
# msci <- as.matrix(msci[,2])
# datos <- cbind(precios.mk,msci)
# msci <- datos[,6]  

r.indice <- diff(log(msci))[-1,]
(rm = mean(r.indice)*12)
(sigmam <- sd(r.indice)*sqrt(12))

n <- length(rx)
betas <- matrix(0,ncol=n)
for(i in 1:n){
    capm <- lm(retornos[,i]~r.indice)
    coef <- capm[["coefficients"]]
    betas[,i] <- coef[2]
}

datos <- round(rbind(rx,sigma,betas),4)
datos # Omit betas FX

## ----------------------------------------------------------------

n <- 3
diagbetas <- diag(c(betas[1:3]),n,n)
B <- matrix(c(diagbetas,
              0,1,0,
              0,0,1),nrow=5,byrow=TRUE)
w <- c(0.4,0.3,0.3)
(Bw <- B%*%w)

(omega <- cov)
(rsp <- t(Bw)%*%rx)
(sigmasp <- sqrt(t(Bw)%*%omega%*%Bw))
round(omega,4)
## ---------------------------------------------------------
## Calculos Retorno esperado y Varianza de Markowitz

mu <- rx[1:3]
cov <- omega[1:3,1:3]
var <- diag(cov)
sigma <- sqrt(var)

(rp <- t(w)%*%mu)
(sigmap <- sqrt(t(w)%*%cov%*%w))

## -----------------------------------------------------

# Optimizacion convexa: CVXR

library(CVXR)
# PMVG de Markowitz
w <- Variable(n)
problem <- Problem(Minimize(quad_form(w,cov)), 
                   constraints = list(w>=0, sum(w)==1))
result <- solve(problem)
(wpmvgCVX <- round(result$getValue(w),4))
rownames(wpmvgCVX) <- activos[1:3]

# ---------------------------------------------------------
# FOrmulacion del portafolio internacional
# Portafolio de minima varianza (sistematica)

w <- Variable(n)
Bwcvx <- B%*%w
problem <- Problem(Minimize(quad_form(Bwcvx,omega)), 
                   constraints = list(w>=0, sum(w)==1))
result <- solve(problem)
wpmvsCVX <- round(result$getValue(w),4)
rownames(wpmvsCVX) <- activos[1:3]

## Comparacion pesos
pesos <- cbind(wpmvgCVX,wpmvsCVX)

windows()
barplot(t(pesos),beside=TRUE,ylim=c(0,0.7),cex.names = 0.85,ylab="Part. (%)")
legend("topright",legend = c("PMVG","PIMV"),
       fill=c("black","gray"),cex=0.8)

## ------------------------------------------------------
## ------------------------------------------------------

## Efecto correlaciones

library(zoo)

# Informacion para los mercados locales
fechai <- '1994-12-01'
fechaf <- '2021-12-31'
periodicidad <- "monthly"   
activos <- c('^FTSE','^GSPC','^GDAXI')

precios.mkhist <- precios(activos,fechai,fechaf,periodicidad)
retornos.hist <- diff(log(precios.mkhist))[-1,]

corrSPFTSE <- na.omit(rollapply(retornos.hist, width=60, function(x) cor(x[,1],x[,2]), by.column=FALSE))[-1]
corrSPDAX <- na.omit(rollapply(retornos.hist, width=60, function(x) cor(x[,2],x[,3]), by.column=FALSE))[-1]
corrSPFTSE <- ts(corrSPFTSE,start=c(2000), frequency=12)
corrSPDAX <- ts(corrSPDAX,start=c(2000), frequency=12)

windows()
plot(corrSPFTSE,ylim=c(0.4,1),type="l",bty="L",ylab="Correlaciones")
lines(corrSPDAX,lty=2)
legend("topright",legend = c("SP&500-FTSE","SP&500-DAX"),
       lty=c(1,2),cex=0.9)



