#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestión Moderna de Portafolio
# Capitulo 8: Modelos Multifactoriales
# Copyright 2022
#-----------------------------------------------------------
#-----------------------------------------------------------

# Ejemplo 8.1
rm(list=ls())
library(zoo)
library(xts)
library(quantmod)

## ----------------------------------------------------------------

# Información de WMT
fechai <- '2015-12-01'
fechaf <- '2020-12-31'
periodicidad <- "monthly"    
precios.wmt <- precios("WMT",fechai,fechaf,periodicidad)
r.wmt <- diff(log(precios.wmt))[-1,]

indice <- c("^GSPC")
indice.hist <- precios(indice,fechai,fechaf,periodicidad)
r.indice <- diff(log(indice.hist))[-1,]
(rm = mean(r.indice))
(sigmam <- sd(r.indice))

# Tasa de interes T-Bills
tnx <- c("^TNX")
tnx.hist <- precios(tnx,fechai,fechaf,periodicidad)/100
r.tnx <- diff(log(tnx.hist))[-1,]
mu.tnx = mean(r.tnx)
sigma.tnx <- sd(r.tnx)

windows()
plot(tnx.hist,type="l")

## Modelo de mercado
modelo1 <- lm(r.wmt~r.indice)
summary(modelo1)
(beta.wmt <- modelo1[["coefficients"]][2])
(r2adj1 <- summary(modelo1)$adj.r.squared)

## Modelo Factorial
modelo2 <- lm(r.wmt~r.indice+r.tnx)
summary(modelo2)

# Retornos estimados 
# CAPM
(betam <- modelo1[["coefficients"]][2])
(re.wmtcapm <- round((betam*rm),4))*12

# Factorial
(betam <- modelo2[["coefficients"]][2])
(betatnx <- modelo2[["coefficients"]][3])
(re.wmtfact <- round(betam*rm+betatnx*mu.tnx,4))*12

# --------------------------------------------------------------
## Graficos
ri <- matrix(r.wmt)
rmk <- matrix(r.indice)
rf <- matrix(r.tnx)

windows()
par(mfrow=c(1,2))
plot(rmk,ri,xlab="Índice", ylab="WMT",pch=19, ylim=c(-0.1,0.12))
abline(lm(ri~rmk), col="black")
plot(rf,ri,xlab="TNX", ylab="WMT",pch=19, ylim=c(-0.1,0.12))
abline(lm(ri~rf), col="black")


## ----------------------------------------------------------------

## Ejemplo 8.2

# Informacion
activos <- c("AAPL","AMZN","GOOG","MSFT")
fechai <- '2015-12-01'
fechaf <- '2020-12-31'
periodicidad <- "monthly"    
precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,]

indice <- c("^GSPC")
indice.hist <- precios(indice,fechai,fechaf,periodicidad)
r.indice <- diff(log(indice.hist))[-1,]

# Factores: INDIPRO, UNEMP
# library(readr)
# factores <- read_csv(".../factores.csv")
# factores <- as.matrix(factores[,2:3])
# ----------------------------------------------
# Se crea una sola matriz con los tres factores
# factores <- cbind(r.indice,factores)
# ----------------------------------------------

ri <- as.matrix(retornos[,'AAPL'])
rmk <- as.matrix(factores[,1])
rprod <- as.matrix(factores[,2])
runemp <- as.matrix(factores[,3])

windows()
par(mfrow=c(1,3))
plot(rmk,ri,xlab="Índice", ylab="AAPL",pch=19, ylim=c(-0.1,0.12))
abline(lm(ri~rmk), col="black")
plot(rprod,ri,xlab="INDIPRO", ylab="WMT",pch=19) #, ylim=c(-0.1,0.12)
abline(lm(ri~rprod), col="black")
plot(runemp,ri,xlab="UNEMP", ylab="WMT",pch=19)
abline(lm(ri~runemp), col="black")

# Regresiones factoriales
modeloAAPL <- lm(retornos[,1]~rmk+rprod+runemp)
summary(modeloAAPL)

modeloAMZN <- lm(retornos[,2]~rmk+rprod+runemp)
summary(modeloAMZN)

modeloGOOG <- lm(retornos[,3]~rmk+rprod+runemp)
summary(modeloGOOG)

modeloMSFT <- lm(retornos[,4]~rmk+rprod+runemp)
summary(modeloMSFT)

# Matriz de betas
(matriz.beta<- round(t(lm(retornos~rmk+rprod+runemp)$coef[-1,]),4))

# Covarianza factores
(Sigma.F <- round(cov(cbind(rmk, rprod, runemp)),6))

# Matriz Sigma errores
Sigma.error <-diag(c(summary(modeloAAPL)$sigma,
                     + summary(modeloAMZN)$sigma, 
                     + summary(modeloGOOG)$sigma,
                     + summary(modeloMSFT)$sigma)^2)
round(Sigma.error,4)

# Covarianza de los activos
(Sigma <- round(matriz.beta%*%Sigma.F%*%t(matriz.beta) + Sigma.error,4))
(round(cov(retornos),4))


# -----------------------------------------------------------------------

## Ejemplo 8.3
## Estimaciones usando la información de Fama-French

# Información de AAPL
fechai <- '2015-12-01'
fechaf <- '2020-12-31'
periodicidad <- "monthly" 
precios.aapl <- precios("AAPL",fechai,fechaf,periodicidad)
r.aapl <- diff(log(precios.aapl))[-1,]

# Importar la info descargada desde la web-page de Keneth French
# library(readr)
# datosff3 <- read_csv(".../datosff3.csv")
# ----------------------------------------------
# Se crea una sola matriz con todos factores y el
# activo
# datos <- cbind(datosff3[,2:5],r.aapl)
# ----------------------------------------------

## Modelo de mercado
ri <- datos[,5]-datos[,4]
rm <- datos[,1]

modelo1 <- lm(ri~rm)
summary(modelo1)
(betam <- modelo1[["coefficients"]][2])
(r2adj1 <- summary(modelo1)$adj.r.squared)

## Modelo FF3
ri 
rm 
smb <- datos[,2]
hml <- datos[,3]

modelo2 <- lm(ri~rm+smb+hml)
summary(modelo2)

## Modelo FF5
# datosff5 <- read_csv(".../datosff5.csv")
# ----------------------------------------------
# Se crea una sola matriz con todos factores y el
# datos <- cbind(datosff5[,2:7],r.aapl)
# ----------------------------------------------
ri 
rm 
smb 
hml 
rmw <- datos[,4]
cma <- datos[,5]

modelo3 <- lm(ri~rm+smb+hml+rmw+cma)
summary(modelo3)

# Retornos estimados 
# CAPM
(betam <- modelo1[["coefficients"]][2])
(re.aaplcapm <- round((betam*mean(rm)),4))

# FF3
(betam <- modelo2[["coefficients"]][2])
(betahml <- modelo2[["coefficients"]][4])
(re.aaplff3 <- round((betam*mean(rm)+betahml*mean(hml)),4))

# FF5
(betam <- modelo3[["coefficients"]][2])
(betarmw <- modelo3[["coefficients"]][5])
(re.aaplff5 <- round((betam*mean(rm)+betarmw*mean(rmw)),4))

# --------------------------------------------------------------
## Graficos

windows()
par(mfrow=c(3,2))
plot(rm,ri,xlab="Índice", ylab="AAPL",pch=19)
abline(lm(ri~rm), col="black")
plot(smb,ri,xlab="SMB", ylab="AAPL",pch=19)
abline(lm(ri~smb), col="black")
plot(hml,ri,xlab="HML", ylab="AAPL",pch=19)
abline(lm(ri~hml), col="black")
plot(rmw,ri,xlab="RMW", ylab="AAPL",pch=19)
abline(lm(ri~rmw), col="black")
plot(cma,ri,xlab="CMA", ylab="AAPL",pch=19)
abline(lm(ri~cma), col="black")

## ----------------------------------------------------------------





