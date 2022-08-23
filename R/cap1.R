
#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestion Moderna de portafolio
# Capitulo 1: Medidas de retorno y riesgo
# Copyright 2022
#-----------------------------------------------------------
#-----------------------------------------------------------
#
# Notas:
# 1. Requiere cargar primero la funcion para importar precios  
#    de Yahoo Finance: "func_precios"
#-----------------------------------------------------------
#-----------------------------------------------------------

## Ejemplos 1.2, 1.3 y 1.4

rm(list=ls())
library(quantmod)
library(xts)
library(zoo)

# Para importar la informacion de precios para las acciones AAPL y AMZN
# se definen los parametros de la función: fechas y periodicidad

fechai <- "2009-12-01"
fechaf <- "2021-12-31"
periodicidad <- "monthly"
activos <- c("AAPL","AMZN")

# Con la funcion precios se obtiene la info historica
precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,] 
#preciosc <- as.zoo(precios[2:145,])

# Se estiman retornos esperados y riesgo (covarianzas, varianzas, volatilidades)
# Ejemplo 1.2: Retornos anualizados
(mu <- colMeans(retornos)*12)

# Medidas de riesgo: covarianzas, varianzaas y volatilidades - anualizadas
# Ejemplo 1.3
cov <- cov(retornos)*12
var <- diag(cov)
(sigma <- sqrt(var))
(corr <- cor(retornos))

n <- length(mu)

# Portafolio inicial: Ejemplo 1.4
# Se construye el portafolio inicial: w1=w2=0.5
wi <- rep(0.5,n)
(rp <- t(wi)%*%mu)
(sigmap2 <- t(wi)%*%cov%*%wi)
(sigmap <- sqrt(sigmap2))

# Grafico de precios
windows()
par(mfrow=c(1,2),bty="L")
plot(precios.hist$AAPL,type="l",xlab="Tiempo",ylab="Precio",main="AAPL",bty="L",lwd=2)
plot(precios.hist$AMZN,type="l",xlab="Tiempo",ylab="Precio",main="AMZN",bty="L",lwd=2)

# Histogramas de frecuencias
windows()
par(mfrow=c(1,2))
hist(retornos$AAPL,breaks=20,main="AAPL",probability=TRUE,
     xlim=c(-0.3,0.3),yaxt='n',xlab="",ylab="")
lines(density(retornos$AAPL),lwd=2)
abline(v=mean(retornos$AAPL),lty=2,lwd=2)
hist(retornos$AMZN,breaks=20,main="AMZN",probability = TRUE,
     xlim=c(-0.3,0.3),yaxt='n',xlab="",ylab="")
lines(density(retornos$AMZN),lwd=2)
abline(v=mean(retornos$AMZN),lty=2,lwd=2)

