#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestioón Moderna de portafolio
# Capitulo 4: Utilidad esperada
# Copyright 2022
#-----------------------------------------------------------
#-----------------------------------------------------------
#
# Notas:
# 1. Requiere cargar primero la funcion para importar precios  
#    de Yahoo Finance: "func_precios"
#-----------------------------------------------------------
#-----------------------------------------------------------

rm(list=ls())

library(zoo)
library(xts)
library(quantmod)
library(quadprog)

# Ejemplos 4.1 y 4.2

## Acciones seleccionadas
activos <- c("AAPL","AMZN")
fechai <- '2009-12-01'
fechaf <- '2021-12-31'
periodicidad <- "monthly"   

precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,]

(mu <- colMeans(retornos)*12)
(cov <- cov(retornos)*12)
var <- diag(cov)
(sigma <- sqrt(var))
(corr <- cor(retornos))
rf <- 0 
n <- length(mu)

# Portafolio tangente
rf <- 0.0
r <- mu -rf 
Z <- solve(cov,r)
sumaZ <- sum(Z)
(wpt <- Z/sumaZ)
rpt <- t(wpt)%*%mu
sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)
(delta <- (rpt -rf)/sigmapt^2)

# Portafolio óptimo
(wpMVu <- solve(c(delta)*cov)%*%r/sum(solve(c(delta)*cov)%*%r))

windows()
barplot(t(wpMVu))

##-----------------------------------------------------------
## ----------------------------------------------------------------

