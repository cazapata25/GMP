

#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestion Moderna de portafolio
# Funciones auxiliares
# Copyright 2022
#-----------------------------------------------------------
#-----------------------------------------------------------

## ---------------------------------------------------------
## Funcion para importar precios de Yahoo Finance

# Librerias necesarias
rm(list=ls())
if (!require("quantmod")) install.packages("quantmod"); library("quantmod")
if (!require("xts")) install.packages("xts"); library("xts")
if (!require("zoo")) install.packages("zoo"); library("zoo")

# O si ya estan instaladas, solo llamarlas como sigue:
library(quantmod)
library(xts)
library(zoo)

# Funcion: 
# Obtiene las series de precios de cierre ajustados paras los activos 
# indicados y segun el periodo de tiempo (fechas inicio y fin) y 
# periodicidad: dias, semanas, meses, años, etc...

precios <- function(activos,fechai,fechaf,periodicidad){
    precios <- xts()
    for(i in 1:length(activos)){
        tmp <- Ad(getSymbols(activos[i],from=fechai,to=fechaf,periodicity=periodicidad,auto.assign=FALSE))
        precios <- cbind(precios,tmp)
        precios <- na.approx(precios,na.rm=FALSE)
    }
    colnames(precios) <- activos
    tclass(precios) <- "Date"
    return(precios)
}

## Ej: Importar precios para las acciones de AAPL y AMZN
## Requiere definir las fechas de inicio, fin y periodicidad

fechai <- "2009-12-01"
fechaf <- "2021-12-31"
periodicidad <- "monthly"
activos <- c("AAPL","AMZN")

precios.hist <- precios(activos,fechai,fechaf,periodicidad)
# Resultado: genera la matriz de precios: "precios.hist"
