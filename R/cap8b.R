#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestión Moderna de Portafolio
# Capitulo 9: Modelos Multifactoriales - estrategias de 
#             inversión
# Copyright 2022
#-----------------------------------------------------------
#-----------------------------------------------------------

# Ejemplo 9.1
rm(list=ls())
library(zoo)
library(xts)
library(quantmod)
library(quadprog)

## ----------------------------------------------------------------

## Ejemplo 9.1: Portafolio factorial

fechai <- '2015-12-01'
fechaf <- '2020-12-31'
periodicidad <- "monthly" 
activos <- c("AAPL","AMZN","MSFT","GOOG")
precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,]
(mu <- colMeans(retornos))

# Factores del modelo de estimacion
# Importar la info descargada desde la web-page de Keneth French
# library(readr)
# datosff3 <- read_csv(".../datosff3.csv")
# ----------------------------------------------
# Se crea una sola matriz con todos factores y el activo
# datos <- cbind(datosff3[,2:5],retornos)
# ----------------------------------------------

# Variables del modelo
rft <- datos[,4]
# Matriz de exesos de retornos de los activos
r.act <- datos[,5:8]-rft
# Factores
rmk <- datos[,1]
smb <- datos[,2]
hml <- datos[,3]

# Matriz de betas
(matriz.beta<- round(t(lm(retornos~rmk+smb+hml)$coef[-1,]),4))

# Covarianza factores
(Sigma.F <- round(cov(cbind(rmk, smb, hml)),6))

# Regresiones factoriales
modeloAAPL <- lm(retornos[,1]~rmk+smb+hml)
modeloAMZN <- lm(retornos[,2]~rmk+smb+hml)
modeloGOOG <- lm(retornos[,3]~rmk+smb+hml)
modeloMSFT <- lm(retornos[,4]~rmk+smb+hml)

# Matriz Sigma errores
Sigma.error <-diag(c(summary(modeloAAPL)$sigma,
                     + summary(modeloAMZN)$sigma, 
                     + summary(modeloGOOG)$sigma,
                     + summary(modeloMSFT)$sigma)^2)
round(Sigma.error,4)

# Covarianza de los activos
(Sigma <- round(matriz.beta%*%Sigma.F%*%t(matriz.beta) + Sigma.error,4))

# Portafolio optimo de minima varianza
n <- length(activos)
Sol <-solve.QP(Dmat=2*Sigma, dvec=rep(0, n), 
               Amat=cbind(rep(1,n), diag(n)), bvec=c(1, rep(0, n)),
               meq=1)$solution

wpmvf <- round(Sol,4)
names(wpmvf) <- activos

windows()
barplot(t(wpmvf),beside=TRUE,ylim=c(0,0.8),cex.names = 0.85,ylab="Part. (%)")


# Portafolio ocn retorno objetivo: 
rep <- mean(mu)
Sol2 <-solve.QP(Dmat=2*Sigma, dvec=rep(0, n), 
               Amat=cbind(mu,rep(1,n), diag(n)), bvec=c(rep,1, rep(0, n)),
               meq=2)$solution
wpmvf2 <- round(Sol2,4)
names(wpmvf2) <- activos

windows()
barplot(t(wpmvf2),beside=TRUE,ylim=c(0,0.8),cex.names = 0.85,ylab="Part. (%)")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Estrategia factorial usando Momentum

library(TTR)
library(FinancialInstrument)
library(PerformanceAnalytics)

symbols=c("AAPL","FE",'PBR',"BAC", "GE")
getSymbols(symbols, from="2016-01-02", to="2022-07-01", 
           periodicity = "monthly") 
for (symbol in symbols){
    x<-get(symbol)
    x<-to.monthly(x,indexAt='lastof',drop.time = TRUE)
    tformat(x)<-'%Y-%m-%d'
    colnames(x)<-gsub("x",symbol,colnames(x))
    x<-x[,6]
    assign(symbol,x)
}

ajustados<-do.call(merge,lapply(symbols,get))
#hago un ROC de 3 periodos
roc<-ROC(ajustados,n=3,type= "discrete")
#aplico la funcion Rank al ROC ( 1 el mas alto ROC)
r<-as.xts(t(apply(-roc,1,rank)))

rankRB<-function(x){
    r<-as.xts(t(apply(-x,1,rank,na.last="keep")))
    return(r)
    
}
MonthlyAd<-function(x){
    sym<-sub("\\..*$", "",names(x)[1]) 
    Ad(to.monthly(x,indexAt='lastof',drop.time=TRUE,name=sym))
    
}
CAGR<-function(x,m){
    x<-na.omit(x)
    cagr<-apply(x,2,function(x,m) prod(1+x)^(1/(length(x) / m))-1,m=m)
    return(cagr)
}
simpleMomentumtest<-function(xts.ret,xts.rank,n=1,ret.fill.na=3){
    lag.rank<-lag(xts.rank, k=1, na.pad=TRUE)
    n2<-nrow(lag.rank[is.na(lag.rank[,1])==TRUE])
    z<-max(n2,ret.fill.na)
    lag.rank<-as.matrix(lag.rank)
    lag.rank[lag.rank>n]<-NA
    lag.rank[lag.rank<=n]<-1
    mat.ret<-as.matrix(xts.ret)*lag.rank
    vec.ret<-rowMeans(mat.ret, na.rm=TRUE)
    vec.ret[1:z]<-NA
    vec.ret<-xts(x=vec.ret,order.by = index(xts.ret))
    f<-list(mat=mat.ret,ret=vec.ret,rank=lag.rank)
    return(f)
}

monthly.returns<-ROC(x=ajustados, n=1,type="discrete",na.pad=TRUE)
roc.three<-ROC(x=ajustados, n=3, type="discrete")
rank.three<-rankRB(roc.three)

roc.six<-ROC(x=ajustados, n=6, type="discrete")
rank.six<-rankRB(roc.six)

roc.nine<-ROC(x=ajustados, n=9, type="discrete")
rank.nine<-rankRB(roc.nine)

roc.twelve<-ROC(x=ajustados, n=12, type="discrete")
rank.twelve<-rankRB(roc.twelve)

num.assets<-4

case1<-simpleMomentumtest(xts.ret=monthly.returns,xts.rank=rank.three,
                          n= num.assets, ret.fill.na= 15)

case2<-simpleMomentumtest(xts.ret=monthly.returns,xts.rank=rank.six,
                          n= num.assets, ret.fill.na= 15)

case3<-simpleMomentumtest(xts.ret=monthly.returns,xts.rank=rank.nine,
                          n= num.assets, ret.fill.na= 15)

case4<-simpleMomentumtest(xts.ret=monthly.returns,xts.rank=rank.twelve,
                          n= num.assets, ret.fill.na= 15)

retornos<-cbind(case1$ret,case2$ret, case3$ret, case4$ret)
colnames(retornos)<-c("3-Month","6-Month","9-Month","12-Month")

windows()
charts.PerformanceSummary(R=retornos,Rf=0, geometric=TRUE,
                          main= "Momentum Cululative Return: 4 assets")

table.Stats(retornos)
CAGR(retornos,m=12)
print("End")






