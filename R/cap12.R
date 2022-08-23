#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestioón Moderna de portafolio
# Capitulo 12: Modelo Black-Litterman
# Copyright 2022
#-----------------------------------------------------------
#-----------------------------------------------------------

rm(list=ls())
library(zoo)
library(xts)
library(quantmod)
library(quadprog)

## Acciones seleccionadas
activos <- c("ADBE","MCD","MSCI","MSFT","NEE","PG","RSG","WMT")
fechai <- '2015-12-01'
fechaf <- '2020-12-31'
periodicidad <- "monthly"   
precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,]

mu <- colMeans(retornos)
cov <- cov(retornos)
var <- diag(cov)
sigma <- sqrt(var)
corr <- cor(retornos)
rf <- 0 
n <- length(mu)

# Indice
indice <- c("^GSPC")
p.indice <- precios(indice,fechai,fechaf,periodicidad)
r.indice <- diff(log(p.indice))[-1,]

mu.indice <- mean(r.indice)
sigma.indice <- sd(r.indice)

# Estimacion de parametros y anualizacion
betas <- matrix(0,ncol=n)
for(i in 1:n){
    capm <- lm(retornos[,i]~r.indice)
    coef <- capm[["coefficients"]]
    betas[,i] <- coef[2]
}

# Retornos de equilibrio
pi <- c(betas*(mu.indice*12-rf))
(datos <- round(rbind(mu*12,sigma*sqrt(12),betas,pi),4))
cov <- cov*12# Covarianzas anualizadas

delta <- (mu.indice-rf)/sigma.indice^2 # E(Rp) -rf / sigmap^2
wpequi <- solve(c(delta)*cov)%*%pi/sum(solve(c(delta)*cov)%*%pi)

windows()
barplot(t(wpequi))


##-----------------------------------------------------------
## ----------------------------------------------------------------

## Estimaciones de Views usando Fama-French

## Importar la info descargada desde la web-page
# library(readr)
# ff3 <- read_csv(".../ff3.csv")

## Acción: MSCI
datos <- cbind(ff3,retornos[,3])

## Modelo de mercado
ri <- datos[,5]-datos[,4]
rm <- datos[,1]
smb <- datos[,2]
hml <- datos[,3]

modelo1 <- lm(ri~rm+smb+hml)
summary(modelo1)

# Retornos estimados 
# FF3
(betam <- modelo1[["coefficients"]][2])
(betahml <- modelo1[["coefficients"]][4])
(re.msciff3 <- round((betam*mean(rm)+betahml*mean(hml))*12,4))

## Acción: NEE
datos <- cbind(ff3,retornos[,5])

## Modelo de mercado
ri <- datos[,5]-datos[,4]
rm <- datos[,1]
smb <- datos[,2]
hml <- datos[,3]

modelo2 <- lm(ri~rm+smb+hml)
summary(modelo2)

# Retornos estimados 
# FF3
(alpha <- modelo2[["coefficients"]][1])
(betam <- modelo2[["coefficients"]][2])
(betasmb <- modelo2[["coefficients"]][3])
(re.neeff3 <- round((alpha+betam*mean(rm)+betasmb*mean(smb))*12,4))

## Acción: RSG
datos <- cbind(ff3,retornos[,7])

## Modelo de mercado
ri <- datos[,5]-datos[,4]
rm <- datos[,1]
smb <- datos[,2]
hml <- datos[,3]

modelo3 <- lm(ri~rm+smb+hml)
summary(modelo3)

# Retornos estimados 
# FF3
(betam <- modelo3[["coefficients"]][2])
(betasmb <- modelo3[["coefficients"]][3])
(re.rsgff3 <- round((betam*mean(rm)+betasmb*mean(smb))*12,4))


## -----------------------------------------------------------
## -----------------------------------------------------------

## Implementacion del Modelo BL
## Parametros

tau <- 0.025 # NIvel de confiabilidad (Ej: al 95% es 0.025)

# Expectativas: views

q <- c(0.035,0.148,0.11)  # MCD-ADBE: 0.035; NEE: 0.148; MSCI-RSG: 0.13 
P <- rbind(c(-1,1,0,0,0,0,0,0),
           c(0,0,0,0,1,0,0,0),
           c(0,0,1,0,0,0,-1,0)) 
k <- dim(P)[1] # No. de views
omega <- matrix(0,k,k)
for(i in 1:k){
    omega[i,i] <- tau*(t(P[i,])%*%cov%*%P[i,])  # tau * P' V P
}

# Formula BL: retornos actualizados
mu.BL <- solve(solve(c(tau)*cov)+t(P)%*%solve(omega)%*%P)%*%(solve(c(tau)*cov)%*%pi+t(P)%*%solve(omega)%*%q)

Er <- t(rbind(t(mu.BL),pi))

windows()
barplot(t(Er),beside=TRUE,main="Grafico de retornos",ylim=c(0,0.18),
        ylab="Retornos esperados")
legend("topright",c("BL","Equilibrio"),fill=c("black","gray"))

# Vector de pesos
w.BL <- solve(c(delta)*cov)%*%mu.BL/sum(solve(c(delta)*cov)%*%mu.BL)

# Comparación de pesos para los modelos BL y MV
windows()
barplot(t(cbind(w.BL,wpequi)),beside=TRUE,main="Grafico de pesos",
        ylim=c(-0.2,0.5))
legend("bottomright",c("BL","Equilibrio"),fill=c("black","gray"))


## ----------------------------------------------------------------
## ----------------------------------------------------------------

# Evaluacion de desempeno 

# Dentro de muestra
# Portafolio BL
# Retorno 
(rp.BL <- t(w.BL)%*%mu.BL)
(rp.Equi <- t(wpequi)%*%pi)
(mu.indice*12)

# Riesgo
(sigma.BL <- sqrt(t(w.BL)%*%cov%*%w.BL))
(sigma.Equi <- sqrt(t(wpequi)%*%cov%*%wpequi))
(sigma.indice*sqrt(12))

# Coef. Sharpe
(S.BL <- (rp.BL-rf)/sigma.BL)
(S.Equi <- (rp.Equi-rf)/sigma.Equi)
(S.Indice <- (mu.indice*12-rf)/(sigma.indice*sqrt(12)))

# Desempeño historico - Valor del portafolio
valor <- 100 
t <- nrow(retornos)
n <- ncol(retornos)

# Retornos historicos

# Portafolio BL
rpbl <- retornos%*%w.BL
# Portafolio Equilibrio
rpequi <- retornos%*%wpequi

# Valor del portafolio BL
port.bl <- matrix(0, nrow=t)
port.bl[1] <- valor
for(i in 2:t){
    port.bl[i] <- port.bl[i-1]*exp(rpbl[i-1])
}

# Valor del portafolio Equilibrio
port.Eq <- matrix(0, nrow=t)
port.Eq[1] <- valor
for(i in 2:t){
    port.Eq[i] <- port.Eq[i-1]*exp(rpequi[i-1])
}

# Valor del benchmark
v.benchmark <- matrix(0, nrow=t)
v.benchmark[1] <- valor
for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.indice[i-1])
}

Performance <- cbind(port.bl,port.Eq,v.benchmark)
Performance <- ts(Performance,start=2016, frequency=12)
colnames(Performance) <- c("BL","PEq","Benchmark")
                     
windows()
plot(Performance[,1], type='l',main="Evaluación de desempeño", 
     ylab="Valor del portafolio", xlab="Tiempo (años)",bty="L")
lines(Performance[,2], lty=2)
lines(Performance[,3], col='darkgray')
legend("topleft",c("BL","PEq","S&P 500"),
       lty =c(1,2,1), col=c("black","black","darkgray"))

## ----------------------------------------------------------------

# Desempeno Out-sample

fechaifm <- '2020-12-01'
fechaffm <- '2021-12-31'

preciosfm <- precios(activos,fechaifm,fechaffm,periodicidad)
retornosfm <- diff(log(preciosfm))[-1,]
indicefm <- precios(indice,fechaifm,fechaffm,periodicidad)
r.indicefm <- diff(log(indicefm))[-1,]
mu.indicefm <- mean(r.indicefm)
sigma.indicefm <- sd(r.indicefm)

valor <- 100 
t <- nrow(retornosfm)+1
n <- ncol(retornosfm)

# Portafolio BL
rpblfm <- retornosfm%*%w.BL
(mu.rpbl <- mean(rpblfm)*12)
(sigma.pbl <- sd(rpblfm)*sqrt(12))
(S.bl <- (mu.rpbl-rf)/sigma.pbl)

# Valor del portafolio BL
port.blfm <- matrix(0, nrow=t)
port.blfm[1] <- valor
for(i in 2:t){
    port.blfm[i] <- port.blfm[i-1]*exp(rpblfm[i-1])
}

# Portafolio de equilibrio
rpEqfm <- retornosfm%*%wpequi
(mu.rpEq <- mean(rpEqfm)*12)
(sigma.pEq <- sd(rpEqfm)*sqrt(12))
(S.Eqfm <- (mu.rpEq-rf)/sigma.pEq)

# Valor del portafolio
port.Eqfm <- matrix(0, nrow=t)
port.Eqfm[1] <- valor
for(i in 2:t){
    port.Eqfm[i] <- port.Eqfm[i-1]*exp(rpEqfm[i-1])
}

# Benchmark
(S.indicefm <- (mu.indicefm*12-rf)/(sigma.indicefm*sqrt(12)))

v.benchmarkfm <- matrix(0, nrow=t)
v.benchmarkfm[1] <- valor
for(i in 2:t){
    v.benchmarkfm[i] <- v.benchmarkfm[i-1]*exp(r.indicefm[i-1])
}

Performancefm <- cbind(port.blfm,port.Eqfm,v.benchmarkfm)
Performancefm <- ts(Performancefm,start=c(2021,1), frequency=12)
colnames(Performancefm) <- c("BL","PEq","Benchmark")

windows()
plot(Performancefm[,1], type='l',main="Evaluacion de desempeno", bty="L",
     ylab="Valor del portafolio", xlab="Tiempo (meses)", ylim=c(min(Performancefm), max(Performancefm)))
lines(Performancefm[,2], lty=2)
lines(Performancefm[,3], col='darkgray')
legend("topleft",c("BL","PEq","S&P 500"),
       lty =c(1,2,1), col=c("black","black","darkgray"))


#----------------------------------------------------------------------
#----------------------------------------------------------------------

## Optimizacion sin cortos
Dmat <- cov*2
dvec <- rep(0,n)
Amat <- matrix(c(r.BL,-r.BL,rep(1,n),-rep(1,n),diag(length(r.BL))),n,n+4)
nport <- 1000
j <- seq(min(r.BL)+0.01,max(r.BL)-0.01,length=nport)
sigmapBL <- matrix(0,nrow=nport)
wpoBL <- matrix(0,nrow=nport,ncol=n)

for(i in 1:nport){
    rp <- j[i]
    bvec <- c(rp,-rp,1,-1,rep(0,n))
    poptimo <- solve.QP(Dmat,dvec,Amat,bvec,meq=2)
    wpoBL[i,] <- poptimo$solution
    sigmapBL[i,] <- sqrt(poptimo$value)
}

# Portafolio tangente del modelo BL

sharpe_port <- (j-rf)/sigmapBL
sharpe <- cbind(sharpe_port,wpoBL)
sharpe.sort <- sharpe[order(-sharpe[,1]),]
sharpe.sel <- cbind(sharpe.sort[1,])
wptBL <- round(cbind(sharpe.sel[2:length(sharpe.sel)]),6)
rownames(wptBL) <- c(activos)

windows()
barplot(t(wptBL),beside=TRUE,ylab= "Part. (%)")

# Comparación de pesos para el modelo BL
windows()
barplot(t(cbind(w.BL,wptBL)),beside=TRUE,main="Grafico de pesos")
legend("topright",c("BL con cortos","BL sin cortos"),fill=c("black","gray"))

# FE del modelo BL
rpSBL <- (t(wptBL)%*%mu.BL)*0.99
sigmaSBL <- sqrt(t(wptBL)%*%cov%*%wptBL)

windows()
plot(sigmapBL,rpBL,type="l",main="Plano Riesgo-Retorno",
     xlim=c(0.08,max(sigmapBL)),
     ylim=c(0.06,max(rpBL)),
     xlab="Riesgo",ylab="Retorno esperado")
points(sigmaSBL,rpSBL,pch=19)
text(sigmaSBL,rpSBL,labels="T",cex=0.8,pos=2)
