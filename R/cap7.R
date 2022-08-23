#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestion Moderna de portafolio
# Capitulo 6: Evalaucion de desempeño
# Copyright 2022
#-----------------------------------------------------------
#------------------------------------------------------------------
#
# Notas:
# 1. Requiere cargar primero la funcion para importar precios  
#    de Yahoo Finance: "func_precios"
#-----------------------------------------------------------
#-----------------------------------------------------------

## Informacion usada para todos los ejemplos

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
n <- length(mu)

# Ejemplo 7.1
# Tabla 7.1
(datos <- round(rbind(mu,sigma),4))

# Informacion del Indice
indice <- c("^GSPC")
p.indice <- precios(indice,fechai,fechaf,periodicidad)
r.indice <- diff(log(p.indice))[-1,]

(mu.indice <- mean(r.indice))
(sigma.indice <- sd(r.indice))

rf <- 0 
n <- length(mu)

# FE y PMVG de Markowitz
ones <- rep(1,n)
x <- t(mu)%*%solve(cov)%*%mu
y <- t(mu)%*%solve(cov)%*%ones
z <- t(ones)%*%solve(cov)%*%ones
d <- x*z - y*y
g <- (solve(cov,ones)%*%x-solve(cov,mu)%*%y)%*%solve(d)
h <- (solve(cov,mu)%*%z-solve(cov,ones)%*%y)%*%solve(d)

# FE
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

# PMVG
cov_inv_1 <- solve(cov, ones) 
(wpmvg <- (1/as.numeric(ones %*% cov_inv_1)) * cov_inv_1)
(rpmvg <- mu%*%wpmvg)
(sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg))

# Plano Riesgo-Retorno
windows()
plot(sigma,mu, ylim=c(min(mu)*0.9,max(mu)*1.15), 
     xlim=c(min(sigma)*0.5,max(sigma)*1.15),
     main="Plano riesgo-retorno",xlab= "Riesgo",
     ylab="Retorno esperado")
points(sigmapo,rpo, type="l") 
points(sigmapmvg,rpmvg, col="black",pch=19)
text(sigmapmvg,rpmvg,labels="PMVG",pos = 2, cex=0.7)
text(sigma,mu,labels = activos,pos = 2,cex=0.6) 

windows()
barplot(t(wpmvg), main="PMVG",axisnames = TRUE, beside = TRUE,ylab= "Part. (%)")

# ---------------------------------

# Evaluacion de desempeno in-sample

# Retorno activo: 
(Ramv = rpmvg - mu.indice)

# In-sample
valor <- 100 
t <- nrow(retornos)
n <- ncol(retornos)

# Retornos historicos
# PMVG
rhpmvg <- retornos%*%wpmvg

# Valor del portafolio
# PMV
port.mv <- matrix(0, nrow=t)
port.mv[1] <- valor
for(i in 2:t){
    port.mv[i] <- port.mv[i-1]*exp(rhpmvg[i-1])
}

# Benchmark
v.benchmark <- matrix(0, nrow=t)
v.benchmark[1] <- valor
for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.indice[i-1])
}

Performance <- cbind(port.mv,v.benchmark)
Performance <- ts(Performance,start=2016, frequency=12)
colnames(Performance) <- c("PMVG","Benchmark")
                     
windows()
plot(Performance[,1], type='l',main="Evaluacion de desempeno", 
     ylab="Valor del portafolio", xlab="Tiempo (anos)",bty="L")
lines(Performance[,2], col='darkgray')
legend("topleft",c("PMVG","S&P 500"),
       lty =c(1,1), col=c("black","darkgray"))

#----------------------------
# Desempeño out-sample

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

# Retornos 
# PMVG
rpmvfm <- retornosfm%*%wpmvg
(mu.rpmvg <- mean(rpmvfm))
(sigma.pmvg <- sd(rpmvfm))

# Valor del portafolio
port.mvfm <- matrix(0, nrow=t)
port.mvfm[1] <- valor
for(i in 2:t){
    port.mvfm[i] <- port.mvfm[i-1]*exp(rpmvfm[i-1])
}

# Benchmark
v.benchmarkfm <- matrix(0, nrow=t)
v.benchmarkfm[1] <- valor
for(i in 2:t){
    v.benchmarkfm[i] <- v.benchmarkfm[i-1]*exp(r.indicefm[i-1])
}

Performancefm <- cbind(port.mvfm,v.benchmarkfm)
Performancefm <- ts(Performancefm,start=c(2021,1), frequency=12)
colnames(Performancefm) <- c("PMVG","Benchmark")

windows()
plot(Performancefm[,1], type='l',main="Evaluacion de desempeno", bty="L",
     ylab="Valor del portafolio", xlab="Tiempo (meses)", ylim=c(min(Performancefm), max(Performancefm)))
lines(Performancefm[,2], col='darkgray')
legend("topleft",c("PMVG","S&P 500"),
       lty =c(1,1), col=c("black","darkgray"))


# ----------------------------------------------------------
# Optimizacion sin cortos

library(quadprog)

nport <- 1000
Dmat <- cov*2
Amat <- cbind(mu,rep(1,n),diag(1,nrow=n))
dvec <- rep(0,n)
j <- seq(min(mu)*1.001,max(mu)*0.999,length=nport)
wpo <- matrix(0,nrow=nport,ncol=n)
sigmapo <- matrix(0,nrow=nport)
rpo <- matrix(0,nrow=nport)

for(i in 1:nport){
    bvec <- c(j[i],1,rep(0,n)) # restricciones de igualdad del sistema
    result <- solve.QP(Dmat, dvec, Amat, bvec, meq=2) #meq=restriccones de igualdad
    wj <- result[["solution"]]
    wpo[i,] <- wj
    rpo[i] <- mu%*%wj
    sigmapo[i] <- sqrt(t(wj)%*%cov%*%wj)
}

# Plano Riesgo-Retorno
windows()
plot(sigma,mu, xlim=c(min(sigma)*0.5,max(sigma)*1.2), ylim=c(min(mu)*0.8,max(mu)*1.2),
     cex=0.6,xlab="Riesgo",ylab="Retorno",bty="L")
points(sigmapo,rpo,type='l')
points(sigmapmvg,rpmvg,pch=19)
points(sigmapt,rpt,pch=19)
text(sigma,mu,labels=activos,pos = 4, cex=0.7)
text(sigmapmvg,rpmvg,labels="PMVG",pos = 2, cex=0.7)

# ----------------------------------------------------------

## Ejemplo 7.2
## Medida de Sharpe

sharpe_port <- (rpo-rf)/sigmapo
tabla <- cbind(sharpe_port,rpo,sigmapo)
colnames(wpo) <- activos
tabla <- cbind(tabla,wpo)
sort_tabla <- tabla[order(-tabla[,1]),]
port_maxsharpe <- round(sort_tabla[1,],4)

(wpt <- port_maxsharpe[4:length(port_maxsharpe)])
(rpt <- port_maxsharpe[2])
(sigmapt <- port_maxsharpe[3])

windows()
barplot(t(wpt),axisnames=TRUE,beside = TRUE,ylab= "Part. (%)",
        ylim=c(0,0.4), main="PT")

# -------------------------------
# Evaluacion de desempeno in-sample

# Retorno activo: 
(Rasharpe = rpt - mu.indice)

# Coef. Sharpe:
(c.sharpe = (rpt-rf)/sigmapt)

# In-sample
valor <- 100 
t <- nrow(retornos)
n <- ncol(retornos)

# Retornos historicos
rpsharpe <- retornos%*%wpt

# Valor del portafolio
port.sharpe <- matrix(0, nrow=t)
port.sharpe[1] <- valor
for(i in 2:t){
    port.sharpe[i] <- port.sharpe[i-1]*exp(rpsharpe[i-1])
}

# Benchmark
v.benchmark <- matrix(0, nrow=t)
v.benchmark[1] <- valor
for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.indice[i-1])
}

Performance <- cbind(port.sharpe,v.benchmark)
Performance <- ts(Performance,start=2016, frequency=12)
colnames(Performance) <- c("Sharpe","Benchmark")

windows()
plot(Performance[,1], type='l',main="Evaluacion de desempeño", 
     ylab="Valor del portafolio", xlab="Tiempo (años)",bty="L")
lines(Performance[,2], col='darkgray')
legend("topleft",c("Sharpe","S&P 500"),
       lty =c(1,1), col=c("black","darkgray"))

#----------------------------
# Desempeno Out-sample

valor <- 100 
t <- nrow(retornosfm)+1
n <- ncol(retornosfm)

# Retornos del portafolio
rpsharpefm <- retornosfm%*%wpt
(mu.sharpefm <- mean(rpsharpefm))
(sigma.sharpefm <- sd(rpsharpefm))

# Valor del portafolio
port.sharpefm <- matrix(0, nrow=t)
port.sharpefm[1] <- valor
for(i in 2:t){
    port.sharpefm[i] <- port.sharpe[i-1]*exp(rpsharpefm[i-1])
}

# Benchmark
v.benchmarkfm <- matrix(0, nrow=t)
v.benchmarkfm[1] <- valor
for(i in 2:t){
    v.benchmarkfm[i] <- v.benchmarkfm[i-1]*exp(r.indicefm[i-1])
}

Performancefm <- cbind(port.sharpefm,v.benchmarkfm)
Performancefm <- ts(Performancefm,start=c(2021,1), frequency=12)
colnames(Performancefm) <- c("Sharpe","Benchmark")

windows()
plot(Performancefm[,1], type='l',main="Evaluacion de desempeño", bty="L",
     ylab="Valor del portafolio", xlab="Tiempo (meses)", ylim=c(min(Performancefm), max(Performancefm)))
lines(Performancefm[,2], col='darkgray')
legend("topleft",c("Sharpe","S&P 500"),
       lty =c(1,1), col=c("black","darkgray"))

# Retorno activo: 
(Rasharpe = mu.sharpefm - mu.indicefm)

# Coef. Sharpe:
(c.sharpefm = (mu.sharpefm -rf)/sigma.sharpefm)
(c.sharpeindicefm = (mu.indicefm -rf)/sigma.indicefm)


## ----------------------------------------------------------------
## ----------------------------------------------------------------

## Ejemplo 7.3
## Medida de Treynor

# Estimación de betas
varm <- sigma.indice^2
betas <- matrix(0,ncol=n)
varerror <- matrix(0,ncol=n)
for(i in 1:n){
    capm <- lm(retornos[,i]~r.indice)
    coef <- capm[["coefficients"]]
    betas[,i] <- coef[2]
    varerror[i] <- var(modelo[["residuals"]])
}

treynor <- (mu-rf)/betas

matrix <- t(rbind(treynor,t(mu),t(sigma),betas,varerror))
sort.matrix <- matrix[order(-matrix[,1]),]
ratio1 <- ((sort.matrix[,2]-rf)*sort.matrix[,4])/sort.matrix[,5]
ratio2 <- sort.matrix[,4]^2/sort.matrix[,5]
sumacu1 <- cumsum(ratio1)
sumacu2 <- cumsum(ratio2)
tasac <- (varm%*%sumacu1)/(1+varm%*%sumacu2)

diff <- sort.matrix[,1]-tasac
cond.diff <- diff[!is.na(diff)&diff>0]
n.optimo <- length(cond.diff)
cuttoff <- max(tasac)

Zi <- (sort.matrix[,4]/sort.matrix[,5])*(sort.matrix[,1]-cuttoff)
Zi <- pmax(Zi,0)
wpot <- round(Zi/sum(Zi),4)
(wpot <- cbind(wpot[activos]))
(rpot <- t(wpot)%*%mu)
(sigmapot <- sqrt(t(wpot)%*%cov%*%wpot))

(betaP <- betas%*%wpot)

windows()
barplot(t(wpot),axisnames=TRUE,beside = TRUE,ylab= "Part. (%)",
        ylim=c(0,0.4))

# Evaluacion de desempeño 

# Retorno activo: 
(Ratreyor = rpot - mu.indice)

# Coef. Treynor:
(c.treynor = (rpot-rf)/betaP)

# In-sample
valor <- 100 
t <- nrow(retornos)
n <- ncol(retornos)

# Retornos historicos
rptreynor <- retornos%*%wpot

# Valor del portafolio
port.treynor <- matrix(0, nrow=t)
port.treynor[1] <- valor
for(i in 2:t){
    port.treynor[i] <- port.treynor[i-1]*exp(rptreynor[i-1])
}

# Benchmark
v.benchmark <- matrix(0, nrow=t)
v.benchmark[1] <- valor
for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.indice[i-1])
}

Performance <- cbind(port.treynor,v.benchmark)
Performance <- ts(Performance,start=2016, frequency=12)
colnames(Performance) <- c("Treynor","Benchmark")

windows()
plot(Performance[,1], type='l',main="Evaluacion de desempeño", 
     ylab="Valor del portafolio", xlab="Tiempo (años)",bty="L")
lines(Performance[,2], col='darkgray')
legend("topleft",c("Treynor","S&P 500"),
       lty =c(1,1), col=c("black","darkgray"))

#----------------------------
# Desempeno Out-sample

valor <- 100 
t <- nrow(retornosfm)+1
n <- ncol(retornosfm)

# Retornos del portafolio
rptreynorfm <- retornosfm%*%wpot
(mu.treynorfm <- mean(rptreynorfm))
(sigma.treynorfm <- sd(rptreynorfm))

# Valor del portafolio
port.treynorfm <- matrix(0, nrow=t)
port.treynorfm[1] <- valor
for(i in 2:t){
    port.treynorfm[i] <- port.treynorfm[i-1]*exp(rptreynorfm[i-1])
}

# Benchmark
v.benchmarkfm <- matrix(0, nrow=t)
v.benchmarkfm[1] <- valor
for(i in 2:t){
    v.benchmarkfm[i] <- v.benchmarkfm[i-1]*exp(r.indicefm[i-1])
}

Performancefm <- cbind(port.treynorfm,v.benchmarkfm)
Performancefm <- ts(Performancefm,start=c(2021,1), frequency=12)
colnames(Performancefm) <- c("Treynor","Benchmark")

windows()
plot(Performancefm[,1], type='l',main="Evaluacion de desempeño", bty="L",
     ylab="Valor del portafolio", xlab="Tiempo (meses)", ylim=c(min(Performancefm), max(Performancefm)))
lines(Performancefm[,2], col='darkgray')
legend("topleft",c("Treynor","S&P 500"),
       lty =c(1,1), col=c("black","darkgray"))

# Retorno activo: 
(Ratreynorfm = mu.treynorfm - mu.indicefm)

# Coef. Treynor
(c.treynorfm = (mu.treynorfm - rf)/betaP)


## ----------------------------------------------------------------
## ----------------------------------------------------------------

## Ejemplo 7.4
## Medida de Sortino
h <- 0
semiret <- pmin(retornos,h)
semicov <- cov(semiret)
semivar <- diag(semicov)
semisigma <- sqrt(semivar)

# FE de Sortino sin cortos
nport <- 1000
Dmat <- semicov*2
Amat <- cbind(mu,rep(1,n),diag(1,nrow=n))
dvec <- rep(0,n)
j <- seq(min(mu)*1.001,max(mu)*0.999,length=nport)
wpoS <- matrix(0,nrow=nport,ncol=n)
sigmapoS <- matrix(0,nrow=nport)
rpoS <- matrix(0,nrow=nport)
for(i in 1:nport){
    bvec <- c(j[i],1,rep(0,n)) 
    result <- solve.QP(Dmat, dvec, Amat, bvec, meq=2) 
    wj <- result[["solution"]]
    wpoS[i,] <- wj
    rpoS[i] <- mu%*%wj
    sigmapoS[i] <- sqrt(t(wj)%*%semicov%*%wj)
}

#--------------------------------------------------------
# Portafolio Tangente de Sortino

sortino_port <- (rpoS-rf)/sigmapoS
tabla <- cbind(sortino_port,rpoS,sigmapoS)
colnames(tabla) <- c("Sortino","Retorno","Riesgo")
colnames(wpo) <- activos
tabla2 <- cbind(tabla,wpoS)
sort_tabla2 <- tabla2[order(-tabla2[,1]),]
port_maxsortino <- round(sort_tabla2[1,],4)
(wpts <- port_maxsortino[4:length(port_maxsortino)])
names(wpts) <- activos
(rpts <- t(wpts)%*%mu)
(sigmapts <- sqrt(t(wpts)%*%semicov%*%wpts))

windows()
barplot(t(wpts),axisnames=TRUE,beside = TRUE,ylab= "Part. (%)",
        ylim=c(0,0.35))

# Evaluacion de desempeño 

# Retorno activo: 
(RaSortino = rpts - mu.indice)

# Coef. Sortino:
(c.Sortino = (rpts-rf)/sigmapts)

# Semi-desiación estándar índice y coef. Sortino:
semiretindice <- pmin(r.indice,h)
(semisigmaindice <- sd(semiretindice))
(c.Sortinoindice = (mu.indice-rf)/semisigmaindice)

# In-sample
valor <- 100 
t <- nrow(retornos)
n <- ncol(retornos)

# Retornos historicos
rpsortino <- retornos%*%wpts

# Valor del portafolio
port.sortino <- matrix(0, nrow=t)
port.sortino[1] <- valor
for(i in 2:t){
    port.sortino[i] <- port.sortino[i-1]*exp(rpsortino[i-1])
}

# Benchmark
v.benchmark <- matrix(0, nrow=t)
v.benchmark[1] <- valor
for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.indice[i-1])
}

Performance <- cbind(port.sortino,v.benchmark)
Performance <- ts(Performance,start=2016, frequency=12)
colnames(Performance) <- c("Sortino","Benchmark")

windows()
plot(Performance[,1], type='l',main="Evaluacion de desempeño", 
     ylab="Valor del portafolio", xlab="Tiempo (años)",bty="L")
lines(Performance[,2], col='darkgray')
legend("topleft",c("Sortino","S&P 500"),
       lty =c(1,1), col=c("black","darkgray"))

#----------------------------
# Desempeno Out-sample

valor <- 100 
t <- nrow(retornosfm)+1
n <- ncol(retornosfm)

# Retornos del portafolio
rpsortinofm <- retornosfm%*%wpts
(mu.sortinofm <- mean(rpsortinofm))
(sigma.sortinofm <- sd(rpsortinofm))

# Valor del portafolio
port.sortinofm <- matrix(0, nrow=t)
port.sortinofm[1] <- valor
for(i in 2:t){
    port.sortinofm[i] <- port.sortinofm[i-1]*exp(rpsortinofm[i-1])
}

# Benchmark
v.benchmarkfm <- matrix(0, nrow=t)
v.benchmarkfm[1] <- valor
for(i in 2:t){
    v.benchmarkfm[i] <- v.benchmarkfm[i-1]*exp(r.indicefm[i-1])
}

Performancefm <- cbind(port.sortinofm,v.benchmarkfm)
Performancefm <- ts(Performancefm,start=c(2021,1), frequency=12)
colnames(Performancefm) <- c("Sortino","Benchmark")

windows()
plot(Performancefm[,1], type='l',main="Evaluacion de desempeño", bty="L",
     ylab="Valor del portafolio", xlab="Tiempo (meses)", ylim=c(min(Performancefm), max(Performancefm)))
lines(Performancefm[,2], col='darkgray')
legend("topleft",c("Sortino","S&P 500"),
       lty =c(1,1), col=c("black","darkgray"))

# Retorno activo: 
(Rasortinofm = mu.sortinofm - mu.indicefm)

# Semi-desviaciones estándar
# POrtafolio
semiretportfm <- pmin(rpsortinofm,h)
(semisigmasortinofm <- sd(semiretportfm))
(c.Sortinoindicefm = (mu.sortinofm-rf)/semisigmaindicefm)

# Indice
semiretindicefm <- pmin(r.indicefm,h)
(semisigmaindicefm <- sd(semiretindicefm))
(c.Sortinoindicefm = (mu.indicefm-rf)/semisigmaindicefm)


# Coef. Sortino:
(c.sortinofm = (mu.sortinofm -rf)/ssigma.sortinofm)
(c.sortinoindicefm = (mu.indicefm -rf)/semisigmaindicefm)


## ----------------------------------------------------------------
## ----------------------------------------------------------------

## Ejemplo 7.5
# Medida Omega

library(quadprog)
library(ROI)

h <- 0
omega <- function(retornos, h){
    s <- nrow(retornos)
    Amat <- rbind(cbind(as.matrix(retornos), diag(s), 0),# u_s >= h - r_s' y
                  c(rep(0, n), rep(1, s),  0), # sum(u) = 1
                  c(rep(1, n), rep(0, s), -1), # sum(y) = z
                  c(mu,        rep(0, s), - h),
                  c(rep(0, n), rep(0, s),  1)) # mu' y  >= k * z
    var.names <- c(paste0("y_omega_aux", seq_len(n)),
                   paste0("u_omega_aux", seq_len(s)), 
                   "z_omega")
    constraint <-  L_constraint(L = Amat,  
                                dir = c(rep(">=", s), "==", "==", ">=", ">"),
                                rhs = c(rep(h, s), 1, 0, 0, 1e-05), 
                                names=var.names)
    objective  <- L_objective(L = c(mu, rep(0, s), -h))
    list(objective = objective, constraint = constraint) 
}

budget_constraint <- function(retornos, dir = "==", rhs = 1) {
    x.names <- activos
    L_constraint(L = rep(1, n), 
                 dir = dir,  rhs = rhs, names = x.names)
}

tmp <- omega(retornos,h)
pomega  <- OP(objective  =  tmp$objective,
              constraints = tmp$constraint, maximum = T)
solom <- ROI_solve(pomega, solver = "glpk")
sol_omega <- solution(solom)
wpomega <- round(sol_omega[1:n]/sol_omega["z_omega"], 4)
names(wpomega) <- activos

(rpomega <- mu%*%wpomega)
(sigmapomega <- sqrt(t(wpomega)%*%cov%*%wpomega))

windows()
barplot(t(wpomega),axisnames=TRUE,beside = TRUE,ylab= "Part. (%)",
        ylim=c(0,0.4))

# Evaluacion de desempeño 

# Retorno activo: 
(RaOmega = rpomega - mu.indice)

# Medida Omega:
rpohist <- retornos%*%wpomega
(omegap <- sum(pmax(rpohist-h,0))/sum(pmax(h-rpohist,0)))
(riskomegap <- sum(pmax(h-rpohist,0)))

# Medida del Indice
(omegaindice <- sum(pmax(r.indice-h,0))/sum(pmax(h-r.indice,0)))
(riskoindice <- sum(pmax(h-r.indice,0)))

# In-sample
valor <- 100 
t <- nrow(retornos)
n <- ncol(retornos)

# Retornos historicos
rpomega <- rpohist

# Valor del portafolio
port.omega <- matrix(0, nrow=t)
port.omega[1] <- valor
for(i in 2:t){
    port.omega[i] <- port.omega[i-1]*exp(rpomega[i-1])
}

# Benchmark
v.benchmark <- matrix(0, nrow=t)
v.benchmark[1] <- valor
for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.indice[i-1])
}

Performance <- cbind(port.omega,v.benchmark)
Performance <- ts(Performance,start=2016, frequency=12)
colnames(Performance) <- c("Omega","Benchmark")

windows()
plot(Performance[,1], type='l',main="Evaluacion de desempeño", 
     ylab="Valor del portafolio", xlab="Tiempo (años)",bty="L")
lines(Performance[,2], col='darkgray')
legend("topleft",c("Omega","S&P 500"),
       lty =c(1,1), col=c("black","darkgray"))

#----------------------------
# Desempeno Out-sample

valor <- 100 
t <- nrow(retornosfm)+1
n <- ncol(retornosfm)

# Retornos del portafolio
rpomegafm <- retornosfm%*%wpomega
(mu.omegafm <- mean(rpomegafm))
(sigma.omegafm <- sd(rpomegafm))

# Valor del portafolio
port.omegafm <- matrix(0, nrow=t)
port.omegafm[1] <- valor
for(i in 2:t){
    port.omegafm[i] <- port.omegafm[i-1]*exp(rpomegafm[i-1])
}

# Benchmark
v.benchmarkfm <- matrix(0, nrow=t)
v.benchmarkfm[1] <- valor
for(i in 2:t){
    v.benchmarkfm[i] <- v.benchmarkfm[i-1]*exp(r.indicefm[i-1])
}

Performancefm <- cbind(port.omegafm,v.benchmarkfm)
Performancefm <- ts(Performancefm,start=c(2021,1), frequency=12)
colnames(Performancefm) <- c("Omega","Benchmark")

windows()
plot(Performancefm[,1], type='l',main="Evaluacion de desempeño", bty="L",
     ylab="Valor del portafolio", xlab="Tiempo (meses)", ylim=c(min(Performancefm), max(Performancefm)))
lines(Performancefm[,2], col='darkgray')
legend("topleft",c("Omega","S&P 500"),
       lty =c(1,1), col=c("black","darkgray"))

# Retorno activo: 
(Raomegafm = mu.omegafm - mu.indicefm)

# Medida Omega:
rpohistfm <- retornosfm%*%wpomega
(omegapfm <- sum(pmax(rpohistfm-h,0))/sum(pmax(h-rpohistfm,0)))
(riskomegapfm <- sum(pmax(h-rpohistfm,0)))

# Indice
(omegaindicefm <- sum(pmax(r.indicefm-h,0))/sum(pmax(h-r.indicefm,0)))
(riskomegap <- sum(pmax(h-r.indicefm,0)))

## ----------------------------------------------------------------
## ----------------------------------------------------------------

# Comparacion de los resultados para todos los portafolios
# Pesos optimos

windows()
par(mfrow = c(3, 2))
barplot(t(wpt), main="Sharpe", axisnames = TRUE, beside = TRUE)
barplot(t(wpmvg), main="PMVG",axisnames = TRUE, beside = TRUE)
barplot(t(wpot), main="Treynor",axisnames = TRUE, beside = TRUE)
barplot(t(wpts), main="Sortino",axisnames = TRUE, beside = TRUE)
barplot(t(wpomega), main="Omega",axisnames = TRUE, beside = TRUE)


#----------------------------------------------------------------------
#----------------------------------------------------------------------

## Ejemplo 7.6
## Coeficiente de Información: IR

## Calculos in-sample
# Retomar series historicas de retornos de los portafolios

rpsharpe 
rptreynor
rpsortino
rpomega
rb <- r.indice

# Retornos activos anualizados
(Rasharpe <- mean(rpsharpe-rb)*12)
(Ratreynor <- mean(rptreynor-rb)*12)
(Rasortino <- mean(rpsortino-rb)*12)
(Raomega <- mean(rpomega-rb)*12)

# TE
(tesharpe <- sd(rpsharpe-rb)*sqrt(12))
(tetreynor <- sd(rptreynor-rb)*sqrt(12))
(tesortino <- sd(rpsortino-rb)*sqrt(12))
(teomega <- sd(rpomega-rb)*sqrt(12))

# IR:
(irsharpe <- Rasharpe/tesharpe)
(irtreynor <- Ratreynor/tetreynor)
(irsortino <- Rasortino/tesortino)
(iromega <- Raomega/teomega)

## ----------------------------------------------

## Calculos Out-sample
# Series historicas de retornos

rpsharpefm
rptreynorfm
rpsortinofm
rpomegafm
rbfm <- r.indicefm

# Retornos activos anualizados
(Rasharpefm <- mean(rpsharpefm-rbfm)*12)
(Ratreynorfm <- mean(rptreynorfm-rbfm)*12)
(Rasortinofm <- mean(rpsortinofm-rbfm)*12)
(Raomegafm <- mean(rpomegafm-rbfm)*12)

# TE
(tesharpefm <- sd(rpsharpefm-rbfm)*sqrt(12))
(tetreynorfm <- sd(rptreynorfm-rbfm)*sqrt(12))
(tesortinofm <- sd(rpsortinofm-rbfm)*sqrt(12))
(teomegafm <- sd(rpomegafm-rbfm)*sqrt(12))

# IR:
(irsharpefm <- Rasharpefm/tesharpefm)
(irtreynorfm <- Ratreynorfm/tetreynorfm)
(irsortinofm <- Rasortinofm/tesortinofm)
(iromegafm <- Raomegafm/teomegafm)

#----------------------------------------------------------------------
#----------------------------------------------------------------------

## Ejemplo 7.7
## Optimización Tracking-Error

# Informacion del Benchmark
(mu <- mu*12)
(sigma <- sigma*sqrt(12))
cov <- cov*12
(mu.indice <- mean(r.indice)*12)
(sigma.indice <- sd(r.indice)*sqrt(12))


# Construcion frontera TE

n <- length(mu)
unos <- rep(1,n)
x <- t(mu)%*%solve(cov)%*%mu
y <- t(mu)%*%solve(cov)%*%unos
z <- t(unos)%*%solve(cov)%*%unos
d <- x - y^2/z

(rpmvg <- y/z)
(sigmapmvg <- sqrt(1/z))

mu.indice <- 0.25
(delta1 <- mu.indice - rpmvg)
(delta2 <- sigma.indice^2 - sigmapmvg^2)

# Calculo del vector de pesos: x
tep <- 0.03
sigmaa2 <- tep^2
x <- as.numeric(sqrt(sigmaa2/d))*solve(cov) %*% (mu - as.numeric(y/z)*unos)
round(x,4)

# Refeencia: PMVG
(wp <- wpmvg+x)
sum(wp)

# Conparación de resultados
(ra <- t(wp)%*%mu)
(sigmapa <- sqrt(t(wp)%*%cov%*%wp))

rpmvg
sigmapmvg

# Frontera FE

nportte <- 100
tepj <- seq(0,0.05,length=nportte)
sigmaa2j <- tepj^2

#wpote <- matrix(c(0), ncol=n, nrow=nportte) 
#wb <- rep(1/n,n)
rb <- t(wb)%*%mu 
sigmab <- sqrt(t(wb)%*%cov%*%wb)

rpote <- matrix(c(0), nrow=nportte)
sigmapote <- matrix(c(0), nrow=nportte)

for(i in 1:nportte){
    xj <- as.numeric(sqrt(sigmaa2j[i]/d))*solve(cov)%*%(mu - as.numeric(y/z)*unos)
    wpote <- wb+xj
    rpote[i] <- t(wpote)%*%mu
    sigmapote[i] <- sqrt(t(wpote)%*%cov%*%wpote)
}

windows()
plot(sigma,mu, ylim=c(0,max(mu)*1.15), 
     xlim=c(0,max(sigma)*1.15),
     main="Plano riesgo-retorno",xlab= "Riesgo",
     ylab="Retorno esperado")
points(sigmapo*sqrt(12),rpo*12, type="l") 
points(sigmapmvg,rpmvg, col="black",pch=19)
text(sigmapmvg,rpmvg,labels="PMVG",pos = 2, cex=0.7)
text(sigma,mu,labels = activos,pos = 2,cex=0.6) 
points(sigmab,rb, col="black",pch=19)
text(sigmab,rb,labels="Bench",pos = 2, cex=0.7)
#points(sigma.indice,mu.indice, col="black",pch=19)
#text(sigma.indice,mu.indice,labels="Bench",pos = 2, cex=0.7)
lines(sigmapote,rpote,col="blue") 

windows()
plot(sigmapote,rpote, type="l",col="blue") 


# Construccion frontera TE Constante
wb <- wb#rep(1/n,n)
tep <- 0.03
sigmaa2 <- tep^2

delta1
delta2
(phi <- sigmapmvg^2 - sigma.indice^2 - sigmaa2)

(lambda3 <- -(delta1/delta2)-(phi/delta2)*sqrt((d*delta2-delta1^2)/(4*sigmaa2*delta2-phi)))
(lambda1 <- -(lambda3+y)/z)
(lambda2 <- -(-2*sqrt((d*delta2-delta1^2)/(4*sigmaa2*delta2-phi))-lambda3))

(x <- -as.numeric((1/(lambda2+lambda3)))*solve(cov)%*%(mu + as.numeric(lambda1)*unos+ as.numeric(lambda3)*cov%*%wb))
round(x,4)

# Refeencia: PMVG
wb <- wpmvg
wp <- wb+x
sum(wp)

# Conparación de resultados
(ra <- t(wp)%*%mu)
(sigmapa <- sqrt(t(wp)%*%cov%*%wp))

rpmvg
sigmapmvg

# --------------------------
# COnstrucción de la elipse
nport <- 100
#wb <- rep(1/n,n)
tep <- 0.03
sigmaa2 <- tep^2
delta1
delta2
sigmapel <- seq(sigmaa2,sigmaa2*1.5,length=nport)
(phi <- sigmapel^2 - sigma.indice^2 - sigmaa2)

rpoel <- matrix(c(0), nrow=nport)
sigmapoel <- matrix(c(0), nrow=nport)

for(i in 1:nportte){
    lambda3 <- -(delta1/delta2)+(phi[i]/delta2)*sqrt((d*delta2-delta1^2)/(4*sigmaa2*delta2-phi[i]))
    lambda1 <- -(lambda3+y)/z
    lambda2 <- (-2*sqrt((d*delta2-delta1^2)/(4*sigmaa2*delta2-phi))-lambda3)
    xj <- as.numeric(sqrt(sigmapel[i]/d))*solve(cov)%*%(mu - as.numeric(y/z)*unos)
    wpote <- wb+xj
    rpoel[i] <- t(wpote)%*%mu
    sigmapoel[i] <- sqrt(t(wpote)%*%cov%*%wpote)
}

windows()
plot(sigma,mu, ylim=c(0,max(mu)),xlim=c(0,max(sigma)),
     main="Plano riesgo-retorno",xlab= "Riesgo",
     ylab="Retorno esperado")
points(sigmapo*sqrt(12),rpo*12, type="l") 
points(sigmapmvg,rpmvg, col="black",pch=19)
text(sigmapmvg,rpmvg,labels="PMVG",pos = 2, cex=0.7)
text(sigma,mu,labels = activos,pos = 2,cex=0.6) 
points(sigmab,rb, col="black",pch=19)
text(sigmab,rb,labels="Bench",pos = 2, cex=0.7)
#lines(sigmapoel,rpoel,col="red") 
points(sigmapa,ra, col="black",pch=19)
lines(sigmapote,rpote,col="blue") 

windows()
plot(sigmapoel,rpoel, type="l",col="blue") 


## -------------------------------------------------------------

