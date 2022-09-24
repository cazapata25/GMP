# -*- coding: utf-8 -*-
"""cap3.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1DXBEftJ_i_gbJ9jB8QNiJJRhB7y2Wn2Q

# Gestión Moderna de portafolio
### Autores Bernardo León y Carlos Zapata
### Copyright 2022

## Capitulo 3: Modelo CAPM

### Ejemplos 3.1 y 3.2. Estimación de parámetros: $\alpha_i$ y $\beta_i$
"""

# Commented out IPython magic to ensure Python compatibility.
#Librerías 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
# %pip install yfinance
import yfinance as yf

## Estimaciones para AAPL
# Información histórica para AAPL 
fechai = "2015-12-01"
fechaf = "2020-12-31"
periodicidad = "1Mo"
aapl = ["AAPL"]
p_aapl = yf.download(aapl,start=fechai,end=fechaf,interval=periodicidad)['Adj Close'].dropna()
r_aapl = np.log(p_aapl/p_aapl.shift(1)).dropna()
aapl = ["AAPL"]
indice = yf.download("^GSPC",start=fechai,end=fechaf,interval=periodicidad)['Adj Close'].dropna()
r_indice = np.log(indice/indice.shift(1)).dropna()
rm = r_indice.mean()
sigmam = r_indice.std()
round(rm,4), round(sigmam,4)

# Estimación lineal
slope, intercept, r, p, std_err = stats.linregress(r_indice, r_aapl)
def linestim(x):
  return intercept + slope * x

model = list(map(linestim, r_indice))

# Gráfico de dispersión
fig = plt.figure(figsize = (6, 5))
ax = fig.add_subplot(111)
plt.plot(r_indice, r_aapl, ".", color = 'gray')
plt.plot(r_indice, model)
plt.xlabel("Retornos del S&P 500")
plt.ylabel("Retornos de AAPL")
plt.show()

# Estimacion betas
beta = round(slope,4)
alpha = round(intercept,4)
print('Beta AAPL: ', beta )
print('Alpha AAPL: ', alpha)

## Estimaciones para AAPL, AMZN, GOOG y MSFT
# Información histórica 
fechai = "2015-12-01"
fechaf = "2020-12-31"
periodicidad = "1Mo"
activos = ["AAPL","AMZN","GOOG","MSFT"]
precios = yf.download(activos,start=fechai,end=fechaf,interval=periodicidad)['Adj Close'].dropna()
retornos = np.log(precios/precios.shift(1)).dropna()
indice = yf.download("^GSPC",start=fechai,end=fechaf,interval=periodicidad)['Adj Close'].dropna()
r_indice = np.log(indice/indice.shift(1)).dropna()

# Estimaciones AAPL
slope, intercept, r, p, se = stats.linregress(r_indice, retornos['AAPL'])
def linestim(x):
  return intercept + slope * x

beta_aapl = slope
alpha_aapl = intercept
print("Beta: ", round(beta_aapl,4))
print("Alpha: ", round(alpha_aapl,4))
print(f"R-squared: {r**2:.6f}")

# Estimación varianza del error
t=60
model = list(map(linestim, r_indice))
error_aapl = (((retornos['AAPL']-model)**2).sum())/(t-1)
round(error_aapl,4)

# Estimaciones AMZN
slope, intercept, r, p, se = stats.linregress(r_indice, retornos['AMZN'])
def linestim(x):
  return intercept + slope * x

beta_amzn = slope
alpha_amzn = intercept
print("Beta: ", round(beta_amzn,4))
print("Alpha: ", round(alpha_amzn,4))
print(f"R-squared: {r**2:.6f}")

# Estimación varianza del error
t=60
model = list(map(linestim, r_indice))
error_amzn = (((retornos['AMZN']-model)**2).sum())/(t-1)
round(error_amzn,4)

# Estimaciones GOOG
slope, intercept, r, p, se = stats.linregress(r_indice, retornos['GOOG'])
def linestim(x):
  return intercept + slope * x

beta_goog = slope
alpha_goog = intercept

print("Beta: ", round(beta_goog,4))
print("Alpha: ", round(alpha_goog,4))
print(f"R-squared: {r**2:.6f}")

# Estimación varianza del error
t=60
model = list(map(linestim, r_indice))
error_goog = (((retornos['GOOG']-model)**2).sum())/(t-1)
round(error_goog,4)

# Estimaciones MSFT
slope, intercept, r, p, se = stats.linregress(r_indice, retornos['MSFT'])
def linestim(x):
  return intercept + slope * x

beta_msft = slope
alpha_msft = intercept
print("Beta: ", round(beta_msft,4))
print("Alpha: ", round(alpha_msft,4))
print(f"R-squared: {r**2:.6f}")

# Estimación varianza del error
t=60
model = list(map(linestim, r_indice))
error_msft = (((retornos['MSFT']-model)**2).sum())/(t-1)
round(error_msft,4)

# Estimación de retornos esperados 
rf = 0.02/12
re_aapl = rf + beta_aapl*(rm-rf)
re_amzn = rf + beta_amzn*(rm-rf)
re_goog = rf + beta_goog*(rm-rf)
re_msft = rf + beta_msft*(rm-rf)

print("E(R)_aapl: ", round(re_aapl,4))
print("E(R)_amzn: ", round(re_amzn,4))
print("E(R)_goog: ", round(re_goog,4))
print("E(R)_msft: ", round(re_msft,4))

"""## Ejemplo 3.3
## Modelo de Sharpe para portafolio óptimo
"""

round(rm,4), round(sigmam,4)

mu = np.array(retornos.mean())
n = len(mu)
unos = np.repeat(1,n)
sigma_error = np.array([error_aapl, error_amzn, error_goog, error_msft, sigmam**2])
betas = np.array([beta_aapl, beta_amzn, beta_goog, beta_msft, -1])
alphas = np.array([alpha_aapl, alpha_amzn, alpha_goog, alpha_msft, rm])
uno = np.array([1,1,1,1,0])
sigmaei = np.diagflat(sigma_error)
alb = np.concatenate(([alphas],[uno],[betas])).T
G = alb.T @ np.linalg.inv(sigmaei) @ alb

# Vector resultados
res = np.array([0.02, 1, 0])

# Calculo pesos optimos
wpo =  np.linalg.inv(sigmaei) @ alb @ np.linalg.inv(G) @ res
wpo[0:4]

# Pesos del PMVG
fig = plt.figure(figsize = (6, 5))
plt.bar(activos, wpo[0:4], width = 0.4)
plt.ylabel("Part. (%)")
plt.show()

# Construcción de la FE
cov = retornos.cov()
sigma = retornos.std()
rpmin = mu.min()
rpmax = mu.max()
nport = 100

j = np.linspace(rpmin,rpmax, nport) 
wpo = np.zeros((nport, n+1))
rpo = np.zeros((nport,1))
sigmapo = np.zeros((nport,1))

for i in range(nport):
  res = np.array([j[i], 1, 0])
  wj =  np.linalg.inv(sigmaei) @ alb @ np.linalg.inv(G) @ res
  wpo[i,:] = wj
  rpo[i] = wj[0:n] @ mu
  sigmapo[i] = np.sqrt( wj[0:n].T @ cov @ wj[0:n])

# Plano Riesgo-Retorno
fig = plt.figure(figsize = (6, 5))
ax = fig.add_subplot(111)
plt.plot(sigmapo, rpo)
plt.xlabel("Riesgo")
plt.ylabel("Retorno esperado")
plt.show()

"""## Ejemplos 3.4 y 3.5
## Clasificación y optimización ussando el modelo de Treynor
"""

## Información histórica 
fechai = "2015-12-01"
fechaf = "2020-12-31"
periodicidad = "1Mo"
activos = ["AAPL","ABT","AMZN","CAT","CSX","META","GOOG","HD","JNJ","MSFT","MCD","V"]
precios = yf.download(activos,start=fechai,end=fechaf,interval=periodicidad)['Adj Close'].dropna()
retornos = np.log(precios/precios.shift(1)).dropna()
indice = yf.download("^GSPC",start=fechai,end=fechaf,interval=periodicidad)['Adj Close'].dropna()
r_indice = np.log(indice/indice.shift(1)).dropna()
rm = r_indice.mean()*12
sigmam = (r_indice.std())*np.sqrt(12)
round(rm,4), round(sigmam,4)

# Retornos esperados y volatilidades de los activos
mu = retornos.mean()*12
cov = retornos.cov()*12
var = np.diag(cov)
sigma = (retornos.std())*np.sqrt(12)
estimaciones = pd.concat([mu, sigma],axis=1).T
estimaciones.index=["Retorno","Riesgo"]
estimaciones

# Estimaciones: betas y varerror
n = len(mu)
t = 60
def linestim(x):
  return intercept + slope * x

betas = np.zeros((n,1))
varerror = np.zeros((n,1))

for i in range(n):
  slope, intercept, r, p, se = stats.linregress(r_indice, retornos.iloc[:,i])
  model = list(map(linestim, r_indice))
  betas[i] = slope
  varerror[i] = ((((retornos.iloc[:,i]-model)**2).sum())/(t-1))*12

# Resultados estimaciones: betas y var error
betas.T, varerror.T

# Cálculo coeficiente de Treynor
rf = 0
mu = np.array(mu).reshape((n,1))
betas = np.array(betas).reshape((n,1))
treynor = (mu-rf)/betas
treynor.T

estimaciones = np.matrix(np.concatenate([np.array(treynor).reshape((n,1)),np.array(mu).reshape((n,1)), 
                                      np.array(sigma).reshape((n,1)), np.array(betas).reshape((n,1)),
                                      np.array(varerror).reshape((n,1))],1))

tabla = pd.DataFrame(estimaciones)
tabla.columns=['Treynor','Retorno','Riesgo','Beta','VarError']
tabla.index=activos 
tabla

# Se ordenan los activos usando el coef. de Treynor
tabla_sort = tabla.sort_values(by=['Treynor'],ascending=False)
tabla_sort

# Se calcula la tasa de corte y se estiman los pesos óptimos
ratio1 = ((tabla_sort['Retorno']-rf)*tabla_sort['Beta'])/tabla_sort['VarError']
ratio2 = tabla_sort['Beta']**2/tabla_sort['VarError']

sumacu1 = ratio1.cumsum()
sumacu2 = ratio2.cumsum()
tasac = (sigmam**2 * sumacu1)/(1+ sigmam**2 * sumacu2)
cuttoff = tasac.max()

Zi = (tabla_sort['Beta']/tabla_sort['VarError'])*(tabla_sort['Treynor']-cuttoff)
Zi = Zi.mask(Zi<0, 0)
wi = round(Zi/sum(Zi),4)
wpot = round(wi,4)
wpot

# Pesos del Portafolio de Treynor
fig = plt.figure(figsize = (6, 5))
plt.bar(activos, wpot, width = 0.4)
plt.ylabel("Part. (%)")
plt.show()