# -*- coding: utf-8 -*-
"""cap10.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/10INpxDFn3cl6mY_tvV4wYI2frvP_lvve

# Gestión Moderna de Portafolio
### Autores Bernardo León y Carlos Zapata
### (C) Copyright 2023

## Capitulo 10: Portafolio internacional
"""

# Commented out IPython magic to ensure Python compatibility.
#Librerías usadas
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
# %pip install --quiet yfinance
import yfinance as yf
import warnings
warnings.filterwarnings("ignore")

"""### Ejemplo 10.1: Modelo Semivarianza
##### Retorno sistematico: BR + X
##### Varianza sistematica: V(BR+X)
"""

# Información histórica para los mercados locales (DJI,N100,N225) y 
# las tasas de cambio (USDEUR=X, USDJPY=X)
fechai = "2015-12-01"
fechaf = "2020-12-31"
periodicidad = "1Mo"
mercados = ['^DJI','^N100','^N225']
precios = yf.download(mercados,start=fechai,end=fechaf,interval=periodicidad)['Adj Close'].dropna()
retornos = np.log(precios/precios.shift(1)).dropna()
tasas_c = ['USDEUR=X','USDJPY=X']
p_tasas = yf.download(tasas_c,start=fechai,end=fechaf,interval=periodicidad)['Adj Close'].dropna()
r_tasas = np.log(p_tasas/p_tasas.shift(1)).dropna()
retornosx = retornos.join(r_tasas)
rx = retornosx.mean()*12
omega = retornosx.cov()*12
corr = retornosx.corr()
print(rx), print(omega), print(corr)

# Load data MSCI
msci= pd.read_csv("msci.csv")[['MSCI World']]
msci.index = precios.index
r_msci = np.log(msci/msci.shift(1)).dropna()
rm = r_msci.mean()*12
sigmam = r_msci.std()*np.sqrt(12)
round(rm,4), round(sigmam,4)

# Estimación de betas - mercados locales
ret = retornos[['^DJI','^N100','^N225']]
n = 3
betas = np.zeros((n,1))
for i in range(n):
    modelo = modelo = sm.OLS(ret.iloc[:,i], sm.add_constant(r_msci)).fit()
    betas[i] = modelo.params[1]
print(betas)

diagbetas = np.diagflat(betas)

# Retorno y riesgo sistemático del portafolio
# Calculos para el portafolio inicial

w = [0.4,0.3,0.3]
f1 = [0,1,0]
f2 = [0,0,1]
B =  np.column_stack([diagbetas, f1, f2]).T
Bw = B @ w
rsp = Bw @ rx
sigmasp = np.sqrt(Bw.T@omega@Bw)
print('E(R+X): ', round(rsp,4)), print('V(R+X): ', round(sigmasp,4))

"""### Calculos Retorno esperado y Varianza de Markowitz"""

# Calculos para el portafolio inicial
w = [0.4,0.3,0.3]
mu = ret.mean()*12
sigma = ret.std()*np.sqrt(12)
datos = pd.concat([mu, sigma], axis=1).T
cov = ret.cov()*12
rp = w @ mu
sigmap = np.sqrt(w @ cov @ w)
print('E(Rp): ', round(rp,4)), print('V(R): ', round(sigmap,4))

"""#### Optimización del portafolio"""

import cvxpy as cp

# Optimización Markowitz
n = 3
w = cp.Variable(n)
risk = cp.quad_form(w, cov)
objective = cp.Minimize(risk)
constraints = [cp.sum(w) == 1, w >= 0]
prob = cp.Problem(objective, constraints)
#solvers = ['SCS']
prob.solve()
wpo = pd.DataFrame(w.value)
round(wpo,4)

# Optimización Markowitz
n = 3
w = cp.Variable(n)
Bw = B @ w
sist_risk = cp.quad_form(Bw, omega)
objective = cp.Minimize(sist_risk)
constraints = [cp.sum(w) == 1, w >= 0]
prob = cp.Problem(objective, constraints)
#solvers = ['SCS']
prob.solve()
wps = pd.DataFrame(w.value)
round(wps,4)