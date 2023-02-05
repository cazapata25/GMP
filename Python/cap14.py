# -*- coding: utf-8 -*-
"""cap14.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1igJxg_EzCLNODArBWo3gwpnHzp1bFEKj

# Gestión Moderna de Portafolio
### Autores Bernardo León y Carlos Zapata
### (C) Copyright 2023

## Capitulo 14: Paridad de Riesgo
"""

# Commented out IPython magic to ensure Python compatibility.
#Librerías usadas
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cp
from scipy import stats
# %pip install --quiet yfinance
import yfinance as yf
import warnings
warnings.filterwarnings("ignore")

"""### Ejemplo 14.2: Paridad de Riesgo Naive"""

# Información histórica dentro de muestra (In-sample) para las acciones
fechai = "2015-12-01"
fechaf = "2020-12-31"
periodicidad = "1Mo"
activos = ["ADBE","MCD","MSCI","MSFT","NEE","PG","RSG","WMT"]
precios = yf.download(activos,start=fechai,end=fechaf,interval=periodicidad)['Adj Close']#.dropna()
retornos = np.log(precios/precios.shift(1)).dropna()
mu = retornos.mean()
cov = retornos.cov()
sigma = retornos.std()

# PMVG y CR
n = len(mu)
def pmvg(cov, ones):  
    covis = np.linalg.inv(cov)  
    w = np.dot(covis, ones)  
    return w / np.sum(np.abs(w)) 
unos = np.repeat(1,n)
wpmvg = pmvg(cov, unos)

crMV = wpmvg * sigma
np.array(crMV)

fig, (ax1, ax2)= plt.subplots(1,2, figsize=(12, 5))
ax1.bar(activos, wpmvg, width = 0.4)
ax1.set_ylabel("Part. (%)")
ax2.bar(activos, crMV, width = 0.4)
ax2.set_ylabel("Part. (%)")
plt.show();

# Paridad de Riesgo Naive
wrpnaive = 1/sigma/((1/sigma).sum())

# Contribución al riesgo (CR):
crnaive = wrpnaive * sigma
crnaive

fig, (ax1, ax2)= plt.subplots(1,2, figsize=(12, 5))
ax1.bar(activos, wrpnaive, width = 0.4)
ax1.set_ylabel("Part. (%)")
ax2.bar(activos, crnaive, width = 0.4)
ax2.set_ylabel("Part. (%)")
plt.show();

"""# Ejemplos 14.3 y 14.4: Risk Parity Vanilla"""

# Forma 1: Root Solve
from scipy.optimize import fsolve

b = np.repeat(1/n,n)
f = lambda x: cov @ x - b/x

x_root = fsolve(f, b)
wRPv = x_root/sum(x_root)
wRPv

# Contribución al riesgo (CR):
crv = wRPv * (cov @ wRPv)

fig, (ax1, ax2)= plt.subplots(1,2, figsize=(12, 5))
ax1.bar(activos, wRPv, width = 0.4)
ax1.set_ylabel("Part. (%)")
ax2.bar(activos, crv, width = 0.4)
ax2.set_ylabel("Part. (%)")
plt.show();

# Forma 2: Metodo de Newton

from scipy.optimize import minimize
n = cov.shape[0]
def fn_convex(x, cov):
    return 0.5 * x @ cov @ x - (1/n)*np.sum(np.log(x))

x0 = np.full((n,), 1/n)
Sigma = cov
result = minimize(fun=fn_convex, x0=x0, args=(Sigma,), method="BFGS")
x_cvx = result.x
wRPv2 = x_cvx/np.sum(x_cvx)
wRPv2

# Contribución al riesgo (CR vanilla):
crvanilla2 = wRPv2 * (Sigma @ wRPv2)
crrvanilla2 = crvanilla2/np.sum(crvanilla2)
np.array(crrvanilla2)