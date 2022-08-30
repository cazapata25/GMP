# -*- coding: utf-8 -*-
"""cap1.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1J1haPugsROmhEQqQhNV9qeTl0iqAPhfy

# Gestión Moderna de portafolio
### Autores Bernardo León y Carlos Zapata
### Copyright 2022

## Capitulo 1: Medidas de retorno y riesgo de los activos financieros y portafolios de inversión

### Ejemplos 1.2, 1.3 y 1.4
"""

#Librerías usadas
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Commented out IPython magic to ensure Python compatibility.
# %pip install yfinance
import yfinance as yf



# Get the data for the stock AAPL
fechai = "2009-12-01"
fechaf = "2021-12-31"
periodicidad = "1Mo"
activos = ["AAPL","AMZN"]
precios = yf.download(activos,start=fechai,end=fechaf,interval=periodicidad)['Adj Close'].dropna()
precios.plot()

# Ejemplo 1.2
# Calculo de los retornos anualizados
retornos = np.log(precios/precios.shift(1)).dropna()
mu = retornos.mean()*12
round(mu, 4)

# Histograma de los retornos
fig, hist = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))

hist[0].hist(x=retornos['AAPL'], bins=30, alpha=0.5)
hist[0].set_title('AAPL')

hist[1].hist(x=retornos['AMZN'], bins=30, alpha=0.5)
hist[1].set_title('AMZN')

# Medidas de riesgo de los activos
# covarianzas anualizadas
# Ejemplo 1.3
cov = retornos.cov()
cov = cov*12
round(cov,4)

# Volatilidad  anualizada
sigma = retornos.std()
sigma = sigma*np.sqrt(12)
round(sigma,4)

cor = retornos['AAPL'].corr(retornos['AMZN'])
round(cor,4)

# Ejemplo 1.4
# Retorno del portafolio
# w(AAPL) = w(AMZN) = 0.5
n = 2
w = np.repeat(0.5,n)
rp = w @ mu
round(rp,4)

# Riesgo del portafolio
# w(AAPL) = w(AMZN) = 0.5
n = 2
w = np.repeat(0.5,n)
# Formula matricial
sigmap = np.sqrt(w.T @ cov @ w)
round(sigmap,4)

# Formula simple
w1 = w2 = 0.5
sigma1 = sigma['AAPL']
sigma2 = sigma['AMZN']
sigmap = np.sqrt(w1**2 * sigma1**2+ w2**2 * sigma2**2 + 2*w1*w2*sigma1*sigma2*cor)
round(sigmap,4)