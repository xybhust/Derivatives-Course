# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 22:47:40 2015

@author: Yibing
"""
import pandas.io.data as web
import datetime
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

def calculate_mean_variance(ticker_list, start, end, basis=252):
    """
    Calculate the expected log return and the standard deviation of return 
    for a given stock. 
    yearly - 1
    monthly - 12
    weekly - 52
    daily - 252
    hourly - tbd   
    """
    price = web.DataReader(ticker_list, 'yahoo', start, end)['Adj Close']
    price = pd.DataFrame(price, columns=ticker_list)
    divide = price / price.shift(1)
    divide = divide.dropna()
    log_return = np.log(divide)
    # annualized expected return
    mu = np.mean(log_return, axis=0) * basis
    mu = pd.Series(mu, name='Return')
    # annualized volatility, which is adjusted to unbiasedness
    sigma = np.std(log_return, axis=0, ddof=1) * np.sqrt(basis)
    sigma = pd.Series(sigma, name='Volatility')
    return mu, sigma

def calculate_correlation(ticker_list, start, end, basis=252):
    """
    correlation between the log return.
    """
    price = web.DataReader(ticker, 'yahoo', start, end)['Adj Close']
    price = pd.DataFrame(price, columns=ticker_list)
    divide = price / price.shift(1)
    divide = divide.dropna()
    log_return = np.log(divide)
    
    # compute covariance and correlation matrix
    cov = np.cov(log_return, rowvar=False)
    corr_mat = np.corrcoef(log_return, rowvar=False)
    return log_return, cov, corr_mat

def sigma_ij(cov, sigma, i, j, basis=252):
    """
    Calculate sigma_ij recursively.
    
    Params:
    cov - covariance matrix.
    sigma - volatility list.
    """
    assert i >= j, "i cannot be less than j"
    if i==1 and j==1:
        return sigma[0] # np index starts from 0.
    elif j==1:
        return cov[i-1, 0]/sigma_ij(cov, sigma, 1, 1)*basis
    elif i>j:
        tmp = 0.
        for p in range(1, j):
            tmp += sigma_ij(cov, sigma, j, p)*sigma_ij(cov, sigma, i, p)
        return (cov[i-1, j-1]*basis - tmp) / sigma_ij(cov, sigma, j, j)
    elif i==j:
        tmp = 0.
        for p in range(1, j):
            tmp += sigma_ij(cov, sigma, i, p)**2
        return np.sqrt(sigma[i-1]**2 - tmp)
        
    
    
def calculate_capm_return(ticker, start, end, risk_free_rate, basis=252):
    """
    CAPM:
    r = rf + b(rm - rf) + e
    r - Log-return of the stock
    rm - Log-return of the benchmark
    rf - Risk-free rate
    b - Beta
    """
    # S&P 500 index as benchmark
    price = web.DataReader((ticker, '^GSPC'), 'yahoo', start, end)['Adj Close']
    plt.figure(1, figsize = (8, 10))
    plt.subplot(2, 1, 1)
    plt.plot(price.index, price[ticker], color = 'blue')
    plt.xlabel("Time")
    plt.ylabel("Price")
    plt.title("Price of %s" % ticker)
    plt.subplot(2, 1, 2)
    plt.plot(price.index, price['^GSPC'], color = 'green')
    plt.xlabel("Time")
    plt.ylabel("Price")
    plt.title("Price of S&P 500")    
    plt.grid(True)
    plt.savefig("CAPM.jpg", dpi = 600)
    plt.close(1)
    
    divide = price / price.shift(1)
    divide = divide.dropna()
    log_return = np.log(divide)
    excess_return = log_return - risk_free_rate
    # Run OLS without intercept
    y = excess_return[ticker]
    x = excess_return['^GSPC']
    reg = sm.OLS(y, x).fit() 
    beta = reg.params[0]
    sigma = np.sqrt(reg.mse_total * basis) 
    mu = (risk_free_rate + beta * (np.mean(log_return['^GSPC']) - 
        risk_free_rate)) * basis
    return mu, sigma, beta
    
if __name__ == '__main__':
#    ticker = ['^GSPC', 'GE', 'AAPL', 'XOM']
#    start = datetime.date(2015, 9, 3)
#    end = datetime.date(2016, 2, 20)
#    mu, sigma = calculate_mean_variance(ticker, start, end)
#    log_return, cov, corr_mat = calculate_correlation(ticker, start, end)
#    
#    vol_mat = np.zeros((4, 4)) # sigma_ij
#    r_mat = np.zeros((4, 4)) # r_ij = sigma_ij / sigma_i for j <= i
#    for i in range(4):
#        for j in range(i+1):
#            vol_mat[i, j] = sigma_ij(cov, sigma, i+1, j+1)
#            r_mat[i, j] = vol_mat[i, j] / sigma[i]
            
    ticker = ['^GSPC', 'GE', 'XOM']
    start = datetime.date(2015, 9, 3)
    end = datetime.date(2016, 2, 20)
    mu, sigma = calculate_mean_variance(ticker, start, end)
    log_return, cov, corr_mat = calculate_correlation(ticker, start, end)
    
    vol_mat = np.zeros((3, 3)) # sigma_ij
    r_mat = np.zeros((3, 3)) # r_ij = sigma_ij / sigma_i for j <= i
    for i in range(3):
        for j in range(i+1):
            vol_mat[i, j] = sigma_ij(cov, sigma, i+1, j+1)
            r_mat[i, j] = vol_mat[i, j] / sigma[i]
    