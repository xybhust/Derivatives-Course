# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 23:18:29 2015

@author: Yibing
"""
from abc import ABCMeta, abstractmethod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import data_handler

class MonteCarlo(object):
    """
    Abstract class for General Monte Carlo method.
    """  
    metaclass = ABCMeta
    
    def __init__(self):
        pass
    
    @abstractmethod
    def _simulate_once(self):
        pass
    
    @abstractmethod
    def _run_simulation(self):
        pass
    
    @abstractmethod
    def _output_results(self):
        pass
    
    @abstractmethod
    def do_simulation(self):
        """
        Execute Monte Carlo Method, output results.
        """
        pass    
    
class MonteCarloStrategy(MonteCarlo):
    """
    Abstract class for Monte Carlo Option Strategy.
    """
    def __init__(self, strategy_type, ticker, initial_stock_price, exp_return,
                 volatility, risk_free_rate, expiration, num_of_sims, basis):
        """
        Parameters:
        exp_return - Expected annual log return.
        volatility - Annual volatility.    
        expiration - Maturity in years      
        num_of_periods - Number of small intervals for one path. 
        all_paths - A list of list containing initial_price plus subsequent 
            prices. (num_of_periods + 1) in total.  
        all_payoffs - A list ot all payoffs at expiration.
        all_future_stock_prices - A list of all stock prices at expiration.
        """      
        self.strategy_type = strategy_type
        self.ticker = ticker
        self.initial_stock_price = initial_stock_price
        self.exp_return = exp_return
        self.volatility = volatility
        self.risk_free_rate = risk_free_rate
        self.expiration = expiration   
        self.num_of_sims = num_of_sims
        self.num_of_periods = int(np.round(self.expiration * basis))
        self.delta_t = self.expiration / self.num_of_periods # Delta t
        self.is_finished = False
        self.all_paths = []  
        self.all_payoffs = []
        self.all_future_stock_prices = []
        
    @abstractmethod
    def _simulate_once(self):
        """
        This method will be overridden by the child classes.
        """
        pass
    
    def _generate_stock_path(self):
        """
        Generate one path of stock evolution.
        """
        stock_path = [0] * (self.num_of_periods + 1)
        stock_path[0] = self.initial_stock_price
        for i in range(1, self.num_of_periods + 1):
            noise = np.random.normal(0.0, 1.0)
            stock_path[i] = stock_path[i - 1] * np.exp((self.exp_return - 
                0.5 * self.volatility ** 2) * self.delta_t + \
                self.volatility * np.sqrt(self.delta_t) * noise)
        self.all_paths.append(stock_path)
        self.all_future_stock_prices.append(stock_path[-1])   
        
    
    def _run_simulation(self):
        """
        Run the entire simulation.
        """
        for i in range(self.num_of_sims):
            self._simulate_once()
        self.is_finished = True
    
    def _output_results(self):
        """
        Plot figure, save figure to the file.
        """
        if self.is_finished == False:
            self.do_simulation()
        mean_payoff = np.mean(self.all_payoffs)
        mean_stock_price = np.mean(self.all_future_stock_prices)
        discounted_payoff = np.exp(-self.risk_free_rate * self.expiration) * \
            mean_payoff
        print("Expected Future Stock Price := $%.2f" % mean_stock_price)
        print("Expected Payoff of %s Strategy := $%.2f" % (self.strategy_type, 
                                                          mean_payoff))
        print("Present Value of %s Strategy := $%.2f" % (self.strategy_type, 
                                                        discounted_payoff))                                                  
        index = [x * self.delta_t for x in range(self.num_of_periods + 1)]
        # Plot stock price evolution
        plt.figure(1, figsize = (8, 10)) 
        plt.subplot(2, 1, 1)
        plt.title("Stock Price Evolution of %s (Expected: $%.2f)" % (self.ticker,
                                                             mean_stock_price))
        plt.xlabel("Time in Years")
        plt.ylabel("Price")
        plt.subplot(2, 1, 2)
        plt.title("Payoff of %s Strategy ($%.2f Per Share)" % (
                  self.strategy_type, discounted_payoff))
        plt.xlabel("Price at Maturity")
        plt.ylabel("Payoff")
        plt.grid(True)
        for i, curve in enumerate(self.all_paths):
            plt.subplot(2, 1, 1)
            plt.plot(index, curve)
            plt.subplot(2, 1, 2)
            plt.scatter(curve[-1], self.all_payoffs[i])
        plt.savefig("%s.jpg" % self.strategy_type, dpi = 600)
        plt.show()
        plt.close(1)        
        
    def do_simulation(self):
        """
        Execute Monte Carlo Method, output results.
        """
        self._run_simulation()
        self._output_results()     
        result = pd.DataFrame(np.array([self.all_future_stock_prices,
                                       self.all_payoffs]).transpose(), 
                                columns = ['Stock Price', 'Payoff'])
        result.to_csv('%s.csv' % self.strategy_type)
        
class MonteCarloNaiveOption(MonteCarloStrategy):
    """
    This is simpy the option pricing for European call and put.
    """
    def __init__(self, option_type, ticker, initial_stock_price, exp_return, 
                 volatility, risk_free_rate, expiration, num_of_sims, strike,
                 basis = 252):
        """
        Paramaters:
        option_type - Could be 'Call' or 'Put'
        """
        super(MonteCarloNaiveOption, self).__init__('Option', ticker, initial_stock_price, exp_return, 
                        volatility, risk_free_rate, expiration, num_of_sims,
                        basis)
        self.option_type = option_type
        self.strike = strike

    def _simulate_once(self):
        self._generate_stock_path()
        payoff = None
        if self.option_type == 'Call':
            payoff = max(0.0, self.all_future_stock_prices[-1] - self.strike)
        elif self.option_type == 'Put':
            payoff = max(0.0, self.strike - self.all_future_stock_prices[-1])
        self.all_payoffs.append(payoff)


class MonteCarloBullSpread(MonteCarloStrategy):
    """
    Bull spread strategy by call options:
    Buy one European call with K1 and sell one European call
    with the same maturity with K2, where K1 < K2.
    """
    def __init__(self, ticker, initial_stock_price, exp_return, 
                 volatility, risk_free_rate, expiration, num_of_sims, strike_1,
                 strike_2, basis = 252):
        """
        Paramaters:
        """
        super(MonteCarloBullSpread, self).__init__('Bull Spread', ticker, initial_stock_price, exp_return, 
                        volatility, risk_free_rate, expiration, num_of_sims,
                        basis)
        self.strike_1 = strike_1
        self.strike_2 = strike_2
   
    def _simulate_once(self):
        self._generate_stock_path()       
        payoff = max(0.0, self.all_future_stock_prices[-1] - self.strike_1) - \
            max(0.0, self.all_future_stock_prices[-1] - self.strike_2)
        self.all_payoffs.append(payoff)


class MonteCarloBearSpread(MonteCarloStrategy):
    """    
    Bear spread strategy by put options:
    Sell one European put with K1 and buy one European put  
    with the same maturity with K2, where K1 < K2.
    """
    def __init__(self, ticker, initial_stock_price, exp_return, 
                 volatility, risk_free_rate, expiration, num_of_sims, strike_1,
                 strike_2, basis = 252):
        """
        Paramaters:
        """
        super(MonteCarloBearSpread, self).__init__('Bear Spread', ticker, initial_stock_price, exp_return, 
                        volatility, risk_free_rate, expiration, num_of_sims,
                        basis)
        self.strike_1 = strike_1
        self.strike_2 = strike_2
        
    def _simulate_once(self):
        self._generate_stock_path()
        payoff = -max(0.0, self.strike_1 - self.all_future_stock_prices[-1]) + \
            max(0.0, self.strike_2 - self.all_future_stock_prices[-1])
        self.all_payoffs.append(payoff)
       
        
class MonteCarloButterflySpread(MonteCarloStrategy):
    """        
    Butterfly spread strategy by call options:
    Buy two European calls with K1 and K3 respectively, and sell two 
    European calls with the same maturity with K2, where K1 < K2 < K3.
    """
    def __init__(self, ticker, initial_stock_price, exp_return, 
                 volatility, risk_free_rate, expiration, num_of_sims, strike_1,
                 strike_2, strike_3, basis = 252):
        """
        Paramaters:
        """
        super(MonteCarloButterflySpread, self).__init__('Butterfly Spread', ticker, initial_stock_price, 
                        exp_return, volatility, risk_free_rate, expiration, 
                        num_of_sims, basis)
        self.strike_1 = strike_1
        self.strike_2 = strike_2
        self.strike_3 = strike_3
        
    def _simulate_once(self):
        self._generate_stock_path()
        payoff = max(0.0, self.all_future_stock_prices[-1] - self.strike_1) + \
            max(0.0, self.all_future_stock_prices[-1] - self.strike_3) - \
            2 * max(0.0, self.all_future_stock_prices[-1] - self.strike_2)
        self.all_payoffs.append(payoff)  


class MonteCarloStraddle(MonteCarloStrategy):
    """            
    Straddle:
    Buy one European call with K and one European put with
    the same strike price and maturity.
    """
    def __init__(self, ticker, initial_stock_price, exp_return, 
                 volatility, risk_free_rate, expiration, num_of_sims, strike,
                 basis = 252):
        """
        Paramaters:
        """
        super(MonteCarloStraddle, self).__init__('Straddle', ticker, initial_stock_price, 
                        exp_return, volatility, risk_free_rate, expiration, 
                        num_of_sims, basis)
        self.strike = strike

    def _simulate_once(self):
        self._generate_stock_path()
        payoff = max(0.0, self.all_future_stock_prices[-1] - self.strike) + \
            max(0.0, self.strike - self.all_future_stock_prices[-1])
        self.all_payoffs.append(payoff)   


class MonteCarloStrip(MonteCarloStrategy):
    """            
    Strip:
    Buy one European call and two European puts with the K and the same
    and maturity.
    """
    def __init__(self, ticker, initial_stock_price, exp_return, 
                 volatility, risk_free_rate, expiration, num_of_sims, strike,
                 basis = 252):
        """
        Paramaters:
        """
        super(MonteCarloStrip, self).__init__('Strip', ticker, initial_stock_price, 
                        exp_return, volatility, risk_free_rate, expiration, 
                        num_of_sims, basis)
        self.strike = strike

    def _simulate_once(self):
        self._generate_stock_path()
        payoff = max(0.0, self.all_future_stock_prices[-1] - self.strike) + \
            2 * max(0.0, self.strike - self.all_future_stock_prices[-1])
        self.all_payoffs.append(payoff)  


class MonteCarloCollar(MonteCarloStrategy):
    """            
    Collar:
    Buy one stock and one European put with K1, and sell one European call
    with K2, where K1 < S < K2. (Both call and put are OTM.)
    """
    def __init__(self, ticker, initial_stock_price, exp_return, 
                 volatility, risk_free_rate, expiration, num_of_sims, strike_1,
                 strike_2, basis = 252):
        """
        Paramaters:
        """
        super(MonteCarloCollar, self).__init__('Collar', ticker, initial_stock_price, 
                        exp_return, volatility, risk_free_rate, expiration, 
                        num_of_sims, basis)
        self.strike_1 = strike_1
        self.strike_2 = strike_2
        
    def _simulate_once(self):
        self._generate_stock_path()
        payoff = -max(0.0, self.all_future_stock_prices[-1] - self.strike_2) + \
            max(0.0, self.strike_1 - self.all_future_stock_prices[-1]) + \
            self.all_future_stock_prices[-1]
        self.all_payoffs.append(payoff)          
        
if __name__ == '__main__':
    ticker = 'XOM'
    start = datetime.date(2015, 5, 6)
    end = datetime.date(2015, 11, 6)
    S = 84.475
    r = 0.003414
    sigma = 0.2141
    K1 = 80
    K2 = 85
    K3 = 90
    T = (50.0 / 252)  
    N = 1000

    strategy = []
    strategy.append(MonteCarloBullSpread(ticker, S, r, 
                 sigma, r, T, N, K1, K3))            
    strategy.append(MonteCarloBearSpread(ticker, S, r, 
                 sigma, r, T, N, K1, K3)) 
            
    strategy.append(MonteCarloButterflySpread(ticker, S, r, 
                 sigma, r, T, N, K1, K2, K3))      
                
    strategy.append(MonteCarloStraddle(ticker, S, r, 
                 sigma, r, T, N, K2))    
        
    strategy.append(MonteCarloStrip(ticker, S, r, 
                 sigma, r, T, N, K2))    
 
    strategy.append(MonteCarloCollar(ticker, S, r, 
                 sigma, r, T, N, K1, K3))
    for s in strategy[:]:
        s.do_simulation()
#        
    quote = [4.83, 6.02, 1.83, 5.90, 9.77, 85.63]
    names = ['Bull Spread', 'Bear Spread', 'Butterfly Spread', 'Straddle', 
        'Strip', 'Collar']
    plt.figure(1, figsize = (12, 14))
    for i, s in enumerate(names):
        stock = pd.read_csv("%s.csv" % s)['Stock Price'] 
        payoff = pd.read_csv("%s.csv" % s)['Payoff'] - quote[i]
        plt.subplot(3, 2, i + 1)
        plt.title(names[i])
        plt.xlabel("Price at Maturity")
        plt.ylabel("Net Payoff")
        plt.scatter(stock, payoff, s = 0.3)
        plt.grid(True)
    plt.savefig("Net Payoff.jpg", dpi = 600)
    plt.close(1)
    
    plt.figure(2, figsize = (12, 14))
    for i, s in enumerate(names):
        stock = pd.read_csv("%s.csv" % s)['Stock Price'] 
        payoff = pd.read_csv("%s.csv" % s)['Payoff']
        plt.subplot(3, 2, i + 1)
        plt.title(names[i])
        plt.xlabel("Price at Maturity")
        plt.ylabel("Payoff")
        plt.scatter(stock, payoff, s = 0.3)
        plt.grid(True)
    plt.savefig("Payoff.jpg", dpi = 600)
    plt.close(2)