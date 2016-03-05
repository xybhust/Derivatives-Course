/*
b = r is the Black¨CScholes stock option model.
b = r - q is the Morton model with continuous dividend yield q.
b = 0 is the Black¨CScholes futures option model.
b = r - R is the Garman and Kohlhagen currency option model, where R is the foreign risk-free interest rate.
*/

#ifndef BLACKSCHOLESANALYTICAL_H
#define BLACKSCHOLESANALYTICAL_H
#include <vector>
#include <map>

class Black_Scholes_Analytical {
public:
	static double price_call(double S, double K, double T, double r, double sig, double b);
	static double price_put(double S, double K, double T, double r, double sig, double b);

	static double delta_call(double S, double K, double T, double r, double sig, double b);
	static double delta_put(double S, double K, double T, double r, double sig, double b);

	static double gamma_call(double S, double K, double T, double r, double sig, double b); 
	static double gamma_put(double S, double K, double T, double r, double sig, double b); 

	static double theta_call(double S, double K, double T, double r, double sig, double b);
	static double theta_put(double S, double K, double T, double r, double sig, double b);

	static double vega_call(double S, double K, double T, double r, double sig, double b);
	static double vega_put(double S, double K, double T, double r, double sig, double b);  

	static double rho_call(double S, double K, double T, double r, double sig, double b);
	static double rho_put(double S, double K, double T, double r, double sig, double b);

	static double implied_volatility_call(double S_mkt, double S, double K, double T, double r, double b);
	static double implied_volatility_put(double S_mkt, double S, double K, double T, double r, double b);

	static std::map<double, double> get_volatility_smile_call(double S_mkt, double S, double T, double r, double b, const std::vector<double>& range_k);
	static std::map<double, double> get_volatility_smile_put(double S_mkt, double S, double T, double r, double b, const std::vector<double>& range_k);
};

#endif