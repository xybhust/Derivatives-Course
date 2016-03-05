
#include "BlackScholesAnlytical.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/tools/roots.hpp>
#include <functional>

namespace
{
	double n(double x)
	{
		boost::math::normal_distribution<double> myNormal(0.0, 1.0);
		return boost::math::pdf(myNormal, x);
	}

	double N(double x)
	{
		boost::math::normal_distribution<> myNormal(0.0, 1.0);
		return boost::math::cdf(myNormal, x);
	}
}

double Black_Scholes_Analytical::price_call(double S, double K, double T, double r, double sig, double b)
{
	double tmp = sig * sqrt(T);
	double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
	double d2 = d1 - tmp;

	return (S * exp((b - r) * T) * N(d1)) - (K * exp(-r * T)* N(d2)); 
}


double Black_Scholes_Analytical::price_put(double S, double K, double T, double r, double sig, double b)
{
	double tmp = sig * sqrt(T);
	double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
	double d2 = d1 - tmp;

	return (K * exp(-r * T)* N(-d2)) - (S * exp((b - r) * T) * N(-d1));
}

double Black_Scholes_Analytical::delta_call(double S, double K, double T, double r, double sig, double b)
{
	double tmp = sig * sqrt(T);
	double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;

	return exp((b - r) * T) * N(d1);
}

double Black_Scholes_Analytical::delta_put(double S, double K, double T, double r, double sig, double b)
{
	double tmp = sig * sqrt(T);
	double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;

	return exp((b - r) * T) * (N(d1) - 1.0);
}

double Black_Scholes_Analytical::gamma_call(double S, double K, double T, double r, double sig, double b)
{
	double tmp = sig * sqrt(T);
	double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
	return n(d1) * exp((b - r) * T) / (S * tmp);
}
double Black_Scholes_Analytical::gamma_put(double S, double K, double T, double r, double sig, double b)
{
	return Black_Scholes_Analytical::gamma_call(S, K, T, r, sig, b);
}

double Black_Scholes_Analytical::theta_call(double S, double K, double T, double r, double sig, double b)
{
	double d1 = (log(S / K) + (b + (sig * sig)*0.5) * T) / (sig * sqrt(T));
	double d2 = d1 - sig * sqrt(T);
	double tmp = -S * sig * exp((b - r) * T)*n(d1) / (2 * sqrt(T));
	return (tmp - (b - r) * S * exp((b - r) * T) * N(d1) - r * K * exp(-r * T) * N(d2));
}
double Black_Scholes_Analytical::theta_put(double S, double K, double T, double r, double sig, double b)
{
	double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / (sig * sqrt(T));
	double d2 = d1 - sig * sqrt(T);
	double tmp = -S * sig * exp((b - r) * T) * n(d1) / (2 * sqrt(T));
	return (tmp + (b - r) * S * exp((b - r)*T) * N(-d1) + r * K * exp(-r * T) * N(-d2));
}

double Black_Scholes_Analytical::vega_call(double S, double K, double T, double r, double sig, double b)
{
	double tmp = sig * sqrt(T);
	double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
	return S * sqrt(T) * n(d1) * exp((b - r) * T);
}


double Black_Scholes_Analytical::vega_put(double S, double K, double T, double r, double sig, double b)
{
	return Black_Scholes_Analytical::vega_call(S, K, T, r, sig, b);
}

double Black_Scholes_Analytical::rho_call(double S, double K, double T, double r, double sig, double b)
{
	double d2 = (log(S / K) + (b - (sig*sig)*0.5) * T) / (sig * sqrt(T));
	return (K * T * exp(-r * T)*N(d2));
}
double Black_Scholes_Analytical::rho_put(double S, double K, double T, double r, double sig, double b)
{
	double d2 = (log(S / K) + (b - (sig*sig)*0.5) * T) / (sig * sqrt(T));
	return (-K * T * exp(-r * T) * N(-d2));
}

double Black_Scholes_Analytical::implied_volatility_call(double S_mkt, double S, double K, double T, double r, double b)
{
	auto implied_call = [=](double vol)
	{	
		double tmp = vol * sqrt(T);
		double d1 = (log(S / K) + (b + (vol * vol) * 0.5) * T) / tmp;
		double d2 = d1 - tmp;
		double bs_call = S * exp((b - r) * T) * N(d1) - (K * exp(-r * T)* N(d2));
		double dd1dvol = (vol * T * vol * sqrt(T) - (log(S / K) + (b - vol * vol / 2) * T) * sqrt(T)) / (vol * vol * T);
		double deri_call = S * exp((b - r) * T) * n(d1) * dd1dvol - (K * exp(-r * T)* n(d2)) * (dd1dvol - sqrt(T));
		return std::tr1::make_tuple(bs_call - S_mkt, deri_call);
	};
	int digits = std::numeric_limits<double>::digits;
	double guess = 0.3;
	double min = 0;
	double max = 10;
	return boost::math::tools::newton_raphson_iterate(implied_call, guess, min, max, digits);
}

double Black_Scholes_Analytical::implied_volatility_put(double S_mkt, double S, double K, double T, double r, double b)
{
	auto implied_put = [=](double vol)
	{
		double tmp = vol * sqrt(T);
		double d1 = (log(S / K) + (b + (vol * vol) * 0.5) * T) / tmp;
		double d2 = d1 - tmp;
		double bsrisk_neutral_prob_up_ut = - S * exp((b - r) * T) * N(-d1) + (K * exp(-r * T)* N(-d2));
		double dd1dvol = (vol * T * vol * sqrt(T) - (log(S / K) + (b - vol * vol / 2) * T) * sqrt(T)) / (vol * vol * T);
		double deririsk_neutral_prob_up_ut = S * exp((b - r) * T) * n(-d1) * dd1dvol - (K * exp(-r * T)* n(-d2)) * (dd1dvol - sqrt(T));
		return std::tr1::make_tuple(bsrisk_neutral_prob_up_ut - S_mkt, deririsk_neutral_prob_up_ut);
	};
	int digits = std::numeric_limits<double>::digits;
	double guess = 0.3;
	double min = 0;
	double max = 10;
	return boost::math::tools::newton_raphson_iterate(implied_put, guess, min, max, digits);
}

std::map<double, double> Black_Scholes_Analytical::get_volatility_smile_call(double S_mkt, double S, double T, double r, double b, const std::vector<double>& range_k)
{
	std::map<double, double> volatility_smile;

	auto bind_implied_call = std::bind(Black_Scholes_Analytical::implied_volatility_call, S_mkt, S, std::placeholders::_1, T, r, b);

	for (double K : range_k)
	{
		volatility_smile.insert(std::pair<double, double>(K, bind_implied_call(K)));
	}
	return volatility_smile;
}


std::map<double, double> Black_Scholes_Analytical::get_volatility_smile_put(double S_mkt, double S, double T, double r, double b, const std::vector<double>& range_k)
{
	std::map<double, double> volatility_smile;

	auto bind_implied_put = std::bind(Black_Scholes_Analytical::implied_volatility_put, S_mkt, S, std::placeholders::_1, T, r, b);

	for (double K : range_k)
	{
		volatility_smile.insert(std::pair<double, double>(K, bind_implied_put(K)));
	}
	return volatility_smile;
}