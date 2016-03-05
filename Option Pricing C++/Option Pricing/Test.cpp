
#include <iostream>
#include <Windows.h>
#include "BinomialMethod.h"
#include "BlackScholesAnlytical.h"
#include "MonteCarlo.h"
#include "ExplicitFDM.h"


using namespace std;

void main()
{
	ULONGLONG t1, t2;
	t1 = GetTickCount64();

	// IBM Stock on 12/03/2015.
	//double r = 0.0024;
	//double T = 29.0 / 252;
	//double sigma = 0.1681;
	//double S = 140.815;
	//double K = 140;
	//double y = 0.0373;
	//size_t mc_num = 100000;
	//unsigned long binomial_period = 30; // Binomial periods.
	//double delta_t = T / binomial_period;
	//double u = exp(sigma * sqrt(delta_t));
	//double d = exp(-sigma * sqrt(delta_t));
	//double discount = exp(-r * delta_t);
	//double q = (exp((r - y) * delta_t) - d) / (u - d);
	//double delta_s = S / 60;
	//double SMAX = delta_s * 1000;
	//size_t FDM_M = 1000;
	//size_t FDM_N = 30;
	//AmericanBinomialMethod b("Put", "IBM", binomial_period, S, K, discount, u, d, q);
	//b.do_binomial_pricing();
	//b.output_results();
	////MonteCarloNaiveOption mc("Put", "XOM", S, r, sigma, r, y, T, N, K);
	////mc.run_simulation();


	//AmericanExplicitFDM fdm("Put", "IBM", SMAX, K, T, r, sigma, y, FDM_M, FDM_N);
	//fdm.do_exiplicit_fdm();
	//fdm.output_results();
// ============================================================================
	// VRSN Stock on 12/03/2015.
	//double r = 0.0024;
	//double T = 29.0 / 252;
	//double sigma = 0.2897;
	//double S = 92.04;
	//double K = 90;
	//double y = 0.0;
	//size_t mc_num = 100000;
	//unsigned long binomial_period = 30; // Binomial periods.
	//double delta_t = T / binomial_period;
	//double u = exp(sigma * sqrt(delta_t));
	//double d = exp(-sigma * sqrt(delta_t));
	//double discount = exp(-r * delta_t);
	//double q = (exp((r - y) * delta_t) - d) / (u - d);
	//double delta_s = S / 40;
	//double SMAX = delta_s * 5000;
	//size_t FDM_M = 5000;
	//size_t FDM_N = 30;

	//EuropeanBinomialMethod b("Call", "VRSN", binomial_period, S, K, discount, u, d, q);
	//b.do_binomial_pricing();
	//b.output_results();

	// NIKE
	double r = 0.0013;
	double T = 30.0 / 252;
	double S = 131.36;
	double K = 130;
	double y = 0.0099;
	unsigned long binomial_period = 30; // Binomial periods.
	double delta_t = T / binomial_period;
	double u = 1.017239;
	double d = 1 /u;
	double discount = exp(-r * 1.0 / 252);
	double q = (exp((r - y) * delta_t) - d) / (u - d);

	AmericanBinomialMethod b("Put", "NIKE", binomial_period, S, K, discount, u, d, q);
	b.do_binomial_pricing();
	b.output_results();

	//EuropeanExplicitFDM fdm("Call", "VRSN", SMAX, K, T, r, sigma, y, FDM_M, FDM_N);
	//fdm.do_exiplicit_fdm();
	//fdm.output_results();

	//
	//cout << "Analytical formula: " << Black_Scholes_Analytical::price_call(S, K, T, r, sigma, r - y) << endl;
	//cout << "Delta" << Black_Scholes_Analytical::delta_call(S, K, T, r, sigma, r - y) << endl;
	//cout << "Gamma" << Black_Scholes_Analytical::gamma_call(S, K, T, r, sigma, r - y) << endl;
	//cout << endl;

	t2 = GetTickCount64();
	printf("Use Time: %f seconds\n", (t2 - t1)*1.0 / 1000);
	getchar();
}