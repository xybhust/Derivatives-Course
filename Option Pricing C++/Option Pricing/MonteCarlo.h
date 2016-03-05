#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <string>
#include <vector>
#include <memory>
#include "PayOff.h"

using namespace std;

class MonteCarlo { // Abstract class for General Monte Carlo method.
public:
	MonteCarlo() = default;
    virtual ~MonteCarlo();

	virtual void run_simulation() = 0;
protected:
	virtual void simulate_once(size_t seed) = 0;
};

class MonteCarloStrategy : public MonteCarlo { // ABC for MC option strategy.
public:
	MonteCarloStrategy(
		string type,  // Strategy type.
		string ticker,
		double S, // Initial stock price.
		double mu,  // Expected log return.
		double sigma,
		double r,
		double y,
		double T,  // Expiration in years.
		size_t N, // Number of simulations.
		size_t basis  // 252 - Daily period.
		);
	virtual ~MonteCarloStrategy();

	virtual void run_simulation() override;
protected:
	virtual void generate_stock_path(size_t seed);
	virtual void simulate_once(size_t seed) = 0;

	string strategy_type_;
	string ticker_;
	double spot_;  
	double exp_return_;  
	double volatility_; 
	double risk_free_rate_;
	double yield_;
	double expiration_;  
	size_t num_of_sims_;  
	size_t num_of_periods_;  // Number of time intervals.
	double delta_t_;  // Length of each time interval.
	vector<vector<double>> all_paths_;  // Stock evolution paths. Each path has
	                                    // size: (num_of_periods_ + 1).
	vector<double> all_payoffs_;  // Payoffs of each simulation.
	vector<double> all_future_stock_prices_;  // Future prices at expiration.
};

class MonteCarloNaiveOption : public MonteCarloStrategy {
public:
	MonteCarloNaiveOption(
		string type,  // "Call" or "Put".
		string ticker,
		double S, // Initial stock price.
		double mu,  // Expected log return.
		double sigma,
		double r,
		double y,
		double T,  // Expiration in years.
		size_t N, // Number of simulations.
		double K,
		size_t basis = 252  // 252 - Daily period.
		);
	virtual ~MonteCarloNaiveOption();

protected:
	virtual void simulate_once(size_t seed) override;
private:
	double strike_price_;
	unique_ptr<PayOff> payoff_ptr_;
};
#endif