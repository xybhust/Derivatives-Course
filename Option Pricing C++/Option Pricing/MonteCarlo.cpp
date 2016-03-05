#include <random>
#include <numeric>
#include <chrono>
#include <ctime>
#include "MonteCarlo.h"

MonteCarlo::~MonteCarlo() { }

MonteCarloStrategy::MonteCarloStrategy(
	string type,  
	string ticker,
	double S, 
	double mu, 
	double sigma,
	double r,
	double y,
	double T,  
	size_t N, 
	size_t basis  
	) : 
	strategy_type_(type),
	ticker_(ticker),
	spot_(S),
	exp_return_(mu),
	volatility_(sigma),
	risk_free_rate_(r),
	yield_(y),
	expiration_(T),
	num_of_sims_(N)
{
	num_of_periods_ = static_cast<size_t>(round(expiration_ * basis));  
	delta_t_ = expiration_ / num_of_periods_;  // Length of each time interval. 
}

MonteCarloStrategy::~MonteCarloStrategy() { }

void MonteCarloStrategy::run_simulation() {
	size_t seed = chrono::system_clock::now().time_since_epoch().count() % 100;
	for (size_t i = 0; i < num_of_sims_; ++i) {
		simulate_once(seed + i);
	}
	double sum = accumulate(all_payoffs_.begin(), all_payoffs_.end(), 0.0);
	double mean = sum / num_of_sims_ * exp(-risk_free_rate_ * expiration_);
	printf("%f\n", mean);
}

void MonteCarloStrategy::generate_stock_path(size_t seed) {
	vector<double> stock_path(num_of_periods_ + 1);
	stock_path[0] = spot_;
	default_random_engine generator(seed);  // minstd_rand0 is a standard linear_congruential_engine
	normal_distribution<double> norm(0, 1);

	for (size_t i = 1; i <= num_of_periods_; ++i) {
		double noise = norm(generator);
		stock_path[i] = stock_path[i - 1] * exp((exp_return_ - yield_ -
			0.5 * volatility_ * volatility_) * delta_t_ + \
			volatility_ * sqrt(delta_t_) * noise);
	}
	all_future_stock_prices_.push_back(stock_path[num_of_periods_]);
	all_paths_.push_back(move(stock_path));

}

// =================
MonteCarloNaiveOption::MonteCarloNaiveOption(
	string type,  
	string ticker,
	double S, 
	double mu, 
	double sigma,
	double r,
	double y,
	double T, 
	size_t N, 
	double K,
	size_t basis 
	) :
	MonteCarloStrategy(type, ticker, S, mu, sigma, r, y, T, N, basis), 
	strike_price_(K) { 
	if (type == "Call")
		payoff_ptr_.reset(new PayOffCall(K));
	else if (type == "Put")
		payoff_ptr_.reset(new PayOffPut(K));
}

MonteCarloNaiveOption::~MonteCarloNaiveOption() { }

void MonteCarloNaiveOption::simulate_once(size_t seed) {
	generate_stock_path(seed);
	double payoff = payoff_ptr_->calculate_payoff(
		all_future_stock_prices_.back());
	all_payoffs_.push_back(payoff);
}





