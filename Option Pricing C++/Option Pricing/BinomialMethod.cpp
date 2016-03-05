
#include <algorithm>
#include <fstream>
#include "BinomialMethod.h"
using namespace std;

BinomialMethod::BinomialMethod(
	string type,
	string ticker,
	size_t num_of_periods,
	double S,
	double K,
	double discount,
	double up,
	double down,
	double risk_neutral_q
    ) :
	type_(type),
	ticker_(ticker),
	num_of_periods_(num_of_periods),
	spot_(S),
	strike_price_(K),
	discount_(discount),
	up_(up),
	down_(down),
	risk_neutral_q_(risk_neutral_q),
	lattice_(Lattice<double, 2>(num_of_periods)),
	payoff_ptr_(nullptr)
{ 
	if (type == "Call") {
		payoff_ptr_.reset(new PayOffCall(K));
	}
	else if (type == "Put") {
		payoff_ptr_.reset(new PayOffPut(K));
	}
}

BinomialMethod::~BinomialMethod() { }

void BinomialMethod::fill_lattice_forward() {
	for (size_t i = 0; i <= num_of_periods_; ++i) {

		get<0>(lattice_[i][0]) = spot_ * pow(down_, i);
		if (i > 0) {
			size_t length = lattice_[i].size();
			for (int j = 1; j < length; ++j) {
				// Update the binomial tree
				// =================================
				get<0>(lattice_[i][j]) = get<0>(lattice_[i][j - 1])
					/ down_ * up_;
			}
		}
	}
}

void BinomialMethod::output_results() const {
	// Write into csv file.
	string file_name = "Binomial" + ticker_  + style_ + type_ + ".csv";
	ofstream output(file_name);
	double option_price = std::get<1>(lattice_[0][0]);
	printf("Option Price := %.2f\n", option_price);
	output << "Ticker:," << ticker_ << ",\n";
	output << "Option Type:," << type_ << ",\n";
	output << "Stock Price:," << spot_ << ",\n";
	output << "Strike Price:," << strike_price_ << ",\n";
	output << "Option Price:," << option_price << ",\n";
	if (up_ == 1 / down_) {
		double delta = (get<1>(lattice_[1][1]) - get<1>(lattice_[1][0])) /
			(get<0>(lattice_[1][1]) - get<0>(lattice_[1][0]));
		double gamma = ((get<1>(lattice_[2][2]) - get<1>(lattice_[2][1])) /
			(get<0>(lattice_[2][2]) - get<0>(lattice_[2][1])) -
			(get<1>(lattice_[2][1]) - get<1>(lattice_[2][0])) /
			(get<0>(lattice_[2][1]) - get<0>(lattice_[2][0]))) /
			(0.5 * (get<0>(lattice_[2][2]) - get<0>(lattice_[2][0])));
		printf("Delta := %.4f\nGamma := %.4f\n", delta, gamma);
		output << "Delta:," << delta << ",\n";
		output << "Gamma:," << gamma << ",\n";
	}
	if (num_of_periods_ <= 50) {
		output << "Binomial Tree:" << ",\n";
		for (size_t i = 0; i < num_of_periods_; ++i) {
			for (size_t j = 0; j <= i; ++j) {
				output << i << " " << j << ",Option Value:," <<
					get<1>(lattice_[i][j]) << ",Early Exercise:," << 
					get<2>(lattice_[i][j]) << ",\n";
			}
		}
	}
	output.close();


	
}

// Definition of European
// ======================
EuropeanBinomialMethod::EuropeanBinomialMethod(
	string type,
	string ticker,
	size_t num_of_periods,
	double S,
	double K,
	double r,
	double up,
	double down,
	double risk_neutral_q) :
    BinomialMethod(type, ticker, num_of_periods, S, K, r, up, 
					down, risk_neutral_q) 
{
	style_ = "European";
}

EuropeanBinomialMethod::~EuropeanBinomialMethod() { }

void EuropeanBinomialMethod::do_binomial_pricing() {
	fill_lattice_forward(); // Update the selected lattice

    for (int i = lattice_.size() - 1; i >= 0; --i) {
        size_t length = lattice_[i].size();
		if (i == lattice_.size() - 1)
            for (size_t j = 0; j < length; ++j){
				get<1>(lattice_[i][j]) = payoff_ptr_->calculate_payoff(
					get<0>(lattice_[i][j]));
            }
        else
            for (size_t j = 0; j < length; ++j) {
				get<1>(lattice_[i][j]) =
					(risk_neutral_q_ * get<1>(lattice_[i + 1][j + 1]) +
					(1 - risk_neutral_q_) * get<1>(lattice_[i + 1][j])) *
					discount_;
				get<2>(lattice_[i][j]) = "No";
			}
    }
}

// Definition of American
// ======================
AmericanBinomialMethod::AmericanBinomialMethod(
	string type,
	string ticker,
	size_t num_of_periods,
	double S,
	double K,
	double r,
	double up,
	double down,
	double risk_neutral_q) :
	BinomialMethod(type, ticker, num_of_periods, S, K, r, up,
	down, risk_neutral_q)
{ 
	style_ = "American";
}

AmericanBinomialMethod::~AmericanBinomialMethod() { }

void AmericanBinomialMethod::do_binomial_pricing()
{
	fill_lattice_forward(); // Update the selected lattice

	for (int i = lattice_.size() - 1; i >= 0; --i) {
		size_t length = lattice_[i].size();
		if (i == lattice_.size() - 1)
			for (size_t j = 0; j < length; ++j){
				get<1>(lattice_[i][j]) = payoff_ptr_->calculate_payoff(
					get<0>(lattice_[i][j]));
			}
		else
			for (size_t j = 0; j < length; ++j) {
				 double wait =
					 (risk_neutral_q_ * get<1>(lattice_[i + 1][j + 1]) +
					 (1 - risk_neutral_q_) * get<1>(lattice_[i + 1][j])) *
					discount_;
				 double strike_now = payoff_ptr_->calculate_payoff(
					 get<0>(lattice_[i][j]));
				 get<1>(lattice_[i][j]) = (max)(strike_now, wait);
				 get<2>(lattice_[i][j]) = strike_now > wait ? "Yes" : "No";
			}
	}
}
