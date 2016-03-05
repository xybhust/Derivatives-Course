#include <algorithm>
#include <fstream>
#include "ExplicitFDM.h"
using namespace std;


ExplicitFDM::ExplicitFDM(
	string type,
	string ticker,
	double SMAX,
	double K,
	double T,
	double r,
	double sig,
	double q,
	size_t M,  
	size_t N) :
	type_(type), ticker_(ticker), max_price_(SMAX), strike_price_(K), 
	time_to_maturity_(T), risk_free_rate_(r), volatility_(sig), yield_(q), 
	num_of_price_intervals_(M),	num_of_periods_(N)
{
	delta_t_ = T / N;
	delta_s_ = SMAX / M;
	grid_.resize(num_of_periods_ + 1);
	for (auto& x : grid_) {
		x.resize(num_of_price_intervals_ + 1);
	}
}

ExplicitFDM::~ExplicitFDM() { }

void ExplicitFDM::set_boundary() {
	if (type_ == "Call") {
		for (size_t i = 0; i <= num_of_periods_; ++i) {
			get<0>(grid_[i][0]) = 0.0;
			get<0>(grid_[i][num_of_price_intervals_]) = max_price_;
		}
		for (size_t j = 0; j <= num_of_price_intervals_; ++j) {
			get<0>(grid_[num_of_periods_][j]) = 
				max(0.0, j * delta_s_ - strike_price_);
		}
	}
	else if (type_ == "Put") {
		for (size_t i = 0; i <= num_of_periods_; ++i) {
			if (style_ == "American") {
				get<0>(grid_[i][0]) = strike_price_;
			} 
			else if (style_ == "European") {
				get<0>(grid_[i][0]) = strike_price_ *
					exp(-risk_free_rate_ * i * delta_t_);
			}
			get<0>(grid_[i][num_of_price_intervals_]) = 0.0;
		}
		for (size_t j = 0; j <= num_of_price_intervals_; ++j) {
			get<0>(grid_[num_of_periods_][j]) = 
				max(0.0, strike_price_ - j * delta_s_);
		}
	}
}

void ExplicitFDM::output_results() {
	string file = "Explicit" + style_ + type_ + ticker_ + ".csv";
	ofstream output(file);
	output << ticker_ << "," << style_ << "," << type_ << ",\n";
	output << ",Time to Maturity,";
	for (size_t i = 0; i <= num_of_periods_; ++i) {
		output << (num_of_periods_ - i) * delta_t_ <<
			",,,";  
	}
	output << "\nStock Price";
	for (size_t j = 0; j <= num_of_price_intervals_; ++j) {
		output << "\n" << j * delta_s_ << ",,";
		for (size_t i = 0; i <= num_of_periods_; ++i) {
			output << get<0>(grid_[i][j]) << "," << get<1>(grid_[i][j]) << ",,";
		}

	}
	for (size_t j = 0; j <= num_of_price_intervals_; ++j) {
		if (j < 100) {
			printf("Stock Price := %.2f    Option Price := %.2f\n", j * delta_s_,
				get<0>(grid_[0][j]));
		}
		else
			break;
	}
	output.close();
}

// European Method
EuropeanExplicitFDM::EuropeanExplicitFDM(
	string type,
	string ticker,
	double SMAX,
	double K,
	double T,
	double r,
	double sig,
	double q,  // Dividend yield.
	size_t M,  // Number of intervals in price.
	size_t N   // Number of periods in time.
	) : ExplicitFDM(type, ticker, SMAX, K, T, r, sig, q, M, N) { 
	style_ = "European";
}

EuropeanExplicitFDM::~EuropeanExplicitFDM() { }

void EuropeanExplicitFDM::do_exiplicit_fdm() {
	set_boundary();
	for (int i = num_of_periods_ - 1; i >= 0; --i) {
		for (int j = 1; j < num_of_price_intervals_; ++j) {
			// Initial coefficients a, b, c.
			double a = (-0.5 * (risk_free_rate_ - yield_) * j * delta_t_ +
				0.5 * volatility_ * volatility_ * j * j * delta_t_) /
				(1 + risk_free_rate_ * delta_t_);
			double b = (1 - volatility_ * volatility_ * j * j * delta_t_) /
				(1 + risk_free_rate_ * delta_t_);
			double c = (0.5 * (risk_free_rate_ - yield_) * j * delta_t_ +
				0.5 * volatility_ * volatility_ * j * j * delta_t_) /
				(1 + risk_free_rate_ * delta_t_);
			get<0>(grid_[i][j]) = a * get<0>(grid_[i + 1][j - 1]) + 
				b * get<0>(grid_[i + 1][j]) + c * get<0>(grid_[i + 1][j + 1]);
			get<1>(grid_[i][j]) = "No";
		}
	}
}

// American
AmericanExplicitFDM::AmericanExplicitFDM(
	string type,
	string ticker,
	double SMAX,
	double K,
	double T,
	double r,
	double sig,
	double q,  // Dividend yield.
	size_t M,  // Number of intervals in price.
	size_t N   // Number of periods in time.
	) : ExplicitFDM(type, ticker, SMAX, K, T, r, sig, q, M, N) {
	style_ = "American";
}

AmericanExplicitFDM::~AmericanExplicitFDM() { }

void AmericanExplicitFDM::do_exiplicit_fdm() {
	set_boundary();
	for (int i = num_of_periods_ - 1; i >= 0; --i) {
		for (int j = 1; j < num_of_price_intervals_; ++j) {
			// Initial coefficients a, b, c.
			double a = (-0.5 * (risk_free_rate_ - yield_) * j * delta_t_ +
				0.5 * volatility_ * volatility_ * j * j * delta_t_) /
				(1 + risk_free_rate_ * delta_t_);
			double b = (1 - volatility_ * volatility_ * j * j * delta_t_) /
				(1 + risk_free_rate_ * delta_t_);
			double c = (0.5 * (risk_free_rate_ - yield_) * j * delta_t_ +
				0.5 * volatility_ * volatility_ * j * j * delta_t_) /
				(1 + risk_free_rate_ * delta_t_);
			double wait = a * get<0>(grid_[i + 1][j - 1]) + 
				b * get<0>(grid_[i + 1][j]) + c * get<0>(grid_[i + 1][j + 1]);
			double early = (type_ == "Call" ? j * delta_s_ - strike_price_ :
				strike_price_ - j * delta_s_);
			get<0>(grid_[i][j]) = max(wait, early);
			get<1>(grid_[i][j]) = early > wait ? "Yes" : "No";
		}
	}
}
