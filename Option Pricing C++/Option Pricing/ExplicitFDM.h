#pragma once
#include <vector>
#include <tuple>
#include <string>
using namespace std;

class ExplicitFDM { // Abstract Base class.
	typedef vector<vector<tuple<double, string>>> Grid;
public:
	ExplicitFDM(
		string type, // "Call" or "Put"
		string ticker,
		double SMAX,
		double K,
		double T,
		double r,
		double sig,
		double q,  // Dividend yield.
		size_t M,  // Number of intervals in price.
		size_t N   // Number of periods in time.
		);
	virtual ~ExplicitFDM();

	void set_boundary();
	virtual void do_exiplicit_fdm() = 0;

	void output_results();

protected:
	string type_;
	string style_; // "European" or "American"
	string ticker_;
	double max_price_;
	double strike_price_;
	double time_to_maturity_;
	double risk_free_rate_;
	double volatility_;
	double yield_;
	size_t num_of_price_intervals_;
	size_t num_of_periods_;
	double delta_t_;
	double delta_s_; 
	Grid grid_; // (N + 1) * (M + 1) matrix. Access: (i, j)
	            // 2 elements: option price; early exercise "Yes" or "No"
};

// ===================== European ==============================
class EuropeanExplicitFDM : public ExplicitFDM {
public:
	EuropeanExplicitFDM(
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
		);
	virtual ~EuropeanExplicitFDM();

	void do_exiplicit_fdm() override;
};

// ======================= American =========================================
class AmericanExplicitFDM : public ExplicitFDM {
public:
	AmericanExplicitFDM(
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
		);
	virtual ~AmericanExplicitFDM();

	void do_exiplicit_fdm() override;
};


