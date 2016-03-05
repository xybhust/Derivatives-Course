
#ifndef BINOMIAL_METHOD_H
#define BINOMIAL_METHOD_H

#include <memory>
#include <string>
#include "PayOff.h"
#include "Lattice.h"

using namespace std;

class BinomialMethod {
public:
	BinomialMethod(
		string type,  // "Call" or "Put".
		string ticker, 
		size_t num_of_periods,
		double S,
		double K,
		double discount,  // It changes with different partition.
		double up,
		double down,
		double risk_neutral_q
		); 
	virtual ~BinomialMethod();
	
	void fill_lattice_forward();
	virtual void do_binomial_pricing() = 0;
	void output_results() const;

protected:
	string type_;
	string style_;
	string ticker_;
	size_t num_of_periods_;
	double spot_;
	double strike_price_;
	double discount_;
	double up_;
	double down_;
	double risk_neutral_q_;
	Lattice<double, 2> lattice_;
	unique_ptr<PayOff> payoff_ptr_;
};

// European method, it override do_binomial_pricing(S) and clone() 
// ==============================================================
class EuropeanBinomialMethod : public BinomialMethod
{
public:
	EuropeanBinomialMethod(
		string type,
		string ticker,
		size_t num_of_periods,
		double S,
		double K,
		double discount,
		double up,
		double down,
		double risk_neutral_q
		);
	virtual ~EuropeanBinomialMethod();

	virtual void do_binomial_pricing() override;

};

// American method, it override do_binomial_pricing(S) and clone() 
// ===============================================================
class AmericanBinomialMethod : public BinomialMethod
{
public:
	AmericanBinomialMethod(
		string type,
		string ticker,
		size_t num_of_periods,
		double S,
		double K,
		double discount,
		double up,
		double down,
		double risk_neutral_q
		);
	virtual ~AmericanBinomialMethod();

	virtual void do_binomial_pricing() override;
};

#endif