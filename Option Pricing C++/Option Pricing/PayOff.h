
#ifndef PAYOFF_H
#define PAYOFF_H
#include <algorithm>

// =======================================================
// Abstract base class 
// Override () operator to make it behave like a function
// =======================================================
class PayOff {
public:
	explicit PayOff(double K);
	virtual ~PayOff();

	virtual double calculate_payoff(double S) const = 0; 

protected:
	double strike_price_;
};

// Call option PayOff
// ==================
class PayOffCall : public PayOff {
public:
	explicit PayOffCall(double K);
	virtual ~PayOffCall();

	virtual double calculate_payoff(double S) const override {
		return (std::max)(0.0, S - strike_price_);
	}
};

// Put option PayOff
// =================
class PayOffPut : public PayOff {
public:
	explicit PayOffPut(double K);
	virtual ~PayOffPut();

	virtual double calculate_payoff(double S) const override {
		return (std::max)(0.0, strike_price_ - S);
	}
};

#endif