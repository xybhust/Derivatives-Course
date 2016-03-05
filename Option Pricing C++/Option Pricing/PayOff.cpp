
#include "PayOff.h"

// Definition of base class
// ========================
PayOff::PayOff(double K) : strike_price_(K) { }

PayOff::~PayOff() { }

// Definition of PayOffCall class
// ==============================
PayOffCall::PayOffCall(double K) : PayOff(K) { }

PayOffCall::~PayOffCall() { }

// Definition of PayOffPut class
// =============================
PayOffPut::PayOffPut(double K) : PayOff(K) { }

PayOffPut::~PayOffPut() { }

