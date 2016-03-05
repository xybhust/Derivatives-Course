
#ifndef LATTICE_H
#define LATTICE_H
#include <vector>
#include <tuple>
#include <string>
#include <iostream>

using namespace std;

/* =============================================================================
 This is a template class for non-path-dependent option pricing.
 NumberNodes = 2 represents binomial, 3 trinomial
 Data member:
 tree_ - the outler vector is the horizontal periods while the inner is the 
         vertical nodes at each period. 

 Index starting from 0, from left to the right, from bottom to the top
 Each node stores a tuple of three elements:
0: stock price, 1: option price, 2: "Yes" or "No" indicate early exercise.
=============================================================================*/
template <typename T, int NumberNodes> 
class Lattice {
public:
	Lattice() = default;
	explicit Lattice(size_t num_of_periods); 
	Lattice(const Lattice<T, NumberNodes>& source);
	virtual ~Lattice();

	Lattice<T, NumberNodes>& operator = (const Lattice<T, NumberNodes>& source);
	
	vector<tuple<T, T, string>>& operator [] (size_t n);
	const vector<tuple<T, T, string>>& operator [] (size_t n) const;

	size_t min_index() const;
	size_t max_index() const;
	size_t size() const;
	size_t num_of_periods() const;
private:
	vector<vector<tuple<T, T, string>>> tree_;

};
template <typename T, int NumberNodes>
inline vector<tuple<T, T, string>>& Lattice<T, NumberNodes>::operator [] (size_t n) { 
	return tree_[n]; 
}


template <typename T, int NumberNodes>
inline const vector<tuple<T, T, string>>& Lattice<T, NumberNodes>::operator [] (size_t n) 
const { 
	return tree_[n]; 
}

#include "Lattice.cpp"
#endif

