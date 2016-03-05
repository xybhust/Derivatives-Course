
#ifndef LATTICE_CPP
#define LATTICE_CPP

#include "Lattice.h"
#include <algorithm>

template <typename T, int NumberNodes>
Lattice<T, NumberNodes>::Lattice(size_t num_of_periods) : 
tree_(num_of_periods + 1) {
	tree_[0].resize(1); // There is only one root
	if (num_of_periods > 0)	{
		unsigned long length = tree_.size();
		unsigned long index = 1; // start from 1
		std::for_each(tree_.begin() + 1, tree_.end(), 
			[&index](vector<tuple<T, T, string>>& x) -> void {
			x.resize(index + NumberNodes - 1); index += NumberNodes - 1; }
		);
	}
}

template <typename T, int NumberNodes>
Lattice<T, NumberNodes>::Lattice(const Lattice<T, NumberNodes>& source) :
tree_(source.tree_) { }

template <typename T, int NumberNodes>
Lattice<T, NumberNodes>::~Lattice() { }

template <typename T, int NumberNodes>
Lattice<T, NumberNodes>& Lattice<T, NumberNodes>::operator = (
	const Lattice<T, NumberNodes>& source) {
	if (this != &source)
		tree_ = source.tree_;
	return *this;
}

template <typename T, int NumberNodes>
size_t Lattice<T, NumberNodes>::min_index() const { return 0; }

template <typename T, int NumberNodes>
size_t Lattice<T, NumberNodes>::max_index() const { return tree_.size() - 1; }

template <typename T, int NumberNodes>
size_t Lattice<T, NumberNodes>::size() const { return tree_.size(); }

template <typename T, int NumberNodes>
size_t Lattice<T, NumberNodes>::num_of_periods() const { 
	return tree_.size() - 1; 
}


#endif