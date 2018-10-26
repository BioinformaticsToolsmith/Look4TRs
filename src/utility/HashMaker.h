/*
 * HashMaker.h
 *
 *  Created on: Jun 26, 2015
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef SRC_UTILITY_HASHMAKER_H_
#define SRC_UTILITY_HASHMAKER_H_

#include <vector>
#include <math.h>
#include "../exception/InvalidInputException.h"
#include "../exception/InvalidStateException.h"
#include "../utility/Util.h"

using namespace std;
using namespace exception;

namespace utility {

template<class I>
class HashMaker {
private:
	int k;
	// [4^0, 4^1, ... , 4^(k-1)]
	I * bases;
	I * mMinusOne;

public:
	HashMaker(int);
	virtual ~HashMaker();
	I hash(const char *);
	I hash(const char *, int);
	void hash(const char *, int, int, vector<I> *);

};

} /* namespace utility */

#include "HashMaker.cpp"

#endif /* SRC_UTILITY_HASHMAKER_H_ */
