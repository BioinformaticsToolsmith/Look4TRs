/*
 * ScorerSat.h
 *
 *  Created on: Jan 6, 2017
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef SATELLITES_SCORERSAT_H_
#define SATELLITES_SCORERSAT_H_

#include <iostream>
#include <vector>
#include <numeric>

#include "ScorerAdjusted.h"
#include "../nonltr/ChromosomeOneDigit.h"
#include "../utility/Util.h"

using namespace std;
using namespace nonltr;
using namespace utility;

namespace satellites {

class ScorerSat {

private:
	ChromosomeOneDigit& chrom;
	bool flatten;
	double base;
	int minK;
	int maxK;
	int halfW;
	void processSegment(int, int);
	void make_flattened();
	vector<double>& compList;
	vector<ScorerAdjusted *> * scorerList;
	vector<int> * adjustedList;
	vector<char> * bestKList;
	vector<int> * adjustedList_flattened;
	vector<int> * adjustedList_zeroed;
	bool canDeleteZeroed;

	inline double round(double number) {
		return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
	}

public:
	ScorerSat(ChromosomeOneDigit&, int, int, int, vector<double>&,
			vector<ScorerAdjusted *> *, double);
	virtual ~ScorerSat();
	vector<int>* getScores() const;
	vector<int>* getFlatScores() const;
	vector<char>* getBestKList() const;
	void make_zeroed(double);
	vector<int>* getZeroScores() const;
	ChromosomeOneDigit* get_chrom() const;
	vector<double>* get_compList() const;

	int getMax();
	double standardDeviation(double);
	double avgScore();

	void printScores(string, bool);
	const bool getCanDeleteZeroed();
};

}

#endif /* SATELLITES_SCORERSAT_H_ */
