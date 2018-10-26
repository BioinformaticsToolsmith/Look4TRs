/*
 * ScorerAdjusted.h
 *
 *  Created on: Jan 8, 2017
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef SATELLITES_SCORERADJUSTED_H_
#define SATELLITES_SCORERADJUSTED_H_

#include <string>
#include <vector>
#include "../utility/HashMaker.h"
#include "../nonltr/KmerHashTable.h"
#include "../exception/InvalidStateException.h"
#include "../exception/InvalidInputException.h"

using namespace std;
using namespace exception;

namespace nonltr {

class ScorerAdjusted {

private:
	int k;
	int halfW;
	int center;
	const string* seq;
	// Initial value of each entry in the table
	const int initValue = 0;
	KmerHashTable<unsigned long, int> * table;
	HashMaker<unsigned long> * hasher;
	vector<unsigned long> * hashList;
	vector<double>& compList;
	vector<double> * expectedScoreList;
	int lastIndexInHashList;
	void fillExpectedScoreList(int, int);
	inline double round(double number) {
		return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
	}

public:
	ScorerAdjusted(int, int, vector<double>&);
	virtual ~ScorerAdjusted();
	// Start processing a whole segment of the same sequence.
	// And return the score of nucleotide at zero.
	int processSegment(const string*, int, int);
	// Delete the first nucleotide and add the last nucleotide.
	// And return the score of the nucleotide at the center.
	int moveOneNucleotide();

	// Get the raw score
	int getScoreOfCenter();
	// Get the adjusted score
	int getAdjustedScoreOfCenter();

	// Nothing is stored in the table after calling this method
	void clear();

};

} /* namespace nonltr */

#endif /* SATELLITES_SCORERADJUSTED_H_ */
