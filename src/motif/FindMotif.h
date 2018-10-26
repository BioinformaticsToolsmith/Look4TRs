// Author: Alfredo Velasco II
// The Bioinformatics Toolsmith Laboratory, the University of Tulsa

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <limits>

#include "../align/GlobAlignE.h"
#include "../nonltr/ChromosomeOneDigit.h"
#include "../train/Predictor.h"
#include "../align/RepeatAlignE.h"

#ifndef FINDMOTIF_H_
#define FINDMOTIF_H_

namespace motif {

class FindMotif {
private:
	string sequence;
	string foundMotif;
	bool isFound;
	double identityScore;
	double threshold;
	const int MICRO_MAX_SIZE = 10;
	const int TOP = 10;
	string h1;
	Predictor<int> * pred;
	Point<int> * seqPoint;

	string dequeToString(deque<char>&);
	void greedyConfirmation(vector<string> *);
	void searchMicro();

public:
	FindMotif(string, double, Predictor<int> *);
	virtual ~FindMotif();
	static string makeExact(string, int);

	bool getIsFound();
	string getFoundMotif();
	double getIdentityScore();
};
}
#endif
