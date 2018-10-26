/**
 * Contains the shared code between the trainer and the scanner
 *
 */

#ifndef ICLIENT_H_
#define ICLIENT_H_

#include <vector>
#include "ScorerSat.h"
#include "../nonltr/ChromosomeOneDigit.h"
#include "../nonltr/HMM.h"

using namespace std;
using namespace satellites;

namespace satellites {

class IClient {
public:
	IClient(HMM *, vector<double>&, int, int, int, double);
	virtual ~IClient();
	HMM * getHMM();
	void decode(ChromosomeOneDigit *, ScorerSat *,vector<ILocation*>*);

protected:
	// Used by the scorers
	int minK;
	// Used by the scorers
	int maxK;
	// Used by the scorers
	int halfW;
	// Used by the HMM
	double base;
	// The HMM
	HMM * hmm;

	vector<double>& compList;
	vector<ScorerAdjusted*> * tableList;
	ScorerSat * makeScorer(ChromosomeOneDigit *);
};

}
#endif
