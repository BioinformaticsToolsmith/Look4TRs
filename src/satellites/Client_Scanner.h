#ifndef CLIENT_SCANNER_H_
#define CLIENT_SCANNER_H_

#include <stdio.h>
#include <stdlib.h>

#include "IClient.h"
#include "../motif/FindMotif.h"
#include "../satellites/ScorerSat.h"
#include "../utility/ILocation.h"
#include "../train/Predictor.h"

using namespace std;
using namespace satellites;
using namespace motif;

namespace satellites {
class Client_Scanner: public IClient {

private:
	const vector<ChromosomeOneDigit*>* chromList;
	vector<vector<ILocation*>*> * hmm_sats;
	string oneDigitToNucleotide(const string *, int, int);

	// The identify score used by the motif discovery module and filtering
	double idn;

	// The smoothing window used for finding peaks in the mini and full
	int smoothingWindow;

	int mtf;

	Predictor<int> * pred;
	int minReg;
	int will_merge;

	void decode(ChromosomeOneDigit *, ScorerSat *, vector<ILocation*>*);

public:
	Client_Scanner(const vector<ChromosomeOneDigit*>*, HMM*, vector<double>&,
			int, int, int, double, double, int, int, Predictor<int> *, int, bool);
	virtual ~Client_Scanner();

	void get_hmm_sats(
			vector<
					tuple<ILocation*, ChromosomeOneDigit *, string, string,
							double> >&);
};
}
#endif
