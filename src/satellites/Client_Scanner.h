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


	int mtf;

	Predictor<int> * pred;
	int minReg;

	void decode(ChromosomeOneDigit *, ScorerSat *, vector<ILocation*>*);
	void decode(vector<vector<utility::ILocation*>*>&, vector<ScorerSat*>*,
			const vector<ChromosomeOneDigit*>*);

public:
	Client_Scanner(const vector<ChromosomeOneDigit*>*, HMM*, vector<double>&,
			int, int, int, double, double, int, Predictor<int> *, int);
	virtual ~Client_Scanner();

	// ToDo: convert tuple to a class
	void get_hmm_sats(
			vector<
					tuple<ILocation*, ChromosomeOneDigit *, string, string,
							double> >&);
};
}
#endif
