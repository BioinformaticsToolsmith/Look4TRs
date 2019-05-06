/**
 * Author: Hani Z. Girgis, PhD, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 * Based on code by Vincent Wells, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 *
 * This class uses the mean shift to train the HMM.
 */

#ifndef CLIENT_TRAINER_H_
#define CLIENT_TRAINER_H_

#include "../train/Predictor.h"
#include "IClient.h"
#include "../nonltr/ChromosomeOneDigit.h"
#include "../nonltr/ChromosomeRandom.h"
#include "../nonltr/ChromosomeTR.h"
#include "../motif/FindMotif.h"

#include <algorithm>
#include <chrono>
#include <ctime>

using namespace std;
using namespace satellites;

namespace satellites {
	class Client_Trainer: public IClient {
	public:
		Client_Trainer(ChromosomeTR*, ChromosomeTR*, HMM*, int, int, int, double, int,
			vector<double>&);
		virtual ~Client_Trainer();

		vector<ILocation*>* getChromSats();

		double getPrecision();
		double getSensitivity();
		double getFMeasure();
		void trainPredictor(Predictor<int> *, int);
	private:
		int champminK;
		int champmaxK;
		ChromosomeTR * trainChrom;
		ChromosomeTR * testChrom;
		ScorerSat * train_scorer;
		ScorerSat * test_scorer;
		vector<ILocation*>* chromSats;

		void train();
	};
}
#endif
