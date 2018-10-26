#ifndef EXECUTOR_H_
#define EXECUTOR_H_

#include <algorithm>

#include "Client_Trainer.h"
#include "Client_Scanner.h"
#include "../nonltr/ChromListMaker.h"
#include "../train/Predictor.h"

using namespace std;
using namespace satellites;

namespace satellites {

	class Executor {
	public:
	//remove last two ints after testing filt params
		Executor(string, double, int, string,
			int, string, int, double, int, int, int);
		~Executor();

	private:
		int maxK;
		int minK;


		int champminK;
		int champmaxK;
		int champhalf_win;

		string scanDir;
		double base;
		int min_reg;
		int mtf;
		string outfile;

	string trainFile; 	// The training file
	vector<double> compList;
	int half_win;
	HMM* hmm;

	Predictor<int> * pred;

	// The identify score used by the motif discovery module and filtering
	double idn;

	// Fragment size used while training the mean shift algorithm
	const int trainChromSize = 30000000;
	// // Fragment size used while scanning
	const int scanChromSize = 10000000;

	const int stateNumber = 20;
	// // Fragment size used while training the mean shift algorithm
	// const int trainChromSize = 500000000;
	// // // Fragment size used while scanning
	// const int scanChromSize = 500000000;

	void write_out(
		vector<
		tuple<ILocation*, ChromosomeOneDigit *, string, string,
		double> >&, string&);

	void train();
	void scan();
};

}

#endif
