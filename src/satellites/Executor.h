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
		Executor(string, double, int, int, string,
			int, int, string, int, double, int, int,
			bool, int, bool, int, int, int,
			int, int,
			 string = "",
			string = "", string = "", string = "",
			 string = "./", string = "", string = "");
		~Executor();

	private:
		int maxK;
		int minK;

		int minM;
		int maxM;

		int order;

		int window_loops;
		int champminK;
		int champmaxK;
		int champhalf_win;

		string scanDir;
		double base;
		int min_reg;
		int mtf;
		bool will_merge;
		bool lng_mtf;
		int trainingSize;
		string outfile;

		string trainFile; 	// The training file
		vector<double> compList;
		int half_win;
		HMM* hmm;

		string hmm_file;
		string glm_file;
		string chmp_file;
		string save_file;
		string chromTR_dir;
		string bed_motif_file;
		string fa_motif_file;
		int seg_size;

		Predictor<int> * pred;

		// The identify score used by the motif discovery module and filtering
		double idn;
		// The smoothing window used for finding peaks in the mini and full
		int smoothingWindow;

		// Fragment size used while training the mean shift algorithm
		const int trainChromSize = 30000000;
		// // Fragment size used while scanning
		const int scanChromSize = 10000000;

		int stateNumber;
		int stateNumberL;
		int stateNumberU;

		void write_out(
			vector<
			tuple<ILocation*, ChromosomeOneDigit *, string, string,
			double> >&, string&);

		void fillCompList();
		void train();
		void scan();
	};

}

#endif
