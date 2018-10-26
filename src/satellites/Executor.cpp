/**
 * Author: Vincent Wells, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 *
 * Modified and refactored by Hani Z. Girgis, PhD, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 * Modified by Robert Geraghty, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 * Coordinate the training and the scanning processes
 * Modified by Alfredo Velasco, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 * Changed the parameters so that Executor only needs a training file, not a directory
 */

#include <libgen.h>
#include "Executor.h"
#include "../utility/Util.h"
#include "../cluster/Progress.h"
#include <unordered_map>

 using namespace std;
 using namespace satellites;



 namespace satellites {

 	Executor::Executor(string scanDirIn, double baseIn,
 		int min_regIn,
 		string outfileIn, int minKIn, string trainFileIn,
 		int maxKIn, double idnIn, int mtfIn,
 		int kmerSize, int trainingSize) {


 		scanDir = scanDirIn;
 		half_win = min_regIn / 2;
 		base = baseIn;
 		min_reg = min_regIn;
 		outfile = outfileIn;
 		minK = minKIn;
 		trainFile = trainFileIn;
 		maxK = maxKIn;
 		idn = idnIn;
 		mtf = mtfIn;

 		compList = vector<double>(0);


 		hmm = new HMM(1.000, stateNumber); // The first parameter seems to be unused

 		pred = new Predictor<int>(kmerSize, idn, PRED_MODE_REGR, PRED_FEAT_FAST,
 			trainingSize);

 		train();
 		scan();
 	}

 	Executor::~Executor() {
 		delete hmm;
 		delete pred;
 	}

// ToDo: Handle the case when the chromosome is very short
 	void Executor::train() {
 		cout << "Training on " << trainFile << endl;
 		ChromListMaker chr_list(trainFile, trainChromSize);
 		const vector<nonltr::ChromosomeOneDigit*>* chrs =
 		chr_list.makeChromOneDigitList();

		// Make the composition list
 		vector<int> * baseCount = chrs->at(0)->getBaseCount();

 		if (chrs->at(0)->getEffectiveSize() == 0) {
 			cerr << "The training chromosome starts with too many Ns" << endl;
 			cerr << "Please use another training chromosome" << endl;
 			throw std::exception();
 		}
 		double total = baseCount->at(0) + baseCount->at(1) + baseCount->at(2)
 		+ baseCount->at(3);
 		for (int bc : *baseCount) {
 			compList.push_back(log2((double) bc / total));
 		}

 		const int window_loops = 4;

 		champminK = -1;
 		champmaxK = -1;
 		champhalf_win = -1;
 		double champFMeaure = -1;
 		// Progress bar((maxK - minK + 1) * (maxK - minK + 2) / 2 * 30, "Optimizing " );
 		std::vector<tuple<int, int, int>> parameterList;
 		parameterList.resize((maxK - minK + 1) * (maxK - minK + 2) / 2 * window_loops);
 		int segmentEnd = -1;
 		{
 			int segmentSize = 0;
 			for (int i = 0; i < chrs->at(0)->getSegment()->size() && segmentSize <= 3000000; i++) {
 				segmentEnd = chrs->at(0)->getSegment()->at(i)->at(1);
 				segmentSize += chrs->at(0)->getSegment()->at(i)->at(1)
 				- chrs->at(0)->getSegment()->at(i)->at(0) + 1;
 			}
 		}


 		int index = 0;
 		int corNum = (int) Util::CORE_NUM;

 		// Micro
 		ChromosomeTR * shuffledChrom = 
 		new ChromosomeTR(0, chrs->at(0), 
 			'N', segmentEnd + 1, 1, 10, min_reg);
 		// Full
 		// ChromosomeTR * shuffledChrom = 
 		// new ChromosomeTR(max((minK + maxK) / 4 - 2, 0), chrs->at(0), 
 		// 	'N', segmentEnd + 1, 10, 100, min_reg);
 		// ChromosomeTR * shuffledChrom = 
 		// new ChromosomeTR(0, chrs->at(0), 
 		// 	'N', segmentEnd + 1, 10, 100, min_reg);
 		// ChromosomeTR * shuffledChrom = 
 		// new ChromosomeTR(0, chrs->at(0), 
 		// 	'N', segmentEnd + 1, 10, 20, min_reg);

 		for(int i = minK; i <= maxK; i++){
 			for(int j = i; j <= maxK; j++){
 				for(int w = half_win; w <= window_loops * half_win; w += half_win){
 					std::tuple<int,int, int> mytuple (i, j, w);
 					parameterList.at(index) = mytuple;
 					index++;
 				}
 			}
 		}

 		#pragma omp parallel for schedule(dynamic)
 		for(int i = 0; i < parameterList.size(); i++){

 			HMM * hmmTest = new HMM(1.000, stateNumber);
 				// Client_Trainer * ctTest = new Client_Trainer(chrs->at(0), hmmTest, i, j,
 				// 	half_win, base, max_reg, min_reg, itr_1D, itr_2D, win_2D, kernel,
 				// 	compList, st_spc_1D);
 			int minK = std::get<0>(parameterList.at(i));
 			int maxK = std::get<1>(parameterList.at(i));
 			int win  = std::get<2>(parameterList.at(i));

 			Client_Trainer * ctTest = new Client_Trainer(
 				shuffledChrom, 
 				hmmTest,
 				std::get<0>(parameterList.at(i)), 
 				std::get<1>(parameterList.at(i)),
 				std::get<2>(parameterList.at(i))
 				, base, min_reg,
 				compList);
 			double FMeasure = ctTest->getFMeasure();
 			#pragma omp critical
 			{
 				std::cout << "Mink:" << std::get<0>(parameterList.at(i)) << " MaxK:" << std::get<1>(parameterList.at(i)) << " half-wsize:" << std::get<2>(parameterList.at(i)) << std::endl;
 				std::cout << "Sensitivity:"<< ctTest->getSensitivity() << " Precision:" << ctTest->getPrecision() <<  " FMeasure:" << FMeasure << " ChampionFMeasure:" << champFMeaure << std::endl;
 			}
 			#pragma omp critical
 			{
 				if(FMeasure > champFMeaure){
 					std::cout << "Found new champion!" << std::endl;
 					champminK = std::get<0>(parameterList.at(i));
 					champmaxK = std::get<1>(parameterList.at(i));
 					champhalf_win = std::get<2>(parameterList.at(i));
 					champFMeaure = ctTest->getFMeasure();
 				}
 				// bar++;
 			}
 				// for(int i = 0; i < chromSats->size(); i++){
 				// 	std::cout << chromSats->at(i)->getStart() << " " << chromSats->at(i)->getEnd() << std::endl;
 				// }

 			delete ctTest;
 			delete hmmTest;
 		}


	// Train
 		// HMM * hmmTest = new HMM(1.000, stateNumber);
 		// Client_Trainer* ctTrim = new Client_Trainer(shuffledChrom, hmmTest, champminK, champmaxK,
 		// 	champhalf_win, base, min_reg,
 		// 	compList);

 		// std::vector<ILocation *> * satList = ctTrim->getChromSats();
 		// std::vector<Location *> * falsePositiveList = Util::wholeLocationSubtract(ctTrim->getChromSats(), shuffledChrom->getRegionList());

 		// shuffledChrom->removeLocations(falsePositiveList);


 		// Util::deleteInVector(falsePositiveList);
 		// delete falsePositiveList;
 		// delete ctTrim;


 		std::cout << "The champion parameters are actually " << champminK << " " << champmaxK << " " << champhalf_win << std::endl;
 		Client_Trainer* ct = new Client_Trainer(shuffledChrom, hmm, champminK, champmaxK,
 			champhalf_win, base, min_reg,
 			compList);

 		
 		if(mtf){
 			ct->trainPredictor(pred, idn);
 		}

	// Free memory
 		delete ct;
 	}

/**
 * Scanning all files in a directory
 */
 void Executor::scan() {
	// cout << "Scanning ... " << endl;
 	vector<string> fileList;
 	Util::readChromList(scanDir, &fileList, string("fa"));

 	int corNum = min((int) Util::CORE_NUM, (int) fileList.size());
# pragma omp parallel for schedule(dynamic) num_threads(corNum) //corNum
 // # pragma omp parallel for num_threads(corNum) //corNum
 	// # pragma omp parallel for ordered num_threads(corNum) //corNum
 	for (int k = 0; k < fileList.size(); k++) {

 		string addr = fileList.at(k);

		# pragma omp critical
 		{
 			cout << "Scanning " << addr << endl;
 		}

 		const char* file_name = addr.c_str();
 		string str = outfile + "/" + string(basename((char*) file_name));
 		str = str.substr(0, str.size() - 2) + "bed";

 		ChromListMaker chromListMaker(addr, scanChromSize);
 		const vector<nonltr::ChromosomeOneDigit*>* chromList =
 		chromListMaker.makeChromOneDigitList();

 		HMM * copyHMM = new HMM(*hmm);

 		Client_Scanner * cs;

		// #pragma omp ordered
 		cs = new Client_Scanner(chromList, copyHMM, compList,
 			champminK, champmaxK, champhalf_win, base, idn, mtf, pred,
 			min_reg);


		// Get results
 		vector<tuple<ILocation*, ChromosomeOneDigit *, string, string, double> > good_sats(
 			0);
 		cs->get_hmm_sats(good_sats);

		/*
		 * Fix the coordinates here
		 * For each item in good_sat, get the start and add the start from ChromLIstMakker::getSplitList()
		 * It's a vector of vectors
		 * Each item corresponds to a chromosome and each vector corresponds to the regions/bases
		 */
		 for (int i = 0; i < good_sats.size(); i++) {
		 	ChromosomeOneDigit * oneDigit = std::get<1>(good_sats[i]);
		 	int newStart =
		 	(chromListMaker.getStartOfChromosome(oneDigit)).second;
		 	ILocation * oldLocation = std::get<0>(good_sats[i]);
		 	oldLocation->setEnd(oldLocation->getEnd() + newStart);
		 	oldLocation->setStart(oldLocation->getStart() + newStart);
		 }

		 write_out(good_sats, str);

		// Clean up
		 delete copyHMM;
		 delete cs;
		}

		fileList.clear();
	}

	void Executor::write_out(
		vector<tuple<ILocation*, ChromosomeOneDigit *, string, string, double> >& regs,
		string& addr) {
#pragma omp critical
		{
			cout << "Writing to " << addr << endl;
		}

		ofstream output_p;

		output_p.open((addr).c_str(), fstream::out);
		if (!output_p.good()) {
			cout << "Cannot write to " << addr << endl;
		}

// The end is exclusive
		for (auto& reg : regs) {
		// This length filter MAY CAUSE PROBLEMS for the mini or the full
			if (std::get<0>(reg)->getLength() > (min_reg / 2)
				&& (!mtf || std::get<4>(reg) >= idn)) {
				string header = std::get<1>(reg)->getHeader();
			header = header.substr(1);
			replace(header.begin(), header.end(), ' ', '_');

			output_p << header << "\t";
			output_p << std::get<0>(reg)->getStart() << "\t";
			output_p << std::get<0>(reg)->getEnd() + 1 << "\t";
			output_p << std::get<2>(reg) << "\t";
			output_p << std::get<3>(reg) << "\t";
			output_p << std::get<4>(reg) << endl;
		}
	}
	output_p.close();
}

}
