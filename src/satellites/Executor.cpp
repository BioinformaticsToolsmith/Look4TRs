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
#include "../nonltr/ChromosomeTR.h"
#include "../nonltr/ChromosomeSelfTR.h"
#include "../nonltr/ChromosomeReadTR.h"
#include <unordered_map>

using namespace std;
using namespace satellites;

namespace satellites
{

	Executor::Executor(string scanDirIn, double baseIn,
		int min_regIn, int order_in,
		string outfileIn, int minKIn, int window_loops_in, string trainFileIn,
		int maxKIn, double idnIn, int smoothingWindowIn, int mtfIn,
		bool will_merge_in, int trainingSizeIn,
		bool lng_mtf_in, int minMIn, int maxMIn, int seg_size_in,
		int stateNumberL_in, int stateNumberU_in,
		string chromTR_dir_in,
		string hmm_file_in, string glm_file_in, string chmp_file_in,
		string save_file_in, string bed_motif_file_in, string fa_motif_file_in)
	{

		scanDir = scanDirIn;
		half_win = min_regIn / 2;
		base = baseIn;
		min_reg = min_regIn;
		outfile = outfileIn;
		minK = minKIn;
		trainFile = trainFileIn;
		maxK = maxKIn;
		idn = idnIn;
		smoothingWindow = smoothingWindowIn;
		mtf = mtfIn;
		trainingSize = trainingSizeIn;
		order = order_in;
		minM = minMIn;
		maxM = maxMIn;
		window_loops = window_loops_in;
		will_merge = will_merge_in;
		lng_mtf = lng_mtf_in;
		seg_size = seg_size_in;
		chromTR_dir = chromTR_dir_in;
		stateNumberL = stateNumberL_in;
		stateNumberU = stateNumberU_in;

		compList = vector<double>(0);

		hmm_file = hmm_file_in;
		glm_file = glm_file_in;
		chmp_file = chmp_file_in;
		save_file = save_file_in;
		bed_motif_file = bed_motif_file_in;
		fa_motif_file = fa_motif_file_in;

		train();
		scan();
	}

	Executor::~Executor()
	{
		delete hmm;
		delete pred;
	}

	void Executor::train()
	{
		cout << "Training on " << trainFile << endl;
		ChromListMaker chr_list(trainFile, trainChromSize);
		const vector<nonltr::ChromosomeOneDigit *> *chrs =
		chr_list.makeChromOneDigitList();

	// Make the composition list
		vector<int> *baseCount = chrs->at(0)->getBaseCount();

		if (chrs->at(0)->getEffectiveSize() == 0)
		{
			cerr << "The training chromosome starts with too many Ns" << endl;
			cerr << "Please use another training chromosome" << endl;
			throw std::exception();
		}
		double total = baseCount->at(0) + baseCount->at(1) + baseCount->at(2) + baseCount->at(3);
		for (int bc : *baseCount)
		{
			compList.push_back(log2((double)bc / total));
		}

		if (hmm_file == "" || glm_file == "" || chmp_file == "")
		{
			champminK = -1;
			champmaxK = -1;
			champhalf_win = -1;
			double champFMeaure = -1;
			std::vector<tuple<int, int, int, int>> parameterList;

			parameterList.resize((maxK - minK + 1) * (maxK - minK + 2) / 2 * window_loops * ( (stateNumberU - stateNumberL) / 2 + 1) );

			int index = 0;
			int corNum = (int)Util::CORE_NUM;

			ChromosomeTR * trainingChrom;
			ChromosomeTR * testingChrom;

			if(bed_motif_file == ""){
				trainingChrom = new ChromosomeSelfTR(order, chrs->at(0),
					'N', seg_size, minM, maxM, min_reg, 0);
				testingChrom = new ChromosomeSelfTR(order, chrs->at(0),
					'N', seg_size, minM, maxM, min_reg, 524287);
			}
			else {
				int segmentEnd = -1;
				{
					int segmentSize = 0;
					for (int i = 0; i < chrs->at(0)->getSegment()->size() && segmentSize <= seg_size; i++)
					{
						segmentEnd = chrs->at(0)->getSegment()->at(i)->at(1);
						segmentSize += chrs->at(0)->getSegment()->at(i)->at(1) - chrs->at(0)->getSegment()->at(i)->at(0) + 1;
					}
				}
				if(fa_motif_file == ""){

					trainingChrom = new ChromosomeReadTR(order, chrs->at(0),
						'N', segmentEnd + 1, trainFile, bed_motif_file, 0);
					testingChrom = new ChromosomeReadTR(order, chrs->at(0),
						'N', segmentEnd + 1, trainFile, bed_motif_file, 524287);
				} else {
					ChromosomeOneDigit * read_from_chrom = new ChromosomeOneDigit(fa_motif_file);
					trainingChrom = new ChromosomeReadTR(order, chrs->at(0),
						'N', segmentEnd + 1, trainFile, bed_motif_file, read_from_chrom, 0);
					testingChrom = new ChromosomeReadTR(order, chrs->at(0),
						'N', segmentEnd + 1, trainFile, bed_motif_file, read_from_chrom, 524287);
				}
			}



			if (chromTR_dir.size() != 0)
			{
				trainingChrom->printSequence(chromTR_dir + "/Training.fa");
				trainingChrom->printBedData(chromTR_dir + "/Training.bed");
				testingChrom->printSequence(chromTR_dir + "/Testing.fa");
				testingChrom->printBedData(chromTR_dir + "/Testing.bed");
			}
			

			for (int i = minK; i <= maxK; i++)
			{
				for (int j = i; j <= maxK; j++)
				{
					for (int w = half_win; w <= window_loops * half_win; w += half_win)
					{
						for(int h = stateNumberL; h <= stateNumberU; h += 2){
							std::tuple<int, int, int, int> mytuple(i, j, w, h);
							parameterList.at(index) = mytuple;
							index++;
						}
					}
				}
			}

			#pragma omp parallel for schedule(dynamic) num_threads(Util::CORE_NUM)
			for (int i = 0; i < parameterList.size(); i++)
			{

				int minK = std::get<0>(parameterList.at(i));
				int maxK = std::get<1>(parameterList.at(i));
				int win = std::get<2>(parameterList.at(i));
				int stateNumberLoop = std::get<3>(parameterList.at(i));

				HMM *hmmTest = new HMM(stateNumberLoop);

				Client_Trainer *ctTest = new Client_Trainer(
					trainingChrom, trainingChrom,
					hmmTest,
					std::get<0>(parameterList.at(i)),
					std::get<1>(parameterList.at(i)),
					std::get<2>(parameterList.at(i)), base, min_reg,
					compList);
				double FMeasure = ctTest->getFMeasure();
				#pragma omp critical
				{
					std::cout << "Mink:" << std::get<0>(parameterList.at(i)) << " MaxK:" << std::get<1>(parameterList.at(i)) << " half-wsize:" << std::get<2>(parameterList.at(i))
					<< " hmm-states:" << std::get<3>(parameterList.at(i)) << std::endl;
					std::cout << "Sensitivity:" << ctTest->getSensitivity() << " Precision:" << ctTest->getPrecision() << " FMeasure:" << FMeasure << " ChampionFMeasure:" << champFMeaure << std::endl;
				}
				#pragma omp critical
				{
					if (FMeasure > champFMeaure)
					{
						std::cout << "Found new champion!" << std::endl;
						champminK = std::get<0>(parameterList.at(i));
						champmaxK = std::get<1>(parameterList.at(i));
						champhalf_win = std::get<2>(parameterList.at(i));
						champFMeaure = ctTest->getFMeasure();
						stateNumber = stateNumberLoop;
					}
				}

				delete ctTest;
				delete hmmTest;
			}

			hmm = new HMM(stateNumber);

			std::cout << "The champion parameters are actually " << champminK << " " << champmaxK << " " << champhalf_win << " " << stateNumber << std::endl;
			std::cout << "The K we are using is " << trainingChrom->getK() << std::endl;
			Client_Trainer *ct = new Client_Trainer(trainingChrom,
				testingChrom,
				hmm, champminK, champmaxK,
				champhalf_win, base, min_reg,
				compList);
			hmm->print(save_file + "/hmm.txt");

			{
				ofstream myfile;
				myfile.open(save_file + "/chmp.txt");
				myfile << champminK << std::endl;
				myfile << champmaxK << std::endl;
				myfile << champhalf_win << std::endl;
				myfile.close();
			}

			std::cout << "Training size is " << trainingSize << std::endl;
			pred = new Predictor<int>(trainingChrom->getK(), idn, PRED_MODE_REGR, PRED_FEAT_FAST,
				trainingSize);
			if (mtf)
			{
				
				ct->trainPredictor(pred, idn);
				
				pred->save(save_file + "/glm.txt");
			}

		// Free memory
			delete ct;
		} 
		else
		{
			std::cout << "Reading HMM and GLM" << std::endl;

			pred = new Predictor<int>(glm_file);
			hmm = new HMM(hmm_file);

			string line;
			ifstream myfile(chmp_file);
			if (myfile.is_open())
			{
				getline(myfile, line);
				champminK = std::stoi(line);
				getline(myfile, line);
				champmaxK = std::stoi(line);
				getline(myfile, line);
				champhalf_win = std::stoi(line);
				myfile.close();
			}
			else
			{
				std::cerr << "Unable to open file" << std::endl;
			}
		}
	}

/**
   * Scanning all files in a directory
   */
	void Executor::scan()
	{
		vector<string> fileList;
		Util::readChromList(scanDir, &fileList, string("fa"));

		int corNum = min((int)Util::CORE_NUM, (int)fileList.size());
		#pragma omp parallel for schedule(dynamic) num_threads(corNum)
		for (int k = 0; k < fileList.size(); k++)
		{

			string addr = fileList.at(k);

			#pragma omp critical
			{
				cout << "Scanning " << addr << endl;
			}

			const char *file_name = addr.c_str();
			string str = outfile + "/" + string(basename((char *)file_name));
			str = str.substr(0, str.size() - 2) + "bed";

			ChromListMaker chromListMaker(addr, scanChromSize);
			const vector<nonltr::ChromosomeOneDigit *> *chromList =
			chromListMaker.makeChromOneDigitList();

			HMM *copyHMM = new HMM(*hmm);

			Client_Scanner *cs;

		// #pragma omp ordered
			cs = new Client_Scanner(chromList, copyHMM, compList,
				champminK, champmaxK, champhalf_win, base, idn, smoothingWindow, mtf, pred,
				min_reg, will_merge);

		// Get results
			vector<tuple<ILocation *, ChromosomeOneDigit *, string, string, double>> good_sats(
				0);
			cs->get_hmm_sats(good_sats);

		/*
       * Fix the coordinates here
       * For each item in good_sat, get the start and add the start from ChromLIstMakker::getSplitList()
       * It's a vector of vectors
       * Each item corresponds to a chromosome and each vector corresponds to the regions/bases
       */
			for (int i = 0; i < good_sats.size(); i++)
			{
				ChromosomeOneDigit *oneDigit = std::get<1>(good_sats[i]);
				int newStart =
				(chromListMaker.getStartOfChromosome(oneDigit)).second;
				ILocation *oldLocation = std::get<0>(good_sats[i]);
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
		vector<tuple<ILocation *, ChromosomeOneDigit *, string, string, double>> &regs,
		string &addr)
	{
#pragma omp critical
		{
			cout << "Writing to " << addr << endl;
		}

		ofstream output_p;

		output_p.open((addr).c_str(), fstream::out);
		if (!output_p.good())
		{
			cout << "Cannot write to " << addr << endl;
		}

	// The end is exclusive
		for (auto &reg : regs)
		{
		// This length filter MAY CAUSE PROBLEMS for the mini or the full
			if (std::get<0>(reg)->getLength() > (min_reg / 2) && (!mtf || std::get<4>(reg) >= idn))
			{
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

} // namespace satellites
