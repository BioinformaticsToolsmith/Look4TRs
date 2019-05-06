/**
 * Author: Vincent Wells, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 *
 * Edited by Hani Z. Girgis, PhD, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 * Edited by Robert Geraghty, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 * Edited by Alfredo Velasco, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 *
 * The driver function
 */

#include "Executor.h"
#include "../train/Predictor.h"
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <map>
#include <algorithm>

using namespace std;

int argCount = 35;

bool has_only_digits(const string s)
{
	return s.find_first_not_of("-.0123456789") == string::npos;
}

void varSetter(pair<string, bool> &val, int &var, string varFlag)
{

	if (has_only_digits(val.first.c_str()))
	{
		var = atoi(val.first.c_str());
	}
	else
	{
		cerr << "Error: \"" << varFlag << "\" must be an integer! You passed: "
		<< val.first.c_str() << endl;
		exit(1);
	}
}

void printHelp()
{
	string top =
	"[#]------------------------------------------------------------------------------[#]";
	cout
	<< "||==================================LOOK4TRS-HELP=================================||"
	<< endl
	<< endl;
	cout
	<< "- - - - - - - - - - - - - -ACCEPTABLE COMMAND FORMATS- - - - - - - - - - - - - -"
	<< endl;
	cout << "General Command Format" << endl;
	cout << "   ./Look4TRs --help" << endl;
	cout << "Command Format Using Default Option Flags" << endl;
	cout << "   ./Look4TRs --adr Input --out Out --default" << endl;
	cout << "Command Format Using the Parameter Option Flags" << endl;
	cout << "   ./Look4TRs [required options][(all)parameter options]"
	<< endl
	<< endl;
	cout
	<< "- - - - - - - - - - - - - - - INFORMATIONAL OPTIONS- - - - - - - - - - - - - - -"
	<< endl;
	cout << "--help or -h" << endl;
	cout << "   Prints out useful information on using the tool." << endl
	<< endl;
	cout
	<< "- - - - - - - - - - - - - - - - REQUIRED OPTIONS - - - - - - - - - - - - - - - -"
	<< endl;
	cout << "--adr <string>" << endl;
	cout << "    Takes the address of the directory of the sequences to be\n"
	<< "    scanned. Each file must be a FASTA file and end with .fa" << endl;

	cout << "--out <string>" << endl;
	cout << "   Takes the address of the output prediction directory." << endl
	<< endl;
	cout << top << endl;
	cout
	<< " | NOTE: After supplying the required options, either provide a single default    | "
	<< endl
	<< " | option, or provide all of the parameter options!                               | "
	<< endl;
	cout << top << endl
	<< endl;
	cout
	<< "- - - - - - - - - - - - - - - - -DEFAULT OPTIONS - - - - - - - - - - - - - - - -"
	<< endl;
	cout
	<< "--default                                                                         "
	<< endl;
	cout
	<< "    Runs the tool with the default parameters."
	<< endl;
	cout
	<< endl;
	cout
	<< "- - - - - - - - - - - - - - - - PARAMETER OPTIONS- - - - - - - - - - - - - - - -"
	<< endl;
	cout << "--trn <string>" << endl;
	cout << "    Takes the address of the directory of the training file."
	<< endl;
	cout
	<< "--motif-location-file <file.bed>"
	<< endl;
	cout
	<< "   This option will make Look4TRs train using repeats that are specified by a bed file." << endl
	<< "   Doing this will mean that minm and maxm will not have to be initialized." << std::endl;
	cout <<
	"   WARNING: This file should correspond with the file used for training Look4TRs or the file indicated by --motif-fa-file!"
	<< endl;
	cout << "--motif-fa-file <file.fa>"
	<< endl;
	cout <<
	"   This will determine the fa file used to extract repeats for training." << endl <<
	"   If ommited, the training file (either specifed by the --trn flag or the largest one in the --adr directory) will be used." <<
	endl;
	cout <<
	"   If included, a bed file (initialized by --motif-fa-file) must be included with it."
	<< endl; 
	cout <<
	"   WARNING: This file should correspond with the file indicated by --motif-location-file!"
	<< endl << endl;

	cout
	<< "--min <integer>                                                                 "
	<< endl;
	cout << "    The minimum k-mer length analyzed (Must be less than or equal" << endl;
	cout << "   to the maximum k-mer length and greater than 0)." << endl;
	cout << "   This is used by the scoring module." << endl;
	cout
	<< "--max <integer>                                                                 "
	<< endl;
	cout
	<< "    The maximum k-mer length analyzed (Must be greater than or equal" << endl;
	cout << "   to the minimum k-mer length and greater than 0)." << endl;
	cout << "   This is used by the scoring module." << endl;
	cout
	<< "--int <integer>                                                                 "
	<< endl;
	cout << "   The minimum allowed region size for predictions. (Must be greater than 0)" <<  endl;
	cout << "--win <integer>" 
	<< endl;
	cout << "   The number of multiples of --int to use. (Must be positive)" << endl;
	cout << "   This will allow Look4TRs to auto-callibrate the scoring module." << endl;
	cout << "--hmm-states-upper <integer>" << endl;
	cout << "   the upper bound for the number of hmm states used when auto-callibrating Look4TRs." << endl;
	cout << "   (must be positive, even, and higher than or equal to --hmm-states-lower)" << endl;
	cout << "--hmm-states-lower <integer>" << endl;
	cout << "   The lower bound for the number of HMM states used when auto-callibrating Look4TRs." << endl;
	cout << "   (Must be postive, even, and lower than or equal to --hmm-states-upper)" <<  endl;
	cout << endl << endl;

	cout << "--seg <integer>"
	<< endl;
	cout << "   The size of the training chromosome" << endl;
	cout << "   (Must be positive)" << endl;
	cout << "--minm <integer>                                                                 "
	<< endl;
	cout << "   The minimum size of the motifs used to generate the training Chromosome." << endl;
	cout << "   (Must be at least 1 and less than or equal to maxm)" << endl;

	cout << "--maxm <interger>                                                                 " << endl;
	cout << "   The maximum size of the motifs used to generate the training Chromosome." << endl;
	cout << "   (Must be at least 1 and bigger than or equal to minm)" << endl;
	cout << "--ord <integer>"
	<< endl;
	cout << "   The markov order used to generate the training chromosome."
	<< endl;
	cout << "   (Must not be negative)"
	<< endl << endl;
	
	cout
	<< "--mtf <0|1>                                                                 "
	<< endl;
	cout << "    Enable ('1') or disable ('0') the motif discovery feature."
	<< endl;
	cout << "    If disabled ('0'), the locations of the repeats will be outputted, but Look4TRs will not attempt to discover " << endl;
	cout << "    the motif that comprises the repeat." << endl;
	cout << "--prn <integer>" << endl;
	cout
	<< "   The size of the training set of the motif analyzer. Must be positive."
	<< endl;
	cout
	<< "--idn <integer>                                                                 "
	<< endl;
	cout << "   The minimum allowed identity score between the exact repeat and a candidate"
	<< endl;
	cout << "   (Must be 0 or greater)" << endl << endl;

	cout << endl;
	cout << "--save-data <directory>"
	<< endl;
	cout << "   A directory where the --hmm, --glm, and --chmp files will be stored."
	<< endl;
	cout << "--hmm <file>" << endl;
	cout << "   The hmm file that is produced from a previous run that used --save-data." << endl;
	cout << "   This HMM will be loaded, instead of training one."
	<< endl;
	cout << "   Note: If --hmm, --glm, or --chmp is used, all of them must be used!" << endl;
	
	cout << "--glm <file>" << endl;
	cout << "   The glm file that is produced from a previous run that used --save-data." << endl;
	cout << "   This GLM will be loaded, instead of training one."
	<< endl;
	cout << "   Note: If --hmm, --glm, or --chmp is used, all of them must be used!" << endl;
	
	cout << "--chmp <file>" << endl;
	cout << "   The file that is produced from a previous run that used --save-data." << endl;
	cout << "   This will be used for parameters used from training the HMM and GLM."
	<< endl;
	cout << "   Note: If --hmm, --glm, or --chmp is used, all of them must be used!" << endl;
	cout << endl;

	cout << "--get-testing-data <directory>" << endl;
	cout << "   The directory where the training and testing chromosomes will be stored." << endl;
	cout << "   These chromosomes are the ones used to auto-callibrate Look4TRs." << endl;
	cout << "   This will also store the locations of the inserted tandem repeats." << endl;
	cout << endl;

	cout << "--hmm-states <integer>" << endl;
	cout << "   The numeber of HMM states (Must be positive and even)" << endl;
	cout << "   This should be used when not wanting to see the upper and lower limits on HMM states for auto-callibration." << endl << endl;

	cout << "--thr <integer>                                                                 "
	<< endl;
	cout << "   Number of threads." << endl;

	cout 
	<< "||==================================LOOK4TRS-HELP==================================||"
	<< endl
	<< endl;
}

int main(int argc, char *argv[])
{
	string bottom =
	"====================================================================================";


	if (argc == 1)
	{

		printHelp();
		return 0;
	}

	if (!(argc == argCount || argc == 6 || argc == 2 || argc != argCount - 1))
	{
		cout << "Error: Incorrect number of arguments!: " << argc - 1 << endl;
		cout
		<< "Use \"--help\" or \"-h\" flag for assistance in using this tool."
		<< endl;
		exit(0);
	}

	double base = 2.0;

	string addr;
	string outfile;
	int minK;
	string trainDir;
	string trainFile;
	int maxK;
	int init_reg;
	int idn;
	int smt = 0;
	int mtf;
	int thr;
	int prn;
	int minm;
	int maxm;
	int order;
	int win;
	int will_merge = 0;
	int lng = 0;
	string hmm_file;
	string glm_file;
	string chmp_file;
	string save_file = "./";
	std::string chromTR_dir;
	int seg_size;
	int hmm_states;
	int hmm_state_lower;
	int hmm_state_upper;
	std::string bedMotifFile = "";
	std::string faMotifFile = "";

	string addrFlag = "--adr";
	string outfileFlag = "--out";
	string minKFlag = "--min";
	string maxKFlag = "--max";
	string trainingFlag = "--trn";
	string microFlag = "--default";
	string initFlag = "--int";
	string idnFlag = "--idn";
	string mtfFlag = "--mtf";
	string thrFlag = "--thr";
	string prnFlag = "--prn";
	string minmFlag = "--minm";
	string maxmFlag = "--maxm";
	string ordFlag = "--ord";
	string winFlag = "--win";
	string hmmFlag = "--hmm";
	string glmFlag = "--glm";
	string chmpFlag = "--chmp";
	string saveDataFlag = "--save-data";
	string segSizeFlag = "--seg";
	string chromTRFlag = "--get-testing-data";
	string hmmStateFlag = "--hmm-states";
	string hmmStateLFlag = "--hmm-states-lower";
	string hmmStateUFlag = "--hmm-states-upper";
	string motifFileFlag   = "--motif-location-file";
	string faMotifFileFlag = "--motif-fa-file";

	pair<string, bool> uninitVal("0", false);

	map<string, pair<string, bool>> optTable;

	optTable[addrFlag]        = uninitVal;
	optTable[outfileFlag]     = uninitVal;
	optTable[minKFlag]        = uninitVal;
	optTable[maxKFlag]        = uninitVal;
	optTable[trainingFlag]    = uninitVal;
	optTable[initFlag]        = uninitVal;
	optTable[idnFlag]         = uninitVal;
	optTable[mtfFlag]         = uninitVal;
	optTable[thrFlag]         = uninitVal;
	optTable[prnFlag]         = uninitVal;
	optTable[minmFlag]        = uninitVal;
	optTable[maxmFlag]        = uninitVal;
	optTable[ordFlag]         = uninitVal;
	optTable[winFlag]         = uninitVal;
	optTable[hmmFlag]         = uninitVal;
	optTable[glmFlag]         = uninitVal;
	optTable[chmpFlag]        = uninitVal;
	optTable[saveDataFlag]    = uninitVal;
	optTable[segSizeFlag]     = uninitVal;
	optTable[chromTRFlag]     = uninitVal;
	optTable[hmmStateFlag]    = uninitVal;
	optTable[hmmStateLFlag]   = uninitVal;
	optTable[hmmStateUFlag]   = uninitVal;
	optTable[motifFileFlag]   = uninitVal;
	optTable[faMotifFileFlag] = uninitVal;

	if (argc == 2)
	{
		string arg = argv[1];
		if (arg == string("--help") || arg == string("-h"))
		{
			printHelp();

			exit(0);
		}
		else if (arg == "--micro" || arg == "--full" || optTable.count(arg) == 1)
		{
			cerr << "Error: Incorrect number of arguments (" << argc - 1 << ") for the \"" << arg << "\" flag!" << endl;
			exit(1);
		}
		else
		{
			cerr << "Error: Unrecognized option!: " << arg << endl;
			exit(1);
		}
	}
	else if (argc == 6 || argc == 8)
	{
		for (int i = 1; i < argc; i++)
		{
			string arg = argv[i];
			if (arg == microFlag)
			{
				optTable[initFlag] = pair<string, bool>("20", true);
				optTable[minKFlag] = pair<string, bool>("4", true);
				optTable[maxKFlag] = pair<string, bool>("6", true);
				optTable[idnFlag] = pair<string, bool>("50", true);
				optTable[prnFlag] = pair<string, bool>("2000", true);
				optTable[mtfFlag] = pair<string, bool>("1", true);
				optTable[thrFlag] = pair<string, bool>(
					std::to_string(Util::CORE_NUM), true);
				optTable[minmFlag] = pair<string, bool>("1", true);
				optTable[maxmFlag] = pair<string, bool>("10", true);
				optTable[ordFlag] = pair<string, bool>("0", true);
				optTable[winFlag] = pair<string, bool>("4", true);
				optTable[segSizeFlag] = pair<string, bool>("3000000", true);
				optTable[hmmStateLFlag] = pair<string, bool>("6", true);
				optTable[hmmStateUFlag] = pair<string, bool>("12", true);
			}
		
			else if (arg == outfileFlag)
			{
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == trainingFlag)
			{
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == addrFlag)
			{
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == hmmFlag)
			{
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == glmFlag)
			{
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == chmpFlag)
			{
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == saveDataFlag)
			{
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == chromTRFlag)
			{
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == hmmStateFlag)
			{
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == hmmStateUFlag)
			{
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == hmmStateLFlag)
			{
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == motifFileFlag){
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == faMotifFileFlag){
				optTable[arg] == pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (optTable.count(arg) == 1 || arg == "--help")
			{
				cerr << "Error: Incorrect number of arguments (" << argc - 1 << ") for the \"" << arg << "\" flag!" << endl;
				exit(1);
			}
			else
			{
				cerr << "Error: Unrecognized option!: " << arg << endl;
				exit(1);
			}
		}
	}
	else
	{
		for (int i = 1; i < argc; i++)
		{
			string arg = argv[i];
			if (optTable.count(arg) == 1)
			{
				optTable.at(arg) = pair<string, bool>(argv[i + 1], true);
				i++;
			}
			else if (arg == microFlag || arg == "--full" || arg == "--help")
			{
				cerr << "Error: Incorrect number of arguments (" << argc - 1 << ") for the \"" << arg << "\" flag!" << endl;
				exit(1);
			}
			else
			{
				cerr << "Error: Unrecognized option!: " << arg << endl;
				exit(1);
			}
		}
	}

	for (auto it = optTable.begin(); it != optTable.end(); ++it)
	{
		if (it->first == outfileFlag)
		{
			outfile = it->second.first;
		}
		else if (it->first == minKFlag)
		{
			varSetter(it->second, minK, minKFlag);
		}
		else if (it->first == maxKFlag)
		{
			varSetter(it->second, maxK, maxKFlag);
		}
		else if (it->first == trainingFlag)
		{
			trainDir = it->second.first;
		}
		else if (it->first == addrFlag)
		{
			addr = it->second.first;
		}
		else if (it->first == initFlag)
		{
			varSetter(it->second, init_reg, initFlag);
		}
		else if (it->first == idnFlag)
		{
			varSetter(it->second, idn, idnFlag);
		}
		else if (it->first == mtfFlag)
		{
			varSetter(it->second, mtf, mtfFlag);
		}
		else if (it->first == thrFlag)
		{
			varSetter(it->second, thr, thrFlag);
		}
		else if (it->first == prnFlag)
		{
			varSetter(it->second, prn, prnFlag);
		}
		else if (it->first == minmFlag)
		{
			varSetter(it->second, minm, minmFlag);
		}
		else if (it->first == maxmFlag)
		{
			varSetter(it->second, maxm, maxmFlag);
		}
		else if (it->first == ordFlag)
		{
			varSetter(it->second, order, ordFlag);
		}
		else if (it->first == winFlag)
		{
			varSetter(it->second, win, winFlag);
		}
		else if (it->first == hmmFlag)
		{
			hmm_file = it->second.first;
		}
		else if (it->first == glmFlag)
		{
			glm_file = it->second.first;
		}
		else if (it->first == chmpFlag)
		{
			chmp_file = it->second.first;
		}
		else if (it->first == saveDataFlag)
		{
			save_file = it->second.first;
		}
		else if (it->first == segSizeFlag)
		{
			varSetter(it->second, seg_size, segSizeFlag);
		}
		else if (it->first == chromTRFlag)
		{
			chromTR_dir = it->second.first;
		}
		else if (it->first == hmmStateFlag)
		{
			varSetter(it->second, hmm_states, hmmStateFlag);
		}
		else if (it->first == hmmStateLFlag)
		{
			varSetter(it->second, hmm_state_lower, hmmStateLFlag);
		}
		else if (it->first == hmmStateUFlag)
		{
			varSetter(it->second, hmm_state_upper, hmmStateUFlag);
		}
		else if (it->first == motifFileFlag)
		{
			bedMotifFile = it->second.first;
		}
		else if (it->first == faMotifFileFlag){
			faMotifFile = it->second.first;
		}
		else
		{
			cerr << "OptTable Error: unrecognized OptTable Element!" << endl;
			exit(1);
		}
	}

	//Criteria for minK
	if (!optTable[minKFlag].second)
	{
		cerr << "Error: The minimum k-mer length (" << minKFlag << ") is uninitialized!" << endl;
		exit(1);
	}
	else if (minK <= 0)
	{
		cerr << "Error: The minimum k-mer length (" << minKFlag << ") must be greater than 0!" << endl;
		exit(1);
	}
	else if (minK > 15)
	{
		cerr << "Error: Minimum k-mer length (" << minK << ") must be less than 16!" << endl;
		exit(1);
	}

	//Criteria for maxK
	if (!optTable[maxKFlag].second)
	{
		cerr << "Error: The maximum k-mer length (" << maxKFlag << ") is uninitialized!" << endl;
		exit(1);
	}
	else if (maxK <= 0)
	{
		cerr << "Error: The maximum k-mer length (" << maxKFlag << ") must be greater than 0!" << endl;
		exit(1);
	}
	else if (maxK < minK)
	{
		cerr << "Error: Maximum k-mer length (" << maxK << ") is shorter than the minimum "
		<< "k-mer length (" << minK << ")!" << endl;
		exit(1);
	}
	else if (maxK > 15)
	{
		cerr << "Error: Maximum k-mer length (" << maxK << ") must be less than 16!" << endl;
		exit(1);
	}

	// Hani Z. Girgis added the following
	// Criteria for the identity score
	if (!optTable[idnFlag].second)
	{
		cerr << "Error: The identity score (" << idnFlag << ") is uninitialized!" << endl;
		exit(1);
	}
	else if (idn < 0 || idn > 100)
	{
		cerr << "Error: The identity score (" << idnFlag << ") must be between 0 and 100!" << endl;
		exit(1);
	}

	// End of Hain's modifications

	//Criteria for addr
	if (!optTable[addrFlag].second)
	{
		cerr << "Error: Address of files to be processed (" << addrFlag << ") is "
		<< "uninitialized!" << endl;
		exit(1);
	}

	//Criteria for outfile
	if (!optTable[outfileFlag].second)
	{
		cerr << "Error: Address of output file (" << outfileFlag << ") is uninitialized!" << endl;
		exit(1);
	}
	else
	{
		Util::checkDir(optTable[outfileFlag].first);
	}

	//criteria for trainDir
	if (!optTable[trainingFlag].second)
	{
		trainFile = Util::getLargestFile(addr);
	}
	else
	{
		vector<string> fileList;
		Util::readChromList(trainDir, &fileList, string("fa"));

		if (fileList.size() == 0)
		{
			cerr << "The training folder is empty!" << endl;
			throw std::exception();
		}
		else if (fileList.size() != 1)
		{
			cerr << "The training folder must include ";
			cerr << "exactly one long-enough chromosome";
			cerr << endl;
			throw std::exception();
		}
		trainFile = fileList.at(0);
	}

	//Criteria for init_reg
	if (!optTable[initFlag].second)
	{
		cerr << "Error: The initial region size (" << initFlag << ") is uninitialized!" << endl;
		exit(1);
	}
	else if (init_reg <= 0)
	{
		cerr << "Error: The initial region size (" << initFlag << ") must be greater than 0!" << endl;
		exit(1);
	}

	if (!optTable[mtfFlag].second)
	{
		cerr << "Error: Want motif finding " << mtfFlag << " is uninitialized!"
		<< endl;
		exit(1);
	}
	else if (mtf != 0 && mtf != 1)
	{
		cerr << "Want motif finding should be 0 or 1: " << mtf << endl;
		exit(1);
	}

	if (!optTable[thrFlag].second)
	{
		cerr << "Error: The number of threads " << thrFlag
		<< " is uninitialized!" << endl;
		exit(1);
	}
	else if (thr < 1)
	{
		cerr << "Error: The number of threads " << thrFlag
		<< " must be 1 or more!" << endl;
		exit(1);
	}
	else
	{
		Util::CORE_NUM = min((int)Util::CORE_NUM + 1, thr);
	}

	if (!optTable[prnFlag].second)
	{
		cerr << "Error: The number of times to train motif analyzer " << prnFlag
		<< " is uninitialized!" << endl;
		exit(1);
	}
	else if (prn < 1)
	{
		cerr << "Error: The number of times to train motif analyzer " << prnFlag
		<< " must be 1 or more!" << endl;
		exit(1);
	}

	//Criteria for order
	if (!optTable[ordFlag].second)
	{
		cerr << "Error: the order " << ordFlag
		<< " is uninitialized!" << endl;
		exit(1);
	}
	else if (order < 0)
	{
		cerr << "Error: the order " << order << " must be non-negative!" << endl;
		exit(1);
	}

	if(optTable[minmFlag].second && optTable[motifFileFlag].second){
		cerr << "-minm and --motif-location-file have both been initialized! " <<
		"Only one needs to be initialized!"<<
		endl;
		exit(1);
	}

	if(optTable[maxmFlag].second && optTable[motifFileFlag].second){
		cerr << "-maxm and --motif-location-file have both been initialized! " <<
		"Only one needs to be initialized!"<<
		endl;
		exit(1);
	}


	//Criteria for minm
	if (!optTable[minmFlag].second && !optTable[motifFileFlag].second)
	{
		if( !optTable[hmmFlag].second ||!optTable[glmFlag].second || !optTable[chmpFlag].second )
		{
			cerr << "Error: the minimum motif size " << minmFlag
			<< " is uninitialized!" << endl;
			cerr << "Either initialize minm or initialize --motif-location-file" << endl;
			exit(1);
		}
	}
	else if (minm <= 0 && !optTable[motifFileFlag].second)
	{
		if( !optTable[hmmFlag].second ||!optTable[glmFlag].second || !optTable[chmpFlag].second )
		{
			cerr << "Error: The minimum " << minmFlag << " " << minm << " must be 1 or more!"
			<< endl;
			exit(1);
		}
	}

	//Criteria for maxm
	if (!optTable[maxmFlag].second  && !optTable[motifFileFlag].second)
	{
		if( !optTable[hmmFlag].second ||!optTable[glmFlag].second || !optTable[chmpFlag].second )
		{
			cerr << "Error: the maximum motif size " << maxmFlag
			<< " is uninitialized!" << endl;
			cerr << "Either initialize maxm or initialize --motif-location-file" << endl;
			exit(1);
		}
	}
	else if (maxm <= 0 && !optTable[motifFileFlag].second)
	{
		if( !optTable[hmmFlag].second ||!optTable[glmFlag].second || !optTable[chmpFlag].second )
		{
			cerr << "Error: The maximum motif" << maxmFlag << " " << maxm << " must be 1 or more!"
			<< endl;
			exit(1);
		}
	}
	else if (maxm < minm  && !optTable[motifFileFlag].second)
	{
		cerr << "Error: The maximum motif " << maxm << " must be larger than the minimum motif " << minm
		<< endl;
		exit(1);
	}

	// Criteria for window loops
	if (!optTable[winFlag].second)
	{
		cerr << "Error: the window loops " << winFlag
		<< " is uninitialized!" << endl;
		exit(1);
	}
	else if (win <= 0)
	{
		cerr << "Error: The window loops " << winFlag << " must be 1 or more!" << endl;
		exit(1);
	}

	if (optTable[hmmFlag].second)
	{
		Util::checkFile(hmm_file);
	}
	else
	{
		hmm_file = "";
	}

	if (optTable[glmFlag].second)
	{
		Util::checkFile(glm_file);
	}
	else
	{
		glm_file = "";
	}

	if (optTable[chmpFlag].second)
	{
		Util::checkFile(chmp_file);
	}
	else
	{
		chmp_file = "";
	}

	if (optTable[saveDataFlag].second)
	{
		Util::mkDir(save_file);
	}
	else
	{
		save_file = "./";
	}

	if (!optTable[segSizeFlag].second)
	{
		cerr << "Error: the segment size " << segSizeFlag << " is uninitialized!" << endl;
		exit(1);
	}
	else if (seg_size <= 0)
	{
		cerr << "Error: the segment size " << seg_size << " must be positive!" << std::endl;
	}

	if (optTable[chromTRFlag].second)
	{
		Util::mkDir(chromTR_dir);
	}
	else
	{
		chromTR_dir = "";
	}

	if (!optTable[hmmStateFlag].second && !optTable[hmmStateLFlag].second && !optTable[hmmStateUFlag].second)
	{
		hmm_state_lower = 10;
		hmm_state_upper = 10;
	}
	else if (optTable[hmmStateFlag].second && !optTable[hmmStateLFlag].second && !optTable[hmmStateUFlag].second)
	{
		if (hmm_states <= 0)
		{
			std::cerr << "The number of hmm states cannot be negative(" << hmm_states << ")!" << std::endl;
			exit(1);
		}
		else if (hmm_states % 2 == 1)
		{
			std::cerr << "The number of HMM states must be even(" << hmm_states << ")!" << std::endl;
			exit(1);
		}
		else
		{
			hmm_state_lower = hmm_states;
			hmm_state_upper = hmm_states;
		}
	}
	else if (!optTable[hmmStateFlag].second && optTable[hmmStateLFlag].second && optTable[hmmStateUFlag].second)
	{
		if (hmm_state_lower <= 0)
		{
			std::cerr << "The lower number of hmm states cannot be negative(" << hmm_states << ")!" << std::endl;
			exit(1);
		}
		else if (hmm_state_lower % 2 == 1)
		{
			std::cerr << "The lower number of HMM states must be even(" << hmm_states << ")!" << std::endl;
			exit(1);
		}

		if (hmm_state_upper <= 0)
		{
			std::cerr << "The lower number of hmm states cannot be negative(" << hmm_states << ")!" << std::endl;
			exit(1);
		}
		else if (hmm_state_upper % 2 == 1)
		{
			std::cerr << "The lower number of HMM states must be even(" << hmm_states << ")!" << std::endl;
			exit(1);
		}
		if (hmm_state_lower > hmm_state_upper)
		{
			std::cerr << "The lower number of HMM states(" << hmm_state_lower << ") cannot be greater than the upper number of HMM states(" << hmm_state_upper << ")!" << std::endl;
			exit(1);
		}
	}
	else if (optTable[hmmStateLFlag].second && !optTable[hmmStateUFlag].second)
	{
		std::cerr << "The lower number of HMM states was initialized, but the upper one wasn't!" << std::endl;
		exit(1);
	}
	else if (!optTable[hmmStateLFlag].second && optTable[hmmStateUFlag].second)
	{
		std::cerr << "The upper number of HMM states was initialized, but the lower one wasn't!" << std::endl;
		exit(1);
	}
	else if (optTable[hmmStateFlag].second && optTable[hmmStateLFlag].second && optTable[hmmStateUFlag].second)
	{
		std::cerr << "Cannot use the " << hmmStateFlag << " and the " << hmmStateLFlag << " and " << hmmStateUFlag << " parameters at the same time!" << std::endl;
		exit(1);
	}
	else
	{
		std::cerr << "This HMM State option was not explored" << std::endl;
		exit(1);
	}

	if(!optTable[motifFileFlag].second){
		bedMotifFile = "";
	}
	else if(optTable[motifFileFlag].second && optTable[motifFileFlag].second && optTable[motifFileFlag].second){
		Util::checkFile(bedMotifFile);
	}

	if(!optTable[faMotifFileFlag].second){
		faMotifFile = "";
	} else if (optTable[faMotifFileFlag].second && !optTable[motifFileFlag].second){
		std::cout << "The motif.fa file has been instantiated with --motif-fa-file, but the corresponding bed file to find its motifs has not been!" << std::endl;
		std::cout << "Please use the --motif-location-file flag to declare the bed file to use to obtain the repeats for " << faMotifFile << std::endl;
		exit(1);
	} else if(optTable[faMotifFile].second){

	}

	Executor e(addr, base, init_reg, order,
		outfile, minK, win, trainFile, maxK, (idn / 100.0), smt, mtf,
		will_merge = 0, prn, lng, minm, maxm, seg_size,
		hmm_state_lower, hmm_state_upper,
		chromTR_dir, hmm_file, glm_file, chmp_file, save_file, bedMotifFile, faMotifFile);
}
