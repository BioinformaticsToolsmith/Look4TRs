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

 bool has_only_digits(const string s) {
 	return s.find_first_not_of("-.0123456789") == string::npos;
 }

 void varSetter(pair<string, bool> &val, int &var, string varFlag) {

 	if (has_only_digits(val.first.c_str())) {
 		var = atoi(val.first.c_str());
 	} else {
 		cerr << "Error: \"" << varFlag << "\" must be an integer! You passed: "
 		<< val.first.c_str() << endl;
 		exit(1);
 	}
 }

 void printHelp() {
 	string top =
 	"[#]------------------------------------------------------------------------------[#]";
 	cout
 	<< "||================================LOOK4TRS-HELP================================||"
 	<< endl << endl;
 	cout
 	<< "- - - - - - - - - - - - - -ACCEPTABLE COMMAND FORMATS- - - - - - - - - - - - - -"
 	<< endl;
 	cout << "General Command Format" << endl;
 	cout << "   ./Look4TRs --help" << endl;
 	cout << "Command Format Using Default Option Flags" << endl;
 	cout << "   ./Look4TRs [required options]" << endl;
 	cout << "Command Format Using the Parameter Option Flags" << endl;
 	cout << "   ./Look4TRs [required options][parameter options]"
 	<< endl << endl;
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
 	<< "    predicted.      " << endl;

 	cout << "--out <string>" << endl;
 	cout << "   Takes the address of the output prediction file." << endl
 	<< endl;
 	cout << top << endl;
 	cout
 	<< " | NOTE: After supplying the required options, either provide a single default    | "
 	<< endl
 	<< " | option, or provide all of the parameter options!                               | "
 	<< endl;
 	cout << top << endl << endl;
 	cout
 	<< "- - - - - - - - - - - - - - - - -DEFAULT OPTION - - - - - - - - - - - - - - - -"
 	<< endl;
 	cout
 	<< "--micro                                                                         "
 	<< endl;
 	cout
 	<< "    Runs the tool with the default parameters for micro satellites."
 	<< endl << endl;
 	cout
 	<< "- - - - - - - - - - - - - - - - PARAMETER OPTIONS- - - - - - - - - - - - - - - -"
 	<< endl;
 	cout << top << endl;
 	cout
 	<< " | NOTE: All parameter options must be greater than 0!                            | "
 	<< endl;
 	cout << top << endl;
 	cout << "--trn <string>" << endl;
 	cout << "    Takes the address of the directory of the training files."
 	<< endl;
 	cout
 	<< "--min <integer>                                                                 "
 	<< endl;
	cout << "    The minimum k-mer length analyzed (Must be less than or equal" //
		<< endl;
		cout << "   to the maximum k-mer length)." << endl;
cout
<< "--max <integer>                                                                 "
<< endl;
cout
			<< "    The maximum k-mer length analyzed (Must be greater than or equal" //)
			<< endl;
			cout << "   to the minimum k-mer length)." << endl;
cout
<< "--int <integer>                                                                 "
<< endl;
cout << "   The minimum allowed region size for predictions." << endl;

cout
<< "--idn <integer>                                                                 "
<< endl;
cout
<< "    The minimum allowed identity score between the exact repeat and a candidate"
<< endl;

cout
<< "--mtf <integer>                                                                 "
<< endl;
cout << "    Enable ('1') or disable ('0') the motif discovery feature."
<< endl;

cout
<< "--thr <integer>                                                                 "
<< endl;
cout << "   Number of threads." << endl;
cout << "--prn <integer>" << endl;
cout
<< "     The number of times that the motif analyzer will train on motifs. Must be positive."
<< endl;
cout << "--kmr <integer>" << endl;
cout << "     The k-mer size used for motif analysis. Must be positive." << endl;

cout

<< "||================================LOOK4TRS-HELP=================================||"
<< endl << endl;

}

int main(int argc, char *argv[]) {
	string top =
	"[#]------------------------------------------------------------------------------[#]";
	string bottom =
	"====================================================================================";
	//cout << top << endl;
	//cout << bottom << endl;
	cout << top << endl << endl;

	string author1 =
	" |                                     Authors                                    | ";
	string author2 =
	" |  Alfredo Velasco II*, Benjamin T. James, Vincent D. Wells*, and Hani Z. Girgis | ";
	string author3 =
	" |                 *These authors contributed equally to this work                | ";

	string author4 =
	" |                                    Thanks to                                   | ";
	string author5 =
	" |      Joseph Valencia for coding the global alignment algorithm                 | ";

	cout << author1 << endl;
	cout << author2 << endl;
	cout << author3 << endl;
	cout << author4 << endl;
	cout << author5 << endl;

	string lab =
	" |         The Bioinformatics Toolsmith Laboratory, The University of Tulsa       | ";
	string version =
	" |                                    10/19/2018                                  | ";
	cout << top << endl;
	// cout << author << endl;
	cout << lab << endl;
	cout << version << endl;
	cout << top << endl;
	cout << endl;

	if (argc == 1) {

		printHelp();
		return 0;
	}

	if (!(argc == argCount || argc == 6 || argc == 2 || argc != argCount - 1)) {
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
	int mtf;
	int thr;
	int prn;
	int kmr;

	string addrFlag = "--adr";
	string outfileFlag = "--out";
	string minKFlag = "--min";
	string maxKFlag = "--max";
	string trainingFlag = "--trn";
	string microFlag = "--micro";
	string initFlag = "--int";
	string idnFlag = "--idn";
	string mtfFlag = "--mtf";
	string thrFlag = "--thr";
	string prnFlag = "--prn";
	string kmrFlag = "--kmr";

	pair<string, bool> uninitVal("0", false);

	map<string, pair<string, bool> > optTable;

	optTable[addrFlag] = uninitVal;
	optTable[outfileFlag] = uninitVal;
	optTable[trainingFlag] = uninitVal;

	optTable[initFlag] = pair<string, bool>("20", true);
	optTable[minKFlag] = pair<string, bool>("4", true);
	optTable[maxKFlag] = pair<string, bool>("6", true);
	optTable[idnFlag] = pair<string, bool>("50", true);
	optTable[kmrFlag] = pair<string, bool>("3", true);
	optTable[prnFlag] = pair<string, bool>("2000", true);
	optTable[mtfFlag] = pair<string, bool>("1", true);
	optTable[thrFlag] = pair<string, bool>(
		std::to_string(Util::CORE_NUM), true);

	if (argc == 2) {
		string arg = argv[1];
		if (arg == string("--help") || arg == string("-h")) {
			printHelp();

			exit(0);
		} else if ( optTable.count(arg) == 1) {
			cerr << "Error: Incorrect number of arguments (" << argc - 1 << ") for the \"" << arg << "\" flag!" << endl;
			exit(1);
		} else {
			cerr << "Error: Unrecognized option!: " << arg << endl;
			exit(1);
		}

	} else if (argc == 6) {
		for (int i = 1; i < argc; i++) {
			string arg = argv[i];
			if (arg == outfileFlag) {
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			} else if (arg == trainingFlag) {
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			} else if (arg == addrFlag) {
				optTable[arg] = pair<string, bool>(argv[i + 1], true);
				i++;
			} else if (optTable.count(arg) == 1 || arg == "--help") {
				cerr << "Error: Incorrect number of arguments (" << argc - 1 << ") for the \"" << arg << "\" flag!" << endl;
				exit(1);
			} else {
				cerr << "Error: Unrecognized option!: " << arg << endl;
				exit(1);
			}
		}
	} 
	else {
		for (int i = 1; i < argc; i++) {
			string arg = argv[i];
			if (optTable.count(arg) == 1) {
				optTable.at(arg) = pair<string, bool>(argv[i + 1], true);
				i++;
			} else if (arg == "--help") {
				cerr << "Error: Incorrect number of arguments (" << argc - 1 << ") for the \"" << arg << "\" flag!" << endl;
				exit(1);
			} else {
				cerr << "Error: Unrecognized option!: " << arg << endl;
				exit(1);
			}
		}
	}

	for (auto it = optTable.begin(); it != optTable.end(); ++it) {
		if (it->first == outfileFlag) {
			outfile = it->second.first;
		} else if (it->first == minKFlag) {
			varSetter(it->second, minK, minKFlag);
		} else if (it->first == maxKFlag) {
			varSetter(it->second, maxK, maxKFlag);
		} else if (it->first == trainingFlag) {
			trainDir = it->second.first;
		} else if (it->first == addrFlag) {
			addr = it->second.first;
		} else if (it->first == initFlag) {
			varSetter(it->second, init_reg, initFlag);
			//remove the next two else ifs after testing
		} else if (it->first == idnFlag) {
			varSetter(it->second, idn, idnFlag);
		} else if (it->first == mtfFlag) {
			varSetter(it->second, mtf, mtfFlag);
		} else if (it->first == thrFlag) {
			varSetter(it->second, thr, thrFlag);
		} else if (it->first == prnFlag) {
			varSetter(it->second, prn, prnFlag);
		} else if (it->first == kmrFlag) {
			varSetter(it->second, kmr, kmrFlag);
		} else {
			cerr << "OptTable Error: unrecognized OptTable Element!" << endl;
			exit(1);
		}
	}

	//Criteria for minK
	if (!optTable[minKFlag].second) {
		cerr << "Error: The minimum k-mer length (" << minKFlag << ") is uninitialized!" << endl;
		exit(1);
	} else if (minK <= 0) {
		cerr << "Error: The minimum k-mer length (" << minKFlag << ") must be greater than 0!" << endl;
		exit(1);
	} else if (minK > 15) {
		cerr << "Error: Minimum k-mer length (" << minK << ") must be less than 16!" << endl;
		exit(1);
	}

	//Criteria for maxK
	if (!optTable[maxKFlag].second) {
		cerr << "Error: The maximum k-mer length (" << maxKFlag << ") is uninitialized!" << endl;
		exit(1);
	} 
	else if (maxK <= 0) {
		cerr << "Error: The maximum k-mer length (" << maxKFlag << ") must be greater than 0!" << endl;
		exit(1);
	} else if (maxK < minK) {
		cerr << "Error: Maximum k-mer length (" << maxK << ") is shorter than the minimum " << "k-mer length (" << minK << ")!" << endl;
		exit(1);
	} else if (maxK > 15) {
		cerr << "Error: Maximum k-mer length (" << maxK << ") must be less than 16!" << endl;
		exit(1);
	}

	// Hani Z. Girgis added the following
	// Criteria for the identity score
	if (!optTable[idnFlag].second) {
		cerr << "Error: The identity score (" << idnFlag << ") is uninitialized!" << endl;
		exit(1);
	} else if (idn < 0 || idn > 100) {
		cerr << "Error: The identity score (" << idnFlag << ") must be between 0 and 100!" << endl;
		exit(1);
	}

	// End of Hain's modifications

	//Criteria for addr
	if (!optTable[addrFlag].second) {
		cerr << "Error: Address of files to be processed (" << addrFlag << ") is " << "uninitialized!" << endl;
		exit(1);
	}

	//Criteria for outfile
	if (!optTable[outfileFlag].second) {
		cerr << "Error: Address of output file (" << outfileFlag << ") is uninitialized!" << endl;
		exit(1);
	} else {
		Util::checkDir(optTable[outfileFlag].first);
	}

	//criteria for trainDir
	if (!optTable[trainingFlag].second) {
		trainFile = Util::getLargestFile(addr);
	} else {
		vector<string> fileList;
		Util::readChromList(trainDir, &fileList, string("fa"));

		if (fileList.size() == 0) {
			cerr << "The training folder is empty!" << endl;
			throw std::exception();
		} else if (fileList.size() != 1) {
			cerr << "The training folder must include ";
			cerr << "exactly one long-enough chromosome";
			cerr << endl;
			throw std::exception();
		}
		trainFile = fileList.at(0);
	}

	//Criteria for init_reg
	if (!optTable[initFlag].second) {
		cerr << "Error: The initial region size (" << initFlag << ") is uninitialized!" << endl;
		exit(1);
	} else if (init_reg <= 0) {
		cerr << "Error: The initial region size (" << initFlag << ") must be greater than 0!" << endl;
		exit(1);
	}

	if (!optTable[mtfFlag].second) {
		cerr << "Error: Want motif finding " << mtfFlag << " is uninitialized!"
		<< endl;
		exit(1);
	} else if (mtf != 0 && mtf != 1) {
		cerr << "Want motif finding should be 0 or 1: " << mtf << endl;
		exit(1);
	}

	if (!optTable[thrFlag].second) {
		cerr << "Error: The number of threads " << thrFlag
		<< " is uninitialized!" << endl;
		exit(1);
	} else if (thr < 1) {
		cerr << "Error: The number of threads " << thrFlag
		<< " must be 1 or more!" << endl;
		exit(1);
	} else {
		Util::CORE_NUM = min((int) Util::CORE_NUM + 1, thr);
	}

	if (!optTable[prnFlag].second) {
		cerr << "Error: The number of times to train motif analyzer " << prnFlag
		<< " is uninitialized!" << endl;
		exit(1);
	} else if (prn < 1) {
		cerr << "Error: The number of times to train motif analyzer " << prnFlag
		<< " must be 1 or more!" << endl;
		exit(1);
	}

	if (!optTable[kmrFlag].second && mtf) {
		cerr << "Error: The k-mer size " << kmrFlag << " is uninitialized!"
		<< endl;
		exit(1);
	} else if (kmr < 1  && mtf) {
		cerr << "Error: The k-mer size " << kmrFlag << " must be 1 or more!"
		<< endl;
		exit(1);
	}

	//remove filt1 and filt2 after testing
	Executor e(addr, base, init_reg,
		outfile, minK, trainFile, maxK,
		(idn / 100.0), mtf, kmr, prn);
}
