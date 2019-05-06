#ifndef _CHROMSOMEREADTR_H_
#define _CHROMSOMEREADTR_H_

#include <string>
#include <random>
#include "../utility/Location.h"
#include "ChromosomeTR.h"

namespace nonltr {

	class ChromosomeReadTR : public ChromosomeTR {

	public:
		ChromosomeReadTR(int, ChromosomeOneDigit*, char, int, std::string, std::string, int seed_in = 0);
		ChromosomeReadTR(int, ChromosomeOneDigit*, char, int, std::string, std::string, ChromosomeOneDigit *, int seed_in = 0);
		~ChromosomeReadTR();
		void shuffle();
		std::string getRandTR();
		void printBedData(std::string, std::string = "");
	private:
		std::vector<std::string> * readTRList;
		int seed = 0;
		int TR_file_index = 0;

	};


}
#endif
