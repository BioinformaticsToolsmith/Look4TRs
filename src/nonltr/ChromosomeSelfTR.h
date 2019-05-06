#ifndef _CHROMSOMESELFTR_H_
#define _CHROMSOMESELFTR_H_

#include <string>
#include <random>
#include "../utility/Location.h"
#include "ChromosomeTR.h"

namespace nonltr {

	class ChromosomeSelfTR : public ChromosomeTR {

	public:
		ChromosomeSelfTR(int, ChromosomeOneDigit*, char, int, int, int, int, int seed_in = 0);
		~ChromosomeSelfTR();
		void shuffle();
		std::string getRandTR();
		void printBedData(std::string, std::string = "");
	private:
		std::vector<std::string> * motifList;
		int init_reg;
		int minMotif, maxMotif;
		int mutationRate = 0;

	};


}
#endif
