/*
* Author: Alfredo Velasco II
* Bioinformatics Toolsmith Laboratory, the University of Tulsa
* This class takes in a string and computes
* the score needed to find it's motif and
* align it against itself, but efficiently
* The alignRepeatWithBackTracking method uses 2 rows
* instead of n * (n + 1) / 2
* Resources: DIMACS Educational Module Series Module 09-2: Finding Repeats Within Strings
*     Algorithms in Bioinformatics I, WS’06, ZBIT, C. Dieterich, October 19, 2006
**/

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdlib.h>

#ifndef REPEATALIGNE_H_
#define REPEATALIGNE_H_

namespace motif{
	class RepeatAlignE {
	private:
		std::string seq;
		int* cost;
	public:
		RepeatAlignE(std::string);
		~RepeatAlignE();
		void alignRepeatWithBackTracking();
		std::vector<int> getLastRow();

	};
}

#endif
