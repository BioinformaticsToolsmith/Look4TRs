using namespace std;

#include "FindMotif.h"


namespace motif {

	FindMotif::FindMotif(string sequenceIn, double thresholdIn,
		Predictor<int> * predIn) {
		sequence = sequenceIn;
		threshold = thresholdIn;
		smoothingWindow = -1;
		h1 = ">Sequence";

		isFound = false;
		foundMotif = string("");
		identityScore = -1;
		pred = predIn;
		ChromosomeOneDigit seqCode(sequence, h1);
		seqPoint = pred->get_point(&seqCode);

		searchMicro();
	}

	

	FindMotif::~FindMotif() {
		delete seqPoint;

	}

/*
 *   Turns a deque of characters into a string
 */
 string FindMotif::dequeToString(deque<char>& l) {
 	string result(l.size(), '\0');
 	for (int i = 0; i < l.size(); i++) {
 		result.at(i) = l.at(i);
 	}
 	return result;
 }

/*
 * This makes an exact repeat from the word. The repeat will be of size len
 */
 string FindMotif::makeExact(string word, int len) {
 	string result = word;

 	if (word.size() == 0) {
 		std::cerr << "You cannot make an exact repeat of an empty word!"
 		<< std::endl;
 		throw std::exception();
 	}

 	if (len < word.size()) {
 		cerr << "The len " << len << " requested is smaller than the word "
 		<< word << endl;
 		throw std::exception();
 	}
 	for (int i = 0; result.size() != len; i++) {
 		result += result.at(i);
 	}
 	if (result.size() != len) {
 		std::cerr << "makeExact is not the same length as the input len!"
 		<< std::endl;
 		throw std::exception();
 	}
 	return result;
 }

/**
 * This method takes in a non-empty vector of strings that represent words.
 * These words represent candidate motifs.
 * This method will see which candidate makes the best motif and sets it if one is found.
 */

 void FindMotif::greedyConfirmation(vector<string> * copyList) {
 	static int callCount = 0;
 	if (copyList->empty()) {
 		std::cerr << "The given copyList in the greedyConfirmation is empty!"
 		<< std::endl;
 		throw std::exception();
 	}

	// Find the best exact repeat, i.e. (motif)
 	pair<string, double> result;
 	result.first = string("");
 	result.second = -1;

	// Collects the distances into a list
 	for (int i = 0; i < copyList->size(); i++, callCount++) {

 		string copy = copyList->at(i);


		// Make exact repeat, i.e. (motif)n
 		string w = makeExact(copy, sequence.size());

 		ChromosomeOneDigit wCode(w, copy);

 		Point<int> * wPoint = pred->get_point(&wCode);

 		double similarity = pred->similarity(seqPoint, wPoint);
 		delete wPoint;

 		if (fabs(result.second - similarity) < std::numeric_limits<double>::epsilon() && copy.size() < result.first.size()) {
 			result.first = copy;
 			result.second = similarity;
 		} 
 		else if (similarity > result.second) {
 			result.first = copy;
 			result.second = similarity;
 		}
 	}

 	if (result.second < 0 || result.first.size() == 0) {
 		if(result.second == -1){
 			std::cerr << "The score is -1" << std::endl;
 		}
 		if(result.first == ""){
 			std::cerr << "The result is an empty string" << std::endl;
 		}
 		for(int i = 0; i < copyList->size(); i++){
 			std::cerr << copyList->at(i) << ", ";
 		}
 		std::cerr << std::endl;
 		std::cerr << "The result from greedyConfirmation is an empty string!"
 		<< std::endl;
 		throw std::exception();
 	}

 	foundMotif = result.first;
 	identityScore = result.second;
 	isFound = result.second > threshold;

 }



/*
 This method will create a vector of words that were found in the sequence.
 These words will be compared to the sequence.
 This will call greedyConfirmation.
 Author: Alfredo Velasco
 */
 void FindMotif::searchMicro() {
	map<string, int> wordSet; // Set of all words seen in the sequence
	vector<string> * wordList = new vector<string>();
	for (int size = 1; size <= MICRO_MAX_SIZE && size <= sequence.size();
		size++) {
		deque<char> word;
	for (int i = 0; i < sequence.size(); i++) {
		word.push_back(sequence.at(i));
		if (word.size() > size) {
			word.pop_front();
		} else if (word.size() < size) {
			continue;
		}

		string w = dequeToString(word);

		if (wordSet.count(w) == 0) {
			wordSet.emplace(w, 1);
		} else {
			wordSet.at(w)++;}
		}
	}

	for (auto it = wordSet.begin(); it != wordSet.end(); it++) {
		if (it->second > 1) {
			wordList->push_back(it->first);
		}
	}

	if (!wordList->empty()) {
		greedyConfirmation(wordList);
		wordList->clear();
	} else {
		isFound = false;
		foundMotif = "";
		identityScore = 0;
	}
	delete wordList;
}

bool FindMotif::getIsFound() {
	return isFound;
}

string FindMotif::getFoundMotif() {
	return foundMotif;
}

double FindMotif::getIdentityScore() {
	return identityScore;
}

}
