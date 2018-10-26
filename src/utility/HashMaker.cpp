/*
 * HashMaker.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: Hani Zakaria Girgis, PhD
 */
namespace utility {

template<class I>
HashMaker<I>::HashMaker(int k) {
	this->k = k;

	// Initialize bases
	bases = new I[k];
	for (int i = k - 1; i >= 0; i--) {
		bases[k - 1 - i] = (I) pow(4.0, i);
	}

	// Initialize mMinusOne
	mMinusOne = new I[4];
	for (int i = 0; i < 4; i++) {
		mMinusOne[i] = i * bases[0];
	}

}

template<class I>
HashMaker<I>::~HashMaker() {
	delete[] bases;
	delete[] mMinusOne;
}

template<class I>
I HashMaker<I>::hash(const char * key) {
	return hash(key, 0);
}

template<class I>
I HashMaker<I>::hash(const char * sequence, int keyStart) {
	I index = 0;
	for (int i = 0; i < k; i++) {
		char nucleotide = sequence[keyStart + i];
		if (nucleotide >= 0 && nucleotide <= 3) {
			index += bases[i] * sequence[keyStart + i];
		} else {
			string msg("The value of the char representing the nucleotide ");
			msg.append("must be between 0 and 3.");
			msg.append("The int value is ");
			msg.append(Util::int2string((int) nucleotide));
			msg.append(" of nucleotide at index ");
			msg.append(Util::int2string(keyStart + i));

			for (int h = 0 + keyStart; h < k + keyStart; h++) {
				cerr << (int) sequence[h];
			}
			cerr << endl;

			throw InvalidInputException(msg);
		}
	}
	return index;
}

template<class I>
void HashMaker<I>::hash(const char * sequence, int start, int end,
		vector<I> * hashList) {

	for (int i = start; i <= end; i++) {
		char nucleotide = sequence[i];

		if (!(nucleotide >= 0 && nucleotide <= 3)) {
			string msg("The value of the char representing the nucleotide ");
			msg.append("must be between 0 and 3.");
			msg.append("The int value is ");
			msg.append(Util::int2string((int) nucleotide));
			msg.append(" of nucleotide at index ");
			msg.append(Util::int2string(i));

			throw InvalidInputException(msg);
		}
	}

	I lastHash = hash(sequence, start);
	hashList->push_back(lastHash);

	for (int i = start + 1; i <= end; i++) {
		I s1 = 4 * (lastHash - mMinusOne[(int) sequence[i - 1]])
				+ (int) sequence[i + k - 1];
		hashList->push_back(s1);
		lastHash = s1;
	}
}

} /* namespace utility */
