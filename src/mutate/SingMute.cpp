#include "SingMute.h"
#include <set>
#include <random>

int intRandMod(int max) {
    static thread_local std::mt19937 generator;
    std::uniform_int_distribution<int> distribution(0, max-1);
    return distribution(generator);
}
long longRandMod(long max) {
    static thread_local std::mt19937 generator;
    std::uniform_int_distribution<long> distribution(0, max-1);
    return distribution(generator);
}


void generate_unique_set(size_t cmd_size, std::set<long>& ret, int num_elts, const std::set<long>& bad_set_1, const std::set<long>& bad_set_2, const std::vector<bool> &valid)
{
	while (ret.size() <= num_elts) {
		long idx = longRandMod(cmd_size);
		if (valid[idx] &&
		    ret.find(idx) == ret.end() &&
		    bad_set_1.find(idx) == bad_set_1.end() &&
		    bad_set_2.find(idx) == bad_set_2.end()) {

			ret.insert(idx);
		}
	}
}
char SingMute::randNucl()
{
	char character;
	int value = intRandMod(percAs + percCs + percGs + percTs);
	if (value < percAs) {
		character = (char) 0;
	} else if (value < percAs + percCs) {
		character = (char) 1;
	} else if (value < percAs + percCs + percGs) {
		character = (char) 2;
	} else {
		character = (char) 3;
	}
	return character;
}
void SingMute::init(const std::vector<bool> &valid)
{
	maxInsert = 0;
	maxDel = 0;
	maxSwitch = 0;
	if (num_mut == 1) {
		maxInsert = 1;
		maxDel = 0;
		maxSwitch = 0;
	} else {
		maxSwitch = intRandMod(num_mut);
		num_mut -= maxSwitch;

		if (maxSwitch % 2 == 1 && num_mut >= 1) {
			maxSwitch++;
			num_mut--;
		} else if (num_mut == 0) {
			maxSwitch--;
			num_mut++;
		}
		if (num_mut > 1) {
			maxInsert = intRandMod(num_mut);
			num_mut -= maxInsert;
		} else {
			maxInsert = num_mut;
			num_mut -= maxInsert;
		}
		maxDel = num_mut;
	}
	size_t seq_len = seq->length();

	maxDel *= seq_len / 100.0;
	maxInsert *= seq_len / 100.0;
	maxSwitch *= seq_len / 100.0;
	alignmentLength = maxInsert;
	IBP = maxDel + maxSwitch;



	std::vector<char> command_str(seq_len, 'S');

	std::set<long> s_ins, s_del, s_switch;
	generate_unique_set(command_str.size(), s_ins, maxInsert, s_del, s_switch, valid);
	generate_unique_set(command_str.size(), s_del, maxDel, s_ins, s_switch, valid);
	generate_unique_set(command_str.size(), s_switch, maxSwitch, s_ins, s_del, valid);
	for (auto idx : s_ins) {
		command_str[idx] = 'I';
	}
	for (auto idx : s_del) {
		command_str[idx] = 'D';
	}
	for (auto idx : s_switch) {
		command_str[idx] = 'W';
	}
	out_seq = "";
	out_seq.reserve(maxInsert + seq_len - maxDel + 1);

	for (long i = 0; i < seq_len; i++) {
		auto cmd = command_str.at(i);
		switch (cmd) {
		case 'I': {
			out_seq += randNucl();
			out_seq += seq->at(i);
			break;
		}
		case 'S': {
			out_seq += seq->at(i);
			break;
		}
		case 'D': {
			break;
		}
		case 'W': {
			out_seq += randNucl();
			break;
		}
		}
	}
}
