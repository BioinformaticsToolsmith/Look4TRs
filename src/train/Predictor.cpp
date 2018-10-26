/* -*- C++ -*-
 *
 * Predictor.cpp
 *
 * Author: Benjamin T. James, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 */

#include "Predictor.h"
#include "../matrix/Matrix.h"
#include <algorithm>
#include "../nonltr/KmerHashTable.h"
#include "../cluster/DivergencePoint.h"
#include "string.h"
#include <random>



std::default_random_engine generatorPredictor (0);  // minstd_rand0 is a standard linear_congruential_engine


using namespace nonltr;
template<class V>
void fill_table(KmerHashTable<unsigned long, V> &table,
	ChromosomeOneDigit *chrom, std::vector<V>& values) {

	const int k = table.getK();
	auto segment = chrom->getSegment();

	const char *seg_bases = chrom->getBase()->c_str();

	for (vector<int> *v : *segment) {
		int start = v->at(0);
		int end = v->at(1);
		table.wholesaleIncrement(seg_bases, start, end - k + 1);
	}

	std::vector<std::string> *keys = table.getKeys();
	for (std::string str : *keys) {
		values.push_back(table.valueOf(str.c_str()));
	}
	keys->clear();

	delete keys;
}

template<class T>
Point<T>* Predictor<T>::get_point(nonltr::ChromosomeOneDigit *chrom) {
	if (chrom == NULL) {
		return NULL;
	}

	KmerHashTable<unsigned long, T> table(k, 1);
	KmerHashTable<unsigned long, uint64_t> table_k1(1, 0);
	std::vector<T> values;
	vector<uint64_t> values_k1;

	values.clear();
	fill_table<T>(table, chrom, values);
	fill_table<uint64_t>(table_k1, chrom, values_k1);
	Point<T> *p = new DivergencePoint<T>(values, chrom->size());

	p->set_1mers(values_k1);
	p->set_header(chrom->getHeader());
	p->set_length(chrom->getBase()->length());
	p->set_data_str(*chrom->getBase());

	DivergencePoint<T>* q = dynamic_cast<DivergencePoint<T>*>(p);
	const auto N = q->points.size();
	double aq = (double) q->getPseudoMagnitude() / N;
	double sq = 0;
	for (auto i = 0; i < N; i++) {
		double qdiff = q->points[i] - aq;
		sq += qdiff * qdiff;
	}
	sq = sqrt(sq / N);
	q->set_stddev(sq);
	p->set_id(seq_num++);
	return p;
}

template<class T>
void Predictor<T>::save(std::string file) {
	std::ofstream out(file);
	out << "k: " << k << endl;
	out << "mode: " << (unsigned int) mode << endl;
	out << "max_features: " << max_num_feat << endl;
	out << "ID: " << id << endl;
	if (mode & PRED_MODE_CLASS) {
		write_to(out, feat_c, c_glm);
	}
	if (mode & PRED_MODE_REGR) {
		write_to(out, feat_r, r_glm);
	}
}

template<class T>
Predictor<T>::Predictor(const std::string filename) {
	std::ifstream in(filename);
	std::string buf;
	unsigned mode_ = 0;
	in >> buf >> k;
	cout << buf << k << endl;
	in >> buf >> mode_;
	mode = mode_;
	cout << buf << mode << endl;
	in >> buf >> max_num_feat;
	cout << buf << max_num_feat << endl;
	in >> buf >> id;
	cout << buf << id << endl;
	is_trained = true;
	if (mode & PRED_MODE_CLASS) {
		auto pr = read_from(in, k);
		c_glm = pr.first;
		feat_c = pr.second;
	}
	if (mode & PRED_MODE_REGR) {
		auto pr = read_from(in, k);
		r_glm = pr.first;
		feat_r = pr.second;
	}
}

template<class T>
void Predictor<T>::write_to(std::ofstream &out, Feature<T>* feat,
	matrix::GLM glm) {
	auto combos = feat->get_combos();
	auto lookup = feat->get_lookup();
	auto mins = feat->get_mins();
	auto maxs = feat->get_maxs();
	out << std::endl << "n_combos: " << combos.size() << std::endl;
	out << glm.get_weights().get(0, 0) << endl;
	for (int j = 0; j < combos.size(); j++) {
		auto cmb = combos[j];
		unsigned int val = 0;
		uint64_t flags = 0;
		for (auto i : cmb.second) {
			flags |= lookup[i];
		}
		switch (cmb.first) {
			case Combo::xy:
			val = 0;
			break;
			case Combo::xy2:
			val = 1;
			break;
			case Combo::x2y:
			val = 2;
			break;
			case Combo::x2y2:
			val = 3;
			break;
		}
		out << val << " ";
		out << flags << " ";
		out << glm.get_weights().get(j + 1, 0) << std::endl;
	}
	out << std::endl << "n_singles: " << lookup.size() << std::endl;
	for (int j = 0; j < lookup.size(); j++) {
		out << lookup[j] << " ";
		out << mins[j] << " ";
		out << maxs[j] << std::endl;
	}
}

template<class T>
pair<matrix::GLM, Feature<T>*> Predictor<T>::read_from(std::ifstream& in,
	int k_) {
	matrix::GLM glm;
	int c_num_raw_feat, c_num_combos;
	Feature<T> *feat = new Feature<T>(k_);
	std::string buf;
	in >> buf >> c_num_combos;
	cout << buf << "\"" << c_num_combos << "\"" << endl;
	matrix::Matrix weights(c_num_combos + 1, 1);
	double d_;
	in >> d_;
	weights.set(0, 0, d_);
	for (int i = 0; i < c_num_combos; i++) {
		int cmb;
		in >> cmb;
		cout << (int) cmb << endl;
		uint64_t flags;
		in >> flags;
		cout << flags << endl;
		double d;
		in >> d;
		cout << "[" << 0 << "," << i << "] " << d << endl;
		weights.set(i + 1, 0, d); //push_back(d);
		Combo cmb_ = Combo::xy;
		switch (cmb) {
			case 0:
			cmb_ = Combo::xy;
			break;
			case 1:
			cmb_ = Combo::xy2;
			break;
			case 2:
			cmb_ = Combo::x2y;
			break;
			case 3:
			cmb_ = Combo::x2y2;
			break;
			default:
			cerr << "error reading weights file" << endl;
			break;
		}
		feat->add_feature(flags, cmb_);
	}

	in >> buf >> c_num_raw_feat;
	cout << buf << "\"" << c_num_raw_feat << "\"" << endl;
	for (int i = 0; i < c_num_raw_feat; i++) {
		uint64_t single_flag;
		double min_, max_;
		in >> single_flag;
		cout << single_flag << endl;
		in >> min_;
		cout << min_ << endl;
		in >> max_;
		cout << max_ << endl;
		feat->set_normal(single_flag, min_, max_);
	}
	feat->finalize();
	glm.load(weights);
	return {glm, feat};
}

template<class T>
void Predictor<T>::add_feats(std::vector<std::pair<uint64_t, Combo> >& vec,
	uint64_t feat_flags) {
	for (uint64_t i = 1; i <= feat_flags; i *= 2) {
		if ((i & feat_flags) == 0) {
			continue;
		}
		for (uint64_t j = 1; j <= i; j *= 2) {
			if ((j & feat_flags) == 0) {
				continue;
			}
			vec.emplace_back(i | j, Combo::xy);
			vec.emplace_back(i | j, Combo::x2y2);
			if (i != j) {
				vec.emplace_back(i | j, Combo::x2y);
				vec.emplace_back(i | j, Combo::xy2);
			}
		}
	}
}

template<class T>
double Predictor<T>::similarity(Point<T>* a, Point<T>* b) {
	if (!is_trained) {
		double d = Selector<T>::align(a, b);
		pra<T> pr;
		pr.first = a->clone();
		pr.second = b->clone();
		pr.val = d;
		if (!is_training) {
			omp_set_lock(&lock);
			if (training.size() < testing.size()
				&& training.size() < threshold) {
				training.push_back(pr);
		} else if (training.size() >= testing.size()
			&& testing.size() < threshold) {
			testing.push_back(pr);
		} else {
			is_training = true;
			train();
			is_training = false;
		}
		omp_unset_lock(&lock);
	}
	return d;
} else {
	return predict(a, b);
}
}

template<class T>
bool Predictor<T>::close(Point<T> *a, Point<T> *b) {
	if (!is_trained) {
		double d = Selector<T>::align(a, b);
		pra<T> pr;
		pr.first = a->clone();
		pr.second = b->clone();
		pr.val = d;
		if (!is_training) {
			omp_set_lock(&lock);
			if (training.size() < testing.size()
				&& training.size() < threshold) {
				training.push_back(pr);
		} else if (training.size() >= testing.size()
			&& testing.size() < threshold) {
			testing.push_back(pr);
		} else {
			is_training = true;
			train();
			is_training = false;
		}
		omp_unset_lock(&lock);
	}
	return d;
} else {
	bool val = p_close(a, b);
	if ((mode & PRED_MODE_REGR) && val) {

	}
	return val;
}
}

template<class T>
double Predictor<T>::p_predict(Point<T>* a, Point<T>* b) {
	auto cache = feat_r->compute(*a, *b);
	auto weights = r_glm.get_weights();
	double sum = weights.get(0, 0);
	for (int col = 0; col < feat_r->size(); col++) {
		double val = (*feat_r)(col, cache);
		sum += weights.get(col + 1, 0) * val;
	}
	if (sum < 0) {
		sum = 0;
	} else if (sum > 1) {
		sum = 1;
	}
	return sum;
}
template<class T>
double Predictor<T>::predict(Point<T>* a, Point<T>* b) {
	if ((mode & PRED_MODE_CLASS) && !p_close(a, b)) {
		return 0;
	}
	return p_predict(a, b);
}

template<class T>
bool Predictor<T>::p_close(Point<T>* a, Point<T>* b) {
	auto weights = c_glm.get_weights();
	double sum = weights.get(0, 0);
	auto cache = feat_c->compute(*a, *b);
	for (int col = 1; col < weights.getNumRow(); col++) {
		double d = (*feat_c)(col - 1, cache);
		sum += weights.get(col, 0) * d;
	}
	return sum > 0.0;
}

template<class T>
std::pair<matrix::Matrix, matrix::Matrix> generate_feat_mat(
	const vector<pra<T> > &data, Feature<T>& feat, double cutoff) {
	bool classify = (cutoff >= 0);
	int nrows = data.size();
	int ncols = feat.size() + 1;
	matrix::Matrix feat_mat(nrows, ncols);
	matrix::Matrix labels(nrows, 1);

	for (int row = 0; row < data.size(); row++) {
		auto kv = data.at(row);
		vector<double> cache;
		{
			cache = feat.compute(*kv.first, *kv.second);
		}
		feat_mat.set(row, 0, 1);
		if (classify) {
			labels.set(row, 0, kv.val >= cutoff ? 1 : -1);
		} else {
			labels.set(row, 0, kv.val);
		}
		for (int col = 1; col < ncols; col++) {
			double val = feat(col - 1, cache);
			feat_mat.set(row, col, val);
		}
	}
	return std::make_pair(feat_mat, labels);
}

template<class T>
void Predictor<T>::train(const vector<Point<T> *> &points, size_t sample_size,
	size_t max_pts_from_one) {
	if (is_trained) {
		return;
	}
	selector = new Selector<T>(points, sample_size, max_pts_from_one);
	selector->select(id);
	auto tr = selector->get_training();
	training.insert(training.end(), tr.first.begin(), tr.first.end());
	training.insert(training.end(), tr.second.begin(), tr.second.end());
	auto te = selector->get_testing();
	testing.insert(testing.end(), te.first.begin(), te.first.end());
	testing.insert(testing.end(), te.second.begin(), te.second.end());
	train();
}
template<class T>
std::pair<double, matrix::GLM> regression_train(const vector<pra<T> > &data,
	Feature<T>& feat) {
	auto pr = generate_feat_mat(data, feat, -1);
	matrix::GLM glm;
	glm.train(pr.first, pr.second);
	auto result1 = pr.first * glm.get_weights();
	auto diff1 = result1 - pr.second;
	double sum = 0;
	for (int i = 0; i < diff1.getNumRow(); i++) {
		sum += fabs(diff1.get(i, 0));
	}
	sum /= diff1.getNumRow();
	return {sum, glm};
}

template<class T>
std::pair<std::tuple<double, double, double>, matrix::GLM> class_train(
	vector<pra<T> > &data, Feature<T>& feat, double cutoff) {
	auto pr = generate_feat_mat(data, feat, cutoff);
	matrix::GLM glm;
	glm.train(pr.first, pr.second);
	matrix::Matrix p = glm.predict(pr.first);
	for (int row = 0; row < p.getNumRow(); row++) {
		if (p.get(row, 0) == 0) {
			p.set(row, 0, -1);
		}
	}

	auto stats = glm.accuracy(pr.second, p);
	double acc = get<0>(stats);
	return {stats, glm};
}

template<class T>
double regression_test(const vector<pra<T> >& data, Feature<T>& feat,
	const matrix::GLM& glm) {
	auto pr = generate_feat_mat(data, feat, -1);
	auto result1 = pr.first * glm.get_weights();
	auto diff1 = result1 - pr.second;
	double sum = 0;
	for (int i = 0; i < diff1.getNumRow(); i++) {
		sum += fabs(diff1.get(i, 0));
	}
	sum /= diff1.getNumRow();
	return sum;
}

template<class T>
tuple<double, double, double> class_test(const vector<pra<T> >& data,
	Feature<T>& feat, const matrix::GLM& glm, double cutoff) {
	auto pr = generate_feat_mat(data, feat, cutoff);
	matrix::Matrix p = glm.predict(pr.first);
	for (int row = 0; row < p.getNumRow(); row++) {
		if (p.get(row, 0) == 0) {
			p.set(row, 0, -1);
		}
	}
	return glm.accuracy(pr.second, p);
}

template<class T>
void Predictor<T>::filter() {
	training.clear();
	testing.clear();

	auto tr = selector->get_training();
	size_t tr_size = 0.1 * tr.first.size();
	training.insert(training.end(), tr.first.begin(), tr.first.end());
	std::shuffle(tr.second.begin(), tr.second.end(), generatorPredictor);
	for (size_t i = 0; i < tr_size; i++) {
		training.push_back(tr.second[i]);
	}

	auto te = selector->get_training();
	size_t te_size = 0.1 * te.first.size();
	testing.insert(testing.end(), te.first.begin(), te.first.end());
	std::shuffle(te.second.begin(), te.second.end(), generatorPredictor);
	for (size_t i = 0; i < te_size; i++) {
		testing.push_back(te.second[i]);
	}
}

template<class T>
void balance(std::vector<pra<T> > &vec, double id) {
	std::vector<pra<T> > positive, negative;
	for (auto p : vec) {
		if (p.val > id) {
			positive.push_back(p);
		} else {
			negative.push_back(p);
		}
	}
	uint64_t tr_size = std::min(positive.size(), negative.size());
	vec.clear();
	for (uint64_t i = 0; i < tr_size; i++) {
		vec.push_back(positive[i]);
		vec.push_back(negative[i]);
	}
	// cout << "Positive: " << tr_size << " Negative: " << tr_size << endl;
}
template<class T>
void Predictor<T>::train() {
	uint64_t pos = 0;
	for (auto p : training) {
		if (p.val > id) {
			pos++;
		}
	}

	/*
	cout << "Training: Positive: " << pos << " Negative: "
			<< training.size() - pos << endl;
	*/

			uint64_t tr_size = std::min(pos, training.size() - pos);

			pos = 0;
			for (auto p : testing) {
				if (p.val > id) {
					pos++;
				}
			}
	/*
	cout << "Testing: Positive: " << pos << " Negative: "
			<< testing.size() - pos << endl;
	cout << "Training: ";
	*/

	balance(training, id);

	// cout << "Testing: ";

	balance(testing, id);
	Feature<T> feat(k);
	feat.set_save(true);
	if (mode & PRED_MODE_CLASS) {
		train_class(&feat);
		if (mode & PRED_MODE_REGR) {
			if (selector == NULL) {
				auto func = [&](pra<T> a) {return a.val < id;};
				training.erase(
					std::remove_if(training.begin(), training.end(), func),
					training.end());
				testing.erase(
					std::remove_if(training.begin(), training.end(), func),
					training.end());
			} else {
				filter();
			}
		}
	}
	if (mode & PRED_MODE_REGR) {
		train_regr(&feat);
	}
	feat.set_save(false);
	training.clear();
	testing.clear();
	possible_feats.clear();
	is_trained = true;
}

template<class T>
void Predictor<T>::train_class(Feature<T>* feat) {
	auto c_size = feat->get_combos().size();
	for (int i = 0; i < c_size; i++) {
		feat->remove_feature();
	}
	vector<uintmax_t> used_list;
	std::tuple<double, double, double> abs_best_acc, abs_best_tr_acc;
	get<0>(abs_best_acc) = 0;
	for (auto num_feat = 1; num_feat <= max_num_feat; num_feat++) {
		std::tuple<double, double, double> best_class_acc = abs_best_acc;
		std::tuple<double, double, double> best_class_tr_acc;
		uintmax_t best_idx = -1, cur_idx = 1;
		auto best_class_feat = possible_feats.front();
		for (uint64_t i = 0; i < possible_feats.size(); i++) {
			if (std::find(used_list.begin(), used_list.end(), i)
				!= used_list.end()) {
				continue;
		}
		auto rfeat = possible_feats[i];
		feat->add_feature(rfeat.first, rfeat.second);
		feat->normalize(training);
		feat->finalize();
		auto pr = class_train(training, *feat, id);
		auto class_ac = class_test(testing, *feat, pr.second, id);
		feat->remove_feature();

		if (get<0>(class_ac) > get<0>(best_class_acc)) {
			best_class_acc = class_ac;
			best_class_feat = rfeat;
			best_class_tr_acc = pr.first;
			best_idx = i;
		}
	}
	if (get<0>(best_class_acc) > get<0>(abs_best_acc)) {
		feat->add_feature(best_class_feat.first, best_class_feat.second);
		feat->normalize(training);
		feat->finalize();
		abs_best_acc = best_class_acc;
		abs_best_tr_acc = best_class_tr_acc;
		used_list.push_back(best_idx);
	}
}
feat_c = new Feature<T>(*feat);
feat_c->set_save(false);
auto pr = class_train(training, *feat_c, id);
c_glm = pr.second;
cout << "Classification Training ACC: " << get<0>(pr.first) << " SENS: "
<< get<1>(pr.first) << " SPEC: " << get<2>(pr.first) << endl;
auto stats = class_test(testing, *feat_c, c_glm, id);
cout << "Classification Testing ACC: " << get<0>(stats) << " SENS: "
<< get<1>(stats) << " SPEC: " << get<2>(stats) << endl;
}
template<class T>
void Predictor<T>::train_regr(Feature<T>* feat) {
	auto c_size = feat->get_combos().size();
	for (int i = 0; i < c_size; i++) {
		feat->remove_feature();
	}
	vector<uintmax_t> used_list;
	double abs_best_regr = 1000000;
	for (auto num_feat = 1; num_feat <= max_num_feat; num_feat++) {
		double best_regr_err = abs_best_regr;
		uintmax_t best_idx = -1, cur_idx = 1;
		auto best_regr_feat = possible_feats.front();
		for (uint64_t i = 0; i < possible_feats.size(); i++) {
			if (std::find(used_list.begin(), used_list.end(), i)
				!= used_list.end()) {
				continue;
		}
		auto rfeat = possible_feats[i];
		feat->add_feature(rfeat.first, rfeat.second);
		feat->normalize(training);
		feat->finalize();
		auto pr = regression_train(training, *feat);
		double regr_mse = regression_test(testing, *feat, pr.second);
		feat->remove_feature();

		if (regr_mse < best_regr_err) {
			best_regr_err = regr_mse;
			best_regr_feat = rfeat;
			best_idx = i;
		}
	}
	if (best_regr_err < abs_best_regr) {
		feat->add_feature(best_regr_feat.first, best_regr_feat.second);
		feat->normalize(training);
		feat->finalize();
		abs_best_regr = best_regr_err;
		used_list.push_back(best_idx);

	}
}
feat_r = new Feature<T>(*feat);
feat_r->set_save(false);
auto pr = regression_train(training, *feat_r);
r_glm = pr.second;

double regr_mse = regression_test(testing, *feat_r, r_glm);
}

template class Predictor<uint8_t> ;
template class Predictor<uint16_t> ;
template class Predictor<uint32_t> ;
template class Predictor<uint64_t> ;
template class Predictor<int> ;
template class Predictor<double> ;
