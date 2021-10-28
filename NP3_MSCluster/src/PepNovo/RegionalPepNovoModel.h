#ifndef __REGIONALPEPNOVOMODEL_H__
#define __REGIONALPEPNOVOMODEL_H__

#include "Fragmentation.h"
#include "FragProbTable.h"


// This is the basic model used for pepnovo
// it has uses discrete valued tables to calculate breakage probabilities
// can possibly also use additional models to calculate edge scores (aa bias etc).

struct Edge;

class RegionalPepNovoModel {

public:

	RegionalPepNovoModel() : charge(0), size_idx(-1), region_idx(-1), missed_cleavage_score(NEG_INF) {};

	void read_regional_model(Config *config, istream& is);

	void write_regional_model(ostream& os) const;

	score_t get_missed_cleavage_score() const { return missed_cleavage_score; }

	void calc_missed_cleavage_score();

	score_t calc_breakage_score(Breakage *breakage, bool verbose = false, Config *config=NULL) const;

	void print_table_names(const Config *config, ostream& os = cout) const;

	void convert_to_scores(const vector<double>& q_rand);

private:

	int charge;   
	int size_idx;
	int region_idx;

	vector<int> frag_type_idxs;


	score_t missed_cleavage_score; // score given when all fragments do not have a peak


	vector<FragProbTable> independent_frag_tables; /* basic Dancik type probabilities for fragments
											 to be used in case the parents of this frag
											 are not in the visible area. */

	vector<FragProbTable> single_breakage_tables; /* Probablities from a network that has
												    fragments only from the breakage being scored. */

};




#endif


