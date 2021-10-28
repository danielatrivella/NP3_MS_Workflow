#ifndef __EDGEMODEL_H__
#define __EDGEMODEL_H__

#include "Config.h"
#include "SpectraList.h"
#include "BasicDataStructs.h"


/*
typedef enum EdgeModelFields_EM {

	EM_CONST,

	EM_IND_N_N_TERM, EM_IND_N_C_TERM, EM_IND_N_Gap, EM_IND_N_Xle, EM_IND_N_Ala, EM_IND_N_Arg, EM_IND_N_Asn,
	EM_IND_N_Asp,    EM_IND_N_Cys,    EM_IND_N_Gln, EM_IND_N_Glu, EM_IND_N_Gly, EM_IND_N_His, EM_IND_N_Ile,
	EM_IND_N_Leu,    EM_IND_N_Lys,    EM_IND_N_Met, EM_IND_N_Phe, EM_IND_N_Pro, EM_IND_N_Ser, EM_IND_N_Thr,
	EM_IND_N_Trp,    EM_IND_N_Tyr,    EM_IND_N_Val,

	EM_IND_C_N_TERM, EM_IND_C_C_TERM, EM_IND_C_Gap, EM_IND_C_Xle, EM_IND_C_Ala, EM_IND_C_Arg, EM_IND_C_Asn,
	EM_IND_C_Asp,    EM_IND_C_Cys,    EM_IND_C_Gln, EM_IND_C_Glu, EM_IND_C_Gly, EM_IND_C_His, EM_IND_C_Ile,
	EM_IND_C_Leu,    EM_IND_C_Lys,    EM_IND_C_Met, EM_IND_C_Phe, EM_IND_C_Pro, EM_IND_C_Ser, EM_IND_C_Thr,
	EM_IND_C_Trp,    EM_IND_C_Tyr,    EM_IND_C_Val,

	EM_CAT2, EM_CAT4, EM_CAT6, EM_CAT8, EM_CAT10, EM_CAT12, EM_CAT14, EM_CAT16, EM_CAT18, EM_CAT20,

	EM_PROBLEMATIC_W,
	EM_PROBLEMATIC_Q,
	EM_PROBLEMATIC_N,
		
	EM_NUM_FEATURES

} EdgeModelFields_EM;
*/

class RegionalEdgeModel {
	friend class EdgeModel;
public:
	
	RegionalEdgeModel() : has_values(false) {}

	bool read_model(const char *path);

	bool write_model(const char *path) const;

	void set_params(const Config *config, int c, int s, int r);

	void train_regional_edge_model(void *model_ptr, const SpectraAggregator& sa, SpectraList& sl);
	
	// returns max variant score
	score_t score_edge_variants(MultiEdge& multi_edge) const;


private:

	bool has_values;

	mass_t one_over_tolerance;
	int charge;
	int size_idx;
	int region_idx;

	int num_states; // 2^(number of strong frags)
	vector<int> strong_frag_idxs;
	vector<int> strong_frag_charges;

	vector<float> log_odds_state;	// intesection of N and C states
	vector< vector< float > > log_odds_offset; // per state and offset bin
	
	vector<float> multi_log_odds_state;	// intesection of N and C states
	vector< vector< float > > multi_log_odds_offset; // per state and offset bin

	vector< vector< float > > log_odds_transfer; // between N state and C state

	
	int calc_tol_bin_idx(const Breakage* n_break, const Breakage *c_break,
						 mass_t exp_edge_mass) const;

	int calc_state(const Breakage *breakage) const
	{
		int state=0;
		int p=1;
		int i;
		for (i=0; i<strong_frag_idxs.size(); i++)
		{
			const int frag_idx = strong_frag_idxs[i];
			if (! breakage->is_frag_type_visible(frag_idx) ||
				  breakage->get_position_of_frag_idx(frag_idx)>=0)
				state+=p;
			p+=p;
		}
		return state;
	}

	int calc_intersec_state(int n_state, int c_state) const
	{
		int intersec=0;
		if ((n_state & 0x1) && (c_state & 0x1))
			intersec+=1;
		if ((n_state & 0x2) && (c_state & 0x2))
			intersec+=2;
		if ((n_state & 0x4) && (c_state & 0x4))
			intersec+=4;
		if ((n_state & 0x8) && (c_state & 0x8))
			intersec+=8;
		return intersec;
	}


};


class EdgeModel {
public:
	EdgeModel() : ind_was_initialized(false), config_(NULL) {};

	void init_edge_model_defaults();

	void read_edge_models(const Config *config, const char *model_name, bool silent_ind=false);

	void write_edge_models(const char *model_name);

	void train_all_edge_models(const SpectraAggregator& sa, void *model_ptr, int specific_charge=0);

	// assigns probs and scores to graph edges
	void score_graph_edges(PrmGraph& prm) const;

	static int get_tol_bin_idx(float tol_ratio)
	{
		int i;
		for (i=0; i<tol_limits.size(); i++)
			if (tol_ratio<tol_limits[i])
				return i;
		return (tol_limits.size()-1);
	}

	static int get_num_tol_limits()
	{
		return tol_limits.size();
	}

	bool get_ind_was_initialized() const { return ind_was_initialized; }

private:

	bool ind_was_initialized;

	const Config* config_;

	// parameters for all regional models (magic numbers)
	static float weight_single_state_score;
	static float weight_single_offset_score;
	static float weight_multi_state_score;
	static float weight_multi_offset_score;
	static float weight_transfer_score;
	static float weight_combo_score;
	static float multi_aa_penalty;
	static float bad_pair_penalty;
	static float problematic_pair_penalty;

	static vector<float> tol_limits;

	vector< vector< vector< RegionalEdgeModel * > > > regional_edge_models;

	vector< vector< float > > double_combo_scores; // N-aa / C-aa

	string make_regional_name(const char *model_name,int c, int s, int r) const;
};





#endif




