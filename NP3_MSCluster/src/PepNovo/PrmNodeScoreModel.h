#ifndef __PRMNODESCOREMODEL_H__
#define __PRMNODESCOREMODEL_H__

#include "RegionalPrmNodeScoreModel.h"

class PrmNodeScoreModel {

public:

	PrmNodeScoreModel(): ind_was_initialized(false), config_(NULL) {};

	void read_score_model(Config *config, const char *name, bool silent_ind = false);
	
	void write_score_model(const char *name) const;

//	void train_score_model(const char *name, const FileManager& fm, 
//						   int charge=-1, int size_idx=-1, int region_idx=-1);

	void trainNodeScoreModels(void* allScoreModelsVoidPointer, const char *name, const SpectraAggregator& sa,
							  int specificCharge=-1, int specificSize=-1, int specificRegion=-1);


//	void learn_prm_normalizer_values(const FileManager& fm);

	bool read_prm_normalizer_values();

	void write_prm_normalizer_values() const;

	void normalize_prm_scores(PrmGraph &prm) const;

	// Scores a breakage using a Dancik-like scoring model (independent frag probs vs. random)
	void score_breakage(Spectrum *spec, Breakage *breakage, bool verbose=false) const;
	
	void score_peptide_node_combos(void *ptr, PrmGraph *prm, const Peptide& peptide ) const;

	
	void initial_combos_score(void *model_ptr, PrmGraph *prm) const;

	score_t get_node_combo_score(PrmGraph *prm, int node_idx, 
								 int in_edge_idx, int in_var_idx, 
								 int out_edge_idx, int out_var_idx) const;

	int get_max_score_model_charge() const { return RegionalPrmNodeScoreModels_.size(); }

	bool get_ind_was_initialized() const { return ind_was_initialized; }
private:

	bool ind_was_initialized;

	Config* config_;

	vector< vector< vector< RegionalPrmNodeScoreModel > > > RegionalPrmNodeScoreModels_;

	vector< vector< vector< score_t > > > regional_prm_normalizers;
};




#endif

