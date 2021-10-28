#ifndef __REGIONALPRMNODESCOREMODEL_H__
#define __REGIONALPRMNODESCOREMODEL_H__

#include "RegularFragmentModel.h"
#include "StrongFragmentModel.h"

class AllScoreModels; // fwd dclr

/*****************************************************************************
Scoring model for a specific region (charge / size_idx / region_idx).
Contains a FragModel for each fragment used for scoring in that region.
******************************************************************************/
class RegionalPrmNodeScoreModel
{
	friend class AllScoreModels;
public:
	RegionalPrmNodeScoreModel() :  charge_(NEG_INF), sizeIndex_(NEG_INF), regionIndex_(NEG_INF),
								   num_strong_frags(0), num_regular_frags(0), config_(NULL),
								   was_initialized(false), has_all_breakage_models(false),
								   rand_prob(0), log_random(0), log_one_minus_random(0), 
								   missing_breakage_score(0) {};

	void init(const Config* config, int charge, int sizeIndex, int regionIndex);

	bool read_regional_score_model(const char *name, bool silent_ind);

	bool write_regional_score_model(const char *name) const;

//	bool train_regional_score_model(void *model, const char *name, const FileManager& fm);

	bool trainRegionalScoreModel(void *model, const char *name, const SpectraAggregator& sa);

	void calc_constant_element(Node& node,
							   Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage) const;

	score_t score_a_single_breakage_combo(
							   PrmGraph *prm,
							   Node& node, 
							   const Breakage *breakage,
							   BreakageInfo& info,
							   bool verbose=false) const;

	string make_model_file_name(const char *name) const;



	score_t get_frag_prob(int frag_type_idx) const
	{
		int i;
		for (i=0; i<frag_type_idxs.size(); i++)
			if (frag_type_idxs[i]==frag_type_idx)
				return frag_probs[i];
		return NEG_INF;
	}

	bool get_was_initialized() const { return was_initialized; }
	bool get_has_all_breakage_models() const { return has_all_breakage_models; }
	const vector<int>& get_frag_type_idxs() const { return frag_type_idxs; }
	const vector<score_t>& get_frag_inten_scores() const { return frag_inten_scores; }
	const vector<score_t>& get_frag_no_inten_scores() const { return frag_no_inten_scores; }

private:
	int charge_, sizeIndex_, regionIndex_;
	int num_strong_frags;
	int num_regular_frags;
	const Config *config_;

	bool was_initialized;
	bool has_all_breakage_models;

	vector<StrongFragmentModel>  strong_models;
	vector<RegularFragmentModel> regular_models;

	// weights on how much to emphasize the FragModel vs. the Dancik prob 
	vector<float> strong_inten_weights,  strong_no_inten_weights;
	vector<float> regular_inten_weights, regular_no_inten_weights;

	// holds the (1-weight) * frag prob to be weight*model prob for the acutal probability
	vector<float> strong_inten_danc_part,  strong_no_inten_danc_part;
	vector<float> regular_inten_danc_part, regular_no_inten_danc_part;

	// basic score probs
	vector<int>     frag_type_idxs;
	vector<score_t> frag_probs;
	score_t		    rand_prob;
	score_t log_random; 
	score_t log_one_minus_random; 

	vector<score_t> frag_inten_scores;
	vector<score_t> frag_no_inten_scores;
	score_t			missing_breakage_score;

	// TODO might have to return base class of FragmentModel and have 
	// RegularFragmentModel and StrongFragmentModel be derived from it
	void createTrainingSet(AllScoreModels *model,
							 const FragmentModel& frag_model,
							 const SpectraAggregator& sa,
							 ME_Regression_DataSet& inten_ds,
							 ME_Regression_DataSet& no_inten_ds) const;


};






#endif




