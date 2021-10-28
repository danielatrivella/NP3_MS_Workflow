#ifndef __ALL_SCORE_MODELS_H__
#define __ALL_SCORE_MODELS_H__

#include "ME_REG.h"
#include "BasicDataStructs.h"
#include "Spectrum.h"
#include "Config.h"
#include "BasicDataStructs.h"
#include "Spectrum.h"
#include "PMCSQS.h"
#include "PeptideComp.h"
#include "PrmNodeScoreModel.h"
#include "EdgeModel.h"
#include "AminoAcidProbs.h"
#include "CumulativeSeqProb.h"
#include "PeptideRankScorer.h"





class PrmGraph; // fwd declaration

struct TrainingStage {
	TrainingStage() : name(""), index(-1), indWasInitialized(false) {}
	TrainingStage(const string& s, int i) : name(s), index(i), indWasInitialized(false) {}
	void writeNameHeader() const;

	string name;
	int    index;
	bool   indWasInitialized;
};


/*****************************************************************************
This is a container class that holds all types of scoring models used for
the various denovo functions and stages:
PrmGraphNodeScorer - scores the nodes in the prm spectrum graph
PrmGraphEdgeScorer - scores the edges in the prm spectrum graph
.
.
.

******************************************************************************/
class AllScoreModels  {
public:

	AllScoreModels();

	string getModelName() const { return model_name; }

	void set_model_name(string _name) { model_name = _name; config_.set_model_name(_name); }

	const Config * get_config() const { return &config_; }
	Config* get_config() { return &config_; }

	void setTrainingStageNames();

	// This function performs the entire training process of the model.
	// Though it theoretically can be run sequencially on as a single process,
	// for many of the models the training process is time consuming and is
	// best done in perallel for specific stages/charges/sizes/regions (see 
	// documentation for a more complete explanation). To this end, the function
	// accomadates spefific entry and exit points at different stages (with the
	// optional parameters startTrianingStage and endTrainingStage, and allows
	// for training of specific models with the optional parameters specificChargeToTrain,
	// specificSizeToTrain, and specificRegionToTrain.
	// If the training includes spectrum quality scores (SQS), the function needs to be
	// given a path to the file of negative spectra (i.e., non-peptide containing).
	void trainModelsInStages(
		const char *newModelName, 
		const SpectraAggregator& sa,
		mass_t initialToleranceEstimate, 
		int startTrianingStage = 0,
		int endTrainingStage   = POS_INF,
		int specificChargeToTrain =-1, 
		int specificSizeToTrain   =-1, 
		int specificRegionToTrain =-1,
		const char *pathOfNegativeSetForSqsTraining = NULL);



	

	void test_pmc(char *specs_file, int charge, mass_t min_mass=0,
		mass_t max_mass = POS_INF)
	{
		pmcsqs.test_pmc(&config_,specs_file,charge,min_mass,max_mass);
	}

	

	void benchmark_sqs(char *list, char *anns)
	{
		pmcsqs.benchmark_sqs(&config_, list, anns);
	}


	

	float findBestMzAndCharge(const Config *config, const PeakList& pl, 
						   mass_t* mz1, int* charge1, float *prob1,
						   mass_t* mz2, int* charge2, float *prob2,
						   vector<PmcSqsChargeRes>* all_res = NULL)
	{
		return pmcsqs.computeBestMzAndChargeForSpectrum(&config_,pl,mz1,charge1,prob1,mz2,charge2,prob2,all_res);
	}



	void selectPrecursorMassesAndCharges(const Config *config, 
										 const PeakList& pl,
										 vector<mass_t>& precursorMassesWith19,
										 vector<int>&    charges,
										 vector<PmcSqsChargeRes>* allResults=NULL)
	{
		pmcsqs.selectPrecursorMassesAndCharges(config, pl, precursorMassesWith19, charges, allResults);
	}







	bool get_ind_pmcsqs_was_intialized() const { return (pmcsqs.getIndInitializedPmc() &&
		pmcsqs.getIndInitializedSqs()); }

	const PMCSQS_Scorer* get_pmcsqs_ptr() const { return &pmcsqs; }

	const AminoAcidProbs* get_amino_acid_probs_ptr() const { return &amino_acid_probs; }

	int get_aa_category(int num_aa, int *aas, bool n_term, bool c_term) const
	{
		const vector<int>& org_aa = config_.get_org_aa();
		int conv_aas[6];
		int i;
		for (i=0; i<num_aa; i++)
			conv_aas[i]=org_aa[aas[i]];

		return compAssigner_.get_aa_category(num_aa,conv_aas,n_term,c_term);
	}
	
	
	
	void read_model(const char* model_name, bool silent = false);

	void clone_charge_model(int source_charge, int target_charge);

	void write_model();



	

	// performs scoring on demand (if the combo was not previously scored, calculates
	// values, otherwise returns hashed value
	score_t get_node_combo_score(PrmGraph *prm, int node_idx, 
								 int in_edge_idx, int in_var_idx, 
								 int out_edge_idx, int out_var_idx) const
	{
		return prmNodeScoreModel_.get_node_combo_score(prm, node_idx, 
								 in_edge_idx, in_var_idx, out_edge_idx, out_var_idx);
	}
	// required Model functions
	void init_score_model() { return; }

	void init_model_for_scoring_spectrum(Spectrum *spec) { return; }


//	score_t get_missing_breakage_score(int charge, int size_idx, int region_idx) const
//	{
//		return this->RegionalPrmNodeScoreModels_[charge][size_idx][region_idx].missing_breakage_score;
//	}


	void score_graph_edges(PrmGraph& prm) const;

	

	int get_max_score_model_charge() const;

	PeptideCompAssigner& getPeptideCompositionAssigner() { return compAssigner_; }
	const PeptideCompAssigner& getPeptideCompositionAssigner() const { return compAssigner_; }

	PeakRankModel*& get_peak_prediction_model_ptr(int type)  { return peak_prediction_models[type]; }

	PeptideRankScorer*& get_rank_model_ptr(int type) { return peptideRankModels_[type]; }

	PeptideRankScorer*& get_rank_tag_model_ptr(int length) { return tagRankModels_[length]; }

	CumulativeSeqProbModel *get_cumulative_seq_prob_model_ptr(int tag_length) { return cumulativeSequenceModels_[tag_length]; }

	void read_rank_models(const char *name, bool silent_ind = false);

	void read_cum_seq_prob_models(const char *name, bool silent_ind = false);


	bool get_ind_score_model_was_initialized() const { return ind_was_initialized; }

	const PrmNodeScoreModel& getPrmNodeScoreModel() const { return prmNodeScoreModel_; }
	

protected:

	bool ind_was_initialized;

	string model_name;
	
	Config config_;

	PMCSQS_Scorer pmcsqs;

	PeptideCompAssigner compAssigner_;

	PrmNodeScoreModel prmNodeScoreModel_;

	EdgeModel edgeModel_;

	AminoAcidProbs amino_acid_probs;

	vector<TrainingStage> trainingStages_;

	PeakRankModel       *peak_prediction_models[8]; // for the 4 types of PeptideRankScorer

	PeptideRankScorer* peptideRankModels_[3]; // db score / partial denovo / full denovo
	PeptideRankScorer* tagRankModels_[10]; // model per tag length (typically only 3-6 will have values)

	CumulativeSeqProbModel* cumulativeSequenceModels_[10]; // model 0 is for de novo, the others are for tags


	bool itializeTrainingStages();
	bool checkModelsInitializationStatus(bool report=false);


};





#endif


