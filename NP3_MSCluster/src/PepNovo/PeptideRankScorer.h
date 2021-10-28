#ifndef __PEPTIDERANKSCORER_H__
#define __PEPTIDERANKSCORER_H__

#include "AnnotatedSpectrum.h"
#include "PeptideComp.h"
#include "PeakRankModel.h"
#include "RankBoost.h"
#include "PMCSQS.h"


typedef enum SOL_TYPES {SOL_CORRECT, SOL_INCORRECT_DB, SOL_INCORRECT_DENOVO, SOL_INCORRECT_DB_CROSS} SOL_TYPES;





struct PeptideSet {

	PeptideSet() : ssf(NULL), file_idx(-1), scan(-1), total_set_weight(1.0) {}

	void *ssf;
	int file_idx;
	int scan;
	float total_set_weight;

	PeptideSolution  correct_sol;
	vector<PeptideSolution> incorrect_sols;
};

struct scan_pair {
	scan_pair() : file_idx(-1), scan(-1) {}
	scan_pair(int _f, int _s) : file_idx(_f) , scan(_s) {};

	bool operator< (const scan_pair& other) const
	{
		return ((file_idx<other.file_idx) ||
				 (file_idx==other.file_idx && scan<other.scan));
	}

	int file_idx;
	int scan;
};

typedef map< scan_pair, PeptideSet , less<scan_pair> > PeptideSetMap;





bool compare_cut_lists(const mass_t tolerance,
					   const vector<mass_t>& a_masses,
					   const vector<mass_t>& b_masses);


struct ScalingFactor {
	bool operator< (const ScalingFactor& other) const
	{
		return (max_pm_with_19>other.max_pm_with_19);
	}

	ScalingFactor() : max_pm_with_19(POS_INF), score_shift(0), score_scale(1.0) {}

	mass_t max_pm_with_19;
	float score_shift;
	float score_scale;
};



class DeNovoPartitionModel {
	friend class PeptideRankScorer;
public:

	DeNovoPartitionModel() :charge(NEG_INF), size_idx(NEG_INF),
							ind_was_initialized(false), num_ppp_frags(NEG_INF),
							ppp_start_idx(NEG_INF), use_ppp_features(false), 
							combined_ppp_start_idx(NEG_INF), use_combined_ppp_features(false),
							use_pmc_features(false),  pmc_start_idx(NEG_INF),
							use_comp_features(false), comp_start_idx(NEG_INF),
							use_peak_offset_features(false), peak_offset_start_idx(NEG_INF)
	{}


	bool read_denovo_part_model(const char *path, Config *config);

	void write_denovo_partition_model(int model_type, const char *path);

	void write_denovo_partition_model_header_to_strings(int model_type,
								vector<string>& header_strings) const;

	void set_shifts_and_scales_for_db(Config *config,
									  const RankBoostDataset& train_ds,
									  const vector<string>& peptide_strings);


	void init_features(int model_type, int _charge, int _size_idx, const vector<int>& ppp_frags, Config *config);

	const vector<string>& get_feature_names() const { return feature_names; }

	void simple_print_peak_pairs(
						  Config *config,
						  const vector<idx_weight_pair>& pair_idxs, 
						  vector<string>* peptide_strings,
						  const RankBoostDataset& ds,
						  int max_examples=-1,
						  ostream& os=cout) const;
		
private:

	int charge, size_idx;
	bool ind_was_initialized;

	vector<ScalingFactor> scaling_factors;


	int		num_ppp_frags;
	vector<int> ppp_frag_type_idxs;
	vector<FragmentType> ppp_fragments;

	bool	use_PTM_peak_features;
	int		PTM_peak_start_idx;

	bool	use_tryp_terminal_features;
	int		tryp_terminal_start_idx;

	bool	use_ppp_features;
	int		ppp_start_idx;

	bool    use_combined_ppp_features;
	int		combined_ppp_start_idx;
	
	bool	use_pmc_features;
	int		pmc_start_idx;

	bool	use_comp_features;
	int		comp_start_idx;

	bool	use_prm_features;
	int		prm_start_idx;

	bool	use_peak_offset_features;
	int		peak_offset_start_idx;

	bool	use_ann_peak_features;
	int		ann_peak_start_idx;

	bool    use_inten_balance_features;
	int		inten_balance_start_idx;

	vector<string> feature_names;

	RankBoostModel boost_model;


	void fill_peak_prediction_features(
			const PeptideSolution& sol, 
			const vector< vector<intensity_t> >& intens,
			const PeakRankModel *peak_model,
			RankBoostSample& rbs,
			int specific_size = -1) const;

	void fill_combined_peak_prediction_features(
			const PeptideSolution& sol,
			const vector< vector<intensity_t> >& intens,
			const PeakRankModel *peak_model,
			RankBoostSample& rbs,
			int specific_size = -1) const;

	void fill_pmcsqs_features(
		const PeptideSolution& sol,
		const vector<PmcSqsChargeRes>& res,
		const PMCSQS_Scorer *pmc_model,
		RankBoostSample& rbs) const;

	void fill_composition_features(const PeptideSolution& sol,
								   const Config *config,
								   PeptideCompAssigner *comp_assigner,
								   const SeqPath& path,
								   RankBoostSample& rbs) const;

	void fill_peak_offset_features(const Config *config,
								   const PeptideSolution& sol,
								   const vector< vector<mass_t> >& masses,
								   const vector< vector<intensity_t> >& intens,
								   RankBoostSample& rbs) const;

	void fill_ann_peak_features(const PeptideSolution& sol,
								const vector< vector<mass_t> >& masses,
								const vector< vector<intensity_t> >& intens,
								const AnnotatedSpectrum& as,
								RankBoostSample& rbs) const;

	void fill_inten_balance_features(const Config *config,
									 const PeptideSolution& sol, 
		    						 const SeqPath& path,
									 RankBoostSample& rbs) const;

	void fill_tryp_terminal_features(const PeptideSolution& sol, 
									 const SeqPath& sol_seq_path,
									 RankBoostSample& rbs) const;


	void fill_PTM_peak_features(const Config *config,
								const PeptideSolution& sol,
								const vector< vector<mass_t> >& masses,
								const vector< vector<intensity_t> >& intens,
								const AnnotatedSpectrum& as,
								RankBoostSample& rbs) const;

	void fill_prm_features(const PeptideSolution& sol, const SeqPath& seq_path, 
						   int model_type, RankBoostSample& main_rbs) const;

	const ScalingFactor& get_scaling_factor(mass_t pm_with_19) const
	{
		int i;
		for (i=0; i<scaling_factors.size(); i++)
			if (pm_with_19<scaling_factors[i].max_pm_with_19)
				break;

		if (i==scaling_factors.size())
		{
			cout << "Error: bad pm_with_19: " << pm_with_19 << endl;
			exit(1);
		}
		return scaling_factors[i];
	}
};	



struct weight_pair {
	weight_pair() : idx_corr(NEG_INF), idx_bad(NEG_INF), weight(0) {};
	weight_pair(int ic, int ib, float w) : idx_corr(ic), idx_bad(ib), weight(w) {};

	int idx_corr, idx_bad;
	float weight;
};


class PeptideRankScorer {
public:
	PeptideRankScorer() : model_type(NEG_INF), model_length(0) {}
	
	void set_model_length(int l) { model_length = l; }
	int  get_model_length() const { return model_length; }

	void set_model(void *allScoreModelsPtr) { allScoreModelsPtr_ = allScoreModelsPtr; }

	void set_type(int t) { model_type = t; }

	void read_denovo_rank_scorer_model(const char *path, string type_string, bool silent_ind = false);

	void write_denovo_rank_scorer_model(char *name);

	void rescore_inspect_results(char *spectra_file, char *inspect_res, char *new_res) const;

	void recalibrate_inspect_delta_scores(char *spectra_file, char *inspect_res, char *new_res) const;

	bool get_ind_part_model_was_initialized(int charge, int size_idx) const { 
		return (dnv_part_models[charge][size_idx] && dnv_part_models[charge][size_idx]->ind_was_initialized); }

	const DeNovoPartitionModel * get_dnv_part_model(int charge, int size_idx) const 
		{ return dnv_part_models[charge][size_idx]; }
								 

	// for complete de novo perdictions
	void train_partition_model_for_complete_sequences(
					const string& db_dir,
					const string& correct_dir,
					const string& denovo_dir,
					const string& mgf_list,
					char *report_dir,
					char *name,
					int charge,
					int size_idx,
					int max_num_rounds,
					float max_boost_ratio,
					int   max_num_samples = 200000,
					float ratio_pair_db = 0.3,
					float ratio_pair_denovo = 0.4,
					float ratio_pair_db_cross = 0.3,
					char  *rerank_path=NULL,
					int   rerank_depth = 2000);

	// for de novo predictions that might be incomplete
	void train_partial_denovo_partition_model(
					const string& mgf_list,
					char *report_dir,
					char *name,
					int charge,
					int size_idx,
					int max_num_rounds,
					float max_boost_ratio,
					int   max_num_samples = 200000,
					int   length_limit = 20,
					char  *rerank_path = NULL);


	
	

	void scoreCompleteSequences(const vector<PeptideSolution>& peptide_sols,
								AnnotatedSpectrum& as,
								vector<score_pair>& scores,
								int forced_size_idx=-1) const;



	void scoreDenovoSequences(const vector<SeqPath>& seqPaths,
							  AnnotatedSpectrum& as,
							  vector<score_pair>& scorePairs,
							  int forcedSizeIndex=-1,
							  int maximalIndexForRanking=-1) const;


	void scoreTagSequences(const vector<SeqPath>& seq_paths,
							 AnnotatedSpectrum& as,
							 vector<score_pair>& scores,
							 int forced_size_idx=-1) const;


	void set_model_type(int t) { model_type = t; }

	int  get_model_type() const { return model_type; }

//	PeakRankModel      *get_peak_prediction_model(int type) const { return peak_prediction_models[type]; }

//	AllScoreModels * get_model() const { return model; }

	void give_de_novo_and_peak_match_examples(
				const string& db_dir,
				const string& correct_dir,
				const string& denovo_dir,
				const string& mgf_list,
				const int charge,
				const int size_idx);


private:

	int model_type; // 0 - full de novo (db), 1 - partial de novo, 2- db score , 3 - tag

	int model_length; // 0 all lengths, otherwise specific length

	string dnv_model_name;

	static void* allScoreModelsPtr_;

	vector< vector<DeNovoPartitionModel *> > dnv_part_models;

	void init_tables(bool silent_ind = false);
	
	void create_training_data_for_complete_sequence_ranking(	
				const string& db_dir,
				const string& correct_dir,
				const string& denovo_dir,
				const string& mgf_list,
				const double train_ratio,
				const int max_num_pairs,
				const int charge,
				const int size_idx,
				RankBoostDataset& train_ds, 
				RankBoostDataset&test_ds,
				vector<string>* peptide_strings = NULL,
				char *test_scan_file = NULL,
				float ratio_db       = 0.4,
				float ratio_denovo   = 0.4,
				float ratio_db_cross = 0.2);

	void create_training_data_for_complete_denovo_ranking(	
				const string& db_dir,
				const string& correct_dir,
				const string& mgf_list,
				const double train_ratio,
				const int max_num_pairs,
				const int charge,
				const int size_idx,
				RankBoostDataset& train_ds, 
				RankBoostDataset& test_ds,
				vector<string>* peptide_strings = NULL,
				char *test_scan_file = NULL,
				float ratio_denovo   = 0.8,
				char *rerank_path	 = NULL,
				int   rerank_depth	 = 2000);

	void create_training_data_for_partial_denovo_ranking(
				const string& mgf_list,
				const double train_ratio,
				const int max_num_pairs,
				const int charge,
				const int size_idx,
				const float penalty_for_bad_aa,
				RankBoostDataset& train_ds, 
				RankBoostDataset&test_ds,
				char *test_scan_file = NULL,
				int   length_limit = 20);

	
	void fillCompletePeptideRankboostSample(
							   const PeptideSolution& sol,
							   AnnotatedSpectrum& as,
							   const vector<PmcSqsChargeRes>& res,
							   RankBoostSample& rbs,
							   int size_idx=-1) const;

	void fillDenovoPeptideRankboostSampleWithCombos(PeptideSolution& sol,
								 const SeqPath& path,
								 AnnotatedSpectrum& as,
								 const vector<PmcSqsChargeRes>& res,
								 RankBoostSample&		  main_rbs,
								 vector<RankBoostSample>& peak_prediction_combos,
							     int size_idx=-1) const;

	void fillTagRankboostSample(PeptideSolution& sol,
					  const SeqPath& path,
					  AnnotatedSpectrum& as,
					  RankBoostSample&	main_rbs,
					  int size_idx=-1) const;

	void select_sample_pairs(const vector<SeqPath>& solutions, const vector<int>& corr_idxs,
							 const vector<int>& bad_idxs, vector<weight_pair>& sample_pairs,
							 int num_pairs) const;

};

void create_complete_denovo_set_map(
						Config *config,
						const string& mgf_list,
						const string& db_dir,
						const string& correct_dir,
						const string& denovo_dir,
						int charge,
						int size_idx,
						PeptideSetMap& psm,
						vector<bool>& file_indicators);



/*void find_special_PTM_frags_using_offset_counts(
										   const string& PTM_label,
										   FileManager& fm,
										   const vector<SingleSpectrumFile *>& all_ssfs,
										   void *model,
										   int max_charge);
*/




void run_peak_benchmark(void *model, char *benchmark_file);



#endif




