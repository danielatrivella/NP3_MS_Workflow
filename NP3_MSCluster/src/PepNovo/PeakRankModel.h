#ifndef __PEAKRANKMODEL_H__
#define __PEAKRANKMODEL_H__

#include "AnnotatedSpectrum.h"
#include "BasicDataStructs.h"
#include "RankBoost.h"


#define MOBILE 1
#define PARTIALLYMOBILE 2
#define NONMOBILE 3

struct SeqPath;

void convert_seq_path_to_peptide_soluition(const Config *config, 
										   const SeqPath& seq_path, 
										   PeptideSolution& sol);

void convert_seq_path_to_peptide_soluition_and_fill_in_aas(const Config *config, 
														   const Peptide& correct_pep,
														   const SeqPath& seq_path, 
														   PeptideSolution& sol);

void push_back_all_RHK_pairs(vector<string>& real_feature_names,
							 string prefix_label);

int get_proton_mobility(const Peptide& pep, int charge);

int get_proton_mobility(int charge, int num_arg, int num_his, int num_lys);

int calc_RKH_combo_idx (int r, int k, int h);

void find_ranks(const vector<intensity_t>& intens, vector<int>& ranks);

void normalize_intens(vector<intensity_t>& intens);

void calc_combined_peak_ranks(const vector< vector<float> >& intens, 
							  vector< vector<int> >& peak_ranks);

struct PeakScore {
	PeakScore() : cut_idx(NEG_INF), frag_idx(NEG_INF), rank_score(NEG_INF) {};
	int cut_idx;
	int frag_idx;
	float rank_score;
};

struct PeptidePeakPrediction {
	
	PeptidePeakPrediction() : num_frags(0), most_basic_missing_on_n(0), most_basic_missing_on_c(0) {}

	void print_ranks_vs_intens(const vector< vector<float> >& intens) const;

	void make_rank_tables(const vector< vector<float> >& intens,
		vector< vector<int> >& observed_ranks, vector< vector<int> >& predicted_ranks) const;

	void make_rank_tables_for_combined_peak_predictions(
		const vector< vector<float> >& intens,
		vector< vector<int> >& observed_ranks, 
		vector< vector<int> >& predicted_ranks) const;

	int num_frags;
	int most_basic_missing_on_n, most_basic_missing_on_c;
	vector<int> amino_acids;
	vector<int> frag_idxs;
	vector< vector<float> > rank_scores;
	vector< PeakScore> combined_peak_scores; // used for the combined peak models (feature_type_set >= 3)
};




struct PeakStart {
	PeakStart() : peptide_sample_idx(-1), peak_start_idx(-1), num_peaks(0) {};

	int peptide_sample_idx;
	int peak_start_idx;
	int num_peaks;
};


struct TrainingPeptide;
class  PeakRankModel;


/***************************************************************************
// model for a specific charge/size_idx/mobility
// has models for {set of fragments}
****************************************************************************/
class PartitionModel {
	friend class PeakRankModel;
public:

	PartitionModel() : partition_name("empty"), feature_set_type(NEG_INF), 
		charge(NEG_INF), size_idx(NEG_INF), mobility(NEG_INF), max_frag_idx(NEG_INF),
		num_features_per_frag(0) {};

	const vector<int>& get_fragment_type_idxs() const { return fragment_type_idxs; }

	
	// works for all models (both separate frags and the combined model)
	void calc_peptides_peaks_rank_scores(const PeakRankModel *prank,
						  const PeptideSolution& sol,
						  mass_t min_detected_mass, 
						  mass_t max_detected_mass,
						  PeptidePeakPrediction& ppp,
						  int   feature_set_type = 1,
						  const vector<int>* ptr_frag_type_idxs = NULL) const;

	
	int read_combined_partition_model(const string& path, Config *config,
							 int _charge, int _size_idx, int _mobility);
	

	int write_combined_partition_model(const string& path);



	void train_combined_partition_model(
								const PeakRankModel *prank, 
								char *sample_file_path,
								int	_charge,
								int  _size_idx,
								int  _mobility,
								int  num_frags = 3, 
								char *report_dir = NULL,
								int  max_rounds = -1,
								char *test_set = NULL,
								int	  test_peptide_length=-1,
								char *stop_signal_file = NULL,
								weight_t max_weight_ratio = 5.0);


	void fill_combined_peak_features(
								 const PeakRankModel *prank,
								 const  vector<int>& amino_acids,
								 const int    cut_idx,
								 const mass_t cut_mass,
								 const PeptideSolution& sol,
								 const FragmentType& frag,
								 const int   frag_pos_idx,
								 RankBoostSample& sample) const;

	void fill_combined_simple_peak_features(
								 const PeakRankModel *prank,
								 const  vector<int>& amino_acids,
								 const int    cut_idx,
								 const mass_t cut_mass,
								 const PeptideSolution& sol,
								 const FragmentType& frag,
								 const int   frag_pos_idx,
								 RankBoostSample& sample,
								 bool	verbose = false) const;

	void fill_combined_dnv_peak_features(
								 const PeakRankModel *prank,
								 const mass_t n_mass,   // this is where the possibly partial peptide starts
								 const mass_t c_mass,
								 const  vector<int>& amino_acids,
								 const int    cut_idx,
								 const mass_t cut_mass,
								 const PeptideSolution &sol,
								 const FragmentType& frag,
								 const int position_idx_in_model_fragment_type_idxs,
								 RankBoostSample& sample) const;


	void set_combined_feature_names_in_rankboost_model(const PeakRankModel *prank);
	void set_combined_simple_feature_names_in_rankboost_model(const PeakRankModel *prank);
	void set_combined_dnv_feature_names_in_rankboost_model(const PeakRankModel *prank);

	void train_partition_model(PeakRankModel *prank, 
							char *sample_file_path,
								int	_charge,
								int  _size_idx,
								int  _mobility,
								int frag_idx = -1, 
								char *report_dir = NULL,
								int  max_rounds = -1,
								char *test_set = NULL,
								int	  test_peptide_length=-1,
								char *stop_signal_file = NULL,
								weight_t max_weight_ratio = 5.0);


	int read_partition_model(const string& path, Config *config,
							 int _charge, int _size_idx, int _mobility);

	int write_partition_model(const string& path);





	void print_combined_peak_pairs(const vector<idx_weight_pair>& pair_idxs, 
						  const vector<TrainingPeptide>& tps,
						  const RankBoostDataset& ds,
						  PeakRankModel *prank,
						  int max_examples=-1,
						  ostream& os = cout) const;

	void simple_print_peak_pairs(const vector<idx_weight_pair>& pair_idxs, 
						  const vector<TrainingPeptide>& tps,
						  const RankBoostDataset& ds,
						  PeakRankModel *prank,
						  int frag,
						  int max_examples=-1,
						  ostream& os = cout) const;

	const string& get_partition_name() const { return partition_name; }

	void set_partition_name(const string& peak_rank_model_name, int charge,
							int size_idx, int mobility);

	void set_feature_set_type(int t) { feature_set_type = t; }

	void   print_model_stats() const;

private:

	string partition_name; // name_charge_sizeidx_mobility

	int feature_set_type;

	int charge;
	int size_idx;
	int mobility;

	int max_frag_idx;

	vector<int> fragment_type_idxs;
	vector<RankBoostModel> frag_models;

	int				num_features_per_frag;
	RankBoostModel	combined_frag_boost_model;
};




/***********************************************************************
Container class for all charge/size/mobility partition models
************************************************************************/
class PeakRankModel {
public:

	// 
	bool read_peak_rank_model(Config *_config, const char *name, bool silent_ind=false,
		int specific_charge=-1, int specific_size=-1, int specific_mobility=-1);

	// 
	void write_peak_rank_model(char *name, char *out_dir = NULL);


	void calc_peptide_predicted_scores(const PeptideSolution& sol,
									   PeptidePeakPrediction& ppp,
									   int specific_size = -1,
									   const vector<int>* ptr_frag_type_idxs = NULL) const;

	bool make_peak_prediction_table(
			const PeptideSolution& sol,
			const vector< vector<intensity_t> >& intens,
			int num_peaks) const;


	
	void   set_max_detected_mass(mass_t m) { max_detected_mass = m; }

	mass_t get_max_detected_mass() const { return max_detected_mass; }

	int		get_num_model_aas() const { return model_aa_labels.size(); }

	const vector<string>& get_model_aa_labels() const { return model_aa_labels; }

	mass_t calc_min_detected_mass(mass_t pm_with_19, int charge) const;

	void    set_config(Config *con) { config = con; }

	Config *get_config() const { return config; }

	const vector<string>& get_binary_names() const { return binary_feature_names; }

	const vector<string>& get_real_names()   const { return real_feature_names; }

	const vector<int>&	  get_real_feature_stage_idxs() const { return real_feature_stage_idxs; }

	PartitionModel * get_non_const_model_ptr(int charge, mass_t pm_with_19, int mobility) const
	{
		const int size_idx = get_size_group(charge, pm_with_19);
		return partition_models[charge][size_idx][mobility];
	}

	const PartitionModel * const get_model_ptr(int charge, mass_t pm_with_19, int mobility) const
	{
		const int size_idx = get_size_group(charge, pm_with_19);
		return partition_models[charge][size_idx][mobility];
	}

	const PartitionModel * const get_model_ptr(int charge, int size_idx, int mobility) const
	{
		return partition_models[charge][size_idx][mobility];
	}


	int get_size_group(int charge, mass_t pm_with_19) const
	{
		const vector<mass_t> & thresholds = size_thresholds[charge];
		int size_idx;
		for (size_idx=0; size_idx<thresholds.size(); size_idx++)
			if (pm_with_19<thresholds[size_idx])
				return size_idx;
		return thresholds.size();
	}

	int get_num_size_thresholds(int charge) const { return size_thresholds[charge].size(); }


	void fill_simple_peak_features(
								 const  vector<int>& amino_acids,
								 int    cut_idx,
								 mass_t cut_mass,
								 mass_t pm_with_19,
								 int	spec_charge,
								 const FragmentType& frag,
								 RankBoostSample& sample) const;


	void fill_advanced_peak_features(
								 const  vector<int>& amino_acids,
								 int    cut_idx,
								 mass_t cut_mass,
								 mass_t pm_with_19,
								 int	spec_charge,
								 const FragmentType& frag,
								 RankBoostSample& sample) const;


	void fill_partial_denovo_peak_features(
								 mass_t n_mass,   // this is wehere the possibly partial peptide starts
								 mass_t c_mass,
								 const  vector<int>& amino_acids,
								 int    cut_idx,
								 mass_t cut_mass,
								 mass_t pm_with_19,
								 int	spec_charge,
								 const FragmentType& frag,
								 int most_basic_on_n_side, // if the n side does not reach the terminal
								 int most_basic_on_c_side, // if the c side does not reach the terminal
								 RankBoostSample& sample) const;

	void set_simple_feature_names();

	void set_advanced_feature_names();


	void set_partial_denovo_feature_names();

	void partition_training_samples(const vector< vector<TrainingPeptide> >& all_tps,
									char *file_path_prefix = NULL,
									char *test_path_prefix = NULL,
									int   minimal_size = 4750,
									float prop_ts = 0.25) const;

	void read_training_peptides_into_rank_boost_dataset(
										int frag_type_idx,
										int spec_charge,
										const vector<TrainingPeptide>& sample_tps,
										RankBoostDataset& rank_ds,
										vector<float>& peak_intens,
										vector<PeakStart>& peak_starts,
										vector<float>& max_annotated_intens) const;

	void train_all_partition_models(
								int frag_fill_type,
								char *prefix_path,
								int	  charge=-1,
								int   size_idx=-1,
								int	  mobility=-1,
								int	  frag_type_idx=-1,
								char *report_dir = NULL,
								int   num_rounds = -1,
								char *test_set = NULL,
								int	  test_peptide_length=-1,
								char *stop_signal_file = NULL,
								weight_t max_weight_ratio = 5.0);

	void read_training_peptides_into_combined_rank_boost_dataset(
										int spec_charge,
										int size_idx,
										int mobility,
										const vector<TrainingPeptide>& sample_tps,
										RankBoostDataset& rank_ds,
										vector<float>& peak_intens,
										vector<PeakStart>& peak_starts,
										vector<int>& peak_frag_types) const;

	void train_all_combined_partition_models(
								int frag_fill_type,
								char *prefix_path,
								int	  charge=-1,
								int   size_idx=-1,
								int	  mobility=-1,
								int	  num_frags=3,
								char *report_dir = NULL,
								int   num_rounds = -1,
								char *test_set = NULL,
								int	  test_peptide_length=-1,
								char *stop_signal_file = NULL,
								weight_t max_weight_ratio = 5.0);

	void convert_aas_to_model_aas(const vector<int>& org_aas, vector<int>& model_aas) const;

	void set_binary_feature_names();

	void set_real_feature_names();

	void set_size_thresholds();

	void set_mass_detection_defaults();

	void init_peak_rank_model_with_defaults(Config *_config, char *name, int feature_type = 1);


	const string& get_peak_rank_model_name() const { return peak_rank_model_name; }

	void set_peak_rank_model_name(string name)    { peak_rank_model_name = name; }

	void list_all_model_idxs();

	void print_model_init_stats() const;

	const vector< vector<mass_t> >& get_size_thresholds() const { return size_thresholds; }

	bool has_intialized_models(int charge, int size_idx, int frag_idx) const;

	int get_feature_set_type() const { return feature_set_type; }

	const vector<int>& get_session_aas_to_model_aas() const { return session_aas_to_model_aas; }


private:

	int feature_set_type; // 0 - regular , 1 - advanced, 2 - partial de novo

	vector< vector< vector<PartitionModel *> > > partition_models; // charge / size_group / moility

	vector< vector<mass_t> > size_thresholds; // charge / threshes

	vector<string> binary_feature_names;

	vector<string> real_feature_names;

	vector<int>    real_feature_stage_idxs;

	vector<int>	   session_aas_to_model_aas; // conversion of the session aas to the model aas code

	vector<string> model_aa_labels;

	Config *config;

	mass_t max_detected_mass;

	string peak_rank_model_name;

	vector<mass_t> charge_min_mass_coefficients;


	void init_aa_defaults();
	void convert_session_aas_to_model_aas();

};






struct TrainingPeptide {

	TrainingPeptide() : n_mass(0),  best_n_removed(0), best_c_removed(0),
		length(0), charge(0), mobility(0), pm_with_19(0) {};

	void get_ranks_for_frag_idx(int frag_idx, vector<int>& ranks) const;

	void create_training_peptide(const PeakRankModel& rm, const AnnotatedSpectrum& as);

	void write_to_stream(ofstream& ofs) const;

	bool read_from_stream(ifstream& ifs);

	void print(Config *config, ostream& ofs = cout) const;

	float get_basicity_score() const; // made up score, higher means more basic amino acids

	int   get_frag_idx_pos(int frag_idx) const {
		int i;
		for (i=0; i<frag_idxs.size(); i++)
			if (frag_idxs[i]==frag_idx)
				return i;
		return -1;
	}


	mass_t n_mass;

	int best_n_removed, best_c_removed; // the strongest amino acids R>K>H>Other that was removed from
									    // each end because of the sequence being partial de novo
	int length;
	int charge;
	int mobility;
	mass_t pm_with_19;

	vector<int> amino_acids;
	vector<int> frag_idxs;         // these say what type of fragment's peaks appear in each row, row i
								   // doesn't have fragments of type i, but it has type frag_idxs[i], which
								   // is a bit confusing...
	vector< vector<float> > intens; // frag,cut idx
};


/*
void read_data_into_training_peptides(const FileManager& fm, Config *config, 
									  PeakRankModel& rm, vector<TrainingPeptide>& tps);
									  */


void read_training_peptides_from_file(char *file, vector<TrainingPeptide>& all_tps,
									  int num_tp = -1);

void convert_tps_to_partial_denovo(Config *config, 
								   vector<TrainingPeptide>& all_tps, 
								   int num_to_add = 0,
								   bool verbose = false);


void write_training_peptides_to_file(char *file, const vector<TrainingPeptide>& all_tps);


void convert_list_to_trianing_peptide_file(char *list, char *tp_file,
										  char *model=NULL, char *ptm_line=NULL);


void select_training_peptides(const vector<TrainingPeptide>& all_tps, 
							  vector<int>& selected_idxs,
							  int charge=0, 
							  int mobility=0, 
							  int min_length=-1, 
							  int max_length=100, 
							  mass_t min_pm_with_19=-1, 
							  mass_t max_pm_with_19=10000);


void generate_size_reports();

void create_training_sets();

void train_all();


#endif

