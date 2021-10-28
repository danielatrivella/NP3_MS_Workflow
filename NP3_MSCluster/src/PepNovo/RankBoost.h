#ifndef __RANKBOOST_H__
#define __RANKBOOST_H__

#include "PepNovo_includes.h"


struct idx_weight_pair {
	idx_weight_pair() : idx(-1), weight(0) {};
	idx_weight_pair(int _i, weight_t _w) : idx(_i), weight(_w) {};

	bool operator< (const idx_weight_pair& other) const
	{
		return weight>other.weight;
	}

	int idx;
	weight_t weight;
};


struct peak_rank_stat {
	peak_rank_stat() : true_rank(-1), avg_precited_rank(-1), sd_predicted_rank(-1), precent_predicted_correctly(-1) {};

	int true_rank;
	float avg_precited_rank;
	float sd_predicted_rank;
	float precent_predicted_correctly;
};

class RankBoostModel;

struct RankBoostSample {

	RankBoostSample() { clear(); }

	void clear()
	{
		groupIndex = (int)NEG_INF;
		rank_in_group = (int)NEG_INF;

		tag1=(int)NEG_INF;
		tag2=(int)NEG_INF;
		tag3=(int)NEG_INF;
		float_tag1=NEG_INF;

		binary_non_zero_idxs.clear();
		binary_novote_idxs.clear();
		real_active_idxs.clear();
		real_active_values.clear();
		real_novote_idxs.clear();
	}

	void add_real_feature(int idx, float val) { 
		real_active_idxs.push_back(idx);
		real_active_values.push_back(val);
	}

	bool get_feature_val(int f_idx, float *val) const;

	void print(ostream& os = cout) const;

	int groupIndex;     // for testing purposes
	int rank_in_group; // for testing purposes
	int tag1,tag2,tag3; // extra information that can be stored for debugging
	float float_tag1;
	
	vector<int>	binary_non_zero_idxs; // idxs of all binary variables that have 1
	vector<int> binary_novote_idxs; // idxs for which we should use the default

	vector<int>   real_active_idxs; //
	vector<float> real_active_values; // the real value given to the feature
	vector<int>	  real_novote_idxs; // use the default value with these, should be used
									// only if default value is not 0
};



struct SamplePairWeight {
	SamplePairWeight() : idx1(-1), idx2(-1), tag(-1), weight(0) {};
	SamplePairWeight(int _idx1, int _idx2, weight_t _weight) : idx1(_idx1), idx2(_idx2), weight(_weight) {};
	SamplePairWeight(int _idx1, int _idx2, int _t, weight_t _weight) : idx1(_idx1), idx2(_idx2), tag(_t), 
					weight(_weight) {};

	int   idx1,idx2;
	int	  tag;
	weight_t weight;
};


class RankBoostDataset {

public:

	RankBoostDataset() : num_groups(0), 
						 max_ratio_for_regular_update(10000.0), 
						 total_phi_weight(0) {};

	void     set_max_ratio_for_regular_update(weight_t w) { max_ratio_for_regular_update = w; }
	weight_t get_max_ratio_for_regular_update() const { return max_ratio_for_regular_update; }

	void compute_max_weights_for_normal_update(
		const vector<weight_t>& D0, vector<weight_t>& d) const;

	void compute_total_phi_weight();

	void compute_initial_distribution(vector<weight_t>& d) const;

	void initialize_potenital_lists();

	void initialize_binary_one_lists(int num_binary_features);

	void initialize_binary_ordered_phi_lists(const vector<string>* binary_names = NULL);

	void initialzie_real_feature_table(int num_real_features);

	void initialize_real_vote_lists(const RankBoostModel& rbm);

	void calc_potentials(const vector<weight_t>& D, vector<weight_t>& potentials) const;

	double get_total_phi_weight() const { return total_phi_weight; }

	const vector<RankBoostSample>& get_samples() const { return samples; }
	const RankBoostSample& get_sample(int i) const { return samples[i]; }

	void  add_sample(const RankBoostSample& rbs) { samples.push_back(rbs); }

	void add_to_phi_vector(int x0_idx, int x1_idx, weight_t w = 1.0) 
	{
		phi_support.push_back(SamplePairWeight(x0_idx,x1_idx,w));
	}

	void add_to_phi_vector(int x0_idx, int x1_idx, int tag, weight_t w = 1.0) 
	{
		phi_support.push_back(SamplePairWeight(x0_idx,x1_idx,tag,w));
	}


	int   get_num_samples() const { return samples.size(); }

	const vector<int>& get_binary_one_list(int i) const { return binary_one_lists[i]; }

	const vector< vector<int> >& get_real_vote_lists(int feature_idx) const { return real_vote_lists[feature_idx]; }

	vector<SamplePairWeight>& get_non_const_phi_support() { return phi_support; }

	const vector<SamplePairWeight>& get_phi_support() const { return phi_support; }
	
	const vector<int>& get_binary_pairs_ordered_correctly(int f_idx) const { return binary_pairs_ordered_correctly[f_idx]; }
	
	const vector<int>& get_binary_pairs_ordered_incorrectly(int f_idx) const { return binary_pairs_ordered_incorrectly[f_idx]; }

	double update_dsitribution_according_to_binary_feature(int binary_feature_idx, 
								weight_t alpha, vector<weight_t>& D, bool verbose = false) const;

	double update_distribution_according_to_real_feature( 
		int best_real_idx, float theta_start,  float theta_end, 
		int q_def, weight_t alpha,  vector<weight_t>& D, 
		vector<weight_t> *max_D_for_normal_updates = NULL,
		bool verbose = false) const;

	int  get_num_groups() const { return num_groups; }
	void set_num_groups(int ng) { num_groups = ng; }
	
	void clear();

private:

	int	   num_groups;

	weight_t max_ratio_for_regular_update; // if a sample pair passes X times its intial weight
										   // then the increase in weight is regularized to avoid 
										   // over-training for a small number of samples
	double total_phi_weight;

	vector<RankBoostSample>  samples;

	vector<SamplePairWeight> phi_support; // all pairs of samples with weight > 0 (i.e., idx2 should be ranked higher than idx1)
														   // we perform the regular update when we com to increase the weight
	vector< vector<int> > ahead_lists;    // holds for each sample x, the idxs of samples x' for which x is ahead, that is (x',x) in phi
	vector< vector<int> > behind_lists;   // holds for each sample x, the idxs of samples x' for which x is behind, that is (x,x') in phi
	vector< vector<int> > binary_one_lists; // holds for each feature the idxs of the samples that have a 1 at that position

	vector< vector< vector<int> > > real_vote_lists; // holds the lists of feature idxs that vote on a given feature
	vector<bool>		  ind_all_samples_vote; // holds indicators that state for each real feature if all samples vote on it

	vector< vector<int> > binary_pairs_ordered_correctly;   // indices of phi of pairs (x',x) for which  f(x)=1 and f(x')=0 (or no vote)
	vector< vector<int> > binary_pairs_ordered_incorrectly; // indices of phi of pairs (x',x) for which  f(x)=0 (or no vote) and f(x')=1

	vector< map<int,float> > real_feature_values; // holds a table featuresXsamples of values (NEG_INF == no vote)
};



class RankBoostModel {
public:
	RankBoostModel() : ind_was_initialized(false), 
					   num_binary_features(0), 
					   num_real_features(0),
					   current_round(0),
					   train_error(1),
					   test_error(1),
					   total_default_weight(0) {};

	bool get_ind_was_initialized() const { return ind_was_initialized; }

	void init_rankboost_model_feature_names(
							   const vector<string>& _binary_feature_names,
							   const vector<string>& _real_feature_names);

	void init_rankboost_model_for_training(
							   const RankBoostDataset& training_ds,
							   int	 min_num_active_samples_for_feature = 50,
							   int	 max_num_real_bins_for_real_feature = 50);

	double calc_training_rank_error(const RankBoostDataset& training_ds) const;

	double calc_prediction_error(const RankBoostDataset& rank_ds, 
								 vector<peak_rank_stat>& peak_stats,
								 int peptide_length = -1,
								 int *num_groups = NULL) const;
									
	weight_t calc_rank_score(const RankBoostSample& sample) const;

	weight_t calc_rank_score_with_details(
						const RankBoostSample& sample,
						vector<string>& feature_names,
						vector<float>&	feature_values,
						vector<float>&  feature_scores) const;

	weight_t calc_rank_score_with_details(
						const RankBoostSample& sample,
						vector<int>  &  feature_idxs,
						vector<float>&	feature_values,
						vector<float>&  feature_scores) const;


	weight_t calc_rank_score(const RankBoostSample& sample, 
						  const vector<int>& active_binary_feature_idxs,
						  const vector<int>& active_real_feature_idxs) const;

	bool train_rankboost_model(const RankBoostDataset& training_ds,
							   int   max_num_rounds,
							   vector<idx_weight_pair>* top_misclassified_pairs= NULL,
							   RankBoostDataset* test_set = NULL,
							   int tag3_filter_val = -1,
							   char *report_preifx = NULL,
							   char *stop_signal_file = NULL,
							   const vector<string>* model_header_strings = NULL);

	void print_top_misclassified_pairs(const RankBoostDataset& training_ds,
									   const vector<weight_t>& D,
									   const vector<weight_t>& org_D,
									   int num_top_pairs = 10,
									   ostream& os = cout) const;

	void get_top_misclassified_pairs(
								   const RankBoostDataset& training_ds,
								   const vector<weight_t>& D,
								   const vector<weight_t>& D0,
								   vector<idx_weight_pair>& pair_idxs,
								   int num_top_pairs = -1) const;

	void summarize_features(const vector<RankBoostSample>& samples, ostream& os=cout); 

	void summarize_features_pos_neg(const vector<RankBoostSample>& pos_samples, 
									const vector<RankBoostSample>& neg_samples);

	
	bool read_rankboost_model(istream& is);

	// compresses the weights and limits before writing
	void write_rankboost_model(ostream& os, bool ind_compress = true);

	
	int  get_num_bins_for_real_feature(int real_feature_idx) const { return real_limits[real_feature_idx].size()+1; }

	// returns the index of the value (weight) idx for a given value
	// this is the index j for which our value is less or equal to theta_j
	// since the last theta is assumed to be infinity, the index |limits|
	// is given if our value is larger than all supplied limits
	// normally, the index 0 cannot be given since no value should be less than
	// -infinity which is the value in limits[0]
	int  get_real_bin_idx_for_value(int real_feature_idx, float val) const
	{
		const vector<float>& limits = real_limits[real_feature_idx];
		if (limits.size() == 0)
			return 0;

		int i;
		for (i=1; i<limits.size(); i++)
			if (val<=limits[i])
				break;
		return i;
	}

	int get_num_real_features() const { return num_real_features; }
	int get_num_binary_features() const { return num_binary_features; }
	const vector<weight_t>& get_binary_weights() const { return binary_weights; }
	const vector< vector<weight_t> >& get_real_weights() const { return real_weights; }
	const vector< vector<float> >&    get_real_limits() const  { return real_limits; }
	const vector<weight_t>&           get_real_default_weights() const { return real_default_weights; }
	const vector<string>* get_ptr_to_binary_feature_names() const { return &binary_feature_names; }
	const vector<string>* get_ptr_to_real_feature_names() const { return &real_feature_names; }
	const vector<string>& get_binary_feature_names() const { return binary_feature_names; }
	const vector<string>& get_real_feature_names() const { return real_feature_names; }
	const vector<int>& get_binary_update_counts() const { return binary_update_counts; }
	const vector<int>& get_real_update_counts() const { return real_update_counts; }

	void  set_real_feature_stage_idxs(const vector<int>& stage_idxs) { real_feature_stage_idxs = stage_idxs; }

	void ouput_ranked_feature_list( ostream& os = cout) const;

	void ouput_importance_ranked_feature_list( const RankBoostDataset& training_ds, 
											   ostream& os = cout,
											   int only_fidx = NEG_INF,
											   int round_idx = NEG_INF) ;


	void list_feature_vals_according_to_score(vector<RankBoostSample>& sams) const;


private:

	bool ind_was_initialized;

	int num_binary_features;
	int num_real_features;

	int	   current_round;
	double train_error, test_error;

	weight_t total_default_weight;

	vector<bool>	  ind_active_binary_feature;
	vector<weight_t>  binary_weights;
	vector<string>    binary_feature_names;
	
	vector<bool>	ind_active_real_feature;
	vector< vector<weight_t> > real_weights;  // feature / weights ( |limits|+1 values)
	vector< vector<float> > real_limits;
	vector<weight_t>        real_default_weights; // given to features that don't vote on a sample

	vector<string>			real_feature_names;
	vector<int>				real_feature_stage_idxs; // if set, these values are used to determine
													 // when to enter a feature into the training

	vector<int> non_zero_binary_idxs;
	vector<int> non_zero_real_idxs;

	// used in training
	vector<int> real_first_updates;
	vector<int> real_update_counts;   
	vector<int> binary_update_counts;

	// the values of the best round
	int best_round_idx;
	double best_train_error, best_test_error;
	double best_total_default_weight;
	vector<bool>	  best_ind_active_binary_feature;
	vector<weight_t>  best_binary_weights;
	vector<bool>			   best_ind_active_real_feature;
	vector< vector<weight_t> > best_real_weights;  // feature / weights ( |limits|+1 values)
	vector< vector<float> >	   best_real_limits;
	vector<weight_t>           best_real_default_weights; // given to features that don't vote on a sample
	vector<int> best_real_update_counts; 
	vector<int> best_binary_update_counts;
	vector<int> best_non_zero_binary_idxs;
	vector<int> best_non_zero_real_idxs;


	void select_active_features(const RankBoostDataset& training_ds,
								int min_sample_count = 50);

	void set_real_limits(const vector<RankBoostSample>& samples,
						 int min_bin_size = 50,
						 int max_num_bins = 200,
						 bool clear_weights = true);

	void binary_weak_learn(const vector<weight_t>& potentials,
						   const RankBoostDataset& rank_ds,
						   int&  best_feature_idx,
						   double& r_value,
						   bool only_non_zero_features = false) const;
	
	void real_weak_learn(const vector<weight_t>& potentials,
						 const RankBoostDataset& rank_ds,
						 int&  best_feature_idx,
						 int&  theta_idx,
						 int&  q_def,
						 double& r_value,
						 bool  only_non_zero_features = false) const;

	void real_weak_learn_double_theta(
						 const vector<weight_t>& potentials,
						 const RankBoostDataset& rank_ds,
						 int&  best_feature_idx,
						 int&  theta_start_idx, // x > theta
						 int&  theta_end_idx, // x<= theta
						 int&  q_def,
						 double& r_value,
						 bool  only_non_zero_features = false) const;


	void update_model_weights_for_binary_feature(int best_binary_idx, weight_t alpha);

	void update_model_weights_for_real_feature(weight_t alpha, int best_real_idx, 
								int q_def, int theata_idx_start, int theta_idx_end=-1);

	// compresses the number of limits and weights to represent steps in the weight
	// function without repeatition of the same weights in different limits
	void compress_real_limits_and_weights();

	void remove_default_weights();

	void set_best_model_parameters_to_current_parameters();

	void set_current_model_parameters_to_best_parameters();
};










#endif


