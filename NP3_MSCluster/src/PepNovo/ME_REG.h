#ifndef __ME_REG_H__
#define __ME_REG_H__

#include "PepNovo_includes.h"



typedef float value_t;

struct fval {
	fval() : f_idx(-1), val(0) {};
	fval(int _f_idx, value_t _val) : f_idx(_f_idx), val(_val) {};

	bool operator< (const fval& other) const
	{
		return f_idx<other.f_idx;
	}

	int f_idx;
	value_t val;
};



// The structure for a sample for a two class logistic regression 
// The first class has feature values, the other is always 0
struct ME_Regression_Sample {

	ME_Regression_Sample() : label(-1) , prop(-1), prop_score(NEG_INF), weight(1) { }

	void print(const char **feature_names=NULL) const;

	void remove_feature(int f_idx);

	value_t get_feature_value(int f_idx) const
	{
		int i;
		for (i=0; i<f_vals.size(); i++)
			if (f_vals[i].f_idx==f_idx)
				return f_vals[i].val;
		return 0;
	}

	void add_feature(int idx, float val)
	{
		f_vals.push_back(fval(idx,val));
	}

	void clear()
	{
		label=-1;
		prop=-1;
		prop_score=NEG_INF; 
		weight=1.0;
		f_vals.clear();
	}

	int label;            // the true class of the sample

	int prop;            // optional additional property that is not used in the model

	double prop_score; // optional additional property
	
	double weight;        // the weight of this sample

	vector<fval> f_vals;  // holds feature values
};



class ME_Regression_DataSet {
	public:
	ME_Regression_DataSet() : num_samples(0), num_features(0), total_weight(0), num_classes(2) { samples.clear(); }
	ME_Regression_DataSet(int exp_samples) : num_samples(0) , num_features(0), num_classes(2),
								total_weight(0)  { samples.reserve(exp_samples); }

	void clear(int num_classes =2);

	void add_sample(const ME_Regression_Sample& sam);

	void add_samples(const vector<ME_Regression_Sample>& samples);

	void randomly_remove_samples_with_activated_feature(int label, int feature_idx, float prob_remove = 0.9);

	void tally_samples(bool print = false);

	// checks that features are within bounds and values are ok (not nan)
	int check_samples(bool remove_bad);

	void print() const;

	void print_summary() const;

	// prints info on features (num non zero and p~(f) )
	void print_feature_summary(ostream& os = cout, const char ** feature_names = NULL) const;

	// Splits the dataset into two new mutually exclusive set, must supply
	// a vector with the indices of the test set, all the other samples
	// are sent to the training set.
	void split_set(ME_Regression_DataSet& training,ME_Regression_DataSet& test, 
					vector<int> test_idxs) const;

	// output dataset, if null goes to screnn
	void write_data_to_file(const char *file=NULL) const;

	// calibrates class weights so class 0 has the desired weight
	void calibrate_class_weights(double class_0_weight);

	// calculates for each feature the proportion of samples (weight)
	// for which the feature has a non-zero value
	void calc_feature_non_zero_weights(vector< vector<double> >& ratios,
									   vector< vector<double> >& avg_nz) const;

	// gives various statistics about a certain feature's values
	void report_feature_statistics(int f_idx, const char *name = NULL) const;


	// calibrates the weight of the entire samples so the new weight adds
	// up to *total_weight*
	void rescale_dataset_weights(double total_weight);

	// extracts all the samples of the given class and puts them in a new
	// dataset
	void extract_class_samples(int label, ME_Regression_DataSet& extract) const;

	// exctract samples that have a non-zero value for the given feature
	void extract_samples_with_activated_feature(int feature_idx,ME_Regression_DataSet& extract) const; 

	// adds the samples from the other dataset, and adjust weights
	void add_other_dataset_samples(const ME_Regression_DataSet& other);


	// changes the weights of the data points so according to the weights
	// given in ratios (must have k ratios for values of prop that
	// are assumed to be 0,1,..,k-1
	void rescale_dataset_according_to_prop(const vector<double>& ratios);

	// return all samples in the datatset that have a desired label
	void get_samples_with_label(int label, vector<int>& idxs) const;


	// removes features that have low occurrence counts (to avoid over-fitting
	// and model stability issues)
	void purge_low_count_features(int min_num_occurrences);

	// returns the relative weight in the class of a certain feature
	double get_relative_weight_of_feature(int label, int feature_idx) const;

	// sets the weights of samples in the class in such a way that the relative weight of samples
	// with non-zero values for the given feature is given in the relative_weight
	void scale_samples_to_feature_relative_weight(int label, int feature_idx, double relative_weight);


	void serial_scale(const vector<int>& feature_idxs);


	int max_label;     // maximum index of a class label, labels should be 0,..,k-1
	int num_classes;  // number of classes k in the data = max_label+1
	int num_samples;
	int num_features;
	double total_weight;

	vector<double> class_weights;
	

	vector<ME_Regression_Sample> samples;
};





class ME_Regression_Model {
public:
	ME_Regression_Model() : num_features(-1), num_classes(0), t_iterations(0),
				has_weights(false) { f_weights.clear(); }
	

	// trains model using CG - Logistic Regression
	// returns false if there was no convergence (numerical stability issues...)
	bool train_cg(const ME_Regression_DataSet& ds, int max_interations=100,
		double epsilon=1E-4, int reset_rounds =0 );


	// takes all the samples in the dataset with the right label and calculates the probs
	// sorts them. Let x be the probability at the desired percentile, and t be the
	// target probability. The function retutns y, s.t. x^y=t
	float calc_log_scaling_constant(int label, 
									const ME_Regression_DataSet& ds, 
									float target_prob) const;


	void report_exp_sums(const ME_Regression_Sample& sam) const;


	double p_y_given_x(int label,const ME_Regression_Sample& sam) const
	{
		int i;
		double e,sum_exp=0;

		for (i=0; i<sam.f_vals.size(); i++)
		{
			const int& f_idx = sam.f_vals[i].f_idx;
			sum_exp+=f_weights[f_idx]*sam.f_vals[i].val;
		}

		if (sum_exp<-20)
			sum_exp=-20;
		
		if (sum_exp>20)
			sum_exp=20;
			
		e=exp(sum_exp);
		if (label == 0)
			return (e/(1.0 + e));
		
		return (1.0 / (1.0 + e));
	}

	float get_sum_exp(const ME_Regression_Sample& sam) const
	{
		int i;
		float sum_exp=0;

		for (i=0; i<sam.f_vals.size(); i++)
		{
			const fval& f_val = sam.f_vals[i];
			if (f_val.f_idx<f_weights.size())
				sum_exp+=(f_weights[f_val.f_idx]*f_val.val);
		}
		return sum_exp;
	}



	// sets the weights to the model returns a constant probability p for the class 0 samples
	void set_weigts_for_const_prob(float p); 




	void print_ds_probs(const ME_Regression_DataSet& ds) const;


	void print_ds_histogram(const ME_Regression_DataSet& ds) const;


	double log_likelihood(const ME_Regression_DataSet& ds) const;

	void read_weights_line(char *buff);

	void write_regression_model(ostream& os = cout) const;

	void read_regression_model(istream& is);

	int get_num_features() const { return num_features; }
	void set_num_features(int n) { num_features = n; f_weights.resize(n); };
	int get_num_classes() const { return num_classes; }
	double get_weight(int w) const { return f_weights[w]; }
	void set_num_classes(int nc) { num_classes=nc; }

	bool get_has_weights() const { return has_weights; }


private:




	// LOGISTIC REGRESSION FUNCTIONS

		/************************************************************************
	// calculates p(y|x) for all x,y
	// first index in p is class index 0..k-1
	// second index in p is sample number 0..n-1
	*************************************************************************/
	void calc_p_for_all_samples(const ME_Regression_DataSet& ds, vector<double>& p) const;


	/*********************************************************************
	Calculates the vector for the first derivatives with the current lambdas
	This function assumes that the values of p and fp are already computed.
	**********************************************************************/
	void calc_neg_first_derivative(const ME_Regression_DataSet& ds, const vector<double>& p,
						           vector<double>& vals) const;

	/******************************************************************
	Calculates the matrix for the Hesian of the current lambdas.
	This function assumes that the values of p and fp are already computed.
	*******************************************************************/
	void calc_Hesian2(const ME_Regression_DataSet& ds, const vector<double>& p,
						   vector< vector<double> >& hes) const;


	int t_iterations; // number of interations in training cycle
	int num_features;
	int num_classes;

	vector<double> f_weights;  //  the weights lambda_i of the features

	bool has_weights;

};


#endif


