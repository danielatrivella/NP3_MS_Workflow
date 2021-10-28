#include "RankBoost.h"
#include "../Common/auxfun.h"


bool RankBoostModel::read_rankboost_model(istream& is)
{
	num_binary_features = 0;
	num_real_features =0;
	int num_non_zero_real_features=0;

	is >> num_binary_features;
	if (num_binary_features>0)
	{
		binary_feature_names.resize(num_binary_features);
		binary_weights.resize(num_binary_features,0);
		ind_active_binary_feature.resize(num_binary_features,true);

		int i;
		for (i=0; i<num_binary_features; i++)
		{
			int idx=-1;
			is >> idx >> binary_weights[i] >> binary_feature_names[i];
			if (idx<0 || binary_weights[i] == 0)
			{
				cout << "Error: reading binary weight line " << i << endl;
				exit(1);
			}
		}
	}

	is >> num_real_features;
	is >> num_non_zero_real_features;
	if (num_real_features>0)
	{
		real_feature_names.resize(num_real_features);
		real_default_weights.resize(num_real_features,0);
		real_weights.resize(num_real_features);
		real_limits.resize(num_real_features);
		ind_active_real_feature.resize(num_real_features,false);
		non_zero_real_idxs.clear();
	
		int i;
		for (i=0; i<num_non_zero_real_features; i++)
		{
			int num_weights=0;
			int num_limits=0;
			int j;

			char buff[1024];

			is.getline(buff,1024);
			is.getline(buff,1024);
			int feature_idx;
			double weight;
			if (sscanf(buff,"%d\t%lf",&feature_idx,&weight) != 2)
			{
				cout << "Error: bad line in model file:" << endl << buff << endl;
				exit(1);
			}

			real_default_weights[feature_idx]=weight;
			
			const int len = strlen(buff);
			int tab_count=0;
			for (j=0; j<len && tab_count<2; j++)
				if (buff[j]=='\t')
					tab_count++;

			if (tab_count<2)
			{
				cout << "Error: bad tab count!" << endl;
				exit(1);
			}

			real_feature_names[feature_idx] = string(buff+j);

			ind_active_real_feature[feature_idx] = true;
			non_zero_real_idxs.push_back(feature_idx);

			string &fn = real_feature_names[feature_idx];

			is >> num_limits;
			real_limits[feature_idx].resize(num_limits,NEG_INF);

			for (j=0; j<num_limits; j++)
				is >> real_limits[feature_idx][j];

			is >> num_weights;
			real_weights[feature_idx].resize(num_weights,0);
			for (j=0; j<num_weights; j++)
				is >> real_weights[feature_idx][j];
		}

		// this last part is done to move to the next line (so the next getline gets a new line)
		char buff[16];
		is.getline(buff,10);
		if (strlen(buff)>1)
		{
			cout << "Error: bad parsing of RankBoostModel!" << endl;
			exit(1);
		}
	}

	sort(non_zero_real_idxs.begin(),non_zero_real_idxs.end());

	this->total_default_weight=0;
	int i;
	for (i=0; i<num_real_features; i++)
		total_default_weight+=real_default_weights[i];

	ind_was_initialized=true;

	return true;
}


void RankBoostModel::write_rankboost_model(ostream& os, bool ind_compress)
{
	// compresses the weights and limits before writing
	if (ind_compress)
		compress_real_limits_and_weights();

	if (! ind_was_initialized)
	{
		cout << "Error: trying to write an unitialized RankBoost model!" << endl;
		exit(1);
	}

	int i;
	int num_non_zero_bin = 0;
	for (i=0; i<num_binary_features; i++)
		if (binary_weights[i] != 0)
			num_non_zero_bin++;

	os << num_non_zero_bin << endl;
	os << setprecision(8);

	for (i=0; i<num_binary_features; i++)
		if (binary_weights[i] != 0)
			os << i << "\t" << binary_weights[i] << "\t" << binary_feature_names[i] << endl;

	// assumes the limits and weights were already compressed, otherwise the model
	// will be very large
	int num_non_zero_real = 0;
	for (i=0; i<num_real_features; i++)
		if (real_weights[i].size()>0)
		{
			num_non_zero_real++;
		}

	os << num_real_features << endl;
	os << num_non_zero_real << endl;
	for (i=0; i<num_real_features; i++)
	{
		if ( real_weights[i].size() == 0 )
			continue;

		os << i << "\t" << real_default_weights[i] << "\t" << real_feature_names[i] << endl;
		os << real_limits[i].size();
		int j;
		for (j=0; j<real_limits[i].size(); j++)
			os << "\t" << real_limits[i][j];
		os << endl;
		os << real_weights[i].size();
		for (j=0; j<real_weights[i].size(); j++)
			os << "\t" << real_weights[i][j];
		os << endl;
	}
}


/***********************************************************************
// compresses the number of limits and weights to represent steps in the weight
// function without repeatition of the same weights in different limits
************************************************************************/
void RankBoostModel::compress_real_limits_and_weights()
{
	int i;

	for (i=0; i<num_real_features; i++)
	{
		if (real_limits[i].size() < 1)
		{
			real_weights[i].clear();
			real_limits[i].clear();
			continue;
		}
		else
		{
			int j;
			if (real_default_weights[i]==0)
			{
				for (j=0; j<real_weights[i].size(); j++)
					if (real_weights[i][j] != 0)
						break;
				if (j== real_weights[i].size())
				{
					real_weights[i].clear();
					real_limits[i].clear();
					continue;
				}
			}
		}

		vector<float> new_limits;
		vector<weight_t> new_weights;

		new_weights.push_back(real_weights[i][0]);
		new_limits.push_back(real_limits[i][0]);

		int j;
		for (j=1; j<real_limits[i].size(); j++)
		{
			if (real_weights[i][j] != real_weights[i][j+1])
			{
				new_weights.push_back(real_weights[i][j]);
				new_limits.push_back(real_limits[i][j]);
			}
		}
		new_weights.push_back(real_weights[i][j]); // last weight is for values larger than the last limit

		real_limits[i]=new_limits;
		real_weights[i]=new_weights;
	}

	remove_default_weights();
}

/*************************************************************************
Removes the default weight by subratcing them from everybody else
Also moves the min/max values towards 0 if they are all above or below it
**************************************************************************/
void RankBoostModel::remove_default_weights()
{
	int i;
	for (i=0; i<num_real_features; i++)
	{
		if (this->real_default_weights[i] != 0)
		{
			int j;
			for (j=0 ;j<real_weights[i].size(); j++)
				real_weights[i][j] -= real_default_weights[i];
			real_default_weights[i]=0;
		}
	}
}

/*************************************************************************
**************************************************************************/
double RankBoostModel::calc_training_rank_error(const RankBoostDataset& training_ds) const
{
	const int num_samples = training_ds.get_num_samples();
	const vector<RankBoostSample>& samples = training_ds.get_samples();

	vector<weight_t> rank_scores;

	rank_scores.resize(num_samples);

	int i;
	for (i=0; i<num_samples; i++)
		rank_scores[i]=calc_rank_score(samples[i]);

	const vector<SamplePairWeight>& phi_support = training_ds.get_phi_support();
	double train_error=0;
	for (i=0; i<phi_support.size(); i++)
	{
		if (rank_scores[phi_support[i].idx1]>=rank_scores[phi_support[i].idx2])
			train_error+=phi_support[i].weight;
	}

	return train_error/training_ds.get_total_phi_weight();
}


struct rank_score_pair {
	rank_score_pair() : rank(999), score(NEG_INF) {};
	rank_score_pair(int _r, float _s) : rank(_r), score(_s) {};

	bool operator< (const rank_score_pair& other) const
	{
		return (score>other.score);
	}
	int rank;
	float score;
};

/***************************************************************************
if tag3 != -1, it is used to filter out samples that don't have the same tag value.
Calculates the weighted ranking error for a given dataset.
Also calculates the error in peak ranking for all peptides:
Avg predicted rank for strongest peak and sd of prediction for strongest peaks,
Avg predicted rank fro 2nd strongest peak...
****************************************************************************/
double RankBoostModel::calc_prediction_error(const RankBoostDataset& rank_ds, 
								 vector<peak_rank_stat>& peak_stats,
								 int tag3_filter_val,
								 int *ptr_num_groups_tested) const
{
	const int num_samples = rank_ds.get_num_samples();
	const vector<RankBoostSample>& samples = rank_ds.get_samples();
	const int num_groups = rank_ds.get_num_groups();

	
	vector<weight_t> rank_scores;

	rank_scores.resize(num_samples,NEG_INF);

	if (ptr_num_groups_tested)
	{
		vector<bool> group_inds;
		group_inds.resize(num_groups,false);
		int i;
		for (i=0; i<num_samples; i++)
			if (tag3_filter_val<0 || samples[i].tag3 == tag3_filter_val)
			{
				rank_scores[i]=calc_rank_score(samples[i]);
				group_inds[samples[i].groupIndex]=true;
			}
		*ptr_num_groups_tested=0;
		for (i=0; i<group_inds.size(); i++)
			if (group_inds[i])
				(*ptr_num_groups_tested)++; 
	}
	else
	{
		int i;
		for (i=0; i<num_samples; i++)
			if (tag3_filter_val<0 || samples[i].tag3 == tag3_filter_val)
				rank_scores[i]=calc_rank_score(samples[i]);
	}

	int i;
	const vector<SamplePairWeight>& phi_support = rank_ds.get_phi_support();
	double error=0;
	double phi_weight=0;
	for (i=0; i<phi_support.size(); i++)
	{
		const int idx1 = phi_support[i].idx1;
		const int idx2 = phi_support[i].idx2;

		if (rank_scores[idx1]<=NEG_INF || rank_scores[idx2]<=NEG_INF)
			continue;
		
		if (rank_scores[idx1]>=rank_scores[idx2])
			error+=phi_support[i].weight;

		phi_weight += phi_support[i].weight;
	}
	double weighted_error = error/phi_weight;

//	cout << ">>>" << tag3_filter_val << "\t" << error << "\t" << phi_weight << endl;

	vector< vector<rank_score_pair> > peptide_peak_ranks;
	peptide_peak_ranks.clear();
	peptide_peak_ranks.resize(num_groups);

	vector< vector<int> > peak_predictions; // org rank / predcited rank
	peak_predictions.resize(100); // upto 100 strong peaks...
		for (i=0; i<100; i++)
		peak_predictions[i].resize(100,0);

	for (i=0; i<num_samples; i++)
	{
		if (rank_scores[i]<=NEG_INF)
			continue;

		const RankBoostSample& sam = samples[i];
		peptide_peak_ranks[sam.groupIndex].push_back(rank_score_pair(sam.rank_in_group,rank_scores[i]));
	}

	for (i=0; i<peptide_peak_ranks.size(); i++)
	{
		if (peptide_peak_ranks[i].size() == 0)
			continue;

		sort(peptide_peak_ranks[i].begin(),peptide_peak_ranks[i].end());

		int j;
		for (j=0; j<peptide_peak_ranks[i].size(); j++)
		{
			const int org_rank = peptide_peak_ranks[i][j].rank;
			const int predicted_rank = (j > 99 ? 99 : j);

			if (org_rank<100)
				peak_predictions[org_rank][predicted_rank]++;
		}
	}


	peak_stats.clear();

	for (i=0; i<100; i++)
	{
		vector<int> vals,counts;
		int peak_count=0;
		int j;
		for (j=0; j<peak_predictions[i].size(); j++)
			peak_count += peak_predictions[i][j];

		if (peak_count < 5)
			break;

		double mean=0,sd=-1; 
		for (j=0; j<peak_predictions[i].size(); j++)
		{
			if (peak_predictions[i][j]>0)
			{
				vals.push_back(j);
				counts.push_back(peak_predictions[i][j]);
			}
		}
		calc_mean_sd_from_counts(vals,counts,&mean,&sd);

		peak_rank_stat prs;

		prs.true_rank=i;
		prs.avg_precited_rank = mean;
		prs.sd_predicted_rank = sd;
		prs.precent_predicted_correctly = peak_predictions[i][i]/(float)peak_count;
		
		peak_stats.push_back(prs);
	}
	
	return weighted_error;
}

/*************************************************************************
**************************************************************************/
weight_t RankBoostModel::calc_rank_score(const RankBoostSample& sample) const
{
	const vector<int>& binary_idxs = sample.binary_non_zero_idxs;
	weight_t score=0;

	if (! ind_was_initialized)
	{
		cout << "Error: RankBoostModel not initialized!" << endl;
		exit(1);
	}

	int i;
	for (i=0; i<binary_idxs.size(); i++)
		score += binary_weights[binary_idxs[i]];

	// first add all no vote scores
	score += total_default_weight;

	for (i=0; i<sample.real_active_idxs.size(); i++)
	{
		const int feature_idx = sample.real_active_idxs[i];
		
		if (real_weights[feature_idx].size()==0) // inactive feature
			continue;

		const int bin_idx = get_real_bin_idx_for_value(feature_idx,sample.real_active_values[i]);
		score += real_weights[feature_idx][bin_idx];
		score -= real_default_weights[feature_idx]; // remove novote scores that shouldn't have been added
	}

	return score;
}


/*************************************************************************
**************************************************************************/
weight_t RankBoostModel::calc_rank_score_with_details(
						const RankBoostSample& sample,
						vector<string>& feature_names,
						vector<float>&	feature_values,
						vector<float>&  feature_scores) const
{
	const vector<int>& binary_idxs = sample.binary_non_zero_idxs;

	feature_names.clear();
	feature_scores.clear();
	weight_t score=0;

	if (! ind_was_initialized)
	{
		cout << "Error: RankBoostModel not initialized!" << endl;
		exit(1);
	}

	int i;
	for (i=0; i<binary_idxs.size(); i++)
		score += binary_weights[binary_idxs[i]];

	// first add all no vote scores
	score += total_default_weight;

	feature_names.push_back("defaults");
	feature_values.push_back(0);
	feature_scores.push_back(score);

	for (i=0; i<sample.real_active_idxs.size(); i++)
	{
		const int feature_idx = sample.real_active_idxs[i];
		
		if (real_weights[feature_idx].size()==0) // inactive feature
			continue;

		const int bin_idx = get_real_bin_idx_for_value(feature_idx,sample.real_active_values[i]);
		score += real_weights[feature_idx][bin_idx];
		score -= real_default_weights[feature_idx]; // remove novote scores that shouldn't have been added

		score_t net_score = real_weights[feature_idx][bin_idx] - real_default_weights[feature_idx];

		if (net_score != 0)
		{
			feature_names.push_back(real_feature_names[feature_idx]);
			feature_values.push_back(sample.real_active_values[i]);
			feature_scores.push_back(net_score);
		//	cout << feature_idx << "\t" << real_feature_names[feature_idx] << "\t" << net_score << endl;
		}
	//	cout << endl;
	}

	return score;
}


/*************************************************************************
**************************************************************************/
weight_t RankBoostModel::calc_rank_score_with_details(
						const RankBoostSample& sample,
						vector<int>  &	feature_idxs,
						vector<float>&	feature_values,
						vector<float>&  feature_scores) const
{
	const vector<int>& binary_idxs = sample.binary_non_zero_idxs;

	feature_idxs.clear();
	feature_values.clear();
	feature_scores.clear();
	weight_t score=0;

	if (! ind_was_initialized)
	{
		cout << "Error: RankBoostModel not initialized!" << endl;
		exit(1);
	}

	int i;
	for (i=0; i<binary_idxs.size(); i++)
		score += binary_weights[binary_idxs[i]];

	// first add all no vote scores
	score += total_default_weight;

	for (i=0; i<sample.real_active_idxs.size(); i++)
	{
		const int feature_idx = sample.real_active_idxs[i];
		
		if (real_weights[feature_idx].size()==0) // inactive feature
			continue;

		const int bin_idx = get_real_bin_idx_for_value(feature_idx,sample.real_active_values[i]);
		score += real_weights[feature_idx][bin_idx];
		score -= real_default_weights[feature_idx]; // remove novote scores that shouldn't have been added

		score_t net_score = real_weights[feature_idx][bin_idx] - real_default_weights[feature_idx];

		if (net_score != 0)
		{
			feature_idxs.push_back(feature_idx);
			feature_values.push_back(sample.real_active_values[i]);
			feature_scores.push_back(net_score);
		}
	}

	return score;
}



/*************************************************************************
Initializes a rank boost model 
**************************************************************************/
void RankBoostModel::init_rankboost_model_feature_names(	
							   const vector<string>& _binary_feature_names,
							   const vector<string>& _real_feature_names)
{
	binary_feature_names = _binary_feature_names;
	real_feature_names   = _real_feature_names;

	num_binary_features = binary_feature_names.size();
	num_real_features   = real_feature_names.size();

	real_default_weights.resize(num_real_features,0);
	real_weights.resize(num_real_features);
	real_limits.clear();
	
	binary_weights.clear();
	binary_weights.resize(num_binary_features,0);

	non_zero_binary_idxs.clear();
	non_zero_real_idxs.clear();

}


/***********************************************************************
Initializes the list of active features and limits for the real valued
features according to the data provided by the training_set.
Each active binary feature must have the given minimum number of active
samples per feature. Each bin in the real valued feature must also have
at least the same number of samples in it.
************************************************************************/
void RankBoostModel::init_rankboost_model_for_training(
						const RankBoostDataset& training_ds,
						int	 min_num_active_samples_for_feature,
						int	 max_num_real_bins_for_real_feature)
{

//	cout << "set limits..." << endl;
	set_real_limits(training_ds.get_samples(), 
					min_num_active_samples_for_feature, 
					max_num_real_bins_for_real_feature);

//	cout << "select active..." << endl;
	select_active_features(training_ds, min_num_active_samples_for_feature);
}




void RankBoostModel::summarize_features(const vector<RankBoostSample>& samples,
										ostream& os)
{
	const int num_samples = samples.size();
	vector<int> inactive_binary_count,  non_zero_binary_count;
	vector<int> inactive_real_count, non_zero_real_count;
	vector< vector<float> > real_vals;

	inactive_binary_count.resize(num_binary_features,0);
	non_zero_binary_count.resize(num_binary_features,0);

	inactive_real_count.resize(num_real_features,0);
	non_zero_real_count.resize(num_real_features,0);

	real_vals.resize(num_real_features);

	int i;
	for (i=0; i<num_samples; i++)
	{
		const RankBoostSample& sam = samples[i];

		int j;
		for (j=0; j<sam.binary_novote_idxs.size(); j++)
		{
			sam.print(os);
			inactive_binary_count[sam.binary_novote_idxs[j]]++;
		}

		for (j=0; j<sam.binary_non_zero_idxs.size(); j++)
			non_zero_binary_count[sam.binary_non_zero_idxs[j]]++;

		for (j=0; j<sam.real_novote_idxs.size(); j++)
		{
			inactive_real_count[sam.real_novote_idxs[j]]++;
			sam.print(os);
		}

		for (j=0; j<sam.real_active_idxs.size(); j++)
		{
			non_zero_real_count[sam.real_active_idxs[j]]++;
			real_vals[sam.real_active_idxs[j]].push_back(sam.real_active_values[j]);
		}
	}

	if (binary_feature_names.size()>0)
	{
		os << "REPORT ON BINARY FEATURES:" << endl;
		os << "--------------------------" << endl;
		for (i=0; i<binary_feature_names.size(); i++)
		{
			if (inactive_binary_count[i]<num_samples)
			{
				os << i << "\t" << 1.0 - inactive_binary_count[i]/(double)num_samples << "\t";
				os << non_zero_binary_count[i] << "\t" << non_zero_binary_count[i]/(double)(num_samples-inactive_binary_count[i]) << "\t";
				os << binary_feature_names[i] << endl;
			}
		}
		os << endl;
	}
	os << "REPORT ON REAL FEATURES:" << endl;
	os << "------------------------" << endl;
	for (i=0; i<real_feature_names.size(); i++)
	{
		if (inactive_real_count[i]<num_samples-10)
		{
			os << i << "\t" << 1.0 - inactive_real_count[i]/(double)num_samples << "\t";
			os << non_zero_real_count[i] << "\t" << non_zero_real_count[i]/(double)(num_samples-inactive_real_count[i]) << "\t";
			os << real_feature_names[i] << endl;
			if (real_vals[i].size()>20)
			{
				sort(real_vals[i].begin(),real_vals[i].end());
				int s_inc = real_vals[i].size()/10;
				int j;
				for (j=0; j<10; j++)
					os << real_vals[i][j*s_inc] << " ";
				os << endl;
			}
			os << endl;

	/*		if (i == 1 )
			{
				cout << "FEATURE " << i << endl;
				int j;
				for (j=0; j<real_vals[i].size(); j++)
					cout << real_vals[i][j] << "\t";
				cout << endl;
			}*/
		}
	}
	cout << endl;
}


/******************************************************************************
*******************************************************************************/
void RankBoostModel::summarize_features_pos_neg(
									const vector<RankBoostSample>& pos_samples, 
									const vector<RankBoostSample>& neg_samples)
{
	const int num_pos_samples = pos_samples.size();
	const int num_neg_samples = neg_samples.size();

	vector<int> pos_inactive_real_count, pos_non_zero_real_count;
	vector<int> neg_inactive_real_count, neg_non_zero_real_count;
	vector< vector<float> > pos_real_vals, neg_real_vals;


	pos_inactive_real_count.resize(num_real_features,0);
	pos_non_zero_real_count.resize(num_real_features,0);
	neg_inactive_real_count.resize(num_real_features,0);
	neg_non_zero_real_count.resize(num_real_features,0);

	pos_real_vals.resize(num_real_features);
	neg_real_vals.resize(num_real_features);

	int i,j;
	for (i=0; i<num_pos_samples; i++)
	{
		const RankBoostSample& pos_sam = pos_samples[i];

		for (j=0; j<pos_sam.real_active_idxs.size(); j++)
		{
			pos_non_zero_real_count[pos_sam.real_active_idxs[j]]++;
			pos_real_vals[pos_sam.real_active_idxs[j]].push_back(pos_sam.real_active_values[j]);
		}
	}

	for (i=0; i<num_neg_samples; i++)
	{
		const RankBoostSample& neg_sam = neg_samples[i];

		for (j=0; j<neg_sam.real_active_idxs.size(); j++)
		{
			neg_non_zero_real_count[neg_sam.real_active_idxs[j]]++;
			neg_real_vals[neg_sam.real_active_idxs[j]].push_back(neg_sam.real_active_values[j]);
		}
	}


	cout << "REPORT ON REAL FEATURES:" << endl;
	cout << "------------------------" << endl;
	for (i=0; i<real_feature_names.size(); i++)
	{
		if (pos_inactive_real_count[i]>=num_pos_samples-10 &&
			neg_inactive_real_count[i]>=num_neg_samples-10)
			continue;

		if (pos_inactive_real_count[i]<num_pos_samples-10)
		{
			cout << i <<  "\t";
			cout << pos_non_zero_real_count[i] << "\t" << pos_non_zero_real_count[i]/(double)(num_pos_samples-pos_inactive_real_count[i]) << "\t";
			cout << neg_non_zero_real_count[i] << "\t" << neg_non_zero_real_count[i]/(double)(num_neg_samples-neg_inactive_real_count[i]) << "\t";
			cout << real_feature_names[i] << endl;
			if (pos_real_vals[i].size()>20)
			{
				cout << "POS: ";
				sort(pos_real_vals[i].begin(),pos_real_vals[i].end());
				int s_inc = pos_real_vals[i].size()/10;
				int j;
				for (j=0; j<10; j++)
					cout << pos_real_vals[i][j*s_inc] << " ";
				cout << endl;
			}
		}
		else
			cout << i<< " not enough pos" << endl;

		if (neg_inactive_real_count[i]<num_neg_samples-10)
		{
			if (neg_real_vals[i].size()>20)
			{
				cout << "NEG: ";
				sort(neg_real_vals[i].begin(),neg_real_vals[i].end());
				int s_inc = neg_real_vals[i].size()/10;
				int j;
				for (j=0; j<10; j++)
					cout << neg_real_vals[i][j*s_inc] << " ";
				cout << endl;
			}
		}
		else
			cout << "Not enough neg" << endl;
		cout << endl;
	}
	cout << endl;
}



/******************************************************************
Chooses what values to use as limits for the real valued features.
If there is only one dominant value, or the number of samples for
which the feature is active is too small, then the feature gets
deactivated.
*******************************************************************/
void RankBoostModel::set_real_limits(const vector<RankBoostSample>& samples,
									 int min_num_samples_for_feature,
									 int max_num_bins,
									 bool clear_weights)
{
	const int num_samples = samples.size();
	vector< vector<float> > real_vals;

	real_vals.clear();
	real_limits.clear();

	real_vals.resize(num_real_features);
	real_limits.resize(num_real_features);

	int i;
	for (i=0; i<num_real_features; i++)
		real_vals[i].clear();

	for (i=0; i<num_samples; i++)
	{
		const RankBoostSample& sam = samples[i];

		int j;
	
		for (j=0; j<sam.real_active_idxs.size(); j++)
		{
			const int f_idx = sam.real_active_idxs[j];
			real_vals[f_idx].push_back(sam.real_active_values[j]);
		}
	}


	for (i=0; i<num_real_features; i++)
	{
		real_limits[i].clear();
		if (real_vals.size() <min_num_samples_for_feature )
			continue;

		sort(real_vals[i].begin(),real_vals[i].end());

		// find number of unique values. If it is below max_num_bins, then 
		// make sure there are at least min_bin_size samples per bin. If there are
		// more than that number, split into equal sized bins

		vector<float> unique_vals;
		vector<int>   counts;
		const int num_real_vals = real_vals[i].size();

		int j;
		for (j=0; j<num_real_vals; j++)
		{
			const float val = real_vals[i][j];
			int k;
			for (k=0; k<unique_vals.size(); k++)
				if (val == unique_vals[k])
				{
					counts[k]++;
					break;
				}

			if (k<unique_vals.size())
				continue;
			
			unique_vals.push_back(val);
			counts.push_back(1);

			if (unique_vals.size()>max_num_bins)
				break;
		}

		// end bin size and mid bin size are the sizes that can be given to bins
		// of real values. The end bin values are larger to reduce incedents of overfitting
		// the end values.
		int bin_size = (min_num_samples_for_feature/2)+1;
		if (num_real_vals/max_num_bins>bin_size)
			bin_size = num_real_vals/max_num_bins;

		real_limits[i].push_back(NEG_INF);

		// go according to unique vales
		if (j==num_real_vals)
		{	
			int k;
			for (k=0; k<unique_vals.size(); k++)
			{
				if (k<unique_vals.size()-1 &&
					counts[k]<bin_size)
				{
					counts[k+1]+=counts[k];
					continue;
				}

				if (counts[k]>=bin_size)
					real_limits[i].push_back(unique_vals[k]);
			}

			if (real_limits[i].size()>1)
				real_limits[i].pop_back(); // remove last limit so it includes all values larger than it too
		}
		else
		{
			if (num_real_vals > min_num_samples_for_feature)
			{
				const int max_idx = num_real_vals - bin_size;
				int idx = bin_size;

				float last_limit = NEG_INF;
				while (idx<max_idx)
				{
					float new_limit = real_vals[i][idx];
					if (new_limit != last_limit)
					{
						real_limits[i].push_back(new_limit);
						last_limit = new_limit;
					}
					idx+=bin_size;
				}
				if (real_limits[i].size()>1)
					real_limits[i].pop_back(); // remove last limit so it includes all values larger than it too
			}
		}
	}


	if (clear_weights)
	{
		for (i=0; i<num_real_features; i++)
		{
			real_weights[i].clear();
			real_weights[i].resize(real_limits[i].size()+1,0);
		}
	}

	cout << "Set real limits for " << num_real_features << " real features..." << endl;
}


/************************************************************************
Chooses what features should be considered active (i.e., have their weights
updated in the boosting rounds). Binary features must have at least 
min_sample_count active samples for each feature
*************************************************************************/
void RankBoostModel::select_active_features(const RankBoostDataset& training_ds,
											int min_sample_count)
{

	ind_active_binary_feature.resize(num_binary_features,true);
	ind_active_real_feature.resize(num_real_features,true);

	const vector<RankBoostSample>& samples = training_ds.get_samples();
	const int num_samples                  = samples.size();

	vector<int> binary_one_counts,  binary_zero_counts, binary_inactive_counts;
	
	if (real_limits.size() != real_feature_names.size())
	{
		cout << "Error: must first set real feature limits!" << endl;
		exit(1);
	}

	binary_one_counts.resize(num_binary_features,0);
	binary_zero_counts.resize(num_binary_features,0);
	binary_inactive_counts.resize(num_binary_features,0);

	


	int i;
	for (i=0; i<num_samples; i++)
	{
		const RankBoostSample& sam = samples[i];
		vector<bool> used_inds;

		used_inds.resize(num_binary_features,false);

		int j;
		for (j=0; j< sam.binary_novote_idxs.size(); j++)
		{
			binary_inactive_counts[sam.binary_novote_idxs[j]]++;
			used_inds[sam.binary_novote_idxs[j]]=true;
		}
		
		for (j=0; j<sam.binary_non_zero_idxs.size(); j++)
		{
			binary_one_counts[sam.binary_non_zero_idxs[j]]++;
			used_inds[sam.binary_non_zero_idxs[j]]=true;
		}

		for (j=0; j<num_binary_features; j++)
			if (! used_inds[j])
				binary_zero_counts[j]++;
	}

	// inactiveate binary features
	for (i=0; i<num_binary_features; i++)
	{
		if (binary_one_counts[i]<min_sample_count ||
			binary_zero_counts[i]<min_sample_count)
		{
			ind_active_binary_feature[i]=false;
			cout << "Inactivating BINARY feature " << i << " " << this->binary_feature_names[i];
			cout << "  (counts 0: " << binary_zero_counts[i] << " 1: " <<
				binary_one_counts[i] << " nv: " << binary_inactive_counts[i] << ")" << endl;

		}
	}

	// inactivate real features with no ranges
	for (i=0; i<num_real_features; i++)
	{
		if (real_limits[i].size()==0)
		{
			ind_active_real_feature[i]=false;
			cout << "Inactivating REAL feature " << i << " " << real_feature_names[i] << endl;
		}
	}

}



/************************************************************************
Performs the weak learn function on the binary features.
Assumes that q_def is always 0 (so its selection is ignored).
The function assumes that the potentials have already been computed.
The function examines all features and returns the one for which |r| is
the largest.
*************************************************************************/
void RankBoostModel::binary_weak_learn(const vector<weight_t>& potentials,
						   const RankBoostDataset& rank_ds,
						   int&  best_feature_idx,
						   double& best_r_value,
						   bool only_non_zero_features) const
{
	best_feature_idx=-1;
	best_r_value = 0;
	
	if (only_non_zero_features && non_zero_binary_idxs.size() == 0)
		return;

	int nz_idx=0;
	int i;
	for (i=0; i<num_binary_features; i++)
	{
		if (only_non_zero_features)
		{
			if (nz_idx==non_zero_binary_idxs.size())
				return;

			if (i<non_zero_binary_idxs[nz_idx])
			{
				i=non_zero_binary_idxs[nz_idx];
				nz_idx++;
			}
		}

		if (! ind_active_binary_feature[i])
			continue;

		const vector<int>& active_idxs = rank_ds.get_binary_one_list(i);
		double r=0;
		int j;

		for (j=0; j<active_idxs.size(); j++)
			r+=potentials[active_idxs[j]];

		if (fabs(r) > fabs(best_r_value))
		{
			best_feature_idx=i;
			best_r_value = r;
		}
	}
}


/************************************************************************
Performs the weak learn function on the real functions.
he function assumes that the potentials have already been computed.
The function finds the feature and theta threshold for which |r| is maximized.
The function also returns the q_def (0 or 1) which should be given to samples
that did not vote on this feature (if q_def is 1 than the default weights 
for this feature should be updated accordingly).
*************************************************************************/
void RankBoostModel::real_weak_learn(const vector<weight_t>& potentials,
						 const RankBoostDataset& rank_ds,
						 int&  best_feature_idx,
						 int&  theta_idx,
						 int&  q_def,
						 double& best_r_value,
						 bool	only_non_zero_features) const
{
	const int num_samples = rank_ds.get_num_samples();

	best_feature_idx=0;
	theta_idx=-1;
	q_def=-1;
	best_r_value=0;

	int nz_idx =0;

	int i;
	for (i=0; i<num_real_features; i++)
	{
		if (only_non_zero_features)
		{
			if (nz_idx==non_zero_real_idxs.size())
				return;
			
			i=non_zero_real_idxs[nz_idx];
			nz_idx++;
		}

		if (! ind_active_real_feature[i])
			continue;

		const int num_thetas = real_limits[i].size()+1;
		double R=0;
		vector<double> L_bin_values;

		L_bin_values.resize(num_thetas,0);

		const vector< vector<int> >& lists = rank_ds.get_real_vote_lists(i);
		int j;
		for (j=0; j<lists.size(); j++)
		{
			const vector<int>& bin_list = lists[j];
			int k;
			for (k=0; k<bin_list.size(); k++)
				L_bin_values[j]+=potentials[bin_list[k]];

			R+=L_bin_values[j];
		}
		
		// the values in bin j are all greater than the limit of index j-1
		// therefore we return theta index j, and update all weights in positions
		// j and up. The first bin (j=0) is always empty (it has all values lower than
		// -infinity (empty), so it is not considered.

		double L=0;
		for (j=L_bin_values.size()-1; j>0; j--)
		{
			int q=0; // default for novote

			L+=L_bin_values[j];

			if (fabs(L)<=fabs(L-R))
				q=1;

			if (fabs(L-(double)q*R)>fabs(best_r_value))
			{
				best_feature_idx = i;
				best_r_value = L-(double)q*R;
				theta_idx = j;
				q_def = q;

		/*		if (fabs(best_r_value)>1)
				{
					cout << "Feature: " << i << "," << j << "  best_r_value: " << best_r_value << endl;
					int k;
					for (k=0; k<real_limits[i].size(); k++)
						cout << real_limits[i][k] << "\t";
					cout << endl;
					for (k=0; k<L_bin_values.size(); k++)
						cout << L_bin_values[k] << "\t";
					cout << endl;


					cout << endl << endl;
				}*/
			}
		}
	}
}


/************************************************************************
Performs the weak learn function on the real functions.
he function assumes that the potentials have already been computed.
The function finds the feature and theta threshold for which |r| is maximized.
The function also returns the q_def (0 or 1) which should be given to samples
that did not vote on this feature (if q_def is 1 than the default weights 
for this feature should be updated accordingly).
*************************************************************************/
void RankBoostModel::real_weak_learn_double_theta(const vector<weight_t>& potentials,
						 const RankBoostDataset& rank_ds,
						 int&  best_feature_idx,
						 int&  theta_start_idx,
						 int&  theta_end_idx,
						 int&  q_def,
						 double& best_r_value,
						 bool	only_non_zero_features) const
{
	const int num_samples = rank_ds.get_num_samples();

	best_feature_idx=0;
	theta_start_idx=-1;
	theta_end_idx =-1;
	q_def=-1;
	best_r_value=0;

	int nz_idx =0;

	int i;
	for (i=0; i<num_real_features; i++)
	{
		if (only_non_zero_features)
		{
			if (nz_idx==non_zero_real_idxs.size())
				return;
			
			i=non_zero_real_idxs[nz_idx];
			nz_idx++;
		}

		if (! ind_active_real_feature[i])
			continue;

		const int num_thetas = real_limits[i].size()+1;
		const int last_bin_idx = num_thetas - 1;
		double R=0;
		vector<double> L_bin_values;

		L_bin_values.resize(num_thetas,0);

		const vector< vector<int> >& lists = rank_ds.get_real_vote_lists(i);
		int j;
		for (j=0; j<lists.size(); j++)
		{
			const vector<int>& bin_list = lists[j];
			int k;
			for (k=0; k<bin_list.size(); k++)
				L_bin_values[j]+=potentials[bin_list[k]];

			R+=L_bin_values[j];
		}
		
		// the values in bin j are all greater than the limit of index j-1
		// therefore we return theta index j, and update all weights in positions
		// j and up. The first bin (j=0) is always empty (it has all values lower than
		// -infinity (empty), so it is not considered.

		double L=0;
		for (j=last_bin_idx; j>0; j--)
		{
			int q=0; // default for novote

			L+=L_bin_values[j];

			if (fabs(L)<=fabs(L-R))
				q=1;

			if (fabs(L-(double)q*R)>fabs(best_r_value))
			{
				best_feature_idx = i;
				best_r_value = L-(double)q*R;
				theta_start_idx = j;
				theta_end_idx   = last_bin_idx;
				q_def = q;
			}

			double L_prime = L;
			int k;
			for (k=last_bin_idx; k>j; k--)
			{
				L_prime -= L_bin_values[k];
				if (fabs(L_prime)<=fabs(L_prime-R))
					q=1;

				if (fabs(L_prime-(double)q*R)>fabs(best_r_value))
				{
					best_feature_idx = i;
					best_r_value = L_prime-(double)q*R;
					theta_start_idx = j;
					theta_end_idx   = k-1;
					q_def = q;
				}
			}
		}
	}
}


/************************************************************************
Performs the training of the rank boost model.
This function assumes that the training_ds was already intialized (has
phi_support computed already), and that the model has been initialized
for this dataset (feature names set, active features selected).
The algorithm performs a designated number of rounds (where each round is
consdiered an update of all active featrues).
If a test_set is given, reports about the test error are given for each round.
*************************************************************************/
bool RankBoostModel::train_rankboost_model(
							   const RankBoostDataset& training_ds,
							   int   max_num_rounds,
							   vector<idx_weight_pair>* top_misclassified_pairs,
							   RankBoostDataset* test_ds,
							   int   test_tag3_filter_val,
							   char *report_prefix,
							   char *stop_signal_file,
							   const vector<string>* model_header_strings) // these appear at the top of the boost model file
																		   // in case there is extra info not in the ordinary
																		   // boost models.
{
	fstream out_file_stream, stat_file_stream, train_res, test_res, flist_stream;
	
	bool use_cout = true;
	int running_feature_report_idx = -1; // if set to a feature idx, will give a report on the feature's progress

	// open streams for reporting
	if (report_prefix)
	{
		char buff[512];
		sprintf(buff,"%s_progress.txt",report_prefix);
		out_file_stream.open(buff,ios::out);
		if (! out_file_stream.is_open() || ! out_file_stream.good())
		{
			cout << "Error: coudln't open out_file for wrtiting: " << buff << endl;
			exit(1);
		}

		sprintf(buff,"%s_stats.txt",report_prefix);
		stat_file_stream.open(buff,ios::out);
		if (! stat_file_stream.is_open() || ! stat_file_stream.good())
		{
			cout << "Error: coudln't open test stats for wrtiting: " << buff << endl;
			exit(1);
		}
	
		sprintf(buff,"%s_train_res.txt",report_prefix);
		train_res.open(buff,ios::out);
		if (! train_res.is_open() || ! train_res.good())
		{
			cout << "Error: coudln't open train res for wrtiing: " << buff << endl;
			exit(1);
		}

		if (test_ds)
		{
			sprintf(buff,"%s_test_res.txt",report_prefix);
			test_res.open(buff,ios::out);
			if (! test_res.is_open() || ! test_res.good())
			{
				cout << "Error: coudln't open test res for writing: " << buff << endl;
				exit(1);
			}
		}
		use_cout = false;
	}

	ostream&  out_stream  = (use_cout ? cout : out_file_stream);
	ostream&  stat_stream = (use_cout ? cout : stat_file_stream);

	vector<weight_t> D,D0, max_D_for_normal_updates;
	vector<weight_t> *p_max_D_for_normal_updates = NULL;

	
	training_ds.compute_initial_distribution(D0);

	if (training_ds.get_max_ratio_for_regular_update()<9999)
	{
		training_ds.compute_max_weights_for_normal_update(D0,max_D_for_normal_updates);
		p_max_D_for_normal_updates = &max_D_for_normal_updates;
	}

	real_first_updates.resize(num_real_features,NEG_INF);
	real_update_counts.resize(num_real_features,0);
	binary_update_counts.resize(num_binary_features,0);
	D=D0;

	double Z_prod = 1.0;

	total_default_weight =0;

	best_test_error  = 1.0;
	best_train_error = 1.0;
	
	this->ind_was_initialized = true;

	time_t start_time = time(NULL);

	cout << "Running boosting for at most " << max_num_rounds << " iterations..." << endl;
	int t;
	for (t=1; t<=max_num_rounds; t++)
	{
		current_round = t;

		// determine how often we report progress and test training/test error
		int report_freq = 1; 
		if (t>10)     report_freq=10;
		if (t>500)    report_freq=50;
		if (t>1000)   report_freq=100;
		if (t>5000)   report_freq=500;
		if (t>10000)  report_freq=1000;
		if (t>100000) report_freq=5000;
		if (t>500000) report_freq=10000;

		const int feature_report_rounds = 10000;

		bool report_this_round = ((t % report_freq) == 0);
		
		bool use_double_theta = true;

		bool ind_only_non_zero_features = false;
		if (t>100    && t%5)   ind_only_non_zero_features = true;
		if (t>1000   && t%10)  ind_only_non_zero_features = true;
		if (t>10000  && t%25)  ind_only_non_zero_features = true;
		if (t>100000 && t%100) ind_only_non_zero_features = true;

		
		int best_binary_idx=-1;
		double best_binary_r = 0;

		int real_feature_idx=-1;
		int real_theta_start_idx = -1;
		int real_theta_end_idx   = -1;
		int real_q_def = -1;
		double real_r = 0;

		vector<weight_t> potentials;

		training_ds.calc_potentials(D,potentials);

		if (num_binary_features>0)
			binary_weak_learn(potentials,training_ds,best_binary_idx,best_binary_r);

		if (num_real_features>0)
		{
			if (use_double_theta)
			{
				real_weak_learn_double_theta(potentials, training_ds, real_feature_idx,
					real_theta_start_idx, real_theta_end_idx,
					real_q_def, real_r, ind_only_non_zero_features );
			}
			else
				real_weak_learn(potentials, training_ds, real_feature_idx,
					real_theta_start_idx,  real_q_def, real_r, ind_only_non_zero_features );
		}

		if (real_r == 0 && best_binary_r == 0)
		{
			if (t<2)
			{
				cout << "Error: model converged to quickly, there is a problem with the feature values!" << endl;
				exit(1);
			}
			break;
		}
		
		double Z=1.0;
		if (fabs(real_r)>fabs(best_binary_r))
		{
			const weight_t alpha = 0.5 * log((1+real_r)/(1-real_r));
			const float theta_start = real_limits[real_feature_idx][real_theta_start_idx-1];
			const float theta_end = ( (real_theta_end_idx>0 &&
									  real_theta_end_idx<real_limits[real_feature_idx].size()) 
									  ? real_limits[real_feature_idx][real_theta_end_idx] :
									  POS_INF);
			
			Z=training_ds.update_distribution_according_to_real_feature(real_feature_idx, 
				theta_start, theta_end, real_q_def,alpha, D , p_max_D_for_normal_updates, false);

			update_model_weights_for_real_feature(alpha, real_feature_idx, real_q_def,
				real_theta_start_idx, real_theta_end_idx);

			if (real_first_updates[real_feature_idx]<0)
				real_first_updates[real_feature_idx]=t;
		}	
		else
		{
			const weight_t alpha = 0.5 * log((1+best_binary_r)/(1-best_binary_r));

			Z=training_ds.update_dsitribution_according_to_binary_feature(best_binary_idx, alpha, D);
			
			update_model_weights_for_binary_feature(best_binary_idx, alpha);
		}
		Z_prod *= Z;

		// This part is for the feature reports / not really part of the trianing
		if (running_feature_report_idx>=0 && 
			real_feature_idx == running_feature_report_idx)
		{
			const int num_updates = real_update_counts[running_feature_report_idx];
			int sqr=1;
			while (sqr<num_updates)
				sqr*=2;

			if (sqr == num_updates ) // output only for 1,2,4,8,... updates
			{
				if (! flist_stream.is_open())
				{
					char buff[512];
					sprintf(buff,"%s_fprog_%d.txt",report_prefix,running_feature_report_idx);
					flist_stream.open(buff,ios::out);
					if (! flist_stream.is_open() || ! flist_stream.good())
					{
						cout << "Error: coudln't open flist for writing: " << buff << endl;
						exit(1);
					}
				}
				ouput_importance_ranked_feature_list(training_ds,flist_stream,running_feature_report_idx,t);
				cout << t << "\t" << "feature report: " << num_updates << endl;
			}
		}


		// all output and test should be done in this area
		if (report_this_round)
		{
			time_t current_time = time(NULL);
			double total_time = current_time - start_time;
			clock_t start_test = clock();

			int num_tested_peptides_in_train=0;
			int num_tested_peptides_in_test=0;

			int *ptr_train_num = (t<=1 ? &num_tested_peptides_in_train : NULL);
			int *ptr_test_num  = (t<=1 ? &num_tested_peptides_in_test : NULL);

			vector<peak_rank_stat> train_peak_stats, test_peak_stats;

			train_error = calc_prediction_error(training_ds, train_peak_stats, 
								test_tag3_filter_val, ptr_train_num);
			test_error  = 1.0;
			if (test_ds)
				test_error = calc_prediction_error(*test_ds, test_peak_stats, 
									test_tag3_filter_val, ptr_test_num);

			cout << "Round\t" << t << "\t#rf " << non_zero_real_idxs.size() << "\ttime: " << total_time << " secs.\t" << setprecision(5) << "[ " <<
				setprecision(5) << train_error << " " << test_error << "]" << endl;

			// these should only be outputted for the first round
			if (ptr_train_num)
			{
				out_stream << "ERRORS MEASURED FOR TAG3 VAL " << test_tag3_filter_val << endl;
				out_stream << "TRAIN ERRORS FROM " << *ptr_train_num << endl;
				if (report_prefix)
				{
					stat_stream << "ERRORS MEASURED FOR TAG3 VAL " << test_tag3_filter_val << endl;
					stat_stream << "TRAIN ERRORS FROM " << *ptr_train_num << endl;
				}
			}

			if (ptr_test_num)
			{
				out_stream << "TEST ERRORS FROM " << *ptr_test_num << endl;
				if (report_prefix)
					stat_stream << "TEST ERRORS FROM " << *ptr_test_num << endl;
			}
		
			out_stream << setprecision(7);
			out_stream << "Round " << t << "\t" << setprecision(7) << fixed  << (int)total_time << "\t";

			
			if (num_binary_features>0)
				out_stream << "Act bin " << non_zero_binary_idxs.size() << "/" << num_binary_features ;
			if (num_real_features>0)
				out_stream << " Act real " << non_zero_real_idxs.size() << "/" << num_real_features << endl;

			if (num_binary_features>0)
			{
				out_stream << "Best BINARY feature: " << best_binary_idx << " " << binary_feature_names[best_binary_idx] <<
					"    r: " << best_binary_r << endl;
			}

			if (num_real_features>0)
			{
				out_stream << "Best REAL feature  : " << real_feature_idx << " " << real_feature_names[real_feature_idx] << 
				    "   theta: " << real_theta_start_idx << "-" <<
					real_theta_end_idx << "  r: " <<  real_r << endl;
			}

			clock_t end_test = clock();
			double test_time = (end_test-start_test)/(double)CLOCKS_PER_SEC;

			out_stream << setprecision(6);
			out_stream << "train: " << train_error;
			if (test_ds)
				out_stream << "\ttest: " << test_error;
			out_stream << "\tZ_prod = " << Z_prod << "\t" << "(" << 
					test_time << ")" << endl;

			out_stream << endl;

			// full stats
			if (report_prefix)
			{
			
				stat_stream << fixed << setprecision(6);
				stat_stream << t << "\t" << (int)total_time << "\t" << non_zero_binary_idxs.size() <<
					"\t" << non_zero_real_idxs.size() << "\t" << train_error << "\t";

				int stat_size = train_peak_stats.size();
				if (test_ds)
				{
					stat_stream << test_error << "\t";
					if (test_peak_stats.size()<train_peak_stats.size())
						stat_size = test_peak_stats.size();
				}
				stat_stream << endl;

				// also to cout
				out_stream << setprecision(4) << t << "\t" << fixed << (int)total_time << "\t" << 
					non_zero_binary_idxs.size() << "\t" << non_zero_real_idxs.size() << 
					"\t" << train_error << "\t";

				if (test_ds)
					out_stream << test_error << "\t";
				out_stream << endl << endl;

				// peak stats
				train_res << t <<"\t" << setprecision(6) << train_error << "\t" << stat_size;
				int j;
				for (j=0; j<stat_size; j++)
					train_res << "\t" << setprecision(4) << train_peak_stats[j].precent_predicted_correctly;

				for (j=0; j<stat_size; j++)
					train_res << "\t" << setprecision(4) << train_peak_stats[j].avg_precited_rank + 1.0 << "\t"
							 <<	setprecision(5) << train_peak_stats[j].sd_predicted_rank;
				train_res << endl;

				if (test_ds)
				{
					test_res << t <<"\t" << setprecision(6) << test_error << "\t" << stat_size;
					int j;
					for (j=0; j<stat_size; j++)
						test_res << "\t" << setprecision(4) << test_peak_stats[j].precent_predicted_correctly;

					for (j=0; j<stat_size; j++)
						test_res << "\t" << setprecision(4) << test_peak_stats[j].avg_precited_rank + 1.0 << "\t"
								 <<	setprecision(5) << test_peak_stats[j].sd_predicted_rank;
					test_res << endl;
				}

				if (t % 100 == 0)
				{
					stat_stream.flush();
					out_stream.flush();	
				}

				if (0 && top_misclassified_pairs)
				{
					get_top_misclassified_pairs(training_ds,D,D0,*top_misclassified_pairs);
					print_top_misclassified_pairs(training_ds,D,D0,10,out_stream);
				}
			}

			// check if this is the best round
			// output the current model file
			if (test_error<best_test_error && ((t % 100 == 0) || (t<100 && t % 10 == 0)))
			{
				set_best_model_parameters_to_current_parameters();

				// write the model. Since compressing messes it up, we will copy the model
				// to a new one and write that one
				if (t % 100 == 0 && report_prefix)
				{
					RankBoostModel copy_of_model = *this;

					char name_buff[512];
					sprintf(name_buff,"%s_model.txt",report_prefix);
					ofstream model_stream(name_buff);
					if (! model_stream.is_open() || ! model_stream.good())
					{
						cout << "Error: couldn't feature_stream file for writing:" << name_buff << endl;
						exit(1);
					}

					if (model_header_strings)
					{
						int i;
						for (i=0; i<model_header_strings->size(); i++)
							model_stream << model_header_strings->at(i);
					}

					copy_of_model.write_rankboost_model(model_stream);
					model_stream.close();

				
					sprintf(name_buff,"%s_feature_list.txt",report_prefix);
					ofstream feature_stream(name_buff);
					if (! feature_stream.is_open() || ! feature_stream.good())
					{
						cout << "Error: couldn't feature_stream file for writing:" << name_buff << endl;
						exit(1);
					}
					copy_of_model.ouput_importance_ranked_feature_list(training_ds,feature_stream);
					feature_stream.close();
				}
			}
			
			// if for the last 20% of the rounds we had no improvement in the test error
			// we can stop!
			if (t>100)
			{
				double stop_ratio = (t>100000 ? 1.3 : 2.0);
				if (t<1000)
					stop_ratio = 3.0;
				if(t<20000)
					stop_ratio = 2.5;

				if (max_num_rounds<1000000 &&  // if it is 1000000 then let it run...
					t/(double)best_round_idx>=stop_ratio && 
					best_test_error<test_error)
				{
					out_stream << "TERMINATING AT ROUND " << t << ", NO PROGRESS IN TEST ERROR SINCE ROUND " << 
						best_round_idx << endl << endl;
					out_stream << fixed << setprecision(6) << "Current test error: " << test_error << ", best test error " << best_test_error << endl;
					break;
				}
			}

			
			// check signal
			if (stop_signal_file)
			{
				ifstream signal_stream(stop_signal_file);
				if (signal_stream.is_open())
				{
					out_stream << endl << "TERMINATED BECAUSE STOP FILE WAS DETECTED!" << endl;
					out_stream << "( " << stop_signal_file << " )" << endl;
					break;
				}
			}
		}
	}

	set_current_model_parameters_to_best_parameters();

	vector<peak_rank_stat> dummy_stats;
	train_error = calc_prediction_error(training_ds, dummy_stats, test_tag3_filter_val);
	test_error  = 1.0;
	if (test_ds)
		test_error = calc_prediction_error(*test_ds, dummy_stats, test_tag3_filter_val);

	out_stream << setprecision(6) << fixed << "FINAL ERRORS:" << endl;
	out_stream << "train:\t" << train_error << endl;
	out_stream << "test:\t" << test_error << endl;

	if (top_misclassified_pairs)
		get_top_misclassified_pairs(training_ds,D,D0,*top_misclassified_pairs);

	if (out_file_stream.is_open())
		out_file_stream.close();

	if (stat_file_stream.is_open())
		stat_file_stream.close();

	if (train_res.is_open())
		train_res.close();

	if (test_res.is_open())
		test_res.close();

	if (flist_stream.is_open())
		flist_stream.close();

	ind_was_initialized = true;

	cout << "Finished training boost model (" << t-1 << " rounds)" << endl;

	return (t== max_num_rounds); // normal termination
}



void RankBoostModel::set_best_model_parameters_to_current_parameters()
{
	best_round_idx   = current_round;
	best_train_error = train_error;
	best_test_error  = test_error;
	best_total_default_weight = total_default_weight;

	best_ind_active_binary_feature = ind_active_binary_feature;
	best_binary_weights = binary_weights;

	best_ind_active_real_feature = ind_active_real_feature;
	best_real_weights = real_weights; 
	best_real_limits = real_limits;
	best_real_default_weights = real_default_weights; 
	best_real_update_counts = real_update_counts; 
	best_binary_update_counts = binary_update_counts;

	best_non_zero_binary_idxs = non_zero_binary_idxs;
	best_non_zero_real_idxs = non_zero_real_idxs;
}

void RankBoostModel::set_current_model_parameters_to_best_parameters()
{
	current_round = best_round_idx;
	train_error = best_train_error;
	test_error = best_test_error;
	total_default_weight = best_total_default_weight;

	ind_active_binary_feature = best_ind_active_binary_feature;
	binary_weights = best_binary_weights;

	ind_active_real_feature = best_ind_active_real_feature;
	real_weights = best_real_weights; 
	real_limits = best_real_limits;
	real_default_weights = best_real_default_weights; 
	real_update_counts = best_real_update_counts; 
	binary_update_counts = best_binary_update_counts;

	non_zero_binary_idxs = best_non_zero_binary_idxs;
	non_zero_real_idxs = best_non_zero_real_idxs;
}


/********************************************************************************
Changes the current weights in the model.
Since it is a binary variable, only need to add weight
*********************************************************************************/
void RankBoostModel::update_model_weights_for_binary_feature(int best_binary_idx, 
															 weight_t alpha)
{
	binary_weights[best_binary_idx] += alpha;
	binary_update_counts[best_binary_idx]++;
	if (binary_update_counts[best_binary_idx] == 1)
	{
		non_zero_binary_idxs.push_back(best_binary_idx);
		sort(non_zero_binary_idxs.begin(),non_zero_binary_idxs.end());
	}
}



/********************************************************************************
Changes the current weights in the model.
Since this is a real theta thresholded variable, all weights above theta should
be affected. If the q_def is 1 then also the default (no vote) weights need to be
updated.
*********************************************************************************/
void RankBoostModel::update_model_weights_for_real_feature(weight_t alpha, 
			int best_real_idx, int q_def, int theata_idx_start, int theta_idx_end)
{
	if (theta_idx_end<0)
		theta_idx_end= real_weights[best_real_idx].size()-1;
	int i;
	for (i=theata_idx_start; i<=theta_idx_end; i++)
		real_weights[best_real_idx][i] += alpha;

	real_update_counts[best_real_idx]++;
	if (real_update_counts[best_real_idx] == 1)
	{
		non_zero_real_idxs.push_back(best_real_idx);
		sort(non_zero_real_idxs.begin(),non_zero_real_idxs.end());
	}

	if (q_def>0)
	{
		real_default_weights[best_real_idx] += alpha;
		total_default_weight += alpha;
	}
}


struct feature_pair {
	feature_pair() : idx(-1), score(NEG_INF) {};
	feature_pair(int _i, float _s) : idx(_i), score(_s) {};
	bool operator< (const feature_pair& other) const
	{
		return score>other.score;
	}
	int idx;
	float score;
};

/***************************************************************************

****************************************************************************/
void RankBoostModel::ouput_ranked_feature_list( ostream& os) const
{
	os << "FEATURE LIST FOR ROUND " << current_round << endl;
	
	if (binary_weights.size()>0)
	{
		os << "BINARY FEATURE WEIGHTS: " << endl;
		vector<feature_pair> bin_pairs;
		int i;
		for (i=0; i<binary_weights.size(); i++)
			if (binary_weights[i] != 0)
				bin_pairs.push_back(feature_pair(i,fabs(binary_weights[i])));

		sort(bin_pairs.begin(),bin_pairs.end());
		os << setprecision(7);
		for (i=0; i<bin_pairs.size(); i++)
		{
			os << i+1<< ")\t" << binary_weights[bin_pairs[i].idx] << "\t" << bin_pairs[i].idx << "\t" <<
				binary_feature_names[bin_pairs[i].idx] << " (" << binary_update_counts[bin_pairs[i].idx] << 
				" updates)" << endl;
		}
		os << endl;
	}


	if (real_weights.size()>0)
	{
		int i;
		os << "REAL FEATURE WEIGHTS: " << endl;
		vector<feature_pair> real_pairs;
		for (i=0; i<real_weights.size(); i++)
		{
			float max=0;
			int j;
			for (j=0; j<real_weights[i].size(); j++)
				if (real_weights[i][j] != 0.0)
					if (fabs(real_weights[i][j])>max)
						max=fabs(real_weights[i][j]);
			
			if (fabs(real_default_weights[i])>max)
				max = fabs(real_default_weights[i]);

			if (max == 0)
				continue;
					
			real_pairs.push_back(feature_pair(i,max));
		}

		sort(real_pairs.begin(),real_pairs.end());
		for (i=0; i<real_pairs.size(); i++)
		{
			int idx = real_pairs[i].idx;
			os << i+1 << ")\t" << idx << "\t" << real_feature_names[idx] << "  (" << 
				setprecision(5) << real_weights[idx].size()-1 << 
				" bins, " << real_update_counts[idx] << "  updates)" << endl;
			
			int j;
			for (j=1; j<real_weights[idx].size()-1; j++)
				if (real_weights[idx][j] != real_weights[idx][j+1])
					os << "  " << j << ":" << real_limits[idx][j] << "," << setprecision(4) << real_weights[idx][j];
				
			os << " >  " << "," << setprecision(4) << real_weights[idx][j] << endl;

			if (real_default_weights[idx] != 0)
				os << "default: " << real_default_weights[idx] << endl;
			
			os << endl;
		}
	}
}


struct FeatureStats {
	FeatureStats() : idx(NEG_INF), global_weight(0), local_weight(0), percent_active(0) {};

	bool operator< (const FeatureStats& other) const
	{
		return local_weight>other.local_weight;
	}

	int	   idx;
	double global_weight;
	double local_weight;
	double percent_active;
};



/******************************************************************************
Measures the "weight" of a feature, globally (how much at adds to all examples)
and locally (how much it adds to the samples for which it is applicable)
outputs a list ranked according to the local importanc)e
*******************************************************************************/
void RankBoostModel::ouput_importance_ranked_feature_list( const RankBoostDataset& training_ds, 
														  ostream& os,
														  int only_fidx,
														  int round_idx)
{
	remove_default_weights();

	if (real_weights.size()>0)
	{
		int i;

		// calc sample "weights"
		const vector<SamplePairWeight>& phi = training_ds.get_phi_support();
		const vector<RankBoostSample>& samples = training_ds.get_samples();
		vector<double> sam_weights;
		sam_weights.resize(samples.size(),0);

		for (i=0; i<phi.size(); i++)
		{
			sam_weights[phi[i].idx1]+=phi[i].weight;
			sam_weights[phi[i].idx2]+=phi[i].weight;
		}

		double total_weight =0;
		for (i=0; i<sam_weights.size(); i++)
			total_weight += sam_weights[i];
	
		// sum the weights
		const int num_real = real_weights.size();
		vector<FeatureStats> feature_stats;
		feature_stats.resize(num_real);
		for (i=0; i<samples.size(); i++)
		{
			const RankBoostSample& sam = samples[i];
			int j;
			for (j=0; j<sam.real_active_idxs.size(); j++)
			{
				const int f_idx   = sam.real_active_idxs[j];

				if (real_weights[f_idx].size() == 0)
					continue;

				const int bin_idx = get_real_bin_idx_for_value(f_idx,sam.real_active_values[j]);
				const double w    = real_weights[f_idx][bin_idx] - real_default_weights[f_idx];

				feature_stats[f_idx].local_weight   += (fabs(w) * sam_weights[i]);
				feature_stats[f_idx].percent_active +=  sam_weights[i];
			}
		}

		// remove default weight from all features
		bool changed_a_default = false;
		for (i=0; i<feature_stats.size(); i++)
		{
			if (real_default_weights[i] != 0)
			{
				int j;
				for (j=0 ;j<real_weights[i].size(); j++)
					real_weights[i][j] -= real_default_weights[i];
				real_default_weights[i]=0;
				changed_a_default = true;
			}
		}


		// recompute weight for all if needed
		if (changed_a_default)
		{
			feature_stats.clear();
			feature_stats.resize(num_real);

			int sam_idx;
			for (sam_idx=0; sam_idx<samples.size(); sam_idx++)
			{
				const RankBoostSample& sam = samples[sam_idx];
				int j;
				for (j=0; j<sam.real_active_idxs.size(); j++)
				{
					const int f_idx   = sam.real_active_idxs[j];
					const int bin_idx = get_real_bin_idx_for_value(f_idx,sam.real_active_values[j]);
					const double w    = real_weights[f_idx][bin_idx] - real_default_weights[f_idx];

					feature_stats[f_idx].local_weight += (fabs(w) * sam_weights[sam_idx]);
					feature_stats[f_idx].percent_active+= sam_weights[sam_idx];
				}
			}
		}

		// create global weighting and weighted active percent
		for (i=0; i<feature_stats.size(); i++)
		{

	//		cout << i << "\t" << fixed << setprecision(4) << feature_stats[i].local_weight << "\t" << 
	//		feature_stats[i].percent_active << endl;
			feature_stats[i].idx = i;
			if (feature_stats[i].local_weight>0)
			{
				feature_stats[i].global_weight = feature_stats[i].local_weight;
				feature_stats[i].global_weight /= total_weight;
				feature_stats[i].local_weight  /= feature_stats[i].percent_active;
				feature_stats[i].percent_active /= total_weight;
			}
		}

		sort(feature_stats.begin(),feature_stats.end());
		while (feature_stats.size()>0 && feature_stats[feature_stats.size()-1].local_weight == 0 )
			feature_stats.pop_back();

		os << "REAL FEATURE WEIGHTS: " << endl;
			
		for (i=0; i<feature_stats.size(); i++)
		{
			int idx = feature_stats[i].idx;

			if (only_fidx>=0 && idx != only_fidx)
				continue;

			os << i+1 << ")\t" << idx << "\t" << real_feature_names[idx] << "  (" << 
				setprecision(5) << real_weights[idx].size()-1 << 
				" bins, " << real_update_counts[idx];
			if (round_idx<0)
			{
				os << "  updates [ " << real_first_updates[idx] << " ]  )" << endl;
			}
			else
				os << " updates [ " << real_first_updates[idx] << " ],  " << round_idx << " rounds)" << endl;

			os << setprecision(4) << "LW: " <<  feature_stats[i].local_weight << "\tGW:" <<
				feature_stats[i].global_weight << "\t%ACT: " << feature_stats[i].percent_active << endl;
			
			int j;
			for (j=1; j<real_weights[idx].size()-1; j++)
				if (real_weights[idx][j] != real_weights[idx][j+1])
					os << "  " << j << ":" << real_limits[idx][j] << "," << setprecision(4) << real_weights[idx][j];
				
			os << " >  " << "," << setprecision(4) << real_weights[idx][j] << endl;

			if (real_default_weights[idx] != 0)
				os << "default: " << real_default_weights[idx] << endl;
			
			os << endl;
		}
	}
}




void RankBoostModel::get_top_misclassified_pairs(
								   const RankBoostDataset& training_ds,
								   const vector<weight_t>& D,
								   const vector<weight_t>& D0,
								   vector<idx_weight_pair>& pair_idxs,
								   int num_top_pairs ) const
{
	const int num_samples = training_ds.get_num_samples();
	const vector<RankBoostSample>& samples = training_ds.get_samples();

	vector<weight_t> rank_scores;

	rank_scores.resize(num_samples);

	int i;
	for (i=0; i<num_samples; i++)
		rank_scores[i]=calc_rank_score(samples[i]);

	vector<idx_weight_pair> pairs;
	const vector<SamplePairWeight>& phi_support = training_ds.get_phi_support();
	double train_error=0;
	for (i=0; i<phi_support.size(); i++)
	{
		if (rank_scores[phi_support[i].idx1]>=rank_scores[phi_support[i].idx2])
			pairs.push_back(idx_weight_pair(i,D[i]/D0[i]));
	}

	sort(pairs.begin(),pairs.end());


	pair_idxs.clear();
	for (i=0; i<pairs.size(); i++)
	{
		if (num_top_pairs>0 && i>=num_top_pairs)
			break;
		pair_idxs.push_back(pairs[i]);
	}
}

/**********************************************************************
***********************************************************************/
void RankBoostModel::print_top_misclassified_pairs(
								   const RankBoostDataset& training_ds,
								   const vector<weight_t>& D,
								   const vector<weight_t>& org_D,
								   int num_top_pairs,
								   ostream& os) const
{
	const int num_samples = training_ds.get_num_samples();
	const vector<RankBoostSample>& samples = training_ds.get_samples();

	vector<weight_t> rank_scores;

	rank_scores.resize(num_samples);

	int i;
	for (i=0; i<num_samples; i++)
		rank_scores[i]=calc_rank_score(samples[i]);

	vector<idx_weight_pair> pairs;
	const vector<SamplePairWeight>& phi_support = training_ds.get_phi_support();
	double train_error=0;
	for (i=0; i<phi_support.size(); i++)
	{
		if (rank_scores[phi_support[i].idx1]>=rank_scores[phi_support[i].idx2])
			pairs.push_back(idx_weight_pair(i,D[i]/org_D[i]));
	}

	sort(pairs.begin(),pairs.end());

	os << "Top miscalssified pairs: " << endl;
	for (i=0; i<num_top_pairs && i<pairs.size(); i++)
	{
		os << i << " " << pairs[i].idx << " " << pairs[i].weight << endl;
		const int idx1 = phi_support[pairs[i].idx].idx1;
		const int idx2 = phi_support[pairs[i].idx].idx2;

		os << idx1 << " > " << idx2 << endl;
	}
	os << endl;
}


/***************************************************************************

****************************************************************************/

struct score_pair {
	score_pair() : idx(int(NEG_INF)), score(NEG_INF) {};
	score_pair(int _i, float _n) : idx(_i), score(_n) {};
	bool operator< (const score_pair& other) const
	{
		return score>other.score;
	}
	int idx;
	float score;
};


void RankBoostModel::list_feature_vals_according_to_score(vector<RankBoostSample>& sams) const
{
	const int num_sams = sams.size();
	vector<bool> ind_printed;
	vector< float > total_scores;
	vector< vector<int> >   feature_idxs;
	vector< vector<float> > feature_values, feature_scores;

	ind_printed.resize(real_feature_names.size(),false);
	total_scores.resize(num_sams,0);
	feature_idxs.resize(num_sams);
	feature_values.resize(num_sams);
	feature_scores.resize(num_sams);

	vector<score_pair> pairs;

	int i;
	for (i=0; i<num_sams; i++)
		total_scores[i]=calc_rank_score_with_details(sams[i],feature_idxs[i],feature_values[i],feature_scores[i]);

	for (i=0; i<num_sams; i++)
	{
		int j;
		for (j=0; j<feature_idxs[i].size(); j++)
		{
			score_pair p;
			p.idx = feature_idxs[i][j];

			if (num_sams != 2)
			{
				p.score = fabs(feature_scores[i][j]);
			}
			else
			{
				int other_sam=1;
				if (i==1)
					other_sam=0;

				float other_score=0;
				int k;
				for (k=0; k<feature_idxs[other_sam].size(); k++)
					if (feature_idxs[other_sam][k] == feature_idxs[i][j])
					{
						other_score = feature_scores[other_sam][k];
						break;
					}
				p.score = fabs(feature_scores[i][j]-other_score);
			}
			pairs.push_back(p);
		}
		cout << i << setprecision(3) << fixed << "\t" << total_scores[i] << endl;
	}
	sort(pairs.begin(),pairs.end());

	cout << "\tf_idx";
	for (i=0; i<num_sams; i++)
	{
		cout << "\tval" << i << "\tscore" <<i;
	}

	if (num_sams == 2)
		cout << "\tdiff";
	cout << endl;

	int num_printed=0;
	for (i=0; i<pairs.size(); i++)
	{
		if (ind_printed[pairs[i].idx])
			continue;

		int idx = pairs[i].idx;
		ind_printed[idx]=true;

		cout << num_printed++ << "\t" << idx ;
		float v0=0,v1=0;
		int j;
		for (j=0; j<num_sams; j++)
		{
			int k;
			for (k=0; k<feature_idxs[j].size(); k++)
				if (feature_idxs[j][k]==idx)
					break;

			if (k==feature_idxs[j].size())
			{
				cout << "\t  -\t  0";
			}
			else
			{
				cout << "\t" << feature_values[j][k] << "\t" << feature_scores[j][k];
				if (j==0)
				{
					v0=feature_scores[j][k];
				}
				else
					v1=feature_scores[j][k];
			}
		}
		if (num_sams==2)
		{
			cout << "\t" << v0-v1;
		}
		cout << "\t" << real_feature_names[idx] << endl;
		

		if (num_printed>=40)
			break;
	}
	
}

