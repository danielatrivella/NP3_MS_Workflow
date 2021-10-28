#include "RankBoost.h"


bool RankBoostSample::get_feature_val(int f_idx, float *val) const
{
	int f;
	for (f=0; f<real_active_idxs.size(); f++)
	{
		if (real_active_idxs[f]==f_idx)
		{
			*val = real_active_values[f];
			return true;
		}
	}
	return false;
}

void RankBoostSample::print(ostream& os) const
{
	int i;

	cout << "GROUP IDX    : " << this->groupIndex << endl;
	cout << "RANK IN GROUP: " << this->rank_in_group << endl;
	cout << "TAG1: " << tag1 << endl;
	cout << "TAG2: " << tag2 << endl;
	cout << "TAG3: " << tag3 << endl;

	if (binary_novote_idxs.size()>0)
	{
		os << "BINARY NOVOTE IDXS: " << binary_novote_idxs.size() << " : ";
		for (i=0; i<binary_novote_idxs.size(); i++)
			os << " " << binary_novote_idxs[i];
		os << endl;
	}

	if (binary_non_zero_idxs.size())
	{
		os << "BINARY NON_ZERO IDXS: " << binary_non_zero_idxs.size() << " : ";
		for (i=0; i<binary_non_zero_idxs.size(); i++)
			os << " " << binary_non_zero_idxs[i];
		os << endl;
	}

	if (real_novote_idxs.size()>0)
	{
		os << "REAL INACTIVE IDXS  : " << real_novote_idxs.size() << " : ";
		for (i=0; i<real_novote_idxs.size(); i++)
			os << " " << real_novote_idxs[i];
		os << endl;
	}
	
	if (real_active_idxs.size()>0)
	{
		os << "REAL ACTIVE IDXS  : " << real_active_idxs.size() << " : ";
		for (i=0; i<real_active_idxs.size(); i++)
			os << " " << real_active_idxs[i] << "," << real_active_values[i];
		os << endl;
	}
}



void RankBoostDataset::clear()
{
	num_groups=0;
	total_phi_weight=0;
	samples.clear();
	phi_support.clear();
	ahead_lists.clear();
	behind_lists.clear();
}


void RankBoostDataset::compute_total_phi_weight()
{
	total_phi_weight=0;
	int i;

	for (i=0; i<phi_support.size(); i++)
		if (phi_support[i].weight>0)
			total_phi_weight+=phi_support[i].weight;
}

/*****************************************************************************
Creates two list for each feature idx what sample pairs in phi:
- idxs of pairs ordered correctly, all (x',x) for which f(x)=1 and f(x')=0 (or no vote)
- idxs of pairs ordered incorreclty,  all (x',x) for which f(x')=1 and f(x)=0 (or no vote)
******************************************************************************/
void RankBoostDataset::initialize_binary_ordered_phi_lists(const vector<string>* binary_names)
{
	const int num_pairs = phi_support.size();
	const int num_binary_features = binary_one_lists.size();

	vector<double> weight_ordered_correctly, weight_ordered_incorrectly;


	if (num_binary_features == 0)
		return;

	if (num_pairs == 0)
	{
		cout << "Error: must first initialize phi support and binary features!" << endl;
		exit(1);
	}

	binary_pairs_ordered_correctly.resize(num_binary_features);
	binary_pairs_ordered_incorrectly.resize(num_binary_features);

	weight_ordered_correctly.resize(num_binary_features,0);
	weight_ordered_incorrectly.resize(num_binary_features,0);

	int i;
	for (i=0; i<num_pairs; i++)
	{
		const int sam_idx1 = phi_support[i].idx1;
		const int sam_idx2 = phi_support[i].idx2;
		const weight_t pair_weight = phi_support[i].weight;
		const vector<int>& ones_idx1 = samples[sam_idx1].binary_non_zero_idxs;
		const vector<int>& ones_idx2 = samples[sam_idx2].binary_non_zero_idxs;

		if (binary_names)
		{
		//	if (ones_idx1.size() == 0)
		//		cout << "Warning: zero binary feautres with 1 for sample " << sam_idx1 << endl;

		//	if (ones_idx2.size() == 0)
		//		cout << "Warning: zero binary feautres with 1 for sample " << sam_idx2 << endl;


			int j;
			for (j=1; j<ones_idx1.size(); j++)
				if (ones_idx1[j]<=ones_idx1[j-1])
				{
					cout << "Error in the ordereing of feature idxs in phi pair " << i << endl;
					cout << ones_idx1[j] << "\t" << binary_names->at(ones_idx1[j]) << "\t" <<
						ones_idx1[j-1] << "\t" << binary_names->at(ones_idx1[j-1]) << endl;
					
					exit(1);
				}
		}


		int i1=0,i2=0;

		while (i1<ones_idx1.size() || i2<ones_idx2.size())
		{
			if (i1==ones_idx1.size())
			{
				const int feature_idx = ones_idx2[i2];
				binary_pairs_ordered_correctly[feature_idx].push_back(i);
				weight_ordered_correctly[feature_idx] += pair_weight; 
				i2++;
				continue;
			}

			if (i2==ones_idx2.size())
			{
				const int feature_idx = ones_idx1[i1];
				binary_pairs_ordered_incorrectly[feature_idx].push_back(i);
				weight_ordered_incorrectly[feature_idx] += pair_weight; 
				i1++;
				continue;
			}

			const int feature_idx1 = ones_idx1[i1];
			const int feature_idx2 = ones_idx2[i2];

			if (feature_idx1 == feature_idx2)
			{
				i1++;
				i2++;
				continue;
			}

			if (feature_idx1<feature_idx2)
			{
				binary_pairs_ordered_incorrectly[feature_idx1].push_back(i);
				weight_ordered_incorrectly[feature_idx1] += pair_weight; 
				i1++;
				continue;
			}
			else
			{
				binary_pairs_ordered_correctly[feature_idx2].push_back(i);
				weight_ordered_correctly[feature_idx2] += pair_weight; 
				i2++;
				continue;
			}
		}
	}

	if (binary_names)
	{
		int i;
		cout << "BINARY PHI LIST REPORT: " << endl;
		for (i=0; i<num_binary_features; i++)
		{
			cout << i << "\t" << binary_pairs_ordered_correctly[i].size() << "\t" <<
				binary_pairs_ordered_incorrectly[i].size() << "\t" <<
				phi_support.size() -
				binary_pairs_ordered_correctly[i].size() - 
				binary_pairs_ordered_incorrectly[i].size()  << "\t";

			if (weight_ordered_correctly[i]>0 && weight_ordered_incorrectly[i]>0)
			{
				double ratio = weight_ordered_correctly[i]/weight_ordered_incorrectly[i];
				cout << ratio << "\t";
			}
			else
				cout << "-\t";
			cout << binary_names->at(i) << endl;
		}
	}

}


/*****************************************************************************
Creates a large featureXsamples table for all real values (NEG_INF = no vote)
******************************************************************************/
void RankBoostDataset::initialzie_real_feature_table(int num_real_features)
{
	const int num_samples = samples.size(); 

	real_feature_values.clear();
	real_feature_values.resize(num_real_features);

	int i;
	for (i=0; i<num_samples; i++)
	{
		const RankBoostSample& sam = samples[i];
		const int num_active_idxs = sam.real_active_idxs.size();
		int j;
		for (j=0; j<num_active_idxs; j++)
		{
			const int& f_idx = sam.real_active_idxs[j];
			real_feature_values[f_idx][i]=sam.real_active_values[j];
		}
	}
}


/*****************************************************************************
Computes an initial distribution. Probs are set so the sum of all pairs in phi
with positive weight equals 1.
******************************************************************************/
void RankBoostDataset::compute_initial_distribution(vector<weight_t>& d) const
{
	int i;
	weight_t norm_val = 1.0 / total_phi_weight;

	d.resize(phi_support.size(),0);
	for (i=0; i<phi_support.size(); i++)
		if (phi_support[i].weight>0)
			d[i] = phi_support[i].weight * norm_val;

}


/*****************************************************************************
Computes the maximal weights for which each pair will get a normal update.
******************************************************************************/
void RankBoostDataset::compute_max_weights_for_normal_update(
					const vector<weight_t>& D0,
					vector<weight_t>& d) const
{
	int i;
	d.resize(D0.size(),0);
	for (i=0; i<D0.size(); i++)
		d[i] = D0[i] * max_ratio_for_regular_update;
}



/***************************************************************************
Creates the ahead and behind lists for each sample (these can be used to 
calculate the potential quickly). This method is useful if these lists are
much smaller than the number of samples (|phi|<<n^2).
****************************************************************************/
void RankBoostDataset::initialize_potenital_lists()
{
	ahead_lists.clear();
	behind_lists.clear();

	ahead_lists.resize(samples.size());
	behind_lists.resize(samples.size());

	int i;
	for (i=0; i<phi_support.size(); i++)
	{
		ahead_lists[phi_support[i].idx2].push_back(i);
		behind_lists[phi_support[i].idx1].push_back(i);
	}
}


/**************************************************************************
Creates a list for each feature of the idxs of samples which have a 1 at 
a certain binary feature.
***************************************************************************/
void RankBoostDataset::initialize_binary_one_lists(int num_binary_features)
{
	const int num_samples = samples.size();
	int i;

	binary_one_lists.clear();
	binary_one_lists.resize(num_binary_features);

	for (i=0; i<num_samples; i++)
	{
		const vector<int>& idxs = samples[i].binary_non_zero_idxs;
		int j;
		for (j=0; j<idxs.size(); j++)
			binary_one_lists[idxs[j]].push_back(i);
	}
}



/***********************************************************************
Creates lists of samples that vote in each of the real feature's values
************************************************************************/
void RankBoostDataset::initialize_real_vote_lists(const RankBoostModel& rbm)
{
	const int num_samples = samples.size();
	const int num_real_features = rbm.get_num_real_features();
	int i;
	
	real_vote_lists.clear();
	real_vote_lists.resize(num_real_features);
	
	ind_all_samples_vote.resize(num_real_features,false);

	vector< vector<int> > counts;
	counts.resize(num_real_features);

	for (i=0; i<num_real_features; i++)
		counts[i].resize(rbm.get_num_bins_for_real_feature(i),0);

	for (i=0; i<num_samples; i++)
	{
		const vector<int>&   idxs = samples[i].real_active_idxs;
		const vector<float>& vals = samples[i].real_active_values;
		int j;

		for (j=0; j<idxs.size(); j++)
		{
			const int bin_idx = rbm.get_real_bin_idx_for_value(idxs[j],vals[j]);
			counts[idxs[j]][bin_idx]++;
		}

	}

	// resize vectors according to counts
	int f;
	for (f=0; f<num_real_features; f++)
	{
		int total_counts=0;
		int i;
		for (i=0; i<counts[f].size(); i++)
			total_counts+=counts[f][i];

		real_vote_lists[f].resize(counts[f].size());
		
		for (i=0; i<counts[f].size(); i++)
			real_vote_lists[f][i].reserve(counts[f][i]);
	}

	for (i=0; i<num_samples; i++)
	{
		const vector<int>&   idxs = samples[i].real_active_idxs;
		const vector<float>& vals = samples[i].real_active_values;
		int j;

		for (j=0; j<idxs.size(); j++)
		{
			const int feature_idx = idxs[j];
			const int bin_idx = rbm.get_real_bin_idx_for_value(feature_idx,vals[j]);
			if (! ind_all_samples_vote[feature_idx])
				real_vote_lists[feature_idx][bin_idx].push_back(i);
		}
	}
}



/***********************************************************************
Calcs the potential of each sample x (this is the difference between
the weight of the pairs in which x is ahead minus the weight of the
samples in which x is behind).
************************************************************************/
void RankBoostDataset::calc_potentials(const vector<weight_t>& D, vector<weight_t>& potentials) const
{
	const int num_samples = samples.size();
	int i;
	potentials.resize(num_samples);
	for (i=0; i<num_samples; i++)
	{
		const vector<int>& ahead_list  = ahead_lists[i];
		const vector<int>& behind_list = behind_lists[i];
		weight_t potential = 0;
		
		int j;
		for (j=0; j<ahead_list.size(); j++)
			potential+=D[ahead_list[j]];

		for (j=0; j<behind_list.size(); j++)
			potential-=D[behind_list[j]];

		potentials[i]=potential;
	}
}


/***************************************************************************************
Updates the weights of samples affected by the current boosting round.
Weights of pairs correctly ordered should go down, and weights of pairs incorrectly 
ordered should go up.

The function returns the normalizing constant Z.
to avoid numerical stability issues, some of the summations are done in batches.
****************************************************************************************/
double RankBoostDataset::update_dsitribution_according_to_binary_feature(
											int binary_feature_idx,
											weight_t alpha,
											vector<weight_t>& D,
											bool verbose) const
{
	const vector<int>& correct_pairs = binary_pairs_ordered_correctly[binary_feature_idx];
	const vector<int>& incorrect_pairs = binary_pairs_ordered_incorrectly[binary_feature_idx];
	const weight_t norm_correct   = exp(-alpha);
	const weight_t norm_incorrect = exp(alpha);
	const int add_batch_size = 1000;

	double mod_weight_before =0;
	double mod_weight_after  =0;

	double batch_before  =0;
	double batch_after   =0;
	int batch_count =0;

	int i;
	for (i=0; i<correct_pairs.size(); i++)
	{
		weight_t& pair_weight = D[correct_pairs[i]];

		batch_before += pair_weight;
		pair_weight   *= norm_correct;
		batch_after  += pair_weight;

		batch_count++;
		if (batch_count == add_batch_size)
		{
			mod_weight_before += batch_before;
			mod_weight_after  += batch_after;
			batch_count=0;
			batch_before=0;
			batch_after=0;
		}
	}

	mod_weight_before += batch_before;
	mod_weight_after  += batch_after;

	batch_count=0;
	batch_before=0;
	batch_after=0;

	for (i=0; i<incorrect_pairs.size(); i++)
	{
		weight_t& pair_weight = D[incorrect_pairs[i]];

		batch_before += pair_weight;
		pair_weight   *= norm_incorrect;
		batch_after  += pair_weight;

		batch_count++;
		if (batch_count == add_batch_size)
		{
			mod_weight_before += batch_before;
			mod_weight_after  += batch_after;
			batch_count=0;
			batch_before=0;
			batch_after=0;
		}
	}

	mod_weight_before += batch_before;
	mod_weight_after  += batch_after;

	double total_weight_after = (1.0 - mod_weight_before) + mod_weight_after;

	return total_weight_after;
}


/************************************************************************************
Update the weights of samples affected in current boosting round.
Examines all samples (including no vote ones)
to avoid numerical stability issues, some of the summations are done in batches.

use POS_INF if there is no theta end
*************************************************************************************/
double RankBoostDataset::update_distribution_according_to_real_feature(
														int best_real_idx, 
														float theta_start, 
														float theta_end, 
														int q_def, 
														weight_t alpha, 
														vector<weight_t>& D,
														vector<weight_t> *max_D_for_normal_updates,
														bool verbose) const
{
	const map<int,float>& feature_values = real_feature_values[best_real_idx];
	const int add_batch_size = 1000;
	const weight_t norm_correct   = exp(-alpha);
	const weight_t norm_incorrect = exp(alpha);

	double mod_weight_before = 0; // the weight of the modified samples
	double mod_weight_after  = 0; //

	double batch_before  =0;
	double batch_after   =0;
	int num_diff = 0;
	int batch_count =0;
	int i;

	for (i=0; i<phi_support.size(); i++)
	{
		const int x0 = phi_support[i].idx1;
		const int x1 = phi_support[i].idx2;

		float x0_val = NEG_INF, x1_val = NEG_INF;
		bool x0_vote =true, x1_vote=true;
		
		map<int,float>::const_iterator it0 = feature_values.find(x0);
		map<int,float>::const_iterator it1 = feature_values.find(x1);

		if ( it0 != feature_values.end())
		{
			x0_val = it0->second;
		}
		else
			x0_vote = false;
		
		if ( it1 != feature_values.end())
		{
			x1_val = it1->second;
		}
		else
			x1_vote = false;

		int h0 = q_def, h1 = q_def;

		if (x0_vote)
			h0 = (x0_val>theta_start ? 1 : 0);
		
		if (x1_vote)
			h1 = (x1_val>theta_start ? 1 : 0);

		// check if the vals exceed the maximum theta (only if an upper bound is used
		if (x0_val > theta_end)
			h0 = 0;
		if (x1_val > theta_end)
			h1 = 0;
		
		if (h0 == h1)
			continue;

		num_diff++;

		weight_t& pair_weight = D[i];
		batch_before += pair_weight;

		if (! max_D_for_normal_updates || pair_weight < max_D_for_normal_updates->at(i))
		{
			pair_weight  *= (h1>0 ? norm_correct : norm_incorrect);
		//	cout << pair_weight << " <-> " << max_D_for_normal_updates->at(i) << endl;
		}
		else
		{
			double update_weight = (h1>0 ? norm_correct : norm_incorrect);
			
			if (update_weight<=1.0)
			{
				pair_weight *= update_weight;
			}
			else // this update is not as strong
			{
				double delta = update_weight - 1.0;
				delta *= (max_D_for_normal_updates->at(i)/pair_weight);
				pair_weight *= (1.0+delta);
			}
		}
	
		batch_after  += pair_weight;
	
		batch_count++;
		if (batch_count == add_batch_size)
		{
			mod_weight_before += batch_before;
			mod_weight_after  += batch_after;
			batch_count=0;
			batch_before=0;
			batch_after=0;
		}
	}

	mod_weight_before += batch_before;
	mod_weight_after  += batch_after;

	double total_weight_after = (1.0 - mod_weight_before) + mod_weight_after;

	if (verbose)
	{
		cout << setprecision(7);
		cout << "Theta: " << theta_start;
		if (theta_end<POS_INF) 
			cout << "-" << theta_end;
		cout << " #diff " << num_diff << endl;
		cout << "mod weight " << mod_weight_before << " -> " << mod_weight_after << "  Z=" << total_weight_after << endl;

	}

	return total_weight_after;
}


