#include "PeakRankModel.h"

extern const int num_RKH_combos = 6;
extern const int num_RKH_pairs = num_RKH_combos * num_RKH_combos;
extern const string combos[]={"NB","H","K","HK","R","RR"}; 

int calc_RKH_combo_idx (int r, int k, int h)
{
	const int sum = r + k + h;

	if (sum==0)
		return 0; // no basics

	if (sum == 1)
	{
		if (h)
			return 1; // H only
		if (k)
			return 2; // K only
	}

	if (r==0)     // KH, HH, KK ...
		return 3;

	if (r==1)	  // R, RK , RH ...
		return 4;

	return 5; // RR, RRK, RRR....
}

extern const float RKH_pair_matrix[6][6]={
	{  0.0,    1.0,   2.0,    3.0,   4.0,   5.0},
	{ -1.0,    0.0,   2.0,    3.0,   4.0,   5.0},
	{ -2.0,   -2.0,   0.0,    2.0,   4.0,   5.0},
	{ -3.0,   -3.0,  -2.0,    0.0,   3.0,   4.0},
	{ -4.0,   -4.0,  -4.0,   -3.0,   0.0,   4.0},
	{ -5.0,   -5.0,  -5.0,   -4.0,  -4.0,   0.0} };

/**********************************************************************
This fills in all possible sequence features that can be used in any
of the models (these features are good for all fragments/mobility/
size/charge). The models can later choose to ignore some of these
features by asigning the weight 0. All features are derived directly
from the peptide sequence and concern the specified cut idx.
***********************************************************************/
void PeakRankModel::fill_advanced_peak_features(
								 const  vector<int>& org_amino_acids,
								 int    cut_idx,
								 mass_t cut_mass,
								 mass_t pm_with_19,
								 int	spec_charge,
								 const FragmentType& fragment,
								 RankBoostSample& sample) const
{
	const int length = org_amino_acids.size();
	const int num_aas = model_aa_labels.size();
	int rankFeatureIdx=0;
	int i;

	vector<int> amino_acids;
	convert_aas_to_model_aas(org_amino_acids, amino_acids);

	if (amino_acids.size() != org_amino_acids.size())
	{
		cout << "Error: aa size mismatch!" << endl;
		exit(1);
	}

	if (cut_idx<=0 || cut_idx>=amino_acids.size())
	{
		cout << "Error: cut_idx is bad!" << endl;
		exit(1);
	}

	// need to use the special Idx variables and not the regular enumerations
	const int HisIdx = session_aas_to_model_aas[His];
	const int LysIdx = session_aas_to_model_aas[Lys];
	const int ArgIdx = session_aas_to_model_aas[Arg];
	const int SerIdx = session_aas_to_model_aas[Ser];
	const int ThrIdx = session_aas_to_model_aas[Thr];
	const int ProIdx = session_aas_to_model_aas[Pro];
	const int GlyIdx = session_aas_to_model_aas[Gly];
	const int AlaIdx = session_aas_to_model_aas[Ala];
	const int LeuIdx = session_aas_to_model_aas[Leu];
	const int AsnIdx = session_aas_to_model_aas[Asn];
	const int AspIdx = session_aas_to_model_aas[Asp];
	const int GluIdx = session_aas_to_model_aas[Glu];


	sample.clear();

	// special N C side aa indicators
	int num_nH=0, num_cH=0;
	int num_nK=0, num_cK=0;
	int num_nR=0, num_cR=0;
	
	for (i=0; i<cut_idx; i++)
	{
		if (amino_acids[i] == HisIdx)
			num_nH++;

		if (amino_acids[i] == LysIdx)
			num_nK++;

		if (amino_acids[i] == ArgIdx)
			num_nR++;
	}

	for (i=cut_idx; i<length; i++)
	{
		if (amino_acids[i] == HisIdx)
			num_cH++;

		if (amino_acids[i] == LysIdx)
			num_cK++;

		if (amino_acids[i] == ArgIdx)
			num_cR++;
	}

	// MASS / LOCATION FEATURES (REAL + BINARY)
	const mass_t exp_peak_mass = fragment.calc_expected_mass(cut_mass,pm_with_19);
	const mass_t min_obs_mass = calc_min_detected_mass(pm_with_19,spec_charge);
	const mass_t max_obs_mass = (pm_with_19>max_detected_mass ? max_detected_mass : pm_with_19);
	const float peak_mass_prop = ((exp_peak_mass - min_obs_mass)/(max_obs_mass - min_obs_mass));
	const float rounded_peak_prop = 0.1*floor(peak_mass_prop * 10.0);

	// give values within a resolution of 20 Da
	const mass_t dis_from_min = 25.0*floor((exp_peak_mass - min_obs_mass)*0.04);
	const mass_t dis_from_max = 25.0*floor((max_obs_mass  - exp_peak_mass)*0.04);

	const int RKH_n_combo_idx = calc_RKH_combo_idx(num_nR,num_nK,num_nH);
	const int RKH_c_combo_idx = calc_RKH_combo_idx(num_cR,num_cK,num_cH);

	const int RKH_pair_idx = (RKH_n_combo_idx * num_RKH_combos) + RKH_c_combo_idx;

	const float RKH_liniar_pair_idx = RKH_pair_matrix[RKH_n_combo_idx][RKH_c_combo_idx];
	const int n_aa = amino_acids[cut_idx-1];
	const int c_aa = amino_acids[cut_idx];

	// proportion of mass of the N/C fragments (special values are given to the first 3 
	// cuts on each side. If the cut is not in those regions, prop is assigned the
	// index of the fifth in which it falls 
	int side_length=3;
	if (length>=12) side_length=4;
	if (length>=15) side_length=5;

	float cut_prop;
	if (cut_idx<=side_length)
	{
		cut_prop=(float)cut_idx;
	}
	else if (cut_idx>=length-side_length)
	{
		cut_prop=(float)(11+cut_idx-length);
	}
	else
	{
		cut_prop = 5.1+floor(3.0*(cut_mass/pm_with_19))*0.1;
	}


	// fill N RKH and C RKH
	sample.add_real_feature(rankFeatureIdx++,RKH_n_combo_idx);
	sample.add_real_feature(rankFeatureIdx++,RKH_c_combo_idx);

	// peak prop
	sample.add_real_feature(rankFeatureIdx++,rounded_peak_prop);
	sample.add_real_feature(rankFeatureIdx+RKH_pair_idx,rounded_peak_prop);
	rankFeatureIdx+=num_RKH_pairs;
	
	// fill dis features
	if (dis_from_min<dis_from_max)
	{
		sample.add_real_feature(rankFeatureIdx,dis_from_min);
		rankFeatureIdx++;
		sample.add_real_feature(rankFeatureIdx+RKH_pair_idx,dis_from_min);
		rankFeatureIdx+=(2*num_RKH_pairs+1);
	}
	else
	{
		rankFeatureIdx+=(num_RKH_pairs+1);
		sample.add_real_feature(rankFeatureIdx,dis_from_max);
		rankFeatureIdx++;
		sample.add_real_feature(rankFeatureIdx+RKH_pair_idx,dis_from_max);
		rankFeatureIdx+=num_RKH_pairs;
	} 

	//  fill prop features
	sample.add_real_feature(rankFeatureIdx++,cut_prop);
	sample.add_real_feature(rankFeatureIdx+RKH_pair_idx,cut_prop);
	rankFeatureIdx+=num_RKH_pairs;

	// fill prop X dis features
	if (dis_from_min<dis_from_max)
	{
		if (dis_from_min<75.0)
			sample.add_real_feature(rankFeatureIdx,cut_prop);
		rankFeatureIdx++;
		if (dis_from_min<150.0)
			sample.add_real_feature(rankFeatureIdx,cut_prop);
		rankFeatureIdx++;
		if (dis_from_min<250.0)
			sample.add_real_feature(rankFeatureIdx,cut_prop);
		rankFeatureIdx++;
		if (dis_from_min<=400.0)
			sample.add_real_feature(rankFeatureIdx,cut_prop);
		rankFeatureIdx++;

		rankFeatureIdx+=4;
	}
	else
	{
		rankFeatureIdx+=4;
		if (dis_from_max<75.0)
			sample.add_real_feature(rankFeatureIdx,cut_prop);
		rankFeatureIdx++;
		if (dis_from_max<150.0)
			sample.add_real_feature(rankFeatureIdx,cut_prop);
		rankFeatureIdx++;
		if (dis_from_max<250.0)
			sample.add_real_feature(rankFeatureIdx,cut_prop);
		rankFeatureIdx++;
		if (dis_from_max<=400.0)
			sample.add_real_feature(rankFeatureIdx,cut_prop);
		rankFeatureIdx++;
	}

	// fill aa count features (up top 3 aa's away from cut)
	vector<int> n_aa_counts, c_aa_counts;
	n_aa_counts.resize(num_aas+1,0);
	c_aa_counts.resize(num_aas+1,0);

	for (i=0; i<cut_idx-3; i++)
		n_aa_counts[amino_acids[i]]++;

	for (i=cut_idx+3; i<length; i++)
		c_aa_counts[amino_acids[i]]++;

	int a;
	for (a=0; a<num_aas; a++)
		sample.add_real_feature(rankFeatureIdx++,n_aa_counts[a]);
	
	for (a=0; a<num_aas; a++)
		sample.add_real_feature(rankFeatureIdx++,c_aa_counts[a]);

	// including the aas up top the count
	int start_cut = cut_idx-3;
	if (start_cut<0)
		start_cut=0;
	for (i=start_cut; i<cut_idx; i++)
		n_aa_counts[amino_acids[i]]++;

	int end_cut = cut_idx+3;
	if (end_cut>length)
		end_cut = length;
	for (i=cut_idx; i<end_cut; i++)
		c_aa_counts[amino_acids[i]]++;

	for (a=0; a<num_aas; a++)
		sample.add_real_feature(rankFeatureIdx++,n_aa_counts[a]);
	
	for (a=0; a<num_aas; a++)
		sample.add_real_feature(rankFeatureIdx++,c_aa_counts[a]);

	// fill aa flanking features N side
	if (cut_idx>0)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;
	if (cut_idx>1)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],cut_prop);
	rankFeatureIdx+=num_aas;
	if (cut_idx>2)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-3],cut_prop);
	rankFeatureIdx+=num_aas;

	// fill aa flanking features C side
	if (cut_idx<length)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;
	if (cut_idx<length-1)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],cut_prop);
	rankFeatureIdx+=num_aas;
	if (cut_idx<length-2)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+2],cut_prop);
	rankFeatureIdx+=num_aas;

	// fill cut pair features X-Y
	sample.add_real_feature(rankFeatureIdx+(n_aa*num_aas+c_aa),cut_prop);
	rankFeatureIdx+=(num_aas*num_aas);

	// fill cut pair features X-Y
	sample.add_real_feature(rankFeatureIdx+(n_aa*num_aas+c_aa),rounded_peak_prop);
	rankFeatureIdx+=(num_aas*num_aas);

	// fill flanking aa info with RKH_pair data
	// fill aa flanking features N side
	if (cut_idx>0)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],RKH_liniar_pair_idx);
	rankFeatureIdx+=num_aas;
	if (cut_idx>1)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],RKH_liniar_pair_idx);
	rankFeatureIdx+=num_aas;
	if (cut_idx>2)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-3],RKH_liniar_pair_idx);
	rankFeatureIdx+=num_aas;

	// fill aa flanking features C side
	if (cut_idx<length)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],RKH_liniar_pair_idx);
	rankFeatureIdx+=num_aas;
	if (cut_idx<length-1)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],RKH_liniar_pair_idx);
	rankFeatureIdx+=num_aas;
	if (cut_idx<length-2)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+2],RKH_liniar_pair_idx);
	rankFeatureIdx+=num_aas;
	

	// fill flanking aa info with peak prop data
	// fill aa flanking features N side
	if (cut_idx>0)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],rounded_peak_prop);
	rankFeatureIdx+=num_aas;
	if (cut_idx>1)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],rounded_peak_prop);
	rankFeatureIdx+=num_aas;
	if (cut_idx>2)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-3],rounded_peak_prop);
	rankFeatureIdx+=num_aas;

	// fill aa flanking features C side
	if (cut_idx<length)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],rounded_peak_prop);
	rankFeatureIdx+=num_aas;
	if (cut_idx<length-1)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],rounded_peak_prop);
	rankFeatureIdx+=num_aas;
	if (cut_idx<length-2)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+2],rounded_peak_prop);
	rankFeatureIdx+=num_aas;

	// add features for flanking pairs of amino acids
	if (cut_idx>1)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2]*num_aas+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas*num_aas;

	if (cut_idx<length-2)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx]*num_aas+amino_acids[cut_idx+1],cut_prop);
	rankFeatureIdx+=num_aas*num_aas;

	
	// X != R
	// features of the form  |LXK   |LXXK   |LXXXK   |LXXXXK
	if (cut_idx < length-2 &&
		amino_acids[cut_idx+2]==LysIdx &&
		amino_acids[cut_idx] != LysIdx && amino_acids[cut_idx] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-3 &&
		amino_acids[cut_idx+3]==LysIdx &&
		amino_acids[cut_idx] != LysIdx && amino_acids[cut_idx] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-4 &&
		amino_acids[cut_idx+4]==LysIdx &&
		amino_acids[cut_idx] != LysIdx && amino_acids[cut_idx] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-5 &&
		amino_acids[cut_idx+5]==LysIdx &&
		amino_acids[cut_idx] != LysIdx && amino_acids[cut_idx] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx && amino_acids[cut_idx+4] != ArgIdx )
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	// features of the form  |XLXK  |XLXXK  |XLXXXK  |XLXXXXK
	if (cut_idx < length-3 &&
		amino_acids[cut_idx+3]==LysIdx &&
		amino_acids[cut_idx+1] != LysIdx && amino_acids[cut_idx+1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx 
		&& amino_acids[cut_idx] != LysIdx )
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-4 &&
		amino_acids[cut_idx+4]==LysIdx &&
		amino_acids[cut_idx+1] != LysIdx && amino_acids[cut_idx+1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx && amino_acids[cut_idx] != LysIdx )
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-5 &&
		amino_acids[cut_idx+5]==LysIdx &&
		amino_acids[cut_idx+1] != LysIdx && amino_acids[cut_idx+1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx && amino_acids[cut_idx+4] != ArgIdx 
		&& amino_acids[cut_idx] != LysIdx )
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-6 &&
		amino_acids[cut_idx+6]==LysIdx &&
		amino_acids[cut_idx+1] != LysIdx && amino_acids[cut_idx+1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx && amino_acids[cut_idx+4] != ArgIdx &&
		amino_acids[cut_idx+5] != ArgIdx && amino_acids[cut_idx] != LysIdx )
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],cut_prop);
	rankFeatureIdx+=num_aas;


	// features of the form L|XK   L|XXK   L|XXXK   L|XXXXK
	if (cut_idx>0 &&cut_idx < length-1 &&
		amino_acids[cut_idx+1]==LysIdx &&
		amino_acids[cut_idx-1] != LysIdx && amino_acids[cut_idx-1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>0 &&cut_idx < length-2 &&
		amino_acids[cut_idx+2]==LysIdx &&
		amino_acids[cut_idx-1] != LysIdx && amino_acids[cut_idx-1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx && amino_acids[cut_idx+1] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>0 &&cut_idx < length-3 &&
		amino_acids[cut_idx+3]==LysIdx &&
		amino_acids[cut_idx-1] != LysIdx && amino_acids[cut_idx-1] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>0 && cut_idx < length-4 &&
		amino_acids[cut_idx+4]==LysIdx &&
		amino_acids[cut_idx-1] != LysIdx && amino_acids[cut_idx-1] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx && amino_acids[cut_idx] != ArgIdx )
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	
	// features of the form  |LXR   |LXXR   |LXXXR   |LXXXXR
		if (cut_idx < length-2 &&
		amino_acids[cut_idx+2]==ArgIdx &&
		amino_acids[cut_idx] != LysIdx && amino_acids[cut_idx] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-3 &&
		amino_acids[cut_idx+3]==ArgIdx &&
		amino_acids[cut_idx] != LysIdx && amino_acids[cut_idx] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-4 &&
		amino_acids[cut_idx+4]==ArgIdx &&
		amino_acids[cut_idx] != LysIdx && amino_acids[cut_idx] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-5 &&
		amino_acids[cut_idx+5]==ArgIdx &&
		amino_acids[cut_idx] != LysIdx && amino_acids[cut_idx] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx && amino_acids[cut_idx+4] != ArgIdx )
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	// features of the form  |XLXR  |XLXXR  |XLXXXR  |XLXXXXR
		if (cut_idx < length-3 &&
		amino_acids[cut_idx+3]==ArgIdx &&
		amino_acids[cut_idx+1] != LysIdx && amino_acids[cut_idx+1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx && amino_acids[cut_idx] != LysIdx &&
		amino_acids[cut_idx+2] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-4 &&
		amino_acids[cut_idx+4]==ArgIdx &&
		amino_acids[cut_idx+1] != LysIdx && amino_acids[cut_idx+1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx && amino_acids[cut_idx] != LysIdx )
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-5 &&
		amino_acids[cut_idx+5]==ArgIdx &&
		amino_acids[cut_idx+1] != LysIdx && amino_acids[cut_idx+1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx && amino_acids[cut_idx+4] != ArgIdx &&
		amino_acids[cut_idx] != LysIdx )
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx < length-6 &&
		amino_acids[cut_idx+6]==ArgIdx &&
		amino_acids[cut_idx+1] != LysIdx && amino_acids[cut_idx+1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx && amino_acids[cut_idx+4] != ArgIdx &&
		amino_acids[cut_idx+5] != ArgIdx && amino_acids[cut_idx] != LysIdx )
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],cut_prop);
	rankFeatureIdx+=num_aas;

	// features of the form L|XR   L|XXR   L|XXXR   L|XXXXR

	if (cut_idx>0 &&cut_idx < length-1 &&
		amino_acids[cut_idx+1]==ArgIdx &&
		amino_acids[cut_idx-1] != LysIdx && amino_acids[cut_idx-1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>0 &&cut_idx < length-2 &&
		amino_acids[cut_idx+2]==ArgIdx &&
		amino_acids[cut_idx-1] != LysIdx && amino_acids[cut_idx-1] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx && amino_acids[cut_idx+1] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>0 &&cut_idx < length-3 &&
		amino_acids[cut_idx+3]==ArgIdx &&
		amino_acids[cut_idx-1] != LysIdx && amino_acids[cut_idx-1] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx] != ArgIdx)
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>0 && cut_idx < length-4 &&
		amino_acids[cut_idx+4]==ArgIdx &&
		amino_acids[cut_idx-1] != LysIdx && amino_acids[cut_idx-1] != ArgIdx &&
		amino_acids[cut_idx+1] != ArgIdx && amino_acids[cut_idx+2] != ArgIdx &&
		amino_acids[cut_idx+3] != ArgIdx && amino_acids[cut_idx] != ArgIdx )
			sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;
	

	// features of the form KXF KXXF KXXXF KXXXXF
	if (cut_idx>2 && amino_acids[cut_idx-3] == LysIdx &&amino_acids[cut_idx-2] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>3 && amino_acids[cut_idx-4] == LysIdx && amino_acids[cut_idx-3] != ArgIdx &&
		amino_acids[cut_idx-2] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>4 && amino_acids[cut_idx-5] == LysIdx && amino_acids[cut_idx-4] != ArgIdx &&
		amino_acids[cut_idx-3] != ArgIdx && amino_acids[cut_idx-2] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>5 && amino_acids[cut_idx-6] == LysIdx && amino_acids[cut_idx-5] != ArgIdx &&
		amino_acids[cut_idx-4] != ArgIdx && amino_acids[cut_idx-3] != ArgIdx &&
		amino_acids[cut_idx-2] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	// features of the form KXFX KXXFX KXXXFX KXXXXFX
	if (cut_idx>3 && amino_acids[cut_idx-4] == LysIdx && amino_acids[cut_idx-3] != ArgIdx &&
		amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>4 && amino_acids[cut_idx-5] == LysIdx && amino_acids[cut_idx-4] != ArgIdx &&
		amino_acids[cut_idx-3] != ArgIdx && amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>5 && amino_acids[cut_idx-6] == LysIdx && amino_acids[cut_idx-5] != ArgIdx &&
		amino_acids[cut_idx-4] != ArgIdx && amino_acids[cut_idx-3] != ArgIdx && 
		amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>6 && amino_acids[cut_idx-7] == LysIdx && amino_acids[cut_idx-6] != ArgIdx &&
		amino_acids[cut_idx-5] != ArgIdx && amino_acids[cut_idx-4] != ArgIdx && 
		amino_acids[cut_idx-3] != ArgIdx && amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],cut_prop);
	rankFeatureIdx+=num_aas;

	// features of the form KX|F KXX|F KXXX|F KXXXX|F
	if (cut_idx>1 && cut_idx<length && amino_acids[cut_idx-2] == LysIdx && 
		amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>2 && cut_idx<length && amino_acids[cut_idx-3] == LysIdx && 
		amino_acids[cut_idx-2] != ArgIdx && amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>3 && cut_idx<length && amino_acids[cut_idx-4] == LysIdx && 
		amino_acids[cut_idx-3] != ArgIdx && amino_acids[cut_idx-2] != ArgIdx && 
		amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>4 && cut_idx<length && amino_acids[cut_idx-5] == LysIdx && 
		amino_acids[cut_idx-4] != ArgIdx && amino_acids[cut_idx-3] != ArgIdx && 
		amino_acids[cut_idx-2] != ArgIdx && amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

		// features of the form RXF RXXF RXXXF RXXXXF
	if (cut_idx>2 && amino_acids[cut_idx-3] == ArgIdx &&amino_acids[cut_idx-2] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>3 && amino_acids[cut_idx-4] == ArgIdx && amino_acids[cut_idx-3] != ArgIdx &&
		amino_acids[cut_idx-2] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>4 && amino_acids[cut_idx-5] == ArgIdx && amino_acids[cut_idx-4] != ArgIdx &&
		amino_acids[cut_idx-3] != ArgIdx && amino_acids[cut_idx-2] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>5 && amino_acids[cut_idx-6] == ArgIdx && amino_acids[cut_idx-5] != ArgIdx &&
		amino_acids[cut_idx-4] != ArgIdx && amino_acids[cut_idx-3] != ArgIdx &&
		amino_acids[cut_idx-2] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;

	// features of the form RXFX RXXFX RXXXFX RXXXXFX
	if (cut_idx>3 && amino_acids[cut_idx-4] == ArgIdx && amino_acids[cut_idx-3] != ArgIdx &&
		amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>4 && amino_acids[cut_idx-5] == ArgIdx && amino_acids[cut_idx-4] != ArgIdx &&
		amino_acids[cut_idx-3] != ArgIdx && amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>5 && amino_acids[cut_idx-6] == ArgIdx && amino_acids[cut_idx-5] != ArgIdx &&
		amino_acids[cut_idx-4] != ArgIdx && amino_acids[cut_idx-3] != ArgIdx && 
		amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>6 && amino_acids[cut_idx-7] == ArgIdx && amino_acids[cut_idx-6] != ArgIdx &&
		amino_acids[cut_idx-5] != ArgIdx && amino_acids[cut_idx-4] != ArgIdx && 
		amino_acids[cut_idx-3] != ArgIdx && amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],cut_prop);
	rankFeatureIdx+=num_aas;

	// features of the form RX|F RXX|F RXXX|F RXXXX|F
	if (cut_idx>1 && cut_idx<length && amino_acids[cut_idx-2] == ArgIdx && 
		amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>2 && cut_idx<length && amino_acids[cut_idx-3] == ArgIdx && 
		amino_acids[cut_idx-2] != ArgIdx && amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>3 && cut_idx<length && amino_acids[cut_idx-4] == ArgIdx && 
		amino_acids[cut_idx-3] != ArgIdx && amino_acids[cut_idx-2] != ArgIdx && 
		amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;

	if (cut_idx>4 && cut_idx<length && amino_acids[cut_idx-5] == ArgIdx && 
		amino_acids[cut_idx-4] != ArgIdx && amino_acids[cut_idx-3] != ArgIdx && 
		amino_acids[cut_idx-2] != ArgIdx && amino_acids[cut_idx-1] != ArgIdx)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;


	// Add sepcial C-terminal features
	const int c_cut_dis = length - cut_idx;
	if (c_cut_dis<=5 && amino_acids[length-1]==LysIdx)
	{
		if (c_aa_counts[LeuIdx]>0 &&
			(c_aa_counts[AspIdx]+c_aa_counts[GluIdx]>0 &&
			c_aa_counts[LeuIdx] + c_aa_counts[GluIdx] + +c_aa_counts[AspIdx] +
			c_aa_counts[LysIdx] == c_cut_dis) )
			sample.add_real_feature(rankFeatureIdx,c_cut_dis);
		rankFeatureIdx++;

		if (c_aa_counts[LeuIdx]>0 && c_aa_counts[AlaIdx]>0 &&
			c_aa_counts[LeuIdx] + c_aa_counts[AlaIdx] + c_aa_counts[LysIdx] == c_cut_dis)
			sample.add_real_feature(rankFeatureIdx,c_cut_dis);
		rankFeatureIdx++;

		if (c_cut_dis>=3 && c_aa_counts[LeuIdx] + c_aa_counts[LysIdx] + c_aa_counts[GluIdx] + 
			c_aa_counts[AspIdx] == c_cut_dis-1)
			sample.add_real_feature(rankFeatureIdx,c_cut_dis);
		rankFeatureIdx++;
	}
	else
		rankFeatureIdx+=3;

	if (amino_acids[length-1]==ArgIdx && cut_idx<length-1)
	{
		if (c_aa_counts[LeuIdx] + c_aa_counts[ArgIdx] + c_aa_counts[GluIdx] + 
			c_aa_counts[AspIdx] == c_cut_dis)
		sample.add_real_feature(rankFeatureIdx,c_cut_dis);
		
	}
	rankFeatureIdx++;
}



void push_back_all_RHK_pairs(vector<string>& real_feature_names, string prefix_label)
{
	int n;
	for (n=0; n<num_RKH_combos; n++)
	{
		int c;
		for (c=0; c<num_RKH_combos; c++)
		{
			string label = prefix_label + " " + combos[n] + "/" + combos[c];
			real_feature_names.push_back(label);
		}
	}
}


void PeakRankModel::set_advanced_feature_names()
{
	const int num_aas = this->get_num_model_aas();
	int i;

	real_feature_names.clear();
	real_feature_stage_idxs.clear();

	real_feature_names.push_back("N RKH LEVEL");
	real_feature_names.push_back("C RKH LEVEL");

	real_feature_names.push_back("PROP PEAK MASS IN VIS RANGE"); 
	push_back_all_RHK_pairs(real_feature_names,"PROP PEAK MASS IN VIS RANGE");

	real_feature_names.push_back("DIS FROM MIN");
	push_back_all_RHK_pairs(real_feature_names,"DIS FROM MIN");

	
	real_feature_names.push_back("DIS FROM MAX");
	push_back_all_RHK_pairs(real_feature_names,"DIS FROM MAX");

	real_feature_names.push_back("CUT POS");
	push_back_all_RHK_pairs(real_feature_names,"CUT POS");

	real_feature_names.push_back("DIS MIN<75, CUT POS");
	real_feature_names.push_back("DIS MIN<150, CUT POS");
	real_feature_names.push_back("DIS MIN<250, CUT POS");
	real_feature_names.push_back("DIS MIN<=400, CUT POS");

	real_feature_names.push_back("DIS MAX<75, CUT POS");
	real_feature_names.push_back("DIS MAX<150, CUT POS");
	real_feature_names.push_back("DIS MAX<250, CUT POS");
	real_feature_names.push_back("DIS MAX<=400, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("#N-SIDE " + model_aa_labels[i] + "(without flanking 3)");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("#C-SIDE " + model_aa_labels[i] + "(without flanking 3)");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("#N-SIDE " + model_aa_labels[i]);

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("#C-SIDE " + model_aa_labels[i]);

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT -1 = " + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT -2 = " + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT -3 = " + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT +1 = " + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT +2 = " + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT +3 = " + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
	{
		int j;
		for (j=0; j<num_aas; j++)
			real_feature_names.push_back("CUT IS " + model_aa_labels[i] + "-" + model_aa_labels[j] + ", CUT POS");
	}

	for (i=0; i<num_aas; i++)
	{
		int j;
		for (j=0; j<num_aas; j++)
			real_feature_names.push_back("CUT IS " + model_aa_labels[i] + "-" + model_aa_labels[j] + ", PEAK PROP");
	}


	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT -1 = " + model_aa_labels[i] + ", RHK LEVEL");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT -2 = " + model_aa_labels[i] + ", RHK LEVEL");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT -3 = " + model_aa_labels[i] + ", RHK LEVEL");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT +1 = " + model_aa_labels[i] + ", RHK LEVEL");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT +2 = " + model_aa_labels[i] + ", RHK LEVEL");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT +3 = " + model_aa_labels[i] + ", RHK LEVEL");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT -1 = " + model_aa_labels[i] + ", PEAK PROP");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT -2 = " + model_aa_labels[i] + ", PEAK PROP");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT -3 = " + model_aa_labels[i] + ", PEAK PROP");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT +1 = " + model_aa_labels[i] + ", PEAK PROP");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT +2 = " + model_aa_labels[i] + ", PEAK PROP");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT +3 = " + model_aa_labels[i] + ", PEAK PROP");

	for (i=0; i<num_aas; i++)
	{
		int j;
		for (j=0; j<num_aas; j++)
			real_feature_names.push_back(model_aa_labels[i] + model_aa_labels[j] + " BEFORE CUT, CUT POS");
	}

	for (i=0; i<num_aas; i++)
	{
		int j;
		for (j=0; j<num_aas; j++)
			real_feature_names.push_back(model_aa_labels[i] + model_aa_labels[j] + " AFTER CUT, CUT POS");
	}


	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|"+ model_aa_labels[i] + "XK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|"+ model_aa_labels[i] + "XXK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|"+ model_aa_labels[i] + "XXXK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|"+ model_aa_labels[i] + "XXXXK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|X"+ model_aa_labels[i] + "XK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|X"+ model_aa_labels[i] + "XXK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|X"+ model_aa_labels[i] + "XXXK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|X"+ model_aa_labels[i] + "XXXXK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( model_aa_labels[i] + "|XK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( model_aa_labels[i] + "|XXK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( model_aa_labels[i] + "|XXXK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( model_aa_labels[i] + "|XXXXK, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|"+ model_aa_labels[i] + "XR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|"+ model_aa_labels[i] + "XXR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|"+ model_aa_labels[i] + "XXXR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|"+ model_aa_labels[i] + "XXXXR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|X"+ model_aa_labels[i] + "XR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|X"+ model_aa_labels[i] + "XXR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|X"+ model_aa_labels[i] + "XXXR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "|X"+ model_aa_labels[i] + "XXXXR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( model_aa_labels[i] + "|XR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( model_aa_labels[i] + "|XXR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( model_aa_labels[i] + "|XXXR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( model_aa_labels[i] + "|XXXXR, CUT POS");


	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KX" + model_aa_labels[i] + "|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KXX" + model_aa_labels[i] + "|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KXXX" + model_aa_labels[i] + "|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KXXXX" + model_aa_labels[i] + "|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KX" + model_aa_labels[i] + "X|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KXX" + model_aa_labels[i] + "X|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KXXX" + model_aa_labels[i] + "X|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KXXXX" + model_aa_labels[i] + "X|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KX|" + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KXX|" + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KXXX|" + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "KXXXX|" + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RX" + model_aa_labels[i] + "|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RXX" + model_aa_labels[i] + "|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RXXX" + model_aa_labels[i] + "|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RXXXX" + model_aa_labels[i] + "|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RX" + model_aa_labels[i] + "X|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RXX" + model_aa_labels[i] + "X|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RXXX" + model_aa_labels[i] + "X|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RXXXX" + model_aa_labels[i] + "X|, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RX|" + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RXX|" + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RXXX|" + model_aa_labels[i] + ", CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( "RXXXX|" + model_aa_labels[i] + ", CUT POS");

	real_feature_names.push_back("|[LD/EK] - permuted");
	real_feature_names.push_back("|[LAK] - permuted");
	real_feature_names.push_back("|[LDEK] - permuted with one X");
	real_feature_names.push_back("|[LDER] - permuted");
}




