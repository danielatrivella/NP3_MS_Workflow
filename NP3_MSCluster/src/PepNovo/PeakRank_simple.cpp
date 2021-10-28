#include "PeakRankModel.h"

extern const int num_RKH_combos ;
extern const int num_RKH_pairs ;
extern const string combos[]; 


extern const float RKH_pair_matrix[6][6];

/**********************************************************************
This fills in all possible sequence features that can be used in any
of the models (these features are good for all fragments/mobility/
size/charge). The models can later choose to ignore some of these
features by asigning the weight 0. All features are derived directly
from the peptide sequence and concern the specified cut idx.
***********************************************************************/
void PeakRankModel::fill_simple_peak_features(
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
	const float rounded_peak_prop = 0.05*floor(peak_mass_prop * 20.0);

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
//	if (length>=12) side_length=4;
//	if (length>=15) side_length=5;

	float cut_prop;
	if (cut_idx<=side_length)
	{
		cut_prop=(float)cut_idx-4;
	}
	else if (cut_idx>=length-side_length)
	{
		cut_prop=(float)(4+cut_idx-length);
	}
	else
	{
		cut_prop = 0;
	}
	

	// peak prop
	sample.add_real_feature(rankFeatureIdx++,rounded_peak_prop);
	
/*	// fill dis features
	if (dis_from_min<dis_from_max)
	{
		sample.add_real_feature(rankFeatureIdx++,dis_from_min);
		rankFeatureIdx++;
	}
	else
	{
		rankFeatureIdx++;
		sample.add_real_feature(rankFeatureIdx++,dis_from_max);
	} 

	//  fill prop features
	sample.add_real_feature(rankFeatureIdx++,cut_prop);*/
	


	// fill aa count features 
	vector<int> n_aa_counts, c_aa_counts;
	n_aa_counts.resize(num_aas+1,0);
	c_aa_counts.resize(num_aas+1,0);

	for (i=0; i<cut_idx-2; i++)
		n_aa_counts[amino_acids[i]]++;

	for (i=cut_idx+2; i<length; i++)
		c_aa_counts[amino_acids[i]]++;

	int a;
	for (a=0; a<num_aas; a++)
		sample.add_real_feature(rankFeatureIdx++,n_aa_counts[a]);
	
	for (a=0; a<num_aas; a++)
		sample.add_real_feature(rankFeatureIdx++,c_aa_counts[a]);


/*	
	// fill aa flanking features N side
	if (cut_idx>0)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],cut_prop);
	rankFeatureIdx+=num_aas;
	if (cut_idx>1)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],cut_prop);
	rankFeatureIdx+=num_aas;

	// fill aa flanking features C side
	if (cut_idx<length)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],cut_prop);
	rankFeatureIdx+=num_aas;
	if (cut_idx<length-1)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],cut_prop);
	rankFeatureIdx+=num_aas;
	*/

	// fill aa flanking features N side
	if (cut_idx>0)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-1],1.0);
	rankFeatureIdx+=num_aas;
	if (cut_idx>1)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx-2],1.0);
	rankFeatureIdx+=num_aas;

	// fill aa flanking features C side
	if (cut_idx<length)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx],1.0);
	rankFeatureIdx+=num_aas;
	if (cut_idx<length-1)
		sample.add_real_feature(rankFeatureIdx+amino_acids[cut_idx+1],1.0);
	rankFeatureIdx+=num_aas;
}




void PeakRankModel::set_simple_feature_names()
{
	const int num_aas = this->get_num_model_aas();
	int i;

	real_feature_names.clear();
	real_feature_stage_idxs.clear();


	real_feature_names.push_back("PROP PEAK MASS IN VIS RANGE"); 
//	real_feature_names.push_back("DIS FROM MIN");
//	real_feature_names.push_back("DIS FROM MAX");

//	real_feature_names.push_back("CUT POS");


	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("#N-SIDE " + model_aa_labels[i] + " (wihtout closest 2)");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("#C-SIDE " + model_aa_labels[i] + " (wihtout closest 2)");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT -1 = " + model_aa_labels[i]);

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT -2 = " + model_aa_labels[i]);


	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT +1 = " + model_aa_labels[i]);

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("CUT +2 = " + model_aa_labels[i]);

}




