#include "PeakRankModel.h"
#include "PepNovo_auxfun.h"

extern const int num_RKH_combos;
extern const int num_RKH_pairs;
extern const string combos[];
extern const float RKH_pair_matrix[6][6];

/**********************************************************************
This fills in all possible sequence features that can be used in any
of the models (these features are good for all fragments/mobility/
size/charge). The models can later choose to ignore some of these
features by asigning the weight 0. All features are derived directly
from the peptide sequence and concern the specified cut idx.
***********************************************************************/
void PeakRankModel::fill_partial_denovo_peak_features(
								 mass_t n_mass,
								 mass_t c_mass,
								 const  vector<int>& org_amino_acids,
								 int    cut_idx,
								 mass_t cut_mass,
								 mass_t pm_with_19,
								 int	spec_charge,
								 const FragmentType& fragment,
								 int most_basic_on_n_side, // if the n side does not reach the terminal
								 int most_basic_on_c_side, // if the c side does not reach the terminal
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

	// uses regular amino acid codes
	if (most_basic_on_n_side>0)
	{
		if (most_basic_on_n_side == His) 
			num_nH++;
		if (most_basic_on_n_side == Lys) 
			num_nK++;
		if (most_basic_on_n_side == Arg) 
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

	// uses regular amino acid codes
	if (most_basic_on_c_side>0)
	{
		if (most_basic_on_c_side == His) 
			num_cH++;
		if (most_basic_on_c_side == Lys) 
			num_cK++;
		if (most_basic_on_c_side == Arg) 
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
	if (pm_with_19>=1000) side_length=4;
	if (pm_with_19>=1500) side_length=5;

	float cut_prop;
	if (cut_idx<=side_length && n_mass<1.0)
	{
		cut_prop=(float)cut_idx;
	}
	else if (cut_idx>=length-side_length && c_mass > pm_with_19 - 30.0)
	{
		cut_prop=(float)(15+cut_idx-length);
	}
	else
	{
		cut_prop = 6.1+floor(3.0*(cut_mass/pm_with_19))*0.1;
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
		sample.add_real_feature(rankFeatureIdx++,dis_from_min);
		rankFeatureIdx++;
	}
	else
	{
		rankFeatureIdx++;
		sample.add_real_feature(rankFeatureIdx++,dis_from_max);
	} 

	//  fill prop features
	sample.add_real_feature(rankFeatureIdx++,cut_prop);




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
}



void push_back_all_RHK_pairs_dnv(vector<string>& real_feature_names,
							 string prefix_label)
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



void PeakRankModel::set_partial_denovo_feature_names()
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
	
	real_feature_names.push_back("DIS FROM MAX");

	real_feature_names.push_back("CUT POS");

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
		real_feature_names.push_back( model_aa_labels[i] + "|XR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( model_aa_labels[i] + "|XXR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( model_aa_labels[i] + "|XXXR, CUT POS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back( model_aa_labels[i] + "|XXXXR, CUT POS");
}




/****************************************************************************
Reduces a TrainingPeptide in size to form a denovo training peptide.
The shorter peptide might have intensities in cuts 0, and cut |length|.
*****************************************************************************/
void convert_tps_to_partial_denovo(Config *config, 
								   vector<TrainingPeptide>& all_tps, 
								   int num_to_add,
								   bool verbose)
{
	const int y_frag_idx = config->get_frag_idx_from_label("y");
	const int b_frag_idx = config->get_frag_idx_from_label("b");
	const int y2_frag_idx = config->get_frag_idx_from_label("s2+10.2");
	const int b2_frag_idx = config->get_frag_idx_from_label("b2");
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const int org_num_tps = all_tps.size();
	vector<double> rand_threshes;
	
	rand_threshes.resize(100,0);
	int i;
	for (i=0; i<=6; i++)
		rand_threshes[i]=1.0;
	rand_threshes[7] = 0.7;    rand_threshes[8] = 0.65;  rand_threshes[9]=0.6;
	rand_threshes[10]= 0.6;    rand_threshes[11]= 0.55; rand_threshes[12]=0.5;
	rand_threshes[13]= 0.45;   rand_threshes[14]=0.4;   rand_threshes[15]=0.3;
	rand_threshes[16]= 0.2;    rand_threshes[17]=0.15;  rand_threshes[18]=0.1;
	rand_threshes[19]= 0.075;  rand_threshes[20]=0.05;

	if (num_to_add>0)
		all_tps.resize(org_num_tps*(1+num_to_add));

	int tp_idx;
	for (tp_idx=0; tp_idx<org_num_tps; tp_idx++)
	{
		const TrainingPeptide org_tp = all_tps[tp_idx];
		const int org_length = org_tp.length;

		if (verbose)
		{
			cout << tp_idx << "\t";
			org_tp.print(config);
			cout << endl;
		}

		int i;
		for (i=0; i<= num_to_add; i++)
		{
			TrainingPeptide& new_tp = all_tps[i*org_num_tps + tp_idx];

			if (rand_threshes[org_length]>0 && myRandom()<rand_threshes[org_length])
			{
				new_tp = org_tp;
				continue;
			}

			// create a ruff score vector to use for subset selection
			vector<int> ruff_scores;
			ruff_scores.resize(org_tp.length+1,0);
			ruff_scores[0]=2;
			ruff_scores[org_tp.length]=1;
			
			const int last_aa = org_tp.amino_acids[org_tp.length-1];
			if (last_aa==Arg || last_aa==Lys)
				ruff_scores[org_tp.length]=2;
			
			vector<int> frags;
			frags.push_back(y_frag_idx);
			frags.push_back(b_frag_idx);
			if (org_tp.charge == 3 || (org_tp.charge == 2 && org_tp.pm_with_19 > 1800))
			{
				frags.push_back(y2_frag_idx);
				frags.push_back(b2_frag_idx);
			}

			int f,c;
			for (f=0; f<frags.size(); f++)
			{
				const int pos = org_tp.get_frag_idx_pos(frags[f]);
				if (pos<0)
					continue;

				int c;
				for (c=1; c<org_tp.length; c++)
					if (org_tp.intens[pos][c]>0)
						ruff_scores[c]++;
			}
			for (c=1; c<ruff_scores.size(); c++)
				ruff_scores[c]+=ruff_scores[c-1];


			int min_length = 6;
			int max_length = org_tp.length - 3;
			if (org_tp.length == 7)
				max_length = 6;
			if (org_tp.length == 8)
				max_length = 7;
			if (org_tp.length == 9)
				max_length = 7;

			if (max_length>20)
				max_length=20;

			const int selected_length = min_length + (int)(myRandom()*(max_length-min_length+1));
			
			// find best subset of this length
			const int max_start = org_tp.length - selected_length;
			int best_idx=-1;
			int best_ruff_score=0;
			

			for (c=0; c<=max_start; c++)
			{
				int score = ruff_scores[c+selected_length];
				if (c>0)
					score -= ruff_scores[c-1];
				if (score>=best_ruff_score)
				{
					best_idx=c;
					best_ruff_score = score;
				}
				if (verbose)
				{
					cout << c << ":" << score << " ";
				}
			}
			if (verbose)
				cout << endl;

			if (best_idx<0)
			{
				cout << "Error: something is wrong with the ruff scores!" << endl;
				exit(1);
			}

			// create new tp
			const vector<int>& amino_acids = org_tp.amino_acids;
			int j;
			new_tp.amino_acids.resize(selected_length,0);
			mass_t n_mass =0;
			for (j=0; j<best_idx; j++)
				n_mass+=aa2mass[amino_acids[j]];

			new_tp.n_mass = n_mass;
			mass_t total_mass = n_mass;
			for (j=0; j<selected_length; j++)
			{
				new_tp.amino_acids[j]=amino_acids[j+best_idx];
				total_mass += aa2mass[new_tp.amino_acids[j]];
			}

			if (total_mass>org_tp.pm_with_19)
			{
				cout << "Error: mismatch in masses with de novo tp!" << endl;
				exit(1);
			}



			if (best_idx>0)
			{
				int j;
				for (j=0; j<best_idx; j++)
					if (amino_acids[j] == His)
						new_tp.best_n_removed = His;
				for (j=0; j<best_idx; j++)
					if (amino_acids[j] == Lys)
						new_tp.best_n_removed = Lys;
				for (j=0; j<best_idx; j++)
					if (amino_acids[j] == Arg)
						new_tp.best_n_removed = Arg;
			}

			const int c_idx = best_idx + selected_length;
			if (c_idx<org_tp.length)
			{
				int j;
				for (j=c_idx; j<org_tp.length; j++)
					if (amino_acids[j] == His)
						new_tp.best_c_removed = His;
				for (j=c_idx; j<org_tp.length; j++)
					if (amino_acids[j] == Lys)
						new_tp.best_c_removed = Lys;
				for (j=c_idx; j<org_tp.length; j++)
					if (amino_acids[j] == Arg)
						new_tp.best_c_removed = Arg;
			}

			new_tp.length = selected_length;
			new_tp.charge = org_tp.charge;
			new_tp.frag_idxs = org_tp.frag_idxs;
			
			// change mobility according to observed + removed aas
			int num_arg=0,num_lys=0,num_his=0;
			for (j=0; j<amino_acids.size(); j++)
			{
				if (amino_acids[j]==Arg)
					num_arg++;
				if (amino_acids[j]==Lys)
					num_lys++;
				if (amino_acids[j]==His)
					num_his++;
			}

			if (new_tp.best_c_removed>0)
			{
				if (new_tp.best_c_removed == Arg) num_arg++;
				if (new_tp.best_c_removed == Lys) num_lys++;
				if (new_tp.best_c_removed == His) num_his++;
			}

			if (new_tp.best_n_removed>0)
			{
				if (new_tp.best_n_removed == Arg) num_arg++;
				if (new_tp.best_n_removed == Lys) num_lys++;
				if (new_tp.best_n_removed == His) num_his++;
			}

			new_tp.mobility = get_proton_mobility(org_tp.charge,num_arg,num_his,num_lys);
			
			new_tp.pm_with_19 = org_tp.pm_with_19;

			// reduce the sets of intensities
			new_tp.intens.resize(new_tp.frag_idxs.size());
			for (f=0; f<new_tp.frag_idxs.size(); f++)
			{
				new_tp.intens[f].clear();
				if (org_tp.intens[f].size()>0)
				{
					new_tp.intens[f].resize(selected_length+1,NEG_INF);
					int j;
					for (j=0; j<selected_length; j++)
					{
						new_tp.intens[f][j]=org_tp.intens[f][j+best_idx];
					}

				}
			}

			if (verbose)
			{
				cout << "\t";
				new_tp.print(config);
				cout << endl;
			}
		}
		if (verbose)
			cout << endl;
	}
}
