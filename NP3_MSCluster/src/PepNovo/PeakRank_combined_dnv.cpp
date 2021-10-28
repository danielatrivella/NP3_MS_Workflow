#include "PeakRankModel.h"
#include "PepNovo_auxfun.h"

extern const int num_RKH_combos;
extern const int num_RKH_pairs;
extern const string combos[];
extern const float RKH_pair_matrix[6][6];

void PartitionModel::fill_combined_dnv_peak_features(
								 const PeakRankModel *prank,
								 const mass_t n_mass,   // this is where the possibly partial peptide starts
								 const mass_t c_mass,
								 const  vector<int>& org_amino_acids,
								 const int    cut_idx,
								 const mass_t cut_mass,
								 const PeptideSolution& sol,
								 const FragmentType& frag,
								 const int position_idx_in_model_fragment_type_idxs,
								 RankBoostSample& sample) const
{
	const int most_basic_on_n_side = sol.most_basic_aa_removed_from_n;
	const int most_basic_on_c_side = sol.most_basic_aa_removed_from_c;
	const bool reaches_n_terminal = sol.reaches_n_terminal;
	const bool reaches_c_terminal = sol.reaches_c_terminal;
	const mass_t pm_with_19 = sol.pm_with_19;
	const int spec_charge=sol.charge;
	const vector<string>& model_aa_labels = prank->get_model_aa_labels();
	const vector<int>& session_aas_to_model_aas = prank->get_session_aas_to_model_aas();
	const int length = org_amino_acids.size();
	const int num_aas = model_aa_labels.size();
	int rankFeatureIdx=0;
	int i;

	vector<int> amino_acids;
	prank->convert_aas_to_model_aas(org_amino_acids, amino_acids);

	if (amino_acids.size() != org_amino_acids.size())
	{
		cout << "Error: aa size mismatch!" << endl;
		exit(1);
	}


	bool bad_cut = false;
	if ((reaches_n_terminal && cut_idx<=0) || (! reaches_n_terminal && cut_idx<0))
		bad_cut=true;

	if ((reaches_c_terminal && cut_idx>=amino_acids.size()) || (! reaches_c_terminal && cut_idx>amino_acids.size()))
		bad_cut = true;

	if (bad_cut)
	{
		cout << "Error: cut_idx is bad ( " << cut_idx << " , length = " << amino_acids.size() << " )" << endl;
		cout << "c-mass " << fixed << c_mass + MASS_OHHH << " , PM+19=" << pm_with_19 << endl;
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
	const mass_t max_detected_mass = prank->get_max_detected_mass();
	const mass_t exp_peak_mass = frag.calc_expected_mass(cut_mass,pm_with_19);
	const mass_t min_obs_mass = prank->calc_min_detected_mass(pm_with_19,spec_charge);
	const mass_t max_obs_mass = (pm_with_19>max_detected_mass ? max_detected_mass : pm_with_19);
	
	const float peak_mass_prop = ((exp_peak_mass - min_obs_mass)/(max_obs_mass - min_obs_mass));
	const float rounded_peak_prop = 0.02*floor(peak_mass_prop * 50.0);

	const mass_t dis_from_min = 20.0*floor((exp_peak_mass - min_obs_mass)*0.05);
	const mass_t dis_from_max = 20.0*floor((max_obs_mass  - exp_peak_mass)*0.05);
	const mass_t dis_from_minmax = (dis_from_min<dis_from_max ? dis_from_min : prank->get_max_detected_mass() - dis_from_max);

	const int RKH_n_combo_idx = calc_RKH_combo_idx(num_nR,num_nK,num_nH);
	const int RKH_c_combo_idx = calc_RKH_combo_idx(num_cR,num_cK,num_cH);
	const int RKH_pair_idx = (RKH_n_combo_idx * num_RKH_combos) + RKH_c_combo_idx;
	const float RKH_liniar_pair_idx = RKH_pair_matrix[RKH_n_combo_idx][RKH_c_combo_idx];

	const float cut_prop = 0.02 * floor((cut_mass / pm_with_19)*50);
	const float cut_n_mass = 20 * floor(cut_mass*0.05);
	const float cut_c_mass = 20 * floor((pm_with_19-cut_mass)*0.05);
	const int cut_dis_from_n = cut_idx;
	const int cut_dis_from_c = length-cut_idx;


	const int third=length/3;
	const int two_thirds = length - third;
	float cut_pos=NEG_INF;
	if (reaches_n_terminal && (cut_idx<=4 || (cut_idx<=third && cut_idx<8)))
	{
		cut_pos = (float)cut_idx;
	}
	else if (reaches_c_terminal && (cut_idx>=length-4 || (cut_idx>=two_thirds && cut_idx>length-8)))
	{
		cut_pos = (float)(20-length+cut_idx);
	}
	else
		cut_pos = 10.0 + cut_prop;

	const int n_aa = (cut_idx > 0 ?  amino_acids[cut_idx-1] : -1);
	const int c_aa = (cut_idx < length ? amino_acids[cut_idx] : -1);
	const int n_term_aa = amino_acids[0];
	const int c_term_aa = amino_acids[length-1];

	sample.clear();

	// first feature is always 0, tells what fragment this is
	sample.add_real_feature(0,(float)position_idx_in_model_fragment_type_idxs);

	// this is the new starting index for features in this sample. The offset
	// is computaed according to the fragment type
	int f_idx =1 + position_idx_in_model_fragment_type_idxs * num_features_per_frag; // num_features_per_frag_type

	if (num_features_per_frag==0)
	{
		cout << "Error: num_features_per_frag is 0!" << endl;
		exit(1);
	}


	// add general position features
	sample.add_real_feature(f_idx++,dis_from_minmax);
	sample.add_real_feature(f_idx++,rounded_peak_prop);
	sample.add_real_feature(f_idx++,RKH_pair_idx);
	sample.add_real_feature(f_idx++,RKH_liniar_pair_idx);
	sample.add_real_feature(f_idx++,cut_prop);
	sample.add_real_feature(f_idx++,cut_pos);
	sample.add_real_feature(f_idx++,cut_n_mass);
	sample.add_real_feature(f_idx++,cut_c_mass);

	if (reaches_n_terminal)
		sample.add_real_feature(f_idx,(float)cut_dis_from_n);
	f_idx++;

	if (reaches_c_terminal)
		sample.add_real_feature(f_idx,(float)cut_dis_from_c);
	f_idx++;

//	cout << reaches_n_terminal << reaches_c_terminal << endl;
	
	const float position_var = rounded_peak_prop; // this will be used for all amino acid features

	// fill aa flanking features N side
	if (cut_idx>0)
		sample.add_real_feature(f_idx+amino_acids[cut_idx-1],position_var);
	f_idx+=num_aas;
	if (cut_idx>1)
		sample.add_real_feature(f_idx+amino_acids[cut_idx-2],position_var);
	f_idx+=num_aas;
	if (cut_idx>2)
		sample.add_real_feature(f_idx+amino_acids[cut_idx-3],position_var);
	f_idx+=num_aas;

	// fill aa flanking features C side
	if (cut_idx<length)
		sample.add_real_feature(f_idx+amino_acids[cut_idx],position_var);
	f_idx+=num_aas;
	if (cut_idx<length-1)
		sample.add_real_feature(f_idx+amino_acids[cut_idx+1],position_var);
	f_idx+=num_aas;
	if (cut_idx<length-2)
		sample.add_real_feature(f_idx+amino_acids[cut_idx+2],position_var);
	f_idx+=num_aas;

	// fill cut pair features X-Y
	if (cut_idx>0 && cut_idx<length)
		sample.add_real_feature(f_idx+(n_aa*num_aas+c_aa),position_var);
	f_idx+=(num_aas*num_aas);

	// fill aa count features (up top 3 aa's away from cut)
	if (reaches_n_terminal)
	{
		vector<int> n_aa_counts;
		n_aa_counts.resize(num_aas+1,0);
		
		int i;
		for (i=0; i<cut_idx-3; i++)
			n_aa_counts[amino_acids[i]]++;

		
		int a;
		for (a=0; a<num_aas; a++)
			if (n_aa_counts[a]>0)
				sample.add_real_feature(f_idx+a,n_aa_counts[a]);
	}
	f_idx+=num_aas;

	if (reaches_c_terminal)
	{
		vector<int> c_aa_counts;
		c_aa_counts.resize(num_aas+1,0);

		int i;
		for (i=cut_idx+3; i<length; i++)
			c_aa_counts[amino_acids[i]]++;
				
		int a;
		for (a=0; a<num_aas; a++)
			if (c_aa_counts[a]>0)
				sample.add_real_feature(f_idx+a,c_aa_counts[a]);
	}
	f_idx+=num_aas;

	// add Cut\terminal features
	if (reaches_n_terminal)
		sample.add_real_feature(f_idx+n_term_aa,cut_pos);
	f_idx+=num_aas;

	if (reaches_c_terminal)
		sample.add_real_feature(f_idx+c_term_aa,cut_pos);
	f_idx+=num_aas;

	// X|
	if (cut_idx>0)
	{
		sample.add_real_feature(f_idx+amino_acids[cut_idx-1],cut_pos);
		f_idx+=num_aas;
		if ((reaches_c_terminal && c_term_aa == LysIdx) || most_basic_on_c_side == Lys)
			sample.add_real_feature(f_idx+amino_acids[cut_idx-1],cut_pos);
		f_idx+=num_aas;
		if ((reaches_c_terminal && c_term_aa == ArgIdx) || most_basic_on_c_side == Arg)
			sample.add_real_feature(f_idx+amino_acids[cut_idx-1],cut_pos);
		f_idx+=num_aas;
	}
	else
		f_idx+=3*num_aas;

	// X_|
	if (cut_idx>1)
	{
		sample.add_real_feature(f_idx+amino_acids[cut_idx-2],cut_pos);
		f_idx+=num_aas;
		if ((reaches_c_terminal && c_term_aa == LysIdx) || most_basic_on_c_side == Lys)
			sample.add_real_feature(f_idx+amino_acids[cut_idx-2],cut_pos);
		f_idx+=num_aas;
		if ((reaches_c_terminal && c_term_aa == ArgIdx) || most_basic_on_c_side == Arg)
			sample.add_real_feature(f_idx+amino_acids[cut_idx-2],cut_pos);
		f_idx+=num_aas;
	}
	else
		f_idx+=3*num_aas;

	// |X
	if (cut_idx<length)
	{
		sample.add_real_feature(f_idx+amino_acids[cut_idx],cut_pos);
		f_idx+=num_aas;
		if ((reaches_c_terminal && c_term_aa == LysIdx) || most_basic_on_c_side == Lys)
			sample.add_real_feature(f_idx+amino_acids[cut_idx],cut_pos);
		f_idx+=num_aas;
		if ((reaches_c_terminal && c_term_aa == ArgIdx) || most_basic_on_c_side == Arg)
			sample.add_real_feature(f_idx+amino_acids[cut_idx],cut_pos);
		f_idx+=num_aas;
	}
	else
		f_idx+=3*num_aas;

	// |_X
	if (cut_idx<length-1)
	{
		sample.add_real_feature(f_idx+amino_acids[cut_idx+1],cut_pos);
		f_idx+=num_aas;
		if ((reaches_c_terminal && c_term_aa == LysIdx) || most_basic_on_c_side == Lys)
			sample.add_real_feature(f_idx+amino_acids[cut_idx+1],cut_pos);
		f_idx+=num_aas;
		if ((reaches_c_terminal && c_term_aa == ArgIdx) || most_basic_on_c_side == Arg)
			sample.add_real_feature(f_idx+amino_acids[cut_idx+1],cut_pos);
		f_idx+=num_aas;
	}
	else
		f_idx+=3*num_aas;

}


void PartitionModel::set_combined_dnv_feature_names_in_rankboost_model(const PeakRankModel *prank)
{
	const vector<string>& model_aa_labels = prank->get_model_aa_labels();
	const int num_aas = model_aa_labels.size();
	Config *config = prank->get_config();
	vector<string> frag_feature_names;
	vector<string> all_feature_names;

	frag_feature_names.clear();

	frag_feature_names.push_back("Dis Min/Max");
	frag_feature_names.push_back("Peak prop [Min-Max]");
	frag_feature_names.push_back("RHK pair idx");
	frag_feature_names.push_back("RHK liniar pair idx");
	frag_feature_names.push_back("Cut prop [0-M+19]");
	frag_feature_names.push_back("Cut pos");
	frag_feature_names.push_back("Cut N mass");
	frag_feature_names.push_back("Cut C mass");
	frag_feature_names.push_back("Cut idx from N");
	frag_feature_names.push_back("Cut idx from C");
//	push_back_all_RHK_pairs(frag_feature_names,"Cut prop");

	int i;
	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|_");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "_|_");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "__|_");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is _|" + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is _|_" + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is _|__" + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
	{
		int j;
		for (j=0; j<num_aas; j++)
			frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|" + model_aa_labels[j]);
	}

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("# N-side " + model_aa_labels[i]);

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("# C-side " + model_aa_labels[i]);

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("N-term aa is " + model_aa_labels[i] + ", cut pos" );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("C-term aa is " + model_aa_labels[i] + ", cut pos" );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|, cut pos");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|, cut pos, C-term is K");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|, cut pos, C-term is R");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "_|, cut pos");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "_|, cut pos, C-term is K");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "_|, cut pos, C-term is R");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |" + model_aa_labels[i] + ", cut pos");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |" + model_aa_labels[i] + ", cut pos, C-term is K");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |" + model_aa_labels[i] + ", cut pos, C-term is R");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |_" + model_aa_labels[i] + ", cut pos");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |_" + model_aa_labels[i] + ", cut pos, C-term is K");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |_" + model_aa_labels[i] + ", cut pos, C-term is R");

	// create all feature names
	all_feature_names.clear();
	all_feature_names.push_back("FRAG TYPE");
	int f;
	for (f=0; f<this->fragment_type_idxs.size(); f++)
	{
		const string frag_label = config->get_fragment(fragment_type_idxs[f]).label;
		int i;
		for (i=0; i<frag_feature_names.size(); i++)
		{
			string feature_name = frag_label + ": " + frag_feature_names[i];
			all_feature_names.push_back(feature_name);
		}
	}

	num_features_per_frag = frag_feature_names.size();
	vector<string> empty;
	empty.clear();
	combined_frag_boost_model.init_rankboost_model_feature_names(empty, all_feature_names);
}





