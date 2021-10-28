#include "RegularFragmentModel.h"

const char * ScoreModelFields_RI_names[]={
"RI_CONST",	"RI_LOG_LOCAL_RANK",	"RI_LOG_GLOBAL_RANK",	"RI_ISO_LEVEL",	"RI_IND_LOG_INTEN_LESS1",	"RI_LOG_INTEN_LESS1",	"RI_IND_LOG_INTEN_LESS2",	"RI_LOG_INTEN_LESS2",	"RI_IND_LOG_INTEN_LESS3",	"RI_LOG_INTEN_LESS3",	"RI_IND_LOG_INTEN_LESS4",	"RI_LOG_INTEN_LESS4",	"RI_IND_LOG_INTEN_MORE",	"RI_LOG_INTEN_MORE",	"RI_IND_DIS_FROM_MINMAX_LESS_50",	"RI_DIS_FROM_MINMAX0",	"RI_LOG_INTEN_DIS50",	
"RI_IND_DIS_FROM_MINMAX_LESS_150",	"RI_DIS_FROM_MINMAX50",	"RI_LOG_INTEN_DIS150",	"RI_IND_DIS_FROM_MINMAX_LESS_250",	"RI_DIS_FROM_MINMAX150",	"RI_LOG_INTEN_DIS250",	"RI_IND_DIS_FROM_MINMAX_MORE",	"RI_DIS_FROM_MINMAX250",	"RI_LOG_INTEN_DISMORE",	"RI_REL_POS0",	"RI_REL_POS1",	"RI_REL_POS2",	"RI_REL_POS3",	"RI_REL_POS4",	"RI_REL_POS5",	"RI_REL_POS6",	"RI_REL_POS7",	"RI_REL_POS8",	
"RI_REL_POS9",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_0",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_1",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_2",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_3",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_4",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_5",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_6",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_MORE6",	"RI_IND_PARENT_COMBO_0",	"RI_IND_PARENT_COMBO_1",	"RI_IND_PARENT_COMBO_2",	
"RI_IND_PARENT_COMBO_3",	"RI_IND_PARENT_COMBO_4",	"RI_IND_PARENT_COMBO_5",	"RI_IND_PARENT_COMBO_6",	"RI_IND_PARENT_COMBO_7",	"RI_IND_PARENT_COMBO_8",	"RI_IND_PARENT_COMBO_9",	"RI_IND_PARENT_COMBO_10",	"RI_IND_PARENT_COMBO_11",	"RI_IND_PARENT_COMBO_12",	"RI_IND_PARENT_COMBO_13",	"RI_IND_PARENT_COMBO_14",	"RI_IND_PARENT_COMBO_15",	"RI_IND_GOT_BOTH_ORIS",	"RI_IND_GOT_PREFIX",	"RI_IND_GOT_SUFFIX",	
"RI_IND_PARENT1_NOT_VIZ",	"RI_IND_PARENT1_NO_INTEN",	"RI_PARENT1_ISO_LEVEL",	"RI_IND_PARENT1_INTEN_MORE",	"RI_PARENT1_INTEN_DIFF_MORE",	"RI_IND_PARENT1_INTEN_LESS",	"RI_PARENT1_INTEN_DIFF_LESS",	"RI_IND_PARENT2_NOT_VIZ",	"RI_IND_PARENT2_NO_INTEN",	"RI_PARENT2_ISO_LEVEL",	"RI_IND_PARENT2_INTEN_MORE",	"RI_PARENT2_INTEN_DIFF_MORE",	"RI_IND_PARENT2_INTEN_LESS",	"RI_PARENT2_INTEN_DIFF_LESS",	
"RI_IND_PARENT3_NOT_VIZ",	"RI_IND_PARENT3_NO_INTEN",	"RI_PARENT3_ISO_LEVEL",	"RI_IND_PARENT3_INTEN_MORE",	"RI_PARENT3_INTEN_DIFF_MORE",	"RI_IND_PARENT3_INTEN_LESS",	"RI_PARENT3_INTEN_DIFF_LESS",	"RI_IND_PARENT4_NOT_VIZ",	"RI_IND_PARENT4_NO_INTEN",	"RI_PARENT4_ISO_LEVEL",	"RI_IND_PARENT4_INTEN_MORE",	"RI_PARENT4_INTEN_DIFF_MORE",	"RI_IND_PARENT4_INTEN_LESS",	"RI_PARENT4_INTEN_DIFF_LESS",	
"RI_IND_PARENT5_NOT_VIZ",	"RI_IND_PARENT5_NO_INTEN",	"RI_PARENT5_ISO_LEVEL",	"RI_IND_PARENT5_INTEN_MORE",	"RI_PARENT5_INTEN_DIFF_MORE",	"RI_IND_PARENT5_INTEN_LESS",	"RI_PARENT5_INTEN_DIFF_LESS",	"RI_IND_PARENT6_NOT_VIZ",	"RI_IND_PARENT6_NO_INTEN",	"RI_PARENT6_ISO_LEVEL",	"RI_IND_PARENT6_INTEN_MORE",	"RI_PARENT6_INTEN_DIFF_MORE",	"RI_IND_PARENT6_INTEN_LESS",	"RI_PARENT6_INTEN_DIFF_LESS",	
"RI_IND_PARENT7_NOT_VIZ",	"RI_IND_PARENT7_NO_INTEN",	"RI_PARENT7_ISO_LEVEL",	"RI_IND_PARENT7_INTEN_MORE",	"RI_PARENT7_INTEN_DIFF_MORE",	"RI_IND_PARENT7_INTEN_LESS",	"RI_PARENT7_INTEN_DIFF_LESS",	"RI_IND_PARENT8_NOT_VIZ",	"RI_IND_PARENT8_NO_INTEN",	"RI_PARENT8_ISO_LEVEL",	"RI_IND_PARENT8_INTEN_MORE",	"RI_PARENT8_INTEN_DIFF_MORE",	"RI_IND_PARENT8_INTEN_LESS",	"RI_PARENT8_INTEN_DIFF_LESS",	
"RI_IND_N_IS_GAP",	"RI_IND_C_IS_GAP",	"RI_IND_N_HAS_N_TERM",	"RI_IND_N_HAS_C_TERM",	"RI_IND_N_HAS_Gap",	"RI_IND_N_HAS_Xle",	"RI_IND_N_HAS_Ala",	"RI_IND_N_HAS_Arg",	"RI_IND_N_HAS_Asn",	"RI_IND_N_HAS_Asp",	"RI_IND_N_HAS_Cys",	"RI_IND_N_HAS_Gln",	"RI_IND_N_HAS_Glu",	"RI_IND_N_HAS_Gly",	"RI_IND_N_HAS_His",	"RI_IND_N_HAS_Ile",	"RI_IND_N_HAS_Leu",	"RI_IND_N_HAS_Lys",	"RI_IND_N_HAS_Met",	
"RI_IND_N_HAS_Phe",	"RI_IND_N_HAS_Pro",	"RI_IND_N_HAS_Ser",	"RI_IND_N_HAS_Thr",	"RI_IND_N_HAS_Trp",	"RI_IND_N_HAS_Tyr",	"RI_IND_N_HAS_Val",	"RI_N_N_TERM_SELF_INTEN",	"RI_N_C_TERM_SELF_INTEN",	"RI_N_Gap_SELF_INTEN",	"RI_N_Xle_SELF_INTEN",	"RI_N_Ala_SELF_INTEN",	"RI_N_Arg_SELF_INTEN",	"RI_N_Asn_SELF_INTEN",	"RI_N_Asp_SELF_INTEN",	"RI_N_Cys_SELF_INTEN",	"RI_N_Gln_SELF_INTEN",	"RI_N_Glu_SELF_INTEN",	
"RI_N_Gly_SELF_INTEN",	"RI_N_His_SELF_INTEN",	"RI_N_Ile_SELF_INTEN",	"RI_N_Leu_SELF_INTEN",	"RI_N_Lys_SELF_INTEN",	"RI_N_Met_SELF_INTEN",	"RI_N_Phe_SELF_INTEN",	"RI_N_Pro_SELF_INTEN",	"RI_N_Ser_SELF_INTEN",	"RI_N_Thr_SELF_INTEN",	"RI_N_Trp_SELF_INTEN",	"RI_N_Tyr_SELF_INTEN",	"RI_N_Val_SELF_INTEN",	"RI_IND_C_HAS_N_TERM",	"RI_IND_C_HAS_C_TERM",	"RI_IND_C_HAS_Gap",	"RI_IND_C_HAS_Xle",	
"RI_IND_C_HAS_Ala",	"RI_IND_C_HAS_Arg",	"RI_IND_C_HAS_Asn",	"RI_IND_C_HAS_Asp",	"RI_IND_C_HAS_Cys",	"RI_IND_C_HAS_Gln",	"RI_IND_C_HAS_Glu",	"RI_IND_C_HAS_Gly",	"RI_IND_C_HAS_His",	"RI_IND_C_HAS_Ile",	"RI_IND_C_HAS_Leu",	"RI_IND_C_HAS_Lys",	"RI_IND_C_HAS_Met",	"RI_IND_C_HAS_Phe",	"RI_IND_C_HAS_Pro",	"RI_IND_C_HAS_Ser",	"RI_IND_C_HAS_Thr",	"RI_IND_C_HAS_Trp",	"RI_IND_C_HAS_Tyr",	
"RI_IND_C_HAS_Val",	"RI_C_N_TERM_SELF_INTEN",	"RI_C_C_TERM_SELF_INTEN",	"RI_C_Gap_SELF_INTEN",	"RI_C_Xle_SELF_INTEN",	"RI_C_Ala_SELF_INTEN",	"RI_C_Arg_SELF_INTEN",	"RI_C_Asn_SELF_INTEN",	"RI_C_Asp_SELF_INTEN",	"RI_C_Cys_SELF_INTEN",	"RI_C_Gln_SELF_INTEN",	"RI_C_Glu_SELF_INTEN",	"RI_C_Gly_SELF_INTEN",	"RI_C_His_SELF_INTEN",	"RI_C_Ile_SELF_INTEN",	"RI_C_Leu_SELF_INTEN",	"RI_C_Lys_SELF_INTEN",	
"RI_C_Met_SELF_INTEN",	"RI_C_Phe_SELF_INTEN",	"RI_C_Pro_SELF_INTEN",	"RI_C_Ser_SELF_INTEN",	"RI_C_Thr_SELF_INTEN",	"RI_C_Trp_SELF_INTEN",	"RI_C_Tyr_SELF_INTEN",	"RI_C_Val_SELF_INTEN",	"RI_NUM_FEATURES",	"ScoreModelFields_RI"};

const char * ScoreModelFields_RNI_names[]={
"RNI_CONST",	"RNI_IND_DIS_FROM_MINMAX_LESS_50",	"RNI_DIS_FROM_MINMAX0",	"RNI_IND_DIS_FROM_MINMAX_LESS_150",	"RNI_DIS_FROM_MINMAX50",	"RNI_IND_DIS_FROM_MINMAX_LESS_250",	"RNI_DIS_FROM_MINMAX150",	"RNI_IND_DIS_FROM_MINMAX_MORE",	"RNI_DIS_FROM_MINMAX250",	"RNI_REL_POS0",	"RNI_REL_POS1",	"RNI_REL_POS2",	"RNI_REL_POS3",	"RNI_REL_POS4",	"RNI_REL_POS5",	"RNI_REL_POS6",	"RNI_REL_POS7",	
"RNI_REL_POS8",	"RNI_REL_POS9",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_0",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_1",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_2",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_3",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_4",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_5",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_6",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_MORE6",	"RNI_IND_PARENT_COMBO_0",	"RNI_IND_PARENT_COMBO_1",	
"RNI_IND_PARENT_COMBO_2",	"RNI_IND_PARENT_COMBO_3",	"RNI_IND_PARENT_COMBO_4",	"RNI_IND_PARENT_COMBO_5",	"RNI_IND_PARENT_COMBO_6",	"RNI_IND_PARENT_COMBO_7",	"RNI_IND_PARENT_COMBO_8",	"RNI_IND_PARENT_COMBO_9",	"RNI_IND_PARENT_COMBO_10",	"RNI_IND_PARENT_COMBO_11",	"RNI_IND_PARENT_COMBO_12",	"RNI_IND_PARENT_COMBO_13",	"RNI_IND_PARENT_COMBO_14",	"RNI_IND_PARENT_COMBO_15",	"RNI_IND_GOT_BOTH_ORIS",	
"RNI_IND_GOT_PREFIX",	"RNI_IND_GOT_SUFFIX",	"RNI_IND_PARENT1_NOT_VIZ",	"RNI_IND_PARENT1_NO_INTEN",	"RNI_PARENT1_LOG_INTEN",	"RNI_PARENT1_LOG_GLOBAL_RANK",	"RNI_PARENT1_ISO_LEVEL",	"RNI_IND_PARENT2_NOT_VIZ",	"RNI_IND_PARENT2_NO_INTEN",	"RNI_PARENT2_LOG_INTEN",	"RNI_PARENT2_LOG_GLOBAL_RANK",	"RNI_PARENT2_ISO_LEVEL",	"RNI_IND_PARENT3_NOT_VIZ",	"RNI_IND_PARENT3_NO_INTEN",	"RNI_PARENT3_LOG_INTEN",	
"RNI_PARENT3_LOG_GLOBAL_RANK",	"RNI_PARENT3_ISO_LEVEL",	"RNI_IND_PARENT4_NOT_VIZ",	"RNI_IND_PARENT4_NO_INTEN",	"RNI_PARENT4_LOG_INTEN",	"RNI_PARENT4_LOG_GLOBAL_RANK",	"RNI_PARENT4_ISO_LEVEL",	"RNI_IND_PARENT5_NOT_VIZ",	"RNI_IND_PARENT5_NO_INTEN",	"RNI_PARENT5_LOG_INTEN",	"RNI_PARENT5_LOG_GLOBAL_RANK",	"RNI_PARENT5_ISO_LEVEL",	"RNI_IND_PARENT6_NOT_VIZ",	"RNI_IND_PARENT6_NO_INTEN",	
"RNI_PARENT6_LOG_INTEN",	"RNI_PARENT6_LOG_GLOBAL_RANK",	"RNI_PARENT6_ISO_LEVEL",	"RNI_IND_PARENT7_NOT_VIZ",	"RNI_IND_PARENT7_NO_INTEN",	"RNI_PARENT7_LOG_INTEN",	"RNI_PARENT7_LOG_GLOBAL_RANK",	"RNI_PARENT7_ISO_LEVEL",	"RNI_IND_PARENT8_NOT_VIZ",	"RNI_IND_PARENT8_NO_INTEN",	"RNI_PARENT8_LOG_INTEN",	"RNI_PARENT8_LOG_GLOBAL_RANK",	"RNI_PARENT8_ISO_LEVEL",	"RNI_IND_N_IS_GAP",	"RNI_IND_C_IS_GAP",	
"RNI_IND_N_HAS_N_TERM",	"RNI_IND_N_HAS_C_TERM",	"RNI_IND_N_HAS_Gap",	"RNI_IND_N_HAS_Xle",	"RNI_IND_N_HAS_Ala",	"RNI_IND_N_HAS_Arg",	"RNI_IND_N_HAS_Asn",	"RNI_IND_N_HAS_Asp",	"RNI_IND_N_HAS_Cys",	"RNI_IND_N_HAS_Gln",	"RNI_IND_N_HAS_Glu",	"RNI_IND_N_HAS_Gly",	"RNI_IND_N_HAS_His",	"RNI_IND_N_HAS_Ile",	"RNI_IND_N_HAS_Leu",	"RNI_IND_N_HAS_Lys",	"RNI_IND_N_HAS_Met",	"RNI_IND_N_HAS_Phe",	
"RNI_IND_N_HAS_Pro",	"RNI_IND_N_HAS_Ser",	"RNI_IND_N_HAS_Thr",	"RNI_IND_N_HAS_Trp",	"RNI_IND_N_HAS_Tyr",	"RNI_IND_N_HAS_Val",	"RNI_IND_C_HAS_N_TERM",	"RNI_IND_C_HAS_C_TERM",	"RNI_IND_C_HAS_Gap",	"RNI_IND_C_HAS_Xle",	"RNI_IND_C_HAS_Ala",	"RNI_IND_C_HAS_Arg",	"RNI_IND_C_HAS_Asn",	"RNI_IND_C_HAS_Asp",	"RNI_IND_C_HAS_Cys",	"RNI_IND_C_HAS_Gln",	"RNI_IND_C_HAS_Glu",	"RNI_IND_C_HAS_Gly",	
"RNI_IND_C_HAS_His",	"RNI_IND_C_HAS_Ile",	"RNI_IND_C_HAS_Leu",	"RNI_IND_C_HAS_Lys",	"RNI_IND_C_HAS_Met",	"RNI_IND_C_HAS_Phe",	"RNI_IND_C_HAS_Pro",	"RNI_IND_C_HAS_Ser",	"RNI_IND_C_HAS_Thr",	"RNI_IND_C_HAS_Trp",	"RNI_IND_C_HAS_Tyr",	"RNI_IND_C_HAS_Val",	"RNI_NUM_FEATURES",	"ScoreModelFields_RNI"};




bool RegularFragmentModel::write_model(ostream& os) const
{
	if (! ind_has_models)
		return false;

	os << model_frag_idx << " " << model_frag_charge << endl;

	if (fabs(inten_log_scaling_factor)>5 || fabs(no_inten_log_scaling_factor)>5)
	{
		os <<  NEG_INF << " " << NEG_INF << endl;
	}
	else
		os << setprecision(4) << inten_log_scaling_factor << " " << no_inten_log_scaling_factor << endl;

	os << parent_idxs.size() << " " << parent_idx_with_same_charge_ori << endl;
	int i;
	for (i=0; i<parent_idxs.size() && i<num_parents; i++)
		os << parent_idxs[i] << " ";
	os << endl;

	inten_model.write_regression_model(os);
	no_inten_model.write_regression_model(os);

	return true;
}


bool RegularFragmentModel::read_model(const Config* config, istream& is, bool silent_ind)
{
	config_ = config;

	char buff[256];

	is.getline(buff,256);
	istringstream iss(buff);
	iss >> model_frag_idx >> model_frag_charge;
	
	is.getline(buff,256);
	sscanf(buff,"%f %f",&inten_log_scaling_factor,&no_inten_log_scaling_factor);

	is.getline(buff,256);
	iss.str(buff);
	iss >> num_parents >> parent_idx_with_same_charge_ori;

	is.getline(buff,256);
	iss.str(buff);
	int i;
	parent_idxs.resize(num_parents);
	for (i=0; i<num_parents; i++)
		iss >> parent_idxs[i];

	inten_model.read_regression_model(is);
	no_inten_model.read_regression_model(is);

	ind_has_models = true;

	// check if something was wrong with this fragment
	if (fabs(inten_log_scaling_factor)>5 || fabs(no_inten_log_scaling_factor)>5)
	{
		if (! silent_ind)
		{
			cout << "Model for frag " << model_frag_idx << " (" << 
				config_->get_fragment(model_frag_idx).label << ") had scaling problems: " << endl;
			cout << inten_log_scaling_factor << " " << no_inten_log_scaling_factor << endl;
		}
		ind_has_models=false;
	}

	if (! inten_model.get_has_weights() || ! no_inten_model.get_has_weights())
		ind_has_models=false;

	return true;
}



void RegularFragmentModel::fill_constant_vals(
							Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage, 
							vector<fval>& f_vals) const
{
	const int num_parents = parent_idxs.size();
	FragStats frag_stats;
	vector<FragStats> parent_frag_stasts;

	parent_frag_stasts.resize(num_parents);
	frag_stats.fill_from_breakage(breakage,spec,model_frag_idx);
	int i;
	for (i=0; i<num_parents; i++)
		parent_frag_stasts[i].fill_from_breakage(breakage,spec,parent_idxs[i]);
	
	f_vals.clear();

	if (frag_stats.has_intensity) // fill features for visible fragment
	{
		f_vals.push_back(fval(RI_CONST,1.0));

		const float log_inten = frag_stats.log_intensity;
		f_vals.push_back(fval(RI_LOG_LOCAL_RANK,frag_stats.log_local_rank));
		f_vals.push_back(fval(RI_LOG_GLOBAL_RANK,frag_stats.log_global_rank));
		if (frag_stats.iso_level != 0)
			f_vals.push_back(fval(RI_ISO_LEVEL,frag_stats.iso_level));

		if (log_inten<1.0)
		{
			f_vals.push_back(fval(RI_IND_LOG_INTEN_LESS1,1.0));
			f_vals.push_back(fval(RI_LOG_INTEN_LESS1,log_inten));
		} 
		else if (log_inten<2.0)
		{
			f_vals.push_back(fval(RI_IND_LOG_INTEN_LESS2,1.0));
			f_vals.push_back(fval(RI_LOG_INTEN_LESS2,log_inten-1.0));
		}
		else if (log_inten<3.0)
		{
			f_vals.push_back(fval(RI_IND_LOG_INTEN_LESS3,1.0));
			f_vals.push_back(fval(RI_LOG_INTEN_LESS3,log_inten-2.0));
		}
		else if (log_inten<4.0)
		{
			f_vals.push_back(fval(RI_IND_LOG_INTEN_LESS4,1.0));
			f_vals.push_back(fval(RI_LOG_INTEN_LESS4,log_inten-3.0));
		}
		else
		{
			f_vals.push_back(fval(RI_IND_LOG_INTEN_MORE,1.0));
			f_vals.push_back(fval(RI_LOG_INTEN_MORE,log_inten-4.0));
		}


			// self distance
		const mass_t dis_min = frag_stats.mass - spec->get_min_peak_mass();
		const mass_t dis_max = spec->get_max_peak_mass() - frag_stats.mass;
		const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);

		if (dis<50)
		{
			f_vals.push_back(fval(RI_IND_DIS_FROM_MINMAX_LESS_50,1.0));
			f_vals.push_back(fval(RI_DIS_FROM_MINMAX0,dis));
			f_vals.push_back(fval(RI_LOG_INTEN_DIS50,log_inten));
		}
		else if (dis<150)
		{
			f_vals.push_back(fval(RI_IND_DIS_FROM_MINMAX_LESS_150,1.0));
			f_vals.push_back(fval(RI_DIS_FROM_MINMAX50,dis-50.0));
			f_vals.push_back(fval(RI_LOG_INTEN_DIS150,log_inten));
		}
		else if (dis<250)
		{
			f_vals.push_back(fval(RI_IND_DIS_FROM_MINMAX_LESS_250,1.0));
			f_vals.push_back(fval(RI_DIS_FROM_MINMAX150,dis-150.0));
			f_vals.push_back(fval(RI_LOG_INTEN_DIS250,log_inten));
		}
		else
		{
			f_vals.push_back(fval(RI_IND_DIS_FROM_MINMAX_MORE,1.0));
			f_vals.push_back(fval(RI_DIS_FROM_MINMAX250,dis-250.0));
			f_vals.push_back(fval(RI_LOG_INTEN_DISMORE,log_inten));
		}

		const int rel_pos = int(10*breakage->mass/pm_with_19);
		f_vals.push_back(fval(RI_REL_POS0+rel_pos,1.0));

		bool got_prefix = false;
		bool got_suffix = false;
		int combo_idx=0;
		int pow=1;
		int num_parents_with_intensity=0;
		int i;
		for (i=0; i<num_parents; i++)
		{
			if (parent_frag_stasts[i].has_intensity)
			{
				num_parents_with_intensity++;

				if (config_->get_fragment(parent_idxs[i]).orientation == PREFIX)
				{
					got_prefix = true;
				}
				else
					got_suffix = true;

				if (i<4)
					combo_idx += pow;
			}
			pow*=2;
		}

		if (num_parents_with_intensity>6)
			num_parents_with_intensity=7;

		f_vals.push_back(fval(RI_IND_NUM_PARENTS_WITH_INTEN_IS_0+num_parents_with_intensity,1.0));

		f_vals.push_back(fval(RI_IND_PARENT_COMBO_0+combo_idx,1.0));

		if (got_prefix && got_suffix)
		{
			f_vals.push_back(fval(RI_IND_GOT_BOTH_ORIS,1.0));	   
		}
		else if (got_prefix)
		{
			f_vals.push_back(fval(RI_IND_GOT_PREFIX,1.0)); 
		}
		else if (got_suffix)
			f_vals.push_back(fval(RI_IND_GOT_SUFFIX,1.0));

		for (i=0; i<num_parents; i++)
		{
			const int offset_idx = 7*i;
			if (! parent_frag_stasts[i].is_viz)
			{
				f_vals.push_back(fval(RI_IND_PARENT1_NOT_VIZ+offset_idx,1.0));
			}
			else
			{
				if (! parent_frag_stasts[i].has_intensity)
				{
					f_vals.push_back(fval(RI_IND_PARENT1_NO_INTEN+offset_idx,1.0));
				}
				else
				{
					if (parent_frag_stasts[i].iso_level != 0)
						f_vals.push_back(fval(RI_PARENT1_ISO_LEVEL+offset_idx,parent_frag_stasts[i].iso_level));

					const float parent_log = parent_frag_stasts[i].log_intensity;
					if (parent_log>log_inten)
					{
						f_vals.push_back(fval(RI_IND_PARENT1_INTEN_MORE+offset_idx,1.0));
						f_vals.push_back(fval(RI_PARENT1_INTEN_DIFF_MORE+offset_idx,parent_log-log_inten));
					}
					else
					{
						f_vals.push_back(fval(RI_IND_PARENT1_INTEN_LESS+offset_idx,1.0));
						f_vals.push_back(fval(RI_PARENT1_INTEN_DIFF_LESS+offset_idx,parent_log-log_inten));
					}
				}
			}
		}
	}
	else // Fill features for non-visible
	{
		const mass_t expected_mass = config_->get_fragment(model_frag_idx).calc_expected_mass(breakage->mass,pm_with_19);
		const mass_t dis_min = expected_mass - spec->get_min_peak_mass();
		const mass_t dis_max = spec->get_max_peak_mass() - expected_mass;
		const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);

		f_vals.push_back(fval(RNI_CONST,1.0));

		if (dis<50)
		{
			f_vals.push_back(fval(RNI_IND_DIS_FROM_MINMAX_LESS_50,1.0));
			f_vals.push_back(fval(RNI_DIS_FROM_MINMAX0,dis));
		}
		else if (dis<150)
		{
			f_vals.push_back(fval(RNI_IND_DIS_FROM_MINMAX_LESS_150,1.0));
			f_vals.push_back(fval(RNI_DIS_FROM_MINMAX50,dis-50.0));
		}
		else if (dis<250)
		{
			f_vals.push_back(fval(RNI_IND_DIS_FROM_MINMAX_LESS_250,1.0));
			f_vals.push_back(fval(RNI_DIS_FROM_MINMAX150,dis-150.0));
		}
		else
		{
			f_vals.push_back(fval(RNI_IND_DIS_FROM_MINMAX_MORE,1.0));
			f_vals.push_back(fval(RNI_DIS_FROM_MINMAX250,dis-250.0));
		}

		const int rel_pos = int(10*breakage->mass/pm_with_19);
		f_vals.push_back(fval(RNI_REL_POS0+rel_pos,1.0));


		bool got_prefix=false;
		bool got_suffix=false;
		int combo_idx=0;
		int pow=1;
		int num_parents_with_intensity=0;
		int i;
		for (i=0; i<num_parents; i++)
		{
			if (parent_frag_stasts[i].has_intensity)
			{
				num_parents_with_intensity++;

				if (config_->get_fragment(parent_idxs[i]).orientation == PREFIX)
				{
					got_prefix = true;
				}
				else
					got_suffix = true;

				if (i<4)
					combo_idx += pow;
			}
			pow*=2;
		}

		if (num_parents_with_intensity>6)
			num_parents_with_intensity=7;

		f_vals.push_back(fval(RNI_IND_NUM_PARENTS_WITH_INTEN_IS_0+num_parents_with_intensity,1.0));

		f_vals.push_back(fval(RNI_IND_PARENT_COMBO_0+combo_idx,1.0));

		if (got_prefix && got_suffix)
		{
			f_vals.push_back(fval(RI_IND_GOT_BOTH_ORIS,1.0));	   
		}
		else if (got_prefix)
		{
			f_vals.push_back(fval(RI_IND_GOT_PREFIX,1.0)); 
		}
		else if (got_suffix)
			f_vals.push_back(fval(RI_IND_GOT_SUFFIX,1.0));

		for (i=0; i<num_parents; i++)
		{
			const int offset_idx = 5*i;
			if (! parent_frag_stasts[i].is_viz)
			{
				f_vals.push_back(fval(RNI_IND_PARENT1_NOT_VIZ+offset_idx,1.0));
			}
			else
			{
				if (! parent_frag_stasts[i].has_intensity)
				{
					f_vals.push_back(fval(RNI_IND_PARENT1_NO_INTEN+offset_idx,1.0));
				}
				else
				{
					f_vals.push_back(fval(RNI_PARENT1_LOG_INTEN+offset_idx,parent_frag_stasts[i].log_intensity));
					f_vals.push_back(fval(RNI_PARENT1_LOG_GLOBAL_RANK+offset_idx,parent_frag_stasts[i].log_global_rank));
					if (parent_frag_stasts[i].iso_level != 0)
						f_vals.push_back(fval(RNI_PARENT1_ISO_LEVEL+offset_idx,parent_frag_stasts[i].iso_level));
				}
			}
		}	
	}
}






/**********************************************************************************
***********************************************************************************/
void RegularFragmentModel::fill_aa_variable_vals(
							   Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage,
							   const BreakageInfo* info,
							   vector<fval>& f_vals) const
{
	const vector<int>& org_aas = config_->get_org_aa();
	const int n_aa = (info->n_aa>=0 ? org_aas[info->n_aa] : Gap);
	const int c_aa = (info->c_aa>=0 ? org_aas[info->c_aa] : Gap);
	const int pos = breakage->get_position_of_frag_idx(model_frag_idx);
	
	const bool do_n_features = (n_aa != Gap);
	const bool do_c_features = (c_aa != Gap);


	// fill intensity
	if (pos>=0)
	{
		if (! do_n_features)
			f_vals.push_back(fval(RI_IND_N_IS_GAP,1.0));
		
		if (! do_c_features)
			f_vals.push_back(fval(RI_IND_C_IS_GAP,1.0));

		const float log_inten = spec->get_peak_log_intensity(breakage->fragments[pos].peak_idx);

		// add aa presence indicators
		if (do_n_features)
		{
			int *var_ptr = info->n_var_ptr;
			const int num_aa = *var_ptr++;
			int k;
			for (k=0; k<num_aa; k++)
			{
				const int n_aa = org_aas[var_ptr[k]];
				f_vals.push_back(fval(RI_IND_N_HAS_N_TERM,1.0));
				f_vals.push_back(fval(RI_IND_N_HAS_N_TERM+n_aa,1.0));
				f_vals.push_back(fval(RI_N_N_TERM_SELF_INTEN,log_inten));
				f_vals.push_back(fval(RI_N_N_TERM_SELF_INTEN+n_aa,log_inten));
			}
		}

		if (do_c_features)
		{
			int *var_ptr = info->c_var_ptr;
			const int num_aa = *var_ptr++;
			int k;
			for (k=0; k<num_aa; k++)
			{
				const int c_aa = org_aas[var_ptr[k]];
				f_vals.push_back(fval(RI_IND_C_HAS_N_TERM,1.0));
				f_vals.push_back(fval(RI_IND_C_HAS_N_TERM+c_aa,1.0));
				f_vals.push_back(fval(RI_C_N_TERM_SELF_INTEN,log_inten));
				f_vals.push_back(fval(RI_C_N_TERM_SELF_INTEN+c_aa,log_inten));
			}
		}
	}
	else // Fill no frag features
	{
		if (! do_n_features)
			f_vals.push_back(fval(RNI_IND_N_IS_GAP,1.0));
		
		if (! do_c_features)
			f_vals.push_back(fval(RNI_IND_C_IS_GAP,1.0));


		// add aa presence indicators
		if (do_n_features)
		{
			int *var_ptr = info->n_var_ptr;
			const int num_aa = *var_ptr++;
			int k;
			for (k=0; k<num_aa; k++)
			{
				const int n_aa = org_aas[var_ptr[k]];
				f_vals.push_back(fval(RNI_IND_N_HAS_N_TERM,1.0));
				f_vals.push_back(fval(RNI_IND_N_HAS_N_TERM+n_aa,1.0));
			}
		}

		if (do_c_features)
		{
			int *var_ptr = info->c_var_ptr;
			const int num_aa = *var_ptr++;
			int k;
			for (k=0; k<num_aa; k++)
			{
				const int c_aa = org_aas[var_ptr[k]];
				f_vals.push_back(fval(RNI_IND_C_HAS_N_TERM,1.0));
				f_vals.push_back(fval(RNI_IND_C_HAS_N_TERM+c_aa,1.0));
			}
		}
	}
}




void RegularFragmentModel::fill_combo_vectors(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							const vector<BreakageInfo>& infos,
							vector< ME_Regression_Sample > & samples) const
{
	vector<fval> const_vals;
	fill_constant_vals(spec,pm_with_19,breakage,const_vals);
	samples.resize(infos.size());
	int i;
	for (i=0; i<infos.size(); i++)
	{
		vector< fval > var_vals;
		fill_aa_variable_vals(spec,pm_with_19,breakage,&infos[i],var_vals);

		samples[i].f_vals = const_vals;
		int j;
		for (j=0; j<var_vals.size(); j++)
			samples[i].f_vals.push_back(var_vals[j]);
	}
}



void RegularFragmentModel::fill_single_frag_vector(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							BreakageInfo& info,
							vector< fval >& f_vals) const
{

	fill_constant_vals(spec,pm_with_19,breakage,f_vals);
	fill_aa_variable_vals(spec,pm_with_19,breakage,&info,f_vals);
	sort(f_vals.begin(),f_vals.end());
}

