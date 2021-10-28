#include "RegionalPrmNodeScoreModel.h"
#include "PrmGraph.h"
#include "PepNovo_auxfun.h"
#include "AllScoreModels.h"

extern const char* ScoreModelFields_SI_names[];
extern const char* ScoreModelFields_SNI_names[];
extern const char* ScoreModelFields_RI_names[];
extern const char* ScoreModelFields_RNI_names[];


/****************************************************************************
Returns the score of the constant element in the breakage:
strong features tha don't depend on aa, and all weak features
*****************************************************************************/
void RegionalPrmNodeScoreModel::calc_constant_element(
							   Node& node,
							   Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage) const
{
	node.const_strong_exps.clear();
	node.const_strong_exps.resize(strong_models.size(),0);

	int f;
	for (f=0; f<this->strong_models.size(); f++)
	{
		const StrongFragmentModel& strong_model = strong_models[f];
		const int frag_idx = strong_model.model_frag_idx;

		if (! strong_model.ind_has_models)
			continue;

		if (! breakage->is_frag_type_visible(frag_idx))
			continue;

		ME_Regression_Sample sam;	
		strong_model.fill_constant_vals(spec,pm_with_19,breakage,sam.f_vals);

		const int pos = breakage->get_position_of_frag_idx(frag_idx);
		if (pos>=0)
		{
			node.const_strong_exps[f]=strong_model.inten_model.get_sum_exp(sam);
		}
		else
			node.const_strong_exps[f]=strong_model.no_inten_model.get_sum_exp(sam);
	}

	node.const_regular_exps.clear();
	node.const_regular_exps.resize(regular_models.size(),0);

	for (f=0; f<regular_models.size(); f++)
	{
		const RegularFragmentModel& regular_model = regular_models[f];
		const int frag_idx = regular_model.model_frag_idx;

		if (! regular_model.ind_has_models)
			continue;

		if (! breakage->is_frag_type_visible(frag_idx))
			continue;

		ME_Regression_Sample sam;
		regular_model.fill_constant_vals(spec,pm_with_19,breakage,sam.f_vals);

		const int pos = breakage->get_position_of_frag_idx(frag_idx);
		if (pos>=0)
		{
			node.const_regular_exps[f] = regular_model.inten_model.get_sum_exp(sam);
		}
		else
			node.const_regular_exps[f] = regular_model.no_inten_model.get_sum_exp(sam);
	}
}


/****************************************************************************
Returns the score of the variable element in the breakages:
the strong features that are aa dependant.
*****************************************************************************/
score_t RegionalPrmNodeScoreModel::score_a_single_breakage_combo(
							   PrmGraph *prm,  
							   Node& node,
							   const Breakage *breakage,
							   BreakageInfo& info,
							   bool verbose) const
{
	if (node.type == NODE_N_TERM || node.type == NODE_C_TERM)
	{
		const score_t terminal_score=config_->get_terminal_score();
		info.score = terminal_score;
		return terminal_score;
	}
	
	Spectrum *spec = prm->get_source_spectrum();
	const mass_t pm_with_19 = prm->get_pm_with_19();

	if (node.const_strong_exps.size()==0)
		calc_constant_element(node,spec,pm_with_19,breakage);

	const bool print_all = false;

	score_t score=0;

	int f;
	for (f=0; f<this->strong_models.size(); f++)
	{
		const StrongFragmentModel& strong_model = strong_models[f];
		const int frag_idx = strong_model.model_frag_idx;

		if (! strong_model.ind_has_models)
			continue;

		if (! breakage->is_frag_type_visible(frag_idx))
			continue;

		ME_Regression_Sample sam;
		sam.f_vals.clear();
		strong_model.fill_aa_variable_vals(spec,pm_with_19,breakage,&info,sam.f_vals);
	
		score_t prev = score;
		const int pos = breakage->get_position_of_frag_idx(frag_idx);
		if (pos>=0)
		{
			const float var_exp = strong_model.inten_model.get_sum_exp(sam);
			const float e = exp(var_exp + node.const_strong_exps[f]);

			float prob = e/(1.0 + e);
			if (prob>0.99)
				prob=0.99;
			if (prob<0.001)
				prob=0.001;

			const float log_random_peak = spec->getPeakLogRandomProbability(breakage->fragments[pos].peak_idx);
			score += (strong_model.inten_log_scaling_factor + log(prob) - log_random_peak);

			if (print_all)
			{
				cout << "viz >> SCORE " << score << "\t(lrp= " << log_random_peak << " , ilsf=" << strong_model.inten_log_scaling_factor <<
					" , lp=" << log(prob) << " , lrandp=" << log_random_peak << endl;
			}
		}
		else
		{
			const float var_exp = strong_model.no_inten_model.get_sum_exp(sam);
			const float e = exp(var_exp + node.const_strong_exps[f]);

			float prob = e/(1.0 + e);
			if (prob>0.99)
				prob=0.99;
			if (prob<0.01)
				prob=0.01;

			score += (strong_model.no_inten_log_scaling_factor + log(prob) - log_one_minus_random);

			if (print_all)
			{
				cout << "no viz >> SCORE " << score << "\t(nilsf=" << strong_model.inten_log_scaling_factor <<
					" , lp=" << log(prob) << " , loneminusrandp=" << log_one_minus_random << endl;
			}
		}

		if (verbose)
		{
			cout << setprecision(4) << fixed << f << "\t" << score-prev << endl;
		}
	}

	for (f=0; f<regular_models.size(); f++)
	{
		const RegularFragmentModel& regular_model = regular_models[f];
		const int frag_idx = regular_model.model_frag_idx;

		if (! regular_model.ind_has_models)
			continue;

		if (! breakage->is_frag_type_visible(frag_idx))
			continue;

		ME_Regression_Sample sam;
		sam.f_vals.clear();
		regular_model.fill_aa_variable_vals(spec,pm_with_19,breakage,&info,sam.f_vals);
		
		score_t prev= score;

		const int pos = breakage->get_position_of_frag_idx(frag_idx);
		if (pos>=0)
		{
			const float var_exp = regular_model.inten_model.get_sum_exp(sam);
			const float e = exp(var_exp + node.const_regular_exps[f]);

			float prob = e/(1.0 + e);
			if (prob>0.99)
				prob=0.99;
			if (prob<0.01)
				prob=0.01;

			const float log_random_peak = spec->getPeakLogRandomProbability(breakage->fragments[pos].peak_idx);
			score += (regular_model.inten_log_scaling_factor + log(prob) - log_random_peak);
		}
		else
		{
			const float var_exp = regular_model.no_inten_model.get_sum_exp(sam);
			const float e = exp(var_exp + node.const_regular_exps[f]);

			float prob = e/(1.0 + e);
			if (prob>0.99)
				prob=0.99;
			if (prob<0.02)
				prob=0.02;

			score += (regular_model.no_inten_log_scaling_factor + log(prob) - log_one_minus_random);
		}

		if (verbose)
		{
			cout << f << "\t" << score-prev << endl;
		}
	}

	// correct for digest node scores
	if (node.type == NODE_DIGEST)
	{
		const score_t digest_score=config_->get_digest_score();

		if ((info.connects_to_N_term && info.n_edge_idx<0) || 
			(info.connects_to_C_term && info.c_edge_idx<0) )
		{
			score -= digest_score; // penalty for ending here!
		}
			
		if (info.connects_to_N_term && info.preferred_digest_aa_N_term && score<0)
			score = 0;
		if (info.connects_to_C_term && info.preferred_digest_aa_C_term && score<0)
			score = 0;
	}

	return score;
}

string RegionalPrmNodeScoreModel::make_model_file_name(const char *name) const
{
	char dir_path[256];
	char model_name[64];
	strcpy(dir_path,config_->get_resource_dir().c_str());
	strcat(dir_path,"/");
	strcat(dir_path,name);
	strcat(dir_path,"_SCORE/");
	sprintf(model_name,"%s_%d_%d_%d.txt", name, charge_, sizeIndex_, regionIndex_);
	string file_path = dir_path;
	file_path += model_name;

	return file_path;
}






bool RegionalPrmNodeScoreModel::write_regional_score_model(const char *name) const
{
	if (! was_initialized)
		return false;

	string file_path = make_model_file_name(name);
	ofstream os(file_path.c_str());

	if (! os.good() || ! os.is_open())
	{
		cout << "Error: couldn't write model file." << endl;
		cout << "Make sure the following path can be written: " << file_path << endl;
		exit(1);
	}

	int i;
	for (i=0; i<strong_models.size(); i++)
		if (! strong_models[i].write_model(os))
		{
			cout << "Error: no data exists for " << config_->get_fragment(strong_models[i].model_frag_idx).label << endl;
			exit(1);
		}

	for (i=0; i<regular_models.size(); i++)
		if (! regular_models[i].write_model(os))
		{
			cout << "Warning: no data exists for " << config_->get_fragment(regular_models[i].model_frag_idx).label << endl;
			exit(1);
		}

	os.close();

	return true;
}


bool RegionalPrmNodeScoreModel::read_regional_score_model(const char *name, bool silent_ind)
{
	string file_path = make_model_file_name(name);
	ifstream is(file_path.c_str());

	if (! is.good() || ! is.is_open())
	{
		is.close();
		return false;
	}

	int i;
	for (i=0; i<strong_models.size(); i++)
	{
		strong_models[i].read_model(config_, is, silent_ind);
		if (! strong_models[i].ind_has_models)
			return false;
	}
	
	for (i=0; i<regular_models.size(); i++)
	{
		regular_models[i].read_model(config_, is, silent_ind);
	}

	is.close();
	
	was_initialized = true;

	return true;
}

void RegionalPrmNodeScoreModel::init(const Config* config, int charge, int sizeIndex, int regionIndex)
{
	config_ = config;
	charge_ = charge;
	sizeIndex_ = sizeIndex;
	regionIndex_ = regionIndex;

	const RegionalFragments& rf = config->get_regional_fragments(charge_, sizeIndex_, regionIndex_);

	frag_type_idxs = rf.get_frag_type_idxs();
	frag_probs     = rf.get_frag_probs();
	rand_prob	   = rf.get_rand_prob();
	log_random	   = log(rand_prob);
	log_one_minus_random = log(1.0-rand_prob);

	// init strong and regular models
	const vector<int>& strong_idxs = rf.get_strong_frag_type_idxs();
	const vector< vector<int> >& mirror_frag_idxs = rf.get_mirror_frag_idxs();

/*	cout << "CHARGE " << charge << " SIZE " << size_idx << " REGION " << region_idx << endl;
	cout << "frag_idxs: ";
	int f;
	for (f=0; f<frag_type_idxs.size(); f++)
		cout << "\t" << frag_type_idxs[f];
	cout << endl;
	cout << "strong_idxs: ";
	for (f=0; f<strong_idxs.size(); f++)
		cout << "\t" << strong_idxs[f];
	cout << endl << endl;*/
	
	strong_models.resize(strong_idxs.size());
	regular_models.resize(frag_type_idxs.size()-strong_idxs.size());

	num_strong_frags=0;
	num_regular_frags=0;
	int i;
	for (i=0; i<strong_idxs.size(); i++)
	{
		const int frag_idx = strong_idxs[i];

		StrongFragmentModel& strong_model=strong_models[num_strong_frags++];

		strong_model.model_frag_idx = frag_idx;
		strong_model.model_frag_charge = config->get_fragment(frag_idx).charge;
		strong_model.set_config_and_tolerance(config);

		strong_model.parent1_idx=-1;
		strong_model.parent2_idx=-1;
		strong_model.parent1_charge=-1;
		strong_model.parent2_charge=-1;

		// add parents
		int f;
		for (f=0; f<2 && f<frag_type_idxs.size(); f++)
		{
			if (frag_type_idxs[f] == frag_idx)
				break;

			const int parent_idx = frag_type_idxs[f];
			const int parent_charge = config->get_fragment(parent_idx).charge;
			if (f==0)
			{
				strong_model.parent1_idx = parent_idx;
				strong_model.parent1_charge = parent_charge;
			}
			else
			{
				strong_model.parent2_idx = parent_idx;
				strong_model.parent2_charge = parent_charge;
			}
		}

		// add mirrors
		int g;
		for (g=0; g<2 && g<mirror_frag_idxs[f].size(); g++)
		{
			const int mirror_idx = mirror_frag_idxs[f][g];
			const int mirror_charge = config->get_fragment(mirror_idx).charge;
			if (f==0)
			{
				strong_model.mirror1_idx = mirror_idx;
				strong_model.mirror1_charge = mirror_charge;
			}
			else
			{
				strong_model.mirror2_idx = mirror_idx;
				strong_model.mirror2_charge = mirror_charge;
			}
		}
	}

	// add regular models
	for (i=0; i<frag_type_idxs.size(); i++)
	{
		int j;
		for (j=0; j<strong_idxs.size(); j++)
			if (strong_idxs[j] == frag_type_idxs[i])
				break;

		if (j<strong_idxs.size())
			continue;

		const int reg_frag_idx = frag_type_idxs[i];

		RegularFragmentModel& regular_model=regular_models[num_regular_frags++];
		regular_model.model_frag_idx = reg_frag_idx;
		regular_model.model_frag_charge = config->get_fragment(reg_frag_idx).charge;
		regular_model.set_config_and_tolerance(config);
		regular_model.parent_idxs.clear();
		regular_model.parent_idx_with_same_charge_ori = -1;

		// add parents
		int f;
		for (f=0; f<frag_type_idxs.size() && f<8; f++)
		{
			if (frag_type_idxs[f] == reg_frag_idx)
				break;

			regular_model.parent_idxs.push_back(frag_type_idxs[f]);

			if (regular_model.parent_idx_with_same_charge_ori<0)
			{
				const FragmentType& parent_frag = config->get_fragment(frag_type_idxs[f]);
				if (parent_frag.charge == regular_model.model_frag_charge &&
					parent_frag.orientation == config->get_fragment(reg_frag_idx).orientation)
				{
					regular_model.parent_idx_with_same_charge_ori=frag_type_idxs[f];
				}
			}
		}
		regular_model.num_parents = regular_model.parent_idxs.size();
	}

	// set scores
	const int num_all_frags = config->get_all_fragments().size();
	frag_inten_scores.resize(num_all_frags,0);
	frag_no_inten_scores.resize(num_all_frags,0);
	missing_breakage_score=0;
	
	for (i=0; i<frag_type_idxs.size(); i++)
	{
		const int frag_idx=frag_type_idxs[i];
		const score_t missing_score = log(1.0 - frag_probs[i]) - log_one_minus_random;
		missing_breakage_score += missing_score;

		frag_no_inten_scores[frag_idx] = missing_score;
		frag_inten_scores[frag_idx] = log(frag_probs[i])- log_random;
	}

	//  init the weights 
//	vector<float> strong_inten_weights,  strong_no_inten_weights;
//	vector<float> regular_inten_weights, regular_no_inten_weights;

	const int num_strong = strong_models.size();
	strong_inten_weights.resize(num_strong,0.99);
	strong_no_inten_weights.resize(num_strong,0.98);

	strong_inten_danc_part.resize(num_strong);
	strong_no_inten_danc_part.resize(num_strong);

	for (i=0; i<num_strong; i++)
	{
		const float frag_prob = get_frag_prob(strong_models[i].model_frag_idx);
		strong_inten_danc_part[i]=(1.0-strong_inten_weights[i])*frag_prob;
		strong_no_inten_danc_part[i]=(1.0-strong_no_inten_weights[i])*(1.0-frag_prob);

//		if (size_idx == 1 && region_idx == 1)
//			cout << i << "\t" << frag_prob << "\t" << strong_inten_danc_part[i] << "\t" << strong_no_inten_danc_part[i] << endl;
	}


	// holds the (1-weight) * frag prob to be weight*model prob for the acutal probability
//	vector<float> strong_inten_danc_part,  strong_no_inten_danc_part;
//	vector<float> regular_inten_danc_part, regular_no_inten_danc_part;

	const int num_regular = regular_models.size();
	regular_inten_weights.resize(num_regular,0.97);
	regular_no_inten_weights.resize(num_regular,0.99);

	regular_inten_danc_part.resize(num_regular);
	regular_no_inten_danc_part.resize(num_regular);

	for (i=0; i<num_regular; i++)
	{
		regular_inten_weights[i] = 1.0 - 0.01*i;
		regular_no_inten_weights[i] = 0.5;
	}

	for (i=0; i<num_regular; i++)
	{
		const float frag_prob = get_frag_prob(regular_models[i].model_frag_idx);
		regular_inten_danc_part[i]=(1.0-regular_inten_weights[i])*frag_prob;
		regular_no_inten_danc_part[i]=(1.0-regular_no_inten_weights[i])*(1.0-frag_prob);

	//	if (size_idx == 1 && region_idx == 1)
	//		cout << i << "\t" << frag_prob << "\t" << regular_inten_danc_part[i] << "\t" << regular_no_inten_danc_part[i] << endl;
	}
}


bool RegionalPrmNodeScoreModel::trainRegionalScoreModel(void *modelPointer, 
														const char *name, 
														const SpectraAggregator& sa)
{
	const int min_num_of_samples_per_feature=6;
	const int num_me_rounds = 750;
	int i;

	AllScoreModels* model = static_cast<AllScoreModels*>(modelPointer);
	for (i=0; i<strong_models.size(); i++)
	{
		ME_Regression_DataSet inten_ds, no_inten_ds;
		const int frag_idx = strong_models[i].model_frag_idx;
		const float frag_prob = get_frag_prob(frag_idx);

		cout << endl << endl << "TRAINING INTENSITY MODEL FOR CHARGE " <<
			charge_ << " SIZE " << sizeIndex_ << " REGION " << regionIndex_ << " FRAGMENT " << i << " " <<
			config_->get_fragment(frag_idx).label << endl << endl;

		int j;
		for (j=0; j<3; j++)
		{
			cout << "SEED: " << getRandomSeed() << endl;
			// TODO fix the create trianing set function
			createTrainingSet(model, strong_models[i], sa, inten_ds, no_inten_ds);

			inten_ds.purge_low_count_features(min_num_of_samples_per_feature);
			int num_bad=inten_ds.check_samples(true);
			if (num_bad>0)
				cout << "Warning: had " << num_bad << " bad samples removed!" << endl;

			inten_ds.print_summary();
			inten_ds.print_feature_summary(cout, ScoreModelFields_SI_names);
			if (! strong_models[i].inten_model.train_cg(inten_ds,num_me_rounds,2E-5))
			{
				cout << "Coudln't train ME model, setting all weights to 0! (" <<j << ")" << endl;
			}
			else
			{
				cout << endl << "INTENSTY - Charge " << charge_ << " size " << sizeIndex_ << " region " << regionIndex_ << 
					" fragment " << frag_idx << " " << config_->get_fragment(frag_idx).label <<  endl;
				strong_models[i].inten_model.print_ds_probs(inten_ds);
				strong_models[i].inten_log_scaling_factor = 
					log(strong_models[i].inten_model.calc_log_scaling_constant(0,inten_ds,1.2*frag_prob));

				break;
			}
		}

		if (j==3)
		{
			cout << "Problem with models!!!" << endl;
			exit(1);
		}

		cout << endl << endl << "TRAINING NO INTENSITY MODEL FOR CHARGE " <<
			charge_ << " SIZE " << sizeIndex_ << " REGION " << regionIndex_ << " FRAGMENT " << i << " " <<
			config_->get_fragment(frag_idx).label << endl << endl;

		no_inten_ds.purge_low_count_features(min_num_of_samples_per_feature);
		int num_bad=no_inten_ds.check_samples(true);
			if (num_bad>0)
				cout << "Warning: had " << num_bad << " bad samples removed!" << endl;

		no_inten_ds.print_summary();
		no_inten_ds.print_feature_summary(cout, ScoreModelFields_SNI_names);
		if (! strong_models[i].no_inten_model.train_cg(no_inten_ds,num_me_rounds,2E-5))
		{
			cout << "Coudln't train no inten ME model, exiting!" << endl;
			exit(1);
		}
		else
		{
			cout << endl << "NO INTENSITY - Charge " << charge_ << " size " << sizeIndex_ << " region " << regionIndex_ << 
				" fragment " << config_->get_fragment(frag_idx).label << endl;
			strong_models[i].no_inten_model.print_ds_probs(no_inten_ds);
			strong_models[i].no_inten_log_scaling_factor = 
				log(strong_models[i].no_inten_model.calc_log_scaling_constant(0,no_inten_ds,1.1*(1.0-frag_prob)));
		}

		strong_models[i].ind_has_models = true;
	}

	for (i=0; i<regular_models.size(); i++)
	{
		ME_Regression_DataSet inten_ds, no_inten_ds;
		const int frag_idx = regular_models[i].model_frag_idx;
		const float frag_prob = get_frag_prob(frag_idx);

		cout << endl << endl << "TRAINING INTENSITY MODEL FOR CHARGE " <<
			charge_ << " SIZE " << sizeIndex_ << " REGION " << regionIndex_ << " FRAGMENT " << i << " " <<
			config_->get_fragment(frag_idx).label << endl << endl;

		int j;
		for (j=0; j<3; j++)
		{
			cout << "SEED: " << getRandomSeed() << endl;
			//TODO FIX create function
			createTrainingSet(model, regular_models[i], sa, inten_ds, no_inten_ds);

			inten_ds.purge_low_count_features(min_num_of_samples_per_feature);
			int num_bad=inten_ds.check_samples(true);
			if (num_bad>0)
				cout << "Warning: had " << num_bad << " bad samples removed!" << endl;


			inten_ds.print_summary();
			inten_ds.print_feature_summary(cout, ScoreModelFields_RI_names);
			if (! regular_models[i].inten_model.train_cg(inten_ds,num_me_rounds,2E-5))
			{
				cout << "Coudln't train ME model, setting all weights to 0! (" <<j<<")"<< endl;
			}
			else
			{
				cout << endl << "INTENSTY - Charge " << charge_ << " size " << sizeIndex_ << " region " << regionIndex_ << 
					" fragment " << config_->get_fragment(frag_idx).label << endl ;
				regular_models[i].inten_model.print_ds_probs(inten_ds);
				regular_models[i].inten_log_scaling_factor = 
					log(regular_models[i].inten_model.calc_log_scaling_constant(0,inten_ds,1.1*frag_prob));
				break;
			}
		}
		
		cout << endl << endl << "TRAINING NO INTENSITY MODEL FOR CHARGE " <<
				charge_ << " SIZE " << sizeIndex_ << " REGION " << regionIndex_ << " FRAGMENT " << i << " " <<
				config_->get_fragment(frag_idx).label << endl << endl;

		no_inten_ds.purge_low_count_features(min_num_of_samples_per_feature);

		int num_bad=no_inten_ds.check_samples(true);
		if (num_bad>0)
			cout << "Warning: had " << num_bad << " bad samples removed!" << endl;


		no_inten_ds.print_summary();
		no_inten_ds.print_feature_summary(cout, ScoreModelFields_RNI_names);
		if (! regular_models[i].no_inten_model.train_cg(no_inten_ds,num_me_rounds,2E-5))
		{
			cout << "Coudln't train ME model, setting all weights to 0!" << endl;
		}
		else
		{
			cout << endl << "NO INTENSTY - Charge " << charge_ << " size " << sizeIndex_ << " region " << regionIndex_ << 
				" fragment " << config_->get_fragment(frag_idx).label << endl;
				regular_models[i].no_inten_model.print_ds_probs(no_inten_ds);
				regular_models[i].no_inten_log_scaling_factor = 
					log(regular_models[i].no_inten_model.calc_log_scaling_constant(0,no_inten_ds,(1.0-frag_prob)));
		}

		regular_models[i].ind_has_models = true;
	}	


	was_initialized = true;

	write_regional_score_model(name);

	return true;
}

/*
bool RegionalPrmNodeScoreModel::train_regional_score_model(void *model_ptr, const char *name, const FileManager& fm)
{
	const int min_num_of_samples_per_feature=6;
	const int num_me_rounds = 750;
	int i;

	AllScoreModels* model = static_cast<AllScoreModels*>(model_ptr);
	for (i=0; i<strong_models.size(); i++)
	{
		ME_Regression_DataSet inten_ds, no_inten_ds;
		const int frag_idx = strong_models[i].model_frag_idx;
		const float frag_prob = get_frag_prob(frag_idx);

		cout << endl << endl << "TRAINING INTENSITY MODEL FOR CHARGE " <<
			charge_ << " SIZE " << sizeIndex_ << " REGION " << regionIndex_ << " FRAGMENT " << i << " " <<
			config_->get_fragment(frag_idx).label << endl << endl;

		int j;
		for (j=0; j<3; j++)
		{
			cout << "SEED: " << getRandomSeed() << endl;
			// TODO fix the create trianing set function
		//	create_training_set(model, strong_models[i], fm, inten_ds, no_inten_ds);

			inten_ds.purge_low_count_features(min_num_of_samples_per_feature);
			int num_bad=inten_ds.check_samples(true);
			if (num_bad>0)
				cout << "Warning: had " << num_bad << " bad samples removed!" << endl;

			inten_ds.print_summary();
			inten_ds.print_feature_summary(cout, ScoreModelFields_SI_names);
			if (! strong_models[i].inten_model.train_cg(inten_ds,num_me_rounds,2E-5))
			{
				cout << "Coudln't train ME model, setting all weights to 0! (" <<j << ")" << endl;
			}
			else
			{
				cout << endl << "INTENSTY - Charge " << charge_ << " size " << sizeIndex_ << " region " << regionIndex_ << 
					" fragment " << frag_idx << " " << config_->get_fragment(frag_idx).label <<  endl;
				strong_models[i].inten_model.print_ds_probs(inten_ds);
				strong_models[i].inten_log_scaling_factor = 
					log(strong_models[i].inten_model.calc_log_scaling_constant(0,inten_ds,1.2*frag_prob));

				break;
			}
		}

		if (j==3)
		{
			cout << "Problem with models!!!" << endl;
			exit(1);
		}

		cout << endl << endl << "TRAINING NO INTENSITY MODEL FOR CHARGE " <<
			charge_ << " SIZE " << sizeIndex_ << " REGION " << regionIndex_ << " FRAGMENT " << i << " " <<
			config_->get_fragment(frag_idx).label << endl << endl;

		no_inten_ds.purge_low_count_features(min_num_of_samples_per_feature);
		int num_bad=no_inten_ds.check_samples(true);
			if (num_bad>0)
				cout << "Warning: had " << num_bad << " bad samples removed!" << endl;

		no_inten_ds.print_summary();
		no_inten_ds.print_feature_summary(cout, ScoreModelFields_SNI_names);
		if (! strong_models[i].no_inten_model.train_cg(no_inten_ds,num_me_rounds,2E-5))
		{
			cout << "Coudln't train no inten ME model, exiting!" << endl;
			exit(1);
		}
		else
		{
			cout << endl << "NO INTENSITY - Charge " << charge_ << " size " << sizeIndex_ << " region " << regionIndex_ << 
				" fragment " << config_->get_fragment(frag_idx).label << endl;
			strong_models[i].no_inten_model.print_ds_probs(no_inten_ds);
			strong_models[i].no_inten_log_scaling_factor = 
				log(strong_models[i].no_inten_model.calc_log_scaling_constant(0,no_inten_ds,1.1*(1.0-frag_prob)));
		}

		strong_models[i].ind_has_models = true;
	}

	for (i=0; i<regular_models.size(); i++)
	{
		ME_Regression_DataSet inten_ds, no_inten_ds;
		const int frag_idx = regular_models[i].model_frag_idx;
		const float frag_prob = get_frag_prob(frag_idx);

		cout << endl << endl << "TRAINING INTENSITY MODEL FOR CHARGE " <<
			charge_ << " SIZE " << sizeIndex_ << " REGION " << regionIndex_ << " FRAGMENT " << i << " " <<
			config_->get_fragment(frag_idx).label << endl << endl;

		int j;
		for (j=0; j<3; j++)
		{
			cout << "SEED: " << getRandomSeed() << endl;
			//TODO FIX create function
		//	create_training_set(model, regular_models[i], fm, inten_ds, no_inten_ds);

			inten_ds.purge_low_count_features(min_num_of_samples_per_feature);
			int num_bad=inten_ds.check_samples(true);
			if (num_bad>0)
				cout << "Warning: had " << num_bad << " bad samples removed!" << endl;


			inten_ds.print_summary();
			inten_ds.print_feature_summary(cout, ScoreModelFields_RI_names);
			if (! regular_models[i].inten_model.train_cg(inten_ds,num_me_rounds,2E-5))
			{
				cout << "Coudln't train ME model, setting all weights to 0! (" <<j<<")"<< endl;
			}
			else
			{
				cout << endl << "INTENSTY - Charge " << charge_ << " size " << sizeIndex_ << " region " << regionIndex_ << 
					" fragment " << config_->get_fragment(frag_idx).label << endl ;
				regular_models[i].inten_model.print_ds_probs(inten_ds);
				regular_models[i].inten_log_scaling_factor = 
					log(regular_models[i].inten_model.calc_log_scaling_constant(0,inten_ds,1.1*frag_prob));
				break;
			}
		}
		
		cout << endl << endl << "TRAINING NO INTENSITY MODEL FOR CHARGE " <<
				charge_ << " SIZE " << sizeIndex_ << " REGION " << regionIndex_ << " FRAGMENT " << i << " " <<
				config_->get_fragment(frag_idx).label << endl << endl;

		no_inten_ds.purge_low_count_features(min_num_of_samples_per_feature);

		int num_bad=no_inten_ds.check_samples(true);
		if (num_bad>0)
			cout << "Warning: had " << num_bad << " bad samples removed!" << endl;


		no_inten_ds.print_summary();
		no_inten_ds.print_feature_summary(cout, ScoreModelFields_RNI_names);
		if (! regular_models[i].no_inten_model.train_cg(no_inten_ds,num_me_rounds,2E-5))
		{
			cout << "Coudln't train ME model, setting all weights to 0!" << endl;
		}
		else
		{
			cout << endl << "NO INTENSTY - Charge " << charge_ << " size " << sizeIndex_ << " region " << regionIndex_ << 
				" fragment " << config_->get_fragment(frag_idx).label << endl;
				regular_models[i].no_inten_model.print_ds_probs(no_inten_ds);
				regular_models[i].no_inten_log_scaling_factor = 
					log(regular_models[i].no_inten_model.calc_log_scaling_constant(0,no_inten_ds,(1.0-frag_prob)));
		}

		regular_models[i].ind_has_models = true;
	}	


	was_initialized = true;

	write_regional_score_model(name);

	return true;
}
*/


