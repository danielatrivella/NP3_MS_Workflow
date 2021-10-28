#include "AllScoreModels.h"
#include "PepNovo_auxfun.h"
#include "PrmGraph.h"
#include "SpectraList.h"







/***************************************************************************
****************************************************************************/
void RegionalPrmNodeScoreModel::createTrainingSet(AllScoreModels *model,
							 const FragmentModel& frag_model,
							 const SpectraAggregator& sa,
							 ME_Regression_DataSet& inten_ds,
							 ME_Regression_DataSet& no_inten_ds) const
{
	const mass_t min_mass = config_->get_min_mass_for_size_idx(charge_, sizeIndex_);
	const mass_t max_mass = config_->get_max_mass_for_size_idx(charge_, sizeIndex_);
	const mass_t tolerance_diff = config_->getTolerance()*0.5;
	const int frag_idx = frag_model.get_model_frag_idx();

	bool is_strong=false;
	int s;
	for (s=0; s<strong_models.size(); s++)
		if (frag_idx == strong_models[s].model_frag_idx)
			is_strong=true;

	SpectraList sl(sa);
	sl.selectHeaders((min_mass/charge_)-2.0, (max_mass/charge_)-2.0, charge_, charge_);
	
	cout << "Selected " << sl.getNumHeaders() << " headers..." << endl;
	cout << "Min m/z " << min_mass/charge_ << endl;
	cout << "Max m/z " << max_mass/charge_ << endl;

	PMCSQS_Scorer *pmcsqs = (PMCSQS_Scorer *)model->get_pmcsqs_ptr();
	const bool use_pmc = (pmcsqs && pmcsqs->getIndInitializedPmc());

	if (use_pmc)
		cout << "Using parent mass correction model to set PM.." << endl;

	int i;
	for (i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
		if (header->getPeptideStr().length() <= 0)
			continue;

		Peptide peptide;
		peptide.parseFromString(config_, header->getPeptideStr());

		const mass_t true_mass_with_19 = peptide.get_mass_with_19();
		Spectrum s;

		s.readSpectrum(sa, header);
	
		// calc corrected pm_with_19, if it is good use it, otherwise, use a value with +- U[0,toleance/2]
		mass_t pm_with_19=NEG_INF;
	
		if (use_pmc)
		{
			mass_t	mz1,mz2;
			int		charge1,charge2;
			float	prob1,prob2;
					

			// output m/z and prob values for the different charge states

			pmcsqs->computeBestMzAndChargeForSpectrum(config_, s, &mz1,&charge1,&prob1,&mz2,&charge2,&prob2);

			const mass_t corr1_pm_with_19 = mz1*charge1 - MASS_PROTON*(charge1-1);
			const mass_t corr2_pm_with_19 = mz2*charge2 - MASS_PROTON*(charge2-1);
			if (fabs(corr2_pm_with_19-true_mass_with_19)<tolerance_diff)
				pm_with_19 = corr2_pm_with_19;
			if (fabs(corr1_pm_with_19-true_mass_with_19)<tolerance_diff)
				pm_with_19 = corr1_pm_with_19;
		}

		if (pm_with_19<0) // use a random value
		{
			double r=myRandom();
			mass_t offset = r*r*tolerance_diff;
			if (myRandom()<0.5)
				offset *= -1;
			pm_with_19 = true_mass_with_19 + offset;
	
		}

		int spec_size_idx = config_->calc_size_idx(charge_, pm_with_19);
		if (spec_size_idx != sizeIndex_)
			continue;
	
		PrmGraph prm;
		prm.create_graph_from_spectrum(model,&s,pm_with_19, charge_);


		vector<BreakageInfo> good_peak_examples, bad_peak_examples;
		prm.extract_breakage_infos_for_score_training(model,
													  frag_model.get_model_frag_idx(), 
													  regionIndex_,
													  is_strong,
													  good_peak_examples, 
													  bad_peak_examples);

	//	cout << "GOOD:" << endl;
		int j;
		for (j=0; j<good_peak_examples.size(); j++)
		{
		//	good_peak_examples[j].print(config);
			ME_Regression_Sample sam;
			sam.label=0;

			frag_model.fill_single_frag_vector(&s, pm_with_19, good_peak_examples[j].breakage,
				good_peak_examples[j], sam.f_vals);

			if (good_peak_examples[j].breakage->get_position_of_frag_idx(frag_idx)>=0)
			{
				inten_ds.add_sample(sam);
			//	sam.print((is_strong ? ScoreModelFields_SI_names : ScoreModelFields_RI_names));
			}
			else
				no_inten_ds.add_sample(sam);
		}
		
	//	cout << "BAD:" << endl;
		for (j=0; j<bad_peak_examples.size(); j++)
		{
		//	bad_peak_examples[j].print(config);
			ME_Regression_Sample sam;
			sam.label=1;

			frag_model.fill_single_frag_vector(&s, pm_with_19, bad_peak_examples[j].breakage,
				bad_peak_examples[j], sam.f_vals);

			if (bad_peak_examples[j].breakage->get_position_of_frag_idx(frag_idx)>=0)
			{
				inten_ds.add_sample(sam);
			}
			else
				no_inten_ds.add_sample(sam);

	//		sam.print();
		}

		if (i>0 && i %1000 == 0)
		{
			cout << i << "/" << sl.getNumHeaders() << " ..." << endl;
		}
	}

	
	const score_t frag_prob  = this->get_frag_prob(frag_idx);
	cout << "Probability of observeing fragment: " << frag_prob << endl;

	inten_ds.num_classes=2;
	inten_ds.num_features = (is_strong ? (int)SI_NUM_FEATURES : (int)RI_NUM_FEATURES);
	inten_ds.calibrate_class_weights((is_strong ? 0.45 : 0.2));
	
	no_inten_ds.num_classes=2;
	no_inten_ds.num_features = (is_strong ? (int)SNI_NUM_FEATURES : (int)RNI_NUM_FEATURES);
	no_inten_ds.calibrate_class_weights((is_strong ? 0.2 : 0.5));

	if (is_strong)
	{
		vector<int> inten_features;
		inten_features.push_back(SI_IND_CONNECTS_TO_N_TERM);
		inten_features.push_back(SI_IND_N_IS_GAP);
		inten_features.push_back(SI_IND_C_IS_GAP);

		inten_ds.serial_scale(inten_features);

		vector<int> no_inten_features;
		no_inten_features.push_back(SNI_IND_CONNECTS_TO_N_TERM);
		no_inten_features.push_back(SNI_IND_N_IS_GAP);
		no_inten_features.push_back(SNI_IND_C_IS_GAP);

		no_inten_ds.serial_scale(no_inten_features);
	}
	else
	{
		vector<int> inten_features;
		inten_features.push_back(RI_IND_N_IS_GAP);
		inten_features.push_back(RI_IND_C_IS_GAP);

		inten_ds.serial_scale(inten_features);

		vector<int> no_inten_features;
		no_inten_features.push_back(RNI_IND_N_IS_GAP);
		no_inten_features.push_back(RNI_IND_C_IS_GAP);

		no_inten_ds.serial_scale(no_inten_features);
	}
}






void BreakageInfo::print(Config *config) const
{
	const vector<string>& aa2label = config->get_aa2label();
	cout << setw(2) << this->type << "> ";
	cout << (this->connects_to_N_term ? '[' : ' ');
	cout << (this->n_edge_is_single ? '1' : '2');
	cout <<  " " << aa2label[this->n_aa] << " : " << aa2label[this->c_aa] << " ";
	cout << (this->c_edge_is_single ? '1' : '2');
	cout << (this->connects_to_C_term ? ']' : ' ');
	cout << " " << fixed << setprecision(3) << "  # " <<this->node_idx << " ";
	if (this->missed_cleavage)
		cout << "MC!";
	cout << endl;
}

void PrmGraph::fill_breakage_info(const AllScoreModels *model, BreakageInfo *info, int node_idx,
								  int n_edge_idx, int n_variant_idx,
								  int c_edge_idx, int c_variant_idx, int type) const
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const vector<int>& n_term_digest_aas = config->get_n_term_digest_aas();
	const vector<int>& c_term_digest_aas = config->get_c_term_digest_aas();
	
	info->node_idx = node_idx;
	info->breakage = &nodes[node_idx].breakage;
	info->type = type;

	info->n_edge_idx = n_edge_idx;
	info->n_var_idx  = n_variant_idx;
	info->c_edge_idx = c_edge_idx;
	info->c_var_idx  = c_variant_idx;

	int n_second_before_cut = NEG_INF;
	if (n_edge_idx>=0)
	{
		const MultiEdge& n_edge = multi_edges[n_edge_idx];
		const Node& n_node = nodes[n_edge.n_idx];
		int *v_ptr = n_edge.variant_ptrs[n_variant_idx];

		info->n_var_ptr = v_ptr;

		const int num_aa = *v_ptr++;
		int n_aa = v_ptr[num_aa-1];
		mass_t exp_edge_mass=0;
		int j;
		for (j=0; j<num_aa; j++)
			exp_edge_mass+=aa2mass[v_ptr[j]];

		if (n_aa == Ile || n_aa == Xle)
			n_aa = Leu;

		info->n_aa = n_aa;
		info->nn_aa = v_ptr[0];
		
		if (info->nn_aa == Ile || info->nn_aa == Xle)
			info->nn_aa = Leu;

		info->n_break = n_edge.n_break;
		info->exp_n_edge_mass = exp_edge_mass;
		info->n_edge_is_single = (num_aa ==1);

		if (! info->n_edge_is_single)
			n_second_before_cut = v_ptr[num_aa-2];
		if (n_second_before_cut == Ile || n_second_before_cut == Xle)
			n_second_before_cut = Leu;

		if (n_node.type == NODE_N_TERM)
		{
			info->connects_to_N_term=true;
			int j;
			for (j=0; j<n_term_digest_aas.size(); j++)
				if (n_term_digest_aas[j] == *v_ptr)
				{
					info->preferred_digest_aa_N_term=true;
					break;
				}
		}

		for (j=0; j<c_term_digest_aas.size(); j++)
			if (c_term_digest_aas[j] == v_ptr[num_aa-1])
			{
				info->missed_cleavage= true;
				break;
			}

		info->ind_n_edge_overlaps = (n_edge.ind_edge_overlaps);
		info->n_side_cat = model->get_aa_category(num_aa,v_ptr,info->connects_to_N_term,false);
	}
	else
	{
		info->n_aa=Gap;
		info->nn_aa=Gap;
	}

	if (c_edge_idx>=0)
	{
		const MultiEdge& c_edge = multi_edges[c_edge_idx];
		const Node& c_node = nodes[c_edge.c_idx];
		int *v_ptr = c_edge.variant_ptrs[c_variant_idx];

		info->c_var_ptr = v_ptr;

		const int num_aa = *v_ptr++;
		int c_aa = v_ptr[0];
		mass_t exp_edge_mass=0;
		int j;
		for (j=0; j<num_aa; j++)
			exp_edge_mass+=aa2mass[v_ptr[j]];

		if (c_aa == Ile || c_aa == Xle)
			c_aa = Leu;

		info->c_aa = c_aa;
		info->cc_aa = v_ptr[num_aa-1];

		if (info->cc_aa == Ile || info->cc_aa == Xle)
			info->cc_aa = Leu;

		info->c_break = c_edge.c_break;
		info->exp_c_edge_mass = exp_edge_mass;
		info->c_edge_is_single = (num_aa ==1);

		int c_second_after_cut=NEG_INF;
		if (! info->c_edge_is_single)
			c_second_after_cut = v_ptr[1];
		if (c_second_after_cut == Ile || c_second_after_cut == Xle)
			c_second_after_cut = Leu;

		if (c_node.type == NODE_C_TERM)
		{
			info->connects_to_C_term=true;
			int j;
			for (j=0; j<c_term_digest_aas.size(); j++)
				if (c_term_digest_aas[j] == v_ptr[num_aa-1])
				{
					info->preferred_digest_aa_C_term=true;
					break;
				}
		}

		for (j=0; j<n_term_digest_aas.size(); j++)
			if (n_term_digest_aas[j] == v_ptr[0])
			{
				info->missed_cleavage= true;
				break;
			}

		info->ind_c_edge_overlaps = (c_edge.ind_edge_overlaps);
		info->c_side_cat = model->get_aa_category(num_aa,v_ptr,false,info->connects_to_C_term);
		if (n_edge_idx>=0)
		{
			int aas[2]={info->n_aa,info->c_aa};
			info->span_cat = model->get_aa_category(2,aas,info->connects_to_N_term && info->n_edge_is_single,info->connects_to_C_term && info->c_edge_is_single);
			
			if (! info->n_edge_is_single)
			{
				if (n_second_before_cut<0)
				{
					cout << "Error: bad aa 2nd before n-side!" << endl;
					exit(1);
				}
				int aas[3]={n_second_before_cut,info->n_aa,info->c_aa};
				info->n_double_span_cat = model->get_aa_category(3,aas,info->connects_to_N_term ,info->connects_to_C_term && info->c_edge_is_single);
			}

			if (! info->c_edge_is_single)
			{
				if (c_second_after_cut<0)
				{
					cout << "Error: bad aa 2nd after c-side!" << endl;
					exit(1);
				}
				int aas[3]={info->n_aa,info->c_aa,c_second_after_cut};
				info->c_double_span_cat = model->get_aa_category(3,aas,info->connects_to_N_term && info->n_edge_is_single, info->connects_to_C_term);
			}
		}

	}
	else
	{
		info->c_aa=Gap;
		info->cc_aa=Gap;
	}
}

/************************************************************************************************
*************************************************************************************************/
void PrmGraph::extract_breakage_infos_for_score_training(AllScoreModels *model,
														 int frag_idx,
														 int target_region_idx,
														 bool ind_strong_frag,
														 vector<BreakageInfo>& good_examples,
														 vector<BreakageInfo>& bad_examples) const
{

	const double Gap_ratio = 0.05;
	const mass_t tolerance = config->getTolerance();

	const Peptide& true_pep = source_spectrum->getPeptide();
	const vector<int>& aas = true_pep.get_amino_acids();
	vector<mass_t> correct_break_masses;
	vector<int>    node_to_breakages, correct_node_idxs;
	vector<int>    correct_edge_variant_map;
	
	// find minimal and maximal nodes for which the frag is visible
	const mass_t min_mass = source_spectrum->get_min_peak_mass()-1.0;
	const mass_t max_mass = source_spectrum->get_max_peak_mass()+1.0;
	const FragmentType& frag = config->get_fragment(frag_idx);
	int min_viz_idx=nodes.size()+1;
	int max_viz_idx=0;
	int i;
	for (i=0; i<nodes.size(); i++)
	{
		mass_t exp_mass = frag.calc_expected_mass(nodes[i].mass,pm_with_19);
		if (exp_mass>min_mass && exp_mass<max_mass)
		{
			if (i<min_viz_idx)
				min_viz_idx=i;
			if (i>max_viz_idx)
				max_viz_idx=i;
		}
	}

	good_examples.clear();
	bad_examples.clear();

	true_pep.calc_expected_breakage_masses(config,correct_break_masses);
	node_to_breakages.resize(nodes.size(),NEG_INF);
	correct_node_idxs.clear();
	correct_edge_variant_map.resize(multi_edges.size(),NEG_INF);

	
	for (i=0; i<correct_break_masses.size(); i++)
	{
		const int max_node_idx = get_max_score_node(correct_break_masses[i],tolerance);
		if (max_node_idx>=0)
		{
			node_to_breakages[max_node_idx]=i;
			correct_node_idxs.push_back(max_node_idx);
		}
	}

//	cout << endl <<true_pep.as_string(config) << endl;
	for (i=0; i<multi_edges.size(); i++)
	{
		const MultiEdge& edge = multi_edges[i];
		if (node_to_breakages[edge.n_idx]>=0 && node_to_breakages[edge.c_idx]>=0)
		{
			const int brekage_idx = node_to_breakages[edge.n_idx];
			correct_edge_variant_map[i]=multi_edges[i].get_variant_idx(edge.num_aa,&aas[brekage_idx]);

			if (0 && correct_edge_variant_map[i]>=0)
			{
				cout << brekage_idx << "\t" << i << "\t" << correct_edge_variant_map[i] << "\t";
				int *var_ptr = edge.variant_ptrs[correct_edge_variant_map[i]];
				int num_aa = *var_ptr++;
				int j;
				for (j=0; j<num_aa; j++)
					cout << config->get_aa2label()[*var_ptr++];
				cout << endl;
			}
		}
	}

	vector<int> correct_in_edge_idxs;
	vector<int> correct_out_edge_idxs;
	correct_in_edge_idxs.resize(correct_node_idxs.size(),NEG_INF);
	correct_out_edge_idxs.resize(correct_node_idxs.size(),NEG_INF);

	bool had_good_connect_to_n_term = false;
	bool had_good_connect_to_c_term = false;

	// create infos for good peak samples
	for (i=0; i<correct_node_idxs.size(); i++)
	{
		const int node_idx = correct_node_idxs[i];
		const int node_region_idx = nodes[node_idx].breakage.region_idx;
		if (node_region_idx != target_region_idx || node_idx==0 || 
			node_idx==nodes.size()-1 || node_idx<min_viz_idx || node_idx>max_viz_idx)
			continue;

		const Node& node = nodes[node_idx];
		
		int n_edge_idx = NEG_INF;
		int c_edge_idx = NEG_INF;
		int n_varaint_idx = NEG_INF;
		int c_variant_idx = NEG_INF;

		int j;
		for (j=0; j<node.in_edge_idxs.size(); j++)
			if (correct_edge_variant_map[node.in_edge_idxs[j]]>=0)
			{
				n_edge_idx=node.in_edge_idxs[j];
				n_varaint_idx=correct_edge_variant_map[node.in_edge_idxs[j]];
				correct_in_edge_idxs[i]=n_edge_idx;
				break;
			}

		for (j=0; j<node.out_edge_idxs.size(); j++)
			if (correct_edge_variant_map[node.out_edge_idxs[j]]>=0)
			{
				c_edge_idx=node.out_edge_idxs[j];
				c_variant_idx=correct_edge_variant_map[node.out_edge_idxs[j]];
				correct_out_edge_idxs[i]=c_edge_idx;
				break;
			}

		BreakageInfo info;
		fill_breakage_info(model,&info,node_idx,n_edge_idx,n_varaint_idx,c_edge_idx,c_variant_idx,1);

		if (info.connects_to_N_term)
			had_good_connect_to_n_term=true;
		
		if (info.connects_to_C_term)
			had_good_connect_to_c_term=true;

		good_examples.push_back(info);
		if (n_edge_idx>=0 && c_edge_idx>=0 && myRandom()<Gap_ratio*0.33)
		{
			BreakageInfo gap_info;
			fill_breakage_info(model,&gap_info,node_idx,NEG_INF,NEG_INF,c_edge_idx,c_variant_idx,11);
			good_examples.push_back(gap_info);
		}

		if (n_edge_idx>=0 && c_edge_idx>=0 && myRandom()<Gap_ratio*0.33)
		{
			BreakageInfo gap_info;
			fill_breakage_info(model,&gap_info,node_idx,n_edge_idx,n_varaint_idx,NEG_INF,NEG_INF,11);
			good_examples.push_back(gap_info);
		}
	}

	// create for bad peak samples
	const int num_good = good_examples.size();

	// same as good nodes, but using bad edges
	// perform only for strong fragments!
	const int num_half_bad_examples = (ind_strong_frag ? int(num_good*0.333) : 0);
	for (i=0; i<num_half_bad_examples; i++)
	{
		const int node_idx = correct_node_idxs[i];
		const int node_region_idx = nodes[node_idx].breakage.region_idx;

		if (node_region_idx != target_region_idx || node_idx==0 || 
			node_idx==nodes.size()-1 || node_idx<min_viz_idx || node_idx>max_viz_idx)
			continue;

		const Node& node = nodes[node_idx];
		if (node.breakage.get_position_of_frag_idx(frag_idx)<0) // don't use such samples if there is no peak
			continue;

		if (node.in_edge_idxs.size()>1 && correct_in_edge_idxs[i]>=0 &&
			node.out_edge_idxs.size()>1 && correct_out_edge_idxs[i]>=0)
		{
			int j;
			vector<int> in_idxs,out_idxs;
			for (j=0; j<node.in_edge_idxs.size(); j++)
			{
				const int edge_idx = node.in_edge_idxs[j];
				if (correct_edge_variant_map[edge_idx]>=0)
					continue;
				if (multi_edges[edge_idx].num_aa == 1)
				{
					in_idxs.push_back(edge_idx);
				}
				else if (multi_edges[edge_idx].num_aa != 1 && (in_idxs.size()<3 || myRandom()<0.1))
					in_idxs.push_back(edge_idx);
			}
			
			for (j=0; j<node.out_edge_idxs.size(); j++)
			{
				const int edge_idx = node.out_edge_idxs[j];
				if (correct_edge_variant_map[edge_idx]>=0)
					continue;
				if (multi_edges[edge_idx].num_aa == 1)
				{
					out_idxs.push_back(edge_idx);
				}
				else if (multi_edges[edge_idx].num_aa != 1 && (out_idxs.size()<3 || myRandom()<0.1))
					out_idxs.push_back(edge_idx);
			}

			if (in_idxs.size()==0 || out_idxs.size()==0)
				continue;

			const int bad_n_edge = in_idxs[int(myRandom()*in_idxs.size())];
			const int bad_c_edge = out_idxs[int(myRandom()*out_idxs.size())];
			const int n_var_idx = int(multi_edges[bad_n_edge].variant_ptrs.size() * myRandom());
			const int c_var_idx = int(multi_edges[bad_c_edge].variant_ptrs.size() * myRandom());

			BreakageInfo info;
			fill_breakage_info(model,&info,node_idx,bad_n_edge,n_var_idx,bad_c_edge,c_var_idx,2);
			bad_examples.push_back(info);
			if (bad_n_edge>=0 && bad_c_edge>=0 && myRandom()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,NEG_INF,NEG_INF,bad_c_edge,c_var_idx,22);
				bad_examples.push_back(gap_info);
			}

			if (bad_n_edge>=0 && bad_c_edge>=0 && myRandom()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,bad_n_edge,n_var_idx,NEG_INF,NEG_INF,22);
				bad_examples.push_back(gap_info);
			}
		}
	}


	// mostly single edges from the highest scoring incorrect nodes
	// gives advantage for aa combinations with high scoring nodes
	vector<score_pair> pairs;
	for (i=1; i<nodes.size()-1; i++)
		if (node_to_breakages[i]<0 && i>=min_viz_idx && i<=max_viz_idx &&
			nodes[i].breakage.region_idx==target_region_idx )
				pairs.push_back(score_pair(i,nodes[i].score));
	
	sort(pairs.begin(),pairs.end());

	vector<size_t> selected_idxs, random_idxs;
	int num_to_select = int(num_good * 3);
	for (i=0; i<num_to_select && i<pairs.size(); i++)
		selected_idxs.push_back(pairs[i].idx);

	if (pairs.size()< num_to_select)
		num_to_select = pairs.size();

	chooseKFromN(num_to_select,pairs.size(),random_idxs);

	for (i=0; i<random_idxs.size(); i++)
		selected_idxs.push_back(pairs[random_idxs[i]].idx);

	for (i=0; i<selected_idxs.size(); i++)
	{
		const int node_idx = selected_idxs[i];
		const Node& node = nodes[node_idx];
		vector<int> in_idxs,out_idxs;
		score_t max_n_score=NEG_INF;
		score_t max_c_score=NEG_INF;
		int best_n_idx=-1;
		int best_c_idx=-1;
		if (node.in_edge_idxs.size()>0 && node.out_edge_idxs.size()>0)
		{
			int j;
			for (j=0; j<node.in_edge_idxs.size(); j++)
			{
				const int edge_idx = node.in_edge_idxs[j];
				if (correct_edge_variant_map[edge_idx]>=0)
					continue;

				if (multi_edges[edge_idx].num_aa == 1)
				{
					in_idxs.push_back(edge_idx);
					if (nodes[multi_edges[edge_idx].n_idx].score>max_n_score)
					{
						max_n_score=nodes[multi_edges[edge_idx].n_idx].score;
						best_n_idx=edge_idx;
					}
				}
				else if (multi_edges[edge_idx].num_aa != 1 && (in_idxs.size()<2 || myRandom()<0.1))
					in_idxs.push_back(edge_idx);
			}

			if (in_idxs.size()==0)
				continue;
			
			for (j=0; j<node.out_edge_idxs.size(); j++)
			{
				const int edge_idx = node.out_edge_idxs[j];
				if (correct_edge_variant_map[edge_idx]>=0)
					continue;

				if (multi_edges[edge_idx].num_aa == 1)
				{
					out_idxs.push_back(edge_idx);
					if (nodes[multi_edges[edge_idx].c_idx].score>max_c_score)
					{
						max_c_score=nodes[multi_edges[edge_idx].c_idx].score;
						best_c_idx=edge_idx;
					}
				}
				else if (multi_edges[edge_idx].num_aa != 1 && (out_idxs.size()<2 || myRandom()<0.1))
					out_idxs.push_back(edge_idx);
			}
		}
		if (in_idxs.size()==0 || out_idxs.size()==0)
			continue;

		// add entries with the highest scoring nodes to bias the results
		if (best_n_idx>=0)
		{
			int num_to_add = int(0.3 * in_idxs.size()+0.5);
			int k;
			for (k=0; k<num_to_add; k++)
			{
				in_idxs.push_back(best_n_idx);
				in_idxs.push_back(best_n_idx);
			}
		}

		if (best_c_idx>=0)
		{
			int num_to_add = int(0.3 * out_idxs.size()+0.5);
			int k;
			for (k=0; k<num_to_add; k++)
			{
				out_idxs.push_back(best_c_idx);
				out_idxs.push_back(best_c_idx);
			}
		}

		const int bad_n_edge = in_idxs[int(myRandom()*in_idxs.size())];
		const int bad_c_edge = out_idxs[int(myRandom()*out_idxs.size())];
		const int n_var_idx = int(multi_edges[bad_n_edge].variant_ptrs.size() * myRandom());
		const int c_var_idx = int(multi_edges[bad_c_edge].variant_ptrs.size() * myRandom());

		BreakageInfo info;
		fill_breakage_info(model,&info,node_idx,bad_n_edge,n_var_idx,bad_c_edge,c_var_idx,4);
		bad_examples.push_back(info);

		if (myRandom()<Gap_ratio)
		{
			BreakageInfo gap_info;
			fill_breakage_info(model,&gap_info,node_idx,bad_n_edge,n_var_idx,NEG_INF,NEG_INF,44);
			bad_examples.push_back(gap_info);
		}

		if ( myRandom()<Gap_ratio)
		{
			BreakageInfo gap_info;
			fill_breakage_info(model,&gap_info,node_idx,NEG_INF,NEG_INF,bad_c_edge,c_var_idx,44);
			bad_examples.push_back(gap_info);
		}
	}

	// add a few nodes connecting to the N-terminal
	if (had_good_connect_to_n_term && nodes[0].out_edge_idxs.size()>2)
	{
		vector<int> in_idxs;
		int j;
		for (j=0; j<nodes[0].out_edge_idxs.size(); j++)
		{
			const int edge_idx = nodes[0].out_edge_idxs[j];
			if (correct_edge_variant_map[edge_idx]>=0)
				continue;

			if (multi_edges[edge_idx].num_aa == 1)
			{
				in_idxs.push_back(edge_idx);
			}
			else if (multi_edges[edge_idx].num_aa != 1 && (in_idxs.size()<3 || myRandom()<0.1))
				in_idxs.push_back(edge_idx);
		}

		vector<int> selected_in_idxs;
		if (in_idxs.size()>3)
		{
			vector<size_t> positions;
			chooseKFromN(3,in_idxs.size(),positions);
			int j;
			for (j=0; j<3; j++)
				selected_in_idxs.push_back(in_idxs[positions[j]]);
		}
		else
			selected_in_idxs = in_idxs;
		
		int k;
		for (k=0; k<selected_in_idxs.size(); k++)
		{
			const int bad_n_edge = selected_in_idxs[k];
			const MultiEdge& in_edge = multi_edges[bad_n_edge];
			const int node_idx = in_edge.c_idx;
			const Node& node = nodes[node_idx];

			vector<int> out_idxs;
			int j;
			for (j=0; j<node.out_edge_idxs.size(); j++)
			{
				const int edge_idx = node.out_edge_idxs[j];
				if (correct_edge_variant_map[edge_idx]>=0)
					continue;

				if (multi_edges[edge_idx].num_aa == 1)
				{
					out_idxs.push_back(edge_idx);
				}
				else if (multi_edges[edge_idx].num_aa != 1 && (out_idxs.size()<3 || myRandom()<0.1))
					out_idxs.push_back(edge_idx);
			}

			if (out_idxs.size() == 0)
				continue;

			const int bad_c_edge = out_idxs[int(myRandom()*out_idxs.size())];

			const int n_var_idx = int(multi_edges[bad_n_edge].variant_ptrs.size() * myRandom());
			const int c_var_idx = int(multi_edges[bad_c_edge].variant_ptrs.size() * myRandom());

			BreakageInfo info;
			fill_breakage_info(model,&info,node_idx,bad_n_edge,n_var_idx,bad_c_edge,c_var_idx,5);
			bad_examples.push_back(info);

			if (myRandom()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,bad_n_edge,n_var_idx,NEG_INF,NEG_INF,55);
				bad_examples.push_back(gap_info);
			}

			if ( myRandom()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,NEG_INF,NEG_INF,bad_c_edge,c_var_idx,55);
				bad_examples.push_back(gap_info);
			}
		}
	}

	// add an example that connects to the C-terminal
	if (had_good_connect_to_c_term && nodes[nodes.size()-1].in_edge_idxs.size()>1)
	{
		vector<int> out_idxs;
		int j;
		for (j=0; j<nodes[nodes.size()-1].in_edge_idxs.size(); j++)
		{
			const int edge_idx = nodes[nodes.size()-1].in_edge_idxs[j];
			if (correct_edge_variant_map[edge_idx]>=0)
				continue;

			if (multi_edges[edge_idx].num_aa == 1)
			{
				out_idxs.push_back(edge_idx);
			}
			else if (multi_edges[edge_idx].num_aa != 1 && myRandom()<0.025)
				out_idxs.push_back(edge_idx);
		}

		if (out_idxs.size() == 0)
			return;

		const int bad_c_edge = out_idxs[int(myRandom()*out_idxs.size())];
		
		const MultiEdge& out_edge = multi_edges[bad_c_edge];
		const int node_idx = out_edge.n_idx;
		const Node& node = nodes[node_idx];

		vector<int> in_idxs;
		for (j=0; j<node.in_edge_idxs.size(); j++)
		{
			const int edge_idx = node.in_edge_idxs[j];
			if (correct_edge_variant_map[edge_idx]>=0)
				continue;

			if (multi_edges[edge_idx].num_aa == 1)
			{
				in_idxs.push_back(edge_idx);
			}
			else if (multi_edges[edge_idx].num_aa != 1 && (in_idxs.size()<3 || myRandom()<0.1))
				in_idxs.push_back(edge_idx);
		}

		if (in_idxs.size()>0)
		{
			const int bad_n_edge = in_idxs[int(myRandom()*in_idxs.size())];

			const int n_var_idx = int(multi_edges[bad_n_edge].variant_ptrs.size() * myRandom());
			const int c_var_idx = int(multi_edges[bad_c_edge].variant_ptrs.size() * myRandom());

			BreakageInfo info;
			fill_breakage_info(model,&info,node_idx,bad_n_edge,n_var_idx,bad_c_edge,c_var_idx,6);
			bad_examples.push_back(info);

			if (myRandom()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,bad_n_edge,n_var_idx,NEG_INF,NEG_INF,66);
				bad_examples.push_back(gap_info);
			}

			if ( myRandom()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,NEG_INF,NEG_INF,bad_c_edge,c_var_idx,66);
				bad_examples.push_back(gap_info);
			}
		}
		
	}
}

