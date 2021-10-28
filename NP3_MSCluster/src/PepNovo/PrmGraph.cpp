#include "PrmGraph.h"
#include "AnnotatedSpectrum.h"
#include "PepNovo_auxfun.h"


void PrmGraph::clear()
{
	nodes.clear(); 
	multi_edges.clear(); 
	nodeIndexArray_.clear(); 
	has_node_combo_scores=false;
	out_aa_ind.clear(); 
	in_aa_ind.clear(); 
	forbidden_node_idxs.clear(); 
	cummulative_scores.clear();
	n_digest_node_idxs.clear(); 
	c_digest_node_idxs.clear(); 
	source_spectrum=NULL; 
	config=NULL; model=NULL;
}

bool comp_SeqPath_sort_key(const SeqPath& a, const SeqPath& b)
{
	return (a.sort_key>b.sort_key);
}

bool comp_SeqPath_path_score (const SeqPath& a, const SeqPath& b)
{
	return (a.path_score>b.path_score);
}



/********************************************************
*********************************************************/
void PrmGraph::create_graph_from_spectrum(AllScoreModels *_model, 
										  Spectrum *spectrum,
										  mass_t _pm_with_19, 
										  int spec_charge, 
										  bool add_all_pepitde_nodes, 
										  bool only_basic_score)
{
	if (nodes.size()>0)
		this->clear();

	model  = _model;
	config = model->get_config();
	source_spectrum = spectrum;
	pm_with_19 = _pm_with_19;
	charge = spec_charge;

	has_node_combo_scores=false;

	if (charge==0)
		charge = source_spectrum->getCharge();

	int org_size_idx = spectrum->get_size_idx();
	size_idx = config->calc_size_idx(charge,pm_with_19);

	spectrum->set_size_idx(size_idx);

	digest_node_score = config->get_digest_score();
	model->init_model_for_scoring_spectrum(source_spectrum);

//	config->print_session_aas();

	if (pm_with_19<10)
	{
		cout << "Error: supplied negative/low PM for graph!" << endl;
		exit(1);
	}

	create_nodes();
	add_digest_nodes();
	score_nodes(model);

	fill_single_multi_edges();
	fill_double_multi_edges();

	int l;
	for (l=3; l<=config->get_max_edge_length(); l++)
		fill_longer_multi_edges(l);

	if (! only_basic_score)
		model->getPrmNodeScoreModel().initial_combos_score(model,this);

	prune_low_scoring_nodes();

	rank_nodes_according_to_score();

	set_idxs_max_in_out_for_nodes();

	fill_forbidden_idx_indicators_and_cumulative_scores();

	// restore size idx
	spectrum->set_size_idx(org_size_idx);
}

/***********************************************************************
Creates a graph tailored for a peptide.
Contains only nodes and edges that correspond to that peptide.
************************************************************************/
void PrmGraph::create_graph_for_peptide_and_spectrum(AllScoreModels *_model, Spectrum *spectrum, 
					mass_t _pm_with_19, int spec_charge, const Peptide& peptide)
{
	if (nodes.size()>0)
		this->clear();

	model  = _model;
	config = model->get_config();
	source_spectrum = spectrum;
	pm_with_19 = _pm_with_19;
	charge = spec_charge;

	has_node_combo_scores=false;

	if (charge==0)
		charge = source_spectrum->getCharge();

	int org_size_idx = spectrum->get_size_idx();
	size_idx = config->calc_size_idx(charge,pm_with_19);

	spectrum->set_size_idx(size_idx);

	digest_node_score = config->get_digest_score();
	model->init_model_for_scoring_spectrum(source_spectrum);

	if (pm_with_19<10)
	{
		cout << "Error: supplied negative/low PM for graph!" << endl;
		exit(1);
	}

	create_nodes_for_peptide(peptide);

	const vector<int>& pep_aas = peptide.get_amino_acids();
	mass_t n_digest_mass = NEG_INF;
	mass_t c_digest_mass = NEG_INF;
	
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const vector<int>& n_digest_aas = config->get_n_term_digest_aas();
	if (n_digest_aas.size()>0)
	{
		int j;
		for (j=0; j<n_digest_aas.size(); j++)
			if (pep_aas[0] == n_digest_aas[j])
				break;
		if (j<n_digest_aas.size())
			n_digest_mass = aa2mass[n_digest_aas[j]];
	}

	const vector<int>& c_digest_aas = config->get_c_term_digest_aas();
	if (c_digest_aas.size()>0)
	{
		int j;
		for (j=0; j<c_digest_aas.size(); j++)
			if (pep_aas[pep_aas.size()-1] == c_digest_aas[j])
				break;
		if (j<c_digest_aas.size())
			c_digest_mass = peptide.get_mass() - aa2mass[c_digest_aas[j]];
	}

	if (n_digest_mass>0 || c_digest_mass>0)
		add_digest_nodes(n_digest_mass,c_digest_mass);

	score_nodes(model);

	fill_single_multi_edges();
	fill_double_multi_edges();

	int l;
	for (l=3; l<=config->get_max_edge_length(); l++)
		fill_longer_multi_edges(l);

	model->getPrmNodeScoreModel().score_peptide_node_combos(static_cast<void*>(model), this, peptide);

	rank_nodes_according_to_score();

	set_idxs_max_in_out_for_nodes();

	// restore size idx
	spectrum->set_size_idx(org_size_idx);

}





/*********************************************************************
Initializes the index array.
For each rounded off Dalton m, it gives the index of the closest peak i
with mass m_i>= m.
**********************************************************************/
void PrmGraph::init_index_array()
{
	int i,c,size=(int)pm_with_19+2;
	const int max_node_idx = nodes.size()-1;
	
	nodeIndexArray_.clear();
	nodeIndexArray_.resize(size,-1);
	
	i=0;
	int m=(int)nodes[0].mass;
	while (i<m)
		nodeIndexArray_[i++]=0;

	c=0;
	while (c< max_node_idx)
	{
		int curr_m=(int)nodes[c].mass;
		int next_m = curr_m;
		int next_c = c;

		while (next_m == curr_m && next_c<max_node_idx)
			next_m=(int)nodes[++next_c].mass;

		while (i<next_m)
			nodeIndexArray_[i++]=c;
		
		c=next_c;
	}

	while (i<size)
		nodeIndexArray_[i++]=max_node_idx;
}



/***********************************************************
Merges nodes that are close to each other
Gives preference to the prefix fragment, but merge genereally
goes according to the intensities of the source fragments
************************************************************/
void PrmGraph::merge_close_nodes()
{
	int i;
	const vector<FragmentType>& all_fragments = config->get_all_fragments();

	mass_t delta = config->getTolerance();
	if (delta>0.2)
		delta *= 0.8;

	for (i=0; i<nodes.size()-1; i++)
	{
		if (nodes[i+1].mass - nodes[i].mass < delta &&
			nodes[i+1].source_frag_type_idx != nodes[i].source_frag_type_idx)
		{
			const int frag1_pos = nodes[i].breakage.get_position_of_frag_idx(nodes[i].source_frag_type_idx);
			const int frag2_pos = nodes[i+1].breakage.get_position_of_frag_idx(nodes[i+1].source_frag_type_idx);

			if (frag1_pos <0)
			{
				nodes[i].mass=999999;
				continue;
			}

			if (frag2_pos<0)
			{
				nodes[i+1].mass=9999999;
				continue;
			}
			if (frag1_pos <0 || frag2_pos<0)
			{
				Node& node1 = nodes[i];
				Node& node2 = nodes[i+1];


				cout << i << "\t";
				node1.breakage.print_fragments(config);
				cout << endl;
				cout << i+1 << "\t";
				node2.breakage.print_fragments(config);
				cout << endl;
				cout << "Error: source fragments not found in breakage!: " <<
					config->get_fragment(node1.source_frag_type_idx).label << "  " << 
					config->get_fragment(node2.source_frag_type_idx).label << endl;

				print();
				exit(1);
			}
			mass_t new_node_mass;
			intensity_t inten1 = nodes[i].breakage.fragments[frag1_pos].intensity;
			intensity_t inten2 = nodes[i+1].breakage.fragments[frag2_pos].intensity;
			if (all_fragments[nodes[i].source_frag_type_idx].orientation == PREFIX)
			{
				inten1 *= 2;
			}
			else if (all_fragments[nodes[i+1].source_frag_type_idx].orientation == PREFIX)
			{
				inten2 *= 2;
			}

			mass_t mass_times_inten1 = nodes[i].mass * inten1;
			mass_t mass_times_inten2 = nodes[i+1].mass * inten2;

			new_node_mass = (mass_times_inten1 + mass_times_inten2)/ (inten1 + inten2);

			// create new node, move it to the i+1 position
			nodes[i].mass = 99999999;
			nodes[i+1].mass = new_node_mass;

			// transfer fragments from node i to i+1 if they are not there
			int f;
			for (f=0; f<nodes[i].breakage.fragments.size(); f++)
			{
				if (nodes[i+1].breakage.get_position_of_frag_idx(
					 nodes[i].breakage.fragments[f].frag_type_idx) < 0)
				{
//					cout << "Adding fragments! " << endl;
//					nodes[i].breakage.print();
//					nodes[i+1].breakage.print();
//					cout<<endl;

					nodes[i+1].breakage.add_fragment(nodes[i].breakage.fragments[f]);
					
				}
			}
		}
	}
	
	sort(nodes.begin(),nodes.end());
	while (nodes.size()>0 && nodes[nodes.size()-1].mass > 50000)
		nodes.pop_back();

	init_index_array();

}



struct frag_region_list {
	int frag_type_idx;
	vector<int> region_idxs;
};


/*******************************************************
Selectes nodes for PrmGraph.
Selection done in two stages. First every strong peak is
considered, then combos are considered.
********************************************************/
void PrmGraph::create_nodes()
{
	const int num_regions = config->get_num_regions(charge,size_idx);
	const vector<FragmentType>& all_fragments = config->get_all_fragments();
	const Peak* peaks  = source_spectrum->getPeaks();
	const int numPeaks = source_spectrum->getNumPeaks();
	const vector<int>& strong_peak_idxs = source_spectrum->get_strong_peak_idxs();
	const mass_t max_mass_to_create_node = pm_with_19 - 55;
	const mass_t mid_mass = pm_with_19 * 0.5;
	const mass_t min_peak_mass = source_spectrum->get_min_peak_mass();
	const mass_t max_peak_mass = source_spectrum->get_max_peak_mass();

	nodes.clear();
	nodes.resize(1);

	vector< vector<int> > peak_usages; // holds for each peak all the interpertations
									   // that were given to it for creating nodes

	peak_usages.resize(numPeaks);
	
	// add N-TERM
	nodes[0].mass=0;
	nodes[0].type = NODE_N_TERM;
	nodes[0].breakage.mass = 0;
	nodes[0].breakage.region_idx=0;
	nodes[0].breakage.parent_charge=charge;
	nodes[0].breakage.parent_size_idx = size_idx;


	// create list for each strong frag type, what regions it can be used as 
	// a basis for creating a node
	vector<frag_region_list> strong_frags_lists;
	int r,i;
	for (r=0; r<num_regions; r++)
	{
		const RegionalFragments& rf = config->get_regional_fragments(charge,size_idx,r);
		const vector<int>& strong_frag_type_idxs = rf.get_strong_frag_type_idxs();
		int f;
		for (f=0; f<strong_frag_type_idxs.size(); f++)
		{
			int j;
			for (j=0; j<strong_frags_lists.size(); j++)
			{
				if (strong_frags_lists[j].frag_type_idx == strong_frag_type_idxs[f])
				{
					strong_frags_lists[j].region_idxs[r]=1;
					break;
				}
			}

			if (j==strong_frags_lists.size())
			{
				frag_region_list frl;
				frl.frag_type_idx= strong_frag_type_idxs[f];
				frl.region_idxs.resize(num_regions,0);
				frl.region_idxs[r]=1;
				strong_frags_lists.push_back(frl);
			}
		}
	}
	


	// create nodes from strong peaks

	for (i=0; i<strong_frags_lists.size(); i++)
	{
		const int strong_frag_idx = strong_frags_lists[i].frag_type_idx;
		const FragmentType& ft = all_fragments[strong_frag_idx];
		const vector<int>& permitted_regions = strong_frags_lists[i].region_idxs;

		int q;
		for (q=0; q<strong_peak_idxs.size(); q++)
		{
			const int p_idx = strong_peak_idxs[q];
			const mass_t exp_break_mass = ft.calc_breakage_mass(peaks[p_idx].mass,pm_with_19);

			if (exp_break_mass<mid_mass && ! config->is_allowed_prefix_mass(exp_break_mass))
			{
			//	cout << "Bad prefix mass: " << exp_break_mass << endl;
				continue;
			}

			if (exp_break_mass>mid_mass && ! config->is_allowed_suffix_mass(pm_with_19,exp_break_mass))
			{
			//	cout << "Bad Suffix mass: " << exp_break_mass << endl;
				continue;
			} 

			if (exp_break_mass < 3 || exp_break_mass > max_mass_to_create_node)
				continue;

			const int region_idx = config->calc_region_idx(exp_break_mass,pm_with_19, 
				charge, min_peak_mass, max_peak_mass);

			const RegionalFragments & rf = config->get_regional_fragments(charge,size_idx,region_idx);
			const vector<FragmentCombo>&  combos = rf.get_frag_type_combos();

			if (! permitted_regions[region_idx])
				continue;

			Node node;
			node.breakage.region_idx= region_idx;
			node.mass = exp_break_mass;
			node.breakage.mass = exp_break_mass;

			node.type = NODE_REG;
			node.source_frag_type_idx = strong_frag_idx;

			// make sure fragment is present
			BreakageFragment brf;
			brf.expected_mass = peaks[p_idx].mass;
			brf.frag_type_idx = strong_frag_idx;
			brf.mass = peaks[p_idx].mass;
			brf.intensity = peaks[p_idx].intensity;
			brf.peak_idx = p_idx;

			node.breakage.add_fragment(brf);

			annotate_breakage(source_spectrum, pm_with_19, size_idx, node.breakage);
			model->getPrmNodeScoreModel().score_breakage(source_spectrum,&node.breakage);
			nodes.push_back(node);

			// add peak usage
			peak_usages[p_idx].push_back(strong_frag_idx);
		}
	}

	// Create nodes from combos
	// if peaks were already used for a certain fragment, then don't create a node
	// the first idx in the combo is a strong_frag_type_idx

	// first create a list for the first strong types from the combo
	vector<int> strong_frags_in_combos;
	for (r=0; r<num_regions; r++)
	{
		const vector<FragmentCombo>& combos = config->get_regional_fragments(charge,size_idx,r).get_frag_type_combos();
		int c;
		for (c=0; c<combos.size(); c++)
		{
			int j;
			for (j=0; j<strong_frags_in_combos.size(); j++)
				if (strong_frags_in_combos[j]== combos[c].frag_inten_idxs[0])
					break;
			if (j<strong_frags_in_combos.size())
				continue;
			strong_frags_in_combos.push_back(combos[c].frag_inten_idxs[0]);
		}
	}

	for (i=0; i<numPeaks; i++)
	{
		if (peak_usages[i].size()>0 || source_spectrum->get_peak_iso_level(i)>0)
			continue;

		int f;
		for (f=0; f<strong_frags_in_combos.size(); f++)
		{
			const int strong_frag_idx = strong_frags_in_combos[f]; 
			const FragmentType& strong_frag = all_fragments[strong_frag_idx];

			const mass_t exp_breakage_mass = strong_frag.calc_breakage_mass(peaks[i].mass,pm_with_19);

			if (exp_breakage_mass<mid_mass && ! config->is_allowed_prefix_mass(exp_breakage_mass))
				continue;
	
			if (exp_breakage_mass>mid_mass && ! config->is_allowed_suffix_mass(pm_with_19,exp_breakage_mass))
				continue;
		
			if (exp_breakage_mass < 3 || exp_breakage_mass > max_mass_to_create_node)
				continue;

			const int region_idx = config->calc_region_idx(exp_breakage_mass, pm_with_19, charge, 
															min_peak_mass, max_peak_mass);
			const vector<FragmentCombo>& combos = config->get_regional_fragments(charge,size_idx,region_idx).get_frag_type_combos();
			int c;
			for (c=0; c<combos.size(); c++)
			{
				const FragmentCombo& combo = combos[c];
				if (combo.frag_inten_idxs[0] == strong_frag_idx)
				{
					Node node;
					node.breakage.region_idx= region_idx;
					node.mass = exp_breakage_mass;
					node.breakage.mass = exp_breakage_mass;

					node.type = NODE_REG;
					node.source_frag_type_idx = strong_frag_idx;

					annotate_breakage(source_spectrum, pm_with_19, size_idx, node.breakage);
					// check that all frags are present
					int j;
					for (j=0; j<combo.frag_inten_idxs.size(); j++)
						if (node.breakage.get_position_of_frag_idx(combo.frag_inten_idxs[j])<0)
							break;
					
					if (j<combo.frag_inten_idxs.size())
						continue;

					model->getPrmNodeScoreModel().score_breakage(source_spectrum,&node.breakage);

					nodes.push_back(node);

					// add peak usage
					peak_usages[i].push_back(strong_frag_idx);
					break;
				}
			}
		}
	}

	if (nodes.size()>5)
	{
		sort(nodes.begin(),nodes.end());
		init_index_array();

		const int num_nodes_to_add = 8 +int(nodes[nodes.size()-1].mass * 0.00125);
		vector<Node> cand_nodes;

		// add nodes for peaks interperted as strong fragment types that have an amino acid before or an amino acid 
		// after them with offset tolerance/2
		const mass_t tolerance = config->getTolerance();
		const mass_t third_tolerance = (tolerance<0.1 ? tolerance*0.666 : tolerance * 0.333);
		const vector<int>& session_aas = config->get_session_aas();
		const vector<mass_t>& aa2mass = config->get_aa2mass();

		for (i=0; i<numPeaks; i++)
		{
			if (source_spectrum->get_peak_iso_level(i)>0 || peak_usages[i].size()>0)
				continue;

			const mass_t peak_mass = peaks[i].mass;
			int f;
			for (f=0; f<strong_frags_lists.size(); f++)
			{
				const int strong_frag_idx = strong_frags_lists[f].frag_type_idx;
				const FragmentType& strong_frag = all_fragments[strong_frag_idx];
				const mass_t one_over_frag_charge = 1.0/strong_frag.charge;

				const mass_t exp_breakage_mass = strong_frag.calc_breakage_mass(peaks[i].mass,pm_with_19);

				if (exp_breakage_mass < 3 || exp_breakage_mass > max_mass_to_create_node)
					continue;

				if (exp_breakage_mass<mid_mass && ! config->is_allowed_prefix_mass(exp_breakage_mass))
					continue;
		
				if (exp_breakage_mass>mid_mass && ! config->is_allowed_suffix_mass(pm_with_19,exp_breakage_mass))
					continue;

				if (get_max_score_node(exp_breakage_mass+MASS_NH3*one_over_frag_charge,third_tolerance)>0)
					continue;

				if (get_max_score_node(exp_breakage_mass+MASS_H2O*one_over_frag_charge,third_tolerance)>0)
					continue;
			
				// look for N-side connection
				score_t n_score=NEG_INF;
				if (1)
				{
					int a;
					
					for (a=0; a<session_aas.size(); a++)
					{
						const int idx = this->get_max_score_node(exp_breakage_mass - aa2mass[session_aas[a]],third_tolerance);
						if (idx>=0 && nodes[idx].breakage.score>n_score)
							n_score = nodes[idx].breakage.score;
					}
					if (n_score== NEG_INF)
						continue;
				}
				
				// look for C-side connection with same frag
				score_t c_score=NEG_INF;
				if (1)
				{
					int a;
					
					for (a=0; a<session_aas.size(); a++)
					{
						const int idx = this->get_max_score_node(exp_breakage_mass + aa2mass[session_aas[a]],third_tolerance);
						if (idx>=0 && nodes[idx].breakage.score>c_score)
							c_score = nodes[idx].breakage.score;
					}
					if (c_score== NEG_INF)
						continue;
				}

				const int region_idx = config->calc_region_idx(exp_breakage_mass, pm_with_19, charge, 
																min_peak_mass, max_peak_mass);
			
				Node node;
				node.breakage.region_idx= region_idx;
				node.mass = exp_breakage_mass;
				node.breakage.mass = exp_breakage_mass;

				node.type = NODE_REG;
				node.source_frag_type_idx = strong_frag_idx;

				annotate_breakage(source_spectrum, pm_with_19, size_idx, node.breakage);
				model->getPrmNodeScoreModel().score_breakage(source_spectrum,&node.breakage);

				node.tmp_score = n_score + c_score;

				int min_idx=NEG_INF;
				score_t min_score = POS_INF;
				if (cand_nodes.size()<num_nodes_to_add)
				{
					cand_nodes.push_back(node);
				}
				else
				{
					if (min_idx<0)
					{
						int j;
						for (j=0; j<cand_nodes.size(); j++)
							if (cand_nodes[j].tmp_score<min_score)
							{
								min_idx=j;
								min_score = cand_nodes[j].tmp_score;
							}
					}

					if (node.tmp_score>min_score)
					{
						cand_nodes[min_idx]=node;
						score_t min_score = POS_INF;
						int j;
						for (j=0; j<cand_nodes.size(); j++)
							if (cand_nodes[j].tmp_score<min_score)
							{
								min_idx=j;
								min_score = cand_nodes[j].tmp_score;
							}

					}
				}
			}
		}
		// add nodes
		for (i=0; i<cand_nodes.size(); i++)
		{
			nodes.push_back(cand_nodes[i]);
		//	cout << "Created :\t" << cand_nodes[i].mass << "\t" << cand_nodes[i].tmp_score << endl;
		}
	} 

	Node node;
	node.mass = pm_with_19 - MASS_OHHH;
	node.type = NODE_C_TERM;
	node.breakage.mass = node.mass;
	node.breakage.parent_charge=charge;
	node.breakage.parent_size_idx = size_idx;
	node.breakage.region_idx = config->calc_region_idx(node.mass,pm_with_19,charge,
														min_peak_mass,max_peak_mass);

	nodes.push_back(node);

	sort(nodes.begin(),nodes.end());
	
	merge_close_nodes();

	for (i=0; i<nodes.size(); i++)
	{
		nodes[i].breakage.parent_charge = source_spectrum->getCharge();
		nodes[i].breakage.parent_size_idx = source_spectrum->get_size_idx();
	}
	
}

/**********************************************************
Creates and annotates nodes that correspond to a pepitdes breakages
Only nodes that have some strong fragments associated with them
are kept.
***********************************************************/
void PrmGraph::create_nodes_for_peptide(const Peptide& pep)
{
	const vector< vector< vector< RegionalFragments > > >& all_rfs = config->get_regional_fragment_sets();
	const int num_regions = config->get_num_regions(charge,size_idx);
	const vector<FragmentType>& all_fragments = config->get_all_fragments();
	const Peak* const peaks  = source_spectrum->getPeaks();
	const int numPeaks = source_spectrum->getNumPeaks();
	const vector<int>& strong_peak_idxs = source_spectrum->get_strong_peak_idxs();
	const mass_t max_mass_to_create_node = pm_with_19 - 55;
	const mass_t mid_mass = pm_with_19 * 0.5;
	const mass_t min_peak_mass = source_spectrum->get_min_peak_mass();
	const mass_t max_peak_mass = source_spectrum->get_max_peak_mass();

	nodes.clear();
	nodes.resize(1);
	
	// add N-TERM
	nodes[0].mass=0;
	nodes[0].type = NODE_N_TERM;
	nodes[0].breakage.mass = 0;
	nodes[0].breakage.region_idx=0;
	nodes[0].breakage.parent_charge=charge;
	nodes[0].breakage.parent_size_idx = size_idx;

	vector<mass_t> exp_break_masses;
	
	pep.calc_expected_breakage_masses(config,exp_break_masses);
	
	// create nodes for peptide breakages
	int i;
	for (i=1; i<exp_break_masses.size()-1; i++)
	{
		const mass_t exp_break_mass = exp_break_masses[i];

		const int region_idx = config->calc_region_idx(exp_break_mass,pm_with_19,
												charge, min_peak_mass, max_peak_mass);

	
		Node node;
		node.breakage.region_idx= region_idx;
		node.mass = exp_break_mass;
		node.breakage.mass = exp_break_mass;

		node.type = NODE_REG;
		node.source_frag_type_idx = -1;

		annotate_breakage(source_spectrum, pm_with_19, size_idx, node.breakage);


		if (node.breakage.fragments.size() == 0)
			continue;
		
		const vector<int>& strong_frag_idxs = all_rfs[charge][size_idx][region_idx].get_strong_frag_type_idxs();
		int j;
		for (j=0; j<strong_frag_idxs.size(); j++)
			if (node.breakage.get_position_of_frag_idx(strong_frag_idxs[j])>=0)
				break;
		if (j==strong_frag_idxs.size())
			continue;
	
		nodes.push_back(node);
	}

	Node node;
	node.mass = pm_with_19 - MASS_OHHH;
	node.type = NODE_C_TERM;
	node.breakage.mass = node.mass;
	node.breakage.parent_charge=charge;
	node.breakage.parent_size_idx = size_idx;
	node.breakage.region_idx = config->calc_region_idx(node.mass,pm_with_19,charge,
														min_peak_mass,max_peak_mass);

	nodes.push_back(node);


	for (i=0; i<nodes.size(); i++)
	{
		nodes[i].breakage.parent_charge = source_spectrum->getCharge();
		nodes[i].breakage.parent_size_idx = source_spectrum->get_size_idx();
	}
}


/**********************************************************
Adds the digest nodes for the digest amino acids
***********************************************************/
void PrmGraph::add_digest_nodes(mass_t n_digest_mass, mass_t c_digest_mass)
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const vector<int>& n_term_digest_aas = config->get_n_term_digest_aas();
	const vector<int>& c_term_digest_aas = config->get_c_term_digest_aas();
	const mass_t digest_tolerance = 1.5 * config->getTolerance();
	const mass_t min_peak_mass = source_spectrum->get_min_peak_mass();
	const mass_t max_peak_mass = source_spectrum->get_max_peak_mass();

	bool added_nodes = false;

	if (c_term_digest_aas.size()>0)
	{
		int c_term_node=nodes.size();
		while (--c_term_node>=0)
			if (nodes[c_term_node].type == NODE_C_TERM)
				break;

		if (c_term_node<=0)
		{
			cout << "Error: couldn't find regular C-terminal!"<< endl;
			exit(1);
		}
		

		int i;
		for (i=0; i<c_term_digest_aas.size(); i++)
		{
			const int aa = c_term_digest_aas[i];
			mass_t exp_node_mass = nodes[c_term_node].mass - aa2mass[aa];

			int n=c_term_node-1;
			while (n>=0 && nodes[n].mass>exp_node_mass)
				n--;

			const int min_dis_node = fabs(nodes[n].mass - exp_node_mass) < 
				fabs(nodes[n+1].mass - exp_node_mass) ? n : n+1;


			if (c_digest_mass>0 && fabs(c_digest_mass - exp_node_mass)>1.5)
				continue;

			if (fabs(nodes[min_dis_node].mass - exp_node_mass)<digest_tolerance)
			{
				nodes[min_dis_node].type = NODE_DIGEST;
			}
			else // create the new node
			{
				Node node;
				node.mass = exp_node_mass;
				node.type = NODE_DIGEST;
				node.breakage.mass = node.mass;
				node.breakage.parent_charge=charge;
				node.breakage.parent_size_idx = size_idx;
				node.breakage.region_idx = config->calc_region_idx(node.mass,pm_with_19,charge,
																	min_peak_mass,max_peak_mass);

				nodes.push_back(node);
				added_nodes=true;
			}
		}
	}


	if (n_term_digest_aas.size()>0)
	{
		int n_term_node=0;
	
		int i;
		for (i=0; i<n_term_digest_aas.size(); i++)
		{
			const int aa = n_term_digest_aas[i];
			mass_t exp_node_mass =aa2mass[aa];

			int n=1;
			while (n<nodes.size() && nodes[n].mass<exp_node_mass)
				n++;

			const int min_dis_node = fabs(nodes[n].mass - exp_node_mass) < 
				fabs(nodes[n-1].mass - exp_node_mass) ? n : n-1;

			if (n_digest_mass>0 && fabs(n_digest_mass - nodes[min_dis_node].mass)>1.5)
				continue;

			if (fabs(nodes[min_dis_node].mass - exp_node_mass)<digest_tolerance)
			{
				nodes[min_dis_node].type = NODE_DIGEST;
			}
			else // create the new node
			{
				Node node;
				node.mass = exp_node_mass;
				node.type = NODE_DIGEST;
				node.breakage.mass = node.mass;
				node.breakage.parent_charge=charge;
				node.breakage.parent_size_idx = size_idx;
				node.breakage.region_idx = config->calc_region_idx(node.mass,pm_with_19, charge, min_peak_mass,max_peak_mass);

				nodes.push_back(node);
				added_nodes=true;
			}
		}
	}


	
	if (added_nodes)
	{
		sort(nodes.begin(),nodes.end());
		init_index_array(); // redo the index because digest nodes were added
	}

	n_digest_node_idxs.clear();
	c_digest_node_idxs.clear();

	int i;
	for (i=0; i<nodes.size(); i++)
	{
		if (nodes[i].type != NODE_DIGEST)
			continue;
		if (nodes[i].mass < 200.0)
		{
			n_digest_node_idxs.push_back(i);
		}
		else
			c_digest_node_idxs.push_back(i);
	}
}


/**********************************************************
Does initial scoring of nodes (uses only a node's breakage)
***********************************************************/
void PrmGraph::score_nodes(AllScoreModels *model)
{
	int i;

	const int max_node = nodes.size();
	for (i=0; i<max_node; i++)
	{
		bool verbose = false;
	
		if (nodes[i].type == NODE_REG || nodes[i].type == NODE_DIGEST) 
		{
			model->getPrmNodeScoreModel().score_breakage(source_spectrum,&nodes[i].breakage,verbose);

			nodes[i].score = nodes[i].breakage.score;
		}
		else if (nodes[i].type == NODE_N_TERM || nodes[i].type == NODE_C_TERM)
		{
			nodes[i].score = config->get_terminal_score();
			nodes[i].breakage.score = config->get_terminal_score();
		}
	}
}




/********************************************************************
	Fills multi edges.
	Connects between nodes with peaks (non-terminal and non digest).
*********************************************************************/
void PrmGraph::fill_single_multi_edges()
{
	const vector<AA_combo>& aa_edge_combos    = config->get_aa_edge_combos();
	const vector<int>& single_edge_combo_idxs = config->get_combo_idxs_by_length(1);
	const vector<int>& strong_fragment_idxs = config->get_all_strong_fragment_type_idxs();
	const vector<int>& session_aas = config->get_session_aas();
	const mass_t tolerance = config->getTolerance();
	const mass_t pm_tolerance = (config->get_pm_tolerance() > tolerance) ? config->get_pm_tolerance() : tolerance * 0.75;
	const mass_t digest_tolerance = 1.5 * tolerance;
	const mass_t half_tolerance = 0.5 * tolerance;
	const int num_nodes  = nodes.size();
	const int num_combos = single_edge_combo_idxs.size();
	const mass_t max_combos_mass = aa_edge_combos[single_edge_combo_idxs[single_edge_combo_idxs.size()-1]].total_mass;;
	const mass_t max_node_mass = nodes[num_nodes-1].mass;
	const int num_session_aas = config->get_session_aas().size();
	const int max_aa = config->get_session_aas()[num_session_aas-1];

	// check we need to allocate the single edge indicator arrays

	in_aa_ind.resize(max_aa+1);
	out_aa_ind.resize(max_aa+1);

	int i;
	for (i=0; i<session_aas.size(); i++)
	{
		const int aa = session_aas[i];
		in_aa_ind[aa].resize(num_nodes,false);
		out_aa_ind[aa].resize(num_nodes,false); 
	}


	// fill single aa edges
	for (i=0; i<num_nodes; i++)
	{
		if (nodes[i].type == NODE_C_TERM)
			continue;

		const mass_t& node_mass = nodes[i].mass;
		
		int last_node_idx=i;			
		int current_e_idx = -1;
		int current_c_idx = -1;
		
		int c;
		for (c=0; c<num_combos; c++)
		{
			const int& combo_idx = single_edge_combo_idxs[c];
			const AA_combo& combo = aa_edge_combos[combo_idx];
			const mass_t& combo_mass = combo.total_mass;
			const mass_t exp_connect_mass = node_mass + combo_mass;
		
			const mass_t min_connect_mass1 = exp_connect_mass - half_tolerance;
			const mass_t max_connect_mass1 = exp_connect_mass + half_tolerance;

			int n_idx;
			score_t max_score = NEG_INF;
			int     max_idx   = -1;
			for (n_idx = last_node_idx+1; n_idx<num_nodes; n_idx++)
			{
				const mass_t& next_node_mass = nodes[n_idx].mass;
				
				if (next_node_mass<min_connect_mass1)
					continue;

				if (next_node_mass>max_connect_mass1)
					break;

				if (nodes[n_idx].score>max_score)
				{
					max_idx = n_idx;
					max_score = nodes[n_idx].score;
				}
			}

			// if couldn't connect with small tolerance, try larger margin
			if (max_idx<0)
			{
				const mass_t min_connect_mass2 = exp_connect_mass - digest_tolerance;
				const mass_t max_connect_mass2 = exp_connect_mass + digest_tolerance;

				int n_idx;
				for (n_idx = last_node_idx+1; n_idx<num_nodes; n_idx++)
				{
					const mass_t& next_node_mass = nodes[n_idx].mass;
					
					if (next_node_mass<min_connect_mass2)
					{
						last_node_idx++;
						continue;
					}

					if (next_node_mass>max_connect_mass2)
						break;

					if (nodes[n_idx].score>max_score)
					{
						max_idx = n_idx;
						max_score = nodes[n_idx].score;
					}
				}
			}

		

			if (max_score>NEG_INF)
			{
				// need to check if distance between peaks is within tolerance
				mass_t min_diff = fabs(nodes[max_idx].mass-exp_connect_mass);

				// try seeing if there are peaks that can bridge them
				if (min_diff>tolerance)
				{
					const vector<BreakageFragment>& fragments1 = nodes[i].breakage.fragments;
					const vector<BreakageFragment>& fragments2 = nodes[max_idx].breakage.fragments;

					int f;
					for (f=0; f<strong_fragment_idxs.size(); f++)
					{
						const int& strong_frag_idx = strong_fragment_idxs[f];
						const int pos1 = nodes[i].breakage.get_position_of_frag_idx(strong_frag_idx);
						if (pos1<0)
							continue;

						const int pos2 = nodes[max_idx].breakage.get_position_of_frag_idx(strong_frag_idx);
						if (pos2<0)
							continue;

						mass_t mass_diff = fabs(fragments1[pos1].mass - fragments2[pos2].mass);
						
						const int frag_charge = config->get_fragment(strong_frag_idx).charge;
						if (frag_charge>1)
							mass_diff *= frag_charge;

						mass_t diff = fabs(mass_diff - combo_mass);

						if (diff < min_diff)
							min_diff = diff;
					}
				}

				// larger tolerance is allowed for the digest and terminal nodes because they can be
				// created relative to the terminal node
				if (nodes[max_idx].type == NODE_DIGEST || nodes[i].type == NODE_DIGEST)
				{
					if (min_diff>digest_tolerance)
							continue;
				}
				else
				if (nodes[max_idx].type == NODE_C_TERM || nodes[i].type == NODE_N_TERM )
				{
					if (min_diff>digest_tolerance)
						continue;
				}
				else
					if (min_diff>tolerance)
						continue;

				// add combo idx
				if (current_c_idx != max_idx)
				{
					current_e_idx = find_edge_idx_ij(i,max_idx);
					current_c_idx = max_idx;
				}

				if (current_e_idx>=0)
				{
					MultiEdge& edge = multi_edges[current_e_idx];
					
					// add the variant ptr and scores for this combo
					add_and_score_edge_variants(combo,edge);
				}
				else // create new edge
				{
					MultiEdge new_edge;
	
					new_edge.n_idx = i;
					new_edge.c_idx = max_idx;
					new_edge.n_break = &nodes[i].breakage;
					new_edge.c_break = &nodes[max_idx].breakage;
					new_edge.num_aa = 1;
				
					add_and_score_edge_variants(combo,new_edge);

					if (new_edge.variant_ptrs.size() == 0)
						continue;

					current_e_idx = multi_edges.size();
					nodes[i].out_edge_idxs.push_back(current_e_idx);
					nodes[max_idx].in_edge_idxs.push_back(current_e_idx);
					multi_edges.push_back(new_edge);
				}

					// set amino acid indicators
				int edge_aa = combo.amino_acids[0];
				out_aa_ind[edge_aa][i]=true;
				in_aa_ind[edge_aa][max_idx]=true;

			}
		}
	}
}




/********************************************************************
	Fills in double edges.
	Uses a precomputed list (from config) of aa comobos and masses to quickly
	determine if certain edges are present.
	Flag add_overlap_edges controls if to add a double edge when the 
	same edge can be constructed by single edges.
*********************************************************************/
void PrmGraph::fill_double_multi_edges(bool add_overlap_edges)
{
	const vector<AA_combo>& aa_edge_combos    = config->get_aa_edge_combos();
	const vector<int>& double_edge_combo_idxs = config->get_combo_idxs_by_length(2);
	const vector<int>& strong_fragment_idxs = config->get_all_strong_fragment_type_idxs();
	const mass_t tolerance = config->getTolerance();
	const mass_t pm_tolerance = (config->get_pm_tolerance() > tolerance) ? config->get_pm_tolerance() : tolerance * 0.75;
	const mass_t digest_tolerance = 1.5 * tolerance;
	const mass_t half_tolerance = 0.5 * tolerance;
	const int num_nodes  = nodes.size();
	const int num_combos = double_edge_combo_idxs.size();
	const mass_t max_combos_mass = aa_edge_combos[double_edge_combo_idxs[double_edge_combo_idxs.size()-1]].total_mass;;
	const mass_t max_node_mass = nodes[num_nodes-1].mass;

	const score_t overlap_thresh = 0;


	int i;
	for (i=0; i<num_nodes; i++)
	{
		Node& node = nodes[i];
		if (node.type == NODE_C_TERM)
			continue;

		const mass_t& node_mass = node.mass;
		
		int current_e_idx = -1;
		int current_c_idx = -1;

		int last_node_idx=i;
		
		int c;
		for (c=0; c<num_combos; c++)
		{
			const int& combo_idx = double_edge_combo_idxs[c];
			const AA_combo& combo = aa_edge_combos[combo_idx];
			const mass_t& combo_mass = combo.total_mass;
	
			// should this combo be check - if there is a single edge with one of the
			// amino acids in the combo skip this combo
			bool is_overlap_edge = false;
			int k;
			for (k=0; k<combo.num_aa; k++)
				if (out_aa_ind[combo.amino_acids[k]][i])
					break;
			if (k<combo.num_aa)
			{
				const int& overlap_aa = combo.amino_acids[k];
				if (i>0 && ! add_overlap_edges && 
					overlap_aa != Pro && overlap_aa != Gly && overlap_aa != His &&  overlap_aa != Ser)
					continue;
				is_overlap_edge = true;
			}

			const mass_t exp_connect_mass = node_mass + combo_mass;
			const mass_t min_connect_mass1 = exp_connect_mass - half_tolerance;
			const mass_t max_connect_mass1 = exp_connect_mass + half_tolerance;

			int n_idx;
			score_t max_score = NEG_INF;
			int     max_idx   = -1;
			for (n_idx = last_node_idx+1; n_idx<num_nodes; n_idx++)
			{
				const mass_t& next_node_mass = nodes[n_idx].mass;
				
				if (next_node_mass<min_connect_mass1)
					continue;

				if (next_node_mass>max_connect_mass1)
					break;

				if (nodes[n_idx].score>max_score)
				{
					max_idx = n_idx;
					max_score = nodes[n_idx].score;
				}
			}

			// if couldn't connect with small tolerance, try larger margin
			if (max_idx<0)
			{
				const mass_t min_connect_mass2 = exp_connect_mass - digest_tolerance;
				const mass_t max_connect_mass2 = exp_connect_mass + digest_tolerance;

				int n_idx;
				for (n_idx = last_node_idx+1; n_idx<num_nodes; n_idx++)
				{
					const mass_t& next_node_mass = nodes[n_idx].mass;
					
					if (next_node_mass<min_connect_mass2)
					{
						last_node_idx++;
						continue;
					}

					if (next_node_mass>max_connect_mass2)
						break;

					if (nodes[n_idx].score>max_score)
					{
						max_idx = n_idx;
						max_score = nodes[n_idx].score;
					}
				}
			}


			if (max_score>NEG_INF)
			{
				// need to check if distance between peaks is within tolerance
				mass_t min_diff = fabs(nodes[max_idx].mass-exp_connect_mass);

				// try seeing if there are peaks that can bridge them
				if (min_diff>tolerance)
				{
					const vector<BreakageFragment>& fragments1 = nodes[i].breakage.fragments;
					const vector<BreakageFragment>& fragments2 = nodes[max_idx].breakage.fragments;

					int f;
					for (f=0; f<strong_fragment_idxs.size(); f++)
					{
						const int& strong_frag_idx = strong_fragment_idxs[f];
						const int pos1 = nodes[i].breakage.get_position_of_frag_idx(strong_frag_idx);
						if (pos1<0)
							continue;

						const int pos2 = nodes[max_idx].breakage.get_position_of_frag_idx(strong_frag_idx);
						if (pos2<0)
							continue;

						mass_t mass_diff = fabs(fragments1[pos1].mass - fragments2[pos2].mass);
						
						const int frag_charge = config->get_fragment(strong_frag_idx).charge;
						if (frag_charge>1)
							mass_diff *= frag_charge;

						mass_t diff = fabs(mass_diff - combo_mass);

						if (diff < min_diff)
							min_diff = diff;
					}
				}

				// larger tolerance is allowed for the digest and terminal nodes because they can be
				// created relative to the terminal node
				if (nodes[max_idx].type == NODE_DIGEST || nodes[i].type == NODE_DIGEST)
				{
					if (min_diff>digest_tolerance)
							continue;
				}
				else
				if (nodes[max_idx].type == NODE_C_TERM || nodes[i].type == NODE_N_TERM )
				{
					if (min_diff>digest_tolerance)
						continue;
				}
				else
					if (min_diff>tolerance)
						continue;
				
				// add combo idx
				if (current_c_idx != max_idx)
				{
					current_e_idx = find_edge_idx_ij(i,max_idx);
					current_c_idx = max_idx;
				}

				// check that the node in the middle is not too good
				if (is_overlap_edge)
				{
				}

				if (current_e_idx>=0 && multi_edges[current_e_idx].num_aa == 2)
				{
					MultiEdge& edge = multi_edges[current_e_idx];
					
					add_and_score_edge_variants(combo,edge);

					if (is_overlap_edge)
						edge.ind_edge_overlaps = true;

				}
				else // create new edge
				{
					MultiEdge new_edge;
					
					new_edge.n_idx = i;
					new_edge.c_idx = max_idx;
					new_edge.n_break = &nodes[i].breakage;
					new_edge.c_break = &nodes[max_idx].breakage;
					new_edge.num_aa= 2;
					new_edge.ind_edge_overlaps=is_overlap_edge;

					add_and_score_edge_variants(combo,new_edge);
					if (new_edge.variant_ptrs.size() == 0)
						continue;

					current_e_idx = multi_edges.size();
					nodes[i].out_edge_idxs.push_back(current_e_idx);
					nodes[max_idx].in_edge_idxs.push_back(current_e_idx);
					multi_edges.push_back(new_edge);
				}
			}
		}
	}
}


/********************************************************************
	Fills in long edges (similar to the double function, only looks both
	ath the out edges of node i and the in edges of node max_idx to see
	if there is a possible overlap).
	Uses a precomputed list (from config) of aa comobos and masses to quickly
	determine if certain edges are present.
	Flag add_overlap_edges controls if to add a double edge when the 
	same edge can be constructed by single edges.
*********************************************************************/
void PrmGraph::fill_longer_multi_edges(int edge_length, bool add_overlap_edges)
{
	const vector<AA_combo>& aa_edge_combos    = config->get_aa_edge_combos();
	const vector<int>& edge_combo_idxs = config->get_combo_idxs_by_length(edge_length);
	const vector<int>& strong_fragment_idxs = config->get_all_strong_fragment_type_idxs();
	const mass_t tolerance = config->getTolerance();
	const mass_t pm_tolerance = (config->get_pm_tolerance() > tolerance) ? config->get_pm_tolerance() : tolerance * 0.75;
	const mass_t digest_tolerance = 1.5 * tolerance;
	const mass_t half_tolerance = 0.5 * tolerance;
	const int num_nodes  = nodes.size();
	const int num_combos = edge_combo_idxs.size();
	const mass_t max_combos_mass = aa_edge_combos[edge_combo_idxs[edge_combo_idxs.size()-1]].total_mass;;
	const mass_t max_node_mass = nodes[num_nodes-1].mass;
	

	int i;
	for (i=0; i<num_nodes; i++)
	{
		if (nodes[i].type == NODE_C_TERM)
			continue;

		const mass_t& node_mass = nodes[i].mass;
		const int length_minus_1 = edge_length-1;
		
		int current_e_idx = -1;
		int current_c_idx = -1;

		int last_node_idx=i;
		
		int c;
		for (c=0; c<num_combos; c++)
		{
			const int& combo_idx = edge_combo_idxs[c];
			const AA_combo& combo = aa_edge_combos[combo_idx];
			const mass_t& combo_mass = combo.total_mass;

			// should this combo be check - if there is a single edge with one of the
			// amino acids in the combo skip this combo
			bool is_overlap_edge = false;
			
			int a;
			for (a=0; a<edge_length; a++)
			{
				const int& aa = combo.amino_acids[a];
				if (out_aa_ind[aa][i])
					break;
			}

			if (a<edge_length)
			{
				if (! add_overlap_edges)
					continue;

				is_overlap_edge = true;
			}

			const mass_t exp_connect_mass = node_mass + combo_mass;
			const mass_t min_connect_mass1 = exp_connect_mass - half_tolerance;
			const mass_t max_connect_mass1 = exp_connect_mass + half_tolerance;

			int n_idx;
			score_t max_score = NEG_INF;
			int     max_idx   = -1;
			for (n_idx = last_node_idx+1; n_idx<num_nodes; n_idx++)
			{
				const mass_t& next_node_mass = nodes[n_idx].mass;
				
				if (next_node_mass<min_connect_mass1)
					continue;

				if (next_node_mass>max_connect_mass1)
					break;

				if (nodes[n_idx].score>max_score)
				{
					max_idx = n_idx;
					max_score = nodes[n_idx].score;
				}
			}

			// if couldn't connect with small tolerance, try larger margin
			if (max_idx<0)
			{
				const mass_t min_connect_mass2 = exp_connect_mass - digest_tolerance;
				const mass_t max_connect_mass2 = exp_connect_mass + digest_tolerance;

				int n_idx;
				for (n_idx = last_node_idx+1; n_idx<num_nodes; n_idx++)
				{
					const mass_t& next_node_mass = nodes[n_idx].mass;
					
					if (next_node_mass<min_connect_mass2)
					{
						last_node_idx++;
						continue;
					}

					if (next_node_mass>max_connect_mass2)
						break;

					if (nodes[n_idx].score>max_score)
					{
						max_idx = n_idx;
						max_score = nodes[n_idx].score;
					}
				}
			}


			if (max_score>NEG_INF)
			{
				// need to check if distance between peaks is within tolerance
				mass_t min_diff = fabs(nodes[max_idx].mass-exp_connect_mass);

				// try seeing if there are peaks that can bridge them
				if (min_diff>tolerance)
				{
					const vector<BreakageFragment>& fragments1 = nodes[i].breakage.fragments;
					const vector<BreakageFragment>& fragments2 = nodes[max_idx].breakage.fragments;

					int f;
					for (f=0; f<strong_fragment_idxs.size(); f++)
					{
						const int& strong_frag_idx = strong_fragment_idxs[f];
						const int pos1 = nodes[i].breakage.get_position_of_frag_idx(strong_frag_idx);
						if (pos1<0)
							continue;

						const int pos2 = nodes[max_idx].breakage.get_position_of_frag_idx(strong_frag_idx);
						if (pos2<0)
							continue;

						mass_t mass_diff = fabs(fragments1[pos1].mass - fragments2[pos2].mass);
						
						const int frag_charge = config->get_fragment(strong_frag_idx).charge;
						if (frag_charge>1)
							mass_diff *= frag_charge;

						mass_t diff = fabs(mass_diff - combo_mass);

						if (diff < min_diff)
							min_diff = diff;
					}
				}

				// larger tolerance is allowed for the digest and terminal nodes because they can be
				// created relative to the terminal node
				if (nodes[max_idx].type == NODE_DIGEST || nodes[i].type == NODE_DIGEST)
				{
					if (min_diff>digest_tolerance)
							continue;
				}
				else
				if (nodes[max_idx].type == NODE_C_TERM || nodes[i].type == NODE_N_TERM )
				{
					if (min_diff>digest_tolerance)
						continue;
				}
				else
					if (min_diff>tolerance)
						continue;

				// check if this edge overlaps a subpath of shorter edges with similar 
				// amino acids (already check the node i, so only check what comes in
				// node max_idx
				int a;
				for (a=0; a<edge_length; a++)
				{
					const int& aa = combo.amino_acids[a];
					if (in_aa_ind[aa][max_idx])
						break;
				}

				if (a<edge_length)
				{
					if (! add_overlap_edges)
						continue;

					is_overlap_edge = true;
				}
				
				// add combo idx
				if (current_c_idx != max_idx)
				{
					current_e_idx = find_edge_idx_ij(i,max_idx);
					current_c_idx = max_idx;
				}

				if (current_e_idx>=0 && multi_edges[current_e_idx].num_aa == edge_length)
				{
					MultiEdge& edge = multi_edges[current_e_idx];

					add_and_score_edge_variants(combo,edge);

					if (is_overlap_edge)
						edge.ind_edge_overlaps = true;
				}
				else // create new edge
				{
					MultiEdge new_edge;
		
					new_edge.n_idx = i;
					new_edge.c_idx = max_idx;
					new_edge.n_break = &nodes[i].breakage;
					new_edge.c_break = &nodes[max_idx].breakage;
					new_edge.num_aa= edge_length;
					new_edge.ind_edge_overlaps=is_overlap_edge;

					add_and_score_edge_variants(combo,new_edge);
					if (new_edge.variant_ptrs.size() == 0)
						continue;

					current_e_idx = multi_edges.size();
					nodes[i].out_edge_idxs.push_back(current_e_idx);
					nodes[max_idx].in_edge_idxs.push_back(current_e_idx);
					multi_edges.push_back(new_edge);
				}
			}
		}
	}
}


/******************************************************************
Removes nodes that have a score that is too low. Also removes edges
from adjacent nodes that are connected to this node
*******************************************************************/
void PrmGraph::prune_low_scoring_nodes()
{
	score_t min_score = -12.0;
	if (this->charge>2)
		min_score=-4.0;

	if (this->pm_with_19>2000)
		min_score+=2.0;
//	cout << "BEFORE " << nodes.size() << endl;
	int num_pruned=0;
	int i;
	for (i=1; i<nodes.size(); i++)
	{
		if (nodes[i].type != NODE_REG || nodes[i].score>min_score)
			continue;
		
		nodes[i].score = NEG_INF;
		nodes[i].active=0;
		num_pruned++;
	}

	remove_edges_from_inactive_nodes();
//	cout << "AFTER  " << nodes.size()-num_pruned << endl;
}


/***********************************************************
	assigns a value to each node's (log) rank field
	also sets max_node_score.
************************************************************/
void PrmGraph::rank_nodes_according_to_score()
{
	vector<score_pair> pairs;

	int i;
	pairs.resize(nodes.size());

	for (i=0; i<nodes.size(); i++)
	{
		pairs[i].idx=i;
		pairs[i].score = nodes[i].score;
	}

	sort(pairs.begin(),pairs.end());

	for (i=0; i<pairs.size(); i++)
		nodes[pairs[i].idx].log_rank = (float)log(2.0+i);

	max_node_score = NEG_INF;
	for (i=0; i<nodes.size(); i++)
		if (nodes[i].score>max_node_score)
			max_node_score = nodes[i].score;
}


/**********************************************************
	Sets the idx_max_in_score and idx_max_out_score fields
	for the nodes.
***********************************************************/
void PrmGraph::set_idxs_max_in_out_for_nodes()
{
	int i;
	for (i=0; i<nodes.size(); i++)
	{
		int e;
		score_t max_in_score = NEG_INF;
		nodes[i].idx_max_in_score_node=-1;
		
		const vector<int>& in_idxs = nodes[i].in_edge_idxs;
		for (e=0; e<in_idxs.size(); e++)
		{
			const int prev_node_idx = multi_edges[in_idxs[e]].n_idx;
			if (nodes[prev_node_idx].score>max_in_score)
			{
				nodes[i].idx_max_in_score_node = prev_node_idx;
				max_in_score = nodes[prev_node_idx].score; 
			}
		}

		score_t max_out_score = NEG_INF;
		nodes[i].idx_max_out_score_node=-1;

		const vector<int>& out_idxs = nodes[i].out_edge_idxs;
		for (e=0; e<out_idxs.size(); e++)
		{
			const int next_node_idx = multi_edges[out_idxs[e]].c_idx;
			if (nodes[next_node_idx].score>max_out_score)
			{
				nodes[i].idx_max_out_score_node = next_node_idx;
				max_out_score = nodes[next_node_idx].score; 
			}
		}		
	}
}


/*******************************************************************
********************************************************************/
void PrmGraph::fill_forbidden_idx_indicators_and_cumulative_scores()
{
	int i;

	const int num_nodes = nodes.size();
	
	cummulative_scores.clear();
	cummulative_scores.resize(num_nodes,0);
	cummulative_scores[0] = nodes[0].score;

	for (i=1; i<num_nodes; i++)
	{
		cummulative_scores[i]=cummulative_scores[i-1];
		if (nodes[i].score>0)
			cummulative_scores[i]+=nodes[i].score;
	}

	forbidden_node_idxs.clear();
	forbidden_node_idxs.resize(num_nodes,NEG_INF);

	const mass_t tolerance = 1.5 * config->getTolerance();
	const mass_t forbidden_mass = get_pm_with_19() - MASS_PROTON;
	
	int last_idx = num_nodes-1;
	for (i=0; i<num_nodes; i++)
	{
		const int pos = nodes[i].breakage.get_position_of_frag_idx(nodes[i].source_frag_type_idx);
		if (pos<0)
			continue;

		const int peak_idx = nodes[i].breakage.fragments[pos].peak_idx;
		
		int j;
		for (j=last_idx; j>i; j--)
		{
			const mass_t mass_sum = nodes[i].mass + nodes[j].mass;
			
			if (mass_sum - tolerance > forbidden_mass)
				continue;

			if (mass_sum + tolerance < forbidden_mass)
				break;

			int k;
			for (k=0; k<nodes[j].breakage.fragments.size(); k++)
			{
				if (nodes[j].breakage.fragments[k].peak_idx == peak_idx)
				{
					forbidden_node_idxs[i]=j;
					forbidden_node_idxs[j]=i;
					break;
				}
			}
		}
	}

//	for (i=0; i<nodes.size(); i++)
//		cout << i << "\t" << nodes[i].mass << "\t" << nodes[i].score << "\t" << forbidden_node_idxs[i] << "\t"
//			<< cummulative_scores[i] << endl;
}



// this function performs all the scoring operations on edges 
// (amino acid scores, missing cleavage scores etc.)
// *** Most of the scoring is now done through EdgeModel !!!!!
//
score_t PrmGraph::calc_edge_variant_score(const MultiEdge& edge, int num_aa, int *aa) const
{
	const Node& n_node = nodes[edge.n_idx];
	const Node& c_node = nodes[edge.c_idx];


	// give digest score only if there aren't other digest nodes with intensity
	if (num_aa == 1 && nodes[edge.n_idx].type == NODE_DIGEST && nodes[edge.c_idx].type == NODE_C_TERM)
	{
		//cout << config->get_aa2label()[aa[0]] << " ";
		// check for digest edge
		const vector<int>& c_term_digest_aas = config->get_c_term_digest_aas();
		if (c_term_digest_aas.size()>0)
		{
			int i;
			for (i=0; i<c_term_digest_aas.size(); i++)
				if (aa[0] == c_term_digest_aas[i])
					break;

			if (i == c_term_digest_aas.size())
			{
				return digest_node_score*-0.5;
			}

			bool other_digest_is_good=false;
			for (i=0; i<this->c_digest_node_idxs.size(); i++)
			{
				if (c_digest_node_idxs[i] == edge.n_idx)
					continue;
				if (nodes[c_digest_node_idxs[i]].breakage.fragments.size()>0)
					other_digest_is_good=true;
			}
			if (other_digest_is_good)
			{
				return digest_node_score*0.5;
			}
			return digest_node_score;
		}
	}

	if (num_aa == 1 && nodes[edge.c_idx].type == NODE_DIGEST && nodes[edge.n_idx].type == NODE_N_TERM)
	{
		const vector<int>& n_term_digest_aas = config->get_n_term_digest_aas();
		if (n_term_digest_aas.size()>0)
		{
			int i;
			for (i=0; i<n_term_digest_aas.size(); i++)
				if (aa[0] == n_term_digest_aas[i])
					break;

			if (i == n_term_digest_aas.size())
				return digest_node_score*-0.5;
			
			bool other_digest_is_good=false;
			for (i=0; i<this->n_digest_node_idxs.size(); i++)
			{
				if (n_digest_node_idxs[i] == edge.c_idx)
					continue;
				if (nodes[n_digest_node_idxs[i]].score>0)
					other_digest_is_good=true;
			}
			if (other_digest_is_good)
				return 0;
			return digest_node_score;
		}
	}

	if (num_aa >1)
	{
	//	return ((num_aa-1)*model->get_missing_breakage_score(charge,size_idx,edge.c_break->region_idx));
	//	return ((num_aa-1)*-6);
	}

	return 0;

}













// removes all edges to and from nodes with the active flag set to 0
void PrmGraph::remove_edges_from_inactive_nodes()
{
	int i;

	for (i=0; i<nodes.size(); i++)
	{
		if (nodes[i].active)
			continue;

		int j;
		Node& node = nodes[i];

		for (j=0; j<node.in_edge_idxs.size(); j++)
		{
			// remove edge from other node's list
			const int e_idx = node.in_edge_idxs[j];
			Node& other_node = nodes[multi_edges[e_idx].n_idx];
			int k;

			for (k=0; k<other_node.out_edge_idxs.size(); k++)
				if (other_node.out_edge_idxs[k] == e_idx)
					break;

			if (k== other_node.out_edge_idxs.size())
			{
				cout << "Error: missing out edge idx!" << endl;
				exit(1);
			}

			other_node.out_edge_idxs[k] = other_node.out_edge_idxs[other_node.out_edge_idxs.size()-1];
			other_node.out_edge_idxs.pop_back();			
		}

		for (j=0; j<node.out_edge_idxs.size(); j++)
		{
			// remove edge from other node's list
			const int e_idx = node.out_edge_idxs[j];
			Node& other_node = nodes[multi_edges[e_idx].c_idx];
			int k;

			for (k=0; k<other_node.in_edge_idxs.size(); k++)
				if (other_node.in_edge_idxs[k] == e_idx)
					break;

			if (k== other_node.in_edge_idxs.size())
			{
				cout << "Error: missing in edge idx!" << endl;
				exit(1);
			}

			other_node.in_edge_idxs[k] = other_node.in_edge_idxs[other_node.in_edge_idxs.size()-1];
			other_node.in_edge_idxs.pop_back();			
		}

		node.in_edge_idxs.clear();
		node.out_edge_idxs.clear();
	}
}






struct idx_score {
	idx_score() : edge_idx(-1), score(NEG_INF) {}
	bool operator< (const idx_score& other) const
	{
		return score>other.score;
	}
	int edge_idx;
	score_t score;
};



// sorts edges according to the value to which they can possibly lead
// uses the max_gains table which state for each node i and length n, what is the maximal
// score attainable by using i + n amino acids in the graph
void PrmGraph::sort_outgoing_edges_according_to_max_gains(const vector< vector< score_t > >& max_gains)
{
	int i;
	const int last_size_idx = max_gains[0].size()-1;
	for (i=0; i<nodes.size(); i++)
	{
		if (nodes[i].out_edge_idxs.size()==0)
			continue;
		
		int j;
		vector<idx_score> pairs;
		pairs.resize(nodes[i].out_edge_idxs.size());
		for (j=0; j<nodes[i].out_edge_idxs.size(); j++)
		{
			const int edge_idx = nodes[i].out_edge_idxs[j];
			const MultiEdge& edge = multi_edges[edge_idx];
			pairs[j].edge_idx= edge_idx;
			pairs[j].score = max_gains[edge.c_idx][last_size_idx] +
							nodes[edge.c_idx].score + edge.max_variant_score;
		}
		sort(pairs.begin(),pairs.end());

		for (j=0; j<nodes[i].out_edge_idxs.size(); j++)
			nodes[i].out_edge_idxs[j]=pairs[j].edge_idx;
	}
}

/***********************************************************************
// returns the optimal ordering of nodes for the search
************************************************************************/
void PrmGraph::get_node_ordering_according_to_max_gains(
	 vector< vector< score_t > >& max_gains_for_length,
	 vector<int>& node_order) const
{
	vector<idx_score> pairs;
	pairs.resize(nodes.size());
	int i;
	const int last_size_idx = max_gains_for_length[0].size()-1;
	for (i=0; i<nodes.size(); i++)
	{
		pairs[i].edge_idx=i;
		pairs[i].score = nodes[i].score + max_gains_for_length[i][last_size_idx];
	}

	sort(pairs.begin(),pairs.end());

	node_order.resize(nodes.size());
	
	for (i=0; i<pairs.size(); i++)
		node_order[i]=pairs[i].edge_idx;
}


// sorts edges according to the value to which they lead
void PrmGraph::sort_outgoing_edges()
{
	int i;
	for (i=0; i<nodes.size(); i++)
	{
		if (nodes[i].out_edge_idxs.size()==0)
			continue;
		
		int j;
		vector<idx_score> pairs;
		pairs.resize(nodes[i].out_edge_idxs.size());
		for (j=0; j<nodes[i].out_edge_idxs.size(); j++)
		{
			pairs[j].edge_idx=nodes[i].out_edge_idxs[j];
			pairs[j].score = multi_edges[nodes[i].out_edge_idxs[j]].max_variant_score + 
							 nodes[multi_edges[nodes[i].out_edge_idxs[j]].c_idx].score;
		}
		sort(pairs.begin(),pairs.end());

		for (j=0; j<nodes[i].out_edge_idxs.size(); j++)
			nodes[i].out_edge_idxs[j]=pairs[j].edge_idx;
	}
}





struct edge_idx_pair {
	bool operator< (const edge_idx_pair& other) const
	{
		return n_idx<other.n_idx;
	}

	int edge_idx;
	int n_idx;
};


/*******************************************************************
// creates a path object from a collection of edges that are assumed
// to correspond to a path in the graph
********************************************************************/
void PrmGraph::create_path_from_edges(vector<int>& edge_idxs, MultiPath& path) const
{
	int i;
	vector<edge_idx_pair> pairs;

	if (edge_idxs.size()==0)
	{
		path.path_score = 0;
		return;
	}
	
	for (i=0; i<edge_idxs.size(); i++)
	{
		edge_idx_pair p;
		p.edge_idx = edge_idxs[i];
		p.n_idx = multi_edges[edge_idxs[i]].n_idx;

		pairs.push_back(p);
	}
	sort(pairs.begin(),pairs.end());

	path.edge_idxs.resize(pairs.size());
	for (i=0; i<pairs.size(); i++)
		path.edge_idxs[i]=pairs[i].edge_idx;

	path.n_term_mass= nodes[multi_edges[edge_idxs[0]].n_idx].mass;
	path.c_term_mass= nodes[multi_edges[edge_idxs[edge_idxs.size()-1]].c_idx].mass;

	for (i=1; i<path.edge_idxs.size(); i++)
	{
		if (multi_edges[path.edge_idxs[i]].n_idx != multi_edges[path.edge_idxs[i-1]].c_idx)
		{
			cout << "Error: inconsistent edges when creating path!" << endl;
			exit(1);
		}
	}

	// collect breakage info and edges

	path.edge_idxs = edge_idxs;
	path.breakages.clear();
	path.node_idxs.clear();
	for (i=0; i<edge_idxs.size(); i++)
	{
		path.breakages.push_back(multi_edges[edge_idxs[i]].n_break);
		path.node_idxs.push_back(multi_edges[edge_idxs[i]].n_idx);
	}
	path.breakages.push_back(multi_edges[edge_idxs[i-1]].c_break);
	path.node_idxs.push_back(multi_edges[edge_idxs[i-1]].c_idx);
}



// finds the highest scoring continuous subpath for a given peptide in the graph
SeqPath PrmGraph::get_highest_scoring_subpath(const Peptide& peptide, mass_t start_mass) const
{
	SeqPath ret_path;
	const vector<int>& path_aas = peptide.get_amino_acids();
	mass_t pre_mass = pre_mass = start_mass;
	mass_t double_tolerance = config->getTolerance() * 3.0;
	vector<bool> use_as_start_pos;
	int i;
	
	use_as_start_pos.resize(nodes.size(),true);
	ret_path.path_score = 0;

	// give start idx double tolerance
	for (i=0; i<path_aas.size(); i++)
	{
		int j;
		PeakRange nr = this->get_nodes_in_range(pre_mass - double_tolerance, pre_mass + double_tolerance);
		for (j=0; j<nr.num_peaks; j++)
		{
			int node_idx = nr.low_idx+j;
			if (! use_as_start_pos[node_idx])
				continue;

			// find max correct subpath from this node
			SeqPath path;

			path.n_term_mass = nodes[node_idx].mass;
			path.c_term_mass = nodes[node_idx].mass;
			path.path_score = 0;
			path.positions.clear();

			// loop until end is reached
			int k=i;
			
			int lass_good_edge_idx = -1;
			while (k<path_aas.size())
			{
				const Node& node = nodes[node_idx];
				int e;
				for (e=0; e<node.out_edge_idxs.size(); e++)
				{
					const int& e_idx = node.out_edge_idxs[e];
					const MultiEdge& edge = multi_edges[e_idx];

					int var_idx = edge.get_variant_idx(1,(int *)&path_aas[k]);
					if (var_idx<0 && i<path_aas.size()-1)
						var_idx = edge.get_variant_idx(2,(int *)&path_aas[k]);

					if (var_idx>=0)
					{
						path.add_edge_variant(edge,e_idx,var_idx);
						k+=edge.num_aa;
						lass_good_edge_idx = e_idx;
						break;
					}
				}
				if (e == node.out_edge_idxs.size())
					break;
			}

			// add last position
			if (lass_good_edge_idx>=0)
			{
				const MultiEdge& last_edge = multi_edges[lass_good_edge_idx];
				PathPos last_pos;

				last_pos.breakage = last_edge.c_break;
				last_pos.edge_idx =-1;
				last_pos.mass = last_edge.c_break->mass;
				last_pos.node_idx = last_edge.c_idx;
				last_pos.node_score = last_edge.c_break->score;
				path.positions.push_back(last_pos);
				
				path.path_score += last_pos.node_score;
				path.c_term_mass = last_pos.mass;
			}

			// check if this path is better
			if (path.path_score>ret_path.path_score)
				ret_path = path;
		}	

		pre_mass += config->get_aa2mass()[path_aas[i]];
	}
	
	return ret_path;
}





// finds the longest continuous subpath for a given peptide in the graph
SeqPath PrmGraph::get_longest_subpath(const Peptide& peptide, mass_t start_mass, bool verbose)
{
	const int max_edge_length     = config->get_max_edge_length();
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	vector<mass_t> prefix_masses;
	vector<int> path_aas = peptide.get_amino_acids();
	vector<int> path_edges;
	vector<int> path_aa_count;
	int a_start=-1;
	mass_t pre_mass = 0;
	mass_t double_tolerance = config->getTolerance() *  5.0;
	int i;

	path_edges.clear();

	for (i=0; i<path_aas.size(); i++)
		if (path_aas[i]==Ile)
			path_aas[i]=Leu;

	prefix_masses.push_back(start_mass);
	for (i=0; i<path_aas.size(); i++)
		prefix_masses.push_back(prefix_masses[i]+aa2mass[path_aas[i]]);
	
	// give start idx double tolerance
	for (i=0; i<path_aas.size(); i++)
	{
		int j;
		const mass_t pre_mass = prefix_masses[i];
		PeakRange nr = this->get_nodes_in_range(pre_mass - double_tolerance, pre_mass + double_tolerance);
		for (j=0; j<nr.num_peaks; j++)
		{
			int node_idx = nr.low_idx+j;
			vector<int> e_idxs, aa_count;
			e_idxs.clear();
			aa_count.clear();
			
		
			// find max correct subpath from this node

			int k;
			int n=node_idx;
			for (k=i; k<path_aas.size(); k++)
			{
				int q;
				int curr_aa = path_aas[k];
				bool found_edge=false;
				
				// look for aas in edges
				// start with short edges and increase length
				int num_aa;
				for (num_aa = 1; num_aa<=max_edge_length; num_aa++)
				{
					if (num_aa + k > path_aas.size())
						break;

					for (q=0; q<nodes[n].out_edge_idxs.size(); q++)
					{
						int e_idx = nodes[n].out_edge_idxs[q];
						const MultiEdge& edge = multi_edges[e_idx];

						if (edge.has_variant(num_aa,&path_aas[k]) )
						{
							found_edge=true;
							e_idxs.push_back(e_idx);
							aa_count.push_back(num_aa);
							n = edge.c_idx;

							k+= num_aa -1;
							
							break;
						}			
					}

					if (found_edge)
						break;
				}
					
				if (! found_edge)
					break;
			}

			if (e_idxs.size()>path_edges.size())
			{
				path_edges = e_idxs;
				path_aa_count = aa_count;
				a_start=i;
			}
		}
	}



	// create path
	// need to add treatment of modifications to terminals


	SeqPath ret_path;
	vector<int> var_idxs;
	var_idxs.clear();
	int a_pos= a_start;

	if (path_edges.size()==0)
	{
		ret_path.pm_with_19 = this->pm_with_19;
		ret_path.charge = this->charge;
		ret_path.prm_ptr = (PrmGraph *)this;
		ret_path.path_score = NEG_INF;
		return ret_path;
	}

	if (verbose)
		cout << endl << "Starting at aa " << a_pos << endl;

	for (i=0; i<path_edges.size(); i++)
	{
		const MultiEdge& edge = multi_edges[path_edges[i]];

		int variant_idx = edge.get_variant_idx(edge.num_aa,&path_aas[a_pos]);
		if (variant_idx<0)
		{
			cout << "Error: edge does not contain variant idx!" << endl;
			exit(1);
		}
		
		var_idxs.push_back(variant_idx);

		if (verbose)
			cout << i<< " " << edge.n_idx << " (" << path_aa_count[i] << ") " << 
			edge.c_idx << endl;

		PathPos pos;
		
		pos.node_idx = edge.n_idx;
		pos.node_score = nodes[edge.n_idx].score;
		pos.edge_variant_score = edge.variant_scores[variant_idx];
		pos.breakage   = (Breakage *)&nodes[edge.n_idx].breakage;
		pos.edge_idx   = path_edges[i];
		pos.mass       = nodes[edge.n_idx].mass;
		pos.aa         = path_aas[a_pos];

		ret_path.positions.push_back(pos);

		a_pos++;

		// add skipped positions for the rest of the amino acids in the edge
		int j;
		for (j=1; j<path_aa_count[i]; j++)
		{
			PathPos pos;
		
			pos.aa         = path_aas[a_pos++];
			pos.edge_variant_score =0;
			pos.node_score = 0;

			ret_path.positions.push_back(pos);

		}		
	}

	// add last pos
	if (path_edges.size()>0)
	{
		const MultiEdge& last_edge = multi_edges[path_edges[path_edges.size()-1]];
		PathPos pos;

		pos.node_idx = last_edge.c_idx;
		pos.node_score = nodes[last_edge.c_idx].score;
		pos.breakage   = (Breakage *)&nodes[last_edge.c_idx].breakage;
		pos.mass       = nodes[last_edge.c_idx].mass;

		ret_path.positions.push_back(pos);

		ret_path.n_term_mass = ret_path.positions[0].mass;
		ret_path.c_term_mass = pos.mass;

		ret_path.path_score =0;
		int i;

		for (i=0; i<ret_path.positions.size(); i++)
		{
			if (ret_path.positions[i].node_idx>=0)
				ret_path.path_score += ret_path.positions[i].node_score + 
				ret_path.positions[i].edge_variant_score;
		//	cout << ret_path.positions[i].node_score << "  " <<  ret_path.positions[i].edge_varaint_score <<
		//			" t: " << ret_path.path_score << endl;
		}
	}
		

	// rescore according to the combo scores
	if (this->has_node_combo_scores)
	{
		vector<PathPos>& positions = ret_path.positions;
		ret_path.path_score =0;
		ret_path.multi_path_score =0;

		// compute the node scores based on the edge variants
		// N-term node scor

		positions[0].node_score = model->get_node_combo_score(this,positions[0].node_idx,-1,0,
													positions[0].edge_idx,var_idxs[0]);
		score_t first_edge_var_score = positions[0].edge_variant_score;

		ret_path.path_score += positions[0].node_score + first_edge_var_score;
		ret_path.multi_path_score += nodes[positions[0].node_idx].score + first_edge_var_score;
		// middle nodes score
		int prev_edge_idx = positions[0].edge_idx;
		int k;
		int p=0;
		for (k=1; k<var_idxs.size(); k++)
		{
			const int num_edge_aa = multi_edges[prev_edge_idx].num_aa;
			p+= num_edge_aa;
			const PathPos& pos = positions[p];
			positions[p].node_score = model->get_node_combo_score(this,pos.node_idx,prev_edge_idx,var_idxs[k-1],
													pos.edge_idx,var_idxs[k]);
			ret_path.path_score += pos.node_score;
			ret_path.path_score += pos.edge_variant_score;

			ret_path.multi_path_score += nodes[pos.node_idx].score;
			ret_path.multi_path_score += pos.edge_variant_score;
			prev_edge_idx = pos.edge_idx;
		}

		// C-term node score
		PathPos& last_pos = positions[positions.size()-1];
		last_pos.node_score = model->get_node_combo_score(this,last_pos.node_idx,prev_edge_idx,var_idxs[k-1],-1,0);

		ret_path.path_score += last_pos.node_score;
		ret_path.multi_path_score += nodes[last_pos.node_idx].score;

		// check for forbidden idx pairs
		if (forbidden_node_idxs.size()>0)
		{
			int num_forbidden=0;
			for (k=0; k<positions.size(); k++)
			{
				const int node_idx = positions[k].node_idx;
				if (node_idx<0)
					continue;

				const int mirror_node = forbidden_node_idxs[node_idx];
				if (mirror_node<0)
					continue;

				int j;
				for (j=k+1; j<positions.size(); j++)
					if (positions[j].node_idx == mirror_node)
						break;
				
				if (j<positions.size())
				{
//					cout << nodes[node_idx].mass << " <-> " << nodes[mirror_node].mass << endl;
					num_forbidden++;
				}
			}
	//		if (num_forbidden>0)
	//			source_spectrum->print_expected_by();
	//		cout << endl;
			ret_path.path_score -= num_forbidden * config->get_forbidden_pair_penalty();
			ret_path.multi_path_score -= num_forbidden * config->get_forbidden_pair_penalty();
		}
	}

	ret_path.adjust_complete_sequence_penalty(this);

	ret_path.make_seq_str(config);
	
	ret_path.pm_with_19 = this->pm_with_19;
	ret_path.charge = this->charge;
	ret_path.prm_ptr = (PrmGraph *)this;

	return ret_path; 
}


score_t SeqPath::adjust_complete_sequence_penalty(PrmGraph *prm)
{
	Config *config = prm->get_config();
	const mass_t pm_with_19 = prm->get_pm_with_19();
	if (this->n_term_mass>0 || this->c_term_mass+25<pm_with_19)
		return 0;

	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const mass_t tolerance = config->getTolerance();
	vector<int> aas;
	get_amino_acids(aas);
	mass_t mass=0;
	int i;
	for (i=0; i<aas.size(); i++)
		mass+=aa2mass[aas[i]];

	mass_t diff = fabs(mass+MASS_OHHH-pm_with_19);
//	cout << "Diff: " << diff << endl;
	if (diff<tolerance)
		return 0;

	score_t penalty = (diff/tolerance) * config->get_digest_score();
	path_score -= penalty;
	return penalty;
}

// returns a path from a graph that contains only good nodes
// edges might not be contiguous
SeqPath PrmGraph::get_path_from_peptide_prm_graph(const Peptide& peptide) const
{
	const int max_edge_length     = config->get_max_edge_length();
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	vector<mass_t> prefix_masses;
	vector<int> path_aas = peptide.get_amino_acids();
	vector<int> path_edges;
	vector<int> edge_pos_starts;
	vector<int> path_aa_count;
	int a_start=-1;
	mass_t pre_mass = 0;
	mass_t double_tolerance = config->getTolerance() * 1.75;
	int i;

	path_edges.clear();

	for (i=0; i<path_aas.size(); i++)
		if (path_aas[i]==Ile)
			path_aas[i]=Leu;

	prefix_masses.push_back(0);
	for (i=0; i<path_aas.size(); i++)
		prefix_masses.push_back(prefix_masses[i]+aa2mass[path_aas[i]]);
	
	path_edges.clear();
	for (i=0; i<this->multi_edges.size(); i++)
	{
		const MultiEdge& edge = multi_edges[i];
		if (edge.n_idx +1 == edge.c_idx)
		{
			int j;
			for (j=0; j<prefix_masses.size(); j++)
				if (fabs(nodes[edge.n_idx].mass - prefix_masses[j])<0.1)
					break;
			
			if (j== prefix_masses.size())
			{
				if (nodes[edge.n_idx].type == NODE_DIGEST)
					continue;

				cout << "Error: bad path from peptide prm!" << endl;
				cout << "Peptide: " << peptide.as_string(config) << endl;
				cout << "Edge: " << i << "   n: " << edge.n_idx << "  c: " << edge.c_idx << endl;
				cout << "Prefix masses: " ;
				int k;
				for (k=0; k<prefix_masses.size(); k++)
					cout << setprecision(6) << prefix_masses[k] << " ";
				cout << endl;
				print();
				exit(1);
			}
			int variant_idx = edge.get_variant_idx(edge.num_aa,&path_aas[j]);
			if (variant_idx>=0)
			{
				path_edges.push_back(i);
				edge_pos_starts.push_back(j);
			}
		}
	}

	// create path
	// need to add treatment of modifications to terminals


	SeqPath ret_path;

	ret_path.positions.resize(path_aas.size()+1);
	for (i=0; i<path_aas.size(); i++)
	{
		ret_path.positions[i].aa = path_aas[i];
		ret_path.positions[i].edge_variant_score = 0;
		ret_path.positions[i].node_score = 0;
	}

	for (i=0; i<path_edges.size(); i++)
	{
		const MultiEdge& edge = multi_edges[path_edges[i]];
		const int a_pos = edge_pos_starts[i];

		int variant_idx = edge.get_variant_idx(edge.num_aa,&path_aas[a_pos]);
		if (variant_idx<0)
		{
			cout << "Error: edge does not contain variant idx!" << endl;
			cout << "Peptide: " << peptide.as_string(config) << endl;
			cout << "a_pos: " << a_pos << endl;
			cout << "Edge: " << path_edges[i] << endl;
			print();
			exit(1);
		}
		

		PathPos& pos = ret_path.positions[a_pos];
		
		pos.node_idx   = edge.n_idx;
		pos.node_score = nodes[edge.n_idx].score;
		pos.edge_variant_score = edge.variant_scores[variant_idx];
		pos.breakage   = (Breakage *)&nodes[edge.n_idx].breakage;
		pos.edge_idx   = path_edges[i];
		pos.mass       = nodes[edge.n_idx].mass;

		PathPos& next_pos = ret_path.positions[a_pos + edge.num_aa];
		next_pos.node_idx = edge.c_idx;
		next_pos.node_score = nodes[edge.c_idx].score;
		next_pos.breakage   = (Breakage *)&nodes[edge.c_idx].breakage;
		next_pos.mass       = nodes[edge.c_idx].mass;

	}

	ret_path.n_term_mass = 0;
	ret_path.c_term_mass = prefix_masses[prefix_masses.size()-1];

	// first pos
	PathPos& first_pos = ret_path.positions[0];
	first_pos.node_idx =0;
	first_pos.node_score = nodes[0].score;
	first_pos.breakage = (Breakage *)&nodes[0].breakage;
	first_pos.mass	 = nodes[0].mass;

	// last pos
	PathPos& last_pos = ret_path.positions[ret_path.positions.size()-1];
	last_pos.node_idx = nodes.size()-1;
	last_pos.node_score = nodes[last_pos.node_idx].score;
	last_pos.breakage   = (Breakage *)&nodes[last_pos.node_idx].breakage;
	last_pos.mass		= nodes[last_pos.node_idx].mass;

	int num_scores=0;
	ret_path.path_score =0;
	for (i=0; i<ret_path.positions.size(); i++)
		if (ret_path.positions[i].breakage)
		{
			num_scores++;
			ret_path.path_score += ret_path.positions[i].node_score + ret_path.positions[i].edge_variant_score;
		}

	if (num_scores<2)
		ret_path.path_score=-1000;

	ret_path.make_seq_str(config);
	ret_path.adjust_complete_sequence_penalty((PrmGraph *)this);
	
	ret_path.pm_with_19 = this->pm_with_19;
	ret_path.prm_ptr = (PrmGraph *)this;

	return ret_path; 
}







/***********************************************************************
// returns the idxs of nodes correponding to the expected breakages of the peptide
// returns -1 for an idx of a missing node
************************************************************************/
void PrmGraph::get_all_correct_node_idxs(const Peptide& peptide, vector<int>& idxs) const
{
	const mass_t tolerance = config->getTolerance()*1.5;
	vector<mass_t> break_masses;
	int i;

	peptide.calc_expected_breakage_masses(config,break_masses);

	idxs.clear();
	for (i=0; i<break_masses.size(); i++)
	{
		int idx= this->get_max_score_node(break_masses[i],tolerance);
	//	if (idx<0)
	//		continue;
		idxs.push_back(idx);
	}
}


/***********************************************************************
   returns the idxs of nodes correponding to the mirrors of the
   expected breakages of the peptide
   returns -1 for an idx of a missing node
************************************************************************/
void PrmGraph::get_all_mirror_node_idxs(const Peptide& peptide, vector<int>& idxs) const
{
	const mass_t tolerance = config->getTolerance()*1.5;
	vector<mass_t> break_masses;
	int i;

	peptide.calc_expected_breakage_masses(config,break_masses);


	idxs.clear();
	for (i=break_masses.size()-1; i>=0; i--)
	{
		int idx= this->get_max_score_node(pm_with_19 - break_masses[i] + MASS_PROTON,tolerance);
	//	if (idx<0)
	//		continue;

		idxs.push_back(idx);
	}
}

/***********************************************************************
// returns the idxs of nodes correponding to the expected breakages of the peptide
************************************************************************/
void PrmGraph::get_relevant_node_idxs(const Peptide& peptide, vector<int>& idxs) const
{
	const mass_t tolerance = config->getTolerance()*2.0;
	vector<mass_t> break_masses;
	int i;

	peptide.calc_expected_breakage_masses(config,break_masses);

	idxs.clear();
	for (i=0; i<break_masses.size(); i++)
	{
		int idx= this->get_max_score_node(break_masses[i],tolerance);
		if (idx<0)
			continue;

		idxs.push_back(idx);
	}
}





// computes amino acid probabilities according to the AminoAcidProbs models
// if no model exists for the seq_path, aa probabilties of -1 wil be assigned
// and false will be returned
bool PrmGraph::calc_amino_acid_probs(SeqPath& seq_path, int seq_rank)
{
	const AminoAcidProbs *aap_models = model->get_amino_acid_probs_ptr();
	vector<int> path_aas;
	bool ret_val = true;

	seq_path.get_amino_acids(path_aas);
	int i;
	for (i=0; i<seq_path.positions.size()-1; i++)
	{
		PathPos& pos = seq_path.positions[i];

		if (pos.edge_idx<0)
			continue;

		// check if the score was already computed for the specific edge variant
		MultiEdge edge = multi_edges[pos.edge_idx];
		const int edge_var_idx = edge.get_variant_idx(edge.num_aa,&path_aas[i]);
		if (edge.variant_probs.size()>edge_var_idx &&
			edge.variant_probs[edge_var_idx]>=0)
		{
			pos.edge_variant_prob = edge.variant_probs[edge_var_idx];
			continue;
		}

		if (edge.variant_probs.size()<edge.variant_ptrs.size())
			edge.variant_probs.resize(edge.variant_ptrs.size(),-1.0);

		float prob = aap_models->calc_variant_prob(this,pos.edge_idx,pos.variant_ptr,seq_rank);

		edge.variant_probs[edge_var_idx]=prob;
		pos.edge_variant_prob = prob;

		if (prob<0)
			ret_val = false;
	}
	return ret_val;
}



void Node::print(Config *config, ostream& os) const
{
	os << mass << " " << breakage.region_idx 
	   << ","  << source_frag_type_idx << "," << type << " s: " << setw(4) << score << " ";

	breakage.print_fragments(config,os);
	os << endl;
}




void PrmGraph::print_multi_edges(int node_idx, bool print_edge_scores) const
{
	const vector<AA_combo>& aa_edge_combos = config->get_aa_edge_combos();
	int e;
	const vector<string>& aa2label = config->get_aa2label();
	for (e=0; e<nodes[node_idx].out_edge_idxs.size(); e++)
	{
		int a;
		int edge_idx = nodes[node_idx].out_edge_idxs[e];
		
		cout << " (n:" << multi_edges[edge_idx].c_idx << " e:" << edge_idx << " ";

		int j;
		for (j=0; j<multi_edges[edge_idx].variant_ptrs.size(); j++)
		{
			
			int *p = multi_edges[edge_idx].variant_ptrs[j];
			int num_aa = *p++;
			int *aas = p;
			for (a=0; a<num_aa; a++)
				cout << aa2label[aas[a]];
			cout <<" ";
			if (print_edge_scores)
				cout << fixed << setprecision(2) << multi_edges[edge_idx].variant_scores[j] <<" ";
		}
		cout << ") ";
	}
	cout << endl;
}


void PrmGraph::print_only_scores() const
{
	int i;

	cout << setprecision(3) << fixed;
	for (i=0; i<nodes.size(); i++)
	{
		if (nodes[i].score>NEG_INF && nodes[i].score !=0 )
			cout << nodes[i].mass << " " << nodes[i].score << endl;
	}
}


void PrmGraph::print(ostream& os, bool print_edge_scores) const
{
	int i;
	for (i=0; i<nodes.size(); i++)
	{
		os << fixed << setw(4) << left << i << " " << setw(8) << setprecision(3) << nodes[i].mass << " " << nodes[i].breakage.region_idx 
			   << ","  << nodes[i].source_frag_type_idx << "," << nodes[i].type << " s: " << setw(4) << nodes[i].score << " ";

		nodes[i].breakage.print_fragments(config,os);
		os << endl;
		if (nodes[i].out_edge_idxs.size()>0)
			print_multi_edges(i,print_edge_scores);

		os << endl;

	}
}


void PrmGraph::print_with_combo_tables() const
{
	int i;
	for (i=0; i<nodes.size(); i++)
	{
		cout << fixed << setw(4) << left << i << " " << setw(8) << setprecision(3) << nodes[i].mass << " " << nodes[i].breakage.region_idx 
		   << ","  << nodes[i].source_frag_type_idx << "," << nodes[i].type << " s: " << setw(4) << nodes[i].score << " ";

		nodes[i].breakage.print_fragments(config);
		cout << endl;
		if (nodes[i].out_edge_idxs.size()>0)
			print_multi_edges(i);

		cout << endl;

		if (this->has_node_combo_scores)
		{
			nodes[i].print_combo_table(this);
			cout << endl;
		}
	}	
}




void PrmGraph::print_with_multi_edges() const
{
	int i;
	for (i=0; i<nodes.size(); i++)
	{
		cout << setw(4) << left << i << " " << setw(8) << setprecision(3) << nodes[i].mass << " " << nodes[i].breakage.region_idx 
		   << ","  << nodes[i].source_frag_type_idx << "," << nodes[i].type << " s: " << setw(4) << nodes[i].score << " ";

		nodes[i].breakage.print_fragments(config,cout);
		cout << endl;
		if (nodes[i].out_edge_idxs.size()>0)
			print_multi_edges(i,true);

		cout << endl;
	}

}

void Node::print_combo_table(const PrmGraph *prm, ostream& os) const
{
	Config *config = prm->get_config();
	const vector<string>& aa2label = config->get_aa2label();
	int i;

	if (mass<1.0)
		return;

	map<ScoreComboLoc,score_t>::const_iterator it;

	vector< vector< score_t > > max_scores;
	vector< vector< int > > max_in_vars,max_out_vars;

	max_scores.resize(in_edge_idxs.size()+1);
	max_in_vars.resize(in_edge_idxs.size()+1);
	max_out_vars.resize(in_edge_idxs.size()+1);
	for (i=0; i<max_scores.size(); i++)
	{
		max_scores[i].resize(out_edge_idxs.size()+1,NEG_INF);
		max_in_vars[i].resize(out_edge_idxs.size()+1,NEG_INF);
		max_out_vars[i].resize(out_edge_idxs.size()+1,NEG_INF);
	}

	for (it=score_combos.begin(); it != score_combos.end(); it++)
	{
		const ScoreComboLoc& loc = (*it).first;
		int p=0;
		if (loc.in_edge>=0)
		{
			int j;
			for (j=0; j<in_edge_idxs.size(); j++)
				if (in_edge_idxs[j]==loc.in_edge)
					break;
			p=j+1;
		}

		int q=0;
		if (loc.out_edge>=0)
		{
			int j;
			for (j=0; j<out_edge_idxs.size(); j++)
				if (out_edge_idxs[j]==loc.out_edge)
					break;
			q=j+1;
		}

		score_t comb_score = (*it).second;
		if (comb_score>max_scores[p][q])
		{
			max_scores[p][q]=comb_score;
			max_in_vars[p][q]=(p == 0 ? 0 : loc.in_var);
			max_out_vars[p][q]=(q == 0 ? 0 : loc.out_var);
		}
	}

	cout << "\tGap";
	for (i=0; i<out_edge_idxs.size(); i++)
	{
		cout << "\t";
		prm->print_edge_label(out_edge_idxs[i],max_out_vars[0][i+1]);
	//	cout << out_edge_idxs[i];
	}
	cout << endl;
	cout << "Gap";
	for (i=0; i<max_scores[0].size(); i++)
		cout << "\t" << max_scores[0][i];
	cout << endl;
	for (i=1; i<max_scores.size(); i++)
	{
		prm->print_edge_label(in_edge_idxs[i-1],max_in_vars[i][0]);
	//	cout << in_edge_idxs[i-1];
		int j;
		for (j=0; j<max_scores[i].size(); j++)
			cout << "\t" << max_scores[i][j];
		cout << endl;
	}
	cout << endl;
}

void PrmGraph::print_edge_label(int edge_idx, int var_idx, ostream& os) const
{
	const vector<string>& aa2label= config->get_aa2label();
	int *var_ptr=multi_edges[edge_idx].variant_ptrs[var_idx];
	const int num_aa = *var_ptr++;
	int i;
	for (i=0; i<num_aa; i++)
		os << aa2label[var_ptr[i]];
}










