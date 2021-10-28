#include "DeNovoSolutions.h"
#include "PeptideRankScorer.h"
#include "CumulativeSeqProb.h"
#include "SpectraList.h"
#include "PepNovo_auxfun.h"


mass_t SeqPathKey::tolerance;

int num_messages=0;

/**************************************************************************
	Adds a new seqpath to the existing vector.
	If vector is at maximal size, can replace the lowest confidence path.
	If a path exists with the same seq and n,c masses but lower confidence
	score, replaces it.

  return values:
	0 - score not high enough to be intered in heap
	1 - solution already exists with higher score
	2 - added path with no need to remove something in the heap
	3 - added path but removed a different one from the heap
***************************************************************************/
int SeqPathHeap::add_path(SeqPath& new_path, bool verbose)
{
	if (new_path.sort_key<min_value && paths.size() >= max_size)
		return 0;

	SeqPathKey new_key(new_path);
	set<SeqPathKey>::iterator it=path_keys.find(new_key);
	if (it != path_keys.end())
	{
		if (verbose)
		{
			cout << "Same: " << new_path.seq_str << " " << new_path.n_term_mass << " " << new_path.pm_with_19 << " " << new_path.path_score << endl;
			cout << "Same: " << it->pep_str << " " << it->n_mass << " " << it->sort_key << endl;
		}

		if (it->sort_key < new_path.sort_key)
		{
			if (verbose)
			{
				cout << "Replaced " << setprecision(6) << it->sort_key << " with " << new_path.sort_key << endl;
				cout << it->n_mass << " : " << new_path.n_term_mass << endl;
				cout << it->pep_str << " : " << new_path.seq_str << endl;
			}
			new_key.path_pos_idx = it->path_pos_idx;
			paths[it->path_pos_idx]=new_path;
			path_keys.erase(it);
			path_keys.insert(new_key);
		}

		// no update of the heap value
		return 1;
	}

	// add path, no need to remove
	if (paths.size() < max_size)
	{
		new_key.path_pos_idx = paths.size();
		paths.push_back(new_path);
		path_keys.insert(new_key);

		min_score_heap.push_back(idx_score_pair(new_key));
		push_heap(min_score_heap.begin(),min_score_heap.end());

		min_idx = min_score_heap[0].path_pos_idx;
		min_value = min_score_heap[0].sort_key;
		return 2;
	}

	// remove lowest score and then add
	new_key.path_pos_idx = min_idx;
	SeqPathKey remove_key(paths[min_idx]);
	it=path_keys.find(remove_key);
	if (it == path_keys.end())
	{
	//	if (num_messages++<10)
	//		cout << "Warning: couldn't find remove key!" << endl;
	//	return;
	}

	paths[min_idx] = new_path;

	if (it != path_keys.end())
		path_keys.erase(it);

	path_keys.insert(new_key);

	pop_heap(min_score_heap.begin(),min_score_heap.end());
	min_score_heap.pop_back();

	min_score_heap.push_back(idx_score_pair(new_key));
	push_heap(min_score_heap.begin(),min_score_heap.end());

	min_idx = min_score_heap[0].path_pos_idx;
	min_value = min_score_heap[0].sort_key;

	return 3;
}



void print_prm_graph_scores(AllScoreModels *model, Spectrum *spec, 
					 mass_t pm_with_19, int charge, bool prm_norm)
{
	const bool verbose = false;
	PrmGraph prm;
	const Config *config= model->get_config();

	spec->setCharge(charge);
	model->init_model_for_scoring_spectrum(spec);

	
	if (verbose)
	{
		if (! strcmp(model->getModelName().c_str(),"ETD"))
		{
			vector<string> labels;
			labels.push_back("c");
			labels.push_back("z");
			spec->print_expected_fragment_peaks(labels,cout);	
		}
		else
			spec->print_expected_by();
	}

	prm.create_graph_from_spectrum(model,spec,pm_with_19);

	if (prm_norm)
		model->getPrmNodeScoreModel().normalize_prm_scores(prm);

//	spec->print_expected_by();
	prm.print_only_scores();
	cout << endl;

//	prm.print();
}


/**************************************************************************
	Wrapper funciton that generates the desired solutions.
	Combines both local and global solutions (similar to the PepNovoTag
	and LocalTag solutions).
***************************************************************************/
bool generate_denovo_solutions(PrmGraph*& prm,
							   AllScoreModels *model, 
							   Spectrum *spec,
							   bool denovo_mode,
							   mass_t pm_with_19,
							   int charge,
							   int num_sols, 
							   int min_length, 
							   int max_length,
							   score_t min_score_needed,
							   vector<SeqPath>& solutions,
							   bool only_complete,
							   bool only_from_graph_containing_true_pep,
							   bool need_to_create_PrmGraph)
{
	DeNovoDp dndp;
	const Config *config= model->get_config();
	const score_t forbidden_pair_penalty = config->get_forbidden_pair_penalty();
	

	// shorten the generated sequences to improve the accuracy
	if (max_length>10)
	{
		if (pm_with_19>1800)
		{
			if (min_length>=10)
			{
				min_length=8;
				max_length=14;
			}
			else
				max_length=14;
		}

		if (pm_with_19>2200)
		{
			if (min_length>=10)
			{
				min_length=8;
				max_length=12;
			}
			else
				max_length=12;
		}
	}



	if (! prm)
		prm = new PrmGraph;

	if (need_to_create_PrmGraph)
	{
		spec->setCharge(charge);
		model->init_model_for_scoring_spectrum(spec);
		prm->create_graph_from_spectrum(model,spec,pm_with_19,charge);
		model->score_graph_edges(*prm);
	}

//	prm->print_with_multi_edges();
//	exit(0);

	if (only_from_graph_containing_true_pep)
	{
		const Peptide& true_pep = spec->getPeptide();
		if (true_pep.get_num_aas()<3)
		{
			cout << "Error: looking for graphs containing peptide path, but no peptide was given!" << endl;
			exit(1);
		}

		SeqPath best_path = prm->get_longest_subpath(true_pep,0);
		if (best_path.get_num_aa() < true_pep.get_num_aas())
			return false;
	}

	dndp.fill_dp_table(prm,config->get_forbidden_pair_penalty());
	SeqPath best = prm->get_longest_subpath(spec->getPeptide(),0);
//	best.print_full(config);



	if (1)
	{
		SeqPathHeap denovo_path_heap;
		vector<SeqPath> seq_paths;
		vector<MultiPath> multi_paths;

		//	generate global tags from parsing de novo paths
		denovo_path_heap.init(num_sols,config->getTolerance());

		int num_multi_paths=num_sols;
		if (num_sols<2000)
			num_multi_paths=2000;

		dndp.get_top_scoring_antisymetric_paths_with_length_limits(multi_paths, 
			num_multi_paths, min_length , max_length, forbidden_pair_penalty, min_score_needed, 
			denovo_mode, only_complete);

		prm->expand_all_multi_paths(model, multi_paths, seq_paths, forbidden_pair_penalty, num_sols);

		int i;
		for (i=0; i<seq_paths.size(); i++)
		{
			seq_paths[i].pm_with_19 = pm_with_19;
			seq_paths[i].charge = charge;
			seq_paths[i].make_seq_str(config);
			seq_paths[i].sort_key = seq_paths[i].path_score;

			// if peptide is whole, make sure its mass is within the pm tolerance
			if (seq_paths[i].n_term_mass<1.0 && seq_paths[i].c_term_mass + 30.0 > seq_paths[i].pm_with_19)
			{
				const mass_t pep_mass = seq_paths[i].calculate_peptide_mass(config) + MASS_OHHH;
				const mass_t delta = fabs(pep_mass - spec->get_org_pm_with_19());
			//	cout << i << "\t" << pep_mass << "\t" << spec->get_org_pm_with_19() << "\t" << delta << endl;
			//	if (delta > config->get_pm_tolerance())
			//		continue;
			}
			denovo_path_heap.add_path(seq_paths[i]);


		}

		sort(denovo_path_heap.paths.begin(),denovo_path_heap.paths.end(),comp_SeqPath_path_score);

		solutions = denovo_path_heap.paths;
		for (i=0; i<solutions.size(); i++)
		{
			solutions[i].delta_score = solutions[0].path_score - solutions[i].path_score;
			solutions[i].make_seq_str(config);
		}
	}

	

//	cout << "SOLS: " << solutions.size() << " " << min_score_needed << endl;
	return true;
}

/**************************************************************************
	Wrapper funciton that generates the desired solutions.
	Combines both local and global solutions (similar to the PepNovoTag
	and LocalTag solutions).
***************************************************************************/
bool generate_denovo_solutions_with_good_start_end_idxs(
							   PrmGraph*& prm,
							   AllScoreModels *model, 
							   Spectrum *spec,
							   bool denovo_mode,
							   mass_t pm_with_19,
							   int  charge,
							   int num_sols, 
							   int min_length, 
							   int max_length,
							   vector<SeqPath>& solutions)
{
	DeNovoDp dndp;
	
	const Config *config= model->get_config();
	const score_t forbidden_pair_penalty = config->get_forbidden_pair_penalty();
	
	if (charge>3 || charge<1)
	{
		solutions.clear();
		return false;
	}

	

	if (! prm)
		prm = new PrmGraph;

	spec->setCharge(charge);
	model->init_model_for_scoring_spectrum(spec);
	prm->create_graph_from_spectrum(model,spec,pm_with_19);
	model->score_graph_edges(*prm);

	dndp.fill_dp_table(prm,forbidden_pair_penalty);

	// find which nodes are good start/end points
	const vector<Node>& nodes = prm->get_nodes();
	vector<bool> good_idx_inds;

	good_idx_inds.resize(nodes.size(),false);
	int m_idx=0,n_idx=0;
	vector<mass_t> good_masses;
	spec->getPeptide().calc_expected_breakage_masses(config,good_masses);
	const int num_nodes = nodes.size(), num_masses = good_masses.size();

	good_idx_inds[0]=true;
	good_idx_inds[num_nodes-1]=true;
	while (n_idx<num_nodes && m_idx<num_masses)
	{
		const mass_t m_mass = good_masses[m_idx];
		const mass_t n_mass = nodes[n_idx].mass;

		if (fabs(m_mass-n_mass)<1.25)
			good_idx_inds[n_idx]=true;

		if (m_mass<n_mass-2)
		{
			m_idx++;
		}
		else
			n_idx++;
	}

	SeqPathHeap denovo_path_heap;
	vector<MultiPath> multi_paths;

	//	generate global tags from parsing de novo paths
	denovo_path_heap.init(num_sols,config->getTolerance());

	dndp.get_top_scoring_antisymetric_paths_with_specified_start_end_idxs(
		good_idx_inds, multi_paths, num_sols, min_length, max_length, config->get_forbidden_pair_penalty());

	vector<SeqPath> seq_paths;
	prm->expand_all_multi_paths(model, multi_paths,seq_paths, forbidden_pair_penalty, num_sols);

	// parse seq paths

	int i;
	for (i=0; i<seq_paths.size(); i++)
	{
		seq_paths[i].sort_key = seq_paths[i].path_score;
		denovo_path_heap.add_path(seq_paths[i]);
	}
	
	sort(denovo_path_heap.paths.begin(),denovo_path_heap.paths.end(),comp_SeqPath_path_score);

	solutions = denovo_path_heap.paths;
	
	return true;
}





/***************************************************************************
	Wrapper function that generates several solutions according to different 
	precursor masses.
****************************************************************************/
void generate_denovo_solutions_from_several_pms(vector<PrmGraph *>& prm_ptrs,
												AllScoreModels *model,
												Spectrum *spec,
												bool ind_denovo_mode,
												int num_sols, 
												int min_length, 
												int max_length,
												vector<mass_t>&  different_pms_with_19,
												vector<int>& charges,
												vector<SeqPath>& best_solutions,
												bool ind_only_complete)
{
	static SeqPathHeap seq_heap;
	best_solutions.clear();

	seq_heap.init(num_sols, model->get_config()->getTolerance());


	time_t start_time = time(NULL);
	int i;
	for (i=0; i<different_pms_with_19.size(); i++)
	{
		const score_t min_score_needed = ( (! ind_denovo_mode || i==0 || seq_heap.min_score_heap.size()<num_sols) ?
									 NEG_INF :  seq_heap.min_value );
	
		vector<SeqPath> sols;
		generate_denovo_solutions(prm_ptrs[i], model,spec, ind_denovo_mode,
			different_pms_with_19[i],charges[i],num_sols,min_length, max_length, 
			min_score_needed, sols, ind_only_complete);

		// removes duplicated
		int j;
		for (j=0; j<sols.size(); j++)
		{
			seq_heap.add_path(sols[j]);
		}

		// impose time limit
		time_t curr_time = time(NULL);
		if (curr_time-start_time>7.0)
			break;
	}

	best_solutions = seq_heap.paths;


	if (ind_denovo_mode)
	{
		sort(best_solutions.begin(),best_solutions.end(),comp_SeqPath_path_score);
	}
	else
		sort(best_solutions.begin(),best_solutions.end(),comp_SeqPath_sort_key);

	for (i=0; i<best_solutions.size(); i++)
		best_solutions[i].org_rank=i;
}

void generate_denovo_solutions_from_several_pms_with_good_start_end_idxs(
							   vector<PrmGraph *>& prm_ptrs,
							   AllScoreModels *model,
							   Spectrum *spec,
							   bool ind_denovo_mode,
							   int num_sols, int min_length, int max_length,
							   vector<mass_t>& different_pms_with_19,
							   vector<int>& charges,
							   vector<SeqPath>& best_solutions)
{
	static SeqPathHeap seq_heap;
	best_solutions.clear();

	seq_heap.init(num_sols, model->get_config()->getTolerance());

	int i;
	for (i=0; i<different_pms_with_19.size(); i++)
	{
		vector<SeqPath> sols;
		generate_denovo_solutions_with_good_start_end_idxs(prm_ptrs[i], model,spec, ind_denovo_mode,
			different_pms_with_19[i],charges[i],num_sols,min_length,max_length,sols);

		// removes duplicated
		int j;
		for (j=0; j<sols.size(); j++)
			seq_heap.add_path(sols[j]);
	}

	best_solutions = seq_heap.paths;

	if (ind_denovo_mode)
	{
		sort(best_solutions.begin(),best_solutions.end(),comp_SeqPath_path_score);
	}
	else
		sort(best_solutions.begin(),best_solutions.end(),comp_SeqPath_sort_key);

	for (i=0; i<best_solutions.size(); i++)
		best_solutions[i].org_rank=i;
}


/***************************************************************************
Tests each tag to see how many de novo paths it is part of.
Sets the variablesi: multi_path_rank (first rank),percent5, percent20, percent_all
in each tag.
****************************************************************************/
void determine_tag_containment_ratios(vector<SeqPath>& tags, 
									  const vector<SeqPath>& paths,
									  const vector<score_pair>& denovo_scores,
									  int start_idx)
{

	int i;
	for (i=start_idx; i<tags.size(); i++)
	{
		int num_5=0, num_20=0, num_all=0, first=POS_INF;
		const vector<PathPos>& tag_positions = tags[i].positions;
		const int first_node_idx = tag_positions[0].node_idx;

	/*	cout << "Tag: ";
		int q;
		for (q=0; q<tag_positions.size(); q++)
			cout << tag_positions[q].node_idx << " ";
		cout << endl;*/

		int j;
		for (j=0; j<denovo_scores.size(); j++)
		{
			const SeqPath& denovo_path = paths[denovo_scores[j].idx];
			const vector<PathPos>& path_positions = denovo_path.positions;
			const int last_pos = path_positions.size()-1;
		
	/*		cout << "Path " << j << " : ";
			int q;
			for (q=0; q<path_positions.size(); q++)
				cout << path_positions[q].node_idx << " ";
			cout << endl;*/

			
			int k=0;
			while (k<last_pos && path_positions[k].node_idx<first_node_idx)
				k++;

			if (path_positions[k].node_idx != first_node_idx)
				continue;

			int a;
			for (a=1; a<tag_positions.size() && a+k<path_positions.size(); a++)
				if (tag_positions[a].node_idx != path_positions[a+k].node_idx)
					break;

			if (a<tag_positions.size())
				continue;

			if (j<first)
				first=j;

			num_all++;
			if (j<20)
				num_20++;
			if (j<5)
				num_5++;
		}

		tags[i].multi_path_rank=first;
		if (paths.size()>0)
		{
			tags[i].tag_percent_top_5  = num_5*0.2;
			tags[i].tag_percent_top_20 = num_20*0.05;
			tags[i].tag_percent_all	   = num_all/(float)paths.size();
		}

	//	cout << i << "\t" << tags[i].seq_str << "\t" << tags[i].multi_path_rank <<
	//		"\t" << tags[i].tag_percent_top_5 << "\t" <<
	//		tags[i].tag_percent_top_20 << "\t" << setprecision(4) << tags[i].tag_percent_all << endl;
	}
//	exit(0);
}


/*************************************************************************
Generates tags by making a mixture of local/de novo tags
checks which tags appear in the longer denovo sequences sets the boolean
indicators in the tag seq paths
**************************************************************************
void generate_tags(vector<PrmGraph *>& prm_ptrs,
				   AllScoreModels *model,
				   BasicSpectrum& bs,
				   Spectrum *spec,
				   const vector<int>& final_num_tags,
				   int main_tag_length,				 // the length for which we parse de novo sequences
				   const vector<mass_t>& pms_with_19,
				   const vector<int>& charges,
				   vector<SeqPath>& final_tags,
				   bool use_original_num_tags,
				   int  prm_ptr_start_idx)
{
	const int max_tag_length=10;
	const double search_time_limit = 7.0;
	const int denovo_seq_heap_size = 50;
	static SeqPathHeap denovo_seq_heap;
	static vector<SeqPathHeap> tag_heaps;

	const Config *config = model->get_config();
	const mass_t tolerance = config->getTolerance();
	const int charge   = charges[0];
	const int size_idx = config->calc_size_idx(charge,pms_with_19[0]);
	
	int i;

	int min_denovo_length=6;
	for (i=0; i<final_num_tags.size(); i++)
		if (final_num_tags[i]>0)
		{
			min_denovo_length=i;
			break;
		}
	if (min_denovo_length<6)
		min_denovo_length=6;

	final_tags.clear();

	// init models, num tags, etc.
	// don't use de novo where it is time consuming and not accurate
	int num_denovo_solutions = denovo_seq_heap_size;
	bool perform_denovo_rerank=false;
	bool perform_denovo=true;
	if (charge==2 && pms_with_19[0]>=2400 || 
		charge>2 && pms_with_19[0]>2450)
	{
		perform_denovo=false;
		num_denovo_solutions=0;
	}

	PeptideRankScorer * denovo_rank_model = (PeptideRankScorer *)model->get_rank_model_ptr(1);
//	if (0 && perform_denovo && denovo_rank_model && denovo_rank_model->get_ind_part_model_was_initialized(charge,size_idx))
//	{
//		perform_denovo_rerank=true;
//		num_denovo_solutions = denovo_seq_heap_size * 3;
//	}

	vector<bool> perform_tag_reranks;
	vector<int> num_tags;
	vector<PeptideRankScorer *> tag_rank_models;

	perform_tag_reranks.resize(max_tag_length,false);
	num_tags.resize(max_tag_length,0);
	tag_rank_models.resize(max_tag_length,NULL);

	int tag_round=0;
	int tag_length;
	for (tag_length=0; tag_length<final_num_tags.size() && tag_length<max_tag_length; tag_length++)
	{
		if (final_num_tags[tag_length]<=0)
			continue;
	
		num_tags[tag_length]=final_num_tags[tag_length];
		tag_rank_models[tag_length] = (PeptideRankScorer *)model->get_rank_tag_model_ptr(tag_length);
		if (tag_rank_models[tag_length] && tag_rank_models[tag_length]->get_ind_part_model_was_initialized(charge,size_idx))
		{	
			perform_tag_reranks[tag_length]=true;
			if (! use_original_num_tags)
			{
				num_tags[tag_length] *= (4+2*tag_round);	// add more to larger lengths since they are covered by shorter tags
				num_tags[tag_length] += 10;
			}
		}
		tag_round++;						
	}

	denovo_seq_heap.init(denovo_seq_heap_size, tolerance);
	tag_heaps.resize(max_tag_length);
	for (tag_length=0; tag_length<max_tag_length; tag_length++)
		if (num_tags[tag_length]>0)
			tag_heaps[tag_length].init(num_tags[tag_length],tolerance);
	
	const clock_t start_t = clock();

	// collect de novo sequences and generate local tag
	for (i=0; i<pms_with_19.size(); i++)
	{
		// ignore differnet charges, should not be here
		if (charges[i] != charges[0])
			continue;

		const clock_t end_t = clock();
		const double run_time = (end_t - start_t)/(double)CLOCKS_PER_SEC;

		// limit the time spent here if searches run too long
		if (i>0 && run_time>search_time_limit)
		{
			break;
		}

		// First generate de novo solutions 
		const score_t min_seq_score_needed = ( (i==0 || denovo_seq_heap.min_score_heap.size()<num_denovo_solutions) ?
											 NEG_INF :  denovo_seq_heap.min_value );

		int max_length = 13;
		if (pms_with_19[i]>1800)
			max_length = 10;
		if (charges[i]>2 || pms_with_19[i]>2200)
			max_length = 9;

		if (perform_denovo)
		{
			vector<SeqPath> denovo_sols;


			generate_denovo_solutions(prm_ptrs[i+prm_ptr_start_idx],  model, spec, true,
				pms_with_19[i], charges[i], num_denovo_solutions, min_denovo_length, max_length,
				min_seq_score_needed, denovo_sols, false);

			if (i==0 && bs.ssf->peptide.get_num_aas()>2)
			{
			//	cout << endl;
			//	SeqPath best = prm_ptrs[0+prm_ptr_start_idx]->get_longest_subpath(bs.ssf->peptide,0);
			//	spec->print_expected_by();
			//	best.print_full(config);
			//	prm_ptrs[0]->print();
			//	exit(0);

			}

			// just remove duplicates
			int j;
			for (j=0; j<denovo_sols.size(); j++)
				denovo_seq_heap.add_path(denovo_sols[j]);
		}
		
		// generate local tags
		int tag_length;
		for (tag_length=0; tag_length<max_tag_length; tag_length++)
		{
			if (num_tags[tag_length]>0)
			{
				const score_t min_tag_score_needed = ( (i==0 || tag_heaps[tag_length].min_score_heap.size()<num_tags[tag_length]) ?
													   NEG_INF :  tag_heaps[tag_length].min_value );		
				vector<SeqPath> tag_sols;
				generate_denovo_solutions(prm_ptrs[i+prm_ptr_start_idx],  model, spec, false,
					pms_with_19[i], charges[i], num_tags[tag_length], tag_length, tag_length, 
					min_tag_score_needed, tag_sols, false, false, (!perform_denovo));

				int j;
				for (j=0; j<tag_sols.size(); j++)
				{
					tag_heaps[tag_length].add_path(tag_sols[j]);
					if (tag_sols[j].get_num_aa() != tag_length)
					{
						cout << "problem1  " << j << " " << tag_sols[j].get_num_aa() << " : " << tag_length <<endl;
					}
				}
			}
		}
	}

	// sort de novo seqs
	vector<score_pair> denovo_scores;
	denovo_scores.clear();
	if (perform_denovo)
	{
		vector<SeqPath>& denovo_seqs = denovo_seq_heap.paths;
		if (0) // don't rerank sequences
		{
			//denovo_rank_model->score_denovo_sequences(denovo_seqs,bs.ssf,bs.peaks,bs.num_peaks,denovo_scores,size_idx);
			sort(denovo_scores.begin(),denovo_scores.end());
		}
		else
		{
			int i;
			denovo_scores.resize(denovo_seqs.size());
			for (i=0; i<denovo_scores.size(); i++)
			{
				denovo_scores[i].idx=i;
				denovo_scores[i].score=denovo_seqs[i].path_score;
			}
			sort(denovo_scores.begin(),denovo_scores.end());
		}

		// parse de novo paths into tags and add them to the tag heap
		if (num_tags[main_tag_length]>0)
		{
			int max_num_to_force = num_tags[main_tag_length]/4;
			if (max_num_to_force>15)
				max_num_to_force=15;

			int num_forced=0;
			int i;
			for (i=0; i<denovo_scores.size(); i++)
			{
				SeqPath& big_path = denovo_seqs[denovo_scores[i].idx];
				vector<SeqPath> parsed_tags;
				
				PrmGraph *prm_ptr = big_path.prm_ptr;
				big_path.multi_path_rank = i;
				prm_ptr->parse_seq_path_to_smaller_ones(big_path, main_tag_length, main_tag_length, parsed_tags);

				int j;
				for (j=0; j<parsed_tags.size(); j++)
				{
					int add_return_value = tag_heaps[main_tag_length].add_path(parsed_tags[j]);

					// if the tag came from one of the top 5 sequences and it was rejected
					// because of a low score, we'll force it in by adding a high score +10000
					// that will later be removed
					if (add_return_value == 0 && 
						num_forced<max_num_to_force && 
						parsed_tags[j].multi_path_rank<5)
					{
						parsed_tags[j].sort_key += 10000.0;
						num_forced++;

					//	int ret_val = tag_heaps[main_tag_idx].add_path(parsed_tags[j]);
					//	cout << ret_val << " " << num_forced << " Added " << parsed_tags[j].multi_path_rank << "\t" << 
					//		parsed_tags[j].seq_str << "\t" << parsed_tags[j].path_score << "\t"
					//		<< parsed_tags[j].sort_key <<"\t" << tag_heaps[main_tag_idx].min_value << endl;
					}

					if (parsed_tags[j].get_num_aa() != main_tag_length)
					{
						cout << "Problem2 " << j << " " << main_tag_length << " : " << parsed_tags[j].get_num_aa() << endl;
						exit(1);
					}
				}
			}

			// since the heap will not be used as a heap anymore (it is now just a container)
			// we can go back and remove the bonus scores that were added
			for (i=0; i<tag_heaps[main_tag_length].paths.size(); i++)
			{
				if (tag_heaps[main_tag_length].paths[i].sort_key>8000.0)
				{
					tag_heaps[main_tag_length].paths[i].sort_key-=10000.0;
				}
			}
		}
	}

	// sort tag heaps
	tag_round=0;
	for (tag_length=0; tag_length<max_tag_length; tag_length++)
	{	
		if (num_tags[tag_length]<=0)
			continue;

		vector<SeqPath>& bigger_tags = tag_heaps[tag_length].paths;
		sort(bigger_tags.begin(),bigger_tags.end(),comp_SeqPath_path_score);
		
		const float top_score = (bigger_tags.size()>0 ? bigger_tags[0].path_score : NEG_INF);
		int i;
		for (i=0; i<bigger_tags.size(); i++)
		{
			bigger_tags[i].org_rank=i;
			bigger_tags[i].delta_score = top_score - bigger_tags[i].path_score;
		}

		// check that tags are not covered by previous shorter ones
		if (tag_round>0)
		{
			
			const vector<SeqPath>& smaller_tags = final_tags;
			int i;
			for (i=0; i<bigger_tags.size(); i++)
			{
				const vector<PathPos>& big_positions = bigger_tags[i].positions;
				const int bigger_length = big_positions.size()-1;

				int j;
				for (j=0; j<smaller_tags.size(); j++)
				{
					const vector<PathPos>& small_positions = smaller_tags[j].positions;
					const int smaller_length = small_positions.size()-1;
					const int length_diff  = bigger_length - smaller_length;

					int k;
					for (k=0; k<length_diff; k++)
					{
						if (small_positions[0].aa == big_positions[k].aa)
						{
							int a;
							for (a=1; a<smaller_length; a++)
								if (small_positions[a].aa != big_positions[a+k].aa)
									break;
							if (a==smaller_length)
								break;
						}
					}
					if (k<length_diff)
						break;
				}

				// this tag is covered remove it
				if (j<smaller_tags.size())
				{
					bigger_tags[i]=bigger_tags[bigger_tags.size()-1];
					bigger_tags.pop_back();
					i--;
				}
			}
			
		}
	
		if (perform_denovo)
			determine_tag_containment_ratios(bigger_tags,denovo_seq_heap.paths, denovo_scores,0);

		// sort tags
		vector<score_pair> tag_scores;
		if (perform_tag_reranks[tag_length])
		{
		//	tag_rank_models[tag_length]->score_tag_sequences(bigger_tags,bs.ssf,bs.peaks,bs.num_peaks,tag_scores,size_idx);
			sort(tag_scores.begin(),tag_scores.end());
			int i;
			for (i=0; i<tag_scores.size(); i++)
				bigger_tags[tag_scores[i].idx].rerank_score = tag_scores[i].score;
		}
		else
		{
			int i;
			tag_scores.resize(bigger_tags.size());
			for (i=0; i<tag_scores.size(); i++)
			{
				tag_scores[i].idx=i;
				tag_scores[i].score=bigger_tags[i].path_score;
			}
		}

		while (tag_scores.size()>final_num_tags[tag_length])
			tag_scores.pop_back();

		for (i=0; i<tag_scores.size(); i++)
		{
			final_tags.push_back(bigger_tags[tag_scores[i].idx]);
		}
		tag_round++;	
	}
}*/



/*************************************************************************
Generates tags by making a mixture of local/de novo tags
checks which tags appear in the longer denovo sequences sets the boolean
indicators in the tag seq paths
**************************************************************************/
void generate_tags(vector<PrmGraph *>& prm_ptrs,
				   AllScoreModels* model,
				   AnnotatedSpectrum& as,
				   const vector<int>& final_num_tags,
				   int main_tag_length,				 // the length for which we parse de novo sequences
				   const vector<mass_t>& pms_with_19,
				   const vector<int>& charges,
				   vector<SeqPath>& final_tags,
				   bool use_original_num_tags,
				   int  prm_ptr_start_idx)
{
	const int max_tag_length=10;
	const double search_time_limit = 7.0;
	const int denovo_seq_heap_size = 50;
	static SeqPathHeap denovo_seq_heap;
	static vector<SeqPathHeap> tag_heaps;

	const Config *config = model->get_config();
	const mass_t tolerance = config->getTolerance();
	const int charge   = charges[0];
	const int size_idx = config->calc_size_idx(charge,pms_with_19[0]);

	int min_denovo_length=6;
	for (size_t i=0; i<final_num_tags.size(); i++)
		if (final_num_tags[i]>0)
		{
			min_denovo_length=i;
			break;
		}
	if (min_denovo_length<6)
		min_denovo_length=6;

	final_tags.clear();

	// init models, num tags, etc.
	// don't use de novo where it is time consuming and not accurate
	int num_denovo_solutions = denovo_seq_heap_size;
	bool perform_denovo_rerank=false;
	bool perform_denovo=true;
	if (charge==2 && pms_with_19[0]>=2400 || 
		charge>2 && pms_with_19[0]>2450)
	{
		perform_denovo=false;
		num_denovo_solutions=0;
	}

	PeptideRankScorer * denovo_rank_model = (PeptideRankScorer *)model->get_rank_model_ptr(1);


	vector<bool> perform_tag_reranks;
	vector<int> num_tags;
	vector<PeptideRankScorer *> tag_rank_models;

	perform_tag_reranks.resize(max_tag_length,false);
	num_tags.resize(max_tag_length,0);
	tag_rank_models.resize(max_tag_length,NULL);

	int tag_round=0;
	int tag_length;
	for (tag_length=0; tag_length<final_num_tags.size() && tag_length<max_tag_length; tag_length++)
	{
		if (final_num_tags[tag_length]<=0)
			continue;
	
		num_tags[tag_length]=final_num_tags[tag_length];
		tag_rank_models[tag_length] = (PeptideRankScorer *)model->get_rank_tag_model_ptr(tag_length);
		if (tag_rank_models[tag_length] && tag_rank_models[tag_length]->get_ind_part_model_was_initialized(charge,size_idx))
		{	
			perform_tag_reranks[tag_length]=true;
			if (! use_original_num_tags)
			{
				num_tags[tag_length] *= (4+2*tag_round);	// add more to larger lengths since they are covered by shorter tags
				num_tags[tag_length] += 10;
			}
		}
		tag_round++;						
	}

	denovo_seq_heap.init(denovo_seq_heap_size, tolerance);
	tag_heaps.resize(max_tag_length);
	for (tag_length=0; tag_length<max_tag_length; tag_length++)
		if (num_tags[tag_length]>0)
			tag_heaps[tag_length].init(num_tags[tag_length],tolerance);
	
	const clock_t start_t = clock();

	// collect de novo sequences and generate local tag
	for (size_t i=0; i<pms_with_19.size(); i++)
	{
		// ignore differnet charges, should not be here
		if (charges[i] != charges[0])
			continue;

		const clock_t end_t = clock();
		const double run_time = (end_t - start_t)/(double)CLOCKS_PER_SEC;

		// limit the time spent here if searches run too long
		if (i>0 && run_time>search_time_limit)
		{
			break;
		}

		// First generate de novo solutions 
		const score_t min_seq_score_needed = ( (i==0 || denovo_seq_heap.min_score_heap.size()<num_denovo_solutions) ?
											 NEG_INF :  denovo_seq_heap.min_value );

		int max_length = 13;
		if (pms_with_19[i]>1800)
			max_length = 10;
		if (charges[i]>2 || pms_with_19[i]>2200)
			max_length = 9;

		if (perform_denovo)
		{
			vector<SeqPath> denovo_sols;


			generate_denovo_solutions(prm_ptrs[i+prm_ptr_start_idx],  model, &as, true,
				pms_with_19[i], charges[i], num_denovo_solutions, min_denovo_length, max_length,
				min_seq_score_needed, denovo_sols, false);

			// just remove duplicates
			int j;
			for (j=0; j<denovo_sols.size(); j++)
				denovo_seq_heap.add_path(denovo_sols[j]);
		}
		
		// generate local tags
		int tag_length;
		for (tag_length=0; tag_length<max_tag_length; tag_length++)
		{
			if (num_tags[tag_length]>0)
			{
				const score_t min_tag_score_needed = ( (i==0 || tag_heaps[tag_length].min_score_heap.size()<num_tags[tag_length]) ?
													   NEG_INF :  tag_heaps[tag_length].min_value );		
				vector<SeqPath> tag_sols;
				generate_denovo_solutions(prm_ptrs[i+prm_ptr_start_idx],  model, &as, false,
					pms_with_19[i], charges[i], num_tags[tag_length], tag_length, tag_length, 
					min_tag_score_needed, tag_sols, false, false, (!perform_denovo));

				int j;
				for (j=0; j<tag_sols.size(); j++)
				{
					tag_heaps[tag_length].add_path(tag_sols[j]);
					if (tag_sols[j].get_num_aa() != tag_length)
					{
						cout << "problem1  " << j << " " << tag_sols[j].get_num_aa() << " : " << tag_length <<endl;
					}
				}
			}
		}
	}

	// sort de novo seqs
	vector<score_pair> denovo_scores;
	denovo_scores.clear();
	if (perform_denovo)
	{
		vector<SeqPath>& denovo_seqs = denovo_seq_heap.paths;
		denovo_scores.resize(denovo_seqs.size());
		for (size_t i=0; i<denovo_scores.size(); i++)
		{
			denovo_scores[i].idx=i;
			denovo_scores[i].score=denovo_seqs[i].path_score;
		}
		sort(denovo_scores.begin(),denovo_scores.end());
		
		// parse de novo paths into tags and add them to the tag heap
		if (num_tags[main_tag_length]>0)
		{
			int max_num_to_force = num_tags[main_tag_length]/4;
			if (max_num_to_force>15)
				max_num_to_force=15;

			int num_forced=0;
			for (size_t i=0; i<denovo_scores.size(); i++)
			{
				SeqPath& big_path = denovo_seqs[denovo_scores[i].idx];
				vector<SeqPath> parsed_tags;
				
				PrmGraph *prm_ptr = big_path.prm_ptr;
				big_path.multi_path_rank = i;
				prm_ptr->parse_seq_path_to_smaller_ones(big_path, main_tag_length, main_tag_length, parsed_tags);

				int j;
				for (j=0; j<parsed_tags.size(); j++)
				{
					int add_return_value = tag_heaps[main_tag_length].add_path(parsed_tags[j]);

					// if the tag came from one of the top 5 sequences and it was rejected
					// because of a low score, we'll force it in by adding a high score +10000
					// that will later be removed
					if (add_return_value == 0 && 
						num_forced<max_num_to_force && 
						parsed_tags[j].multi_path_rank<5)
					{
						parsed_tags[j].sort_key += 10000.0;
						num_forced++;
					}

					if (parsed_tags[j].get_num_aa() != main_tag_length)
					{
						cout << "Problem2 " << j << " " << main_tag_length << " : " << parsed_tags[j].get_num_aa() << endl;
						exit(1);
					}
				}
			}

			// since the heap will not be used as a heap anymore (it is now just a container)
			// we can go back and remove the bonus scores that were added
			for (size_t i=0; i<tag_heaps[main_tag_length].paths.size(); i++)
			{
				if (tag_heaps[main_tag_length].paths[i].sort_key>8000.0)
				{
					tag_heaps[main_tag_length].paths[i].sort_key-=10000.0;
				}
			}
		}
	}

	// sort tag heaps
	tag_round=0;
	for (tag_length=0; tag_length<max_tag_length; tag_length++)
	{	
		if (num_tags[tag_length]<=0)
			continue;

		vector<SeqPath>& bigger_tags = tag_heaps[tag_length].paths;
		sort(bigger_tags.begin(),bigger_tags.end(),comp_SeqPath_path_score);
		
		const float top_score = (bigger_tags.size()>0 ? bigger_tags[0].path_score : NEG_INF);
		int i;
		for (i=0; i<bigger_tags.size(); i++)
		{
			bigger_tags[i].org_rank=i;
			bigger_tags[i].delta_score = top_score - bigger_tags[i].path_score;
		}

		// check that tags are not covered by previous shorter ones
		if (tag_round>0)
		{
			
			const vector<SeqPath>& smaller_tags = final_tags;
			int i;
			for (i=0; i<bigger_tags.size(); i++)
			{
				const vector<PathPos>& big_positions = bigger_tags[i].positions;
				const int bigger_length = big_positions.size()-1;

				int j;
				for (j=0; j<smaller_tags.size(); j++)
				{
					const vector<PathPos>& small_positions = smaller_tags[j].positions;
					const int smaller_length = small_positions.size()-1;
					const int length_diff  = bigger_length - smaller_length;

					int k;
					for (k=0; k<length_diff; k++)
					{
						if (small_positions[0].aa == big_positions[k].aa)
						{
							int a;
							for (a=1; a<smaller_length; a++)
								if (small_positions[a].aa != big_positions[a+k].aa)
									break;
							if (a==smaller_length)
								break;
						}
					}
					if (k<length_diff)
						break;
				}

				// this tag is covered remove it
				if (j<smaller_tags.size())
				{
					bigger_tags[i]=bigger_tags[bigger_tags.size()-1];
					bigger_tags.pop_back();
					i--;
				}
			}
			
		}
	
		if (perform_denovo)
			determine_tag_containment_ratios(bigger_tags,denovo_seq_heap.paths, denovo_scores,0);

		// sort tags
		vector<score_pair> tag_scores;
		if (perform_tag_reranks[tag_length])
		{
			tag_rank_models[tag_length]->scoreTagSequences(bigger_tags, as, tag_scores, size_idx);
			sort(tag_scores.begin(),tag_scores.end());
			int i;
			for (i=0; i<tag_scores.size(); i++)
				bigger_tags[tag_scores[i].idx].rerank_score = tag_scores[i].score;
		}
		else
		{
			int i;
			tag_scores.resize(bigger_tags.size());
			for (i=0; i<tag_scores.size(); i++)
			{
				tag_scores[i].idx=i;
				tag_scores[i].score=bigger_tags[i].path_score;
			}
		}

		while (tag_scores.size()>final_num_tags[tag_length])
			tag_scores.pop_back();

		for (i=0; i<tag_scores.size(); i++)
		{
			final_tags.push_back(bigger_tags[tag_scores[i].idx]);
		}
		tag_round++;	
	}
}










/*
void output_denovo_solutions(SingleSpectrumFile *ssf, Config *config, ostream& out_stream,
							 const vector<SeqPath>& solutions, int max_num_sols)
{
	if (max_num_sols<0)
		max_num_sols=solutions.size();

	ssf->print_ssf_stats(config,out_stream);

	if (solutions.size() == 0)
	{
		out_stream << "No solutions found." << endl;
	}
	else 
	{
		out_stream << "#Index\tScore\tN-Gap\tC-Gap\t[M+H]\tCharge\tSequence" << endl;
		int i; 	
		for (i=0; i<solutions.size() && i<max_num_sols; i++) 
		{
			mass_t c_gap=solutions[i].pm_with_19 - solutions[i].c_term_mass;
			if (c_gap<24.0)
				c_gap = 0;

			out_stream << setprecision(3) << fixed << i << "\t";
			out_stream << solutions[i].path_score << "\t";
			out_stream << solutions[i].n_term_mass << "\t";
			out_stream << c_gap << "\t";
			out_stream << solutions[i].pm_with_19 << "\t";
			out_stream << solutions[i].charge << "\t";
			out_stream << solutions[i].seq_str;
			out_stream << endl;
		}
	}
	out_stream << endl;
}
*/

/*
void output_tag_solutions(SingleSpectrumFile *ssf, 
						  const Config *config, ostream& out_stream,
						  const vector<SeqPath>& solutions,
						  bool	output_aa_probs)
{
	ssf->print_ssf_stats(config,out_stream);

	if (solutions.size() == 0)
	{
		out_stream << "No solutions found." << endl;
	}
	else 
	{
		out_stream << "#Index\tRnkScr\tPnvScr\tN-Gap\tC-Gap\t[M+H]\tCharge\tSequence" << endl;
		int i; 	
		for (i=0; i<solutions.size(); i++) 
		{
			mass_t c_gap=solutions[i].pm_with_19 - solutions[i].c_term_mass;
			if (c_gap<24.0)
				c_gap = 0;

			out_stream << setprecision(3) << fixed << i << "\t";
			out_stream << solutions[i].rerank_score << "\t";
			out_stream << solutions[i].path_score << "\t";
			out_stream << solutions[i].n_term_mass << "\t";
			out_stream << c_gap << "\t";
			out_stream << solutions[i].pm_with_19 << "\t";
			out_stream << solutions[i].charge << "\t";
			out_stream << solutions[i].seq_str;

			if (output_aa_probs)
			{
				const vector<PathPos>& positions = solutions[i].positions;
				int j;
				for (j=0; j<positions.size()-1; j++)
				{
					out_stream << "\t" << positions[j].edge_variant_prob;
				}
			}

			out_stream << endl;
		}
	}
	out_stream << endl;
}

*/





void parse_tag_string(char *tag_string, int& main_length, vector<int>& num_tags)
{
	num_tags.clear();
	num_tags.resize(10,0);

	string ts(tag_string);

	int i;
	for (i=0; i<ts.length(); i++)
		if (ts[i]==':')
			ts[i]=' ';

	istringstream iss(ts);

	int len=0;
	while (iss>>len)
	{
		if (len<2 || len>10)
		{
			cout << "Error: bad string for tag lengths, should be something like \"4:20:5:30\""<< endl;
			exit(1);
		}

		int n=0;
		iss >> n;
		if (n<0 || n>1000)
		{
			cout << "Error: number of tags for length " << len << " should be 1-1000" << endl;
			exit(1);
		}
		num_tags[len]=n;
	}

	main_length=0;
	for (i=0; i<10; i++)
		if (num_tags[i]>=5)
		{
			main_length=i;
			break;
		}

	if (main_length<=0)
	{
		for (i=0; i<10; i++)
			if (num_tags[i]>0)
			{
				main_length=i;
				break;
			}
	}

	if (main_length<=0)
	{
		cout << "Error parsing tag_string: " << tag_string << endl;
		exit(1);
	}

}

/****************************************************************
The input uses a tag string which has the format
4:20:5:30 - meaning 20 tags of length 4 and 30 tags of length 5
Generate a text file containing tags for a whole spectrum.
The format is:
<scan number> <number tags>
tweak 0   (each tweak is charge pm_with_19)
..
tweak 5
tag 0 (tweak_idx score num_aas n_term mass AASeq)
..
tag n-1
*****************************************************************/
void create_tag_file_for_inspect(AllScoreModels& model,
								 char *spectrum_file,
								 char *tag_string,
								 char *tag_suffix)
{
	const int max_tweak_charge =3;
	const Config *config = model.get_config();

	vector<int> num_tags; // how many tags from each length
	int main_length;

	parse_tag_string(tag_string,main_length,num_tags);
	
	SpectraAggregator sa;
	sa.initializeFromSpectraFilePath(spectrum_file, config);
	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();

	cout << "Generating tags for: " << spectrum_file << endl;
	cout << "Scans: " << sl.getNumHeaders() << endl;

	string file_name,tag_file;
	getFileNameWithoutExtension(spectrum_file,file_name);
	tag_file = file_name + "_" + string(tag_suffix) + ".txt";
	ofstream tag_stream(tag_file.c_str());
	
	time_t start_time;
	start_time = time(NULL);

	vector<PrmGraph *> prm_ptrs;

	int too_few_peaks=0;
	int bad_charge=0;
	int no_tags=0;
	int scan_count=0;

	prm_ptrs.resize(50,NULL);

	///////////////////////////////////////////////
	int sc;
	for (sc=0; sc<sl.getNumHeaders(); sc++)
	{
		vector<int>    tweak_charges;
		vector<mass_t> tweak_pms_with_19;

		tweak_charges.resize(6,0);
		tweak_pms_with_19.resize(6,0);

		AnnotatedSpectrum as;
		as.readSpectrum(sa, sl.getSpectrumHeader(sc));
		// NP3 GOT change few peaks from 10 to 1 @@@
		if (as.getNumPeaks()<1)
		{
			too_few_peaks++;
			continue;
		}

		if ( as.getCharge() > model.get_max_score_model_charge())
		{
			bad_charge++;
			continue;
		}

		vector<SeqPath> solutions;
		vector<mass_t> pms_with_19, alt_pms_with_19;
		vector<int>    charges,     alt_charges;

		pms_with_19.clear();
		charges.clear();		
		

		// output m/z and prob values for the different charge states
		vector<PmcSqsChargeRes> pmcsqs_res;
	//	model.select_pms_and_charges(config,bs,pms_with_19,charges,&pmcsqs_res);
		model.selectPrecursorMassesAndCharges(config, as, pms_with_19, charges, &pmcsqs_res);

		cout << "Scan " << as.getHeader()->getScanNumber() << ",\tch: " << charges[0] << setprecision(1) << "\tpm19: " <<
			pms_with_19[0] << setprecision(3) << " \tSQS: " <<  pmcsqs_res[charges[0]].sqs_prob;
		if (config->get_filter_flag() && 
			! config->get_use_spectrum_charge() &&
			pmcsqs_res[charges[0]].sqs_prob<0.01)
		{
			cout << "  - skipping..." << endl;
			continue;
		}
		
		alt_charges.clear();
		alt_pms_with_19.clear();
		while (charges.size()>0 && charges[0] != charges[charges.size()-1])
		{
			alt_charges.push_back(charges[charges.size()-1]);
			charges.pop_back();
			alt_pms_with_19.push_back(pms_with_19[pms_with_19.size()-1]);
			pms_with_19.pop_back();
		}

		if (pms_with_19.size()>0)
		{
			generate_tags(prm_ptrs, &model, as, num_tags, main_length,
						  pms_with_19,charges,solutions);

			// set tweaks
			int charge = charges[0];
			if (charge>max_tweak_charge)
			{
				cout << "Error: max allowed tweak charge is " << max_tweak_charge << ", selected charge " << charge << endl;
				exit(1);
			}
			tweak_charges[2*(charge-1)]=charge;
			tweak_pms_with_19[2*(charge-1)]=pms_with_19[0];

			if (pms_with_19.size()>1)
			{
				tweak_charges[2*(charge-1)+1]=charge;
				tweak_pms_with_19[2*(charge-1)+1]=pms_with_19[1];
			}
		}

		if (alt_pms_with_19.size()>0)
		{
			vector<SeqPath> alt_solutions;

			generate_tags(prm_ptrs,&model, as, num_tags, main_length,
						  alt_pms_with_19, alt_charges, alt_solutions, false, alt_solutions.size());
			int j;
			for (j=0; j<alt_solutions.size(); j++)
				solutions.push_back(alt_solutions[j]);

			// set tweaks
			int charge = alt_charges[0];
			if (charge>max_tweak_charge)
			{
				cout << "Error: max allowed tweak charge is " << max_tweak_charge << ", selected charge " << charge << endl;
				exit(1);
			}
			tweak_charges[2*(charge-1)]=charge;
			tweak_pms_with_19[2*(charge-1)]=alt_pms_with_19[0];

			if (alt_pms_with_19.size()>1)
			{
				tweak_charges[2*(charge-1)+1]=charge;
				tweak_pms_with_19[2*(charge-1)+1]=alt_pms_with_19[1];
			}
		}
		
		
		if (solutions.size() == 0)
		{
			cout << "   - no tags" << endl;
			no_tags++;
			continue;
		}
		
		tag_stream << fixed << setprecision(3);
		tag_stream << as.getHeader()->getScanNumber() << "\t" << solutions.size() << endl;
		int i; 
		for (i=0; i<tweak_charges.size(); i++)
			tag_stream << tweak_charges[i] << "\t" << tweak_pms_with_19[i] << endl;

		for (i=0; i<solutions.size(); i++) 
		{
			const float score = ( solutions[i].rerank_score>-200 ? solutions[i].rerank_score : solutions[i].path_score);
			const int charge = solutions[i].charge;
			int best_tweak_idx=2*(charge-1);
			if (tweak_charges[2*(charge-1)+1]>0 && 
				fabs(tweak_pms_with_19[2*(charge-1)+1]-solutions[i].pm_with_19)<
				fabs(tweak_pms_with_19[2*(charge-1)]-solutions[i].pm_with_19))
				best_tweak_idx = 2*(charge-1)+1;
			
			tag_stream << best_tweak_idx << "\t" << setprecision(2) << score << "\t" << setprecision(3) << solutions[i].n_term_mass << "\t" << solutions[i].seq_str << endl;
		}

		scan_count++;

		cout << "   " << solutions.size() << " tags." << endl;

		if (sc % 200 == 0)
		{
			time_t curr_time = time(NULL);
			cout << setprecision(1) << fixed;
			cout << "done " << sc << "/" << sl.getNumHeaders() << "  " << curr_time - start_time << " secs.,  "
				<< scan_count << " scans completed." << endl;
		}
		
	}
	tag_stream.close();

	cout << "Done..." << endl;
	cout << "Tag file: " << tag_file << endl;
	cout << "Wrote tags for " << scan_count << " spectra." << endl;
	if (bad_charge>0)
		cout << "Scans with bad charges  : " << bad_charge << endl;
	if (too_few_peaks>0)
		cout << "Scans with too few peaks: " << too_few_peaks << endl;
	if (no_tags>0)
		cout << "Scans for which no tags were found: " << no_tags << endl;

}







void perform_tags_on_list_of_files(AllScoreModels& model, 
									 const vector<string>& list_vector, 
									 int file_start_idx,
									 int num_solutions, 
									 int	tag_length,
									 bool report_progress,
									 float min_filter_prob,
									 bool	output_aa_probs,
									 bool	output_cumulative_seq_probs,
									 ostream& out_stream)
{
	const Config *config = model.get_config();
	

/*	vector<int> num_tags; // how many tags from each length
	int main_length;

	
	num_tags.resize(10,0);
	num_tags[tag_length]=num_solutions;
	main_length = tag_length;

	// Quick read, get all pointers to begining of spectra

	time_t start_time, last_time;
	start_time = time(NULL);
	last_time = start_time;

	vector<PrmGraph *> prm_ptrs;
	BasicSpecReader bsr;
	int too_few_peaks=0;
	int bad_charge=0;
	int no_tags=0;
	int scan_count=0;
	int	correct_benchmark=0;
	int num_benchmark=0;
	int total_number_processed=0;
	int total_number_headers=0;


	prm_ptrs.resize(50,NULL);

	int f_idx;
	for (f_idx=0; f_idx<list_vector.size(); f_idx++)
	{
		FileManager fm;
		FileSet fs;

		const string spectrum_file = list_vector[f_idx];

		if (get_file_extension_type(spectrum_file) != MZXML)
		{
			fm.init_from_file(config,spectrum_file.c_str());
		}
		else // reads peaks 
			fm.init_and_read_single_mzXML(config,spectrum_file.c_str(),0);

		fs.select_all_files(fm);
		
		///////////////////////////////////////////////
		const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

		total_number_headers += all_ssf.size();
		int sc;
		for (sc=0; sc<all_ssf.size(); sc++)
		{
			static vector<QCPeak> peaks;
			SingleSpectrumFile *ssf = all_ssf[sc];
			if (peaks.size()<ssf->num_peaks)
			{
				int new_size = ssf->num_peaks*2;
				if (new_size<2500)
					new_size=2500;
				peaks.resize(new_size); 
			}

		
			time_t curr_time = time(NULL);
			double elapsed_time = (curr_time - last_time);
			if (report_progress && elapsed_time>5)
			{
				last_time = curr_time;
				cout << "Processing file " << f_idx+1 << "/" << list_vector.size();
				if (all_ssf.size() == 1)
				{
					cout << "  spectrum 1/1 of current file." << endl;
				}
				else
					cout << "  spectrum " << sc+1 << "/" << all_ssf.size() << 
					" of current file." << endl;
				last_time=curr_time;
			}

			const int num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0]);
			ssf->file_idx = f_idx+file_start_idx;

			// convert peak list ot a spectrum with charge (if original charge ==0)
			// the spectrum gets charge 2, but the true charge is computed from the data
		
			if (num_peaks<5)
			{
				ssf->print_ssf_stats(config,out_stream);
				out_stream << "# too few peaks..." << endl;
				continue;
			}

	
			Spectrum s;
			s.init_from_QCPeaks(config,&peaks[0],num_peaks,ssf);

			if ( ssf->charge > model.get_max_score_model_charge())
			{
				bad_charge++;
				continue;
			}

			vector<SeqPath> solutions;
			vector<mass_t> pms_with_19, alt_pms_with_19;
			vector<int>    charges,     alt_charges;

			pms_with_19.clear();
			charges.clear();		
			BasicSpectrum bs;
			bs.ssf = ssf;
			bs.peaks = &peaks[0];
			bs.num_peaks = num_peaks;

			// output m/z and prob values for the different charge states
			vector<PmcSqsChargeRes> pmcsqs_res;
			//model.select_pms_and_charges(config,bs,pms_with_19,charges,&pmcsqs_res);

			if (pmcsqs_res.size()>charges[0] && pmcsqs_res[charges[0]].sqs_prob<min_filter_prob)
			{
				ssf->print_ssf_stats(config,out_stream);
				out_stream << "# low quality, skipping: " << pmcsqs_res[charges[0]].sqs_prob << endl << endl;
				continue;
			}

			if (pms_with_19[0]<100)
			{
				ssf->print_ssf_stats(config,out_stream);
				out_stream << "# Could not process spectrum..." << endl << endl;
				continue;
			}
			
			alt_charges.clear();
			alt_pms_with_19.clear();
			while (charges.size()>0 && charges[0] != charges[charges.size()-1])
			{
				alt_charges.push_back(charges[charges.size()-1]);
				charges.pop_back();
				alt_pms_with_19.push_back(pms_with_19[pms_with_19.size()-1]);
				pms_with_19.pop_back();
			}

			if (pms_with_19.size()>0)
			{
				//generate_tags(prm_ptrs,&model,bs,&s,num_tags, main_length,
				//			  pms_with_19,charges,solutions);

			}

			if (alt_pms_with_19.size()>0)
			{
				vector<SeqPath> alt_solutions;

				//generate_tags(prm_ptrs,&model,bs,&s, num_tags, main_length,
				//			  alt_pms_with_19, alt_charges, alt_solutions, false, alt_solutions.size());
				int j;
				for (j=0; j<alt_solutions.size(); j++)
					solutions.push_back(alt_solutions[j]);
			}
			
			total_number_processed++;
			
			if (solutions.size() == 0)
			{
				ssf->print_ssf_stats(config,out_stream);
				out_stream << "   - no tags found" << endl;
				no_tags++;
				continue;
			}
			
			if (output_aa_probs)
			{
				const vector<mass_t>& aa2mass = config->get_aa2mass();
				int q;
				for (q=0; q<solutions.size(); q++)
				{
					SeqPath& seq_path = solutions[q];
					seq_path.prm_ptr->calc_amino_acid_probs(seq_path,q);
					
					mass_t exp_mass = 0;
					vector<int> aas;
					seq_path.get_amino_acids(aas);
					int t;
					for (t=0; t<aas.size(); t++)
						exp_mass+=aa2mass[aas[t]];
					
					mass_t diff = seq_path.c_term_mass - seq_path.n_term_mass - exp_mass;

					if (fabs(diff)>4)
					{
						cout << "DIFF PROBLEM!" << endl;
						seq_path.print_full(config);
						seq_path.prm_ptr->print();
						exit(0);
					}
				}
			}

			if (output_cumulative_seq_probs)
			{

			}

			output_tag_solutions(ssf,config,out_stream,solutions, output_aa_probs);

			if (ssf->peptide.get_num_aas()>2)
			{
				num_benchmark++;
				int idx;
				for (idx=0; idx<solutions.size() && num_solutions; idx++)
					if (solutions[idx].check_if_correct(ssf->peptide.as_string(config),config))
					{
						correct_benchmark++;
						break;
					}
			}


			scan_count++;
			
		}
	}
	

	if (report_progress)
	{
		time_t curr_time = time(NULL);
		double elapsed_time = (curr_time - last_time);
		cout << "Total running time: " << elapsed_time << endl;
	}

	cout << "#Processed spectra " << total_number_processed << "/" << total_number_headers << endl;


	if (num_benchmark>0)
	{
		cout << endl << "#Correct spectra " << correct_benchmark << "/" << num_benchmark << " (" <<
			setprecision(3) << fixed << (float)correct_benchmark/num_benchmark << ")" << endl;
	}*/
}




// makes a FASTA file with the sequences of full denovo sequences (completed
// from the SEQ in the annotated mgf file)
void make_denovo_training_fa(AllScoreModels& model,
								 char *spectrum_file)
{
/*	const Config *config = model.get_config();
	FileManager fm;
	FileSet fs;


	// Quick read, get all pointers to begining of spectra
	if (get_file_extension_type(spectrum_file) != MZXML)
	{
		fm.init_from_file(config,spectrum_file);
	}
	else // reads peaks 
		fm.init_and_read_single_mzXML(config,spectrum_file,0);

	fs.select_all_files(fm);

	cout << "Generating fa for: " << spectrum_file << endl;
	cout << "Scans: " << fs.get_total_spectra() << endl;

	string file_name,fa_file;
	getFileNameWithoutExtension(spectrum_file,file_name);
	fa_file = file_name + ".dnv.fa";
	ofstream fa_stream(fa_file.c_str());
	
	time_t start_time;
	start_time = time(NULL);

	vector<PrmGraph *> prm_ptrs;
	BasicSpecReader bsr;
	int too_few_peaks=0;
	int bad_charge=0;
	int no_sols=0;
	int scan_count=0;

	prm_ptrs.resize(50,NULL);

	///////////////////////////////////////////////
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
	int sc;
	for (sc=0; sc<all_ssf.size(); sc++)
	{
		static vector<QCPeak> peaks;
		SingleSpectrumFile *ssf = all_ssf[sc];
		if (peaks.size()<ssf->num_peaks)
		{
			int new_size = ssf->num_peaks*2;
			if (new_size<2500)
				new_size=2500;
			peaks.resize(new_size); 
		}

	//	if (ssf->get_scan()<11254)
	//		continue;

	
		const int num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0]);
		ssf->file_idx = 0;

		if (num_peaks<10)
		{
			too_few_peaks++;
			continue;
		}

		Spectrum s;
		s.init_from_QCPeaks(config,&peaks[0],num_peaks,ssf);

		if ( ssf->charge > model.get_max_score_model_charge())
		{
			bad_charge++;
			continue;
		}

		vector<SeqPath> seqpath_solutions;
		vector<int> charges;
		vector<mass_t> pms_with_19;

		pms_with_19.clear();
		charges.clear();		
		BasicSpectrum bs;
		bs.ssf = ssf;
		bs.peaks = &peaks[0];
		bs.num_peaks = num_peaks;

		Peptide correct_peptide = ssf->peptide;

		// output m/z and prob values for the different charge states
		vector<PmcSqsChargeRes> pmcsqs_res;
		//model.select_pms_and_charges(config,bs,pms_with_19,charges,&pmcsqs_res);
	
		cout << "Scan " << bs.ssf->get_scan() << ",\tch: " << charges[0] << setprecision(1) << "\tpm19: " <<
			pms_with_19[0] << setprecision(3) << " \tSQS: " <<  pmcsqs_res[charges[0]].sqs_prob;
		if (pmcsqs_res[charges[0]].sqs_prob<0.01 && ! config->get_use_spectrum_charge())
		{
		//	cout << "  - skipping..." << endl;
		//	continue;
		}

		pms_with_19.clear();
		charges.clear();

		pms_with_19.push_back(s.get_true_mass_with_19());
		charges.push_back(s.getCharge());
	
		generate_denovo_solutions_from_several_pms_with_good_start_end_idxs(prm_ptrs,
				&model,&s,true,150,6,14,pms_with_19,charges,seqpath_solutions);
		
		int j;
		for (j=0; j<seqpath_solutions.size(); j++)
		{
			const int num_correct_aas = seqpath_solutions[j].get_num_correct_aas(correct_peptide,config);
			if (num_correct_aas==seqpath_solutions[j].get_num_aa())
			{
				seqpath_solutions[j]=seqpath_solutions[seqpath_solutions.size()-1];
				seqpath_solutions.pop_back();
			}
		}
		
		vector<PeptideSolution> pep_solutions;
		if (seqpath_solutions.size()>0)
		{
			// create peptide solutions
			pep_solutions.resize(seqpath_solutions.size());
			int s_idx;
			for (s_idx=0; s_idx<seqpath_solutions.size(); s_idx++)
			{
				convert_seq_path_to_peptide_soluition_and_fill_in_aas(config,
					correct_peptide, seqpath_solutions[s_idx],pep_solutions[s_idx]);
			}
		}
		else
			continue;
		

		if (seqpath_solutions.size() == 0)
		{
			cout << "   - no sols" << endl;
			no_sols++;
			continue;
		}
		

		fa_stream << ">" << file_name << ":" << ssf->get_scan() << endl;

		// sample de novo sequences
		int i;
		for (i=9; i<pep_solutions.size(); i+=10)
		{
			int idx = int(myRandom()*i);
			pep_solutions[idx].pep.convert_to_org(config);
			fa_stream << pep_solutions[idx].pep.as_string(config);
			pep_solutions[idx].num_correct_aas=pep_solutions[idx].pep.calc_number_of_correct_aas(config,correct_peptide);
			cout << " " << idx << "(" << pep_solutions[idx].num_correct_aas << ")";
		}
		cout << endl;
		fa_stream << endl;


		scan_count++;

	//	cout << "   " << seqpath_solutions.size() << " sols." << endl;

		if (sc % 200 == 0)
		{
			time_t curr_time = time(NULL);
			cout << setprecision(1) << fixed;
			cout << "done " << sc << "/" << all_ssf.size() << "  " << curr_time - start_time << " secs.,  "
				<< scan_count << " scans completed." << endl;
		}
		
	}
	fa_stream.close();

	cout << "Done..." << endl;
	cout << "fa file: " << fa_file << endl;
	cout << "Wrote fa for " << scan_count << " spectra." << endl;
	if (bad_charge>0)
		cout << "Scans with bad charges  : " << bad_charge << endl;
	if (too_few_peaks>0)
		cout << "Scans with too few peaks: " << too_few_peaks << endl;
	if (no_sols>0)
		cout << "Scans for which no tags were found: " << no_sols << endl;*/
}


void benchmark_tags(AllScoreModels& model,
					char *list,
					char *tag_string,
					int num_test_cases)
{
/*	const Config *config = model.get_config();
	FileManager fm;
	FileSet fs;

	vector<int> num_tags; // how many tags from each length
	int main_length;

	parse_tag_string(tag_string,main_length,num_tags);
	fm.init_from_list_file(config,list);
	fs.select_all_files(fm);

	cout << "benchmarking tags for : " << list << endl;
	cout << "Scans                 : " << fs.get_total_spectra() << endl;
	cout << "Tag types             : " << tag_string << endl;
	if (num_test_cases<0)
	{
		num_test_cases=fs.get_total_spectra();
	}
	else
	{
		fs.randomly_reduce_ssfs(num_test_cases);
	}
	cout << "Will try to get " << num_test_cases << " test cases" << endl;
	
	time_t start_time;
	start_time = time(NULL);

	vector<PrmGraph *> prm_ptrs;
	BasicSpecReader bsr;
	int too_few_peaks=0;
	int bad_charge=0;
	int no_tags=0;
	int scan_count=0;
	vector<int> first_correct_len;
	vector<int> correct_ranks;

	int num_spectra_tested=0;
	int had_correct =0;
	first_correct_len.resize(10,0);

	prm_ptrs.resize(50,NULL);

	///////////////////////////////////////////////
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
	int sc;
	for (sc=0; sc<all_ssf.size(); sc++)
	{
		static vector<QCPeak> peaks;
		SingleSpectrumFile *ssf = all_ssf[sc];
		string peptide_str = ssf->peptide.as_string(config);

		if (peaks.size()<ssf->num_peaks)
		{
			int new_size = ssf->num_peaks*2;
			if (new_size<2500)
				new_size=2500;
			peaks.resize(new_size); 
		}

		vector<int>    tweak_charges;
		vector<mass_t> tweak_pms_with_19;

		tweak_charges.resize(6,0);
		tweak_pms_with_19.resize(6,0);

		const int num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0]);
		ssf->file_idx = 0;

		if (num_peaks<10)
		{
			too_few_peaks++;
			continue;
		}

		Spectrum s;
		s.init_from_QCPeaks(config,&peaks[0],num_peaks,ssf);

	//	s.print_spectrum();

		if ( ssf->charge > model.get_max_score_model_charge())
		{
			bad_charge++;
			continue;
		}

		vector<SeqPath> solutions;
		vector<mass_t> pms_with_19, alt_pms_with_19;
		vector<int>    charges,     alt_charges;

		pms_with_19.clear();
		charges.clear();		
		BasicSpectrum bs;
		bs.ssf = ssf;
		bs.peaks = &peaks[0];
		bs.num_peaks = num_peaks;

	
		
		// output m/z and prob values for the different charge states
		vector<PmcSqsChargeRes> pmcsqs_res;
		//model.select_pms_and_charges(config,bs,pms_with_19,charges,&pmcsqs_res);

		if (pmcsqs_res[charges[0]].sqs_prob<0.01)
		{
		//	continue;
		}

		cout << sc << "\t" << num_peaks << "\t" << setprecision(2) << pmcsqs_res[charges[0]].sqs_prob << "\t" <<
			ssf->peptide.as_string(config) << "  (" << ssf->peptide.get_mass_with_19() - s.get_org_pm_with_19() << ")" << endl;

		
		alt_charges.clear();
		alt_pms_with_19.clear();
		while (charges.size()>0 && charges[0] != charges[charges.size()-1])
		{
			alt_charges.push_back(charges[charges.size()-1]);
			charges.pop_back();
			alt_pms_with_19.push_back(pms_with_19[pms_with_19.size()-1]);
			pms_with_19.pop_back();
		}

		if (pms_with_19.size()>0)
		{
			//generate_tags(prm_ptrs,&model,bs,&s,num_tags, main_length,
			//			  pms_with_19,charges,solutions);

	
		}

		if (alt_pms_with_19.size()>0)
		{
			vector<SeqPath> alt_solutions;

			//generate_tags(prm_ptrs,&model,bs,&s, num_tags, main_length,
			//			  alt_pms_with_19, alt_charges, alt_solutions, false, alt_solutions.size());
			int j;
			for (j=0; j<alt_solutions.size(); j++)
				solutions.push_back(alt_solutions[j]);
		}
		
		num_spectra_tested++;

		if (num_test_cases==num_spectra_tested)
			break;
		
		if (solutions.size() == 0)
		{
			no_tags++;
			continue;
		}
		

		int j;
		for (j=0; j<solutions.size(); j++)
			if (solutions[j].check_if_correct(peptide_str,config))
			{
				had_correct++;
				first_correct_len[solutions[j].get_num_aa()]++;
				correct_ranks.push_back(j);
				break;
			}

		if (sc % 200 == 0)
		{
			time_t curr_time = time(NULL);
			cout << setprecision(1) << fixed;
			cout << "done " << sc << "/" << all_ssf.size() << "  " << curr_time - start_time << " secs.,  "
				<< scan_count << " scans completed." << endl;
		}	
	}

	cout << "Tested " << num_spectra_tested << " spectra" << endl;
	if (bad_charge>0)
		cout << "Scans with bad charges  : " << bad_charge << endl;
	if (too_few_peaks>0)
		cout << "Scans with too few peaks: " << too_few_peaks << endl;
	if (no_tags>0)
		cout << "Scans for which no tags were found: " << no_tags << endl;

	cout << "% with correct " << setprecision(3) << (float)had_correct/(float)num_spectra_tested << endl;
	cout << "First length breakdown: " << endl;
	int j;
	for (j=1; j<first_correct_len.size(); j++)
		if (first_correct_len[j]>0)
			cout << j << "\t" << (float)first_correct_len[j]/(float)num_spectra_tested << endl;
	
	cout << endl;

//	1 3 5 10 25 50 100
	vector<int> counts,seps;

	seps.push_back(1);
	seps.push_back(3);
	seps.push_back(5);
	seps.push_back(10);
	seps.push_back(25);
	seps.push_back(50);
	seps.push_back(100);
	seps.push_back(250);
	seps.push_back(500);

	counts.resize(seps.size(),0);
	for (j=0; j<correct_ranks.size(); j++)
	{
		int k;
		for (k=0; k<seps.size(); k++)
			if (correct_ranks[j]<seps[k])
				counts[k]++;
	}
	
	for (j=0; j<seps.size(); j++)
	{
		cout << seps[j] << "\t" << counts[j]/(float)all_ssf.size() << endl;
	}
	cout << endl;*/
}


void perform_denovo_on_list_of_files(AllScoreModels& model, 
									 const vector<string>& list_vector, 
									 int file_start_idx,
									 int num_solutions, 
									 int min_length, 
									 int max_length,
									 bool report_progress,
									 float min_filter_prob,
									 bool	output_aa_probs,
									 bool	output_cumulative_seq_probs,
									 ostream& out_stream)
{

/*	const int num_rerank_sols_per_charge[] = {0,300,300,300,100,100,100,100,100,100,100,100,100,100};
	
	const Config *config = model.get_config();
	int correct_benchmark=0;
	int total_benchmark=0;
	int spec_counter=0;
	int total_number_headers=0;
	int total_number_processed=0;

	PeptideRankScorer *rank_model = (PeptideRankScorer *)model.get_rank_model_ptr(1);
	static PrmGraph *prm_ptr = NULL;
	static vector<PrmGraph *> prm_ptrs;

	time_t start_time,last_time;
	start_time = time(NULL);
	last_time = start_time;

	int f;
	for (f=0; f<list_vector.size(); f++) 
	{
		const char *spectra_file = list_vector[f].c_str();
		FileManager fm;
		FileSet fs;
		BasicSpecReader bsr;



		///////////////////////////////////////////////
		// Quick read, get all pointers to begining of spectra
		if (get_file_extension_type(list_vector[f]) != MZXML)
		{
			fm.init_from_file(config,spectra_file);
		}
		else // reads peaks 
			fm.init_and_read_single_mzXML(config,spectra_file,f);

		fs.select_all_files(fm);
	//	fs.select_files_in_mz_range(fm,555,645);

		const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

		total_number_headers+= all_ssf.size();

		int sc;
		for (sc=1; sc<all_ssf.size(); sc++)
		{
			static vector<QCPeak> peaks;
			SingleSpectrumFile *ssf = all_ssf[sc];
			if (peaks.size()<ssf->num_peaks)
			{
				int new_size = ssf->num_peaks*2;
				if (new_size<2500)
					new_size=2500;
				peaks.resize(new_size); 
			}

		//	if (sc>241)
		//		break;

		//	if (ssf->get_scan() != 241)
		//		continue;

			time_t curr_time = time(NULL);
			double elapsed_time = (curr_time - last_time);
			if (report_progress && elapsed_time>5)
			{
				last_time = curr_time;
				cout << "Processing file " << f+1 << "/" << list_vector.size();
				if (all_ssf.size() == 1)
				{
					cout << "  spectrum 1/1 of current file." << endl;
				}
				else
					cout << "  spectrum " << sc+1 << "/" << all_ssf.size() << 
					" of current file." << endl;
			}

			const int num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0]);
			ssf->file_idx = f+file_start_idx;
			

			// convert peak list ot a spectrum with charge (if original charge ==0)
			// the spectrum gets charge 2, but the true charge is computed from the data
			
			if (num_peaks<5)
			{
				ssf->print_ssf_stats(config,out_stream);
				out_stream << "# too few peaks..." << endl;
				continue;
			}

			spec_counter++;
			Spectrum s;
			s.init_from_QCPeaks(config,&peaks[0],num_peaks,ssf);

			s.print_expected_by();
			s.print_spectrum();
			cout << endl;
		//	exit(0);

			vector<SeqPath> solutions;
			solutions.clear();

			if ( ssf->charge > model.get_max_score_model_charge())
			{
				ssf->print_ssf_stats(config,out_stream);
				out_stream << "# Charge " << s.getCharge() << " not supported yet..." << endl << endl;
				continue;
			}

			
			bool perform_rerank=false;
			int rerank_size_idx = NEG_INF;
			int num_sols_to_generate_before_ranking=num_solutions;
			float spectrum_quality = 0;
		
			if (1)
			{
				vector<mass_t> pms_with_19;
				vector<int>    charges;
				pms_with_19.clear();
				charges.clear();		
				BasicSpectrum bs;
				bs.ssf = ssf;
				bs.peaks = &peaks[0];
				bs.num_peaks = num_peaks;

				// output m/z and prob values for the different charge states
				vector<PmcSqsChargeRes> pmcsqs_res;
				//model.select_pms_and_charges(config,bs,pms_with_19,charges,&pmcsqs_res);
				if (pmcsqs_res.size()>charges[0] && pmcsqs_res[charges[0]].sqs_prob<min_filter_prob)
				{
					ssf->print_ssf_stats(config,out_stream,false);
					out_stream << " (SQS " << spectrum_quality << ")" << endl;
					out_stream << "# low quality, skipping: " << pmcsqs_res[charges[0]].sqs_prob << endl << endl;
					continue;
				}
				spectrum_quality = pmcsqs_res[charges[0]].sqs_prob;

				cout << "QUAL: " << spectrum_quality << endl;
				int j;
				for (j=0; j<pms_with_19.size(); j++)
					cout << j << "\t" << pms_with_19[j] << "\t" << charges[j] << endl;
			//	exit(0);


				if (prm_ptrs.size()<pms_with_19.size())
					prm_ptrs.resize(pms_with_19.size(),NULL);

				if (pms_with_19[0]<100)
				{
					ssf->print_ssf_stats(config,out_stream);
					out_stream << "# Could not process spectrum..." << endl << endl;
					continue;
				}
				
				const int num_rerank_per_charge =  num_rerank_sols_per_charge[charges[0]];
				if (rank_model && num_sols_to_generate_before_ranking<num_rerank_per_charge)
						num_sols_to_generate_before_ranking=num_rerank_per_charge;
				
				generate_denovo_solutions_from_several_pms(
					prm_ptrs,
					&model,
					&s,
					true, 
					num_sols_to_generate_before_ranking,
					min_length,
					max_length,
					pms_with_19,
					charges,
					solutions,
					false);

				// use charge of top scoring solution
				if (solutions.size()>1)
				{
					const int sol_charge = solutions[0].charge;
					int j;
					for (j=1; j<solutions.size(); j++)
					{
						if (solutions[j].charge != sol_charge)
						{
							if (j==solutions.size()-1)
							{
								solutions.pop_back();
							}
							else
							{
								solutions[j]=solutions[solutions.size()-1];
								solutions.pop_back();
								j--;
							}
						}
					}
				}
			}

			if (rank_model && solutions.size()>0)
			{
				rerank_size_idx = config->calc_size_idx(solutions[0].charge,solutions[0].pm_with_19);
				if (rank_model->get_ind_part_model_was_initialized(solutions[0].charge,rerank_size_idx))
				{
					perform_rerank=true;
				}
			}

			
			vector<score_pair> score_idx_pairs;
			if (perform_rerank)
			{
			//	rank_model->score_denovo_sequences(solutions,ssf,&peaks[0],num_peaks,
			//		score_idx_pairs,rerank_size_idx);
				int i;
				for (i=0; i<score_idx_pairs.size(); i++)
					solutions[score_idx_pairs[i].idx].rerank_score = score_idx_pairs[i].score;
				sort(score_idx_pairs.begin(),score_idx_pairs.end());
			}
			else
			{
				score_idx_pairs.resize(solutions.size());
				int i;
				for (i=0; i<solutions.size(); i++)
					score_idx_pairs[i].idx=i;
			}

			////////////////////////////////////////////////////////////
			// if we are here it is only for denovo/tags
			// print results
			////////////////////////////////////////////////////////////

			
			bool had_pep = false;
			bool had_correct = false;

			total_number_processed++;

			ssf->print_ssf_stats(config,out_stream,false);
			out_stream << " (SQS " << spectrum_quality << ")" << endl;

			if (solutions.size() == 0)
			{
				out_stream << "# No solutions found." << endl;
			}
			else 
			{

				out_stream << "#Index\t";
				out_stream << "RnkScr\t";
				if (output_cumulative_seq_probs)
					out_stream << "CumProb\t";

				out_stream << "PnvScr\tN-Gap\tC-Gap\t[M+H]\tCharge\tSequence" << endl;

				if (output_aa_probs || output_cumulative_seq_probs)
				{
					for (int i=0; i<solutions.size() && i<num_solutions; i++)
					{
						const int idx = (perform_rerank ? score_idx_pairs[i].idx : i);
						const vector<PathPos>& positions = solutions[idx].positions;
						solutions[idx].prm_ptr->calc_amino_acid_probs(solutions[idx],i);
					}
				}

				if (output_cumulative_seq_probs)
				{
					const int first_sol_charge = solutions[0].charge;
					const int first_sol_size_idx = config->calc_size_idx(first_sol_charge,solutions[0].pm_with_19);
					CumulativeSeqProbModel* csp_model = (CumulativeSeqProbModel* )model.get_cumulative_seq_prob_model_ptr(0);
			
					csp_model->calc_cumulative_seq_probs(first_sol_charge, first_sol_size_idx, 
						spectrum_quality, score_idx_pairs, solutions); 
				}

				int i; 	
				for (i=0; i<solutions.size() && i<num_solutions; i++) 
				{
					const int idx = (perform_rerank ? score_idx_pairs[i].idx : i);

					mass_t c_gap=solutions[idx].pm_with_19 - solutions[idx].c_term_mass;
					if (c_gap<24.0)
						c_gap = 0;

					out_stream << setprecision(3) << fixed << i << "\t";
				
					if (perform_rerank)
					{
						out_stream << score_idx_pairs[i].score << "\t";
					}
					else
						out_stream << -999 << "\t";

					if (output_cumulative_seq_probs)
						out_stream << solutions[idx].cumulative_seq_prob << "\t";

					out_stream << solutions[idx].path_score << "\t";
					out_stream << solutions[idx].n_term_mass << "\t";
					out_stream << c_gap << "\t";
					out_stream << solutions[idx].pm_with_19 << "\t";
					out_stream << solutions[idx].charge << "\t";
					out_stream << solutions[idx].seq_str;	
		
					if (output_aa_probs)
					{
						const vector<PathPos>& positions = solutions[idx].positions;
						solutions[idx].prm_ptr->calc_amino_acid_probs(solutions[idx],i);
						int j;
						for (j=0; j<positions.size()-1; j++)
							out_stream << "\t" << positions[j].edge_variant_prob;
					}

					if (ssf->peptide.get_num_aas()>2)
					{
						if (solutions[idx].check_if_correct(ssf->peptide.as_string(config),config))
						{

							out_stream << " *";

							if (! had_correct)
							{
								correct_benchmark++;
								had_correct=true;
							}
						}
						had_pep=true;
					}
					out_stream << endl;
				}
			}

			if (had_pep) // for annotated spectra (benchmark)
				total_benchmark++;

			out_stream << endl;

			exit(0);
		}
	}

	/////////////////////////////////////////////////////////////////
	// this part works only if the spectra are annotated (benchmark)
	/////////////////////////////////////////////////////////////////

	cout << "#Processed spectra " << total_number_processed << "/" << total_number_headers << endl;

	if (total_benchmark>0)
	{
		cout << "#Correct spectra " << correct_benchmark << "/" << total_benchmark << " (" <<
			fixed << setprecision(3) << (double)correct_benchmark/(double)total_benchmark << ")" << endl;
	}

	
	if (report_progress)
	{
		time_t curr_time = time(NULL);
		double elapsed_time = (last_time - start_time);
		cout << "Total running time: " << elapsed_time << endl;
		cout << "Processed " << list_vector.size() << " files and " << spec_counter << " spectra." << endl;
	}
	*/

}



void new_perform_denovo_on_list_of_files(AllScoreModels& model, 
									 const vector<string>& list_vector, 
									 int file_start_idx,
									 int num_solutions, 
									 int min_length, 
									 int max_length,
									 bool report_progress,
									 float min_filter_prob,
									 bool	output_aa_probs,
									 bool	output_cumulative_seq_probs,
									 ostream& out_stream)
{

	const int num_rerank_sols_per_charge[] = {0,300,300,300,100,100,100,100,100,100,100,100,100,100};
	
	const Config *config = model.get_config();
	int correct_benchmark=0;
	int total_benchmark=0;
	int spec_counter=0;
	int total_number_headers=0;
	int total_number_processed=0;

	PeptideRankScorer *rank_model = model.get_rank_model_ptr(1);
	PrmGraph *prm_ptr = NULL;
	vector<PrmGraph *> prm_ptrs;

	time_t start_time,last_time;
	start_time = time(NULL);
	last_time = start_time;

	int f;
	for (f=0; f<list_vector.size(); f++) 
	{
		const char *spectra_file = list_vector[f].c_str();
		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(spectra_file, config);

		SpectraList sl(sa);
		sl.selectAllAggregatorHeaders();


		total_number_headers+= sl.getNumHeaders();

		int sc;
		for (sc=0; sc<sl.getNumHeaders(); sc++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(sc);

			time_t curr_time = time(NULL);
			double elapsed_time = (curr_time - last_time);
			if (report_progress && elapsed_time>5)
			{
				last_time = curr_time;
				cout << "Processing file " << f+1 << "/" << list_vector.size();
				if (sl.getNumHeaders() == 1)
				{
					cout << "  spectrum 1/1 of current file." << endl;
				}
				else
					cout << "  spectrum " << sc+1 << "/" << sl.getNumHeaders() << 
					" of current file." << endl;
			}

			AnnotatedSpectrum as;
			as.readSpectrum(sa, header);

			header->printStats(config, out_stream);
			
		//	as.print_expected_by();
		//	as.print_spectrum();
		//	cout << endl;
			
			
			if (as.getNumPeaks()==0)
			{
				header->printStats(config, out_stream);
			//	ssf->print_ssf_stats(config,out_stream);
				out_stream << "# too few peaks..." << endl;
				continue;
			}

			spec_counter++;
			
			vector<SeqPath> solutions;
			solutions.clear();

			if ( as.getCharge() > model.get_max_score_model_charge())
			{
				header->printStats(config, out_stream);
			//	ssf->print_ssf_stats(config,out_stream);
				out_stream << "# Charge " << as.getCharge() << " not supported yet..." << endl << endl;
				continue;
			}

			
			bool perform_rerank=false;
			int rerank_size_idx = NEG_INF;
			int num_sols_to_generate_before_ranking=num_solutions;
			float spectrum_quality = 0;
		
			if (1)
			{
				vector<mass_t> pms_with_19;
				vector<int>    charges;
				pms_with_19.clear();
				charges.clear();		
			

				// output m/z and prob values for the different charge states
				vector<PmcSqsChargeRes> pmcsqs_res;
				// model.select_pms_and_charges(config,bs,pms_with_19,charges,&pmcsqs_res);
				model.selectPrecursorMassesAndCharges(config, as, pms_with_19, charges, &pmcsqs_res);
				
				if (pmcsqs_res.size()>charges[0] && pmcsqs_res[charges[0]].sqs_prob<min_filter_prob)
				{
					header->printStats(config, out_stream);
				//	ssf->print_ssf_stats(config,out_stream,false);
					out_stream << " (SQS " << spectrum_quality << ")" << endl;
					out_stream << "# low quality, skipping: " << pmcsqs_res[charges[0]].sqs_prob << endl << endl;
					continue;
				}
				spectrum_quality = pmcsqs_res[charges[0]].sqs_prob;

			/*	cout << "QUAL: " << spectrum_quality << endl;
				int j;
				for (j=0; j<pms_with_19.size(); j++)
					cout << j << "\t" << pms_with_19[j] << "\t" << charges[j] << endl;
				exit(0); */


				if (prm_ptrs.size()<pms_with_19.size())
					prm_ptrs.resize(pms_with_19.size(),NULL);

				if (pms_with_19[0]<100)
				{
					header->printStats(config, out_stream);
				//	ssf->print_ssf_stats(config,out_stream);
					out_stream << "# Could not process spectrum..." << endl << endl;
					continue;
				}
				
				const int num_rerank_per_charge =  num_rerank_sols_per_charge[charges[0]];
				if (rank_model && num_sols_to_generate_before_ranking<num_rerank_per_charge)
						num_sols_to_generate_before_ranking=num_rerank_per_charge;
				
				generate_denovo_solutions_from_several_pms(
					prm_ptrs,
					&model,
					&as,
					true, 
					num_sols_to_generate_before_ranking,
					min_length,
					max_length,
					pms_with_19,
					charges,
					solutions,
					false);

				// use charge of top scoring solution
				if (solutions.size()>1)
				{
					const int sol_charge = solutions[0].charge;
					int j;
					for (j=1; j<solutions.size(); j++)
					{
						if (solutions[j].charge != sol_charge)
						{
							if (j==solutions.size()-1)
							{
								solutions.pop_back();
							}
							else
							{
								solutions[j]=solutions[solutions.size()-1];
								solutions.pop_back();
								j--;
							}
						}
					}
				}
			}

			if (rank_model && solutions.size()>0)
			{
				rerank_size_idx = config->calc_size_idx(solutions[0].charge,solutions[0].pm_with_19);
				if (rank_model->get_ind_part_model_was_initialized(solutions[0].charge,rerank_size_idx))
				{
					perform_rerank=true;
				}
			}

			
			vector<score_pair> score_idx_pairs;
			if (perform_rerank)
			{
			//	rank_model->score_denovo_sequences(solutions,ssf,&peaks[0],num_peaks,
			//		score_idx_pairs,rerank_size_idx);
				rank_model->scoreDenovoSequences(solutions, as, score_idx_pairs, rerank_size_idx);
				int i;
				for (i=0; i<score_idx_pairs.size(); i++)
					solutions[score_idx_pairs[i].idx].rerank_score = score_idx_pairs[i].score;
				sort(score_idx_pairs.begin(),score_idx_pairs.end());
			}
			else
			{
				score_idx_pairs.resize(solutions.size());
				int i;
				for (i=0; i<solutions.size(); i++)
					score_idx_pairs[i].idx=i;
			}

			////////////////////////////////////////////////////////////
			// if we are here it is only for denovo/tags
			// print results
			////////////////////////////////////////////////////////////

			
			bool had_pep = false;
			bool had_correct = false;

			total_number_processed++;

			//ssf->print_ssf_stats(config,out_stream,false);
			header->printStats(config, out_stream);
			out_stream << " (SQS " << spectrum_quality << ")" << endl;

			if (solutions.size() == 0)
			{
				out_stream << "# No solutions found." << endl;
			}
			else 
			{

				out_stream << "#Index\t";
				out_stream << "RnkScr\t";
				if (output_cumulative_seq_probs)
					out_stream << "CumProb\t";

				out_stream << "PnvScr\tN-Gap\tC-Gap\t[M+H]\tCharge\tSequence" << endl;

				if (output_aa_probs || output_cumulative_seq_probs)
				{
					for (int i=0; i<solutions.size() && i<num_solutions; i++)
					{
						const int idx = (perform_rerank ? score_idx_pairs[i].idx : i);
						const vector<PathPos>& positions = solutions[idx].positions;
						solutions[idx].prm_ptr->calc_amino_acid_probs(solutions[idx],i);
					}
				}

				if (output_cumulative_seq_probs)
				{
					const int first_sol_charge = solutions[0].charge;
					const int first_sol_size_idx = config->calc_size_idx(first_sol_charge,solutions[0].pm_with_19);
					CumulativeSeqProbModel* csp_model = (CumulativeSeqProbModel* )model.get_cumulative_seq_prob_model_ptr(0);
			
					csp_model->calc_cumulative_seq_probs(first_sol_charge, first_sol_size_idx, 
						spectrum_quality, score_idx_pairs, solutions); 
				}

				int i; 	
				for (i=0; i<solutions.size() && i<num_solutions; i++) 
				{
					const int idx = (perform_rerank ? score_idx_pairs[i].idx : i);

					mass_t c_gap=solutions[idx].pm_with_19 - solutions[idx].c_term_mass;
					if (c_gap<24.0)
						c_gap = 0;

					out_stream << setprecision(3) << fixed << i << "\t";
				
					if (perform_rerank)
					{
						out_stream << score_idx_pairs[i].score << "\t";
					}
					else
						out_stream << -999 << "\t";

					if (output_cumulative_seq_probs)
						out_stream << solutions[idx].cumulative_seq_prob << "\t";

					out_stream << solutions[idx].path_score << "\t";
					out_stream << solutions[idx].n_term_mass << "\t";
					out_stream << c_gap << "\t";
					out_stream << solutions[idx].pm_with_19 << "\t";
					out_stream << solutions[idx].charge << "\t";
					out_stream << solutions[idx].seq_str;	
		
					if (output_aa_probs)
					{
						const vector<PathPos>& positions = solutions[idx].positions;
						solutions[idx].prm_ptr->calc_amino_acid_probs(solutions[idx],i);
						int j;
						for (j=0; j<positions.size()-1; j++)
							out_stream << "\t" << positions[j].edge_variant_prob;
					}

					if (header->getPeptideStr().length()>2 )
					{
						if (solutions[idx].check_if_correct(header->getPeptideStr(),config))
						{

							out_stream << " *";

							if (! had_correct)
							{
								correct_benchmark++;
								had_correct=true;
							}
						}
						had_pep=true;
					}
					out_stream << endl;
				}
			}

			if (had_pep) // for annotated spectra (benchmark)
				total_benchmark++;

			out_stream << endl;

		//	exit(0);
		}
	}

	/////////////////////////////////////////////////////////////////
	// this part works only if the spectra are annotated (benchmark)
	/////////////////////////////////////////////////////////////////

	cout << "#Processed spectra " << total_number_processed << "/" << total_number_headers << endl;

	if (total_benchmark>0)
	{
		cout << "#Correct spectra " << correct_benchmark << "/" << total_benchmark << " (" <<
			fixed << setprecision(3) << (double)correct_benchmark/(double)total_benchmark << ")" << endl;
	}

	
	if (report_progress)
	{
		time_t curr_time = time(NULL);
		double elapsed_time = (last_time - start_time);
		cout << "Total running time: " << elapsed_time << endl;
		cout << "Processed " << list_vector.size() << " files and " << spec_counter << " spectra." << endl;
	}

}




void perform_prm_on_list_of_files(AllScoreModels& model, 
								  const vector<string>& list_vector,
								  float sqs_filter_prob,
								  int file_start_idx,
								  bool prm_norm)
{
/*	const Config *config = model.get_config();

	int f;
	for (f=0; f<list_vector.size(); f++) 
	{
		const char *spectra_file = list_vector[f].c_str();
		FileManager fm;
		FileSet fs;
		BasicSpecReader bsr;

		///////////////////////////////////////////////
		// Quick read, get all pointers to begining of spectra
		if (get_file_extension_type(list_vector[f]) != MZXML)
		{
			fm.init_from_file(config,spectra_file);
		}
		else // reads peaks 
			fm.init_and_read_single_mzXML(config,spectra_file,f);

		fs.select_all_files(fm);

		const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
		int sc;
		for (sc=0; sc<all_ssf.size(); sc++)
		{
			static vector<QCPeak> peaks;
			SingleSpectrumFile *ssf = all_ssf[sc];
			if (peaks.size()<ssf->num_peaks)
			{
				int new_size = ssf->num_peaks*2;
				if (new_size<2500)
					new_size=2500;
				peaks.resize(new_size); 
			}

			const int num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0]);
			ssf->file_idx = f+file_start_idx;

		//	if (ssf->get_scan() != 145)
		//		continue;


			// convert peak list ot a spectrum with charge (if original charge ==0)
			// the spectrum gets charge 2, but the true charge is computed from the data
		
			Spectrum s;
			s.init_from_QCPeaks(config,&peaks[0],num_peaks,ssf);

			ssf->print_ssf_stats(config);
			if ( ssf->charge > model.get_max_score_model_charge())
			{
				cout << "# Charge " << s.getCharge() << " not supported yet..." << endl << endl;
				continue;
			}


			if (num_peaks<5)
			{
				cout << "# too few peaks..." << endl;
				continue;
			}

			vector<mass_t> pms_with_19;
			vector<int>    charges;

			pms_with_19.clear();
			charges.clear();		
			
			BasicSpectrum bs;
			bs.ssf = ssf;
			bs.peaks = &peaks[0];
			bs.num_peaks = num_peaks;

			// output m/z and prob values for the different charge states
			//model.select_pms_and_charges(config,bs,pms_with_19,charges);
			if (pms_with_19.size()<=0)
			{
				cout << "# Couldn't find PM or charge..." << endl;
				continue;
			}

			cout << "Charge= " << charges[0] << "  [M+H] = " << pms_with_19[0] << endl;
			print_prm_graph_scores(&model,&s,pms_with_19[0],charges[0], prm_norm);
		}
	}*/
}








/************************************************************************
Creates a new peptide string with correct aas that fill the N/C gaps
*************************************************************************/
bool fill_with_true_peptide_aas(Config *config,
								const Peptide& true_peptide, 
								const SeqPath& solution, 
								string& new_pep_str, 
								int& num_correct, 
								int& num_aas)
{
	const vector<string>& aa2label = config->get_aa2label();
	const vector<int>& true_amino_acids = true_peptide.get_amino_acids();
	vector<mass_t> break_masses;

	true_peptide.calc_expected_breakage_masses(config,break_masses);
	const mass_t total_pred_mass = break_masses[break_masses.size()-1];

	if (solution.n_term_mass<1.0 && solution.c_term_mass > total_pred_mass-5.0)
	{
		new_pep_str = solution.seq_str;
		num_correct = solution.get_num_correct_aas(true_peptide,config);
		num_aas = solution.get_num_aa();
		return true;
	}

	string true_string = true_peptide.as_string(config);
	
	new_pep_str = "";
	num_correct = 0;
	num_aas		= 0;

	// fill in N-term side
	if (solution.n_term_mass>1.0)
	{
		int i;
		for (i=0; i<break_masses.size(); i++)
			if (fabs(solution.n_term_mass-break_masses[i])<1.0)
				break;
		
		if (i==break_masses.size())
			return false;

		int j;
		for (j=0; j<i; j++)
			new_pep_str += aa2label[true_amino_acids[j]];

		num_aas=j;
		num_correct = j;

	//	cout << "I: \t" << new_pep_str << "\t" << new_pep_str.length() << endl ;
	}

	// fill in middle
	new_pep_str += solution.seq_str;
	num_correct += solution.get_num_correct_aas(true_peptide,config);
	num_aas += solution.get_num_aa();

//	cout << "II: \t" << new_pep_str << "\t" << new_pep_str.length() << endl ;

	// fill end
	if (fabs(solution.c_term_mass-total_pred_mass)>5.0)
	{
		int i;
		for (i=0; i<break_masses.size(); i++)
			if (fabs(solution.c_term_mass-break_masses[i])<1.0)
				break;
		
		if (i==break_masses.size())
			return false;

		int j;
		for (j=i; j<break_masses.size()-1; j++)
			new_pep_str += aa2label[true_amino_acids[j]];

		num_aas+=(j-i);
		num_correct += (j-i);
	}

//	cout << "III: \t" << new_pep_str << "\t" << new_pep_str.length() << endl ;

	
	// check the mass of the peptide
	Peptide new_pep;
	new_pep.parseFromString(config,new_pep_str);
	new_pep.calc_mass(config);

	if (fabs(new_pep.get_mass()-total_pred_mass)>5.0)
	{
		cout << "Error constructing training peptide!!!" << endl;
		solution.print_full(config);
		cout << endl << "Created: " << new_pep_str << " "<< new_pep.get_mass() << endl;
		cout << "True: " << true_peptide.as_string(config) << " " << total_pred_mass << endl;
		exit(1);
	}


	return true;
}




bool complete_denovo_seq(Config *config,
							const Peptide& corr_pep,
							const vector<mass_t>& corr_break_masses,
							const SeqPath& path,
							vector<int>& full_aas,
							bool	ind_random_fill)
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const vector<int>&    corr_amino_acids = corr_pep.get_amino_acids();
	int i;
	int first_break=-1;
	int last_break=-1;

	for (i=0; i<corr_break_masses.size(); i++)
		if (fabs(corr_break_masses[i]-path.n_term_mass)<2)
			break;
	if (i==corr_break_masses.size())
		return false;
	first_break = i;

	for ( ; i<corr_break_masses.size(); i++)
		if (fabs(corr_break_masses[i]-path.c_term_mass)<4)
			break;
	if (i==corr_break_masses.size())
		return false;
	last_break=i;

	full_aas.clear();

	if (! ind_random_fill)
	{
		for (i=0; i<first_break; i++)
			full_aas.push_back(corr_amino_acids[i]);

		for (i=0; i<path.positions.size()-1; i++)
			full_aas.push_back(path.positions[i].aa);

		for (i=last_break; i<corr_amino_acids.size(); i++)
			full_aas.push_back(corr_amino_acids[i]);
	}
	else // random fill - permute the correct amino acids
	{
		const vector<string>& aa2labels = config->get_aa2label();
		if (first_break>0)
		{
			full_aas.push_back(corr_amino_acids[0]);
			vector<int> pre_aas;
			pre_aas.clear();
			for (i=1; i<first_break; i++)
				pre_aas.push_back(corr_amino_acids[i]);

			permute_vector(pre_aas);
			for (i=0; i<pre_aas.size(); i++)
				full_aas.push_back(pre_aas[i]);
		}
		
		for (i=0; i<path.positions.size(); i++)
			if (path.positions[i].aa>=0)
				full_aas.push_back(path.positions[i].aa);

		if (last_break<corr_amino_acids.size())
		{
			vector<int> post_aas;
			post_aas.clear();
			for (i=last_break; i<corr_amino_acids.size(); i++)
				post_aas.push_back(corr_amino_acids[i]);

			permute_vector(post_aas);
			for (i=0; i<post_aas.size(); i++)
				full_aas.push_back(post_aas[i]);
		}
	}

	// test peptide 
	mass_t new_mass=MASS_OHHH;
	for (i=0; i<full_aas.size(); i++)
		new_mass+=aa2mass[full_aas[i]];

	const mass_t diff = fabs(new_mass-corr_pep.get_mass_with_19());
	if (diff>15)
	{
		cout << "Error: mismatch in masses when trying to fill training de novo str!" << endl;
		cout << corr_pep.as_string(config) << endl;
		cout << new_mass << " vs. " << corr_pep.get_mass_with_19() << endl;
		exit(1);
	}

	if (diff>6)
		return false;

	return true;
}











