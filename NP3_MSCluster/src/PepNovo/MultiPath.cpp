#include "PrmGraph.h"

struct pos_score_pair {
	pos_score_pair() : pos(NEG_INF), score(NEG_INF) {};
	bool operator< (const pos_score_pair& other) const
	{
		return (score > other.score);
	}
	int pos;
	score_t score;
};


// This struct uses a simple score based comparison
struct PathHeap {
public:
	void init(int size) 
	{ 
		heap_size = size; 
		paths.resize(heap_size); 
		score_pairs.resize(heap_size);
		int i;
		for (i=0; i<heap_size; i++)
		{
			score_pairs[i].pos = i;
			score_pairs[i].score = NEG_INF;
		}
	}

	score_t get_min_score() const { return score_pairs[0].score; }

	int get_num_real_entries() const
	{
		int i,n=0;
		for (i=0; i<score_pairs.size(); i++)
			if (score_pairs[i].score>NEG_INF)
				n++;
		return n;
	}

	void add_path(const SeqPath& path)
	{
		if (path.path_score<=score_pairs[0].score)
			return;
		
		const int pos = score_pairs[0].pos;
		paths[pos]=path;

		pop_heap(score_pairs.begin(),score_pairs.end());
		score_pairs[heap_size-1].pos = pos;
		score_pairs[heap_size-1].score = path.path_score;
		push_heap(score_pairs.begin(),score_pairs.end());
	}

	void sort_paths() { sort(paths.begin(),paths.end(),comp_SeqPath_path_score); }

	vector<SeqPath> get_paths() { return paths; }

private:
	int heap_size;

	vector<pos_score_pair> score_pairs;
	vector<SeqPath> paths;  // this is a min heap smallest score in front!
};






// add the variant ptr and scores for this combo
void PrmGraph::add_and_score_edge_variants(const AA_combo& aa_combo, MultiEdge& edge)
{
	const vector<int>& aa_positions = config->get_aa_positions();
	const bool reaches_n_term = (nodes[edge.n_idx].type == NODE_N_TERM);
	const bool reaches_c_term = (nodes[edge.c_idx].type == NODE_C_TERM);
	int v;
	
	int *variant_ptr = (int *)config->get_variant_ptr(aa_combo.variant_start_idx);

	for (v=0; v<aa_combo.num_variants; v++)
	{
		int num_aa = *variant_ptr;
		int *aas = variant_ptr+1;

		int i;
		for (i=0; i<num_aa; i++)
			if (aa_positions[*(aas+i)])
				break;

		// need to check that the variant is not violating any of the position restrictions
		// such as +1 -1 positions
		if (i<num_aa)
		{
			if (reaches_n_term)
			{
				int a;
				for (a=0; a<num_aa; a++)
				{
					int aa_idx = aas[a];
					if (aa_positions[aa_idx] != 0 && aa_positions[aa_idx] != a+1)
						break; 
				}
				if (a<num_aa)
				{
					variant_ptr+= num_aa +1;
					continue;
				}
			}
			else // check for +1 positions
			{
				int a;
				for (a=0; a<num_aa; a++)
				{
					int aa_idx = aas[a];
					if (aa_positions[aa_idx] == 1)
						break; 
				}
				if (a<num_aa)
				{
					variant_ptr+= num_aa +1;
					continue;
				}
			}

			if (reaches_c_term)
			{
				int a;
				for (a=0; a<num_aa; a++)
				{
					int aa_idx = aas[a];
					if (aa_positions[aa_idx] != 0 && aa_positions[aa_idx] != a-num_aa )
						break; 
				}
				if (a<num_aa) // found a problem with one of the aa positions
				{
					variant_ptr+= num_aa +1;
					continue;
				}
			}
			else // check for -1 positions
			{
				int a;
				for (a=0; a<num_aa; a++)
				{
					int aa_idx = aas[a];
					if (aa_positions[aa_idx] == -1)
						break; 
				}
				if (a<num_aa)
				{
					variant_ptr+= num_aa +1;
					continue;
				}
			}

			

		}
		
		score_t variant_score = calc_edge_variant_score(edge,num_aa,aas);

		if (variant_score>edge.max_variant_score)
			edge.max_variant_score = variant_score;

		edge.variant_ptrs.push_back(variant_ptr);
		edge.variant_scores.push_back(variant_score);

		variant_ptr+= num_aa +1;
	}
}


// adds the relevant PathPos to the path and adjusts the other non-terminal values 
void SeqPath::add_edge_variant(const MultiEdge& edge, int e_idx, int variant_idx)
{
	PathPos new_pos;

	int *variant_ptr = edge.variant_ptrs[variant_idx];
	int num_aa = *variant_ptr++;
	int *aas = variant_ptr;

	if (num_aa != edge.num_aa)
	{
		cout << "Error: edge and variant mixup!" << endl;
		exit(1);
	}
	
	num_aa += edge.num_aa;
	
	new_pos.breakage = edge.n_break;
	new_pos.edge_idx = e_idx;
	new_pos.mass = edge.n_break->mass;
	new_pos.edge_variant_score = edge.variant_scores[variant_idx];
	new_pos.node_score = edge.n_break->score;
	new_pos.node_idx = edge.n_idx;
	new_pos.aa = aas[0];

	path_score += new_pos.edge_variant_score + new_pos.node_score;

	positions.push_back(new_pos);

	if (edge.num_aa == 1)
		return;

	int i;
	for (i=1; i<edge.num_aa; i++)
	{
		PathPos new_pos;

		new_pos.breakage = NULL;
		new_pos.aa = aas[i]; // the rest of the fields are initialized to the default (NULL) values
		positions.push_back(new_pos);
	}
	
}


/************************************************************************

  Expands the single multi path into all possible sequence variants.
  This is done incrementaly by expanding each multi edges

*************************************************************************/
void PrmGraph::expand_multi_path(const MultiPath& multi_path, 
								 vector<SeqPath>& seq_paths,
								 score_t min_score,
								 score_t forbidden_pair_penalty,
								 int max_num_paths) const
{
	const vector<AA_combo>& aa_edge_comobos = config->get_aa_edge_combos();
	vector<score_t> max_attainable_scores;

	const int num_edges = multi_path.edge_idxs.size();
	max_attainable_scores.resize(num_edges+1,0);

	int e;
	for (e=num_edges-1; e>=0; e--)
	{
		const int e_idx = multi_path.edge_idxs[e];
		const MultiEdge& edge = multi_edges[e_idx];

		max_attainable_scores[e]= max_attainable_scores[e+1] + 
			edge.max_variant_score + edge.c_break->score; 
	}

	if (min_score>multi_edges[multi_path.edge_idxs[0]].n_break->score + max_attainable_scores[0])
		return;

	int i;
	seq_paths.resize(1);

	seq_paths[0].n_term_aa = multi_path.n_term_aa;
	seq_paths[0].c_term_aa = multi_path.c_term_aa;
	seq_paths[0].n_term_mass = multi_path.n_term_mass;
	seq_paths[0].c_term_mass = multi_path.c_term_mass;
	seq_paths[0].multi_path_rank = multi_path.original_rank;
	seq_paths[0].path_score = 0;
	seq_paths[0].prm_ptr = (PrmGraph *)this;
	seq_paths[0].positions.clear();

	for (e=0; e<multi_path.edge_idxs.size(); e++)
	{
		const int e_idx = multi_path.edge_idxs[e];
		const MultiEdge& edge = multi_edges[e_idx];

		if (edge.get_num_variants() == 1)
		{
			const int * variant_ptr = edge.variant_ptrs[0];
			const int num_aa = *variant_ptr++;
			const int *aas = variant_ptr;
		
			int i;
			for (i=0; i<seq_paths.size(); i++)
			{
				if (seq_paths[i].path_score>NEG_INF && 
					seq_paths[i].path_score + max_attainable_scores[e] > min_score)
				{
					seq_paths[i].add_edge_variant(edge,e_idx,0);
				}
				else
					seq_paths[i].path_score = NEG_INF;
			}
		}
		else
		{
			vector<SeqPath> old_paths = seq_paths;
			seq_paths.resize(seq_paths.size()*edge.get_num_variants());

			int v;
			int idx=0;
			for (v=0; v<edge.get_num_variants(); v++)
			{
				int i;
				for (i=0; i<old_paths.size(); i++)
					seq_paths[idx++]=old_paths[i];
			}

			
			idx=0;
			for (v=0; v<edge.get_num_variants(); v++)
			{
				const int * variant_ptr =  edge.variant_ptrs[v];
				const int num_aa = *variant_ptr++;
				const int *aas = variant_ptr;

				int i;
				for (i=0; i<old_paths.size(); i++)
				{
					if (seq_paths[idx].path_score>NEG_INF &&
						seq_paths[idx].path_score + max_attainable_scores[e] > min_score)
					{
						seq_paths[idx].add_edge_variant(edge,e_idx,v);
					}
					else
						seq_paths[idx].path_score=NEG_INF;
					idx++;
				}
			}	
		}

		
		// if too many paths are here sort and pop back
		if (max_num_paths>0 && seq_paths.size()>max_num_paths)
		{
			sort(seq_paths.begin(),seq_paths.end(),comp_SeqPath_path_score);
			while (seq_paths.size()>max_num_paths && seq_paths[seq_paths.size()-1].path_score == NEG_INF)
				seq_paths.pop_back();
		}
		else // remove all seq paths with NEG_INF score
		{
			int last_idx = seq_paths.size()-1;
			int i;
			for (i=0; i<last_idx; i++)
			{
				while (last_idx>i && seq_paths[last_idx].path_score == NEG_INF)
					last_idx--;
				
				if (seq_paths[i].path_score == NEG_INF)
				{
					seq_paths[i]=seq_paths[last_idx];
					seq_paths[last_idx].path_score = NEG_INF;
				}
			}

			while (seq_paths.size()>0 && seq_paths[seq_paths.size()-1].path_score == NEG_INF)
				seq_paths.pop_back();
		}
	}

	const MultiEdge& last_edge = multi_edges[multi_path.edge_idxs[multi_path.edge_idxs.size()-1]];
	PathPos last_pos;

	last_pos.breakage = last_edge.c_break;
	last_pos.edge_idx =-1;
	last_pos.mass = last_edge.c_break->mass;
	last_pos.node_idx = last_edge.c_idx;
	last_pos.node_score = last_edge.c_break->score;


	for (i=0; i<seq_paths.size(); i++)
	{
		SeqPath& path= seq_paths[i];
		path.path_score += last_pos.node_score;
		path.num_forbidden_nodes = multi_path.num_forbidden_nodes;
		path.path_score -= (forbidden_pair_penalty * path.num_forbidden_nodes);
		path.positions.push_back(last_pos);
		path.make_seq_str(config);
	}
}

struct var_combo {
	var_combo() : path_score(0) {};
	bool operator< (const var_combo& other) const
	{
		return (path_score>other.path_score);
	}

	vector<int> var_idxs;
	vector<score_t> node_scores; // dim = #var_idxs+1
	score_t path_score;
};






/************************************************************************

  Expands the single multi path into all possible sequence variants.
  Since this turns out to be the time-limiting process for long de nvoo,
  this is implemented using a quick branch and bound that works on
  edge variant indices. The SeqPaths are created only for the final
  set of sequences.

*************************************************************************/
void PrmGraph::fast_expand_multi_path(const MultiPath& multi_path, 
									  vector<SeqPath>& seq_paths,
									  score_t min_score,
									  score_t forbidden_pair_penalty,
									  int max_num_paths) const
{
	const vector<AA_combo>& aa_edge_comobos = config->get_aa_edge_combos();
	vector<score_t> max_attainable_scores;

	const int num_edges = multi_path.edge_idxs.size();
	max_attainable_scores.resize(num_edges+1,0);
	
	int num_path_variants=1;
	int e;
	for (e=num_edges-1; e>=0; e--)
	{
		const int e_idx = multi_path.edge_idxs[e];
		const MultiEdge& edge = multi_edges[e_idx];

		max_attainable_scores[e]= max_attainable_scores[e+1] + 
			edge.max_variant_score + edge.c_break->score; 
		num_path_variants *= edge.variant_ptrs.size();
	}

	const score_t max_path_score =  multi_edges[multi_path.edge_idxs[0]].n_break->score + 
									max_attainable_scores[0];
	const score_t min_delta_allowed = min_score - max_path_score;

	if (min_delta_allowed>0)
		return; // this multipath won't help much

	if (num_path_variants == 0)
	{
		cout << "Error: had an edge with 0 variants!" <<endl;
		exit(1);
	}

	// perform expansion using a heap and condensed path representation
	vector<int> var_positions;				   // holds the edge idxs
	vector<int> num_vars;
	vector<int> var_edge_positions;
	vector< vector< score_t > > var_scores;
	var_scores.clear();
	var_positions.clear();
	num_vars.clear();
	int num_aas_in_multipath = 0;
	for (e=0; e<num_edges; e++)
	{
		const int e_idx = multi_path.edge_idxs[e];
		const MultiEdge& edge = multi_edges[e_idx];
		const int num_vars_in_edge = edge.variant_ptrs.size();
		vector<score_t> scores;

		num_aas_in_multipath+=edge.num_aa;

		if (num_vars_in_edge>1)
		{
			int i;
			for (i=0; i<num_vars_in_edge; i++)
				scores.push_back(edge.variant_scores[i]);
		
			var_scores.push_back(scores);
			var_positions.push_back(e);
			num_vars.push_back(num_vars_in_edge);
			var_edge_positions.push_back(num_aas_in_multipath-edge.num_aa);

		//	cout << e << "\t" << e_idx << "\t" << num_vars_in_edge << "\t" << edge.num_aa << endl;
		}	
	}

		// create the SeqPaths from the edge_combos...
	const int num_positions = num_aas_in_multipath+1;
	SeqPath template_path;				// holds common elements to all paths
	template_path.n_term_aa = multi_path.n_term_aa;
	template_path.c_term_aa = multi_path.c_term_aa;
	template_path.n_term_mass = multi_path.n_term_mass;
	template_path.c_term_mass = multi_path.c_term_mass;
	template_path.multi_path_rank = multi_path.original_rank;
	template_path.path_score = 0;
	template_path.prm_ptr = (PrmGraph *)this;
	template_path.positions.clear();
	template_path.positions.resize(num_positions);

	template_path.num_forbidden_nodes = multi_path.num_forbidden_nodes;
	template_path.path_score -= template_path.num_forbidden_nodes * forbidden_pair_penalty;

	int pos_idx=0;
	for (e=0; e<multi_path.edge_idxs.size(); e++)
	{
		const int e_idx = multi_path.edge_idxs[e];
		const MultiEdge& edge = get_multi_edge(e_idx);
		PathPos& pos = template_path.positions[pos_idx++];

		pos.breakage = edge.n_break;
		pos.edge_idx = e_idx;
		pos.mass =  edge.n_break->mass;
		pos.node_score = edge.n_break->score;
		pos.node_idx = edge.n_idx;
		template_path.path_score += pos.node_score;

		int *variant_ptr = edge.variant_ptrs[0];
		const int num_aa = *variant_ptr++;
		int *aas = variant_ptr;

		if (edge.variant_ptrs.size() == 1)
		{
			pos.edge_variant_score = edge.variant_scores[0];
			pos.aa = aas[0];
			template_path.path_score += pos.edge_variant_score;
		}

		if (edge.num_aa == 1)
			continue;
		
		int j;
		for (j=1; j<edge.num_aa; j++)
		{
			PathPos& pos = template_path.positions[pos_idx++];
			pos.breakage = NULL;
			pos.aa = aas[j];
		}
	}

	if (pos_idx != num_aas_in_multipath)
	{
		cout << "Error: mismatches between positions and multipath length!" << endl;
		exit(1);
	}
	
	PathPos& last_pos = template_path.positions[num_positions-1];
	const int last_e_idx = multi_path.edge_idxs[e-1];
	const MultiEdge& last_edge = get_multi_edge(last_e_idx);
	last_pos.breakage = last_edge.c_break;
	last_pos.edge_idx =-1;
	last_pos.mass = last_edge.c_break->mass;
	last_pos.node_idx = last_edge.c_idx;
	last_pos.node_score = last_edge.c_break->score;

	template_path.path_score += last_pos.node_score;

	const int num_multi_var_edges = var_positions.size();
	if (num_multi_var_edges == 0)
	{
		seq_paths.resize(1);
		seq_paths[0]=template_path;
		return;
	}

	vector<var_combo>      combo_heap;
	vector<int> idxs;
	idxs.resize(num_multi_var_edges,0);
	const int last_var_pos = num_multi_var_edges-1;
	const int last_var_val = num_vars[last_var_pos];
	const score_t needed_var_score = min_score - template_path.path_score;
	while (1)
	{
		score_t idxs_score = 0;
		int k;
		for (k=0; k<idxs.size(); k++)
			idxs_score += var_scores[k][idxs[k]];
		
		if (idxs_score>=needed_var_score)
		{
			if (combo_heap.size()<max_num_paths)
			{
				var_combo v;
				v.path_score = template_path.path_score + idxs_score;
				v.var_idxs = idxs;
				combo_heap.push_back(v);
				if (combo_heap.size()== max_num_paths)
					make_heap(combo_heap.begin(),combo_heap.end());
			}
			else
			{
				const score_t score = template_path.path_score + idxs_score;
				if (score>combo_heap[0].path_score)
				{
					pop_heap(combo_heap.begin(),combo_heap.end());
					combo_heap[max_num_paths-1].path_score = score;
					combo_heap[max_num_paths-1].var_idxs = idxs;
					push_heap(combo_heap.begin(),combo_heap.end());
				}
			}	
		}

		int j=0;
		while (j<num_multi_var_edges)
		{
			idxs[j]++;
			if (idxs[j]==num_vars[j])
			{
				idxs[j++]=0;
			}
			else
				break;
		}
		if (j==num_multi_var_edges)
			break;

	}

	seq_paths.clear();
	seq_paths.resize(combo_heap.size(),template_path);

	int i;
	for (i=0; i<combo_heap.size(); i++)
	{
		const var_combo& combo  = combo_heap[i];
		SeqPath& seq_path = seq_paths[i];

		// fill in info that changes with multi-variant edges
		int j;
		for (j=0; j<num_multi_var_edges; j++)
		{
			const int var_pos = var_positions[j];
			const int e_idx = multi_path.edge_idxs[var_pos];
			const MultiEdge& edge = get_multi_edge(e_idx);

			const int pos_idx = var_edge_positions[j];
		
			if (seq_path.positions[pos_idx].edge_idx !=  e_idx)
			{
				cout << "POS: " << pos_idx << endl;
				cout << "Error: mismatch in pos_idx and e_idx of a multipath!" << endl;
				cout << "looking for " << e_idx << endl;
				cout << "edge idxs:" << endl;
				int k;
				for (k=0; k<seq_path.positions.size(); k++)
					cout << k << "\t" << seq_path.positions[k].edge_idx << "\tnidx: " <<
					seq_path.positions[k].node_idx << endl;
				cout << endl;
				multi_path.print(config);

				this->print();

				exit(1);
			}

			PathPos& pos = seq_path.positions[pos_idx];
			const int variant_idx = combo.var_idxs[j];
			int *variant_ptr = edge.variant_ptrs[variant_idx];
			const int num_aa = *variant_ptr++;
			int *aas = variant_ptr;

			if (num_aa != edge.num_aa)
			{
				cout << "Error: edge and variant mixup!" << endl;
				exit(1);
			}

			pos.edge_variant_score = edge.variant_scores[variant_idx];
			pos.aa = aas[0];
		//	pos.edge_var_idx = variant_idx;

			if (edge.num_aa == 1)
				continue;

			int k;
			for (k=1; k<edge.num_aa; k++)
			{
				PathPos& pos = seq_path.positions[pos_idx+k];
				pos.aa = aas[k];
			}
		}

		//seq_path.path_score = max_path_score + combo.path_delta;
		seq_path.path_score=combo.path_score;
		seq_path.make_seq_str(config);
	}
}


/************************************************************************

  Expands the single multi path into all possible sequence variants.
  Since this turns out to be the time-limiting process for long de nvoo,
  this is implemented using a quick branch and bound that works on
  edge variant indices. The SeqPaths are created only for the final
  set of sequences.

*************************************************************************/
void PrmGraph::fast_expand_multi_path_for_combo_scores(
									  AllScoreModels *model,
									  const MultiPath& multi_path, 
									  vector<SeqPath>& seq_paths,
									  score_t min_score,
									  score_t forbidden_pair_penalty,
									  int max_num_paths)
{
	const vector<AA_combo>& aa_edge_comobos = config->get_aa_edge_combos();
	vector<score_t> max_attainable_scores;

	const int num_edges = multi_path.edge_idxs.size();
	max_attainable_scores.resize(num_edges+1,0);
	
	int num_path_variants=1;
	int e;
	for (e=num_edges-1; e<num_edges; e++)
	{
		const int e_idx = multi_path.edge_idxs[e];
		const MultiEdge& edge = multi_edges[e_idx];
		num_path_variants *= edge.variant_ptrs.size();
	}


	if (num_path_variants == 0)
	{
		cout << "Error: had an edge with 0 variants!" <<endl;
		exit(1);
	}

	// re-check the max attainable scores, this time compute all combo scores along the path
	// this will give a more accurate max attainable score, and also allow quick computation
	// of the expanded scores

	max_attainable_scores.clear();
	max_attainable_scores.resize(num_edges+1,0);
	vector<score_t> max_node_scores;
	max_node_scores.resize(num_edges+1,NEG_INF);

	for (e=0; e<=num_edges; e++)
	{
		const int n_edge_idx = (e>0 ? multi_path.edge_idxs[e-1] : -1);
		const int c_edge_idx = (e<num_edges ? multi_path.edge_idxs[e] : -1);
		const int node_idx = (n_edge_idx>=0 ? multi_edges[n_edge_idx].c_idx : multi_edges[c_edge_idx].n_idx);
		const int num_n_vars = (e>0 ? multi_edges[n_edge_idx].variant_ptrs.size() : 1);
		const int num_c_vars = (e<num_edges ? multi_edges[c_edge_idx].variant_ptrs.size() : 1);

		int j,k;
		for (j=0; j<num_n_vars; j++)
			for (k=0; k<num_c_vars; k++)
			{
				score_t combo_score = model->get_node_combo_score(this,node_idx,n_edge_idx,j,c_edge_idx,k);
				if (max_node_scores[e]<combo_score)
					max_node_scores[e]=combo_score;
			}
	}

	max_attainable_scores[num_edges]=max_node_scores[num_edges];
	for (e=num_edges-1; e>=0; e--)
	{
		const int e_idx = multi_path.edge_idxs[e];
		const MultiEdge& edge = multi_edges[e_idx];

		max_attainable_scores[e]= max_attainable_scores[e+1] + 
			edge.max_variant_score + max_node_scores[e]; 
	}
	score_t max_path_score =  max_attainable_scores[0] - multi_path.num_forbidden_nodes * forbidden_pair_penalty;
	score_t min_delta_allowed = min_score - max_path_score;

	if (min_delta_allowed>0)
		return;

	// in this expansion method, all edges are stored (not only ones with more than one variant)
	// perform expansion using a heap and condensed path representation
	vector<int> var_positions;				   // holds the edge idxs
	vector<int> num_vars;
	vector<int> var_edge_positions;
	vector< vector< score_t > > var_scores;
	var_scores.clear();
	var_positions.clear();
	num_vars.clear();

	int num_aas_in_multipath = 0;
	for (e=0; e<num_edges; e++)
	{
		const int e_idx = multi_path.edge_idxs[e];
		const MultiEdge& edge = multi_edges[e_idx];
		const int num_vars_in_edge = edge.variant_ptrs.size();
		vector<score_t> scores;

		num_aas_in_multipath+=edge.num_aa;

		int i;
		for (i=0; i<num_vars_in_edge; i++)
			scores.push_back(edge.variant_scores[i]);
		
		var_scores.push_back(scores);
		var_positions.push_back(e);
		num_vars.push_back(num_vars_in_edge);
		var_edge_positions.push_back(num_aas_in_multipath-edge.num_aa);
	}

	// create the SeqPaths from the edge_combos...
	const int num_positions = num_aas_in_multipath+1;
	SeqPath template_path;				// holds common elements to all paths
	template_path.n_term_aa = multi_path.n_term_aa;
	template_path.c_term_aa = multi_path.c_term_aa;
	template_path.n_term_mass = multi_path.n_term_mass;
	template_path.c_term_mass = multi_path.c_term_mass;
	template_path.charge = charge;
	template_path.multi_path_rank = multi_path.original_rank;
	template_path.path_score = 0;
	template_path.prm_ptr = (PrmGraph *)this;
	template_path.positions.clear();
	template_path.positions.resize(num_positions);

	template_path.num_forbidden_nodes = multi_path.num_forbidden_nodes;
	template_path.path_score -= template_path.num_forbidden_nodes * forbidden_pair_penalty;

	int pos_idx=0;
	for (e=0; e<multi_path.edge_idxs.size(); e++)
	{
		const int e_idx = multi_path.edge_idxs[e];
		const MultiEdge& edge = get_multi_edge(e_idx);
		PathPos& pos = template_path.positions[pos_idx++];

		pos.breakage = edge.n_break;
		pos.edge_idx = e_idx;
		pos.mass =  edge.n_break->mass;
		pos.node_idx = edge.n_idx;
		
		int *variant_ptr = edge.variant_ptrs[0];
		const int num_aa = *variant_ptr++;
		int *aas = variant_ptr;

		if (edge.variant_ptrs.size() == 1)
		{
			pos.variant_ptr = edge.variant_ptrs[0];
			pos.edge_variant_score = edge.variant_scores[0];
			pos.aa = aas[0];
		}

		if (edge.num_aa == 1)
			continue;
		
		int j;
		for (j=1; j<edge.num_aa; j++)
		{
			PathPos& pos = template_path.positions[pos_idx++];
			pos.breakage = NULL;
			pos.aa = aas[j];
		}
	}

	if (pos_idx != num_aas_in_multipath)
	{
		cout << "Error: mismatches between positions and multipath length!" << endl;
		exit(1);
	}
	
	PathPos& last_pos = template_path.positions[num_positions-1];
	const int last_e_idx = multi_path.edge_idxs[e-1];
	const MultiEdge& last_edge = get_multi_edge(last_e_idx);
	last_pos.breakage = last_edge.c_break;
	last_pos.edge_idx =-1;
	last_pos.mass = last_edge.c_break->mass;
	last_pos.node_idx = last_edge.c_idx;
	

	const int num_multi_var_edges = var_positions.size();
	
	vector<var_combo>      combo_heap;
	vector<int> var_idxs;
	var_idxs.resize(num_multi_var_edges,0);
	const int last_var_pos = num_multi_var_edges-1;
	const int last_var_val = num_vars[last_var_pos];
	const score_t needed_var_score = min_score; // this is only the penalty for forbidden pairs
	while (1)
	{
		score_t idxs_score = 0;
		int k;
		for (k=0; k<var_idxs.size(); k++)
			idxs_score += var_scores[k][var_idxs[k]]; // this is the score for the edge variants (usually 0) 
		
		const vector<PathPos>& positions = template_path.positions;
		score_t path_score = idxs_score + template_path.path_score;
		
		// compute the node scores based on the edge variants
		// N-term node score
		vector<score_t> node_scores;
		node_scores.resize(var_idxs.size()+1,NEG_INF);

		node_scores[0] = model->get_node_combo_score(this,positions[0].node_idx,-1,0,
													positions[0].edge_idx,var_idxs[0]);
		path_score += node_scores[0];
		// middle nodes score
		int prev_edge_idx = positions[0].edge_idx;
		for (k=1; k<var_idxs.size(); k++)
		{
			const int p=var_edge_positions[k];
			const PathPos& pos = positions[p];
			node_scores[k] = model->get_node_combo_score(this,pos.node_idx,prev_edge_idx,var_idxs[k-1],
													pos.edge_idx,var_idxs[k]);
			path_score += node_scores[k];
			prev_edge_idx = pos.edge_idx;
		}

		// C-term node score
		node_scores[k] = model->get_node_combo_score(this,last_pos.node_idx,prev_edge_idx,var_idxs[k-1],-1,0);

		path_score += node_scores[k];
		
		if (path_score>=needed_var_score)
		{
			if (combo_heap.size()<max_num_paths)
			{
				var_combo v;
				v.path_score  = path_score;
				v.var_idxs    = var_idxs;
				v.node_scores = node_scores;
				combo_heap.push_back(v);
				if (combo_heap.size()== max_num_paths)
					make_heap(combo_heap.begin(),combo_heap.end());
			}
			else
			{
				const score_t score = path_score;
				if (score>combo_heap[0].path_score)
				{
					pop_heap(combo_heap.begin(),combo_heap.end());
					var_combo& combo  = combo_heap[max_num_paths-1];
					combo.path_score  = path_score;
					combo.var_idxs    = var_idxs;
					combo.node_scores = node_scores;
					push_heap(combo_heap.begin(),combo_heap.end());
				}
			}	
		}

		int j=0;
		while (j<num_multi_var_edges)
		{
			var_idxs[j]++;
			if (var_idxs[j]==num_vars[j])
			{
				var_idxs[j++]=0;
			}
			else
				break;
		}
		if (j==num_multi_var_edges)
			break;

	}


	seq_paths.clear();
	seq_paths.resize(combo_heap.size(),template_path);

	int i;
	for (i=0; i<combo_heap.size(); i++)
	{
		const var_combo& combo  = combo_heap[i];
		SeqPath& seq_path = seq_paths[i];

		// fill in info that changes with multi-variant edges
		int j;
		for (j=0; j<num_multi_var_edges; j++)
		{
			const int var_pos = var_positions[j];
			const int e_idx = multi_path.edge_idxs[var_pos];
			const MultiEdge& edge = get_multi_edge(e_idx);

			const int pos_idx = var_edge_positions[j];
		
			if (seq_path.positions[pos_idx].edge_idx !=  e_idx)
			{
				cout << "POS: " << pos_idx << endl;
				cout << "Error: mismatch in pos_idx and e_idx of a multipath!" << endl;
				cout << "looking for " << e_idx << endl;
				cout << "edge idxs:" << endl;
				int k;
				for (k=0; k<seq_path.positions.size(); k++)
					cout << k << "\t" << seq_path.positions[k].edge_idx << "\tnidx: " <<
					seq_path.positions[k].node_idx << endl;
				cout << endl;
				multi_path.print(config);

				this->print();

				exit(1);
			}

			PathPos& pos = seq_path.positions[pos_idx];
			const int variant_idx = combo.var_idxs[j];
			int *variant_ptr = edge.variant_ptrs[variant_idx];
			const int num_aa = *variant_ptr++;
			int *aas = variant_ptr;

			if (num_aa != edge.num_aa)
			{
				cout << "Error: edge and variant mixup!" << endl;
				exit(1);
			}

			pos.node_score = combo.node_scores[j];
			pos.edge_variant_score = edge.variant_scores[variant_idx];
			pos.variant_ptr = edge.variant_ptrs[variant_idx];
			pos.aa = aas[0];
	
			if (edge.num_aa == 1)
				continue;

			int k;
			for (k=1; k<edge.num_aa; k++)
			{
				PathPos& pos = seq_path.positions[pos_idx+k];
				pos.aa = aas[k];
			}
		}

		seq_path.positions[seq_path.positions.size()-1].node_score = combo.node_scores[num_multi_var_edges];
		seq_path.path_score=combo.path_score;
		float pen = seq_path.adjust_complete_sequence_penalty(this);

		seq_path.make_seq_str(config);
	}
}


void PrmGraph::expand_all_multi_paths(AllScoreModels *model, 
									  const vector<MultiPath>& multi_paths, 
									  vector<SeqPath>& paths, 
									  score_t forbidden_pair_penalty,
									  int max_num_paths)
{
	// expand all variants
	PathHeap path_heap;

	path_heap.init(max_num_paths);
	paths.clear();


	int mp_idx;
	for (mp_idx=0; mp_idx<multi_paths.size(); mp_idx++ )
	{
		vector<SeqPath> variants;

		if (multi_paths[mp_idx].path_score== NEG_INF)
			continue;

		if (this->has_node_combo_scores)
		{
			fast_expand_multi_path_for_combo_scores(model, 
													multi_paths[mp_idx], 
													variants, 
													path_heap.get_min_score(), 
													forbidden_pair_penalty, 
													max_num_paths);
		}
		else
		{
			if (multi_paths[mp_idx].path_score<=path_heap.get_min_score())
				continue;

			fast_expand_multi_path(multi_paths[mp_idx], variants, path_heap.get_min_score(), 
				forbidden_pair_penalty, max_num_paths);
		}
		
		if (variants.size()<=0)
			continue;

		sort(variants.begin(),variants.end(),comp_SeqPath_path_score);

	/*	cout << mp_idx << "\t" << multi_paths[mp_idx].path_score << "\t" << variants.size() << endl;
		int q;
		for (q=0; q<variants.size(); q++)
			cout <<  "  " << variants[q].seq_str << "\t" << variants[q].path_score << endl;
	*/
		int j;
		for (j=0; j<variants.size(); j++)
		{
			variants[j].multi_path_rank = mp_idx;
			variants[j].multi_path_score = multi_paths[mp_idx].path_score;

			// check that the peptide has correct mass if it spans the entire graph

			int n_idx = variants[j].positions[0].node_idx;
			int c_idx = variants[j].positions[ variants[j].positions.size()-1].node_idx;

			if (nodes[n_idx].type == NODE_N_TERM &&	nodes[c_idx].type == NODE_C_TERM) 
			{
				if (! config->get_need_to_estimate_pm() )
				{
					const vector<mass_t>& aa2mass = config->get_aa2mass();
					vector<int> aas;
					variants[j].get_amino_acids(aas);
					mass_t pep_mass = 0;
					int a;
					for (a=0; a<aas.size(); a++)
						pep_mass += aa2mass[aas[a]];

					pep_mass+=MASS_OHHH;

					if (fabs(pep_mass-source_spectrum->get_org_pm_with_19())>config->get_pm_tolerance())
					{
					//	cout << "Rejected: " << variants[j].seq_str << endl;
					//	cout << "aa_mass: " << fixed << setprecision(3) << pep_mass << " spec_mass: " << source_spectrum->get_org_pm_with_19();
					//	cout << " (" << pep_mass-source_spectrum->get_org_pm_with_19() << ")" << endl;

						continue;
					}
				}
			}

			// check that there are no PTMs that are specific to the +1, -1 positions

	
			path_heap.add_path(variants[j]);
		}
	}

	paths = path_heap.get_paths();
	sort(paths.begin(),paths.end(),comp_SeqPath_path_score);

	while (paths.size()>0)
	{
		int idx = paths.size()-1;
		if (paths[idx].path_score<= NEG_INF || paths[idx].get_num_aa() < 1)
		{
			paths.pop_back();
		}
		else
			break;
	}

//	int i;
//	for (i=0; i<paths.size(); i++)
//		cout << i << "\t" << paths[i].path_score << "\t" << paths[i].seq_str << endl;

//	cout << "Considered: " << num_vars << endl;
}



bool MultiPath::check_if_correct(const Peptide& p, Config *config) const
{
	const mass_t tolerance = config->getTolerance() * 1.25;
	vector<mass_t> break_masses;
	int idx=0;
	int i;

	p.calc_expected_breakage_masses(config,break_masses);

	for (i=0; i<breakages.size(); i++)
	{
		const mass_t& mass = breakages[i]->mass;
		const mass_t max_mass = mass + tolerance;
		const mass_t min_mass = mass - tolerance;

		while (idx < break_masses.size() && break_masses[idx] < min_mass)
			idx++;

		if (break_masses[idx]>max_mass)
			return false;
	}

	return true;
}



int  MultiPath::get_num_correct_aas(const PrmGraph& prm, const Peptide& p, Config *config) const
{
	const mass_t tolerance = config->getTolerance() * 1.25;
	vector<mass_t> break_masses;
	int idx=0;
	int num_correct=0;
	int i;

	p.calc_expected_breakage_masses(config,break_masses);

	for (i=0; i<breakages.size(); i++)
	{
		const mass_t& mass = breakages[i]->mass;
		const mass_t max_mass = mass + tolerance;
		const mass_t min_mass = mass - tolerance;

		while (idx < break_masses.size() && break_masses[idx] < min_mass)
			idx++;

		if (break_masses[idx]>max_mass)
			continue;
		
		if (idx<breakages.size()-1 && edge_idxs[idx]>=0)
			num_correct += prm.get_multi_edge(edge_idxs[idx]).num_aa;
	}
	return num_correct;
}


int  MultiPath::get_num_aas() const
{
	return (breakages.size()-1);
}


// returns the number of b/y ions
int SeqPath::get_num_frags(const vector<int>& frag_idxs) const
{
	int num_frags=0;
	int i;
	for (i=0; i<positions.size(); i++)
	{
		Breakage *bb = positions[i].breakage;

		if (positions[i].breakage && positions[i].breakage->fragments.size()>0)
		{
			int j;
			for (j=0; j<frag_idxs.size(); j++)
			{
				int k;
				for (k=0; k<positions[i].breakage->fragments.size(); k++)
				{
					if (positions[i].breakage->fragments[k].frag_type_idx == frag_idxs[j])
						num_frags++;
				}
			}
		}
	}
	return num_frags;
}


mass_t SeqPath::calculate_peptide_mass(const Config *config) const
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	mass_t pep_mass=0;

	int i;
	for (i=0; i<positions.size()-1; i++)
	{
		pep_mass+=aa2mass[positions[i].aa];
	}
	return pep_mass;
}

int SeqPath::get_num_correct_aas(const Peptide& pep, const Config *config) const
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const vector<int>& pep_aas = pep.get_amino_acids();

	int num_correct=0;
	int i;

	vector<mass_t> pep_masses;
	vector<int> path_aas;
	
	get_amino_acids(path_aas);

	pep_masses.resize(pep_aas.size(),0);
	for (i=1; i<pep_aas.size(); i++)
		pep_masses[i]=pep_masses[i-1]+aa2mass[pep_aas[i-1]];

	mass_t path_mass = n_term_mass;
	for (i=0; i<path_aas.size(); i++)
	{
		const int path_aa = path_aas[i];
		int j;
		for (j=0; j<pep_aas.size(); j++)
		{
			const int pep_aa = pep_aas[j];
			if (fabs(pep_masses[j]-path_mass)<1.0 && pep_aas[j] == path_aas[i])
			{
				num_correct++;
				break;
			}
		}

		path_mass += aa2mass[path_aas[i]];
	}

	return num_correct;
}





void PrmGraph::parse_seq_path_to_smaller_ones(const SeqPath& org_path, 
										  int min_length, 
										  int max_length, 
										  vector<SeqPath>& new_paths)
{
	const vector<PathPos>& org_positions = org_path.positions;
	new_paths.clear();

	
	const int num_org_positions = org_positions.size();
	const int last_org_position = num_org_positions-1;
	if (num_org_positions<=min_length)
		return;

	const int num_iters = org_positions.size()-min_length;
	new_paths.resize(num_iters);

	int np_idx=0;

	int i;
	for (i=0; i<num_iters; i++)
	{
		if (org_positions[i].aa>=0 && org_positions[i].edge_idx>=0)
		{
			int j=0;
			SeqPath& path = new_paths[np_idx];
			path.positions.clear();

			if (i==0)
				path.n_term_aa = org_path.n_term_aa;

			path.n_term_mass = org_positions[i].mass;
			path.path_score = 0;
			path.multi_path_rank = org_path.multi_path_rank;
			
			const int start_pos=i;
			bool score_this_path = false;
			int pos=i;
			while (pos<num_org_positions )
			{
				path.positions.push_back(org_positions[pos]);
				pos++;
			
				while (pos<num_org_positions  && org_positions[pos].node_idx<0)
					path.positions.push_back(org_positions[pos++]);
				
				const int length = pos-start_pos;

				if (pos<num_org_positions && length>=min_length && length<=max_length)
				{
					path.c_term_mass = org_positions[pos].mass;
					if (pos==last_org_position)
						path.c_term_aa= org_path.c_term_aa;

					path.positions.push_back(org_positions[pos]);
					score_this_path = true;
					break;
				}
			}	


			if (path.positions.size()>max_length+1 || ! score_this_path)
			{
				path.positions.clear();
				continue;
			}

			np_idx++; // stores the path

			// compute combo scores
			vector<PathPos>& new_positions = path.positions;
			const int num_new_positions = new_positions.size();
			const int last_new_position = num_new_positions-1;

			path.path_score=0;

			// special treatment for scores at the begining/end of the parsed tag
			// (need to use the Gap score on the terminal ends)
			if (new_positions[0].node_idx == org_positions[0].node_idx)
			{
				new_positions[0].node_score = org_positions[0].node_score;
				path.path_score += new_positions[0].node_score; 
				path.path_score += new_positions[0].edge_variant_score;
			}
			else
			{
				// find edge variant
				const int aa = new_positions[0].aa;
				const MultiEdge& edge = multi_edges[new_positions[0].edge_idx];
				const int var_idx = edge.get_variant_idx(aa);
				if (var_idx<0)
				{
					cout << "Error: bad aa in var idx! 1" << endl;
					exit(1);
				}
				new_positions[0].node_score = model->get_node_combo_score(this,new_positions[0].node_idx,
														-1,0,new_positions[0].edge_idx,var_idx);

				path.path_score += new_positions[0].node_score;
				path.path_score += new_positions[0].edge_variant_score;
			}

			if (new_positions[last_new_position].node_idx == org_positions[last_org_position].node_idx)
			{
				new_positions[last_new_position].node_score = org_positions[last_org_position].node_score;
				path.path_score += new_positions[last_new_position].node_score;
				path.path_score += new_positions[last_new_position].edge_variant_score;
			}
			else
			{
				int k=last_new_position-1;
				while (k>0 && new_positions[k].edge_idx<0)
					k--;

				if (k==0)
				{
					cout << "Error: bad parse of tag!" << endl;
					exit(1);
				}

				const int aa = new_positions[k].aa;
				const int edge_idx = new_positions[k].edge_idx;
				const MultiEdge& edge = multi_edges[edge_idx];
				const int var_idx = edge.get_variant_idx(aa);
				if (var_idx<0)
				{
					cout << "Error: bad aa in var idx! 2" << endl;
					exit(1);
				}

				new_positions[last_new_position].node_score = model->get_node_combo_score(this,
					new_positions[last_new_position].node_idx,edge_idx,var_idx,-1,0);

				path.path_score += new_positions[last_new_position].node_score;
				////	no edge score for last position!
				//path.path_score += new_positions[last_new_position].edge_variant_score;
				new_positions[last_new_position].edge_variant_score =0;

			}

			for (j=1; j<last_new_position; j++)
			{
				path.positions[j].node_score = org_positions[i+j].node_score;
				path.path_score += new_positions[j].node_score; 
				path.path_score += new_positions[j].edge_variant_score;
			}
			
			path.make_seq_str(config);
			path.charge = charge;
			path.pm_with_19 = pm_with_19;
			path.prm_ptr = this;
			path.sort_key = path.path_score;
		}
	}

	while(new_paths.size()>0 && new_paths[new_paths.size()-1].positions.size()==0)
		new_paths.pop_back();
}



bool SeqPath::check_if_correct(const string& str, const Config *config) const
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const char *path_str = seq_str.c_str();
	const char *corr_str = str.c_str();

	int len_path_str = strlen(path_str);
	int len_corr_str = strlen(corr_str);

	if (len_path_str>len_corr_str)
		return false;

	
	int i;
	for (i=0; i<=len_corr_str-len_path_str; i++)
	{
		int j;
		bool correct_seq = true;
		for (j=0; j<len_path_str; j++)
			if (! (path_str[j] == corr_str[i+j] ||
				  (path_str[j] == 'I' && corr_str[i+j]== 'L') ||
				  (path_str[j] == 'L' && corr_str[i+j]== 'I') ||
				  (path_str[j] == 'Q' && corr_str[i+j]== 'K') ||
				  (path_str[j] == 'K' && corr_str[i+j]== 'Q') ) )
			{
				correct_seq = false;
				break;
			}



		if (correct_seq)
		{
			// check prefix mass
			Peptide pep;
			pep.parseFromString(config,corr_str);
			const vector<int>& aas= pep.get_amino_acids();
			mass_t mass=0;
			int j;

			if (n_term_mass == 0 && i==0)
				return true;

			for (j=0; j<aas.size(); j++)
			{
				mass+=aa2mass[aas[j]];
				if (fabs(mass-this->n_term_mass)<6)
					return true;

				if (mass>n_term_mass)
					break;
			}
		}
	}

	
	
	return false;
}



bool SeqPath::check_if_cut_correct(const vector<mass_t>& exp_cut_masses, mass_t tolerance) const
{
	int c_idx=0;
	int i;

	for (i=0; i<positions.size()-1; i++)
	{
		if (positions[i].node_idx<0)
			continue;

		const mass_t pos_mass = positions[i].mass;
		while (c_idx<exp_cut_masses.size() && exp_cut_masses[c_idx]< pos_mass - tolerance)
			c_idx++;

		if (c_idx == exp_cut_masses.size() || exp_cut_masses[c_idx]-tolerance > pos_mass)
			return false;

	}
	return true;
}





void SeqPath::make_seq_str(const Config *config)
{
	const vector<string>& aa2label = config->get_aa2label();
	int i;

	seq_str = "";

	if (n_term_aa>N_TERM)
		seq_str += aa2label[n_term_aa] ;
	
	if (positions.size()>0)
		for (i=0; i<positions.size()-1; i++)
			seq_str += aa2label[positions[i].aa];

	if (c_term_aa>C_TERM)
		seq_str +=  aa2label[c_term_aa];

}

void MultiPath::print(Config *config, ostream& os) const
{
	os << "MultiPath: " << n_term_mass << " - " << c_term_mass << 
		   " score: " << path_score << "  ";
	int i;

	cout << " Nodes:";
	for (i=0; i<node_idxs.size(); i++)
		cout << " " << node_idxs[i];
	cout << "   Edges:";
	for (i=0; i<edge_idxs.size(); i++)
		cout << " " << edge_idxs[i];
	cout << endl;
}


void SeqPath::print(ostream& os) const
{
	os << setprecision(5);
	os << n_term_mass << " " << seq_str << " " << c_term_mass << " (s: " << 
		this->path_score << ")" << endl;
}


void SeqPath::print_with_probs(ostream& os) const
{
	os << setprecision(5);
	os << n_term_mass << " " << seq_str << " (s: " << 
		this->path_score << ")";

	int i;
	os << setprecision(2);
	for (i=0; i<this->positions.size(); i++)
	{
		os << " X " ;
	}
	os << endl;
}

void SeqPath::print_full(const Config *config, ostream &os) const
{
	const vector<string>& aa2label = config->get_aa2label();
	os << n_term_mass << " " << seq_str << " " << c_term_mass << " (s: " << 
		this->path_score << ")  NF: " << this->num_forbidden_nodes << endl;
	int i;

	if (positions.size()==0)
		return;

	for (i=0; i<positions.size()-1; i++)
	{
		cout << left << setw(3) << i;
		cout << setw(5) << left << aa2label[positions[i].aa] << "\t" <<
			positions[i].node_idx << "\t" << positions[i].edge_idx<< "\t" << 
			fixed << setprecision(2) << positions[i].mass << "\t" << 
			positions[i].node_score << "\t" << positions[i].edge_variant_score << endl;
	}
	cout << "Exact score: " << path_score << "\t" << "Approx score: " << setprecision(4) << multi_path_score << endl;
		
}
