#include "DeNovoDp.h"
#include "AnnotatedSpectrum.h"
#include "PepNovo_auxfun.h" 

struct dis_pair  {
	dis_pair() : dis(0), left_idx(-1), right_idx(-1) {};
	
	bool operator< (dis_pair& other)
	{
		return dis<other.dis;
	}

	bool operator< (dis_pair& other) const
	{
		return dis<other.dis;
	} 

	bool operator< (const dis_pair& other)
	{
		return dis<other.dis;
	}

	bool operator< (const dis_pair& other) const
	{
		return dis<other.dis;
	}

	mass_t dis;
	
	int left_idx,right_idx;
};




// fills in all cells in the dp_table according to the PrmGraph
void DeNovoDp::fill_dp_table(const PrmGraph *_prm, score_t sym_penalty)
{
	prm = (PrmGraph *)_prm;
	config=prm->config;
	const vector<Node>& nodes = prm->nodes;
	const vector<MultiEdge>& edges = prm->multi_edges;
	const int num_nodes = nodes.size();

	// forbidden window of 1 Daltons
	const mass_t pm_with_19 = prm->pm_with_19;
	const mass_t min_forbidden_sum=pm_with_19 - MASS_PROTON - config->getTolerance()*1.5;
	const mass_t max_forbidden_sum=pm_with_19 - MASS_PROTON + config->getTolerance()*1.5;
	const mass_t sym_axis = (pm_with_19 - 1.0) *0.5;

	int i,j;
	cells.clear();
	cells.resize(num_nodes);
	for (i=0; i<num_nodes; i++)
		cells[i].resize(num_nodes);

	vector<bool> skip_ind;
	skip_ind.resize(num_nodes,false);
	for (i=0; i<num_nodes; i++)
		if (nodes[i].in_edge_idxs.size()  == 0 && 
			nodes[i].out_edge_idxs.size() == 0)
			skip_ind[i]=true;

	forbidden_idxs.resize(num_nodes,-1);
	vector<dis_pair> pairs;
	for (i=0; i<num_nodes-1; i++)
		for (j=i+1; j<num_nodes; j++)
		{
			dis_pair p;
			p.dis = nodes[j].mass - nodes[i].mass;
			p.left_idx=i;
			p.right_idx=j;

			if (p.dis< 56.0)
				continue;

			pairs.push_back(p);

			// mark forbidden pairs
			mass_t sum = nodes[j].mass + nodes[i].mass;
			if (sum>min_forbidden_sum && sum<max_forbidden_sum)
			{
				cells[i][j].is_forbidden=1;

//				cout << "F: " << i << " " << j << endl;

				if (forbidden_idxs[j]<0)
				{
					forbidden_idxs[j]=i;
					forbidden_idxs[i]=j;
				}
			}
		}

	prm->set_frobidden_node_idxs(forbidden_idxs);

	sort(pairs.begin(),pairs.end());

	// init diagonal
	for (i=0; i<num_nodes; i++)
		cells[i][i].score = nodes[i].score;

	// fill other cells
	for (i=0; i<pairs.size(); i++)
	{
		const int left_idx = pairs[i].left_idx;
		const int right_idx = pairs[i].right_idx;

		if (sym_axis - nodes[left_idx].mass > nodes[right_idx].mass - sym_axis)
		{
			int j;
			for (j=0; j<nodes[left_idx].out_edge_idxs.size(); j++)
			{
				int e_idx = nodes[left_idx].out_edge_idxs[j];
				score_t s = nodes[left_idx].score + edges[e_idx].max_variant_score + 
							cells[edges[e_idx].c_idx][right_idx].score;
				if (s > cells[left_idx][right_idx].score)
				{
					cells[left_idx][right_idx].score = s;
					cells[left_idx][right_idx].prev_edge_idx = e_idx;
				}
			}
		}
		else  // fill the right side
		{
			int j;
			for (j=0; j<nodes[right_idx].in_edge_idxs.size(); j++)
			{
				int e_idx = nodes[right_idx].in_edge_idxs[j];
				score_t s = nodes[right_idx].score + edges[e_idx].max_variant_score+ 
							cells[left_idx][edges[e_idx].n_idx].score;
				if (s > cells[left_idx][right_idx].score)
				{
					cells[left_idx][right_idx].score = s;
					cells[left_idx][right_idx].prev_edge_idx = e_idx;
				}
			}
		}

		if (cells[left_idx][right_idx].is_forbidden)
			cells[left_idx][right_idx].score -= sym_penalty;
	}


	
}



MultiPath DeNovoDp::get_top_scoring_antisymetric_path( score_t sym_penalty) const
{
	int i,j;
	score_t max_score =NEG_INF;
	int best_left=-1, best_right=-1;

	for (i=0; i<cells.size()-1; i++)
		for (j=i+1; j<cells.size(); j++)
			if (cells[i][j].score > max_score)
			{
				best_left=i;
				best_right=j;
				max_score = cells[i][j].score;
			}

	MultiPath ret_path;
	if (max_score<0)
		return ret_path;

	// collect edges and create path
	vector<int> path_edges;
	int l=best_left, r=best_right;
	while (cells[l][r].prev_edge_idx>=0)
	{
		int e_idx = cells[l][r].prev_edge_idx;
		path_edges.push_back(e_idx);

		if (prm->multi_edges[e_idx].c_idx==r)
		{
			r = prm->multi_edges[e_idx].n_idx;
		}
		else
		{
			l = prm->multi_edges[e_idx].c_idx;
		}
	}
	
	prm->create_path_from_edges(path_edges,ret_path);
	ret_path.path_score = max_score;

	return ret_path;
}




struct edge_idx_set {
	edge_idx_set() : score(NEG_INF), length(0), num_aa(0) {};

	bool operator< (const edge_idx_set& other) const
	{
		return score > other.score;
	}

	edge_idx_set& operator= (const edge_idx_set& e)
	{
		score = e.score;
		length = e.length;
		memcpy(edge_idxs,e.edge_idxs,length*sizeof(int));

		return *this;
	}

	bool is_good_prefix(const vector<int>& good_idxs) const
	{
		
		int i;
		for (i=0; i<length; i++)
		{
			int j;
			for (j=0; j<good_idxs.size(); j++)
				if (edge_idxs[i]==good_idxs[j])
					break;
			if (j==good_idxs.size())
				return false;
		}
		return true;
	}

	void print()
	{
		cout << "P " << num_aa << "\t" << score << endl;
		int i;
		for (i=0; i<length; i++)
			cout << i << "\t" << edge_idxs[i] << endl;
		cout << endl;
	}

	score_t score;
	int length;
	int num_aa; // number of aas used to fill these edges
	int edge_idxs[32]; // max peptide size...
};




/****************************************************************************
   Computes for each node and number of amino acids, what is the maximal
   score attainable by making X steps forward from that node.
   Finds these values by performing a BFS search from all nodes.
   Cell i,j holds the vlaue for num steps i , node j.
*****************************************************************************/
void DeNovoDp::find_max_gains_per_length(int max_length, 
										 vector< vector< score_t > >& max_gains, 
										 bool verbose) const
{
	const int num_nodes = prm->get_num_nodes();
	const vector<Node>& nodes = prm->get_nodes();
	const vector<MultiEdge>& edges = prm->get_multi_edges();
	int n,i;

	max_gains.resize(num_nodes);
	for (n=0; n<num_nodes; n++)
	{
		max_gains[n].resize(max_length+1,NEG_INF);
		max_gains[n][0]=0;
	}


	for (i=1; i<=max_length; i++)
	{
		int n;
		for (n=0; n<num_nodes; n++)
		{
			const Node& node = nodes[n];
			const vector<int>& out_idxs = node.out_edge_idxs;
			int j;
			for (j=0; j<out_idxs.size(); j++)
			{
				const MultiEdge& edge = edges[out_idxs[j]];
				const int num_aa = edge.num_aa;
				const int c_idx  = edge.c_idx;
				if (num_aa>i)
					continue;

				const score_t new_score = max_gains[c_idx][i-num_aa] + nodes[c_idx].score + edge.max_variant_score;

				if (new_score > max_gains[n][i])
					max_gains[n][i]=new_score;
			}
		}
	}

	for (n=0; n<num_nodes; n++)
	{
		score_t max=NEG_INF;
		int max_idx = 0;
		for (i=0; i<=max_length; i++)
			if (max_gains[n][i]>max)
			{
				max=max_gains[n][i];
				max_idx=i;
			}

		for (i=max_idx +1; i<=max_length; i++)
			max_gains[n][i]=max;
	}



/*	max_gains.resize(max_length+1);
	for (i=0; i<max_gains.size(); i++)
		max_gains[i].resize(num_nodes,NEG_INF);

	for (n=0; n<num_nodes; n++)
		max_gains[0][n]=0;

	for (n=num_nodes-1; n>=0; n--)
	{
		if (nodes[n].score<=NEG_INF)
			continue;

		const Node& node = nodes[n];
		const vector<int>& out_idxs = node.out_edge_idxs;
		int j;


		for (j=0; j<out_idxs.size(); j++)
		{
			const MultiEdge& edge = edges[out_idxs[j]];

			// loop over different number of amino acids in edge
			const int num_aa = edge.num_aa;
			const int c_idx  = edge.c_idx;
			const score_t add_score = nodes[c_idx].score + edge.max_variant_score;

			int k;
			for (k=0; k<max_length-num_aa; k++)
				if (max_gains[k][c_idx]>NEG_INF && max_gains[k][c_idx]+add_score>max_gains[k+num_aa][n])
					max_gains[k+num_aa][n]=max_gains[k][c_idx]+add_score;
		}

		score_t max=NEG_INF;
		int max_idx = 0;
		for (i=1; i<=max_length; i++)
			if (max_gains[i][n]>max)
			{
				max=max_gains[i][n];
				max_idx=i;
			}

		for (i=max_idx +1; i<=max_length; i++)
			max_gains[i][n]=max;
	} */


	if (verbose)
	{
		cout << "MAX GAINS: " << endl;
		int n;
		for (n=0; n<num_nodes; n++)
		{
			cout << left << setw(4) << n << "  ";
			for (i=0; i<=max_length; i++)	
				cout << setw(5) << (max_gains[n][i]>-1000 ? max_gains[n][i] : -999) << " ";
			cout << endl;
		}
		cout << endl;
	}
}




/*******************************************************************************
 Returns all the top scoring paths, uses a DFS search with branch and bound pruning.
 limits the length of the solution to be between (approximately) supplied bounds.
************************************************************************************/
void DeNovoDp::get_top_scoring_antisymetric_paths_with_length_limits(
					vector<MultiPath>& multi_paths, 
					int required_num_paths, 
					int min_length,
					int max_length,
					score_t sym_penalty,
					score_t min_score_needed,
					bool try_complete_sequences,
					bool only_complete_sequences,
					double half_life_time) const
{
	const vector<Node>& nodes = prm->nodes;
	const vector<MultiEdge>& multi_edges = prm->multi_edges;
	const int num_nodes = nodes.size();
	const int last_node_idx = num_nodes-1;

	int last_heap_pos = required_num_paths - 1;
	int last_alt_heap_pos = required_num_paths - 1;


	const score_t non_complete_penalty = 0; // this is a penalty added for each end of the peptide
											  // that does not reach the terminal. Used to bias 										     // against sequences that do not read ends.
							
	vector<edge_idx_set> heap,		// heap is used to store the paths. If running time takes too long,
									// we start reducing its size so the  prunning becomes more efficient
									// in such a case, we start putting removed paths into the alt_heap
									// so we don't lose good paths
						 alt_heap;
	
	vector< vector<score_t> > max_gains_for_length; // the maximal score attainable from each node
										            // using a given number of amino acids.
												    // length, node_idx

	vector< int > node_ordering; // holds the optimal order of nodes for the BB search

	vector<score_t> added_scores;  // holds for each depth in the tree, the score that was
								   // added by using the edges in the current path          

	vector<int> out_idx_counters; // for each node, what branch are we going down

	vector<bool> used_nodes;  // indicators for each node if it was used in the current path


	vector<int> debug_idxs;
	bool debug_mode=false;
	if (debug_mode)
	{
		debug_idxs.push_back(127);
		debug_idxs.push_back(10);
		debug_idxs.push_back(38);
		debug_idxs.push_back(56);
	//	debug_idxs.push_back(576);
	//	debug_idxs.push_back(2950);
	//	debug_idxs.push_back(855);
	//	debug_idxs.push_back(860);
	}

	int i;
	if (max_length>31)
		max_length=31;
	
	added_scores.resize(num_nodes,NEG_INF);
	out_idx_counters.resize(num_nodes,0);
	used_nodes.resize(num_nodes,false);
	heap.clear();
	alt_heap.clear();
	heap.resize(required_num_paths);

	find_max_gains_per_length(max_length,max_gains_for_length);
	prm->sort_outgoing_edges_according_to_max_gains(max_gains_for_length);
	
	time_t start_t,end_t;
	start_t = time(NULL); 

	if (only_complete_sequences)
	{
		node_ordering.clear();
		node_ordering.push_back(0);
	}
	else
		prm->get_node_ordering_according_to_max_gains(max_gains_for_length, node_ordering);

	int ns;
	for (ns=0; ns<node_ordering.size(); ns++)
	{
		const int start_idx = node_ordering[ns];

		if (max_gains_for_length[start_idx][max_length]<heap[0].score)
			break;

		const int num_first_out_edges = nodes[start_idx].out_edge_idxs.size();
		edge_idx_set current_path;
		int          current_node_idx;
	
		current_node_idx=start_idx;
		current_path.score = nodes[start_idx].score;
		used_nodes[start_idx] = true;

		if (try_complete_sequences && start_idx>0)
			current_path.score += non_complete_penalty;

		while (1)
		{
			
			// check if the search is running too long, if so decrease the heap size
			// and send the excess paths to the alternate heap
			end_t = time(NULL);
			const double iteration_time = (end_t - start_t);
			
			if (iteration_time>=5.0)
				break;

			if (iteration_time>=half_life_time)
			{
				if (alt_heap.size()==0)
					alt_heap.resize(heap.size());

				// reduce heap only if large enough
				if (heap.size()>7)
				{
					int half_size = heap.size()/2;
					while (heap.size()>half_size)
					{
						pop_heap(heap.begin(),heap.end());
						const edge_idx_set& removed_path = heap[heap.size()-1];
						if (removed_path.score>min_score_needed && removed_path.score>alt_heap[0].score)
						{
							pop_heap(alt_heap.begin(),alt_heap.end());
							alt_heap[last_alt_heap_pos] = removed_path;
							push_heap(alt_heap.begin(),alt_heap.end());
						}
						heap.pop_back();
					}
					last_heap_pos = heap.size()-1;
				}
				start_t = end_t;
			}

			if (out_idx_counters[current_node_idx] >= nodes[current_node_idx].out_edge_idxs.size())
			{
				if (current_node_idx == start_idx)
					break; // we've returned all the way back, exhausted this tree
				// store path if necessary

				bool store_anyway = (nodes[current_node_idx].out_edge_idxs.size() == 0 &&
					current_path.num_aa >= min_length &&
					current_path.num_aa <= max_length &&
					current_path.score > heap[0].score);

				if (try_complete_sequences &&
					! store_anyway &&
					current_path.score > heap[0].score && 
					current_node_idx == last_node_idx  && 
					start_idx == 0)
					store_anyway=true;

				if (store_anyway)
				{
					if (only_complete_sequences)
					{
						if (nodes[current_node_idx].type == NODE_C_TERM)
						{
							pop_heap(heap.begin(),heap.end());

							const edge_idx_set& removed_path = heap[last_heap_pos];
							if (alt_heap.size()>0 && 
								removed_path.score>min_score_needed &&
								removed_path.score>alt_heap[0].score)
							{
								pop_heap(alt_heap.begin(),alt_heap.end());
								alt_heap[last_alt_heap_pos] = removed_path;
								push_heap(alt_heap.begin(),alt_heap.end());
							}

							heap[last_heap_pos] = current_path;
							push_heap(heap.begin(),heap.end());
							if (debug_mode && current_path.is_good_prefix(debug_idxs))
							{
								cout << "PUSH0 : ";
								current_path.print();
							}
						}
					}
					else if (! try_complete_sequences || nodes[current_node_idx].type == NODE_C_TERM)
					{
						pop_heap(heap.begin(),heap.end());
						
						const edge_idx_set& removed_path = heap[last_heap_pos];
						if (removed_path.score>min_score_needed &&
							alt_heap.size()>0 &&
							removed_path.score>alt_heap[0].score)
						{
							pop_heap(alt_heap.begin(),alt_heap.end());
							alt_heap[last_alt_heap_pos] = removed_path;
							push_heap(alt_heap.begin(),alt_heap.end());
						}

						heap[last_heap_pos] = current_path;
						push_heap(heap.begin(),heap.end());
						if (debug_mode && current_path.is_good_prefix(debug_idxs))
						{
							cout << "PUSH1 : ";
							current_path.print();
						}
					}
					else
					{
						const score_t score_with_penalty = current_path.score + non_complete_penalty;
						if (score_with_penalty > min_score_needed &&
							score_with_penalty > heap[0].score)
						{
							current_path.score += non_complete_penalty;
							pop_heap(heap.begin(),heap.end());
							const edge_idx_set& removed_path = heap[last_heap_pos];
							if (alt_heap.size()>0 && removed_path.score>alt_heap[0].score)
							{
								pop_heap(alt_heap.begin(),alt_heap.end());
								alt_heap[last_alt_heap_pos] = removed_path;
								push_heap(alt_heap.begin(),alt_heap.end());
							}

							heap[last_heap_pos] = current_path;
							push_heap(heap.begin(),heap.end());
							current_path.score -= non_complete_penalty;
							if (debug_mode && current_path.is_good_prefix(debug_idxs))
							{
								cout << "PUSH2 : ";
								current_path.print();
							}
						}
					}
				}

				// backtrack
				out_idx_counters[current_node_idx] =0;
				used_nodes[current_node_idx]=false;
				current_path.length--;

				const int& path_length     = current_path.length;
				const MultiEdge& back_edge = multi_edges[current_path.edge_idxs[path_length]];
				current_path.num_aa		 -= back_edge.num_aa;
				current_node_idx          = back_edge.n_idx;
				current_path.score       -= added_scores[path_length];
				continue;
			}


			// discard this path if we are using too many edges or 
			// the score will not be able to improve enough
			// since the edges are sorted according to the gain they can bring,
			// none of the rest can help so we skip the rest

			const int remaining_aas = max_length - current_path.num_aa;
			const score_t threshold_score = (min_score_needed > heap[0].score ?  min_score_needed : heap[0].score);
			const score_t maximal_achievable_score = (remaining_aas>=0 ? current_path.score + max_gains_for_length[current_node_idx][remaining_aas] : NEG_INF);
			if (current_path.num_aa > max_length || maximal_achievable_score<threshold_score)
			{
				out_idx_counters[current_node_idx] = nodes[current_node_idx].out_edge_idxs.size();
				continue; 
			}
		
			// advance on the edge
			const int edge_idx = nodes[current_node_idx].out_edge_idxs[out_idx_counters[current_node_idx]];
			const MultiEdge& e = multi_edges[edge_idx];

			out_idx_counters[current_node_idx]++;
			current_node_idx = e.c_idx;
			out_idx_counters[current_node_idx] =0;
			used_nodes[current_node_idx]=true;
			
			added_scores[current_path.length] = e.max_variant_score + nodes[e.c_idx].score;

			// check if forbidden pair is used..
			if (forbidden_idxs[current_node_idx]>=0 && used_nodes[forbidden_idxs[current_node_idx]])
				added_scores[current_path.length] -= sym_penalty; 

			current_path.edge_idxs[current_path.length] = edge_idx;
			current_path.num_aa += multi_edges[edge_idx].num_aa;
			current_path.score += added_scores[current_path.length];
			current_path.length++;

			if (debug_mode && current_path.is_good_prefix(debug_idxs))
			{
				cout << "GP: " << current_path.length << "\t" << current_path.score << endl;
			}
			
			// 
			score_t heap_score = heap[0].score;
			score_t added_score = added_scores[current_path.length-1];

			// check if the path should be stored at this stage, and if we can mark
		

			if (added_scores[current_path.length-1]>-20.0 && 
				nodes[current_node_idx].out_edge_idxs.size()>0 &&
				current_path.num_aa <= max_length &&
				current_path.num_aa >= min_length &&
				current_path.score> heap[0].score )
			{
				if (only_complete_sequences)
				{
					if (nodes[current_node_idx].type == NODE_C_TERM)
					{
						pop_heap(heap.begin(),heap.end());
						const edge_idx_set& removed_path = heap[last_heap_pos];
						if (alt_heap.size()>0 &&
							removed_path.score > min_score_needed &&
							removed_path.score > alt_heap[0].score)
						{
							pop_heap(alt_heap.begin(),alt_heap.end());
							alt_heap[last_alt_heap_pos] = removed_path;
							push_heap(alt_heap.begin(),alt_heap.end());
						}

						heap[last_heap_pos] = current_path;
						push_heap(heap.begin(),heap.end());
						if (debug_mode && current_path.is_good_prefix(debug_idxs))
						{
							cout << "PUSH3 : ";
							current_path.print();
						}
					}
				}
				else if (! try_complete_sequences || nodes[current_node_idx].type == NODE_C_TERM)
				{
					pop_heap(heap.begin(),heap.end());
					const edge_idx_set& removed_path = heap[last_heap_pos];
					if (alt_heap.size()>0 && 
						removed_path.score > min_score_needed &&
						removed_path.score>alt_heap[0].score)
					{
						pop_heap(alt_heap.begin(),alt_heap.end());
						alt_heap[last_alt_heap_pos] = removed_path;
						push_heap(alt_heap.begin(),alt_heap.end());
					}

					heap[last_heap_pos] = current_path;
					push_heap(heap.begin(),heap.end());
					if (debug_mode && current_path.is_good_prefix(debug_idxs))
					{
						cout << "PUSH4 : ";
						current_path.print();
					}
				}
				else
				{

					if (current_path.score + non_complete_penalty > heap[0].score)
					{
						current_path.score += non_complete_penalty;
						pop_heap(heap.begin(),heap.end());
						const edge_idx_set& removed_path = heap[last_heap_pos];
						if (alt_heap.size()>0 &&
							removed_path.score> min_score_needed &&
							removed_path.score>alt_heap[0].score)
						{
							pop_heap(alt_heap.begin(),alt_heap.end());
							alt_heap[last_alt_heap_pos] = removed_path;
							push_heap(alt_heap.begin(),alt_heap.end());
						}

						heap[last_heap_pos] = current_path;
						push_heap(heap.begin(),heap.end());
						current_path.score -= non_complete_penalty;
						if (debug_mode && current_path.is_good_prefix(debug_idxs))
						{
							cout << "PUSH5 : ";
							current_path.print();
						}
					}
				}
			}
		}

		used_nodes[start_idx] = false;
	}

	if (heap.size()==0)
		return;

	if (alt_heap.size()>0) 	// transfer all paths that are in current heap
	{
		int i;
		for (i=0; i<heap.size(); i++)
			if (heap[i].score>alt_heap[0].score)
			{
				pop_heap(alt_heap.begin(),alt_heap.end());
				alt_heap[last_alt_heap_pos]=heap[i];
				push_heap(alt_heap.begin(),alt_heap.end());
			}
	}

	
	// work on the alt_heap if necessary
	vector<edge_idx_set>& final_heap = (alt_heap.size()>0 ? alt_heap : heap);
	
	sort(final_heap.begin(),final_heap.end());

	while (final_heap.size()>0 && final_heap[final_heap.size()-1].score<-40)
		final_heap.pop_back();

	int actual_num_paths = required_num_paths;
	if (final_heap.size()<actual_num_paths)
		actual_num_paths = final_heap.size();

	multi_paths.resize(actual_num_paths);

	for (i=0; i<actual_num_paths; i++)
	{
		int j;
		vector<int> edge_idxs;
		edge_idxs.resize(final_heap[i].length);
		for (j=0; j<final_heap[i].length; j++)
			edge_idxs[j]=final_heap[i].edge_idxs[j];

		if (0 && debug_mode)
		{
			cout << i << " >>> \t";
			final_heap[i].print();
		}

		prm->create_path_from_edges(edge_idxs, multi_paths[i]);

		multi_paths[i].path_score = final_heap[i].score;
		multi_paths[i].original_rank = i;

		const vector<int>& node_idxs = multi_paths[i].node_idxs;
		const int max_idx = node_idxs.size()-1;
		int num_forbidden =0;
		for (j=0; j<max_idx; j++)
		{
			const int n_idx = node_idxs[j];
			if (forbidden_idxs[n_idx]<0)
				continue;
			int k;
			for (k=j+1; k<node_idxs.size(); k++)
				if (node_idxs[k]==forbidden_idxs[n_idx])
					break;
			if (k<node_idxs.size())
				num_forbidden++;
		}
		multi_paths[i].num_forbidden_nodes = num_forbidden;

		if (edge_idxs.size()==0)
		{
			cout << "Error: path with 0 edges!" << endl;
			exit(0);
		}

		if (debug_mode)
		{
			cout << "MP " << i << "\t";
			multi_paths[i].print(config);
		}
	}


	if (debug_mode)
	{

		exit(0);
	}
}



/*****************************************************************************
Same as the above function, but only returns paths that have good start/end idxs
******************************************************************************/
void  DeNovoDp::get_top_scoring_antisymetric_paths_with_specified_start_end_idxs(
		const vector<bool>& ind_allowed_start_end_idxs,
		vector<MultiPath>& multi_paths, 
		int required_num_paths, 
		int min_length, 
		int max_length, 
		score_t sym_penalty,
		score_t min_score_needed,
		bool try_complete_sequences,
		bool only_complete_sequences,
		double half_life_time) const
{
	const vector<Node>& nodes = prm->nodes;
	const vector<MultiEdge>& multi_edges = prm->multi_edges;
	const int num_nodes = nodes.size();
	
	int last_heap_pos = required_num_paths - 1;
	int last_alt_heap_pos = required_num_paths - 1;

	const score_t non_complete_penalty = 0; // this is a penalty added for each end of the peptide
											  // that does not reach the terminal. Used to bias 										     // against sequences that do not read ends.
							
	vector<edge_idx_set> heap,		// heap is used to store the paths. If running time takes too long,
									// we start reducing its size so the  prunning becomes more efficient
									// in such a case, we start putting removed paths into the alt_heap
									// so we don't lose good paths
						 alt_heap;
	
	vector< vector<score_t> > max_gains_for_length; // the maximal score attainable from each node
										            // using a given number of amino acids.
												    // length, node_idx

	vector< int > node_ordering; // holds the optimal order of nodes for the BB search

	vector<score_t> added_scores;  // holds for each depth in the tree, the score that was
								   // added by using the edges in the current path          

	vector<int> out_idx_counters; // for each node, what branch are we going down

	vector<bool> used_nodes;  // indicators for each node if it was used in the current path


	int i;
	if (max_length>31)
		max_length=31;
	
	added_scores.resize(num_nodes,NEG_INF);
	out_idx_counters.resize(num_nodes,0);
	used_nodes.resize(num_nodes,false);
	heap.clear();
	alt_heap.clear();
	heap.resize(required_num_paths);

	find_max_gains_per_length(max_length, max_gains_for_length);
	prm->sort_outgoing_edges_according_to_max_gains(max_gains_for_length);
	
	clock_t start_t,end_t;
	start_t = clock();

	if (only_complete_sequences)
	{
		node_ordering.clear();
		node_ordering.push_back(0);
	}
	else
		prm->get_node_ordering_according_to_max_gains(max_gains_for_length, node_ordering);

	int ns;
	for (ns=0; ns<node_ordering.size(); ns++)
	{
		const int start_idx = node_ordering[ns];

		if (! ind_allowed_start_end_idxs[start_idx])
			continue;

		if (max_gains_for_length[start_idx][max_length]<heap[0].score)
			break;

		const int num_first_out_edges = nodes[start_idx].out_edge_idxs.size();
		edge_idx_set current_path;
		int          current_node_idx;
	
		current_node_idx=start_idx;
		current_path.score = nodes[start_idx].score;
		used_nodes[start_idx] = true;

		if (try_complete_sequences && start_idx>0)
			current_path.score += non_complete_penalty;

		while (1)
		{
			
			// check if the search is running too long, if so decrease the heap size
			// and send the excess paths to the alternate heap
			end_t = clock();
			const double iteration_time = (end_t - start_t)/(double)CLOCKS_PER_SEC;
			if (0 && iteration_time>half_life_time)
			{
				if (alt_heap.size()==0)
					alt_heap.resize(heap.size());

				// reduce heap only if large enough
				if (heap.size()>7)
				{
					int half_size = heap.size()/2;
					while (heap.size()>half_size)
					{
						pop_heap(heap.begin(),heap.end());
						const edge_idx_set& removed_path = heap[heap.size()-1];
						if (removed_path.score>min_score_needed && removed_path.score>alt_heap[0].score)
						{
							pop_heap(alt_heap.begin(),alt_heap.end());
							alt_heap[last_alt_heap_pos] = removed_path;
							push_heap(alt_heap.begin(),alt_heap.end());
						}
						heap.pop_back();
					}
					last_heap_pos = heap.size()-1;
				//	cout << "Reduced heap to : " << heap.size() << endl;
				}
				start_t = end_t;
			}

			if (out_idx_counters[current_node_idx] >= nodes[current_node_idx].out_edge_idxs.size())
			{
				if (current_node_idx == start_idx)
					break; // we've returned all the way back, exhausted this tree
				// store path if necessary
				if (ind_allowed_start_end_idxs[current_node_idx] &&
					nodes[current_node_idx].out_edge_idxs.size() == 0 &&
					current_path.num_aa >= min_length &&
					current_path.num_aa <= max_length &&
					current_path.score > heap[0].score)
				{
					if (only_complete_sequences)
					{
						if (nodes[current_node_idx].type == NODE_C_TERM)
						{
							pop_heap(heap.begin(),heap.end());

							const edge_idx_set& removed_path = heap[last_heap_pos];
							if (alt_heap.size()>0 && 
								removed_path.score>min_score_needed &&
								removed_path.score>alt_heap[0].score)
							{
								pop_heap(alt_heap.begin(),alt_heap.end());
								alt_heap[last_alt_heap_pos] = removed_path;
								push_heap(alt_heap.begin(),alt_heap.end());
							}

							heap[last_heap_pos] = current_path;
							push_heap(heap.begin(),heap.end());
						}
					}
					else if (! try_complete_sequences || nodes[current_node_idx].type == NODE_C_TERM)
					{
						pop_heap(heap.begin(),heap.end());
						
						const edge_idx_set& removed_path = heap[last_heap_pos];
						if (removed_path.score>min_score_needed &&
							alt_heap.size()>0 &&
							removed_path.score>alt_heap[0].score)
						{
							pop_heap(alt_heap.begin(),alt_heap.end());
							alt_heap[last_alt_heap_pos] = removed_path;
							push_heap(alt_heap.begin(),alt_heap.end());
						}

						heap[last_heap_pos] = current_path;
						push_heap(heap.begin(),heap.end());
					}
					else
					{
						const score_t score_with_penalty = current_path.score + non_complete_penalty;
						if (score_with_penalty > min_score_needed &&
							score_with_penalty > heap[0].score)
						{
							current_path.score += non_complete_penalty;
							pop_heap(heap.begin(),heap.end());
							const edge_idx_set& removed_path = heap[last_heap_pos];
							if (alt_heap.size()>0 && removed_path.score>alt_heap[0].score)
							{
								pop_heap(alt_heap.begin(),alt_heap.end());
								alt_heap[last_alt_heap_pos] = removed_path;
								push_heap(alt_heap.begin(),alt_heap.end());
							}

							heap[last_heap_pos] = current_path;
							push_heap(heap.begin(),heap.end());
							current_path.score -= non_complete_penalty;
						}
					}
				}

				// backtrack
				out_idx_counters[current_node_idx] =0;
				used_nodes[current_node_idx]=false;
				current_path.length--;

				const int& path_length     = current_path.length;
				const MultiEdge& back_edge = multi_edges[current_path.edge_idxs[path_length]];
				current_path.num_aa		 -= back_edge.num_aa;
				current_node_idx          = back_edge.n_idx;
				current_path.score       -= added_scores[path_length];
				continue;
			}


			// discard this path if we are using too many edges or 
			// the score will not be able to improve enough
			// since the edges are sorted according to the gain they can bring,
			// none of the rest can help so we skip the rest

			const int remaining_aas = max_length - current_path.num_aa;
			const score_t threshold_score = (min_score_needed > heap[0].score ?  min_score_needed : heap[0].score);
			const score_t maximal_achievable_score = (remaining_aas>=0 ? current_path.score + max_gains_for_length[current_node_idx][remaining_aas] : NEG_INF);
			if (current_path.num_aa > max_length || maximal_achievable_score<threshold_score)
			{
				out_idx_counters[current_node_idx] = nodes[current_node_idx].out_edge_idxs.size();
				continue; 
			}
		
			// advance on the edge
			const int edge_idx = nodes[current_node_idx].out_edge_idxs[out_idx_counters[current_node_idx]];
			const MultiEdge& e = multi_edges[edge_idx];

			out_idx_counters[current_node_idx]++;
			current_node_idx = e.c_idx;
			out_idx_counters[current_node_idx] =0;
			used_nodes[current_node_idx]=true;
			
			added_scores[current_path.length] = e.max_variant_score + nodes[e.c_idx].score;

			// check if forbidden pair is used..
			if (forbidden_idxs[current_node_idx]>=0 && used_nodes[forbidden_idxs[current_node_idx]])
				added_scores[current_path.length] -= sym_penalty; 

			current_path.edge_idxs[current_path.length] = edge_idx;
			current_path.num_aa += multi_edges[edge_idx].num_aa;
			current_path.score += added_scores[current_path.length];
			current_path.length++;

	
			// 
			score_t heap_score = heap[0].score;
			score_t added_score = added_scores[current_path.length-1];

			// check if the path should be stored at this stage, and if we can mark
			if (ind_allowed_start_end_idxs[current_node_idx] &&
				added_scores[current_path.length-1]>-5 && 
				nodes[current_node_idx].out_edge_idxs.size()>0 &&
				current_path.num_aa <= max_length &&
				current_path.num_aa >= min_length &&
				current_path.score> heap[0].score )
			{
				if (only_complete_sequences)
				{
					if (nodes[current_node_idx].type == NODE_C_TERM)
					{
						pop_heap(heap.begin(),heap.end());
						const edge_idx_set& removed_path = heap[last_heap_pos];
						if (alt_heap.size()>0 &&
							removed_path.score > min_score_needed &&
							removed_path.score > alt_heap[0].score)
						{
							pop_heap(alt_heap.begin(),alt_heap.end());
							alt_heap[last_alt_heap_pos] = removed_path;
							push_heap(alt_heap.begin(),alt_heap.end());
						}

						heap[last_heap_pos] = current_path;
						push_heap(heap.begin(),heap.end());
					}
				}
				else if (! try_complete_sequences || nodes[current_node_idx].type == NODE_C_TERM)
				{
					pop_heap(heap.begin(),heap.end());
					const edge_idx_set& removed_path = heap[last_heap_pos];
					if (alt_heap.size()>0 && 
						removed_path.score > min_score_needed &&
						removed_path.score>alt_heap[0].score)
					{
						pop_heap(alt_heap.begin(),alt_heap.end());
						alt_heap[last_alt_heap_pos] = removed_path;
						push_heap(alt_heap.begin(),alt_heap.end());
					}

					heap[last_heap_pos] = current_path;
				}
				else
				{

					if (current_path.score + non_complete_penalty > heap[0].score)
					{
						current_path.score += non_complete_penalty;
						pop_heap(heap.begin(),heap.end());
						const edge_idx_set& removed_path = heap[last_heap_pos];
						if (alt_heap.size()>0 &&
							removed_path.score> min_score_needed &&
							removed_path.score>alt_heap[0].score)
						{
							pop_heap(alt_heap.begin(),alt_heap.end());
							alt_heap[last_alt_heap_pos] = removed_path;
							push_heap(alt_heap.begin(),alt_heap.end());
						}

						heap[last_heap_pos] = current_path;
						push_heap(heap.begin(),heap.end());
						current_path.score -= non_complete_penalty;
					}
				}
			}
		}

		used_nodes[start_idx] = false;
	}


	if (heap.size()==0)
		return;

	if (alt_heap.size()>0) 	// transfer all paths that are in current heap
	{
		int i;
		for (i=0; i<heap.size(); i++)
			if (heap[i].score>alt_heap[0].score)
			{
				pop_heap(alt_heap.begin(),alt_heap.end());
				alt_heap[last_alt_heap_pos]=heap[i];
				push_heap(alt_heap.begin(),alt_heap.end());
			}
	}

	
	// work on the alt_heap if necessary
	vector<edge_idx_set>& final_heap = (alt_heap.size()>0 ? alt_heap : heap);
	
	sort(final_heap.begin(),final_heap.end());

	while (final_heap.size()>0 && final_heap[final_heap.size()-1].score<-40)
		final_heap.pop_back();

	int actual_num_paths = required_num_paths;
	if (final_heap.size()<actual_num_paths)
		actual_num_paths = final_heap.size();

	multi_paths.resize(actual_num_paths);

	for (i=0; i<actual_num_paths; i++)
	{
		int j;
		vector<int> edge_idxs;
		edge_idxs.resize(final_heap[i].length);
		for (j=0; j<final_heap[i].length; j++)
			edge_idxs[j]=final_heap[i].edge_idxs[j];

		prm->create_path_from_edges(edge_idxs, multi_paths[i]);

		multi_paths[i].path_score = final_heap[i].score;
		multi_paths[i].original_rank = i;

		const vector<int>& node_idxs = multi_paths[i].node_idxs;
		const int max_idx = node_idxs.size()-1;
		int num_forbidden =0;
		for (j=0; j<max_idx; j++)
		{
			const int n_idx = node_idxs[j];
			if (forbidden_idxs[n_idx]<0)
				continue;
			int k;
			for (k=j+1; k<node_idxs.size(); k++)
				if (node_idxs[k]==forbidden_idxs[n_idx])
					break;
			if (k<node_idxs.size())
				num_forbidden++;
		}
		multi_paths[i].num_forbidden_nodes = num_forbidden;
	}
}











