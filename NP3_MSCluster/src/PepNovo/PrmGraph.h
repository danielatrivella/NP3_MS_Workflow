#ifndef __PRMGRAPH_H__
#define __PRMGRAPH_H__

#include "Spectrum.h"
#include "AllScoreModels.h"


class PrmGraph;


struct SeqPath;

// This structure contains a path that can be expanded to many variants
struct MultiPath {

	MultiPath() : original_rank(-1), n_term_mass(-1), c_term_mass(-1),
		n_term_aa(N_TERM), c_term_aa(C_TERM), ind_compete_path(0), num_forbidden_nodes(0), path_score(0) {
		edge_idxs.clear(); node_idxs.clear(); breakages.clear(); }

	MultiPath(const struct Multipath& other);

	bool check_if_correct(const Peptide& p, Config *config) const;

	int  get_num_correct_aas(const PrmGraph& prm, const Peptide& p, Config *config) const;

	int  get_num_aas() const;

	void print(Config *config, ostream& os=cout) const;

	int original_rank;
	mass_t n_term_mass;
	mass_t c_term_mass;

	int n_term_aa;
	int c_term_aa;

	int ind_compete_path;   // set to 1 if the path is from the N-terminal to the C-terminal

	int num_forbidden_nodes;

	score_t  path_score;

	vector<int> edge_idxs; // for the multi edges
	vector<int> node_idxs;
	vector<Breakage *> breakages;
};





// The structure for a single path corresponsing to one amino acid sequence
struct SeqPath {

	SeqPath() : multi_path_rank(POS_INF), n_term_mass(-1), c_term_mass(-1),
		n_term_aa(N_TERM), c_term_aa(C_TERM), pm_with_19(-1), charge(0),
		num_forbidden_nodes (0), path_score(NEG_INF), multi_path_score(NEG_INF), 
		org_rank(NEG_INF), delta_score(NEG_INF),
		rerank_score(NEG_INF), cumulative_seq_prob(0), sort_key(0), prm_ptr(NULL), 
		tag_percent_top_5(0), tag_percent_top_20(0), tag_percent_all(0), 
		is_correct(false), is_cut_correct(false) { positions.clear(); };

	void print(ostream &os = cout) const;

	void print_with_probs(ostream &os = cout) const;

	void print_full(const Config *config, ostream &os = cout) const;

	void get_amino_acids(vector<int>& aas) const
	{
		int i;
		const int num_aa = positions.size()-1;

		if (num_aa<=0)
		{
			aas.clear();
			return;
		}

		aas.resize(num_aa);
		for (i=0; i<num_aa; i++)
			aas[i]=positions[i].aa;
	}

	void convert_QK()
	{
		int i;
		for (i=0; i<positions.size()-2; i++)
			if (positions[i].aa == Lys)
				positions[i].aa = Gln;

		if (c_term_mass<pm_with_19-30 && positions[i].aa == Lys)
			positions[i].aa = Gln;
	}


	// gives a penalty of digest score * # tolerance differences if the
	// peptide is complete and the mass is different form the one
	// used for the prm
	score_t adjust_complete_sequence_penalty(PrmGraph *prm);

	// calculates score according to the aa probabilties
	// gives premium if the path belonged to a denovo solution 
	void calc_confidence_score(int denovo_rank=-1);

	int get_num_correct_aas(const Peptide& pep, const Config *config) const;

	mass_t calculate_peptide_mass(const Config *config) const;

	
	bool check_if_correct(const string& str, const Config *config) const;

	bool check_if_cut_correct(const vector<mass_t>& exp_cut_masses, mass_t tolerance) const;

	// returns the number of b/y ions
	int get_num_frags(const vector<int>& frag_idxs) const;


	void parse_path_to_smaller_ones(const Config *config,
		int min_length, int max_length, vector<SeqPath>& new_paths) const;


	int get_num_aa() const { return positions.size()-1; }


	// adds the relevant PathPos to the path and adjusts the other non-terminal values 
	void add_edge_variant(const MultiEdge& edge, int e_idx, int varaint_idx);

	void make_seq_str(const Config *config);

	vector<PathPos> positions; // there are n+1 positions for n amino acids

	int multi_path_rank;
		
	mass_t n_term_mass;
	mass_t c_term_mass;

	int n_term_aa;
	int c_term_aa;

	mass_t pm_with_19;
	int    charge;

	int num_forbidden_nodes;

	score_t path_score;

	score_t multi_path_score;

	int		org_rank;

	float   delta_score;

	float	rerank_score;
	
	float	cumulative_seq_prob;

	float	sort_key; // according to the seq heap is sorted, higher value gets first position in heap

	PrmGraph *prm_ptr;

	string seq_str;

	float	tag_percent_top_5;
	float	tag_percent_top_20;
	float   tag_percent_all;

	bool    is_correct;
	bool	is_cut_correct;
};

bool comp_SeqPath_sort_key(const SeqPath& a, const SeqPath& b);
bool comp_SeqPath_path_score (const SeqPath& a, const SeqPath& b);



class PrmGraph {
	friend class DeNovoDp;
	friend class EdgeModel;
public:

	void create_graph_from_spectrum(AllScoreModels *model, Spectrum *spectrum, 
					mass_t pm_with_19, int spec_charge=0, bool add_all_pepitde_nodes=false,
					bool only_basic_score = false);

	void create_graph_for_peptide_and_spectrum(AllScoreModels *model, Spectrum *spectrum, 
					mass_t pm_with_19, int spec_charge, const Peptide& peptide);

	void extract_breakage_infos_for_score_training(AllScoreModels *model,
												   int frag_idx,
												   int target_region_idx,
												   bool ind_strong_frag,
												   vector<BreakageInfo>& good_examples,
												   vector<BreakageInfo>& bad_examples) const;


	// removes all edges to and from nodes with the active flag set to 0
	void remove_edges_from_inactive_nodes();

	// finds the highest scoring continuous subpath for a given peptide in the graph
	SeqPath get_highest_scoring_subpath(const Peptide& peptide, mass_t start_mass =0) const;

	SeqPath get_longest_subpath(const Peptide& peptide, mass_t start_mass,
						bool verbose = false);

	SeqPath get_path_from_peptide_prm_graph(const Peptide& peptide) const;


	int get_num_nodes() const { return nodes.size(); }
	Spectrum *get_source_spectrum() const { return source_spectrum; }
	Config   *get_config() const { return config; }
	AllScoreModels* get_model()  const { return model;  }

	// computes amino acid probabilities according to the AminoAcidProbs models
	// if no model exists for the seq_path, aa probabilties of -1 wil be assigned
	// and false will be returned
	bool calc_amino_acid_probs(SeqPath& seq_path, int seq_rank = 5);


	
	void expand_multi_path(const MultiPath& multi_path, 
		                   vector<SeqPath>& seq_paths,
						   score_t min_score,
						   score_t forbidden_pair_penalty,
						   int max_num_paths = -1) const;

	void fast_expand_multi_path(const MultiPath& multi_path, 
		                   vector<SeqPath>& seq_paths,
						   score_t min_score,
						   score_t forbidden_pair_penalty,
						   int max_num_paths = -1) const;

	void fast_expand_multi_path_for_combo_scores(
						   AllScoreModels *model,
						   const MultiPath& multi_path, 
		                   vector<SeqPath>& seq_paths,
						   score_t min_score,
						   score_t forbidden_pair_penalty,
						   int max_num_paths = -1);

	void parse_seq_path_to_smaller_ones(const SeqPath& org_path, 
										  int min_length, 
										  int max_length, 
										  vector<SeqPath>& new_paths);

	void expand_all_multi_paths(AllScoreModels *model, const vector<MultiPath>& multi_paths, 
			vector<SeqPath>& paths, score_t forbidden_pair_penalty = 25, int max_num_paths = 100);

	void get_all_correct_node_idxs(const Peptide& peptide, vector<int>& idxs) const;

	void get_all_mirror_node_idxs(const Peptide& peptide, vector<int>& idxs) const;

	// returns the idxs of nodes correponding to the expected breakages of the peptide
	void get_relevant_node_idxs(const Peptide& peptide, vector<int>& idxs) const;

	void fill_breakage_info(const AllScoreModels *model, BreakageInfo *info, int node_idx,
								  int n_edge_idx, int n_variant_idx,
								  int c_edge_idx, int c_variant_idx, int type =0) const;


	const MultiEdge& get_multi_edge(int e_idx) const { return multi_edges[e_idx]; }
	const Node& get_node(int n_idx) const { return nodes[n_idx]; }
	Node& get_non_const_node(int n_idx) { return nodes[n_idx]; }
	const vector<MultiEdge>& get_multi_edges() const { return multi_edges; }
	const vector<Node>& get_nodes() const { return nodes; }
	const vector<int>& get_forbidden_node_idxs() const { return forbidden_node_idxs; }
	void set_frobidden_node_idxs(const vector<int>& idxs) { forbidden_node_idxs = idxs; }
	const vector<score_t>& get_cummulative_scores() const { return cummulative_scores; }

	mass_t get_max_node_mass() const { return nodes[nodes.size()-1].mass; }


	score_t    get_max_node_score() const { return max_node_score; }

	bool get_has_node_combo_scores() const { return has_node_combo_scores; }
	void set_has_node_combo_scores(bool b) { has_node_combo_scores = b; }

	// for quick access to nodes array
	PeakRange get_nodes_in_range(mass_t min_r, mass_t max_r) const;
	int get_max_score_node(mass_t exp_mass, mass_t tolerance) const;
	int get_min_dis_node(mass_t mass, mass_t tolerance) const;

	void print(ostream& os = cout, bool print_edge_scores = false) const;

	void print_only_scores() const;

	void clear();

	void print_with_multi_edges() const;

	void print_with_combo_tables() const;

	void print_edge_label(int edge_idx, int var_idx, ostream& os = cout) const;

	int	get_charge() const { return charge; }
	int get_size_idx() const { return size_idx; }
	mass_t get_pm_with_19() const { return pm_with_19; }


	
	

private:
	vector<Node>	  nodes;
	vector<MultiEdge> multi_edges;

	vector<int>  nodeIndexArray_;

	Spectrum * source_spectrum;
	Config   * config;
	AllScoreModels    * model;

	mass_t   pm_with_19; // the pm_with_19 used to create the graph
	int      size_idx;   // 
	int      charge;

	mass_t min_significant_mass;  // the mass of the highest scoring node in the below 10% mass range
	mass_t max_significant_mass;  // the mass of the highest scoring node in the above 90% mass range

	int    min_significant_idx, max_significant_idx;

	score_t max_node_score;

	score_t digest_node_score;

	bool has_node_combo_scores;

	vector<int>     forbidden_node_idxs;
	vector<score_t> cummulative_scores;

	vector<int> n_digest_node_idxs, c_digest_node_idxs;
	
	vector< vector< bool > > out_aa_ind, in_aa_ind; // holds for each amino acid and node
												 // whether there is a single aa leaving that node


	// functions
	void init_index_array();
	void merge_close_nodes();
	void create_nodes();
	void add_digest_nodes(mass_t n_digest_mass=NEG_INF, mass_t c_digest_mass=NEG_INF);
	void score_nodes(AllScoreModels *model);

	void create_nodes_for_peptide(const Peptide& pep);
	

	// this function performs all the scoring operations on edges 
	// (amino acid scores, missing cleavage scores etc.)
	score_t calc_edge_variant_score(const MultiEdge& egde, int num_aa, int *aa) const;

	// add the variant ptr and scores for this combo
	void add_and_score_edge_variants(const AA_combo& combo, MultiEdge& edge);


	// fill edges
	void fill_single_multi_edges();
	
	void fill_double_multi_edges(bool add_overlap_edges=false);

	void fill_longer_multi_edges(int max_edge_length, bool add_overlap_edges=false);

	// assigns a value to each node's rank field
	void rank_nodes_according_to_score();

	void prune_low_scoring_nodes();

	void set_idxs_max_in_out_for_nodes();

	void fill_forbidden_idx_indicators_and_cumulative_scores();


	// finds the edge idx of the edge from node i to j, returns -1 if no such node exists
	int find_edge_idx_ij(int i_idx, int j_idx) const
	{
		const vector<int>& out_idxs = nodes[i_idx].out_edge_idxs;
		int e;
		for (e=out_idxs.size()-1; e>=0; e--)
			if (multi_edges[out_idxs[e]].c_idx == j_idx)
				return out_idxs[e];
		return -1;
	}



	// sorts edges according to the value to which they lead
	void sort_outgoing_edges();

	void sort_outgoing_edges_according_to_max_gains(const vector< vector< score_t > >& max_gains);

	void get_node_ordering_according_to_max_gains(
	 vector< vector< score_t > >& max_gains_for_length,
	 vector<int>& node_order) const;


	// creates a path object from a collection of edges that are assumed
	// to correspond to a path in the graph
	void create_path_from_edges(vector<int>& edge_idxs, MultiPath& path) const;


	void print_multi_edges(int node_idx, bool print_edge_scores = false) const;
};


// more or less the same functions as the map with the spectrum
// why not use template? good question...
// this seems easier...

/*****************************************************************
Returns the indices of all peaks that are within the mass range
******************************************************************/
inline PeakRange  PrmGraph::get_nodes_in_range(mass_t min_r, mass_t max_r) const
{
	PeakRange pr;

	int max_node_idx = nodes.size()-1;

	if (max_r>= pm_with_19-3)
		max_r = pm_with_19-3;

	if (min_r<0)
		min_r=0;

	if (max_r<min_r)
		return pr;

	if (min_r>1)
		min_r--;

	int i_min=nodeIndexArray_[(int)(min_r)];
	int i_max=nodeIndexArray_[(int)(max_r+1.0)];
		
	if (i_max<max_node_idx)
		i_max++;

	while (nodes[i_min].mass<min_r && i_min<i_max)
		i_min++;

	while (nodes[i_max].mass>max_r && i_max>0)
		i_max--;

	if (nodes[i_min].mass > max_r || nodes[i_max].mass<min_r)
		return pr;

	pr.num_peaks = i_max - i_min+1;
	pr.low_idx =   i_min;
	pr.high_idx =  i_max;

	return pr;
}


/***********************************************************
returns the idx of the peak that is closes to the mass
(-1 is returned if no peak is found)
************************************************************/
inline int PrmGraph::get_min_dis_node(mass_t mass, mass_t tolerance) const
{
	PeakRange pr = get_nodes_in_range(mass-tolerance,mass+tolerance);
	if (pr.num_peaks==0)
		return -1;
	mass_t min_dis=fabs(mass - nodes[pr.low_idx].mass);
	int min_idx=pr.low_idx;
	int j;

	// find closest peak to exp_pos
	for (j=1; j<pr.num_peaks; j++)
	{
		mass_t dis = fabs(mass - nodes[pr.low_idx+j].mass);
		if (dis<min_dis)
		{
			min_dis=dis;
			min_idx = pr.low_idx+j;
		}
	}
	return min_idx;
}

/***********************************************************
returns the idx of the peak that has the highest intensity
in the range (-1 is returned if no peak is found)
Balances between intensity and proximity to the expected position
A peak at the edge needs to be at least 2 times stronger than
a peak exactly at the middle to be chosen.
************************************************************/
inline int PrmGraph::get_max_score_node(mass_t exp_mass, mass_t tolerance) const
{
	PeakRange pr = get_nodes_in_range(exp_mass-tolerance,exp_mass+tolerance);
	if (pr.num_peaks==0)
		return -1;

	if (pr.num_peaks ==1)
		return pr.low_idx;

	mass_t max_val = (2*tolerance - fabs(nodes[pr.low_idx].mass - exp_mass)) *
								nodes[pr.low_idx].score;

	int max_idx=pr.low_idx;
	int j;
	// find closest peak to exp_pos
	for (j=1; j<pr.num_peaks; j++)
	{
		mass_t peak_val = (2*tolerance - fabs(nodes[pr.low_idx+j].mass - exp_mass)) *
								nodes[pr.low_idx+j].score;
		if (peak_val>max_val)
		{
			max_val=peak_val;
			max_idx = pr.low_idx+j;
		}
		else
			break;
	}
	return max_idx;
}




#endif

