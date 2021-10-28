#ifndef __DENOVODP_H__
#define __DENOVODP_H__

#include "PrmGraph.h"


void expand_path_multi_edges(const PrmGraph& prm, const vector<MultiPath>& multi_paths, 
							 vector<SeqPath>& paths, int max_num_paths = 20000);


// holds info for path between  nodes i,j in the graph
struct DPCell {
	DPCell() : score(NEG_INF), prev_edge_idx(-1), is_forbidden(0) {}

	score_t score; // the maximal socre attainable between the nodes i,j
	int prev_edge_idx; // the edge idx used to get to this cell
	int is_forbidden; // indicator if this cell is a forbidden pair
};


// This class holds the dynamic programming table on which de novo
// searches and tag generation is done.
class DeNovoDp {
public:
	DeNovoDp() : prm(NULL), config(NULL) {}

	MultiPath get_top_scoring_antisymetric_path(score_t sym_penalty = 15) const;

	void get_top_scoring_antisymetric_paths_with_length_limits(
		vector<MultiPath>& multi_paths, 
		int num_paths, 
		int min_length, 
		int max_length, 
		score_t sym_penalty = 15,
		score_t min_score_needed = NEG_INF,
		bool try_complete_sequences = false,
		bool only_complete_sequences = false,
		double half_life_time = 2.0) const;

	void get_top_scoring_antisymetric_paths_with_specified_start_end_idxs(
		const vector<bool>& ind_allowed_start_end_idxs,
		vector<MultiPath>& multi_paths, 
		int num_paths, 
		int min_length, 
		int max_length, 
		score_t sym_penalty = 15,
		score_t min_score_needed = NEG_INF,
		bool try_complete_sequences = false,
		bool only_complete_sequences = false,
		double half_life_time = 2.0) const;


	// fills in all cells in the dp_table according to the PrmGraph
	void fill_dp_table(const PrmGraph *_prm, score_t sym_penalty);

private:
	vector< vector<DPCell> > cells; 
	vector<int> forbidden_idxs; // for each node holds the idx of its closet forbidden pair
							   // (-1  if it doesn't have one)
	PrmGraph *prm;
	Config *config;

	// functions

	void find_max_gains_per_length(int max_length, 
		vector< vector< score_t > >&  max_gains, bool verbose = false) const;

};


#endif


