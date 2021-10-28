#ifndef __DENOVOSOLUTIONS_H__
#define __DENOVOSOLUTIONS_H__


#include "Spectrum.h"
#include "DeNovoDp.h"
#include "AllScoreModels.h"
#include <set>


struct SeqPathKey {
	SeqPathKey() : pep_str(""), n_mass(NEG_INF), sort_key(NEG_INF), path_pos_idx(NEG_INF) {}
	SeqPathKey(const SeqPath& sp) : pep_str(sp.seq_str), n_mass(sp.n_term_mass),
									sort_key(sp.sort_key){}

	bool operator < (const SeqPathKey& other) const
	{
		return (pep_str < other.pep_str ||
			(pep_str == other.pep_str && other.n_mass - n_mass > tolerance));
	}

	static void set_tolerance(mass_t t) { tolerance = t; }
	static mass_t tolerance;
	string pep_str;
	mass_t n_mass;
	float  sort_key;
	int	   path_pos_idx;
};

struct idx_score_pair {
	idx_score_pair() : path_pos_idx(NEG_INF), sort_key((float)NEG_INF) {}
	idx_score_pair(const SeqPathKey& key) : path_pos_idx(key.path_pos_idx),
		sort_key(key.sort_key) {}

	bool operator< (const idx_score_pair& other) const
	{
		return sort_key>other.sort_key;
	}

	int path_pos_idx;
	float sort_key;
};


// well it's not really a heap
struct SeqPathHeap {
	void init(int _max_size, mass_t _tolerance) {
		paths.clear(); 
		paths.reserve(_max_size);
		min_score_heap.clear();
		path_keys.clear();

		max_size = _max_size; 
		tolerance = 8 * _tolerance;
		min_idx=-1; 
		min_value=99999;
		SeqPathKey::set_tolerance(tolerance);
	}

	int add_path(SeqPath& new_path, bool verbose = false);


	int max_size;
	mass_t tolerance;
	int min_idx;
	float min_value;
	vector<SeqPath> paths;
	set<SeqPathKey> path_keys;
	vector<idx_score_pair> min_score_heap;
};



/**************************************************************************
	Wrapper funciton that generates the desired solutions.
	Combines both local and global solutions (similar to the PepNovoTag
	and LocalTag solutions).
***************************************************************************/
bool generate_denovo_solutions(PrmGraph * & prm,
							   AllScoreModels *model, 
							   Spectrum *spec,
							   bool denovo_mode,
							   mass_t pm_with_19,
							   int  charge,
							   int num_sols, 
							   int min_length, 
							   int max_length,
							   score_t min_score_needed,
							   vector<SeqPath>& solutions,
							   bool only_complete = false,
							   bool only_from_graph_containing_true_pep = false,
							   bool need_to_create_PrmGraph = true);






/***************************************************************************
	Wrapper function that generates several solutions according to different 
	precursor masses.
****************************************************************************/
void generate_denovo_solutions_from_several_pms(
							   vector<PrmGraph *>& prm_ptrs,
							   AllScoreModels *model, 
							   Spectrum *spec,
							   bool denovo_mode,
							   int num_sols, 
							   int min_length, 
							   int max_length,
							   vector<mass_t>& different_pms_with_19,
							   vector<int>& charges,
							   vector<SeqPath>& solutions,
							   bool ind_only_complete = false);

void generate_denovo_solutions_from_several_pms_with_good_start_end_idxs(
							   vector<PrmGraph *>& prm_ptrs,
							   AllScoreModels *model,
							   Spectrum *spec,
							   bool denovo_mode,
							   int num_sols, int min_length, int max_length,
							   vector<mass_t>& different_pms_with_19,
							   vector<int>& charges,
							   vector<SeqPath>& solutions);


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
							   vector<SeqPath>& solutions);

/**************************************************************************
Generates an prints the PRM graph
***************************************************************************/
void print_prm_graph_scores(AllScoreModels *model, Spectrum *spec, 
					 mass_t pm_with_19, int charge, bool prm_norm);







/*************************************************************************
Generates tags by making a mixture of local/de novo tags
checks which tags appear in the longer denovo sequences sets the boolean
indicators in the tag seq paths
**************************************************************************/
void generate_tags(vector<PrmGraph *>& prm_ptrs,
				   AllScoreModels *model,
				   AnnotatedSpectrum& as,
				   const vector<int>& max_num_tags,
				   int main_tag_length,				 // the length for which we parse de novo sequences
				   const vector<mass_t>& pms_with_19,
				   const vector<int>& charges,
				   vector<SeqPath>& final_tags,
				   bool use_original_num_tags=false,
				   int  prm_ptr_start_idx=0);



//void output_denovo_solutions(SingleSpectrumFile *ssf, const Config *config, ostream& out_stream,
//							 const vector<SeqPath>& solutions, int max_num_sols = -1);


//void output_tag_solutions(SingleSpectrumFile *ssf, const Config *config, ostream& out_stream,
//							 const vector<SeqPath>& solutions, bool	output_aa_probs);

void perform_denovo_on_list_of_files(AllScoreModels& model, 
									 const vector<string>& list_vector, 
									 int file_start_idx,
									 int num_solutions, 
									 int min_length, 
									 int max_length,
									 bool report_progress,
									 float min_filter_prob,
									 bool  output_aa_probs,
									 bool  output_cumulative_seq_probs,
									 ostream& out_stream = cout);


void perform_tags_on_list_of_files(AllScoreModels& model, 
									 const vector<string>& list_vector, 
									 int file_start_idx,
									 int num_solutions, 
									 int	tag_length,
									 bool report_progress,
									 float min_filter_prob,
									 bool	output_aa_probs,
									 bool	output_cumulative_seq_probs,
									 ostream& out_stream = cout);


void new_perform_denovo_on_list_of_files(AllScoreModels& model, 
									 const vector<string>& list_vector, 
									 int file_start_idx,
									 int num_solutions, 
									 int min_length, 
									 int max_length,
									 bool report_progress,
									 float min_filter_prob,
									 bool  output_aa_probs,
									 bool  output_cumulative_seq_probs,
									 ostream& out_stream = cout);


void create_tag_file_for_inspect(AllScoreModels& model,
								 char *spectrum_file,
								 char *tag_string,
								 char *tag_suffix);

// makes a FASTA file with the sequences of full denovo sequences (completed
// from the SEQ in the annotated mgf file)
void make_denovo_training_fa(AllScoreModels& model,
								 char *mgf);

void benchmark_tags(AllScoreModels& model,
					char *spectrum_file,
					char *tag_string,
					int num_test_cases=-1);
								 

void perform_prm_on_list_of_files(AllScoreModels& model, 
								  const vector<string>& list_vector,
								  float sqs_filter_prob,
								  int file_start_idx,
								  bool prm_norm);

void prm_benchmark(AllScoreModels& model, 
				   const vector<string>& list_vector,
				   float sqs_filter_prob,
				   int file_start_idx);






#endif



