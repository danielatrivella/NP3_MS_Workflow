#ifndef __BASICDATASTRUCTS_H__
#define __BASICDATASTRUCTS_H__

#include "Config.h"

/// \struct Peak
/// Holds the basic information for a single peak as read from a spectrum file
struct Peak {
	Peak() : mass(-1.0),  intensity(0.0), count(1), maxPossible(1), charge(0) {};
	
	bool operator< ( const Peak& rhs) const
	{
		return (mass<rhs.mass);
	}

	void print(ostream& os = cout) const
	{
		os << mass << "\t" << intensity << "\t" << static_cast<int>(count) << "\t" 
			<< static_cast<int>(maxPossible) << "\t" << static_cast<int>(charge) << endl;
	}

	mass_t         mass;         /// in Daltons
	intensity_t    intensity;
	unsigned char count;        // used in clustering (multiplicity observed)
	unsigned char maxPossible;  // used in clustering (maximal multiplicity possible)
	unsigned short charge;      // useful for high resolution, currently ignored
};

/// \struct IndexIntensityPair
/// used for sorting peaks with sort()
/// The peak index with the strongest intensity will be in position 0
struct IndexIntensityPair {
	IndexIntensityPair() : index(-1), intensity(NEG_INF) {};
	IndexIntensityPair(int i, intensity_t f) : index(i), intensity(f) {};

	bool operator< (const IndexIntensityPair& rhs) const
	{
		return (intensity > rhs.intensity);
	}

	int index;
	intensity_t intensity;
};


class Peptide {
public:
	Peptide() : mass(0), n_gap(0), aa_before(Gap), aa_after(Gap), peptideStr_("") {};

	void clear() { mass=0; amino_acids.clear(); n_gap=0; aa_before=Gap; aa_after=Gap;}

	bool operator== (const Peptide& other) const
	{
		if (amino_acids.size() != other.amino_acids.size())
			return false;
		int i;
		for (i=0; i<amino_acids.size(); i++)
			if (amino_acids[i] != other.amino_acids[i])
				return false;
		return true;
	}

	mass_t get_mass() const { return mass; } // mass of aas + terminals, doesn't include 19
	mass_t get_mass_with_19() const { return (n_gap + mass + MASS_OHHH); }

	int calc_charge(mass_t& m_over_z) const;

	int get_length() const { return amino_acids.size(); } // gaps count as amino acids
	int get_num_aas() const { return amino_acids.size(); }

	const string& getPeptideStr() const { return peptideStr_; }
	void setPeptideStr(const string& str) { peptideStr_ = str; }

	const vector<int>& get_amino_acids() const { return amino_acids; }

	void generate_random_peptide(const Config *config, int peptide_length);

	/*********************************************************************
	Returns the global edit distance between two peptides.
	**********************************************************************/
	float peptide_edit_distance(Config *config, Peptide& other) const;

	// changes the amino acids I->L
	// an Q->K if not at terminal and tolerance > 0.1
	void convert_ILQK(const Config *config);

	// changes all the amino acids to their original form (without PTMs)
	void convert_to_org(const Config *config);

	void convert_IL();

	void reverse();

	string as_string(const Config* config) const;

	bool parseFromString(const Config* config, const string& str);

	void set_peptide(vector<int>& aas, mass_t mass, mass_t n_gap,
					 int n_term_aa = N_TERM, int c_term_aa = C_TERM);

	void set_peptide_aas(const vector<int>& aas) { amino_acids = aas; }

	void set_aa_before(int a) { aa_before=a; }
	void set_aa_after(int a)  { aa_after=a; }
	int get_aa_before() const { return aa_before; }
	int get_aa_after() const  { return aa_after; }

	mass_t get_n_gap() const { return n_gap; }
	void   set_n_gap(mass_t gap) { n_gap = gap; }

	void calc_mass(const Config *config);

	void calc_expected_breakage_masses(const Config *config, vector<mass_t>& break_masses) const;

	int calc_number_of_correct_aas(const Config *config, const Peptide& other) const;


private:

	mass_t mass; // mass of aas + terminals, doesn't include 19

	mass_t n_gap; // if partial sequence

	int aa_before, aa_after; // before and after cleavage

	string peptideStr_;

	vector<int>    amino_acids;
	
};


/*****************************************************************************
Used to represent a peptide when predicting peaks or rank socring.
******************************************************************************/
struct PeptideSolution {
	PeptideSolution() : num_correct_aas(0), charge(0), type(-1), weight(0), MQScore(NEG_INF),
		reaches_n_terminal(true), reaches_c_terminal(true),
		pm_with_19(NEG_INF), most_basic_aa_removed_from_n(0), most_basic_aa_removed_from_c(0) {}

	int	   num_correct_aas;
	int	   charge;
	int	   type; // 0 correct , 1 db , 2 full de novo, 3 cross db

	float  weight;
	float  MQScore; 

	bool reaches_n_terminal;
	bool reaches_c_terminal;

	mass_t pm_with_19; // spectrum value / corrected / peptide value (which ever is available)
					   // if possible, use the peptide value, otherwise use the corrected value,
					   // if none are avialble then this is the spectrum's pm_with_19 (least accurate)
	int most_basic_aa_removed_from_n;
	int most_basic_aa_removed_from_c;

	Peptide pep;

	int  calc_mobility() const;
};



struct PeakRange {
	PeakRange() : num_peaks(0), low_idx(-1), high_idx(-1) {};
	int num_peaks;
	int low_idx;
	int high_idx;
};
typedef struct PeakRange PeakRange;




struct BreakageFragment {
	BreakageFragment() : frag_type_idx(-1), peak_idx(-1), peak_level(-1), 
						 is_strong_fragment(0), mass(-1), expected_mass(-1), intensity(0) {};

	int frag_type_idx;    // idx in theFragmenTypeSet
	int peak_idx;

	int peak_level;          // used in models that have discrete rank/intensity levels
	int is_strong_fragment;  // 1 if this is a strong fragment according to the region
							 // in which the breakage occurs, 0 otherwise

	mass_t mass;             // the mass of the peak in the spectrum, -1 if no peak
	mass_t expected_mass; 
	intensity_t intensity;
};


struct Breakage {
	Breakage() : region_idx(-1), mass(-1), 
			     total_intensity(0),  num_frags_detected(0),
				 parent_charge(0),    parent_size_idx(0), score(NEG_INF) {};


	// returns the position in the fragments vector of the given frag_type_idx
	// returns -1 if this frag_type_idx is not found
	int get_position_of_frag_idx(int frag_type_idx) const
	{
		int i;
		for (i=0; i<fragments.size(); i++)
			if (fragments[i].frag_type_idx == frag_type_idx)
				return i;
		return -1;
	}


	bool are_all_frag_types_visible() const
	{
		return (frag_type_idxs_not_visible.size() == 0);
	}

	bool is_frag_type_visible(int f_idx) const
	{
		int i;
		if (frag_type_idxs_not_visible.size() == 0)
			return true;

		for (i=0; i<frag_type_idxs_not_visible.size(); i++)
			if (frag_type_idxs_not_visible[i] == f_idx)
				return false;
		return true;
	}


	void add_fragment(BreakageFragment& brf)
	{
		fragments.push_back(brf);
		total_intensity += brf.intensity;
		num_frags_detected++;
	}


	void remove_fragment(int frag_type_idx)
	{
		const int pos = get_position_of_frag_idx(frag_type_idx);
		if (pos<0)
		{
			cout << "Error: removing fragment type that isn't there:" << frag_type_idx << endl;
			exit(1);
		}

		if (fragments[pos].mass>0)
		{
			total_intensity -= fragments[pos].intensity;
			num_frags_detected--;

			vector<BreakageFragment> new_fragments;
			int j;
			for (j=0; j<fragments.size(); j++)
			{
				if (j==pos)
					continue;
				new_fragments.push_back(fragments[j]);
			}
			fragments = new_fragments;
		}
	}

	void clear()
	{
		region_idx=-1;
		mass=-1;
		total_intensity=0;
		num_frags_detected=0;
		score=0;
		frag_type_idxs_not_visible.clear();
		fragments.clear();
	}

	
	void print() const;

	void print(Config *config, ostream& os = cout) const;
	void print_fragments(Config *config, ostream& os = cout) const;


	int region_idx;  // model region
	mass_t mass;
	intensity_t total_intensity;
	int num_frags_detected;

	int parent_charge;   // for model_size
	int parent_size_idx; // for model_size

	score_t score;

	vector<int> frag_type_idxs_not_visible; // idxs of all frag types of the breakages region
											// that are not visible
	vector<BreakageFragment> fragments;
};


struct BreakageInfo {
	BreakageInfo() : breakage(NULL), n_break(NULL), c_break(NULL), n_var_ptr(NULL), c_var_ptr(NULL),
					 n_edge_idx(-1), n_var_idx(-1), c_edge_idx(-1), c_var_idx(-1), n_aa(Gap), c_aa(Gap), nn_aa(-1), cc_aa(-1),
					 n_edge_is_single(false), c_edge_is_single(false), connects_to_N_term(false), 
					 connects_to_C_term(false), preferred_digest_aa_C_term(false), preferred_digest_aa_N_term(false),
					 missed_cleavage(false), ind_n_edge_overlaps(false), ind_c_edge_overlaps(false), exp_n_edge_mass(0), exp_c_edge_mass(0), 
					 n_side_cat(NEG_INF), c_side_cat(NEG_INF), span_cat(NEG_INF), n_double_span_cat(NEG_INF),
					 c_double_span_cat(NEG_INF), score(NEG_INF),
					type(-1), node_idx(-1) {};

	const Breakage *breakage;
	const Breakage *n_break,*c_break;
	int *n_var_ptr, *c_var_ptr;
	int n_edge_idx, n_var_idx;
	int c_edge_idx, c_var_idx;
	int n_aa, c_aa;
	int nn_aa, cc_aa;
	bool n_edge_is_single, c_edge_is_single;
	bool connects_to_N_term;
	bool connects_to_C_term;
	bool preferred_digest_aa_C_term;
	bool preferred_digest_aa_N_term;
	bool missed_cleavage;
	bool ind_n_edge_overlaps;
	bool ind_c_edge_overlaps;

	mass_t exp_n_edge_mass;
	mass_t exp_c_edge_mass;

	int  n_side_cat, c_side_cat, span_cat;
	int  n_double_span_cat, c_double_span_cat;

	float score; // the score for the n_aa, c_aa combination

	int	  type;
	int	  node_idx;

	void print(Config *config) const;
};



typedef enum { NODE_REG, NODE_N_TERM, NODE_C_TERM, NODE_DIGEST} NODE_TYPES;

typedef enum { EDGE_REG, EDGE_FROM_N_TERM, EDGE_TO_C_TERM, EDGE_DIGEST } EDGE_TYPES;

struct ScoreComboLoc {
	ScoreComboLoc(int _in_edge, int _out_edge, int _in_var, int _out_var) : in_edge(_in_edge), out_edge(_out_edge), 
													in_var(_in_var), out_var(_out_var) {}
	ScoreComboLoc(const BreakageInfo& info) : in_edge(info.n_edge_idx), out_edge(info.c_edge_idx),
								in_var(info.n_var_idx), out_var(info.c_var_idx) {}

	bool operator< (const ScoreComboLoc& other) const
	{
		if (in_edge>other.in_edge)
			return false;
		return  ((in_edge<other.in_edge) ||
				 (in_edge == other.in_edge && out_edge< other.out_edge) ||
				 (in_edge == other.in_edge && out_edge == other.out_edge && in_var<other.in_var) ||
				 (in_edge == other.in_edge && out_edge == other.out_edge && in_var==other.in_var && out_var<other.out_var));
	}

	int in_edge, out_edge;
	int in_var, out_var;
};

#define MAX_EDGE_SIZE 3

class PrmGraph;

struct Node {
	Node() : mass(-1.0), score(0), const_score(NEG_INF), tmp_score(NEG_INF), log_rank(999999), idx_max_in_score_node(-1), 
			 idx_max_out_score_node(-1), type(-99), active(1),source_frag_type_idx(-1) {};

	bool operator< (const Node& other) const
	{
		return mass<other.mass;
	}

	void remove_in_edge_idx(int idx)
	{
		int i;
		for (i=0; i<in_edge_idxs.size(); i++)
			if (in_edge_idxs[i] == idx)
				break;
		if (i== in_edge_idxs.size())
			return;

		in_edge_idxs[i]=in_edge_idxs[in_edge_idxs.size()-1];
		in_edge_idxs.pop_back();
	}

	void remove_out_edge_idx(int idx)
	{
		int i;
		for (i=0; i<out_edge_idxs.size(); i++)
			if (out_edge_idxs[i] == idx)
				break;
		if (i== out_edge_idxs.size())
			return;

		out_edge_idxs[i]=out_edge_idxs[out_edge_idxs.size()-1];
		out_edge_idxs.pop_back();
	}


	void print(Config *config, ostream& os = cout) const;

	void print_combo_table(const PrmGraph *prm, ostream& os = cout) const;

	mass_t mass;
	score_t score;
	score_t const_score; // constant element in the scoring of combos

	score_t tmp_score;

	float log_rank;
	int   idx_max_in_score_node;  // the idx of the highest scoring previous node that connects to this node
	int   idx_max_out_score_node; // the idx of the highest scoring next node theat connects to this node
	int type;
	int active;                             // should this node be used
	int source_frag_type_idx;  // the fragment according to which the node was created
	Breakage breakage;
	vector<int> in_edge_idxs, out_edge_idxs;

	map<ScoreComboLoc,score_t> score_combos;
	vector<score_t> const_strong_exps;		// temporary values stored for combo scores
	vector<score_t> const_regular_exps;		// temporary values stored for combo scores
};




// Contains all the amino acid combos of the same length between two nodes
struct MultiEdge {
	MultiEdge() :  num_aa(0), type(0), max_variant_score(NEG_INF),  n_idx(-1), c_idx(-1), 
				  n_break(NULL), c_break(NULL), ind_edge_overlaps(false) {};

	int get_num_variants() const { return variant_ptrs.size(); }

	bool has_variant(int num_aa, const int *aas) const
	{
		int v;
		for (v=0; v<variant_ptrs.size(); v++)
		{
			if (*variant_ptrs[v] == num_aa )
			{
				int *var_aas = variant_ptrs[v]+1;
				int j;
				for (j=0; j<num_aa; j++)
					if (aas[j] !=  var_aas[j])
						break;

				if (j==num_aa)
					return true;
			}
		}
		return false;
	}

	int get_variant_idx(int num_aa, const int *aas) const
	{
		int v;
		for (v=0; v<variant_ptrs.size(); v++)
		{
			if (*variant_ptrs[v] == num_aa )
			{
				int *var_aas = variant_ptrs[v]+1;
				int j;
				for (j=0; j<num_aa; j++)
					if (aas[j] !=  var_aas[j])
						break;

				if (j==num_aa)
					return v;
			}
		}
		return -1;
	}


	int get_variant_idx(int first_aa) const
	{
		int v;
		for (v=0; v<variant_ptrs.size(); v++)
		{
			int *var_aas = variant_ptrs[v]+1;
			if (first_aa == var_aas[0])
				return v;
		}
		return -1;
	}

	int num_aa;
	int type;  // 
	
	score_t max_variant_score;

	int n_idx, c_idx;
	Breakage *n_break, *c_break;

	bool ind_edge_overlaps;  // if true there is a subpath with shorter edges that 
							 // has the same amino acids as this edge

	vector<int *>     variant_ptrs;
	vector<score_t>   variant_scores; // these are the deltas to the edge scores
	vector<float>	  variant_probs;  // these are filled using the AminoAcidProbs models
};


struct PathPos {
	PathPos() : mass(-1), breakage(NULL), node_idx(-1), edge_idx(-1), variant_ptr(NULL),
			aa(-1), node_score(0), edge_variant_score(0), edge_variant_prob(-1) {};

	void print() const
	{
		cout << "mass\t" << mass << endl;
		cout << "node\t" << node_idx << endl;
		cout << "edge\t" << edge_idx << endl;
		cout << "var_ptr\t" << variant_ptr << endl;
		cout << "aa\t" << aa << endl;
		cout << "score\t" << node_score << endl;
		cout << "prob\t" << edge_variant_prob << endl;
	}

	mass_t mass;

	Breakage *breakage; 

	int node_idx; // the index of the node N-terminal to the aa (-1 means no node is used
				  // as in the case where we are in the middle of a multiple aa edge

	int edge_idx; // the idx of the edge leaving the node (-1 if there is no edge)

	int* variant_ptr;

	int aa; // -1 means the is no aa (should only happen for the last PathPos in the path

	score_t node_score; // this is derived from the edge score + adjustments due to the 
						// specific amino combos in this edge

	score_t edge_variant_score;

	float	edge_variant_prob; // the probability assigned to the specific amino acid combo of the
							   // edge (AminoAcidProbs models calculate this number)
};






struct AA_combo {
	AA_combo() : total_mass(0), num_aa(0) {};

	bool operator< (const AA_combo& other) const
	{
		return (total_mass<other.total_mass);
	}

	void print(ostream& os = cout) const
	{
		int i;
		os << setprecision(4) << total_mass << " ";
		for (i=0; i<num_aa; i++)
			os << amino_acids[i] << " ";
		os << endl;
	}

	mass_t total_mass;
	int num_aa;
	int num_variants;
	int variant_start_idx; // the index in the variant_vector where this edges
						   // permutations are listed, the permutations are all concatanted
						   // in the same vector of int, first comes the number of amino acids
						   // then the amino acids themselves are listed
	int amino_acids[MAX_EDGE_SIZE];
};



struct score_pair {
	score_pair() : idx(int(NEG_INF)), score(NEG_INF) {};
	score_pair(int _i, float _n) : idx(_i), score(_n) {};
	bool operator< (const score_pair& other) const
	{
		return score>other.score;
	}
	int idx;
	float score;
};




#endif




