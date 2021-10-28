#ifndef __FRAGMENTATION_H__
#define __FRAGMENTATION_H__

#include "PepNovo_includes.h"


typedef enum { PREFIX, SUFFIX } FRAG_ORI;


// This class stores combinations of fragments that are later used in the
// model's scoring tables

class Config;

struct FragmentCombo {

	FragmentCombo() { frag_inten_idxs.clear(); frag_no_inten_idxs.clear(); }
	
	void print_combo(const Config *config, ostream& os = cout) const;
	void read_combo(Config *config, istream& is);

	vector<int> frag_inten_idxs;
	vector<int> frag_no_inten_idxs;
};



/*****************************************************************************

  This is the basic class for describing a fragment ion type (e.g. b , y, b2, etc.)
******************************************************************************/
struct FragmentType {

	FragmentType() : orientation(-1), charge(0), offset(0), parent_frag_idx(-1), 
					 offset_from_parent_frag(0),   prob(0), spec_count(0), 
					 ind_is_a_strong_fragment(0), label("") {};

	
	bool operator< (const FragmentType& other) const
	{
		return prob > other.prob;
	}

	bool operator< (FragmentType& other)
	{
		return prob > other.prob;
	}

	// calcs the epxected position of the peak
	mass_t calc_expected_mass(mass_t breakage_mass, mass_t pm_with_19) const
	{
		const mass_t ps_mass = (orientation == PREFIX) ? breakage_mass :
						pm_with_19 - breakage_mass - MASS_OHHH;

		if (charge == 1)
			return ps_mass + offset;

		return (ps_mass/ charge) + offset;
	}

	// calcs the expected mass of the breakage, given the mass of the fragment's peak
	mass_t calc_breakage_mass(mass_t peak_mass, mass_t pm_with_19) const
	{
		const mass_t ps_mass = (peak_mass - offset) * charge;

		if (orientation == PREFIX)
			return ps_mass;

		return (pm_with_19 - ps_mass - MASS_OHHH);
	}


	void read_fragment(istream& is);

	void write_fragment(ostream& os = cout) const;

	void make_frag_label(mass_t tolerance);			

	
	int		 orientation;   // PREFIX / SUFFIX

	int		 charge;       // +1, +2, +3, etc.

	mass_t   offset;      // from prefix/suffix mass

	int		 parent_frag_idx;  // the frag with the highest prob that has
							   // the same charge and orientation

	mass_t   offset_from_parent_frag; // the mass difference between the parent frag
									  // (takes into account the charge)

	score_t  prob;  // the probability of observing this fragment in a true breakage
					  // (this is in the context of the a specific FragmentSet that has
					  // a specific charge, size_idx, region_idx

	int		spec_count; // used when learning

	int     ind_is_a_strong_fragment;

	string   label;
};





// Holds a set of fragments
// It is assumed, but not enforced, that the fragments are entered in 
// some topological order, i.e., parent fragments appear before their
// derivatives (i.e., no loss before losses, etc.)
class FragmentTypeSet {
	friend class Config;
public:

	int get_num_fragments() const { return fragments.size(); }

	void add_fragment_type(const FragmentType& ft) 
	{ 
		if (find_idx_from_label(ft.label)<0)
		{
			fragments.push_back(ft); 
			label2idx.insert(STRING2INT_MAP::value_type(ft.label,fragments.size()-1));
		}
	}

	void output_fragment_labels(ostream& os = cout) const;

	void print() const;

	void clear_set() { fragments.clear(); base_frag_idxs.clear(); label2idx.clear(); }

	// returns idx that has the same label
	int find_idx_from_label(const string& label) const
	{
		STRING2INT_MAP::const_iterator iter = label2idx.find(label);

		if (iter == label2idx.end())
			return -1;

		return (*iter).second;
	}

	void sort_fragments_according_to_probs();

	// removes fragments that appear to be isotopic peaks of previously
	// selected fragments such as b+1, y+2 etc.
	void remove_isotopic_fragments(mass_t tolerance, bool vebose = true);

	void set_parent_frag_idxs();

	// recreates the fragments vector by selecting only fragments
	// whose idxs appear in the frag_idxs vector
	void select_fragment_types(const vector<int>& frag_idxs);

	const FragmentType& get_fragment(int frag_idx) const { return fragments[frag_idx]; }

	FragmentType& get_non_const_fragment(int frag_idx) { return fragments[frag_idx]; }

	const vector<FragmentType>& get_fragments() const    { return fragments; }

	
private:

	vector<FragmentType> fragments;
	vector<int>          base_frag_idxs;
	STRING2INT_MAP label2idx;
};



// holds the fragments that are applicable to a specific region
class RegionalFragments {
	friend class Config;
public:

	void init_pepnovo_types(int charge, Config* config);

	void init_with_all_types(int charge, Config *config);

	int get_num_fragments() const { return frag_type_idxs.size(); }

	// returns the position in the frag_idxs vecctor of the desired fragment
	// return -1 if the frag_idx is not found in frag_idxs
	int get_position_of_frag_type_idx(int frag_type_idx) const
	{
		int i;
		for (i=0; i<frag_type_idxs.size(); i++)
			if (frag_type_idxs[i] != frag_type_idx)
				return i;
		return -1;
	}

	int get_frag_type_idx_with_highest_position(const vector<int>& frag_idxs) const
	{
		int best_pos = 9999999;
		int best_frag = -1;
		int i;
		for (i=0; i<frag_idxs.size(); i++)
		{
			int pos = get_position_of_frag_type_idx(frag_idxs[i]);
			if (pos>=0 && pos<best_pos)
			{
				best_pos = pos;
				best_frag = frag_idxs[i];
			}
		}
		return best_frag;
	}


	bool is_a_strong_frag_type(int f_idx) const
	{
		int i;
		for (i=0; i<strong_frag_type_idxs.size(); i++)
			if (f_idx == strong_frag_type_idxs[i])
				return true;
		return false;
	}

	const vector<int>& get_frag_type_idxs() const { return frag_type_idxs; }
	const vector<int>& get_strong_frag_type_idxs() const { return strong_frag_type_idxs; }
	const vector<FragmentCombo>& get_frag_type_combos() const { return frag_type_combos; }
	const vector< vector<int> >& get_mirror_frag_idxs() const { return mirrorFragmentIndexes_; }
	const vector< vector<int> >& get_parent_idxs() const { return parentFragmentIndexes_; }
	const vector< vector<int> >& get_parents_with_same_charge_ori_idxs() const { return parentFragmentsWithSameChargeAndOrientation_; }
	const vector<score_t>& get_frag_probs() const { return frag_probs; }
	score_t get_frag_prob(int frag_type_idx) const
	{
		int i;
		for (i=0; i<frag_type_idxs.size(); i++)
			if (frag_type_idxs[i]==frag_type_idx)
				return frag_probs[i];
		return NEG_INF;
	}
	float   get_rand_prob() const { return rand_prob; }

	void    set_rand_prob(float p) { rand_prob = p; }

	void select_fragments_with_minimum_prob(score_t min_prob, int max_num_frags=0);

	// this changes the indexes order and invalidates other lists
	// such as strong_frag_idxs and fragment_combos
	void sort_by_prob();

	// checks if the fragments have a parent fragment in the frag_idxs and 
	// give its idx, otherwise -1
	void set_parent_frag_idxs(Config *config);

	void set_frag_probs(vector<score_t>& probs) { frag_probs=probs; }

	// makes the order of the fragments in the frag_idxs and frag_probs vectors
	// be in descending probability order
	void sort_according_to_frag_probs();


	// chooses all fragments that have high enough probability to be strong
	void select_strong_fragments(Config *config, score_t min_prob, int max_num_strong);

	void add_combo(FragmentCombo& combo) { frag_type_combos.push_back(combo); }

	
	// sets the idxs of the different relatioships between the fragments
	// (parent, alt_charge, mirror frag, alt_charge mirror frag)
	void set_fragment_relationships(Config *config);

	void print_fragment_relationships(const Config *config) const;

	void printFragments(ostream& os) const;
	

private:

	float   rand_prob;

	vector<int> frag_type_idxs; // the indexes of the fragments that are applicable in this region
	vector<score_t> frag_probs; // holds the probabilities of observing each fragment in this region
	vector<int> strong_frag_type_idxs; // the frag_idxs that are considerd strong in this region

	vector< vector<int> > mirrorFragmentIndexes_;
	vector< vector<int> > parentFragmentIndexes_;
	vector< vector<int> > parentFragmentsWithSameChargeAndOrientation_;
										
	vector<FragmentCombo>  frag_type_combos; // if any set of peaks meets the restrictions for inten / no_inten
											// in one of these combos, then a node should be created

};


#endif




