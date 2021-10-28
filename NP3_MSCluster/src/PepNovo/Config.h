#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "ConversionTables.h"
#include "Fragmentation.h"
#include "BasicDataStructs.h"
#include "SpectraAggregator.h"

// the Config class holds all configuration variables that are used
// by the models and programs: aa's, PTMs, thresholds etc


typedef enum PTM_REGIONS { PTM_ALL, PTM_N_TERMINAL, PTM_C_TERMINAL, PTM_POSITION } PTM_REGIONS;

typedef enum PTM_TYPES   { PTM_FIXED, PTM_OPTIONAL } PTM_TYPES;

#define NON_SPECIFIC_DIGEST 0
#define TRYPSIN_DIGEST 1

#define MINIMAL_NUMBER_SPECTRA_FOR_FRAGMENT_SELECTION 100


struct AA_combo;   // decleration in BasicDataStructs.h
class Peptide;


// label for fixed PTM stays the same as the orginal amino acid
struct PTM {
	PTM(int _org_aa = Gap, double _delta = 0.0,  int _type = -1, 
		int _region = PTM_ALL, int _position=0, string _label = "", string _name = "") : 
			org_aa(_org_aa), delta(_delta), type(_type), region(_region), position(_position),
			label(_label), name(_name) {};
	int org_aa;
	double delta;

	int type; //
	int region;

	int position; // used only for specific position like Q-17 at the +1 position
	
	string label;
	string name;
};

ostream& operator << (ostream& os, const PTM& ptm);



struct PTM_list {
	friend class Config;
public:

	void clear() { list.clear(); }

	int get_num_PTMs() const { return list.size(); }


	// returns the idx of the label
	// -1 is returned if label is not found
	int get_PTM_idx(const string& label) const;

private:
	
	// adds a PTM, makes sure that it wasn't already added
	// (if already added, returns false
	bool add_PTM(const PTM& ptm);

	vector<PTM> list;
};


// A class that holds a map of mass ranges
class MassRangeMap {
public:
	MassRangeMap() : was_initialized(0) { clear(); }

	void clear(mass_t max_mass = 400); 

	// adds a mass range to the map, if already covered, does nothinh
	void add_range(mass_t min_mass, mass_t range);

	// adds a new set of ranges *shift_size* away from the current existing ranges
	void add_shifted_ranges(mass_t shift_size);

	// checks if the given mass is in one of the maps ranges
	bool is_covered(mass_t mass) const
	{
		if (mass>max_map_mass)
			return true;

		MASS_T_MAP::const_iterator it;
		it = ranges.lower_bound(mass);
		if (it != ranges.end() && it != ranges.begin())
		{
			if (it->first == mass)
				return true;

			it--;
			if (it->first<= mass && it->second>= mass)
				return true;
		}
		return false;
	}

	void print(ostream& os = cout) const;

	void read_ranges(istream& is);
	void write_ranges(ostream& os = cout) const;

	void set_was_initialized(int flag) { was_initialized = flag; }
	int  get_was_initialized() const { return was_initialized; }

private:
	int was_initialized; 
	mass_t max_map_mass;
	MASS_T_MAP ranges;
};




class Config {
public:
	
	Config() : indWasInitialized_(false), keepOriginalDatasetIdx_(false), resource_dir("Models"), trainingDir_("Training"), 
		       min_exclude_range(9999999), max_exclude_range(NEG_INF) {};

	// sets the values for all defined PTMs
	// any additional PTMs can only be user defined

	bool getIndWasInitialized() const { return indWasInitialized_; }

	bool getKeepOriginalDatasetIdx() const { return keepOriginalDatasetIdx_; }
	void setKeepOriginalDatasetIdx(bool b) { keepOriginalDatasetIdx_ = b; }

	void init_with_defaults();

	void init_regional_fragment_set_defaults(int type_set = 0, int max_charge  =5);

	void print_supported_PTMs() const;

	void apply_selected_PTMs(const char *ptm_line);

	void apply_site_input_PTMs(const vector<string>& ptm_lines);

	void read_PTM_file(char *file);

	void fill_allowed_double_edges(bool allow_all = false);

	// returns the idx of an aa from its label
	// -1 if label is not found
	int get_aa_from_label(const string& label) const;

	int get_aa_with_position_from_label(const string& label, int position) const;

	int get_max_session_aa_idx() const;

	int get_frag_idx_from_label(const string& label) const { return (all_fragments.find_idx_from_label(label)); }

	const ConversionTables& get_session_tables() const { return session_tables; }

	const vector<int>& get_session_aas() const { return session_aas; }

	const vector<int>& get_char2aa() const { return session_tables.get_char2aa(); }

	const vector<char>&   get_aa2char() const { return session_tables.get_aa2char(); }
	const vector<string>& get_aa2label() const { return session_tables.get_aa2label(); }
	const vector<mass_t>& get_aa2mass() const { return session_tables.get_aa2mass(); }
	const vector<mass_t>& get_char2mass() const { return session_tables.get_char2mass(); }
	const vector<int>&    get_org_aa() const { return session_tables.get_org_aa(); }
	const vector<int>& get_aa_positions() const { return session_tables.get_aa_positions(); }


	int    get_digest_type() const { return digest_type; }
	void   set_digest_type(int type);
	
	const vector<int>& get_n_term_digest_aas() const { return n_term_digest_aas; }
	const vector<int>& get_c_term_digest_aas() const { return c_term_digest_aas; }
	int   get_num_n_term_digest_aas() const { return n_term_digest_aas.size(); }
	int   get_num_c_term_digest_aas() const { return c_term_digest_aas.size(); }
	bool  is_n_digest_aa(int aa) const { int i; for (i=0; i<n_term_digest_aas.size(); i++) if (aa==n_term_digest_aas[i]) break; return (i<n_term_digest_aas.size()); }
	bool  is_c_digest_aa(int aa) const { int i; for (i=0; i<c_term_digest_aas.size(); i++) if (aa==c_term_digest_aas[i]) break; return (i<c_term_digest_aas.size()); }


	int    get_need_to_estimate_pm() const { return need_to_estimate_pm; }
	void   set_need_to_estimate_pm(int val) { need_to_estimate_pm = val; }
	
	int    get_use_spectrum_charge() const { return use_spectrum_charge; }
	void   set_use_spectrum_charge(int val) { use_spectrum_charge = val; }
	
	int    get_use_spectrum_mz() const { return use_spectrum_mz; }
	void   set_use_spectrum_mz(int val) { use_spectrum_mz = val; }
	
	int    get_need_to_normalize() const { return need_to_normalize; }
	void   set_need_to_normalize(int val) { need_to_normalize = val; }

	int	   get_filter_flag() const { return filter_flag; }
	void   set_filter_flag(int val) { filter_flag = val; }

	int    get_mass_spec_type() const { return mass_spec_type; }

	int    get_itraq_mode() const { return itraq_mode; }
	void   set_itraq_mode(int mode) { itraq_mode = mode; }

	score_t get_terminal_score() const { return terminal_score; }
	void    set_terminal_score(score_t val)  { terminal_score=val; }

	score_t get_digest_score() const { return digest_score; }
	void    set_digest_score(score_t val)  { digest_score=val; }

	score_t get_forbidden_pair_penalty() const { return forbidden_pair_penalty; }
	void    set_forbidden_pair_pennalty(score_t val) { forbidden_pair_penalty=val; }

	int		get_max_edge_length() const { return max_edge_length; }
	void    set_max_edge_length(int val) { max_edge_length = val; }
	
	mass_t getTolerance() const { return tolerance; }
	mass_t get_pm_tolerance() const { return pm_tolerance; }

	void set_tolerance(mass_t t) { tolerance = t; }
	void setPrecursorMassTolerance(mass_t t) { pm_tolerance = t; }
	void setTolerances(mass_t t) { tolerance = t; pm_tolerance = t; }

	// NP3 GOT set and get rt tolerance
	float get_rt_tolerance() const { return rt_tolerance; }
	void set_rt_tolerance(float t) { rt_tolerance = t;}
	// NP3 scale factor
	float get_scale_factor() const { return scale_factor; }
	void set_scale_factor(float f) { scale_factor = f;}

	mass_t get_max_n_term_mod() const { return max_n_term_mod; }
	mass_t get_max_c_term_mod() const { return max_c_term_mod; }
	mass_t get_min_n_term_mod() const { return min_n_term_mod; }
	mass_t get_min_c_term_mod() const { return min_c_term_mod; }

	mass_t get_local_window_size() const { return local_window_size; }
	int    get_max_number_peaks_per_local_window() const { return max_number_peaks_per_local_window; }
	void   set_max_number_peaks_per_local_window(int n) { max_number_peaks_per_local_window = n; }
	int    get_number_of_strong_peaks_per_local_window() const { return number_of_strong_peaks_per_local_window; }
	int    get_max_charge_for_size() const { return max_charge_for_size; }
	void   set_max_charge_for_size(int max_charge) { max_charge_for_size = max_charge; }
	const vector< vector< mass_t > >& get_size_thresholds() const { return massThresholdsForSizes_; }
	const vector< mass_t >& get_region_thresholds() const { return region_thresholds; }
	int get_num_sizes(int charge) const { return regionalFragmentSets_[charge].size(); }
	int get_num_regions(int charge, int size_idx) const { return regionalFragmentSets_[charge][size_idx].size(); }

	int determine_charge(Peptide& pep, mass_t m_over_z) const;

	string get_resource_dir() const { return resource_dir; }
	string getTrainingDir() const { return trainingDir_; }
	string get_config_file() const { return config_file; }
	string getModelName() const { return model_name; }
	string get_fragments_file() const { return fragments_file; }

	string get_regional_fragment_sets_file() const { return regional_fragment_sets_file; }
	string get_aa_combo_file() const { return aa_combo_file; }
	void set_resource_dir(string _resource_dir) { resource_dir = _resource_dir; }
	void setTrainingDir(string dir) { trainingDir_ = dir; }
	void set_config_file(string _config_file) { config_file = _config_file; }
	void set_model_name(string name) { model_name = name; }
	string get_model_name() const { return model_name; }
	void set_fragments_file(string file) { fragments_file = file; }
	void set_regional_fragment_sets_file(string& file) { regional_fragment_sets_file = file; }
	void set_aa_combo_file(string file) { aa_combo_file = file; }

	const vector< vector<int> >& get_aa_variants() const { return aa_variants; }

	bool is_allowed_double_edge(int aa1, int aa2) const { return allowed_double_edge[aa1][aa2]; }

	const vector< vector<bool> >& get_allowed_double_edge() const { return allowed_double_edge; }
	const vector< vector<bool> >& get_double_edge_with_same_mass_as_single() const { return double_edge_with_same_mass_as_single; }

	void computeSizeThresholds(const SpectraAggregator& sa);

	bool selectFragmentIonTypes(const SpectraAggregator& sa,
								int maxNumberFragmentsPerRegion =0, 
								float minFragmentProbability = 0.05);

	void set_size_thresholds_according_to_set_of_masses(int charge,
											vector<mass_t>& spectra_masses);

	mass_t get_min_mass_for_size_idx(int charge, int size_idx) const
	{
		mass_t min=0;
		if (massThresholdsForSizes_[charge].size()>0 && size_idx>0)
			min=massThresholdsForSizes_[charge][size_idx-1];
		return min;
	}

	mass_t get_max_mass_for_size_idx(int charge, int size_idx) const
	{
		mass_t max=POS_INF;
		if (massThresholdsForSizes_[charge].size()>0 && size_idx<massThresholdsForSizes_[charge].size())
			max=massThresholdsForSizes_[charge][size_idx];
		return max;
	}
	
	// returns the size that depends on the charge and total mass
	int    calc_size_idx(int charge, mass_t pm_mass_with_19) const
	{
		const int charge_idx = (charge>max_charge_for_size) ? max_charge_for_size : charge;
		if (charge>=massThresholdsForSizes_.size())
			return 0;
		int i;
		for (i=0; i<massThresholdsForSizes_[charge_idx].size(); i++)
			if (pm_mass_with_19 < massThresholdsForSizes_[charge_idx][i])
				break;
		return i;
	}

	void print_size_thresholds() const;

	// returns the region idx for a (breakge) mass
	int	calc_region_idx(mass_t break_mass, mass_t pm_with_19, int charge,
					 mass_t min_peak_mass, mass_t max_peak_mass) const;


	void set_fragment_type_set(const FragmentTypeSet& fts) { all_fragments = fts; }

	void addFragmentTypes(const FragmentTypeSet& fts);

	const vector<FragmentType>& get_all_fragments() const { return all_fragments.get_fragments(); }


	// accessing the fragment sets
	// these are the fragments that are applicable for a given charge/size_idx/region_idx
	const vector<int>& get_regional_fragment_type_idxs(int charge, int size_idx, int region_idx) const
	{
		if (charge > max_charge_for_size)
			charge = max_charge_for_size;

		if (size_idx>= regionalFragmentSets_[charge].size())
			size_idx = regionalFragmentSets_[charge].size()-1;

		if (region_idx >= regionalFragmentSets_[charge][size_idx].size())
			region_idx = regionalFragmentSets_[charge][size_idx].size()-1;

		return regionalFragmentSets_[charge][size_idx][region_idx].get_frag_type_idxs();
	}

	void set_regional_random_probability(int charge, int size_idx, int region_idx, float p)
	{
		regionalFragmentSets_[charge][size_idx][region_idx].set_rand_prob(p);
	}

	float get_regional_random_probability(int charge, int size_idx, int region_idx)
	{
		return regionalFragmentSets_[charge][size_idx][region_idx].get_rand_prob();
	}

	const FragmentType& get_fragment(int frag_idx) const
	{
		return all_fragments.get_fragment(frag_idx);
	}

	const vector < vector< vector< RegionalFragments> > >& get_regional_fragment_sets() const
	{ return regionalFragmentSets_; }

	const RegionalFragments& get_regional_fragments(int charge, int size_idx, int region_idx) const
	{
		if (charge > max_charge_for_size)
			charge = max_charge_for_size;

		if (size_idx >= regionalFragmentSets_[charge].size())
			size_idx = regionalFragmentSets_[charge].size()-1;

		if (region_idx>= regionalFragmentSets_[charge][size_idx].size())
			region_idx = regionalFragmentSets_[charge][size_idx].size()-1;

		return regionalFragmentSets_[charge][size_idx][region_idx];
	}

	RegionalFragments& get_non_const_regional_fragments(int charge, int size_idx, int region_idx)
	{
		if (charge > max_charge_for_size)
			charge = max_charge_for_size;

		if (size_idx >= regionalFragmentSets_[charge].size())
			size_idx = regionalFragmentSets_[charge].size()-1;

		if (region_idx>= regionalFragmentSets_[charge][size_idx].size())
			region_idx = regionalFragmentSets_[charge][size_idx].size()-1;

		return regionalFragmentSets_[charge][size_idx][region_idx];
	}

	void clear_combos(int charge, int size_idx, int region_idx)
	{
		regionalFragmentSets_[charge][size_idx][region_idx].frag_type_combos.clear();
	}

	void set_all_regional_fragment_relationships();
	void print_all_regional_fragment_relationships() const;


	// all strong fragments (that are strong in at least one regional fragment set)
	const vector<int>& get_all_strong_fragment_type_idxs() const
	{
		return all_strong_fragment_type_idxs;
	}

	// selects the fragments to be used, uses a cutoff that is X times random prob)
	void select_fragments_in_sets(score_t X=1.0, int max_num_frags=0);

	// selects the fragments to be used, uses a cutoff that is X times random prob)
	void applyCuttoffsToRegionalSets(score_t X=1.0, int max_num_frags=0);

	void learnTolerancesFromData(const SpectraAggregator& sa, mass_t initalToleranceEstimate);


	// For each regional fragments selects all fragments that have a minimal probability
	// to be strong.
	void select_strong_fragments(int charge,
						score_t min_prob = 0.5, int max_num_strong = 3, bool verbose = true);

	void selectTwoOverallStrongFragments();

	int get_strong_type1_idx() const { return strong_type1_idx; }
	int get_strong_type2_idx() const { return strong_type2_idx; }


	void sort_accoriding_to_fragment_probs(vector<score_t>& frag_probs, int charge, 
										   int size_idx, int region_idx)
	{
		regionalFragmentSets_[charge][size_idx][region_idx].set_frag_probs(frag_probs);
		regionalFragmentSets_[charge][size_idx][region_idx].sort_according_to_frag_probs();
	}

	void print_session_aas() const;
	void print_fragments(ostream &os) const;
	void print_regional_fragment_sets(ostream& os = cout) const;
	void read_fragments(istream& is);
	void read_regional_fragment_sets(istream& is);
	void clone_regional_fragment_sets(int source_charge, int target_charge);


	void print_all_fragments() const { all_fragments.print(); }


	void read_config(const char* file_name);

	void write_config();

	void print_config_parameters(ostream& os = cout) const;


	// parses a line that is assumed to be from a config file
	// all parameters are assumed to start with
	// #CONF <PARAMETER_NAME> VALUES
	void parse_config_parameter(char *buff);


	// checks if the given mass falls within the allowed regions for prefix masses
	bool is_allowed_prefix_mass(mass_t mass) const
	{
		return  (allowed_node_masses.is_covered(mass));
	}

	// checks if the given mass falls within the allowed regions for suffix masses
	bool is_allowed_suffix_mass(mass_t pm_with_19, mass_t mass) const
	{
		return  (allowed_node_masses.is_covered(pm_with_19 - mass - MASS_OHHH));
	}

	// checks if the given mass falls within the allowed regions for amino acid masses
	bool is_allowed_aa_combo_mass(mass_t mass) const
	{
		return  (allowed_node_masses.is_covered(mass));
	}


	// initializes the allowed_prefix_masses map
	// and the allowed suffix masses map
	void init_allowed_node_masses(mass_t max_mass = 400.0);

	// calclates the masses of the different aa_combos
	// combos are sorted aa lists.
	void calc_aa_combo_masses();
	
	const vector<AA_combo>& get_aa_edge_combos() const { return aa_edge_combos; }
	const vector<int>& get_combo_idxs_by_length(int length) const { return combo_idxs_by_length[length]; }
	const int *get_first_variant_ptr(int combo_idx) const;
	const int *get_variant_ptr(int var_idx) const { return &variant_vector[var_idx]; }

	// returns an index and number of combos that exist in the given mass range
	int get_ptrs_for_combos_in_mass_range(mass_t min_mass, mass_t max_mass, 
												   int& num_combos) const;
		
	// returns true if there is a combo that contains the ordered variant of the given aas
	bool combos_have_variant(const vector<int>& combos, int num_aa, int *var_aas) const;


	// This is the maximal mass of a combo of amino acids that can be in an edge
	// (equals the maximal combo mass + tolerance)
	mass_t get_max_combo_mass() const { return max_combo_mass; }



	mass_t get_min_exclude_range() const { return min_exclude_range; }
	mass_t get_max_exclude_range() const { return max_exclude_range; }
	
	void add_exclude_range(mass_t min_range, mass_t max_range);

	bool check_if_mass_is_in_exclude_range(mass_t m) const
	{
		if (min_ranges.size() == 0)
			return false;

		int i;
		for (i=0; i<min_ranges.size(); i++)
			if (m>=min_ranges[i] && m<=max_ranges[i])
				return true;
		return false;
	}
	
	// calculates the aa_variants vectors (for terms & A-V)
	void set_aa_variants();

	void print_aa_variants() const;

	void print_size_partitions() const;
	

private:
	bool indWasInitialized_;

	bool keepOriginalDatasetIdx_;

	// These conversion tables represent the parameters after PTM modifications
	// All tables have the same aa indices however the actual masses might be
	// different (due to terminal additions for instance)
	ConversionTables session_tables; 

	vector<int> standard_aas; // the 20 unmodified amino acids
	vector<int> session_aas; //  all aas that can be part of a peptide (including terminal mods)


	vector< vector<int> > aa_variants; // holds for each of the amino acids A-V and the terminals
									   // all the variants that can be used for each amino acid (e.g. M-> M,M+16,M+32)
	// maps all labels to their aa_idx
	STRING2INT_MAP label2aa;


	// PTM vectors (hold info on all supported PTMs)
	bool       ind_read_PTM_file;
	
	PTM_list   all_fixed_PTMs;    // region PTM_ALL
	PTM_list   all_optional_PTMs; // region PTM_ALL
	PTM_list   all_terminal_PTMs; // must be either PTM_N_TERMINAL, PTM_C_TERMINAL



	// MASS SPEC TYPE
	int		   mass_spec_type; // type of machine, might influence the way things are done
							   // this parameter should be partof the model file, and not changed.

	int		    digest_type; // 0 - nothing , 1 Trypsin

	vector<int> n_term_digest_aas;
	vector<int> c_term_digest_aas;

	int        need_to_estimate_pm; // 0 - used file pm, 1 - use original pm

	int		   need_to_normalize; // 0 - don't normalize intensities, 1 - do normalize (sum of intensities = m_over_z)

	int		   itraq_mode; // 0 - no itraq, 1 - this is itraq data

	int		   use_spectrum_charge;

	int		   use_spectrum_mz;

	int		   filter_flag;

	score_t    terminal_score;
	score_t	   digest_score;
	score_t	   forbidden_pair_penalty;

	string     resource_dir;  // path to direcotry where resource files can be found
	string	   trainingDir_;
	string     fragments_file;
	string     regional_fragment_sets_file;
	string	   aa_combo_file;
	string	   model_name;   // the name of the model that uses this config
	string     config_file;

	mass_t max_n_term_mod; // for the current session_aas
	mass_t max_c_term_mod; // for the current session_aas
	mass_t min_n_term_mod;
	mass_t min_c_term_mod;


	// Tolerance variables
	mass_t tolerance;       // tolerance for consecutive peaks
	mass_t pm_tolerance;

	// NP3 GOT add rt tolerance and scale factor
	float rt_tolerance;
	float scale_factor;

	mass_t local_window_size;
	int	 max_number_peaks_per_local_window;
	int  number_of_strong_peaks_per_local_window;

//	score_t random_prob; // the random prob of observing a peak in the region (bins/area);

	// Idxs for selected PTMs to be used
	vector<int>  selected_fixed_PTM_idxs;
	vector<int>  selected_optional_PTM_idxs;
	vector<int>  selected_terminal_PTM_idxs;

	// Vectors that hold a sorted list of amino acid combos and masses
	vector< AA_combo >    aa_edge_combos;       // all lengths combined
	vector< vector<int> > combo_idxs_by_length; // the idxs in the aa_edge_combos_vector

	vector< int > variant_vector; // holds the variants for the various combos
								  // each variant is listed with the number of amino acids
								  // followed by the amino acids themselves

	vector< int > combo_start_idxs; // has the index of the first combo with a given integer mass

	mass_t max_combo_mass;
	int  max_edge_length;

	// Thresholds for determining what model to use
	int max_charge_for_size;       // the size of size_thresholds
	vector< vector< mass_t > > massThresholdsForSizes_; // charge / size 
	vector< mass_t >           region_thresholds;

	// The sets of fragments that are to be observed in the different models
	// these depend on the parent charge, the peptide size, and the region
	// in which the breakage is observed
	FragmentTypeSet all_fragments;    // all possible fragment types
	vector < vector< vector< RegionalFragments> > > regionalFragmentSets_; // charge, size, region_idx
	vector<int> all_strong_fragment_type_idxs;
	int strong_type1_idx; // the two strongest (highest probability fragments)
	int strong_type2_idx;

	// mass exclude ranges (when doing clustering, for example ignore iTraq peaks)
	vector< mass_t > min_ranges;
	vector< mass_t > max_ranges;
	mass_t min_exclude_range;
	mass_t max_exclude_range;


	
	// These maps contain the ranges of permitted masses for amino acids combos,
	// prefix node masses, and suffix node masses, repectively.
	// Separate range maps are used because terminal PTMs can change the allowed ranges
	MassRangeMap allowed_node_masses;


	// these are used in the selection of edges for the PrmGraph
	vector< vector<bool> > allowed_double_edge, double_edge_with_same_mass_as_single;


	void init_standard_aas();

	void init_model_size_and_region_thresholds();

	void print_table_aas(const ConversionTables& table, 
							 const vector<int>& aas) const;

	// determines the parent mass tolerance for which *cuttoff_prob* of the abundant fragments
	// are caught
	mass_t calculatePrecursorMassTolerance(const SpectraAggregator& sa, 
										   float cutoffProbability=0.98) const;

	// determines the tolerance for which *cuttoff_prob* of the abundant fragments
	// are caught
	mass_t calculateFragmentMassTolerance(const SpectraAggregator& sa, 
										  mass_t maxTolerance,
									      float cutoffProbability=0.96) const;

};
#endif


