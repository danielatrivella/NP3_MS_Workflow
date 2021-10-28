#ifndef __FRAGPROBTABLE_H__
#define __FRAGPROBTABLE_H__

#include "BasicDataStructs.h"
#include "Config.h"

#define NUM_TABLE_FIELDS 5

typedef int table_entry[NUM_TABLE_FIELDS];


// contains the probability table for a fragment.
// After all calculations are done, the probabilities
// are converted to scores by dividing by the random probability
// Might consider changing to a new object..
class FragProbTable {
	friend class RegionalPepNovoModel;
public:
	FragProbTable() :
	  model_tolerance(0.5), max_table_size(-1), score_probs_type(0) {};
	
	score_t get_score(Breakage *breakagre, Breakage *previous_breakage = NULL,
						mass_t mass_diff =0) const;

	score_t get_score_prob_from_idx(int idx) const { return score_probs[idx]; }

	const vector<int>& get_fields() const { return fields; }
	
	void read_table(Config *config, istream& is);
	void write_table(ostream& os) const;

	void print_pretty(Config *config, ostream& os = cout) const;

	

	void make_table_name(const Config *config, string& name) const;

	void init_fields(Config * config, vector<int>& fields, vector<int>& num_vals);

	void init_counts(double init_val = 0.5);

	void add_instance(Breakage *breakagre, Breakage *previous_breakage = NULL,
					  mass_t exp_mass_diff =0, bool verbose  = false);

	void calc_probs();

	double calc_dkl_sum(const vector<score_t>& ind_probs) const; // assumes probs and counts are set

	bool are_relevant_fragments_visible(Breakage *breakage) const
	{
		int f;
		for (f=0; f<fields.size(); f++)
			if (fields[f]>=0 && ! breakage->is_frag_type_visible(fields[f]))
				return false;
		return true;
	}

	// changes the values in the score_probs vector from probabilities to scores
	void convert_to_score(const vector<double>& rand_probs);

	void set_model_tolerance(mass_t tolerance) { model_tolerance = tolerance; }

private:

	Config *config;

	mass_t model_tolerance; // used for the peak offset values

	int max_table_size;   // the number of entries in the table

	int score_probs_type; // 0 - probabilities, 1 - scores

	int charge_of_frag; // the charge of frag 0

	vector<int> fields;   // 0 - frag_type_idx of table
						  // 1 - frag_type_idx of parent 1
						  // 2 - frag_type_idx of parent 2
						

	
	vector<int> field_multipliers;   // used for calculating an entry's index in the table

	vector<int> num_field_vals;    // the number of values possible for each field in the table
								  // should be either (2, 0 or 1 for all vals>0), or the number 
								 // of peak levels.

	vector<score_t> score_probs;  // the probabilities / scores

	vector<double> counts; // used for training





	// set the multipliers used to calculate the index for an entry
	// multipliers are based on the number of field values for each
	// participating field in the table
	void set_field_multipliers();


	// puts integer values into the entry based on the values of the 
	// tables fragments. If num_field vals == 2, then a binary value is
	// given, where 1 represents all values>0
	// the exp_offset is the expected amino acid mass difference between the breakages
	void fill_table_entry(Breakage *breakage, Breakage *previous_breakage,
						  table_entry& entry, mass_t exp_offset = 0) const;
	

	int calc_table_idx(table_entry& entry) const
	{
		int i;
		int idx;
		
		idx=entry[0];
		for (i=1; i<NUM_TABLE_FIELDS; i++)
			if (num_field_vals[i]>0)
				idx += entry[i] * field_multipliers[i];
		
		return idx;
	}
	

};


#endif

