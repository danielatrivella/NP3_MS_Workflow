#ifndef __PEPTIDE_COMP_H__
#define __PEPTIDE_COMP_H__

#include "BasicDataStructs.h"
#include "Config.h"


#define MAX_COMP_CAT 20 // 1-20
#define MAX_COMP_LEN 3  // 1-3

#define MAX_NUM_M16s 2


struct PeptideCompStats {
	
	PeptideCompStats() : num_aa(0) {};
	void clear_pcs();
	void print_pcs(ostream& os = cout) const;

	int num_aa;
	
	int start_comp[MAX_COMP_LEN+1];
	int end_comp[MAX_COMP_LEN+1];
	int cat_counts[MAX_COMP_LEN+1][MAX_COMP_CAT+1];
};



class PeptideCompAssigner {
public:
	PeptideCompAssigner() : was_initialized(false), config(NULL) {}

	void fill_peptide_stats(const Peptide& peptide, PeptideCompStats& stats) const;

	int get_aa_category(int num_aa, const int *aas, bool n_term, bool c_term) const;

	void read_and_init_from_tables(Config *_config,  const char *name);

	const string& getModelName() const { return model_name; }

	bool get_ind_was_initialized() const { return was_initialized; }

	void init_aa_translations();

private:
	bool was_initialized;
	Config *config;
	string model_name;
	vector<int> aa_translation; // converts aa to digit 1-19 (I=L)
	vector<vector<int> > start_assigns, end_assigns, mid_assigns; // length, aa_code

	

	int calc_aa_code(const int *pos, int length) const
	{
		if (length<1 || length>10)
			return NEG_INF;

		if (length==1)
			return aa_translation[*pos];

		if (length==2)
			return aa_translation[*pos]+20*aa_translation[*(pos+1)];

		return aa_translation[*pos]+20*aa_translation[*(pos+1)]+400*aa_translation[*(pos+2)];
	}

	
	void read_table_to_vector(char *file_path, int num_aa, vector<int>& vec);

};

#endif

