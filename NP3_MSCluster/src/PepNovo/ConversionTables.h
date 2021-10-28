#ifndef __CONVERSIONTABLES_H__
#define __CONVERSIONTABLES_H__

#include "PepNovo_includes.h"

class ConversionTables {
public:
	void init_for_standard_aas();

	bool add_optional_PTM_aa(int aa, const string& label, mass_t delta, int position);

	bool make_fixed_mod(int aa, mass_t delta);

	bool add_optional_PTM_terminal_aa(mass_t mass_t, int position, const string& ptm_label);


	// access to 
	int	   get_char2aa(char c) const { return char2aa[c]; }
	mass_t get_char2mass(char c) const { return char2mass[c]; }
	mass_t get_aa2mass(int a) const { return aa2mass[a]; }
	char   get_aa2char(int a) const { return aa2char[a]; }
	string get_aa2label(int a) const { return aa2label[a]; }
	int    get_org_aa(int a) const { return org_aa[a]; }
	int    get_aa_position(int a) const { return aa_positions[a]; }

	const vector<int>& get_char2aa() const { return char2aa; }
	const vector<mass_t>& get_char2mass() const { return char2mass; }
	const vector<mass_t>& get_aa2mass() const { return aa2mass; }
	const vector<char>& get_aa2char() const { return aa2char; }
	const vector<string>& get_aa2label() const { return aa2label; }
	const vector<int>& get_org_aa() const { return org_aa; }
	const vector<int>& get_aa_positions() const { return aa_positions; }

private:
	vector<int>       char2aa;
	vector<mass_t>    char2mass;
	vector<mass_t>    aa2mass;
	vector<char>      aa2char;
	vector<string>    aa2label;
	vector<int>       org_aa;
	vector<int>       aa_positions; // for each session aa idx holds 0 if the aa can appear anywhere,
							       // otherwise holds the position where we can have that aa for instance
							      // +1 for N-terminal mod aas or Q-17, -1 for C-terminal mods


};


#endif




