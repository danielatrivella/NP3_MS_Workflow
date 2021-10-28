#include "RegionalPepNovoModel.h"


score_t RegionalPepNovoModel::calc_breakage_score(Breakage *breakage, 
												  bool verbose, Config *config) const
{
	if (breakage->are_all_frag_types_visible())
	{
		int f;
		score_t score =0;
		for (f=0; f<single_breakage_tables.size(); f++)
		{
			if (verbose)
			{
		//		single_breakage_tables[f].print_pretty(config);
			}
			score_t frag_score = single_breakage_tables[f].get_score(breakage);
			score += frag_score;
			if (verbose)
			{
				cout << single_breakage_tables[f].get_fields()[0] << " " << frag_score << " " << score << endl;
			}
		}
		if (verbose)
			cout << endl;

		return score;
	}
	
	if (verbose)
	{
		int i;
		cout << "Not all frags visible:" << endl;
		for (i=0; i<breakage->frag_type_idxs_not_visible.size(); i++)
			cout << breakage->frag_type_idxs_not_visible[i] << " ";
		cout << endl;
	}
	//might need to use the independent peak tables for fragments whose parents 
	// are not visible
	score_t score =0;
	int f;
	for (f=0; f<single_breakage_tables.size(); f++)
	{
		const int frag_type_idx =  frag_type_idxs[f];

		// don't score fragments that are beyond visible range
		if (! breakage->is_frag_type_visible(frag_type_idx))
			continue;

		// check if all fragment's parent are visible
		int p;
		bool all_parents_visible = true;
		const vector<int>& fields = single_breakage_tables[f].get_fields();
		for (p=1; p<3; p++) // the first 3 fields are assumed to be from the current breakage
		{
			if (p<0)
				continue;

			if (! breakage->is_frag_type_visible(fields[p]))
			{
				all_parents_visible = false;
				break;
			}
		}

		score_t frag_score = all_parents_visible ? 
			single_breakage_tables[f].get_score(breakage) :
			independent_frag_tables[f].get_score(breakage);

		score += frag_score;
		if (verbose)
		{
			cout << frag_type_idx << " " << frag_score << (all_parents_visible? "viz" : "no_viz") << endl;
		}
		
		if (verbose)
			cout << endl;
	}
	return score;	
}





void RegionalPepNovoModel::convert_to_scores(const vector<double>& q_rand)
{
	int f;
	
	missed_cleavage_score = 0;
	for (f=0; f<independent_frag_tables.size(); f++)
	{
		independent_frag_tables[f].convert_to_score(q_rand);
		missed_cleavage_score += independent_frag_tables[f].get_score_prob_from_idx(0);
	}

	for (f=0; f<single_breakage_tables.size(); f++)
		single_breakage_tables[f].convert_to_score(q_rand);
}


void RegionalPepNovoModel::read_regional_model(Config *config, istream& is)
{
	int f,num_frags;
	char buff[128];

	is.getline(buff,128);
	if (sscanf(buff,"#REGIONAL_MODEL %d %d %d %d",&charge,&size_idx,&region_idx,&num_frags) != 4)
	{
		cout << "Error: bad line in regional model: "<< buff << endl;
		exit(1);
	}

	frag_type_idxs.resize(num_frags,-1);
	independent_frag_tables.resize(num_frags);
	single_breakage_tables.resize(num_frags);

	missed_cleavage_score =0;
	for (f=0; f<num_frags; f++)
	{
		independent_frag_tables[f].read_table(config,is);
		single_breakage_tables[f].read_table(config,is);
		
		frag_type_idxs[f] = independent_frag_tables[f].get_fields()[0];
	}

	calc_missed_cleavage_score();
	
}

void RegionalPepNovoModel::calc_missed_cleavage_score()
{
	int f;

	missed_cleavage_score=0;

	for (f=0; f<frag_type_idxs.size(); f++)
		missed_cleavage_score += independent_frag_tables[f].get_score_prob_from_idx(0);
}


void RegionalPepNovoModel::write_regional_model(ostream& os) const
{
	os << "#REGIONAL_MODEL " << charge << " " << size_idx << " " << region_idx <<
		" " << frag_type_idxs.size() << endl;

	int f;
	for (f=0; f<frag_type_idxs.size(); f++)
	{
		independent_frag_tables[f].write_table(os);
		single_breakage_tables[f].write_table(os);
	}

}


void RegionalPepNovoModel::print_table_names(const Config *config, ostream& os) const
{
	os << "Charge: " << charge << ", size " << size_idx << ", region " << region_idx << endl;
		int f;
	for (f=0; f<frag_type_idxs.size(); f++)
	{
		string name;
		single_breakage_tables[f].make_table_name(config,name);	
		os << setw(3) << left << f+1 << name << endl;
	}
	os << endl;
}




