#include "FragProbTable.h"


void FragProbTable::set_field_multipliers()
{
	int i;
	field_multipliers.resize(NUM_TABLE_FIELDS,0);
	field_multipliers[0]=1;

	for (i=1; i<NUM_TABLE_FIELDS; i++)
		field_multipliers[i] = (this->num_field_vals[i-1]>0) ? 
			field_multipliers[i-1] * num_field_vals[i-1] : field_multipliers[i-1];

	max_table_size = 0;
	for (i=0; i<NUM_TABLE_FIELDS; i++)
		if (num_field_vals[i]>0)
			max_table_size = field_multipliers[i] * num_field_vals[i];
			
}


// 
void FragProbTable::write_table(ostream& os) const
{
	int i;

	if (score_probs_type == 0)
	{
		cout << "Error: must first convert probs to scores before writing table!" << endl;
		exit(1);
	}

	// fields and number of vals
	for (i=0; i<NUM_TABLE_FIELDS; i++)
		os << fields[i] << " " << num_field_vals[i] << " ";
	os << endl;
	
	os << max_table_size << endl;

	for (i=0; i< max_table_size; i++)
		os << setprecision(5) << score_probs[i] << endl;

}


void FragProbTable::print_pretty(Config *config, ostream& os) const
{
	int i,t_lengths=0;
	vector<int> lengths;
	
	lengths.resize(NUM_TABLE_FIELDS,0);
	for (i=0; i<NUM_TABLE_FIELDS; i++)
	{
		if (fields[i]>=0)
		{
			lengths[i] = config->get_fragment(fields[i]).label.length() + 1;
			if (lengths[i]<4)
				lengths[i]=4;
			t_lengths += lengths[i];
		}
	}

	t_lengths+=11;
	os <<"Idx Prob   ";
	for (i=0; i<NUM_TABLE_FIELDS; i++)
	{
		if (num_field_vals[i]>0)
			os << " " << setw(lengths[i]) <<  config->get_fragment(fields[i]).label;
	}
	os << endl;
	os << setw(t_lengths) << setfill('-') << "-" << endl << setfill(' ');

	vector<int> f;
	f.resize(NUM_TABLE_FIELDS,0);

	for (f[4]=0; f[4]<=num_field_vals[4]; f[4]++)
		for (f[3]=0; f[3]<=num_field_vals[3]; f[3]++)
			for (f[2]=0; f[2]<=num_field_vals[2]; f[2]++)
				for (f[1]=0; f[1]<=num_field_vals[1]; f[1]++)
				{
					bool print_line = false;
					for (f[0]=0; f[0]<num_field_vals[0]; f[0]++)
					{
						int i;

						for (i=0; i<=4; i++)
							if (f[i] == num_field_vals[i] && num_field_vals[i]>0)
								break;
						if (i<=4)
							continue;

						print_line = true;

						table_entry e;
						for (i=0; i<NUM_TABLE_FIELDS; i++)
							e[i]=f[i];

						int idx = this->calc_table_idx(e);
						os << setw(4) << left << idx;
						os << setw(6) << setprecision(4) << left << score_probs[idx] << " ";
						for (i=0; i< NUM_TABLE_FIELDS; i++)
						{
							if (num_field_vals[i]>0)
								os << setw(lengths[i]) << left << f[i] << " ";
						}
						os << endl;
					}
					if (print_line)
						os << setw(t_lengths) << setfill('-') << "-" << endl << setfill(' ');
				}
}



void FragProbTable::read_table(Config *_config, istream& is)
{
	int i;
	char buff[256];
	istringstream iss;

	config = _config;

	is.getline(buff,256);
	iss.str(buff);

	fields.resize(NUM_TABLE_FIELDS,-1);
	num_field_vals.resize(NUM_TABLE_FIELDS,0);

	// read fields and num_vals
	for (i=0; i<NUM_TABLE_FIELDS; i++)
	{
		iss >> fields[i];
		iss >> num_field_vals[i];
	}

	charge_of_frag = config->get_fragment(fields[0]).charge; // used to correct previous peak offset

	set_field_multipliers();

	is.getline(buff,256);
	iss.str(buff);
	iss >> max_table_size;
	score_probs.resize(max_table_size);
	for (i=0; i<max_table_size; i++)
	{
		char buff[32];
		is.getline(buff,32);
		istringstream iss(buff);
		iss >> score_probs[i];
	}
	
	score_probs_type = 1;
}



void FragProbTable::init_fields(Config *_config, vector<int>& _fields, vector<int>& num_vals)
{
	config = _config;
	fields = _fields;
	num_field_vals = num_vals;
	set_field_multipliers();
	charge_of_frag = config->get_fragment(fields[0]).charge; // used to correct previous peak offset
}


void FragProbTable::init_counts(double init_val)
{
	counts.resize(max_table_size,init_val);
}



void FragProbTable::add_instance(Breakage *breakage, Breakage *previous_breakage,
								 mass_t exp_mass_diff, bool verbose)
{
	table_entry e;

	fill_table_entry(breakage,previous_breakage,e,exp_mass_diff);
	int idx = this->calc_table_idx(e);
	if (verbose)
	{
		string name;
		this->make_table_name(config,name);
		cout << name << endl;
		breakage->print(config);
		cout << endl;
		previous_breakage->print(config);
		cout << endl;
		int i;
		for (i=0; i<5; i++)
			cout << e[i] << " " << "("<<field_multipliers[i]<<")" << "  ";
		cout << " => " << idx << endl << endl;
	}
	counts[idx]++;
}


void FragProbTable::calc_probs()
{
	int i;
	int nv = num_field_vals[0];
	score_probs.resize(max_table_size,0);

	for (i=0; i<max_table_size; i+= nv)
	{
		int j;
		double total=0;
		for (j=0; j<nv; j++)
			total+=counts[i+j];

		for (j=0; j<nv; j++)
			score_probs[i+j] = (score_t)(counts[i+j] / total);
	}
}


double FragProbTable::calc_dkl_sum(const vector<score_t>& ind_probs) const
{
	double dkl_sum =0, total_c=0;
	int num_vals = this->num_field_vals[0];
	int idx;

	if (ind_probs.size() != num_vals)
	{
		cout << "Error: ind_probs not same number as num_vals: " << ind_probs.size() << " vs. "
			<< num_vals << endl;
		exit(1);
	}

	for (idx=0; idx<max_table_size; idx+= num_vals)
	{
		double c=0;
		double dkl=0;
		int i;
		for (i=0; i<num_vals; i++)
		{
			int bin_idx = idx + i;
			dkl += score_probs[bin_idx] * log(score_probs[bin_idx]/ind_probs[i]); 
			c += counts[bin_idx];
		}
		
		dkl_sum += c * dkl;
		total_c += c;
	}
	return dkl_sum / total_c;
}


// changes the values in the score_probs vector from
// probabilities to scores
void FragProbTable::convert_to_score(const vector<double>& rand_probs)
{
	const int num_vals = num_field_vals[0];

	if (score_probs_type == 1)
		return;

	if (rand_probs.size() != num_field_vals[0])
	{
		cout << "Error: number of random probs != number of values for first field: " <<
			rand_probs.size() << " vs. " << num_vals << endl;
		exit(1);
	}

	int i;
	for (i=0; i<score_probs.size(); i+= num_vals)
	{
		int j;
		for (j=0; j<num_vals; j++)
			score_probs[i+j]=log(score_probs[i+j]/rand_probs[j]);
	}

	score_probs_type=1;
}




// puts integer values into the entry based on the values of the 
// tables fragments. If num_field vals == 2, then a binary value is
// given, where 1 represents all values>0
void FragProbTable::fill_table_entry(Breakage *breakage, Breakage *previous_breakage,
						  table_entry& entry, mass_t exp_offset) const
{
	int i;
	

	for (i=0; i< NUM_TABLE_FIELDS; i++)
		entry[i]=0;

	int frag_pos = breakage->get_position_of_frag_idx(fields[0]);
	entry[0] = (frag_pos<0) ? 0 : breakage->fragments[frag_pos].peak_level;

	if (num_field_vals[1]>0 && fields[1]>=0)
	{
		int parent_pos1 = breakage->get_position_of_frag_idx(fields[1]);
		entry[1] = (parent_pos1<0) ? 0 : breakage->fragments[parent_pos1].peak_level;
	}

	if (num_field_vals[2]>0 && fields[2]>=0)
	{
		int parent_pos2 = breakage->get_position_of_frag_idx(fields[2]);
		entry[2] = (parent_pos2<0) ? 0 : breakage->fragments[parent_pos2].peak_level;
	}

	if (previous_breakage)
	{
		if (num_field_vals[3]>0 && fields[3]>=0)
		{
			int previous_pos1 = previous_breakage->get_position_of_frag_idx(fields[3]);
			entry[3] = (previous_pos1<0) ? 0 : 1;

			// change value to one that considers offsets between peaks
			if (num_field_vals[3]>2 && entry[0]>0 && entry[3]>0)
			{
				if (charge_of_frag>1)
					exp_offset /= charge_of_frag;
				
				mass_t dis = (breakage->fragments[frag_pos].mass >
							  previous_breakage->fragments[previous_pos1].mass) ?

								fabs(breakage->fragments[frag_pos].mass - 
								   previous_breakage->fragments[previous_pos1].mass -
								   exp_offset) 
								   :
								fabs(previous_breakage->fragments[previous_pos1].mass -
									breakage->fragments[frag_pos].mass - 
								    exp_offset);
				

				if (dis < 0.4*model_tolerance)
				{
					entry[3]=1;
				}
				else if (dis < model_tolerance)
				{
					entry[3]=2;
				}
				else if (dis < 4*model_tolerance)
				{
					entry[3]=3;
				}
				else 
				{
					cout << "Error: offsets too large for peaks:" << dis << endl;
					exit(1);
				}
			}
		}

		if (num_field_vals[4]>0 && fields[4]>=0)
		{
			int previous_pos2 = previous_breakage->get_position_of_frag_idx(fields[4]);
			entry[4] = (previous_pos2<0) ? 0 : 1;
		}
	}

	// correct the entry vals for binary case
	for (i=0; i<NUM_TABLE_FIELDS; i++)
		if (num_field_vals[i] == 2)
			if (entry[i]>1)
				entry[i]=1;
}
	



void FragProbTable::make_table_name(const Config *config, string& name) const
{
	name="P( ";
	name += config->get_fragment(this->fields[0]).label;
	int i;
	for (i=1; i<fields.size(); i++)
		if (fields[i]>=0)
			break;
	if (i==fields.size())
	{
		name += " )";
		return;
	}

	name += " | ";
	for (i=1; i<3 && i < fields.size(); i++)
		if (fields[i]>=0)
			name+= config->get_fragment(this->fields[i]).label + " ";

	bool put_comma=false;
	for (i=3; i<fields.size(); i++)
		if (fields[i]>=0)
		{
			if (! put_comma)
			{
				put_comma=true;
				name += " , ";
			}
			name+= config->get_fragment(this->fields[i]).label + " ";
		}
	name += ")";
}



score_t FragProbTable::get_score(Breakage *breakage, Breakage *previous_breakage,
								 mass_t exp_mass_diff) const
{
	table_entry entry;
	fill_table_entry(breakage,previous_breakage,entry,exp_mass_diff);
	const int idx  = calc_table_idx(entry);
	return score_probs[idx];
}




