#include "Config.h"


/**************************************************************
Reads a list of PTMS from a file
C    57.021464   FIXED   ALL    C+57    Carbamidomethyl
***************************************************************/
void Config::read_PTM_file(char *file)
{
	if (ind_read_PTM_file)
		return; 

	char buff[256];
	FILE *stream = fopen(file,"r");
	if (! stream)
	{
		cout << "Error: couldn't open PTM file: " << file << endl;
		exit(1);
	}

	all_fixed_PTMs.clear(); 
	all_optional_PTMs.clear();
	all_terminal_PTMs.clear();

	while (fgets(buff,256,stream))
	{
		char aa_label[64];
		char offset_str[64];
		char type_str[64];
		char region_str[64];
		char symbol[64];
		char name[64];

		if (sscanf(buff,"%s %s %s %s %s %s",aa_label,offset_str,type_str,region_str,symbol,name) != 6)
			continue;

		if (aa_label[0]== '#')
			continue;

		int aa = get_aa_from_label(aa_label);

		if (! strcmp("N_TERM",aa_label))
			aa = N_TERM;

		if (! strcmp("C_TERM",aa_label))
			aa = C_TERM;

		if (aa<0)
		{
			cout << "Unkown aa "<< aa_label << "  in PTM:" << endl << buff << endl;
			exit(1);
		}

		mass_t offset = (mass_t)atof(offset_str);
		int type = (! strcmp(type_str,"OPTIONAL")) ? PTM_OPTIONAL : PTM_FIXED;

		int region=-1;
		int position = 0;
		if (! strcmp("ALL",region_str))
		{
			region = PTM_ALL;
		}
		else if (! strcmp("C_TERM",region_str))
		{
			region = PTM_C_TERMINAL;
		}
		else if (! strcmp("N_TERM",region_str))
		{
			region = PTM_N_TERMINAL;
		}
		else
		{
			position = atoi(region_str);
			if (position<10 && position>-10)
				region = PTM_POSITION;
		}

		if (region<0)
		{
			cout << "Error: bad PTM region : " << region_str << endl;
			exit(1);
		}
		
		// add PTMs
		PTM ptm = PTM(aa, offset , type , region ,position, symbol, name);
		if (region == PTM_C_TERMINAL || region == PTM_N_TERMINAL)
		{
			ptm.position = (region == PTM_N_TERMINAL ? 1 : -1);
			if (type == PTM_FIXED)
			{
				cout << "Error: terminal PTMs must be of type OPTIONAL!" << endl;
				exit(1);
			}
			all_terminal_PTMs.add_PTM(ptm);
		}
		else if (type == PTM_OPTIONAL)
		{
			all_optional_PTMs.add_PTM(ptm);
		}
		else if (type == PTM_FIXED)
		{
			all_fixed_PTMs.add_PTM(ptm);
		}

	
	}

	// sanity check (no two PTMs can have the same label)
	vector<string> used_labels;

	int i;
	for (i=0; i<all_fixed_PTMs.get_num_PTMs(); i++)
	{
		int j;
		for (j=0; j<used_labels.size(); j++)
		{
			if (used_labels[j] == all_fixed_PTMs.list[i].label)
			{
				cout << "Error: " << used_labels[j] << " is already used!" << endl;
				exit(1);
			}
		}
		used_labels.push_back(all_fixed_PTMs.list[i].label);
	}

	for (i=0; i<all_optional_PTMs.get_num_PTMs(); i++)
	{
		int j;
		for (j=0; j<used_labels.size(); j++)
		{
			if (used_labels[j] == all_optional_PTMs.list[i].label)
			{
				cout << "Error: " << used_labels[j] << " is already used!" << endl;
				exit(1);
			}
		}
		used_labels.push_back(all_optional_PTMs.list[i].label);
	}

	for (i=0; i<all_terminal_PTMs.get_num_PTMs(); i++)
	{
		int j;
		const PTM& ptm = all_terminal_PTMs.list[i];

		if (ptm.region == PTM_ALL)
		{
			cout << "Error:  region must be N_TERM or C_TERM for terminal mod!" 
				<< endl << ptm << endl;
			exit(1);
		}

		if (ptm.org_aa != N_TERM && ptm.org_aa != C_TERM)
		{
			cout << "Error: original aa must be N_TERM or C_TERM for terminal modification!" 
				<< endl << ptm << endl;
			exit(1);
		}

		for (j=0; j<used_labels.size(); j++)
		{
			if (used_labels[j] == ptm.label)
			{
				cout << "Error: " << used_labels[j] << " is already used!" << endl;
				exit(1);
			}
		}
		used_labels.push_back(ptm.label);
	}

	// update the conversion vectors according to the PTMS

	session_tables.init_for_standard_aas();

	// add the optional AAs
	for (i=0; i<all_optional_PTMs.get_num_PTMs(); i++)
	{
		const PTM& ptm = all_optional_PTMs.list[i];
		session_tables.add_optional_PTM_aa(ptm.org_aa,ptm.label,ptm.delta,ptm.position);
	}

	// add optional terminal AAs
	for (i=0; i<all_terminal_PTMs.get_num_PTMs(); i++)
	{
		const PTM& ptm = all_terminal_PTMs.list[i];

		if (ptm.type != PTM_OPTIONAL)
			continue;

		session_tables.add_optional_PTM_terminal_aa(ptm.delta,ptm.position,ptm.label);
	}


	// create a mapping between label and aa
	const vector<string>& aa2label = session_tables.get_aa2label();
	for (i=0; i<aa2label.size(); i++)
	{
		label2aa.insert(STRING2INT_MAP::value_type(aa2label[i],i));
	}

	ind_read_PTM_file = true;
	
//	this->print_supported_PTMs();
}



/**************************************************************
// returns the idx of the label
// -1 is returned if label is not found
***************************************************************/
int PTM_list::get_PTM_idx(const string& label) const
{
	int i;

	for (i=0; i<list.size(); i++)
	{
		const string & l_label = list[i].label;
		if (list[i].label == label)
			return i;
	}

	return -1;
}


/**************************************************************
// adds a PTM, makes sure that it wasn't already added
// (if already added, returns false
***************************************************************/
bool PTM_list::add_PTM(const PTM& ptm)
{
	int i;
	for (i=0; i<list.size(); i++)
	{
		if (ptm.label == list[i].label)
			return false;
	}

	if (ptm.label.length()<2)
	{
//		cout << "Error: PTM label must be longer (form should be X-dd or X+dd)\n";
//		exit(1);
	}
	list.push_back(ptm);
	return true;
}







/********************************************************************
Gets a text line with all the PTM labels that are to be used in this
run. Updates the conversion tables and PTM lists appropriately
Only fixed mods cause a change to the masses (the optional mod are
assumed to already have the correction).
This function erases the effects of any previous PTMs selected!
*********************************************************************/
void Config::apply_selected_PTMs(const char *ptm_line)
{
	// replace all colons with white space

	int ptm_line_length = strlen(ptm_line);
	char *new_ptm_line = new char[strlen(ptm_line)+1];

	strcpy(new_ptm_line,ptm_line);

	int i;
	for (i=0; i<ptm_line_length; i++)
		if (new_ptm_line[i] == ':')
			new_ptm_line[i] =  ' ';

	istringstream ptm_stream(new_ptm_line);
	string ptm_label;

	// create lists of AAs for different regions
	session_aas = standard_aas;
	selected_optional_PTM_idxs.clear();
	selected_fixed_PTM_idxs.clear();
	selected_terminal_PTM_idxs.clear();
	
	int num_fixed_C_TERM_PTMs = 0 , num_fixed_N_TERM_PTMs = 0;

	while (ptm_stream >> ptm_label)
	{
		int idx;

		idx=all_optional_PTMs.get_PTM_idx(ptm_label);
		if (idx>=0)
		{
			selected_optional_PTM_idxs.push_back(idx);
			continue;
		}

		idx=all_fixed_PTMs.get_PTM_idx(ptm_label);
		if (idx>=0)
		{
			const PTM& ptm = all_fixed_PTMs.list[idx];
			selected_fixed_PTM_idxs.push_back(idx);

			if (ptm.region == PTM_ALL)
			{
				if (! session_tables.make_fixed_mod(ptm.org_aa,ptm.delta))
				{
					cout << "Warning: could not make fixed PTM: " << ptm.label << " (mass shift too negative)" << endl;
				}
				continue;
			}

			cout << "ERROR: Fixed PTMs must be in region PTM_ALL!" << endl;
			exit(1);
		}  

		idx=all_terminal_PTMs.get_PTM_idx(ptm_label);
		if (idx>=0)
		{
			const PTM& ptm = all_terminal_PTMs.list[idx];
			selected_terminal_PTM_idxs.push_back(idx);

			if (ptm.type != PTM_FIXED)
				continue;


			cout << "Error: terminal modifications must be of type OPTIONAL!" << endl;
			exit(1);


			if (ptm.region == PTM_N_TERMINAL)
			{
				if (num_fixed_N_TERM_PTMs>0)
				{
					cout << "Error: only one fixed terminal PTM allowed for N terminal!"
						<< endl << ptm << endl;
					exit(1);
				}
				session_tables.make_fixed_mod(N_TERM,ptm.delta);
				num_fixed_N_TERM_PTMs++;
				continue;
			}

			if (ptm.region == PTM_C_TERMINAL)
			{
				if (num_fixed_C_TERM_PTMs>0)
				{
					cout << "Error: only one fixed terminal PTM allowed for C terminal!"
						<< endl << ptm << endl;
					exit(1);
				}
				session_tables.make_fixed_mod(C_TERM,ptm.delta);
				num_fixed_C_TERM_PTMs++;
				continue;
			}

			cout << "Error : fixed terminal PTM must have a terminal region!"
				<< endl << ptm << endl;
			exit(1);
		}

		cout << "Error: No support for PTM: " << ptm_label << endl;
		exit(1);
	}


	// add optional PTMs
	for (i=0; i<selected_optional_PTM_idxs.size(); i++)
	{
		const PTM& ptm = all_optional_PTMs.list[selected_optional_PTM_idxs[i]];

		int aa_idx = get_aa_from_label(ptm.label);
		if (aa_idx<0)
		{
			cout << "Error: unknown PTM label " << ptm.label << endl;
			exit(1);
		}


		if (ptm.region == PTM_ALL || ptm.region == PTM_POSITION)
		{
			session_aas.push_back(aa_idx);
		}
	}

	// add optional terminal PTMs
	for (i=0; i<selected_terminal_PTM_idxs.size(); i++)
	{
		const PTM& ptm = all_terminal_PTMs.list[selected_terminal_PTM_idxs[i]];

		if (ptm.type != PTM_OPTIONAL)
			continue;

		if (ptm.region == PTM_N_TERMINAL)
		{
			// create new amino acids that have the optional terminal PTM
			// but can only be applied to the position +1
			const vector<char>& aa2char = get_aa2char();
			int aa;
			for (aa=Ala; aa<=Val; aa++)
			{
				string aa_label = aa2char[aa] + ptm.label;
				int aa_idx = get_aa_with_position_from_label(aa_label,1);
				if (aa_idx<0)
				{
					cout << "Error: unknown PTM label " << aa_label << endl;
					exit(1);
				}
				session_aas.push_back(aa_idx);
			}
			continue;
		}

		if (ptm.region == PTM_C_TERMINAL)
		{
				// create new amino acids that have the optional terminal PTM
			// but can only be applied to the position -1
			const vector<char>& aa2char = get_aa2char();
			int aa;
			for (aa=Ala; aa<=Val; aa++)
			{
				string aa_label = aa2char[aa] + ptm.label;
				int aa_idx = get_aa_with_position_from_label(aa_label,-1);
				if (aa_idx<0)
				{
					cout << "Error: unknown PTM label " << aa_label << endl;
					exit(1);
				}
				session_aas.push_back(aa_idx);
			}
			continue;
		}

		cout << "Error: bad PTM region for terminal PTM : " << ptm.region << endl;
		exit(1);
	}

	sort(session_aas.begin(),session_aas.end());

	calc_aa_combo_masses();
	set_aa_variants();

	// set the maximal n and c terminal PTM values
	const vector<int>& org_aa = this->get_org_aa();
	const vector<mass_t>& aa2mass = this->get_aa2mass();

	max_n_term_mod = aa2mass[N_TERM];
	max_c_term_mod = aa2mass[C_TERM];
	min_n_term_mod = aa2mass[N_TERM];
	min_c_term_mod = aa2mass[C_TERM];

	for (i=0; i<session_aas.size(); i++)
	{
		const int& aa = session_aas[i];

		if (org_aa[aa] == N_TERM)
		{
			if (aa2mass[aa]>max_n_term_mod)
				max_n_term_mod = aa2mass[aa];
			if (aa2mass[aa]<min_n_term_mod)
				min_n_term_mod = aa2mass[aa];
		}
		else if (org_aa[aa] == C_TERM)
		{
			if (aa2mass[aa]>max_c_term_mod)
				max_c_term_mod = aa2mass[aa];
			if (aa2mass[aa]<min_c_term_mod)
				min_c_term_mod = aa2mass[aa];
		} 
	}

	init_allowed_node_masses(400);

//	this->print_aa_variants();

	delete [] new_ptm_line;
}


void Config::apply_site_input_PTMs(const vector<string>& ptm_lines)
{
	all_fixed_PTMs.clear();    // region PTM_ALL
	all_optional_PTMs.clear(); // region PTM_ALL
	all_terminal_PTMs.clear(); // must be either PTM_N_TERMINAL, PTM_C_TERMINAL

	session_tables.init_for_standard_aas();

	int next_aa_idx = Val+1;

	{
		int i;
		const vector<string>& aa2label = session_tables.get_aa2label();
		label2aa.clear();
		for (i=0; i<aa2label.size(); i++)
			label2aa.insert(STRING2INT_MAP::value_type(aa2label[i],i));
	}


	bool read_valid_PTMs = false;
	
	// parse the PTM lines, create PTMs for each line
	
	int i;
	for (i=0; i<ptm_lines.size(); i++)
	{
		char offset_str[32];
		char aa_labels_str[32];
		char type_str[32];
		char name_str[32];

		vector<int> ptm_aas;

		ptm_aas.clear();

	//	mod,[MASS],[RESIDUES],[TYPE],[NAME]
		string line = ptm_lines[i];
		int j;
		for (j=0; j<line.length(); j++)
			if (line[j]==',')
				line[j]='\t';

		int num_parameters = sscanf(line.c_str(),"mod\t%s\t%s\t%s\t%s",
											offset_str,aa_labels_str,type_str,name_str);
		if ( num_parameters<2)
		{
			cout << "Warning: couldn't parse PTM line: " << ptm_lines[i] << " ... ignoring. " << endl;
			continue;
		}

		if (! strcmp(aa_labels_str,"*") )
		{
			ptm_aas = get_session_aas();
		}
		else
		{
			int j;
			for (j=0; j<strlen(aa_labels_str); j++)
			{
				char aa_buff[2];
				aa_buff[0]=aa_labels_str[j];
				aa_buff[1]='\0';
				int aa = get_aa_from_label(aa_buff);
				if (aa>=0)
				{
					ptm_aas.push_back(aa);
				}
				else
				{
					cout << "Warning: Bad aa letter: " << aa_labels_str[j] << " in PTM line: " <<
						ptm_lines[i] << " ...ignoring." << endl;
					continue;
				}
			}
		}

		if (ptm_aas.size()==0)
			continue;


		mass_t offset = (mass_t)atof(offset_str);
		int type     = PTM_OPTIONAL;
		int region   = PTM_ALL;
		string name  = "";
		

		if (num_parameters>2)
		{
			if (! strcmp(type_str,"opt"))
			{
				type = PTM_OPTIONAL;	
			}
			else if (! strcmp(type_str,"fix"))
			{
				type = PTM_FIXED;
			}
			else if (! strcmp(type_str,"nterminal"))
			{
				type = PTM_OPTIONAL;
				region = PTM_N_TERMINAL;
			}
			else if (! strcmp(type_str,"cterminal"))
			{
				type = PTM_OPTIONAL;
				region = PTM_C_TERMINAL;
			}
		}

		// add PTMs
		
		for (j=0; j<ptm_aas.size(); j++)
		{
			// make label
			char symbol_str[32];
			char plus_minus = '+';

			if (offset<0)
				plus_minus='-';

			int integer_offset = (int)((offset<0) ? ceil(offset-0.5) : floor(offset+0.5));
			if (integer_offset<0)
				integer_offset *= -1;
		
			sprintf(symbol_str,"%c%c%d",this->get_aa2char()[ptm_aas[j]],plus_minus,integer_offset);

			name =(num_parameters>3) ? name_str : symbol_str;
			
			PTM ptm = PTM(ptm_aas[j], offset , type , region , 0, symbol_str, name);

			if (region == PTM_C_TERMINAL || region == PTM_N_TERMINAL)
			{
				ptm.position = (region == PTM_N_TERMINAL ? 1 : -1);
				if (type == PTM_FIXED)
				{
					cout << "Error: terminal PTMs must be of type OPTIONAL!" << endl;
					exit(1);
				}
				all_terminal_PTMs.add_PTM(ptm);
			}
			else if (type == PTM_OPTIONAL)
			{
				all_optional_PTMs.add_PTM(ptm);
			}
			else if (type == PTM_FIXED)
			{
				all_fixed_PTMs.add_PTM(ptm);
				if (! session_tables.make_fixed_mod(ptm.org_aa,ptm.delta))
				{
					cout << "Warning: could not make fixed PTM: " << ptm.label << " (mass shift too negative)" << endl;
				}
			}
		}
	}


	// Update the tables according to these PTMs


	// add the optional AAs
	for (i=0; i<all_optional_PTMs.get_num_PTMs(); i++)
	{
		const PTM& ptm = all_optional_PTMs.list[i];
		session_tables.add_optional_PTM_aa(ptm.org_aa,ptm.label,ptm.delta,ptm.position);
		label2aa.insert(STRING2INT_MAP::value_type(ptm.label,next_aa_idx++));

		int aa_idx = get_aa_from_label(ptm.label);
		if (aa_idx<0)
		{
			cout << "Error: unknown PTM label " << ptm.label << endl;
			exit(1);
		}

		if (ptm.region == PTM_ALL || ptm.region == PTM_POSITION)
		{
			session_aas.push_back(aa_idx);
			
		}
	}

	// add optional terminal AAs
	for (i=0; i<all_terminal_PTMs.get_num_PTMs(); i++)
	{
		const PTM& ptm = all_terminal_PTMs.list[i];

		if (ptm.type != PTM_OPTIONAL)
			continue;

		session_tables.add_optional_PTM_aa(ptm.org_aa,ptm.label,ptm.delta,ptm.position);

		label2aa.insert(STRING2INT_MAP::value_type(ptm.label,next_aa_idx++));

		int aa_idx = get_aa_with_position_from_label(ptm.label,ptm.position);
		if (aa_idx<0)
		{
			cout << "Error: unknown PTM label " << ptm.label << endl;
			exit(1);
		}
		session_aas.push_back(aa_idx);
	}

	// update combos
	sort(session_aas.begin(),session_aas.end());
	calc_aa_combo_masses();
	set_aa_variants();

	// set the maximal n and c terminal PTM values
	const vector<int>& org_aa = this->get_org_aa();
	const vector<mass_t>& aa2mass = this->get_aa2mass();

	max_n_term_mod = aa2mass[N_TERM];
	max_c_term_mod = aa2mass[C_TERM];
	min_n_term_mod = aa2mass[N_TERM];
	min_c_term_mod = aa2mass[C_TERM];

	for (i=0; i<session_aas.size(); i++)
	{
		const int& aa = session_aas[i];

		if (org_aa[aa] == N_TERM)
		{
			if (aa2mass[aa]>max_n_term_mod)
				max_n_term_mod = aa2mass[aa];
			if (aa2mass[aa]<min_n_term_mod)
				min_n_term_mod = aa2mass[aa];
		}
		else if (org_aa[aa] == C_TERM)
		{
			if (aa2mass[aa]>max_c_term_mod)
				max_c_term_mod = aa2mass[aa];
			if (aa2mass[aa]<min_c_term_mod)
				min_c_term_mod = aa2mass[aa];
		} 
	}

	init_allowed_node_masses(400);


	// create a mapping between label and aa
	const vector<string>& aa2label = session_tables.get_aa2label();
	for (i=0; i<aa2label.size(); i++)
	{
		label2aa.insert(STRING2INT_MAP::value_type(aa2label[i],i));
	}


	ind_read_PTM_file = true;
}




ostream& operator << (ostream& os, const PTM& ptm)
{
	os << setw(8) << ptm.label << " "<< setw(12) << fixed << right << ptm.delta << "  ";
	os << setw(12) << left;
	if (ptm.region == PTM_ALL)
	{
		os << "ALL";
	}
	else if (ptm.region == PTM_N_TERMINAL)
	{
		os << "N_TERM";
	}
	else if (ptm.region == PTM_C_TERMINAL)
	{
		os << "C_TERM";
	}
	else if (ptm.region == PTM_POSITION)
	{
		os << ptm.position;
	}

	os << setw(10) << left << ((ptm.type == PTM_OPTIONAL) ? "OPTIONAL " : "FIXED ");

	if (! ptm.name.empty())
		os<< " (" << ptm.name << ")";

	return os;
}


void Config::print_supported_PTMs() const
{
	int i;

	cout << endl << "FIXED MODIFICATIONS ( " << all_fixed_PTMs.list.size() <<" )" << endl;
	for (i=0; i<all_fixed_PTMs.list.size(); i++)
		cout << left << setw(3) << i+1 << " " << all_fixed_PTMs.list[i] << endl;


	cout << endl <<"OPTIONAL MODIFICATIONS ( " << all_optional_PTMs.list.size() <<" )" << endl;
	for (i=0; i<all_optional_PTMs.list.size(); i++)
		cout <<  left << setw(3) << i+1 << " " << all_optional_PTMs.list[i] << endl;

	cout << endl << "TERMINAL MODIFICATIONS ( " << all_terminal_PTMs.list.size() <<" )" << endl;
	for (i=0; i<all_terminal_PTMs.list.size(); i++)
		cout <<  left << setw(3) << i+1 << " " << all_terminal_PTMs.list[i] << endl;

	cout << endl << endl; 
}

