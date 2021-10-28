#include "PeptideRankScorer.h"
#include "AllScoreModels.h"
#include "DeNovoSolutions.h"
#include "PepNovo_auxfun.h"


bool add_peptide_set_to_psm(PeptideSetMap& psm, PeptideSet& pep_set)
{
	PeptideSetMap::iterator it;

	scan_pair key(pep_set.file_idx,pep_set.scan);

	it=psm.find(key);
	if (it == psm.end())
	{
		psm.insert(pair<scan_pair,PeptideSet>(key,pep_set));
		return false;
	}
	
	PeptideSet& set_in_map = (*it).second;
	const int corr_aa_before = set_in_map.correct_sol.pep.get_aa_before();
	const int corr_aa_after  = set_in_map.correct_sol.pep.get_aa_after();
	int i;
	for (i=0; i<pep_set.incorrect_sols.size(); i++)
	{
		PeptideSolution& sol = pep_set.incorrect_sols[i];
	
		// check if sol is already there
		int j;
		for (j=0; j<set_in_map.incorrect_sols.size(); j++)
			if (sol.pep == set_in_map.incorrect_sols[j].pep)
				break;

	//	if (j<set_in_map.incorrect_sols.size())
	//		continue;
		
		if (sol.pep.get_aa_before()<0 && sol.pep.get_aa_after()<0)
		{
			sol.pep.set_aa_before(corr_aa_before);
			sol.pep.set_aa_after(corr_aa_after);
		}
		set_in_map.incorrect_sols.push_back(sol);
	}
	
	return true;
}

void read_correct_seqs(const string& path, 
					   vector<string>& seqs, 
					   vector<char>& n_aas,
					   vector<char>& c_aas,
					   vector<float>& mqscores)
{
	ifstream ifs(path.c_str());
	if (! ifs.is_open())
	{
		cout << "Error: couldn't find file with correct sequences: " << path << endl;
		exit(1);
	}

	seqs.clear();
	n_aas.clear();
	c_aas.clear();

	while (! ifs.eof())
	{
		char buff[128];
		string pep_seq;
		char n_aa='$',c_aa='$';
		float mqscore=NEG_INF;
		ifs.getline(buff,128);
		if (ifs.gcount()<5)
			break;

		istringstream iss(buff);
		iss >> pep_seq >> n_aa >> c_aa >> mqscore;

		if (n_aa == '$' || c_aa == '$' || pep_seq.length()<4 || pep_seq.length()>100)
		{
			cout << "Error: bad line:|" << buff << "|" << endl;
			exit(1);
		}
		
		seqs.push_back(pep_seq);
		n_aas.push_back(n_aa);
		c_aas.push_back(c_aa);
		mqscores.push_back(mqscore);
	}
	ifs.close();
}

/************************************************************************
*************************************************************************/
void add_db_hits_to_psm(Config *config,
						const vector<string>& names,
						const string& db_dir,
						const string& correct_seq_dir,
						const string& file_suffix,
						int   file_charge,
						PeptideSetMap& psm)
{
	const vector<int>& char2aa = config->get_char2aa();
	int num_read=0;
	int total_num_sets=0;
	int total_num_peptides=0;
	
	int file_idx;
	for (file_idx=0; file_idx<names.size(); file_idx++)
	{
		string db_file = db_dir + "/" + names[file_idx] + file_suffix;
		ifstream ifs(db_file.c_str());

		if (! ifs.is_open())
			continue;

		cout << file_idx << "\t" << db_file << endl;

		num_read++;
		int num_sets=0;
		int num_peptides=0;
		int num_with_correct=0;
		int num_first_correct=0;

		vector<string> correct_seqs;
		vector<char>   correct_n_aas, correct_c_aas;
		vector<float>  mqscores;

		string seq_path = correct_seq_dir + "/" + names[file_idx] + ".txt";
		read_correct_seqs(seq_path,correct_seqs, correct_n_aas, correct_c_aas, mqscores);
		
		while (! ifs.eof())
		{
			char buff[2048];
			ifs.getline(buff,2048);
			if (ifs.gcount()>2046)
				cout << "Warning: buffer size?" << endl;

			if (ifs.gcount()<5)
				continue;
			
			istringstream iss(buff);
			int scan=-1;
			mass_t mass_with_19 = -1;
			int num_seqs = -1;

			iss >> scan >> mass_with_19 >> num_seqs;

			if (scan<0 || 
				mass_with_19<0 || 
				mass_with_19>10000 || 
				num_seqs<0 || 
				num_seqs>100)
			{
				cout << scan << "\t" <<mass_with_19 << "\t" << num_seqs << endl;
				cout << "Error parsing line in file " << db_file << endl << "LINE:" <<buff << endl;
				exit(1);
			}

			if (scan>correct_seqs.size())
			{
				cout << "Error: mismatch in correct seqs file!" << endl;
				exit(1);
			}

			const string& correct_pep_str = correct_seqs[scan];
			bool has_correct = false;
			bool correct_first = false;

			PeptideSet set;
			set.file_idx=file_idx;
			set.scan = scan;
			set.correct_sol.charge = file_charge;
			set.correct_sol.type = SOL_CORRECT;
			set.correct_sol.pep.parseFromString(config,correct_pep_str);
			set.correct_sol.pep.set_aa_before(char2aa[correct_n_aas[scan]]);
			set.correct_sol.pep.set_aa_after(char2aa[correct_c_aas[scan]]);
			set.correct_sol.num_correct_aas = set.correct_sol.pep.get_num_aas();
			set.correct_sol.pm_with_19 = set.correct_sol.pep.get_mass_with_19();
			set.correct_sol.reaches_n_terminal = true;
			set.correct_sol.reaches_c_terminal = true;
			set.correct_sol.MQScore= mqscores[scan];
		
			vector<string> inserted_strings;
			inserted_strings.clear();
			int i;
			for (i=0; i<num_seqs; i++)
			{
				char char_before, char_after;
				int charge=0;
				float mqscore;
				string peptide_str;
				iss >> char_before >> char_after >> charge >> mqscore >> peptide_str;

				const int aa_before = char2aa[char_before];
				const int aa_after  = char2aa[char_after];

				if (charge<=0 || charge>5)
				{
					cout << "Error: bad charge: " << charge << " in line:" << endl << buff << endl;
					exit(1);
				}


				int j;
				for (j=0; j<peptide_str.length(); j++)
					if (peptide_str[j]=='I')
						peptide_str[j]='L' ;

				// check if peptide is same as correct don't add it
				if (peptide_str == correct_pep_str)
				{
					has_correct   = true;
					if (i==0)
						correct_first = true;
					continue;
				}

				for (j=0; j<inserted_strings.size(); j++)
					if (inserted_strings[j] == peptide_str)
						break;
				if (j<inserted_strings.size())
					continue;

				inserted_strings.push_back(peptide_str);

				PeptideSolution sol;
				sol.charge = charge;
				sol.pep.parseFromString(config,peptide_str);
				sol.pep.set_aa_before(aa_before);
				sol.pep.set_aa_after(aa_after);
				sol.MQScore = mqscore;
				sol.num_correct_aas = set.correct_sol.pep.calc_number_of_correct_aas(config,sol.pep);
				sol.type= SOL_INCORRECT_DB;
				sol.pm_with_19 = sol.pep.get_mass_with_19();
				sol.reaches_c_terminal=true;
				sol.reaches_n_terminal=true;
				set.incorrect_sols.push_back(sol);
				num_peptides++;
			}

			if (add_peptide_set_to_psm(psm,set))
				total_num_sets++;

			num_sets++;
			if (has_correct)
				num_with_correct++;
			if (correct_first)
				num_first_correct++;
		}
		ifs.close();

		cout << "\tsets added:   " << num_sets << endl;
		cout << "\thad correct:  " << float(num_with_correct)/num_sets << endl;
		cout << "\tcorrect first:" << float(num_first_correct)/num_sets << endl;
		cout << endl;
	}

	cout << "Processed " << num_read << "/" << names.size() << " dbh files." << endl;
}


void add_denovo_dbh_to_psm(Config *config,
						const vector<string>& names,
						const string& denovo_dir,
						const string& correct_seq_dir,
						const string& file_suffix,
						int file_charge,
						PeptideSetMap& psm)
{
	const vector<int>& char2aa = config->get_char2aa();
	int num_read=0;
	int total_num_peptides=0;
	
	int file_idx;
	for (file_idx=0; file_idx<names.size(); file_idx++)
	{
		string dnv_file = denovo_dir + "/" + names[file_idx] + file_suffix;
		ifstream ifs(dnv_file.c_str());

		if (! ifs.is_open())
			continue;

		cout << file_idx << "\t" << dnv_file << endl;

		num_read++;
		int num_new_sets=0;
		int num_sets=0;
		int num_with_correct=0;
		int num_first_correct=0;

		vector<string> correct_seqs;
		vector<char>   correct_n_aas, correct_c_aas;
		vector<float>  mqscores;

		string seq_path = correct_seq_dir + "/" + names[file_idx] + ".txt";
		read_correct_seqs(seq_path,correct_seqs, correct_n_aas, correct_c_aas, mqscores);
		
		while (! ifs.eof())
		{
			char buff[2048];
			ifs.getline(buff,2048);
			if (ifs.gcount()>2046)
				cout << "Warning: buffer size?" << endl;

			if (ifs.gcount()<5)
				continue;
			
			istringstream iss(buff);
			int scan=-1;
			mass_t mass_with_19 = -1;
			int num_seqs = -1;

			iss >> scan >> mass_with_19 >> num_seqs;

			if (scan<0 || 
				mass_with_19<0 || 
				mass_with_19>10000 || 
				num_seqs<0 || 
				num_seqs>100)
			{
				cout << scan << "\t" <<mass_with_19 << "\t" << num_seqs << endl;
				cout << "Error parsing line in file " << dnv_file << endl << "LINE:" <<buff << endl;
				exit(1);
			}

			if (scan>correct_seqs.size())
			{
				cout << "Error: mismatch in correct seqs file!" << endl;
				exit(1);
			}

			const string& correct_pep_str = correct_seqs[scan];
			bool has_correct = false;
			bool correct_first = false;

			PeptideSet set;
			set.file_idx=file_idx;
			set.scan = scan;
			set.correct_sol.charge = file_charge;
			set.correct_sol.type = SOL_CORRECT;
			set.correct_sol.pep.parseFromString(config,correct_pep_str);
			set.correct_sol.pep.set_aa_before(char2aa[correct_n_aas[scan]]);
			set.correct_sol.pep.set_aa_after(char2aa[correct_c_aas[scan]]);
			set.correct_sol.num_correct_aas = set.correct_sol.pep.get_num_aas();
			set.correct_sol.pm_with_19 = set.correct_sol.pep.get_mass_with_19();
			set.correct_sol.reaches_n_terminal = true;
			set.correct_sol.reaches_c_terminal = true;
			set.correct_sol.MQScore= mqscores[scan];
		
			vector<string> inserted_strings;
			inserted_strings.clear();
			int i;
			for (i=0; i<num_seqs; i++)
			{
				char char_before, char_after;
				int charge=0;
				float mqscore;
				string peptide_str;
				iss >> char_before >> char_after >> charge >> mqscore >> peptide_str;

				const int aa_before = char2aa[char_before];
				const int aa_after  = char2aa[char_after];

				if (charge<=0 || charge>5)
				{
					cout << "Error: bad charge: " << charge << " in line:" << endl << buff << endl;
					exit(1);
				}


				int j;
				for (j=0; j<peptide_str.length(); j++)
					if (peptide_str[j]=='I')
						peptide_str[j]='L' ;

				// check if peptide is same as correct don't add it
				if (peptide_str == correct_pep_str)
				{
					has_correct   = true;
					if (i==0)
						correct_first = true;
					continue;
				}

				for (j=0; j<inserted_strings.size(); j++)
					if (inserted_strings[j] == peptide_str)
						break;
				if (j<inserted_strings.size())
					continue;

				inserted_strings.push_back(peptide_str);

				PeptideSolution sol;
				sol.charge = charge;
				sol.pep.parseFromString(config,peptide_str);
				sol.pep.set_aa_before(aa_before);
				sol.pep.set_aa_after(aa_after);
				sol.MQScore = mqscore;
				sol.num_correct_aas = set.correct_sol.pep.calc_number_of_correct_aas(config,sol.pep);
				sol.type= SOL_INCORRECT_DENOVO;
				sol.pm_with_19 = sol.pep.get_mass_with_19();
				sol.reaches_c_terminal=true;
				sol.reaches_n_terminal=true;
				set.incorrect_sols.push_back(sol);
			}

			if (! add_peptide_set_to_psm(psm,set))
				num_new_sets++;

			num_sets++;
			if (has_correct)
				num_with_correct++;
			if (correct_first)
				num_first_correct++;
		}
		ifs.close();

		cout << "\tnew sets added : " << num_new_sets << endl;
		cout << "\tsets added :     " << num_sets << endl;
		cout << "\thad correct :    " << float(num_with_correct)/num_sets << endl;
		cout << "\tcorrect first :  " << float(num_first_correct)/num_sets << endl;
		cout << endl;
	}

	cout << "Processed " << num_read << "/" << names.size() << " dbh files." << endl;
}



void add_denovo_paths_to_psm(Config *config,
							const vector<string>& names,
							const string& denovo_dir,
							const string& file_suffix,
							int   file_charge,
							PeptideSetMap& psm)
{
	const vector<int>& char2aa = config->get_char2aa();
	int num_read=0;
	int total_num_sets=0;
	int total_num_peptides=0;
	
	int file_idx;
	for (file_idx=0; file_idx<names.size(); file_idx++)
	{
		string denovo_file = denovo_dir + "/" + names[file_idx] + file_suffix;
		ifstream ifs(denovo_file.c_str());

		if (! ifs.is_open())
			continue;

		cout << file_idx << "\t" << denovo_file << endl;

		num_read++;
		int num_sets=0;
		int num_with_correct=0;
		int num_first_correct=0;

		while (! ifs.eof())
		{
			char buff[2048];
			ifs.getline(buff,2048);
			if (ifs.gcount()>2046)
				cout << "Warning: buffer size?" << endl;
			
			if (ifs.gcount()<5)
				continue;

			istringstream iss(buff);

		//	>> 92 0 6_6260 EGSSLLGSDAGELAGAGK 1618.79
			string dummy;
			int dummy_file;
			int scan=-1;
			string title;
			string correct_pep_str;
			mass_t mass_with_19 = -1;
		
			iss >> dummy >> dummy_file >> scan >> title >> correct_pep_str >> mass_with_19;
			if (scan<0 || mass_with_19<0 || mass_with_19>10000)
			{
				cout << scan << "\t" << mass_with_19 << endl;
				cout << "Error parsing line in file " << denovo_file << endl << "LINE:" <<buff << endl;
				exit(1);
			}

			int num_seqs=0;
			ifs.getline(buff,2048);
			sscanf(buff,"%d",&num_seqs);

			if (num_seqs<=0 || num_seqs>500)
			{
				cout << "Error parsing line in file " << denovo_file << endl << "LINE:" <<buff << endl;
				exit(1);
			}

		
			bool has_correct = false;
			bool correct_first = false;

			PeptideSet set;
			set.file_idx=file_idx;
			set.scan = scan;
			set.correct_sol.charge = file_charge;
			set.correct_sol.type = SOL_CORRECT;
			set.correct_sol.pep.parseFromString(config,correct_pep_str);
			set.correct_sol.pep.set_aa_before(-1);
			set.correct_sol.pep.set_aa_after(-1);
			set.correct_sol.num_correct_aas = set.correct_sol.pep.get_num_aas();
			set.correct_sol.pm_with_19 = set.correct_sol.pep.get_mass_with_19();
			set.correct_sol.reaches_n_terminal = true;
			set.correct_sol.reaches_c_terminal = true;

			vector<string> inserted_strings;
			inserted_strings.clear();

			int i;
			for (i=0; i<num_seqs; i++)
			{
				ifs.getline(buff,2048);
				mass_t start_mass=-1;
				int num_correct_aas =-1;
				int num_aas = -1;
				string denovo_seq="";

				istringstream iss(buff);
				iss >> start_mass >> num_correct_aas >> num_aas >> denovo_seq;

				if (num_correct_aas==0 && num_aas==0)
					continue;

				if (start_mass<0 || start_mass > 10000 || num_correct_aas<0 ||
					num_aas<0 || num_aas>100 || num_correct_aas> num_aas || 
					denovo_seq.length()<3)
				{
					cout << "Error parsing " << denovo_file << endl;
					cout << "LINE: " << buff << endl;
					exit(1);
				}

				// check if peptide is same as correct don't add it
				if (denovo_seq == correct_pep_str || num_correct_aas == num_aas)
				{
					has_correct   = true;
					if (i==0)
						correct_first = true;
					continue;
				}

				int j;
				for (j=0; j<inserted_strings.size(); j++)
					if (inserted_strings[j] == denovo_seq)
						break;
				if (j<inserted_strings.size())
					continue;

				// check how similar the cuts are

				inserted_strings.push_back(denovo_seq);

				PeptideSolution sol;
				sol.charge = file_charge;
				sol.pep.parseFromString(config,denovo_seq);
				sol.pep.set_aa_before(-1);
				sol.pep.set_aa_after(-1);
				sol.num_correct_aas = num_correct_aas;
				sol.type=SOL_INCORRECT_DENOVO;
				sol.pm_with_19 = sol.pep.get_mass_with_19();
				sol.reaches_n_terminal = true;
				sol.reaches_c_terminal = true;

				set.incorrect_sols.push_back(sol);
			}

			if (add_peptide_set_to_psm(psm,set))
				total_num_sets++;

			num_sets++;
			if (has_correct)
				num_with_correct++;
			if (correct_first)
				num_first_correct++;
				
		}

		ifs.close();

		cout << "\tsets added:   " << num_sets << endl;
		cout << "\thad correct:  " << float(num_with_correct)/num_sets << endl;
		cout << "\tcorrect first:" << float(num_first_correct)/num_sets << endl;
		cout << endl;
	}

	cout << "Processed " << num_read << "/" << names.size() << " denovo files" << endl;
}




/*******************************************************************************
Reads the sets of predicted peptides:
list of db hits amd of complete de novo sequences for the spectrum.
********************************************************************************/
void create_complete_denovo_set_map(
						Config *config,
						const string& mgf_list,
						const string& db_dir,
						const string& correct_dir,
						const string& denovo_dir,
						int charge,
						int size_idx,
						PeptideSetMap& psm,
						vector<bool>& file_indicators)
{
	vector<string> base_names,list;
	readListOfPaths(mgf_list.c_str(), list);
	base_names.resize(list.size());

	int i;
	for (i=0; i<list.size(); i++)
		getFileNameWithoutExtension(list[i].c_str(),base_names[i]);

	file_indicators.clear();
	file_indicators.resize(list.size(),false);

	// create map
	psm.clear();

	char suf_buf[64];
	sprintf(suf_buf,"_dbh_%d_%d.txt",charge,size_idx);
	string db_suffix = string(suf_buf);
	add_db_hits_to_psm(config,base_names,db_dir,correct_dir,db_suffix,charge,psm);

	
	sprintf(suf_buf,"_full_dnv_%d_%d.txt",charge,size_idx);
	string dnv_full_suffix = string(suf_buf);
//	add_denovo_paths_to_psm(config,base_names,denovo_dir,db_suffix,charge,psm);
	
	add_denovo_dbh_to_psm(config,base_names,denovo_dir,correct_dir,db_suffix,charge,psm);
	

	
	//
	cout << "Read info for charge " << charge << " size " << size_idx << endl;
	vector<int> type_sum;
	type_sum.resize(10,0);

	int num_sets=0;
	int num_incorrect=0;

	for (PeptideSetMap::const_iterator it = psm.begin(); it != psm.end(); it++)
	{
		num_sets++;
		file_indicators[(*it).first.file_idx]=true;
		num_incorrect += (*it).second.incorrect_sols.size();
		int j;
		for (j=0; j<(*it).second.incorrect_sols.size(); j++)
			type_sum[(*it).second.incorrect_sols[j].type]++;
	}

	cout << "\t" << num_sets << "\tsets" << endl;
	cout << "\t" << num_incorrect <<"\tnegtaive samples" << endl;

	for (i=0; i<type_sum.size(); i++)
		if (type_sum[i]>0)
			cout << "\ttype " << i << " " << type_sum[i] << endl;

}


/*******************************************************************************
Reads the sets of predicted peptides:
list of db hits amd of complete de novo sequences for the spectrum.
********************************************************************************/
void create_complete_dbh_set_map(
						Config *config,
						const string& mgf_list,
						const string& db_dir,
						const string& correct_dir,
						int charge,
						int size_idx,
						PeptideSetMap& psm,
						vector<bool>& file_indicators)
{
	vector<string> base_names,list;
	readListOfPaths(mgf_list.c_str(), list);
	base_names.resize(list.size());

	int i;
	for (i=0; i<list.size(); i++)
		getFileNameWithoutExtension(list[i].c_str(),base_names[i]);

	file_indicators.clear();
	file_indicators.resize(list.size(),false);

	// create map
	psm.clear();

	char suf_buf[64];
	sprintf(suf_buf,"_dbh_%d_%d.txt",charge,size_idx);
	string db_suffix = string(suf_buf);
	add_db_hits_to_psm(config,base_names,db_dir,correct_dir,db_suffix,charge,psm);
	
	//
	cout << "Read info for charge " << charge << " size " << size_idx << endl;
	vector<int> type_sum;
	type_sum.resize(10,0);

	int num_sets=0;
	int num_incorrect=0;

	for (PeptideSetMap::const_iterator it = psm.begin(); it != psm.end(); it++)
	{
		num_sets++;
		file_indicators[(*it).first.file_idx]=true;
		num_incorrect += (*it).second.incorrect_sols.size();
		int j;
		for (j=0; j<(*it).second.incorrect_sols.size(); j++)
			type_sum[(*it).second.incorrect_sols[j].type]++;
	}

	cout << "\t" << num_sets << "\tsets" << endl;
	cout << "\t" << num_incorrect <<"\tnegtaive samples" << endl;

	for (i=0; i<type_sum.size(); i++)
		if (type_sum[i]>0)
			cout << "\ttype " << i << " " << type_sum[i] << endl;

}


// returns the total number of pairs in the psm
// that had weights assigned
// reassigns weights to sets according to weights of db / full de novo
int assign_denovo_weights_to_sets(PeptideSetMap& psm, 
								  float ratio_db, 
								  float ratio_full_denovo)
{
	const float sum = ratio_db + ratio_full_denovo;
	const float req_type1_weight = ratio_db/sum;
	const float req_type2_weight = 1 - req_type1_weight;
	int total_num_pairs =0;

	for (PeptideSetMap::iterator it = psm.begin(); it != psm.end(); it++)
	{
		PeptideSet& set = (*it).second;

		int num_type1=0;
		int num_type2=0;
		float total_type1_weight=0;
		float total_type2_weight=0;

		// set weight according to proportion of correct amino acids

		int i;
		for (i=0; i<set.incorrect_sols.size(); i++)
		{
			PeptideSolution& sol = set.incorrect_sols[i];
			const float correct_aa_ratio = ((float)sol.num_correct_aas/(float)sol.pep.get_num_aas());
			sol.weight = 1.0 - correct_aa_ratio;
			if (sol.type==1)
			{
				num_type1++;
				total_type1_weight+=sol.weight;
			}
			else if (sol.type==2)
			{
				num_type2++;
				total_type2_weight+=sol.weight;
			}
			else
			{
				cout << "Error: unknown de novo sample type: " << sol.type << endl;
				exit(1);
			}
			
		}
		total_num_pairs+=set.incorrect_sols.size();

		const float mul1= (num_type1>0 ? req_type1_weight/total_type1_weight : 0);
		const float mul2= (num_type2>0 ? req_type2_weight/total_type2_weight : 0);

		set.total_set_weight=0;
		for (i=0; i<set.incorrect_sols.size(); i++)
		{
			PeptideSolution& sol = set.incorrect_sols[i];

			sol.weight *= (sol.type == 1 ? mul1 : mul2);
			set.total_set_weight += sol.weight;
		}
	}
	return total_num_pairs;
}


// also removes pairs with weight 0
void re_weight_complete_denovo_phi_support(vector<SamplePairWeight>& phi_support,	
							float ratio_db, float ratio_denovo, float ratio_db_cross)
{
	const int num_types = 4;
	const double ratio_sum = (double)(ratio_db + ratio_denovo + ratio_db_cross);
	const double per_weights[num_types] = {0, (double)ratio_db / ratio_sum, (double)ratio_denovo/ratio_sum,
		(double)ratio_db_cross / ratio_sum};
	double total_weights[num_types]={0};
	int	   num_samples[num_types]={0};

	int i;
	for (i=0; i<phi_support.size(); i++)
	{
		const int type = phi_support[i].tag;
		total_weights[type]+=phi_support[i].weight;
		num_samples[type]++;
	}

	double total_weight=0;
	for (i=0; i<num_types; i++)
		total_weight+=total_weights[i];

	double mult_vals[num_types]={0};
	for (i=0; i<num_types; i++)
		if (per_weights[i]>0)
			mult_vals[i]= (per_weights[i]*total_weight) / total_weights[i];

	for (i=0; i<phi_support.size(); i++)
		phi_support[i].weight *= mult_vals[phi_support[i].tag];

	for (i=0; i<phi_support.size(); i++)
	{
		if (phi_support[i].weight<=0)
		{
			phi_support[i]=phi_support[phi_support.size()-1];
			phi_support.pop_back();
		}
	}


	for (i=0; i<num_types; i++)
	{
		total_weights[i]=0;
		num_samples[i]=0;
	}

	for (i=0; i<phi_support.size(); i++)
	{
		const int type = phi_support[i].tag;
		total_weights[type]+=phi_support[i].weight;
		num_samples[type]++;
	}

	total_weight=0;
	for (i=0; i<num_types; i++)
		total_weight+=total_weights[i];

	cout << endl << "Reweighted samples in set: " << endl;
	for (i=0; i<num_types; i++)
		if (num_samples[i]>0)
			cout << "Type " << i << " " << num_samples[i] << "\t" << per_weights[i] << "\t" << " : " <<
				total_weights[i] / total_weight << endl;
	cout << endl;
}


bool are_same_upto_NDILQK(const Peptide& pep1, const Peptide& pep2)
{
	if (pep1.get_num_aas() != pep2.get_num_aas())
		return false;

	const vector<int>& aas1 = pep1.get_amino_acids();
	const vector<int>& aas2 = pep2.get_amino_acids();
	int i;

	for (i=0; i<aas1.size()-1; i++)
	{
		if (aas1[i] == aas2[i] )
			continue;

		if ((aas1[i] == Ile && aas2[i] == Leu) ||
			(aas1[i] == Leu && aas2[i] == Ile) ||
			(aas1[i] == Asn && aas2[i] == Asp) ||
			(aas1[i] == Asp && aas2[i] == Asn) ||
			(aas1[i] == Gln && aas2[i] == Lys) ||
			(aas1[i] == Lys && aas2[i] == Gln) )
			continue;
		return false;
	}
	
	return (aas1[i]==aas2[i]);
}



void PeptideRankScorer::create_training_data_for_complete_denovo_ranking(	
				const string& db_dir,
				const string& correct_dir,
				const string& mgf_list,
				const double train_ratio,
				const int max_num_pairs,
				const int charge,
				const int size_idx,
				RankBoostDataset& train_ds, 
				RankBoostDataset& test_ds,
				vector<string>* peptide_strings,
				char *test_scan_file,
				float ratio_denovo,
				char *rerank_path,
				int   rerank_depth)
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(this->allScoreModelsPtr_);
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	Config *config = allScoreModels->get_config();
	const mass_t tolerance_diff = config->getTolerance()*0.5;
	const int max_num_ppp_frags = 4;

	bool print_sams=false;

	if (model_type != 2)
	{
		cout << "Error: this training function is only intended for full de novo sequence samples!" << endl;
		cout << "Need to set model type to 2, not " << model_type << endl;
		exit(1);
	}

	if (peak_model->get_feature_set_type() <= 2)
	{
		cout << "Error: training function intended for combined peak model!" << endl;
		exit(1);
	}

/*	vector<bool>   file_indicators;
	PeptideSetMap   psm;

	if (peptide_strings)
		peptide_strings->clear();

	double start_time = time(NULL);

	// read sample peptide sequences
	create_complete_dbh_set_map(config,mgf_list,db_dir,correct_dir,
		charge,size_idx,psm, file_indicators);

	// Read previous rerank path
	PeptideRankScorer previous_rerank_model;
	if (rerank_path)
	{

		cout << "Read previous rerank model, type " << previous_rerank_model.get_model_type() << endl;
		cout << "Option not supported anymore!" << endl;
		exit(1);
	}

	const float ratio_db = 1.0 - ratio_denovo;
	const int num_db_pairs     = (int)(max_num_pairs * ratio_db);
	const int num_denovo_pairs = (int)(max_num_pairs * ratio_denovo);

	// read spectra
	FileManager fm;
	FileSet	    fs;
	BasicSpecReader bsr;
	QCPeak		peaks[2000];

	fm.init_from_list_file(allScoreModels->get_config(),mgf_list.c_str(),file_indicators);
	fs.select_all_files(fm);

	const vector<SingleSpectrumFile *>& all_ssfs = fs.get_ssf_pointers();

	int i;
	int num_spectra_read=0;

	cout << endl << "Read " << all_ssfs.size() << " headers" << endl << endl;
	cout << "Creating RankBoostDatasets... proportion for training " << train_ratio << endl;
	cout << "Using following proportions for training data:" << endl;
	cout << ratio_denovo   << " pairs of correct, incorrect full dnv (same spectrum)" << endl;
	cout << ratio_db	   << " pairs of correct, incorrect db hits (same spectrum)" << endl;
	
	// first select most abundant fragments, use first 1000 ssfs for statistics
	vector<int> frag_counts;
	frag_counts.resize(config->get_all_fragments().size(),0);
	for (i=0; i<all_ssfs.size() && i<1000; i++)
	{
		MGF_single *ssf = (MGF_single *)all_ssfs[i];
		const int num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		BasicSpectrum     bs;
		bs.peaks = peaks;
		bs.num_peaks = num_peaks;
		bs.ssf = ssf;

		vector< vector<float> > ann_intens;
		vector< vector<mass_t> > ann_masses;
		AnnotatedSpectrum as;
		as.init_from_QCPeaks(config,peaks,num_peaks,ssf);
		as.set_peptide(ssf->peptide);
		as.annotate_spectrum(ssf->peptide.get_mass_with_19(),true);
		as.extract_annotated_intens_and_masses(ann_intens,ann_masses);

		int f;
		for (f=0; f<ann_intens.size(); f++)
		{
			int j;
			for (j=0; j<ann_intens[f].size(); j++)
				if (ann_intens[f][j]>0)
					frag_counts[f]++;
		}
	}


	vector<int> ppp_frag_type_idxs;
	ppp_frag_type_idxs.clear();

	cout << "Using a combined peak model" << endl;
	cout << "Selecting peaks of the paritally mobile model" << endl;

	int mobility;
	for (mobility=MOBILE+1; mobility<NONMOBILE; mobility++)
	{
		const vector<int>& frag_idxs = peak_model->get_model_ptr(charge,size_idx,mobility)->get_fragment_type_idxs();
		int f;
		for (f=0; f<frag_idxs.size(); f++)
		{
			if (ppp_frag_type_idxs.size() == max_num_ppp_frags)
				break;

			int j;
			for (j=0; j<ppp_frag_type_idxs.size(); j++)
				if (ppp_frag_type_idxs[j] == frag_idxs[f])
					break;
			if (j==ppp_frag_type_idxs.size())
				ppp_frag_type_idxs.push_back(frag_idxs[f]);
		}
	}
	int f;
	for (f=0; f<ppp_frag_type_idxs.size(); f++)
		cout << f+1 << "\t" << config->get_fragment(ppp_frag_type_idxs[f]).label << endl;
	cout << endl;
	

	sort(ppp_frag_type_idxs.begin(),ppp_frag_type_idxs.end());

	DeNovoPartitionModel *&part_model = dnv_part_models[charge][size_idx];
	part_model->init_features(model_type,charge,size_idx,ppp_frag_type_idxs,config);
		
	int num_groups_in_train=0;
	int num_groups_in_test=0;
	int num_train_pairs=0;

	ofstream test_scan_stream;
	if (test_scan_file)
	{
		test_scan_stream.open(test_scan_file);
		if (! test_scan_stream.is_open())
		{
			cout << "Error: couldn't open test stream!" << endl;
			exit(1);
		}
	}
	
	int set_counter=0;
	for (i=0; i<all_ssfs.size(); i++)
	{
		PeptideSetMap::const_iterator it;
		MGF_single *ssf = (MGF_single *)all_ssfs[i];
		scan_pair key(ssf->file_idx,ssf->idx_in_file);

		it = psm.find(key);
		if (it == psm.end())
			continue;

		set_counter++;
	}

	cout << set_counter << " sets available." << endl;
	int num_db_peptides_per_set = (int)(0.5 + (float)num_db_pairs/set_counter);
	if (num_db_peptides_per_set>4)
		num_db_peptides_per_set=4;
	if (num_db_peptides_per_set<1)
		num_db_peptides_per_set=1;

	int db_pairs = num_db_peptides_per_set * set_counter;

	int num_denovo_peptides_per_set = (int)(0.5 + (float)(max_num_pairs-db_pairs)/set_counter);
	if (num_denovo_peptides_per_set<1)
		num_denovo_peptides_per_set=1;
	if (num_denovo_peptides_per_set>25)
		num_denovo_peptides_per_set=25;

	cout << "NUM DB PER SET     : " << num_db_peptides_per_set << endl;
	cout << "NUM DE NOVO PER SET: " << num_denovo_peptides_per_set << endl;


	static PrmGraph *prm_ptr=NULL;
	static vector<PrmGraph *> prm_ptrs;
	static vector<SeqPath>    seqpath_solutions;

	int num_rand=0;
	int num_first=0;

	// Generate various types of samples from spectra
	int ssf_idx;
	for (ssf_idx=0; ssf_idx<all_ssfs.size(); ssf_idx++)
	{
		MGF_single *ssf = (MGF_single *)all_ssfs[ssf_idx];
		const mass_t true_mass_with_19 = ssf->peptide.get_mass_with_19();
		const Peptide& correct_peptide = ssf->peptide;
		PeptideSetMap::const_iterator it;
		
		scan_pair key(ssf->file_idx,ssf->idx_in_file);

	//	if (num_groups_in_train>800)
	//		break;

		it = psm.find(key);
		if (it == psm.end())
			continue;

		const vector<int>& aas = ssf->peptide.get_amino_acids();
		int a;
		for (a=0; a<aas.size(); a++)
			if (aas[a]>Val)
				break;
		if (a<aas.size())
			continue;

		const PeptideSet& set = (*it).second;
		BasicSpectrum     bs;
		AnnotatedSpectrum as;
		vector<PmcSqsChargeRes> pmc_sqs_res;

		if (fabs(set.correct_sol.pep.get_mass_with_19() - ssf->peptide.get_mass_with_19())>0.001)
		{
			cout << "Error: mismatch between set peptide and peptide in file:" << endl;
			cout << "Set : " << set.correct_sol.pep.as_string(config) << "\t" << set.correct_sol.pep.get_mass_with_19() <<endl;
			cout << "Spec: " << ssf->peptide.as_string(config) << "\t" << ssf->peptide.get_mass_with_19() << endl;
		}
		
		int charge1=0,charge2=0;
		mass_t mz1=0,mz2=0;
		float prob1=0,prob2=0;

		const int num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		bs.peaks = peaks;
		bs.num_peaks = num_peaks;
		bs.ssf = ssf;

		// calc corrected pm_with_19, if it is good use it, otherwise, use a value with +- U[0,toleance/2]
	//	allScoreModels->get_best_mz_charge(config,bs,&mz1,&charge1,&prob1,&mz2,&charge2,&prob2,&pmc_sqs_res);
		mass_t pm_with_19=NEG_INF;
		bool good_first = true;
		const mass_t corr1_pm_with_19 = mz1*charge1 - MASS_PROTON*(charge1-1);
		const mass_t corr2_pm_with_19 = mz2*charge2 - MASS_PROTON*(charge2-1);
		if (fabs(corr2_pm_with_19-true_mass_with_19)<tolerance_diff)
		{
			pm_with_19 = corr2_pm_with_19;
			good_first=false;
		}

		if (fabs(corr1_pm_with_19-true_mass_with_19)<tolerance_diff)
		{
			pm_with_19 = corr1_pm_with_19;
		}
		else
			good_first=false;
		
		bool bad_pm=false;
	
		if (pm_with_19<0) // use a random value
		{
			double r=myRandom();
			mass_t offset = r*r*tolerance_diff;
			if (myRandom()<0.5)
				offset *= -1;
			pm_with_19 = true_mass_with_19 + offset;
			bad_pm=true;
		}


		as.init_from_QCPeaks(config,peaks,num_peaks,ssf);
		as.set_corrected_pm_with_19(pm_with_19);
	
		const bool add_to_train = (myRandom() <= train_ratio);
		RankBoostDataset& rds = (add_to_train ? train_ds : test_ds);

		if (test_scan_file && ! add_to_train)
			test_scan_stream << ssf->file_idx << " " << ssf->idx_in_file << " " <<
			ssf->peptide.get_num_aas() << endl;
	
		bool first=true;
		int corr_sam_idx=-1;
		int groupIndex=-1;

		if (add_to_train)
		{
			groupIndex = num_groups_in_train++;
		}
		else
		{
			groupIndex = num_groups_in_test++;
		}

		RankBoostSample corr_rbs;
		//fill_complete_peptide_rbs(set.correct_sol, peaks, num_peaks, as, pmc_sqs_res, corr_rbs, size_idx);
		corr_sam_idx = rds.get_num_samples();
		corr_rbs.groupIndex=groupIndex;
		corr_rbs.tag1=NEG_INF;
		corr_rbs.float_tag1 = set.correct_sol.pm_with_19;

		if (peptide_strings && add_to_train)
		{
			corr_rbs.tag1 = peptide_strings->size();
			peptide_strings->push_back(set.correct_sol.pep.as_string(config));
		}

		corr_rbs.rank_in_group=0;
		rds.add_sample(corr_rbs);

		if (print_sams)
		{
			float score;
			corr_rbs.get_feature_val(117,&score);
			cout << endl << ssf_idx << "\t"<< set.correct_sol.pep.as_string(config) << "\t" << score << endl;
		}

		// add db samples
		int j;
		for (j=0; j<set.incorrect_sols.size() && j<num_db_peptides_per_set; j++)
		{
			const int type = set.incorrect_sols[j].type;
			if (type != 1)
				break;
			
			if (are_same_upto_NDILQK(set.correct_sol.pep,set.incorrect_sols[j].pep))
				continue;

			RankBoostSample bad_rbs;
			//fill_complete_peptide_rbs(set.incorrect_sols[j], peaks, num_peaks, as, pmc_sqs_res, 
			//						bad_rbs, size_idx);

			int bad_sam_idx=rds.get_num_samples();

			bad_rbs.groupIndex = groupIndex;
			bad_rbs.rank_in_group=(j == 0 ? 1 : 2);

			if (peptide_strings && add_to_train)
			{		
				bad_rbs.tag1=peptide_strings->size();	
				peptide_strings->push_back(set.incorrect_sols[j].pep.as_string(config));
			}

			bad_rbs.tag2=set.incorrect_sols[j].type;
			bad_rbs.tag3=set.incorrect_sols[j].num_correct_aas;
			bad_rbs.float_tag1 = set.correct_sol.pm_with_19;

			rds.add_sample(bad_rbs);
			rds.add_to_phi_vector(bad_sam_idx, corr_sam_idx, set.incorrect_sols[j].type); 
			if (add_to_train)
				num_train_pairs++;

			if (print_sams)
			{
				float score;
				bad_rbs.get_feature_val(117,&score);
				cout << j << "\t" <<  
					set.incorrect_sols[j].pep.as_string(config) << "\t" << score << endl;
			}
		}
	
		// add denovo samples (if the rerank model exists, rerank them)
		if (num_denovo_peptides_per_set>0)
		{
			vector<mass_t> pms_with_19;
			vector<int> charges;
			pms_with_19.push_back(pm_with_19);
			charges.push_back(charge);

	

			if (prm_ptrs.size()<pms_with_19.size())
				prm_ptrs.resize(pms_with_19.size(),NULL);
			
			generate_denovo_solutions_from_several_pms_with_good_start_end_idxs(prm_ptrs,
				allScoreModels,&as,true,300,6,14,pms_with_19,charges,seqpath_solutions);

		//	generate_denovo_solutions_with_good_start_end_idxs(prm_ptr,model,&as,true,
		//								pm_with_19,charge,300,6,14,seqpath_solutions);

	
			int j;
			for (j=0; j<seqpath_solutions.size(); j++)
			{
				const int num_correct_aas = seqpath_solutions[j].get_num_correct_aas(correct_peptide,config);
				if (num_correct_aas==seqpath_solutions[j].get_num_aa())
				{
					seqpath_solutions[j]=seqpath_solutions[seqpath_solutions.size()-1];
					seqpath_solutions.pop_back();
				}
			}
		
			vector<PeptideSolution> pep_solutions;
			if (seqpath_solutions.size()>0)
			{
				// create peptide solutions
				pep_solutions.resize(seqpath_solutions.size());
				int s_idx;
				for (s_idx=0; s_idx<seqpath_solutions.size(); s_idx++)
				{
					convert_seq_path_to_peptide_soluition_and_fill_in_aas(config,
						correct_peptide, seqpath_solutions[s_idx],pep_solutions[s_idx]);
				}
			}
			else
				continue;


			vector<int> pep_sol_idxs;
			pep_sol_idxs.clear();
		
			vector<score_pair> rerank_scores;
			// re-order the de novo solutions according to previous model
			if (0 || rerank_path)
			{
				
			//	previous_rerank_model.score_complete_sequences(pep_solutions,ssf,peaks,
			//						  bs.num_peaks, rerank_scores, size_idx);
				
				sort(rerank_scores.begin(),rerank_scores.end());
				pep_sol_idxs.push_back(0); // include top-scoring de novo in any case
				pep_sol_idxs.push_back(1);

				const int num_top = num_denovo_peptides_per_set/2+1;
				int j;
				for (j=0; j<num_top; j++)
					pep_sol_idxs.push_back(rerank_scores[j].idx);
				
				const int num_left = num_denovo_peptides_per_set - pep_sol_idxs.size();

				for (j=0; j<num_left; j++)
				{
					int p = int(myRandom()*(rerank_scores.size()-pep_sol_idxs.size()))+pep_sol_idxs.size();
					pep_sol_idxs.push_back(rerank_scores[p].idx);
				}
				sort(pep_sol_idxs.begin(),pep_sol_idxs.end());
			}
			else // use original order
			{
				const int num_top = int(num_denovo_peptides_per_set*0.6)+1;
				int j;
				for (j=0; j<num_top; j++)
					pep_sol_idxs.push_back(j);
				
				const int num_left = num_denovo_peptides_per_set - pep_sol_idxs.size();
				for (j=0; j<num_left; j++)
				{
					int p = int(myRandom()*(pep_solutions.size()-pep_sol_idxs.size()))+pep_sol_idxs.size();
					pep_sol_idxs.push_back(p);
				}
			}


			sort(pep_sol_idxs.begin(),pep_sol_idxs.end());

			// add samples
			for (j=0; j<pep_sol_idxs.size(); j++)
			{
				const int pep_sol_idx = pep_sol_idxs[j];

				if (are_same_upto_NDILQK(set.correct_sol.pep,pep_solutions[pep_sol_idx].pep))
					continue;
				
				RankBoostSample bad_rbs;
				//fill_complete_peptide_rbs(pep_solutions[pep_sol_idx], peaks, num_peaks, as, 
				//		pmc_sqs_res, bad_rbs, size_idx);

				const int bad_sam_idx=rds.get_num_samples();

				bad_rbs.groupIndex = groupIndex;
				bad_rbs.rank_in_group=2;

				if (peptide_strings && add_to_train)
				{		
					bad_rbs.tag1=peptide_strings->size();	
					peptide_strings->push_back(pep_solutions[pep_sol_idx].pep.as_string(config));
				}

				// give higher ranked db mismatches a stronger weight;
				bad_rbs.tag2= 2; // full de novo
				bad_rbs.float_tag1 = set.correct_sol.pm_with_19;
				
				rds.add_sample(bad_rbs);
				rds.add_to_phi_vector(bad_sam_idx, corr_sam_idx, 2);
	
				if (add_to_train)
					num_train_pairs++;

				if (print_sams)
				{
					float score;
					bad_rbs.get_feature_val(117,&score);
					cout << j << "\t" << pep_sol_idx << "\t" << 
						pep_solutions[pep_sol_idx].pep.as_string(config) << "\t" << score << endl;
				}
			}
		}

		if (bad_pm)
			num_rand++;

		if (good_first)
			num_first++;
		
		int num_denovo_to_add=0;
	
		num_spectra_read++;

		if ((model_type == 2 && num_spectra_read % 200 == 0) ||
		    (model_type == 0 && num_spectra_read % 1000 == 0) )
		{
			cout << "processed " << num_spectra_read << "/" << set_counter << " spectra (";
			double curr_time = time(NULL);
			cout << curr_time - start_time << " secs.)\t" << setprecision(3) << fixed << " pm acc: " <<
				1.0 - num_rand/(float)num_spectra_read <<"\tfirst: " << num_first/(float)num_spectra_read<<  endl;
		}
	}

	if (test_scan_file)
		test_scan_stream.close();

	train_ds.set_num_groups(num_groups_in_train);

	test_ds.set_num_groups(num_groups_in_test);
			
	cout << "computing phi weights..." << endl;
	train_ds.compute_total_phi_weight();
	
	cout << "initializing potential lists..." << endl;
	train_ds.initialize_potenital_lists();
	
	cout << "initializing real feature table..." << endl;
	train_ds.initialzie_real_feature_table(part_model->get_feature_names().size());

	cout << "Processed " << num_spectra_read << " spectra" << endl;
	cout << num_groups_in_train << " sets in training set (" << train_ds.get_phi_support().size() << " pairs)" << endl;
	cout << num_groups_in_test  << " sets in test set     (" << test_ds.get_phi_support().size() << " pairs)" << endl << endl;
	*/
}


/*************************************************************************************
Creates the training data for models of complete predictions, i.e., it is assumed that
all peptides being scores start at mass 0 and end at the c-terminal. Therefore, we can
use the peptides mass as the spectrum's pm_with_19.
**************************************************************************************/
void PeptideRankScorer::create_training_data_for_complete_sequence_ranking(
				const string& db_dir,
				const string& correct_dir,
				const string& denovo_dir,
				const string& mgf_list,
				const double train_ratio,
				const int max_num_pairs,
				const int charge,
				const int size_idx,
				RankBoostDataset& train_ds, 
				RankBoostDataset& test_ds,
				vector<string>* peptide_strings,
				char *test_scan_file,
				float ratio_db,
				float ratio_denovo,
				float ratio_db_cross)
{	
/*	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(this->allScoreModelsPtr_);
	const int max_num_ppp_frags = 4;
	if (model_type !=0 && model_type != 2)
	{
		cout << "Error: this training function is only intended for full sequence samples!" << endl;
		cout << "Need to set model type to 0  or 2, not " << model_type << endl;
		exit(1);
	}

	if (ratio_db<=0)
		ratio_db=0;

	if (ratio_denovo<0)
		ratio_denovo=0;

	if (ratio_db_cross<0)
		ratio_db_cross=0;

	vector<bool>   file_indicators;
	PeptideSetMap	 psm;

	Config *config = allScoreModels->get_config();

	if (peptide_strings)
		peptide_strings->clear();

	double start_time = time(NULL);

	// read sample peptide sequences
	create_complete_denovo_set_map(config,mgf_list,db_dir,correct_dir,
		denovo_dir,charge,size_idx,psm, file_indicators);

	// Read previous rerank path
	if (model_type == 2)
		ratio_db_cross = 0;


	const int total_num_pairs=assign_denovo_weights_to_sets(psm, ratio_db, ratio_denovo);
	const float total_ratio = ratio_db + ratio_denovo+ratio_db_cross;
	const int num_db_pairs     = (int)(max_num_pairs * (ratio_db/total_ratio));
	const int num_denovo_pairs = (int)(max_num_pairs * (ratio_denovo/total_ratio));

	// read spectra
	FileManager fm;
	FileSet	    fs;
	BasicSpecReader bsr;
	QCPeak		peaks[2000];

	fm.init_from_list_file(allScoreModels->get_config(),mgf_list.c_str(),file_indicators);
	fs.select_all_files(fm);

	const vector<SingleSpectrumFile *>& all_ssfs = fs.get_ssf_pointers();

	int i;
	int num_spectra_read=0;

	cout << endl << "Read " << all_ssfs.size() << " headers" << endl << endl;
	cout << "Creating RankBoostDatasets... proportion for training " << train_ratio << endl;
	cout << "Using following proportions for training data:" << endl;
	cout << ratio_db << " pairs of correct, incorrect db hits (same spectrum)" << endl;
	cout << ratio_denovo   << " pairs of correct, incorrect full dnv (same spectrum)" << endl;
	cout << ratio_db_cross << " pairs of correct, incorrect db hits (different spectra)" << endl;
	
	// first select most abundant fragments, use first 1000 ssfs for statistics
	vector<int> frag_counts;
	frag_counts.resize(config->get_all_fragments().size(),0);
	for (i=0; i<all_ssfs.size() && i<1000; i++)
	{
		MGF_single *ssf = (MGF_single *)all_ssfs[i];
		const int num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		BasicSpectrum     bs;
		bs.peaks = peaks;
		bs.num_peaks = num_peaks;
		bs.ssf = ssf;

		vector< vector<float> > ann_intens;
		vector< vector<mass_t> > ann_masses;
		AnnotatedSpectrum as;
		as.init_from_QCPeaks(config,peaks,num_peaks,ssf);
		as.set_peptide(ssf->peptide);
		as.annotate_spectrum(ssf->peptide.get_mass_with_19(),true);
		as.extract_annotated_intens_and_masses(ann_intens,ann_masses);

		int f;
		for (f=0; f<ann_intens.size(); f++)
		{
			int j;
			for (j=0; j<ann_intens[f].size(); j++)
				if (ann_intens[f][j]>0)
					frag_counts[f]++;
		}
	}


	vector<int> ppp_frag_type_idxs;
	ppp_frag_type_idxs.clear();

	// The peak prediction fragments are used elsewhere too
	// select them here unless they are provided by the PeakRankModel
	// (if it is a combined model)
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	if (peak_model->get_feature_set_type() <= 2)
	{
		cout << "Selecting frag models:" << endl;
		
		for (i=0; i<max_num_ppp_frags; i++)
		{
			int max_frag=-1;
			int max_count=-1;
			int f;
			for (f=0; f<frag_counts.size(); f++)
				if (frag_counts[f]>max_count && 
					peak_model->has_intialized_models(charge,size_idx,f))
				{
					max_count = frag_counts[f];
					max_frag=f;
				}

			if (max_count<0)
			{
				cout << "Warning: not enough frag models initialized for charge " << charge <<
					"  size " << size_idx << endl;
				break;
			}
			ppp_frag_type_idxs.push_back(max_frag);
			cout << max_frag << "\t" << config->get_fragment(max_frag).label << "\t" << frag_counts[max_frag] << endl;
			frag_counts[max_frag]=-1;
		}
	}
	else
	{
		cout << "Using a combined peak model" << endl;
		cout << "Selecting peaks of the paritally mobile model" << endl;

		int mobility;
		for (mobility=MOBILE+1; mobility<NONMOBILE; mobility++)
		{
			const vector<int>& frag_idxs = peak_model->get_model_ptr(charge,size_idx,mobility)->get_fragment_type_idxs();
			int f;
			for (f=0; f<frag_idxs.size(); f++)
			{
				if (ppp_frag_type_idxs.size() == max_num_ppp_frags)
					break;

				int j;
				for (j=0; j<ppp_frag_type_idxs.size(); j++)
					if (ppp_frag_type_idxs[j] == frag_idxs[f])
						break;
				if (j==ppp_frag_type_idxs.size())
					ppp_frag_type_idxs.push_back(frag_idxs[f]);
			}
		}
		int f;
		for (f=0; f<ppp_frag_type_idxs.size(); f++)
			cout << f+1 << "\t" << config->get_fragment(ppp_frag_type_idxs[f]).label << endl;
		cout << endl;
	}

	sort(ppp_frag_type_idxs.begin(),ppp_frag_type_idxs.end());

	DeNovoPartitionModel *&part_model = dnv_part_models[charge][size_idx];
	part_model->init_features(model_type,charge,size_idx,ppp_frag_type_idxs,config);
		
	int num_groups_in_train=0;
	int num_groups_in_test=0;
	int num_train_pairs=0;

	ofstream test_scan_stream;
	if (test_scan_file)
	{
		test_scan_stream.open(test_scan_file);
		if (! test_scan_stream.is_open())
		{
			cout << "Error: couldn't open test stream!" << endl;
			exit(1);
		}
	}
	
	int set_counter=0;
	for (i=0; i<all_ssfs.size(); i++)
	{
		PeptideSetMap::const_iterator it;
		MGF_single *ssf = (MGF_single *)all_ssfs[i];
		scan_pair key(ssf->file_idx,ssf->idx_in_file);

		it = psm.find(key);
		if (it == psm.end())
			continue;

		set_counter++;
	}

	cout << set_counter << " sets available." << endl;
	int num_db_peptides_per_set = (int)(0.5 + (float)num_db_pairs/set_counter);
	if (num_db_peptides_per_set>7)
		num_db_peptides_per_set=7;
	if (num_db_peptides_per_set<1)
		num_db_peptides_per_set=1;
	if (ratio_db_cross>0 && num_db_peptides_per_set<4)
		num_db_peptides_per_set=4;

	int num_denovo_peptides_per_set = (int)(0.75 + (float)num_denovo_pairs/set_counter);
	if (num_denovo_peptides_per_set<1)
		num_denovo_peptides_per_set=1;
	if (num_denovo_peptides_per_set>10)
		num_denovo_peptides_per_set=10;

	const int org_num_denovo_peptides_per_set = num_denovo_peptides_per_set;

	cout << "NUM DB PER SET     : " << num_db_peptides_per_set << endl;
	cout << "NUM DE NOVO PER SET: " << num_denovo_peptides_per_set << endl;

	static vector<PrmGraph *> prm_ptrs;
	static vector<SeqPath>    seqpath_solutions;

	// Generate various types of samples from spectra
	for (i=0; i<all_ssfs.size(); i++)
	{
		PeptideSetMap::const_iterator it;
		MGF_single *ssf = (MGF_single *)all_ssfs[i];

		// remove samples with PTMs
		const vector<int>& aas = ssf->peptide.get_amino_acids();
		int j;
		for (j=0; j<aas.size(); j++)
			if (aas[j]>Val)
				break;
		if (j<aas.size())
			continue;

		scan_pair key(ssf->file_idx,ssf->idx_in_file);

		it = psm.find(key);
		if (it == psm.end())
			continue;

		const PeptideSet& set = (*it).second;
		BasicSpectrum     bs;
		AnnotatedSpectrum as;
		vector<PmcSqsChargeRes> pmc_sqs_res;

		int charge1=0,charge2=0;
		mass_t mz1=0,mz2=0;
		float prob1=0,prob2=0;

		const int num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		bs.peaks = peaks;
		bs.num_peaks = num_peaks;
		bs.ssf = ssf;

	//	allScoreModels->get_best_mz_charge(config,bs, &mz1, &charge1, &prob1,
	//							  &mz2, &charge2, &prob2, &pmc_sqs_res);

		const mass_t mass_with_19 = mz1*charge1 - (charge1-1.0);
		as.init_from_QCPeaks(config,peaks,num_peaks,ssf);
		as.set_corrected_pm_with_19(mass_with_19);
	
		const bool add_to_train = (myRandom() <= train_ratio);
		RankBoostDataset& rds = (add_to_train ? train_ds : test_ds);

		if (test_scan_file && ! add_to_train)
			test_scan_stream << ssf->file_idx << " " << ssf->idx_in_file << " " <<
			ssf->peptide.get_num_aas() << endl;
	
		bool first=true;
		int corr_sam_idx=-1;
		int groupIndex=-1;

		if (add_to_train)
		{
			groupIndex = num_groups_in_train++;
		}
		else
		{
			groupIndex = num_groups_in_test++;
		}

		RankBoostSample corr_rbs;
		//fill_complete_peptide_rbs(set.correct_sol, peaks, num_peaks, as, pmc_sqs_res, corr_rbs, size_idx);
		corr_sam_idx = rds.get_num_samples();
		corr_rbs.groupIndex=groupIndex;
		corr_rbs.tag1=NEG_INF;
		corr_rbs.float_tag1 = set.correct_sol.pm_with_19;

		if (peptide_strings && add_to_train)
		{
			corr_rbs.tag1 = peptide_strings->size();
			peptide_strings->push_back(set.correct_sol.pep.as_string(config));
		}

		corr_rbs.rank_in_group=0;
		rds.add_sample(corr_rbs);

		// add db samples
		for (j=0; j<set.incorrect_sols.size() && j<num_db_peptides_per_set; j++)
		{
			const int type = set.incorrect_sols[j].type;
			if (type != 1)
				break;
			
			RankBoostSample bad_rbs;
			//fill_complete_peptide_rbs(set.incorrect_sols[j], peaks, num_peaks, as, pmc_sqs_res, 
			//						bad_rbs, size_idx);

			int bad_sam_idx=rds.get_num_samples();

			bad_rbs.groupIndex = groupIndex;
			bad_rbs.rank_in_group=(j == 0 ? 1 : 2);

			if (peptide_strings && add_to_train)
			{		
				bad_rbs.tag1=peptide_strings->size();	
				peptide_strings->push_back(set.incorrect_sols[j].pep.as_string(config));
			}

			// give higher ranked db mismatches a stronger weight;

			bad_rbs.tag2=set.incorrect_sols[j].type;
			bad_rbs.tag3=set.incorrect_sols[j].num_correct_aas;
			bad_rbs.float_tag1 = set.correct_sol.pm_with_19;

			rds.add_sample(bad_rbs);
			rds.add_to_phi_vector(bad_sam_idx, corr_sam_idx, set.incorrect_sols[j].type); 
			
			if (add_to_train)
				num_train_pairs++;
		}
		int num_db_added=j;

		while (j<set.incorrect_sols.size() && set.incorrect_sols[j].type==1)
			j++;

		const int den_start_idx = j;
		
		// add the rerank samples
		int num_denovo_to_add = org_num_denovo_peptides_per_set;

		// add precomputed denovo samples
		int num_den_available = set.incorrect_sols.size() - den_start_idx;
		vector<size_t> den_idxs;
		den_idxs.clear();
		if (num_den_available>num_denovo_to_add)
		{
			chooseKFromN(num_denovo_to_add,num_den_available,den_idxs);
			int k;
			for (k=0; k<num_denovo_to_add; k++)
				den_idxs[k]+=den_start_idx;
		}
		else
		{
			int k;
			for (k=0; k<num_den_available; k++)
				den_idxs.push_back(den_start_idx+k);
		}

		for (j=0; j<den_idxs.size(); j++)
		{
			const int sam_idx = den_idxs[j];
			const int type = set.incorrect_sols[sam_idx].type;
			if (type != 2)
			{
				cout << "Error: adding non denovo type: " << set.incorrect_sols[sam_idx].type << endl;
				cout << "sam idx: " << sam_idx << "  out of " << set.incorrect_sols.size() << endl;
				int k;
				for (k=0; k<set.incorrect_sols.size(); k++)
				{
					cout << k << "\t" << set.incorrect_sols[k].type << endl;
				}
				exit(1);
			}

			RankBoostSample bad_rbs;
			//fill_complete_peptide_rbs(set.incorrect_sols[sam_idx], peaks, num_peaks, as, 
			//	pmc_sqs_res, bad_rbs, size_idx);

			const int bad_sam_idx=rds.get_num_samples();

			bad_rbs.groupIndex = groupIndex;
			bad_rbs.rank_in_group=2;

			if (peptide_strings && add_to_train)
			{		
				bad_rbs.tag1=peptide_strings->size();	
				peptide_strings->push_back(set.incorrect_sols[sam_idx].pep.as_string(config));

			}

			// give higher ranked db mismatches a stronger weight;
			bad_rbs.tag2=set.incorrect_sols[sam_idx].type;
			bad_rbs.tag3=set.incorrect_sols[sam_idx].num_correct_aas;

			rds.add_sample(bad_rbs);
			rds.add_to_phi_vector(bad_sam_idx, corr_sam_idx, set.incorrect_sols[sam_idx].type);

			if (add_to_train)
				num_train_pairs++;
		}
		num_spectra_read++;

		if ((model_type == 2 && num_spectra_read % 200 == 0) ||
		    (model_type == 0 && num_spectra_read % 1000 == 0) )
		{
			cout << "processed " << num_spectra_read << "/" << set_counter << " spectra (";
			double curr_time = time(NULL);
			cout << curr_time - start_time << " secs.)" << endl;
		}
	}

	if (test_scan_file)
		test_scan_stream.close();

	
	const int total_num_cross = max_num_pairs - train_ds.get_phi_support().size() - test_ds.get_phi_support().size();
	const int num_train_db_cross = (int)(total_num_cross * train_ratio);
	const int num_test_db_cross  = total_num_cross - num_train_db_cross;

	if (ratio_db_cross>0)
	{
		cout  << "Adding " << num_train_db_cross << " db cross samples to train" <<endl;
		// add cross db samples to train
		if (num_train_db_cross>0) 
		{
			const vector<RankBoostSample>& samples = train_ds.get_samples();
			vector<SamplePairWeight>& phi_support  = train_ds.get_non_const_phi_support();
			
			int i;
			for (i=0; i<num_train_db_cross; i++)
			{
				const int corr_idx = (int)(num_train_pairs * myRandom());
				int badIndex=-1;
				
				while (badIndex<0 || phi_support[badIndex].tag != SOL_INCORRECT_DB) 
				//	samples[phi_support[badIndex].idx1].rank_in_group != 1) // cross db with strongest bad db hit
					badIndex = (int)(num_train_pairs * myRandom());
			
				SamplePairWeight new_pair = phi_support[badIndex];
				new_pair.idx2 = phi_support[corr_idx].idx2;
				new_pair.tag = SOL_INCORRECT_DB_CROSS;
				new_pair.weight=1.0;
				phi_support.push_back(new_pair);
			}
			re_weight_complete_denovo_phi_support(phi_support, ratio_db, ratio_denovo, ratio_db_cross);
		}

		cout  << "Adding " << num_test_db_cross << " db cross samples to test" <<endl;
		// add cross db samples to test
		if (num_test_db_cross>0)
		{
			const vector<RankBoostSample>& samples = test_ds.get_samples();
			vector<SamplePairWeight>& phi_support  = test_ds.get_non_const_phi_support();
			const int org_support_size = phi_support.size();
			for (i=0; i<num_test_db_cross; i++)
			{
				const int corr_idx = (int)(org_support_size * myRandom());
				int badIndex=-1;
				
				while (badIndex<0 || phi_support[badIndex].tag != SOL_INCORRECT_DB )
				//	samples[phi_support[badIndex].idx1].rank_in_group != 1) // cross db with strongest bad db hit
					badIndex = (int)(org_support_size * myRandom());

				SamplePairWeight new_pair = phi_support[badIndex];
				new_pair.idx2 = phi_support[corr_idx].idx2;
				new_pair.tag = SOL_INCORRECT_DB_CROSS;
				new_pair.weight=1.0;
				phi_support.push_back(new_pair);
			}
			re_weight_complete_denovo_phi_support(phi_support,ratio_db,ratio_denovo,ratio_db_cross);
		}
	}


	train_ds.set_num_groups(num_groups_in_train);

	test_ds.set_num_groups(num_groups_in_test);
			
	cout << "computing phi weights..." << endl;
	train_ds.compute_total_phi_weight();
	
	cout << "initializing potential lists..." << endl;
	train_ds.initialize_potenital_lists();
	
	cout << "initializing real feature table..." << endl;
	train_ds.initialzie_real_feature_table(part_model->get_feature_names().size());

	cout << "Processed " << num_spectra_read << " spectra" << endl;
	cout << num_groups_in_train << " sets in training set (" << train_ds.get_phi_support().size() << " pairs)" << endl;
	cout << num_groups_in_test  << " sets in test set     (" << test_ds.get_phi_support().size() << " pairs)" << endl << endl;
	*/
}



bool cmp_rerank_score(const SeqPath& a, const SeqPath& b)
{
	return (a.rerank_score>b.rerank_score);
}


void PeptideRankScorer::select_sample_pairs(const vector<SeqPath>& solutions, 
										   const vector<int>& corr_idxs,
										   const vector<int>& bad_idxs, 
										   vector<weight_pair>& sample_pairs, 
										   int num_pairs) const
{
	const int num_iters=3*num_pairs;
	const double log_bad_size = log((double)bad_idxs.size());
	sample_pairs.clear();
	float total_w = 0;

	int i;
	for (i=0; i<num_iters && sample_pairs.size()<num_pairs; i++)
	{
		int corr_pos = 0;
		int bad_pos = i;

		if (i>2 && corr_idxs.size()>0)
		{
			double r = myRandom();
			if ((r<0.4 && solutions[0].org_rank>0) || r<0.15)
			{
				r=0;
			}
			else
			{
				r = myRandom();
			}
			corr_pos = int(r*(double)corr_idxs.size()); 
		}

		const int corr_idx = corr_idxs[corr_pos];
		const int rank_corr = solutions[corr_idx].org_rank;
		
		if (! (i<3 && rank_corr < 25))
		{
			double log_corr_rank = log(1.0 + rank_corr);
			double min_log = 0;
			if (rank_corr>150)
				min_log = log_corr_rank - 3.5;
			
			double max_u = exp(min_log+(log_bad_size+1-min_log) * myRandom());
			if (max_u>(double)bad_idxs.size())
				max_u=(double)bad_idxs.size();

			double min_u = exp(min_log)-1;
			
			bad_pos = (int)(min_u + myRandom()*(max_u-min_u));
		}

	/*	double min_log,max_log;
		if (rank_corr<=5)
		{
			min_log=0;
			max_log=2.75;
		} 
		else if (rank_corr<=25)
		{
			min_log=0;
			max_log=3.25;
		}
		else if (rank_corr<500)
		{
			min_log=log_corr_rank-3.0;
			max_log=log_corr_rank;
		}
		else
		{
			min_log=log_corr_rank-3.0;
			max_log=log_corr_rank-1.0;
		}

		
		if (max_log>log_bad_size)
			max_log = log_bad_size;

		if (min_log<0)
			min_log=0;

		int bad_pos=i;
		if (rank_corr<20 && i<4)
		{
			bad_pos = i;
			if (i==3)
				bad_pos+=(int)(myRandom()*3);
		}
		else
		{
			if (myRandom()<0.8)
			{
				bad_pos = (int)(exp(min_log + (myRandom()*(max_log-min_log))))-1;
			}
			else
			{
				double m = (log_corr_rank < log_bad_size-1.5 ? log_corr_rank : log_bad_size -1.5);
				bad_pos = (int)(exp(m+ (myRandom()*(log_bad_size-m))))-1;
			}
		} */
		
		if (bad_pos<0)
			bad_pos=0;
		if (bad_pos>=bad_idxs.size())
			continue;

		const int badIndex  = bad_idxs[bad_pos];
		int j;
		for (j=0; j<sample_pairs.size(); j++)
			if (sample_pairs[j].idx_bad == badIndex && sample_pairs[j].idx_corr == corr_idx)
				break;

		if (j<sample_pairs.size())
			continue;

		
		int rank_bad  = solutions[badIndex].org_rank;
	//	float w = (rank_bad<rank_corr ? 1.0 : 1.0 + 0.333*(log((double)(1+rank_bad)) - log((double)(1+rank_corr))));
		float w=1.0;
		sample_pairs.push_back(weight_pair(corr_idx,badIndex,w));
		total_w += w;
	//	cout << "Added " << sample_pairs.size() << " :\t" << corr_idx << "\t" << badIndex << "\t" << w << endl;
	}
//	cout << endl;

	if (total_w>0)
	{
		float mult_val = (float)(sample_pairs.size())/total_w;
		for (i=0; i<sample_pairs.size(); i++)
			sample_pairs[i].weight*=mult_val;
	}
}






/*************************************************************************************
Creates the training data for models of de novo predictions (6-14).
These sequences do not have to reach the N- or C-terminals
**************************************************************************************/
void PeptideRankScorer::create_training_data_for_partial_denovo_ranking(
				const string& mgf_list,
				const double train_ratio,
				const int max_num_pairs,
				const int charge,
				const int size_idx,
				const float penalty_for_bad_aa,
				RankBoostDataset& train_ds, 
				RankBoostDataset& test_ds,
				char *test_scan_file,
				int   length_limit)
{	
/*	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(this->allScoreModelsPtr_);
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	const int max_num_ppp_frags = 4;
	const int max_num_pairs_per_spec = 80;
	const mass_t tolerance_diff = 0.275;

	if (model_type !=1 && model_type !=3)
	{
		cout << "Error: this training function is only intended for partial de novo samples!" << endl;
		cout << "Need to set model type to 1  or 3, not " << model_type << endl;
		exit(1);
	} 
	

	train_ds.clear();
	test_ds.clear();

	Config *config = allScoreModels->get_config();

	config->set_use_spectrum_charge(1);

	const vector< vector<mass_t> >& size_thresholds = peak_model->get_size_thresholds();

	mass_t min_mz=0;
	mass_t max_mz = 10000;
	if (size_idx>0)
		min_mz = (size_thresholds[charge][size_idx-1]+charge-1)/(mass_t)charge;
	if (size_idx< size_thresholds[charge].size())
		max_mz = (size_thresholds[charge][size_idx]+charge-1)/(mass_t)charge;

	int num_solutions= (length_limit<=8 ? 1000 :  2000);

	if (model_type == 3 && model_length<10)
	{
		if (model_length<= 4)
			num_solutions= 100;
		if (model_length == 5)
			num_solutions = 250;
		if (model_length > 5)
			num_solutions = 750;
	}

	if (length_limit>15)
		length_limit=15;

	if (charge >2 && length_limit>10)
		length_limit=10;


	cout << endl << "Model type " << model_type << ", model length " << model_length << endl;
	cout << "Length limit : " << length_limit << endl;
	cout << "Num solutions : " << num_solutions << endl;
//	max_mz = min_mz + 30;
	
	cout << "Charge " << charge << " size " << size_idx << endl;
	cout << "Min m/z " << setprecision(2) << fixed << min_mz  << "  Max m/z " << max_mz << endl;

	double start_time = time(NULL);

	// read spectra
	FileManager fm;
	FileSet	    fs;
	BasicSpecReader bsr;
	QCPeak		peaks[4000];

	fm.init_from_list_file(allScoreModels->get_config(),mgf_list.c_str());
	fs.select_files_in_mz_range(fm,min_mz,max_mz,charge);
	fs.randomly_reduce_ssfs(36000);
	
	int num_pairs_per_spec = (int)(max_num_pairs/(float)fs.get_total_spectra() + 0.5);
	if (num_pairs_per_spec>max_num_pairs_per_spec)
		num_pairs_per_spec=max_num_pairs_per_spec;

	cout << "Read " << fs.get_total_spectra() << " spectra headers" << endl;
	if (fs.get_total_spectra()<10)
	{
		cout << "Error: not enough spectra to train!" << endl;
		exit(1);
	}

	
	cout << "Training with " << num_pairs_per_spec << " pairs per specturm.." << endl;
	cout << "Creating RankBoostDatasets... proportion for training " << train_ratio << endl;

	vector<int> ppp_frag_type_idxs;
	ppp_frag_type_idxs.clear();

	cout << "Using a combined peak model" << endl;
	cout << "Selecting peaks of the paritally mobile model" << endl;

	ppp_frag_type_idxs = peak_model->get_model_ptr(charge,size_idx,PARTIALLYMOBILE)->get_fragment_type_idxs();
	int f;
	for (f=0; f<ppp_frag_type_idxs.size(); f++)
		cout << f+1 << "\t" << config->get_fragment(ppp_frag_type_idxs[f]).label << endl;
	cout << endl;

	sort(ppp_frag_type_idxs.begin(),ppp_frag_type_idxs.end());
	DeNovoPartitionModel *&part_model = dnv_part_models[charge][size_idx];
	part_model->init_features(model_type, charge,size_idx,ppp_frag_type_idxs,config);
		
	int num_groups_in_train=0;
	int num_groups_in_test=0;
	int num_train_pairs=0;

	ofstream test_scan_stream;
	if (test_scan_file)
	{
		test_scan_stream.open(test_scan_file);
		if (! test_scan_stream.is_open())
		{
			cout << "Error: couldn't open test stream!" << endl;
			exit(1);
		}
	}
	
	
	static vector<PrmGraph *> prm_ptrs;
	static vector<SeqPath> solutions;
	const vector<SingleSpectrumFile *>& all_ssfs = fs.get_ssf_pointers();
	int spec_idx;
	int num_spectra_read=0;

	int num_without_sol_in_graph=0;
	int num_without_correct_sol =0;
	int num_possible=0;

	config->set_use_spectrum_charge(1);

	for (spec_idx=0; spec_idx<all_ssfs.size(); spec_idx++)
	{
		MGF_single *ssf = (MGF_single *)all_ssfs[spec_idx];
		const Peptide& full_pep = ssf->peptide;
		const mass_t true_mass_with_19 = full_pep.get_mass_with_19();
		BasicSpectrum     bs;
		AnnotatedSpectrum as;
		PrmGraph prm;

		const vector<int>& aas = ssf->peptide.get_amino_acids();
		int a;
		for (a=0; a<aas.size(); a++)
			if (aas[a]>Val)
				break;
		if (a<aas.size())
			continue;

		num_possible++;
//
//		if (spec_idx>2000)
//			break;

		const int num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		bs.peaks = peaks;
		bs.num_peaks = num_peaks;
		bs.ssf = ssf;
		as.init_from_QCPeaks(config,peaks,num_peaks,ssf);
		
		prm.create_graph_from_spectrum(allScoreModels,&as,true_mass_with_19,charge);
		SeqPath longest_corr = prm.get_longest_subpath(full_pep,0);

	//	longest_corr.print_full(config);

		if (longest_corr.get_num_aa()<5 || 
			longest_corr.get_num_aa()<model_length)
			continue;

		const int num_aa_in_corr = longest_corr.get_num_aa();
		int max_length = num_aa_in_corr+3;
		int min_length = 6;

		if (model_type == 1)
		{
			if (max_length > 15)
				max_length = 15;
			if (max_length > length_limit)
				max_length = length_limit;
		}
		else
		{
			min_length=model_length;
			max_length=model_length;
		}

		const bool add_to_train = (myRandom() <= train_ratio);
		RankBoostDataset& rds = (add_to_train ? train_ds : test_ds);
		if (test_scan_file && ! add_to_train)
			test_scan_stream << ssf->file_idx << " " << ssf->idx_in_file << " " <<
			ssf->peptide.get_num_aas() << endl;

		clock_t c_start = clock();

		// generate de novo solutions
		vector<PmcSqsChargeRes> pmc_sqs_res;
		vector<mass_t> pms_with_19;
		vector<int>    charges;
		
	//	model->get_best_mz_charge(config,bs,&mz1,&charge1,&prob1,&mz2,&charge2,&prob2,&pmc_sqs_res);
	//	allScoreModels->select_pms_and_charges(config,bs,pms_with_19,charges,&pmc_sqs_res);

		SeqPath best_solution_in_a_real_graph;
		if (1)
		{
			const int iters = pms_with_19.size()-1;
			int j;
			(j=0; j<iters; j++)
				if (myRandom()<0.8)
				{
					pms_with_19.pop_back();
					charges.pop_back();
				}


		//	for (j=0; j<pms_with_19.size(); j++)
		//		cout << pms_with_19[j] << " " << charges[j] << endl;


			for (j=0; j<pms_with_19.size(); j++)
				if (fabs(pms_with_19[j]-true_mass_with_19)<0.275)
					break;

			if (j==pms_with_19.size())
			{
				pms_with_19.push_back(true_mass_with_19+(myRandom()-0.5)*0.3);
				charges.push_back(charge);
			//	cout << pms_with_19[pms_with_19.size()-1] << " " << charge << " **" <<endl;
			}
		
			if (prm_ptrs.size()<pms_with_19.size())
				prm_ptrs.resize(pms_with_19.size(),NULL);
		
			if (model_type == 1)
			{
				generate_denovo_solutions_from_several_pms(
					prm_ptrs,
					allScoreModels,
					&as,
					true, 
					num_solutions,
					min_length,
					max_length,
					pms_with_19,
					charges,
					solutions);
			}
			else
			{
				vector<int> num_tags;
				num_tags.resize(10,0);
				num_tags[model_length]=num_solutions;
				//generate_tags(prm_ptrs,allScoreModels,bs,&as,num_tags,model_length,
				//	pms_with_19,charges,solutions,true);
			}


			// find longest correct from the prm with the closest mass to the true mass
			int longest_aa=0;
			for (j=0; j<pms_with_19.size(); j++)
			{
				SeqPath longest = prm_ptrs[j]->get_longest_subpath(full_pep,0);
				const int num_best_aa = longest.get_num_aa();
				if (num_best_aa<min_length )
					continue;

				if (num_best_aa>longest_aa)
				{
					longest_aa = num_best_aa;
					best_solution_in_a_real_graph = longest;
				}
			}

			if (longest_aa==0)
			{
				num_without_sol_in_graph++;
				continue;
			}
		}

		vector<int> correct_solution_idxs;
		vector<int> incorrect_solution_idxs;
	
		correct_solution_idxs.clear();
		incorrect_solution_idxs.clear();
		int i;

		
		if (spec_idx % 100 == 0)
		{
			cout << "processed " << num_spectra_read << "/" << spec_idx << " : " << all_ssfs.size() << " spectra (";
			double curr_time = time(NULL);
			cout << curr_time - start_time << " secs.)" << endl;
		}
		

	//	cout << num_spectra_read << " (" << min_length << "-" << max_length << ")\t st:" << sam_time<< " \t";
	//	ssf->print_ssf_stats(config);
	//	output_denovo_solutions(ssf,config,cout,solutions,20);

		longest_corr.make_seq_str(config);
		for (i=0; i<solutions.size(); i++)
		{
			const SeqPath& path = solutions[i];
			const int num_corr_aa = path.get_num_correct_aas(full_pep,config);
			const int num_aa_in_path = path.get_num_aa();

			if (num_corr_aa == num_aa_in_path)
			{
				if (correct_solution_idxs.size() > 0)
					continue;

				correct_solution_idxs.push_back(i);

			
			}
			else
				incorrect_solution_idxs.push_back(i);
		}

		// don't use spectra without a good sequence in the returned solutions!
		if (correct_solution_idxs.size() == 0|| incorrect_solution_idxs.size()<5)
		{
			num_without_correct_sol++;
			continue;
		}
	
		vector<int> sol_sample_idxs_in_rds;
		sol_sample_idxs_in_rds.resize(solutions.size(),NEG_INF); // idxs of samples in the dataset

		int groupIndex;
		if (add_to_train) 
		{
			groupIndex = num_groups_in_train++;
		}
		else 
		{
			groupIndex = num_groups_in_test++;
		}

		// randomly select the pairs of samples

		vector<weight_pair> weight_pairs;
		select_sample_pairs(solutions,correct_solution_idxs,incorrect_solution_idxs,
			weight_pairs,num_pairs_per_spec);

		for (i=0; i<weight_pairs.size(); i++)
		{
			const int corr_sol_idx = weight_pairs[i].idx_corr;
			const int bad_sol_idx  = weight_pairs[i].idx_bad;
		
			// add good de novo sample if needed
			if (sol_sample_idxs_in_rds[corr_sol_idx]<0)
			{
				const SeqPath& seq_path = solutions[corr_sol_idx];
				PeptideSolution sol;
				vector<int> amino_acids;
				seq_path.get_amino_acids(amino_acids);
				sol.charge = charge;
				sol.pep.set_peptide_aas(amino_acids);
				sol.pep.set_n_gap(seq_path.n_term_mass);
				sol.reaches_n_terminal = (sol.pep.get_n_gap()<1.0);
				sol.pep.calc_mass(config);
				sol.reaches_c_terminal = (sol.pep.get_mass_with_19()>seq_path.pm_with_19 - 9);
				sol.pm_with_19 = (sol.reaches_n_terminal && sol.reaches_c_terminal ? 
					sol.pep.get_mass_with_19() : seq_path.pm_with_19 );
				sol.type = -1;

		
				RankBoostSample corr_rbs;

				if (model_type == 1)
				{
					//fill_denovo_peptide_rbs(sol, seq_path, peaks, num_peaks, as, pmc_sqs_res, corr_rbs, size_idx);
				}
				else
				{
					//fill_tag_rbs(sol,seq_path,peaks,num_peaks,as,corr_rbs,size_idx);
				}

				sol_sample_idxs_in_rds[corr_sol_idx] = rds.get_num_samples();

				corr_rbs.groupIndex=groupIndex;
				corr_rbs.tag1=NEG_INF;
				corr_rbs.rank_in_group=0;
				rds.add_sample(corr_rbs);
			}

			// add bad de novo sample if needed
			if (sol_sample_idxs_in_rds[bad_sol_idx]<0)
			{
				RankBoostSample bad_rbs;
				const SeqPath& seq_path = solutions[bad_sol_idx];
				PeptideSolution sol;
				vector<int> amino_acids;
				seq_path.get_amino_acids(amino_acids);
				sol.charge = charge;
				sol.pep.set_peptide_aas(amino_acids);
				sol.pep.set_n_gap(seq_path.n_term_mass);
				sol.reaches_n_terminal = (sol.pep.get_n_gap()<1.0);
				sol.pep.calc_mass(config);
				sol.reaches_c_terminal = (sol.pep.get_mass_with_19()>seq_path.pm_with_19 - 9);

				sol.pm_with_19 = (sol.reaches_n_terminal && sol.reaches_c_terminal ? 
					sol.pep.get_mass_with_19() : seq_path.pm_with_19 );

				sol.type = -1;

				if (model_type == 1)
				{
					//fill_denovo_peptide_rbs(sol, seq_path, peaks, num_peaks, as, pmc_sqs_res, bad_rbs, size_idx);
				}
				else
				{
					//fill_tag_rbs(sol,seq_path,peaks,num_peaks,as,bad_rbs,size_idx);
				}
				
				sol_sample_idxs_in_rds[bad_sol_idx] = rds.get_num_samples();
				
				bad_rbs.groupIndex=groupIndex;
				bad_rbs.tag1=NEG_INF;
				bad_rbs.rank_in_group=1;
				rds.add_sample(bad_rbs);
			}
			else
				continue;


			if (sol_sample_idxs_in_rds[bad_sol_idx]<0 || sol_sample_idxs_in_rds[corr_sol_idx]<0)
			{
				int qq=1;
				cout << "Bad:  " << bad_sol_idx << " " << sol_sample_idxs_in_rds[bad_sol_idx] << endl;
				cout << "Good: " << corr_sol_idx << " " << sol_sample_idxs_in_rds[corr_sol_idx] << endl;
				cout << "Oy.." << endl;
			}
		
			rds.add_to_phi_vector(sol_sample_idxs_in_rds[bad_sol_idx],
								  sol_sample_idxs_in_rds[corr_sol_idx],0,weight_pairs[i].weight);

		}


		clock_t c_end = clock();
		double sam_time = (c_end - c_start)/(double)CLOCKS_PER_SEC;

		if (sam_time>15.0)
		{
			ssf->print_ssf_stats(config);
			cout << "Warning: spectrum required " << sam_time << " to process (skipping 20)." << endl;
			spec_idx+=20;
		}

		num_spectra_read++;
	}

	cout << endl << "collected results from " << num_spectra_read << " spectra ( num possible " << 
		num_possible << " )" << endl;

	cout << "# without solution in graph = " << num_without_sol_in_graph << " (" << fixed << setprecision(3) <<
		num_without_sol_in_graph/(float)num_possible << ") "<< endl;
	cout << "# without good solution in path set = " << num_without_correct_sol << " (" <<
		num_without_correct_sol/(float)num_possible << ") " << endl;

	train_ds.set_num_groups(num_groups_in_train);
	test_ds.set_num_groups(num_groups_in_test);
			
	cout << "computing phi weights..." << endl;
	train_ds.compute_total_phi_weight();
	
	cout << "initializing potential lists..." << endl;
	train_ds.initialize_potenital_lists();
	
	cout << "initializing real feature table..." << endl;
	train_ds.initialzie_real_feature_table(part_model->get_feature_names().size());

	cout << "Processed " << num_spectra_read << " spectra" << endl;
	cout << num_groups_in_train << " sets in training set (" << train_ds.get_phi_support().size() << " pairs)" << endl;
	cout << num_groups_in_test  << " sets in test set     (" << test_ds.get_phi_support().size() << " pairs)" << endl << endl;
	*/
}




/*************************************************************************
**************************************************************************/
void DeNovoPartitionModel::simple_print_peak_pairs(
						  Config *config,
						  const vector<idx_weight_pair>& pair_idxs, 
						  vector<string>* peptide_strings,
						  const RankBoostDataset& ds,
						  int max_examples,
						  ostream& os) const
{

	const int max_examples_for_set = 3;
	const vector<SamplePairWeight>& phi_support= ds.get_phi_support();
	vector<int> set_idx_counts;

	set_idx_counts.resize(pair_idxs.size(),0);

	int p_idx;
	int counter = 0;
	for (p_idx=0; p_idx<pair_idxs.size(); p_idx++)
	{
		const int pair_idx = pair_idxs[p_idx].idx;
		const int x0_idx = phi_support[pair_idx].idx1;
		const int x1_idx = phi_support[pair_idx].idx2;
		const int pair_type = phi_support[pair_idx].tag;
		const RankBoostSample sam_x0 = ds.get_sample(x0_idx);
		const RankBoostSample sam_x1 = ds.get_sample(x1_idx);
		const int set_idx = sam_x0.groupIndex;

		if (set_idx>=set_idx_counts.size())
			set_idx_counts.resize(set_idx*2+10000,0);

		if (set_idx_counts[set_idx]>=max_examples_for_set)
				continue;
	
		set_idx_counts[set_idx]++;
		counter++;
	
		os << counter << "\t" << phi_support[pair_idx].weight << "\t" << pair_idxs[p_idx].weight << endl;

		// print example
		const RankBoostSample& corr_sam =  ds.get_sample(x0_idx);
		const RankBoostSample& bad_sam  =  ds.get_sample(x1_idx);

		if (corr_sam.groupIndex!=bad_sam.groupIndex && pair_type != SOL_INCORRECT_DB_CROSS)
		{
			cout << "Error: mismatch in correct and bad sample group idxs!" << endl;
			exit(1);
		}
		
		os << "TRUE:\t" << (*peptide_strings)[sam_x1.tag1] << endl;
		os << "MISS:\t" << (*peptide_strings)[sam_x0.tag1] << "\t"
			<< sam_x0.tag2 << "\t" << "(" << sam_x0.tag3 << ")" << endl;


		os << endl << endl;

		if (max_examples>0 && counter>=max_examples)
			break;
	}
}


/***************************************************************************************
The main model training function for complete sequences (model type 0 - db results or 
model tpye 2 - full denovo results).
****************************************************************************************/
void PeptideRankScorer::train_partition_model_for_complete_sequences(
					const string& db_dir,
					const string& correct_dir,
					const string& denovo_dir,
					const string& mgf_list,
					char *report_dir,
					char *name,
					int charge,
					int size_idx,
					int max_num_rounds,
					float max_boost_ratio,
					int	  max_num_samples,
					float ratio_pair_db,
					float ratio_pair_denovo,
					float ratio_pair_db_cross,
					char  *rerank_path,
					int	  rerank_depth)
{
	char report_buff[512],test_scan_buff[512];
	char *report_prefix=NULL;
	if (report_dir)
	{
		sprintf(report_buff,"%s/%s_%d_%d",report_dir, name ,charge, size_idx);
		report_prefix = report_buff;
	}

	bool only_set_shifts = false;

	// init models

	AllScoreModels* allScoreModels = new AllScoreModels;
	allScoreModelsPtr_ = static_cast<void*>(allScoreModels);

	PeptideCompAssigner& comp_assigner = allScoreModels->getPeptideCompositionAssigner();

	allScoreModels->read_model("CID_IT_TRYP");
	Config *config = allScoreModels->get_config();
	config->apply_selected_PTMs("M+16:Q-17:C+57");
	allScoreModels->read_rank_models("CID_IT_TRYP");

	PeakRankModel*& peakPeredictionModel = allScoreModels->get_peak_prediction_model_ptr(model_type);

	if (! peakPeredictionModel)
	{
		peakPeredictionModel = new PeakRankModel;
		if (! peakPeredictionModel->read_peak_rank_model(config,"DBC4_PEAK/DBC4",true,
			charge,size_idx))
		{
			cout << "Error: couldn't read peak prediction model!" << endl;
			exit(1);
		}
	}

	PeakRankModel *&peak_model = peakPeredictionModel;

	comp_assigner.read_and_init_from_tables(config,"LTQ_COMP/IT_TRYP");

	if (only_set_shifts)
	{
		PeptideRankScorer *drs = (PeptideRankScorer *)allScoreModels->get_rank_model_ptr(0);
		dnv_part_models = drs->dnv_part_models;
	}
	else
	{
		if (dnv_part_models.size()<=charge)
			init_tables();
	}

	RankBoostDataset train_ds, test_ds;

	if (! dnv_part_models[charge][size_idx])
		dnv_part_models[charge][size_idx] = new DeNovoPartitionModel;

	vector<string> peptide_strings;
	char *test_scan_file = NULL;
	if (report_prefix)
	{
		sprintf(test_scan_buff,"%s_test_scans.txt",report_prefix);
		test_scan_file = test_scan_buff;
	}

	if (model_type == 0)
	{
		create_training_data_for_complete_sequence_ranking(db_dir,correct_dir,
			denovo_dir, mgf_list, 0.7, max_num_samples, charge, size_idx, train_ds, test_ds,
			&peptide_strings, test_scan_file,0.2,0.0,0.8);
	}
	else if (model_type == 2)
	{
		create_training_data_for_complete_denovo_ranking(db_dir, correct_dir, mgf_list, 0.7, max_num_samples,
			charge, size_idx, train_ds, test_ds, &peptide_strings, test_scan_file, 0.8, rerank_path, rerank_depth);
	}

	

	cout << "READ " << peptide_strings.size() << " peptide examples" << endl;

	if (train_ds.get_num_samples() == 0 || peptide_strings.size() == 0)
	{
		cout << "Error: no trianing samples created for model " << report_prefix << endl;
		exit(1);
	}

	// init boost model for training
	RankBoostModel& boost = dnv_part_models[charge][size_idx]->boost_model;
	const vector<string>& feature_names = dnv_part_models[charge][size_idx]->get_feature_names();

	vector<string> empty;
	boost.init_rankboost_model_feature_names(empty,feature_names);

	train_ds.initialzie_real_feature_table(feature_names.size());

	train_ds.set_max_ratio_for_regular_update(max_boost_ratio);
				
	boost.init_rankboost_model_for_training(train_ds,150,25);
		
	train_ds.initialize_real_vote_lists(boost);

	if (report_prefix)
	{
		char name_buff[512];
		sprintf(name_buff,"%s_feature_summary.txt",report_prefix);
		ofstream sum_stream(name_buff);
		boost.summarize_features(train_ds.get_samples(),sum_stream);
		sum_stream.close();
	}
	else
		boost.summarize_features(train_ds.get_samples());

	vector<idx_weight_pair> miss_pairs;
	vector<string> header_strings;
	dnv_part_models[charge][size_idx]->write_denovo_partition_model_header_to_strings(model_type, header_strings);

	if (! only_set_shifts)
	{
		boost.train_rankboost_model(train_ds, max_num_rounds, &miss_pairs, &test_ds,
			-1, report_prefix, NULL, &header_strings);
	}

	// normalize score if needed
	if (model_type == 0)
		dnv_part_models[charge][size_idx]->set_shifts_and_scales_for_db(config,train_ds,peptide_strings);
	

	// final report
	if (report_dir)
	{
		if (! only_set_shifts)
		{
			char name_buff[512];
			sprintf(name_buff,"%s_train_miss_pairs.txt",report_prefix);
			ofstream report_stream(name_buff);
			if (! report_stream.is_open() || ! report_stream.good())
			{
				cout << "Error: couldn't open pairs report file for writing:" << name_buff << endl;
				exit(1);
			}

			dnv_part_models[charge][size_idx]->simple_print_peak_pairs(config,miss_pairs,
				&peptide_strings, train_ds,200, report_stream);
			report_stream.close();
		}

		string path = string(report_prefix) + "_model.txt";
		dnv_part_models[charge][size_idx]->write_denovo_partition_model(this->model_type, path.c_str());
	}
	else
		dnv_part_models[charge][size_idx]->simple_print_peak_pairs(config,miss_pairs,
			&peptide_strings, train_ds,100);

	cout << "DONE..." << endl;
}


/***********************************************************************
// for de novo predictions that might be incomplete
************************************************************************/
void PeptideRankScorer::train_partial_denovo_partition_model(
					const string& mgf_list,
					char *report_dir,
					char *name,
					int charge,
					int size_idx,
					int max_num_rounds,
					float max_boost_ratio,
					int   max_num_samples,
					int   length_limit,
					char  *rerank_path)
{
	char report_buff[512],test_scan_buff[512];
	char *report_prefix=NULL;
	if (report_dir)
	{
		sprintf(report_buff,"%s/%s_%d_%d",report_dir, name ,charge, size_idx);
		report_prefix = report_buff;
	}


	// init models
	AllScoreModels* allScoreModels =	new AllScoreModels;
	allScoreModelsPtr_ = static_cast<AllScoreModels*>(allScoreModels);
	
	PeptideCompAssigner& compAssigner = allScoreModels->getPeptideCompositionAssigner();

	allScoreModels->read_model("CID_IT_TRYP");
	Config *config = allScoreModels->get_config();
	config->apply_selected_PTMs("C+57");
	allScoreModels->read_rank_models("CID_IT_TRYP");
	
	PeakRankModel*& peakPeredictionModel = allScoreModels->get_peak_prediction_model_ptr(model_type);

	if (! peakPeredictionModel)
	{
		peakPeredictionModel	  = new PeakRankModel;
		if (! peakPeredictionModel->read_peak_rank_model(config,"ITDNV_PEAK/ITDNV4",
			true,charge,size_idx))
		{
			cout << "Error: couldn't read peak model " << "ITDNV_PEAK/ITDNV4" << endl;
			exit(1);

		}
	}
	PeakRankModel *&peak_model = peakPeredictionModel;


	compAssigner.read_and_init_from_tables(config,"LTQ_COMP/IT_TRYP");

	if (dnv_part_models.size()<=charge)
		init_tables();

	RankBoostDataset train_ds, test_ds;

	if (! dnv_part_models[charge][size_idx])
		  dnv_part_models[charge][size_idx] = new DeNovoPartitionModel;

	char *test_scan_file = NULL;
	if (report_prefix)
	{
		sprintf(test_scan_buff,"%s_test_scans.txt",report_prefix);
		test_scan_file = test_scan_buff;
	}

	create_training_data_for_partial_denovo_ranking(mgf_list, 0.667, max_num_samples, 
		charge, size_idx, -1.0,train_ds, test_ds, test_scan_file, length_limit);

	
	// init boost model for training
	RankBoostModel& boost = dnv_part_models[charge][size_idx]->boost_model;
	const vector<string>& feature_names = dnv_part_models[charge][size_idx]->get_feature_names();

	vector<string> empty;

	boost.init_rankboost_model_feature_names(empty,feature_names);

	train_ds.initialzie_real_feature_table(feature_names.size());

	train_ds.set_max_ratio_for_regular_update(max_boost_ratio);
				
	boost.init_rankboost_model_for_training(train_ds,100,50);
		
	train_ds.initialize_real_vote_lists(boost);

	if (report_prefix)
	{
		char name_buff[512];
		sprintf(name_buff,"%s_feature_summary.txt",report_prefix);
		ofstream sum_stream(name_buff);
		boost.summarize_features(train_ds.get_samples(),sum_stream);
		sum_stream.close();
	}
	else
		boost.summarize_features(train_ds.get_samples());

	vector<idx_weight_pair> miss_pairs;
	vector<string> header_strings;
	dnv_part_models[charge][size_idx]->write_denovo_partition_model_header_to_strings(model_type, header_strings);

	boost.train_rankboost_model(train_ds, max_num_rounds, &miss_pairs, &test_ds,
		-1, report_prefix, NULL, &header_strings);
	

	// final report
	if (report_dir)
	{
		char name_buff[512];
	
		sprintf(name_buff,"%s_feature_list.txt",report_prefix);
		ofstream feature_stream(name_buff);
		if (! feature_stream.is_open() || ! feature_stream.good())
		{
			cout << "Error: couldn't feature_stream file for writing:" << name_buff << endl;
			exit(1);
		}
		boost.ouput_importance_ranked_feature_list(train_ds,feature_stream);
		feature_stream.close();

		string path = string(report_prefix) + "_model.txt";
		dnv_part_models[charge][size_idx]->write_denovo_partition_model(this->model_type, path.c_str());
	}


	cout << "DONE..." << endl;
}




