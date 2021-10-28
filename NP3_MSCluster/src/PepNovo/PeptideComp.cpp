#include "PeptideComp.h"


void PeptideCompStats::clear_pcs()
{
	num_aa=0;
	int i;
	
	for (i=0; i<=MAX_COMP_LEN; i++)
	{
		start_comp[i]=0;
		end_comp[i]=0;

		int j;
		for (j=0; j<=MAX_COMP_CAT; j++)
			cat_counts[i][j]=0;
	}
}

void PeptideCompStats::print_pcs(ostream& os ) const
{
	os << "NAA: " << num_aa << endl;
	os << "START:" ;
	int i;
	for (i=1; i<= MAX_COMP_LEN; i++)
		os << "\t" << start_comp[i];
	os << endl;
	os << "END:" ;
	for (i=1; i<= MAX_COMP_LEN; i++)
		os << "\t" << end_comp[i];
	os << endl;
	for (i=1; i<= MAX_COMP_LEN; i++)
	{
		os << i << " -\t";
		int j;
		for (j=1; j<=MAX_COMP_CAT; j++)
			os << " " << cat_counts[i][j];
		os << endl;
	}
}

void PeptideCompAssigner::init_aa_translations()
{
	const vector<int>& session_aas = config->get_session_aas();
	const vector<int>& org_aas = config->get_org_aa();

	aa_translation.resize(config->get_max_session_aa_idx()+1,NEG_INF);
	aa_translation[Ala]=1;
	aa_translation[Arg]=2;
	aa_translation[Asn]=3;
	aa_translation[Asp]=4;
	aa_translation[Cys]=5;
	aa_translation[Gln]=6;
	aa_translation[Glu]=7;
	aa_translation[Gly]=8;
	aa_translation[His]=9;
	aa_translation[Ile]=10;
	aa_translation[Leu]=10;
	aa_translation[Lys]=11;
	aa_translation[Met]=12;
	aa_translation[Phe]=13;
	aa_translation[Pro]=14;
	aa_translation[Ser]=15;
	aa_translation[Thr]=16;
	aa_translation[Trp]=17;
	aa_translation[Tyr]=18;
	aa_translation[Val]=19;

	const int max_aa_idx=config->get_max_session_aa_idx();
	int i;
	for (i=0; i<session_aas.size(); i++)
	{
		aa_translation[session_aas[i]]=aa_translation[org_aas[session_aas[i]]];
	}

//	for (i=0; i<session_aas.size(); i++)
//		cout << config->get_aa2label()[session_aas[i]] << "\t" << aa_translation[session_aas[i]] << endl;
}

void PeptideCompAssigner:: read_table_to_vector(char *file_path, int naa, vector<int>& vec)
{
	const vector<int>& char2aa = config->get_char2aa();
	const vector<int>& org_aas = config->get_org_aa();
	if (naa<1 || naa>3)
	{
		cout << "Error: number of aas should be 1-3!" << endl;
		exit(1);
	}
	int max_size = 20;
	if (naa==2)
		max_size = 400;
	if (naa==3)
		max_size = 8000;

	vec.resize(max_size,NEG_INF);

	ifstream ifs(file_path);
	if (! ifs.is_open())
	{
		cout << "Error: couldn't open file for reading: " << file_path << endl;
		exit(1);
	}

	char buff[64];
	while (! ifs.eof())
	{
		ifs.getline(buff,64);
		if (ifs.gcount()<5)
			continue;

		char sym[4];
		int  trans[3];
		int category=NEG_INF;
		int count=NEG_INF;
		if (sscanf(buff,"%s\t%d\t%d",sym,&category,&count) != 3)
		{
			cout << "Error: bad line in file " << file_path << " :" << endl << buff << endl;
			exit(1);
		}

		int i;
		for (i=0; i<naa; i++)
		{
			const char org_char = sym[i];
			const int  org_aa   = org_aas[char2aa[org_char]];
			trans[i]=org_aa;
		}

		const int code = this->calc_aa_code(trans,naa);
		vec[code]=category;
	}
	ifs.close();
}




void PeptideCompAssigner::read_and_init_from_tables(Config *_config, const char *name)
{
	config = _config;

	model_name = name;

	init_aa_translations();

	start_assigns.resize(MAX_COMP_LEN+1);
	end_assigns.resize(MAX_COMP_LEN+1);
	mid_assigns.resize(MAX_COMP_LEN+1);

	int i;
	for (i=1; i<=MAX_COMP_LEN; i++)
	{
		char file_path[256];

		sprintf(file_path,"%s/%s_start_freq_%d.txt",config->get_resource_dir().c_str(),name,i);
		read_table_to_vector(file_path,i,start_assigns[i]);

		sprintf(file_path,"%s/%s_end_freq_%d.txt",config->get_resource_dir().c_str(),name,i);
		read_table_to_vector(file_path,i,end_assigns[i]);

		sprintf(file_path,"%s/%s_mid_freq_%d.txt",config->get_resource_dir().c_str(),name,i);
		read_table_to_vector(file_path,i,mid_assigns[i]);
	}

	was_initialized = true;
}


void PeptideCompAssigner::fill_peptide_stats(const Peptide& peptide, 
											 PeptideCompStats& stats) const
{

	const vector<int>& org_aa = config->get_org_aa();
	vector<int> aas = peptide.get_amino_acids();
	const int num_aa = aas.size();

	int i;
	for (i=0; i<aas.size(); i++)
		aas[i]=org_aa[aas[i]];

	stats.clear_pcs();
	stats.num_aa = num_aa;

	if (num_aa<3)
	{
		cout << "Error: trying to fill peptide with " << num_aa << " amino acids!" << endl;
		exit(1);
	}

	
	for (i=1; i<=MAX_COMP_LEN; i++)
	{
		stats.start_comp[i]=start_assigns[i][calc_aa_code(&aas[0],i)];
		stats.end_comp[i]=end_assigns[i][calc_aa_code(&aas[num_aa-i],i)];
		int j;
		for (j=1; j<num_aa-i; j++)
		{
			const int cat = mid_assigns[i][calc_aa_code(&aas[j],i)];
			stats.cat_counts[i][cat]++;
		}
	}
}

int PeptideCompAssigner::get_aa_category(int num_aa, const int *aas, bool n_term, bool c_term) const
{
	if (! this->was_initialized)
	{
		cout << "Error: using an uninitialized peptide composition asssigner!" << endl;
		exit(1);
	}

	if (num_aa<1)
		return 0;

	if (num_aa>3)
	{
		if (c_term)
			aas+=num_aa-3;
		num_aa=3;
	}
	int code=calc_aa_code(aas,num_aa);

	if (n_term)
		return start_assigns[num_aa][code];

	if (c_term)
		return end_assigns[num_aa][code];

	return mid_assigns[num_aa][code];
}


