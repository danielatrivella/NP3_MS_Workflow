#include "PeakRankModel.h"
#include "AllScoreModels.h"
#include "PrmGraph.h"
#include "PepNovo_auxfun.h"

const char* mobility_labels[]={"","Mobile","Partially Mobile","Non-Mobile"};

extern const int num_RKH_combos;
extern const int num_RKH_pairs;
extern const string combos[];






/**********************************************************************
Categorizes the peptides as mobile,partially mobile, and non-mobile
according to Wysocki 2000, and Kapp 2003.
***********************************************************************/
int get_proton_mobility(const Peptide& pep, int charge)
{
	const vector<int>& amino_acids = pep.get_amino_acids();
	size_t i;

	int num_arg=0;
	int num_his=0;
	int num_lys=0;
	for (i=0; i<amino_acids.size(); i++)
	{
		if (amino_acids[i]==Arg)
			num_arg++;
		if (amino_acids[i]==Lys)
			num_lys++;
		if (amino_acids[i]==His)
			num_his++;
	}
	if (num_arg+num_his+num_lys<charge)
		return MOBILE;
	if (num_arg>=charge)
		return NONMOBILE;
	return PARTIALLYMOBILE;
}


/**************************************************************
***************************************************************/
int PeptideSolution::calc_mobility() const
{
	const vector<int>& amino_acids = pep.get_amino_acids();
	int numH=0,numK=0,numR=0;
	int i;
	for (i=0; i<amino_acids.size(); i++)
	{
		if (amino_acids[i] == His)
			numH++;

		if (amino_acids[i] == Lys)
			numK++;

		if (amino_acids[i] == Arg)
			numR++;
	}
	if (this->most_basic_aa_removed_from_n == His)
		numH++;
	if (this->most_basic_aa_removed_from_n == Lys)
		numK++;
	if (this->most_basic_aa_removed_from_n == Arg)
		numR++;
	if (this->most_basic_aa_removed_from_c == His)
		numH++;
	if (this->most_basic_aa_removed_from_c == Lys)
		numK++;
	if (this->most_basic_aa_removed_from_c == Arg)
		numR++;

	return get_proton_mobility(this->charge,numR,numH,numK);
	
}

int get_proton_mobility(int charge, int num_arg, int num_his, int num_lys)
{
	if (num_arg+num_his+num_lys<charge)
		return MOBILE;
	if (num_arg>=charge)
		return NONMOBILE;
	return PARTIALLYMOBILE;
}


struct inten_pair {
	inten_pair() : idx(int(NEG_INF)), inten(int(NEG_INF)) {};
	inten_pair(int _i, float _n) : idx(_i), inten(_n) {};
	bool operator< (const inten_pair& other) const
	{
		return inten>other.inten;
	}
	int idx;
	float inten;
};

/***********************************************************************
calculates the relative rank of the cuts. If intensity is 0, the rank 
999 is given (rank[0] also is 999)
************************************************************************/
void calc_cut_ranks(const vector<float>& intens, vector<int>& cut_ranks)
{
	vector<inten_pair> pairs;
	int i;
	for (i=1; i<intens.size(); i++)
		if (intens[i]>0)
			pairs.push_back(inten_pair(i,intens[i]));
	sort(pairs.begin(),pairs.end());
	
	cut_ranks.clear();
	cut_ranks.resize(intens.size(),999);
	for (i=0; i<pairs.size(); i++)
		cut_ranks[pairs[i].idx]=i;
}

void convert_seq_path_to_peptide_soluition(Config *config, const SeqPath& seq_path, PeptideSolution& sol)
{
	vector<int> aas;
	seq_path.get_amino_acids(aas);

	sol.pep.set_peptide_aas(aas);
	sol.pep.set_n_gap(seq_path.n_term_mass);
	sol.pep.calc_mass(config);

	sol.pm_with_19 = seq_path.pm_with_19;
	sol.charge = seq_path.charge;
	if (sol.charge==0)
	{
		cout << "Error: must supply seq path with charge other than 0!" << endl;
		exit(1);
	}
	sol.reaches_n_terminal = (sol.pep.get_n_gap()<1.0);
	sol.reaches_c_terminal = (sol.pep.get_mass_with_19() + 15.0 > sol.pm_with_19);

	if (sol.reaches_n_terminal && sol.reaches_c_terminal)
		sol.pm_with_19 = sol.pep.get_mass_with_19();
}

void convert_seq_path_to_peptide_soluition_and_fill_in_aas(const Config *config, 
														   const Peptide& correct_pep,
														   const SeqPath& seq_path, 
														   PeptideSolution& sol)
{
	const mass_t true_mass_with_19 = correct_pep.get_mass_with_19();
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const vector<int>& correct_aas = correct_pep.get_amino_acids();
	vector<int> full_set_of_aas;
	vector<int> aas;
	seq_path.get_amino_acids(aas);

	full_set_of_aas.clear();
	if (seq_path.n_term_mass>1.0)
	{
		mass_t n_mass=0;
		int i;
		for (i=0; i<correct_aas.size(); i++)
		{
			if (fabs(n_mass-seq_path.n_term_mass)<3)
				break;
			n_mass+=aa2mass[correct_aas[i]];
			full_set_of_aas.push_back(correct_aas[i]);
		}
		if (i==correct_aas.size())
		{
			cout << "Error: seq_path doesn't have correct starting node!" << endl;
			cout << correct_pep.as_string(config) << endl;
			seq_path.print_full(config);
			exit(1);
		}
	}

	int i;
	for (i=0; i<aas.size(); i++)
		full_set_of_aas.push_back(aas[i]);

	if (fabs(seq_path.c_term_mass+19-true_mass_with_19)>5)
	{
		mass_t n_mass=0;
		int i;
		for (i=0; i<correct_aas.size(); i++)
		{
			if (fabs(n_mass-seq_path.c_term_mass)<5)
				break;
			n_mass+=aa2mass[correct_aas[i]];
		}
		if (i==correct_aas.size())
		{
			cout << "Error: seq_path doesn't have correct ending node!" << endl;
			exit(1);
		}

		for ( ; i<correct_aas.size(); i++)
			full_set_of_aas.push_back(correct_aas[i]);
	}
	

	sol.pep.set_peptide_aas(full_set_of_aas);
	sol.pep.set_n_gap(0);
	sol.pep.calc_mass(config);


	sol.charge = seq_path.charge;
	if (sol.charge==0)
	{
		cout << "Error: must supply seq path with charge other than 0!" << endl;
		exit(1);
	}
	sol.reaches_n_terminal = true;
	sol.reaches_c_terminal = true;
	sol.pm_with_19 = sol.pep.get_mass_with_19();
}


void  PeakRankModel::init_peak_rank_model_with_defaults(Config *_config, char *name,
														int feature_type)
{
	feature_set_type = feature_type;

	config = _config;

	set_peak_rank_model_name((name ? string(name) : "Rank"));

	set_size_thresholds();

	set_mass_detection_defaults();

	init_aa_defaults();

	convert_session_aas_to_model_aas();

	if (feature_set_type == 0)
	{
	//	set_binary_feature_names();

	//	set_real_feature_names();
		binary_feature_names.clear();
		set_simple_feature_names();
	}
	else if (feature_set_type == 1)
	{
		binary_feature_names.clear();
		this->set_advanced_feature_names();
	}
	else if (feature_set_type == 2)
	{
		binary_feature_names.clear();
		set_partial_denovo_feature_names();
	}
	else if (feature_set_type == 3)
	{
		binary_feature_names.clear();
		real_feature_names.clear();
	}
	else if (feature_set_type == 4)
	{
		binary_feature_names.clear();
		real_feature_names.clear();
	}
	else if (feature_set_type == 5)
	{
		binary_feature_names.clear();
		real_feature_names.clear();
	}
	else
	{
		cout << "Error: feture set type " << feature_set_type << " not supported!" << endl;
		exit(1);
	}

}

/******************************************************************************
*******************************************************************************/
void PeakRankModel::calc_peptide_predicted_scores(
							  const PeptideSolution& sol,
							  PeptidePeakPrediction& ppp,
							  int specific_size,
							  const vector<int>* ptr_frag_type_idxs) const
{
	const Peptide& pep = sol.pep;
	const mass_t pm_with_19 = sol.pm_with_19;
	const int spec_charge = sol.charge;
	const int mobility = get_proton_mobility(pep,spec_charge);
	const int size_idx = (specific_size>=0 ? specific_size : get_size_group(spec_charge,pm_with_19));
	
	if (! partition_models[spec_charge][size_idx][mobility])
	{
		cout << "Error: no rank partition model for " <<
			spec_charge << " " << size_idx << " " << mobility << endl;
		exit(1);
	}

	const mass_t min_detected_mass = calc_min_detected_mass(pm_with_19, spec_charge);
	const mass_t max_detected_mass = get_max_detected_mass();

	partition_models[spec_charge][size_idx][mobility]->calc_peptides_peaks_rank_scores(this, sol, 
		min_detected_mass, max_detected_mass, ppp, feature_set_type, ptr_frag_type_idxs);
}




bool PeakRankModel::has_intialized_models(int charge, int size_idx, int frag_idx) const
{
	int i;
	for (i=MOBILE; i<=NONMOBILE; i++)
	{
		if (! partition_models[charge][size_idx][i] ||
			! partition_models[charge][size_idx][i]->frag_models[frag_idx].get_ind_was_initialized())
			return false;
	}
	return true;
}



void PeptidePeakPrediction::make_rank_tables(const vector< vector<float> >& intens,
		vector< vector<int> >& observed_ranks, vector< vector<int> >& predicted_ranks) const
{
	const int max_num = rank_scores.size();
	predicted_ranks.resize(max_num);
	observed_ranks.resize(max_num);
	int f;
	for (f=0; f<max_num; f++)
	{
		if (rank_scores[f].size()>0)
		{
			find_ranks(rank_scores[f],predicted_ranks[f]);
			find_ranks(intens[f],observed_ranks[f]);
		}
		else
		{
			predicted_ranks[f].clear();
			observed_ranks[f].clear();
		}
	}
}



void PeptidePeakPrediction::print_ranks_vs_intens(const vector< vector<float> >& intens) const
{
	cout << setprecision(1) << fixed;
	
	cout << "Observed intensities:" << endl;
	int f,i;
	cout <<" ";
	for (i=1; i<intens[0].size(); i++)
		cout << "\t" << i;
	cout << endl;
	for (f=0; f<intens.size(); f++)
	{
		if (intens[f].size()==0)
			continue;
		cout << f << "\t";
		int i;
		for (i=1; i<intens[f].size(); i++)
		{
			if (intens[f][i]>0)
				cout << intens[f][i];
			cout << "\t";
		}
		cout << endl;
	}
	cout << endl;

	cout << "Predicted scores:" << endl;
	cout <<" ";
	for (i=1; i<rank_scores[0].size(); i++)
		cout << "\t" << i;
	cout << endl;
	for (f=0; f<rank_scores.size(); f++)
	{
		if (rank_scores[f].size()==0)
			continue;
		cout << f << "\t";
		int i;
		for (i=1; i<rank_scores[f].size(); i++)
		{
			if (rank_scores[f][i]>NEG_INF)
				cout << rank_scores[f][i];
			cout << "\t";
		}
		cout << endl;
	}
	cout << endl;


	vector< vector<int> > obs_table, pred_table;
	make_rank_tables(intens,obs_table,pred_table);

	cout << "Rank comparison:" << endl;
	cout <<" ";
	for (i=1; i<intens[0].size()-3; i++)
		cout << "\t" << i;
	cout << endl;
	for (f=0; f<pred_table.size(); f++)
	{
		if (intens[f].size()==0)
			continue;

		cout << f << "\t";
		int i;
		for (i=1; i<pred_table[f].size(); i++)
		{
			if (pred_table[f][i]>=0)
			{
				cout << obs_table[f][i] << ":" << pred_table[f][i] << "\t";
			}
			else
				cout << "\t";
		}
		cout << endl;
	}
	cout << endl;
	
}


/********************************************************************
Sets the default aa labels that are used in the model.
The aas that are included are the typical 20 amino acids minus Ile.
There are also M+16 and Q-17 that are included.
This list can be overriden if read from the model file.
*********************************************************************/
void PeakRankModel::init_aa_defaults()
{
	int aa;
	const vector<string>& aa2label =  config->get_aa2label();
	
	model_aa_labels.clear();

	for (aa=Ala; aa<=Val; aa++)
	{
		if (aa==Ile)
			continue;

		model_aa_labels.push_back(aa2label[aa]);
	}
	model_aa_labels.push_back("M+16");
	model_aa_labels.push_back("Q-17");
}



/********************************************************************
Creates a conversion table for converting session aas to the model aa
idxs. If there is a session aa that is not present it gets the org_aa
*********************************************************************/
void PeakRankModel::convert_session_aas_to_model_aas()
{
	const vector<string>& aa2label = config->get_aa2label();
	const vector<int>&    org_aas  = config->get_org_aa();
	const int max_session_aa_idx   = config->get_max_session_aa_idx();
	const vector<int>& session_aas = config->get_session_aas();

	session_aas_to_model_aas.clear();
	session_aas_to_model_aas.resize(max_session_aa_idx+1,-999);

	
	int i;
	for (i=0; i<session_aas.size(); i++)
	{
		const int aa = session_aas[i];
		const string label = aa2label[aa];

		// look for match
		int j;
		for (j=0; j<model_aa_labels.size(); j++)
			if (model_aa_labels[j] == label)
				break;

		if (j<model_aa_labels.size())
		{
			session_aas_to_model_aas[aa]=j;
			continue;
		}

		// try the org aa for this one
		const int org_aa = org_aas[aa];
		const string org_label = aa2label[org_aa];
		
		for (j=0; j<model_aa_labels.size(); j++)
			if (model_aa_labels[j] == org_label)
				break;

		if (j<model_aa_labels.size())
		{
			session_aas_to_model_aas[aa]=j;
			continue;
		}

		// check if we can use I/L
		if (org_label == "I")
		{
			const string l_label = "L";
			for (j=0; j<model_aa_labels.size(); j++)
			if (model_aa_labels[j] == l_label)
				break;

			if (j<model_aa_labels.size())
			{
				session_aas_to_model_aas[aa]=j;
				continue;
			}
		}

		// Error
		cout << "Error: couldn't find aa for " << label << " (" << aa << ")" << endl;
		exit(1);
	}


/*	cout << "++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "AA conversion table:" << endl;
	for (i=0; i<session_aas_to_model_aas.size(); i++)
	{
		cout << i << " " << aa2label[i] << " --> " << session_aas_to_model_aas[i];
		if (session_aas_to_model_aas[i]>=0)
			cout << "  " << model_aa_labels[session_aas_to_model_aas[i]];
		cout << endl;
	}
	cout << "-----------------------------------------" << endl;*/
}


/***********************************************************************
Changes the original sequence aas (with their int encoding) to the model's
own amino acid int codes. This enables the model to be applied to peptides
with PTMs that were not in the original training set and also takes care
of discrepancies that could result due to changes in amino acid encoding.
************************************************************************/
void PeakRankModel::convert_aas_to_model_aas(const vector<int>& org_aas, vector<int>& model_aas) const
{
	const int max_session_aa = config->get_max_session_aa_idx();
	int i;
	model_aas.resize(org_aas.size());
	for (i=0; i<org_aas.size(); i++)
	{
		if (org_aas[i]>max_session_aa)
		{
			cout << "Error: session AAs don't include: " << org_aas[i] << endl;
			exit(1);
		}
		model_aas[i]=session_aas_to_model_aas[org_aas[i]];
		if (model_aas[i]<0)
		{
			cout << "Error: couldn't find a rank model conversion for amino acid " << org_aas[i] << endl;
			exit(1);
		}
	}
}









void PeakRankModel::set_binary_feature_names()
{
	const int num_aas = this->get_num_model_aas();
	const string dis_labels[]={"-3","-2","-1","+1","+2","+3"};
	const int num_dis_labels = 6;

	int i;

//	cout << "Initialzing binary features for " << num_aas << " amino acids " << endl
//	for (i=0; i<num_aas; i++)
//		cout << i << "\t" << model_aa_labels[i] << endl;

	this->binary_feature_names.clear();
	
	for (i=0; i<num_dis_labels; i++)
	{
		const string dis_label = dis_labels[i];
		int aa;
		for (aa =0; aa<num_aas; aa++)
		{
			string fname = "IND CUT HAS " + model_aa_labels[aa] + " AT " + dis_label;
		//	binary_feature_names.push_back(fname);
		}
	}

/*	for (i=0; i<num_aas; i++)
	{
		int j;
		for (j=0; j<num_aas; j++)
		{
			string fname = "IND CUT IS BETWEEN " + model_aa_labels[i] + " - " + model_aa_labels[j];
			binary_feature_names.push_back(fname);
		}
	} */



	binary_feature_names.push_back("IND CUT AT N+1");
	binary_feature_names.push_back("IND CUT AT C-1");

/*	binary_feature_names.push_back("IND CUT AT N+1");
	binary_feature_names.push_back("IND CUT AT N+2");
	binary_feature_names.push_back("IND CUT AT N+3");
	binary_feature_names.push_back("IND CUT AT N+4");
	binary_feature_names.push_back("IND CUT AT N+5");

	binary_feature_names.push_back("IND CUT AT C-1");
	binary_feature_names.push_back("IND CUT AT C-2");
	binary_feature_names.push_back("IND CUT AT C-3");
	binary_feature_names.push_back("IND CUT AT C-4");
	binary_feature_names.push_back("IND CUT AT C-5");*/

	for (i=0; i<num_aas; i++)
	{
		string fname = "IND CUT AT +1, N AA IS " + model_aa_labels[i];
		binary_feature_names.push_back(fname);
	}

	for (i=0; i<num_aas; i++)
	{
		string fname = "IND CUT AT -1, C AA IS " + model_aa_labels[i];
		binary_feature_names.push_back(fname);
	}



/*	binary_feature_names.push_back("IND CUT LESS THAN MID PEPTIDE");
	binary_feature_names.push_back("IND CUT MORE THAN MID PEPTIDE");
	binary_feature_names.push_back("IND CUT OVER CHARGE LESS THAN MID OBS");
	binary_feature_names.push_back("IND CUT OVER CHARGE MORE THAN MID OBS");*/
}

void PeakRankModel::set_real_feature_names()
{
	const int num_aas = this->get_num_model_aas();
	int i;

	real_feature_names.clear();

	
	const string prefixes[]={"BEFORE MID, N SIDE RKH ","AFTER MID, N SIDE RKH ",
							 "BEFORE MID, C SIDE RKH ","AFTER MID, C SIDE RKH "};
	

	int c;
	for (c=0 ; c< num_RKH_combos; c++)
		real_feature_names.push_back("DIS MIN, N SIDE RHK " + combos[c]);
	for (c=0 ; c< num_RKH_combos; c++)
		real_feature_names.push_back("DIS MIN, C SIDE RHK " + combos[c]);
	for (c=0 ; c< num_RKH_combos; c++)
		real_feature_names.push_back("DIS MAX, N SIDE RHK " + combos[c]);
	for (c=0 ; c< num_RKH_combos; c++)
		real_feature_names.push_back("DIS MAX, C SIDE RHK " + combos[c]);

	int p;
	for (p=0; p<4; p++)
	{
		int c;
		for (c=0 ; c< num_RKH_combos; c++)
		{
			int aa;
			for (aa=0; aa<num_aas; aa++)
			{
				string fname = prefixes[p] + combos[c] + ", AA N of cut " + this->model_aa_labels[aa];
				real_feature_names.push_back(fname);
			}
		}

		for (c=0 ; c< num_RKH_combos; c++)
		{
			int aa;
			for (aa=0; aa<num_aas; aa++)
			{
				string fname = prefixes[p] + combos[c] + ", AA C of cut " + this->model_aa_labels[aa];
				real_feature_names.push_back(fname);
			}
		}
	}

/*	real_feature_names.push_back("CUT MASS PROPORTION (LESS THAN MID)");
	real_feature_names.push_back("CUT MASS PROPORTION (MORE THAN MID)");

	real_feature_names.push_back("CUT MASS OVER CHARGE PROPORTION (LESS THAN MID OBS)");
	real_feature_names.push_back("CUT MASS OVER CHARGE ABS DIFF (LESS THAN MID OBS)");

	real_feature_names.push_back("CUT MASS OVER CHARGE PROPORTION (MORE THAN MID OBS)");
	real_feature_names.push_back("CUT MASS OVER CHARGE ABS DIFF (MORE THAN MID OBS)"); */

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("# N SIDE " + model_aa_labels[i]);
	
	real_feature_names.push_back("# N SIDE BASIC AMINO ACIDS");

	for (i=0; i<num_aas; i++)
		real_feature_names.push_back("# C SIDE " + model_aa_labels[i]);

	real_feature_names.push_back("# C SIDE BASIC AMINO ACIDS");

	
}


struct basicity_pair {
	bool operator< (const basicity_pair& other) const
	{
		return basicity_score<other.basicity_score;
	}

	float basicity_score;
	int   idx;
};

/****************************************************************************
// reoprts sizes of partitioned training sets
// partitioned according to charge / size_idx / mobility
// if path is given, creates the training files
*****************************************************************************/
void PeakRankModel::partition_training_samples(
						const vector< vector<TrainingPeptide> >& all_tps,
						char *file_path_prefix,
						char *test_path_prefix,
						int   minimal_size,
						float prop_ts) const
{
	if (file_path_prefix)
		cout << "TRAIN: " << file_path_prefix << endl;
	if (test_path_prefix)
		cout << "TEST : " << test_path_prefix << endl;

	cout << endl << endl;
	cout << "Training partiotioning: " << endl;
	cout << "----------------------- " << endl << endl;

	int charge;
	for (charge=1; charge<all_tps.size(); charge++)
	{
		if (all_tps[charge].size()==0)
			continue;

		cout << endl << "**************************************" << endl;
		cout << endl << "CHARGE " << charge << " total " << all_tps[charge].size() << " peptides..." << endl;

		const int num_size_idxs = this->size_thresholds[charge].size()+1;
		vector< vector<int> > pep_counts;
		vector< vector< vector<int> > > pep_idxs;

		pep_counts.resize(num_size_idxs);
		pep_idxs.resize(num_size_idxs);
		int size_idx;
		for (size_idx=0; size_idx<num_size_idxs; size_idx++)
		{
			pep_counts[size_idx].resize(4,0);
			pep_idxs[size_idx].resize(4);
		}


		int i;
		for (i=0; i<all_tps[charge].size(); i++)
		{
			const TrainingPeptide& tp = all_tps[charge][i];
			const int size_idx = this->get_size_group(charge,tp.pm_with_19);
			pep_counts[size_idx][0]++;
			pep_counts[size_idx][tp.mobility]++;
			pep_idxs[size_idx][tp.mobility].push_back(i);
		}

		for (size_idx=0; size_idx<num_size_idxs; size_idx++)
		{
			cout << size_idx << "\t";
			if (pep_counts[size_idx][0]>0)
			{
				cout << pep_counts[size_idx][0] << "\t" << fixed << setprecision(4) << 
					pep_counts[size_idx][0] / (double)all_tps[charge].size() << "\t";
			}
			else
				cout << "  -\t  -\t";

			int m;
			for (m=MOBILE; m<=NONMOBILE; m++)
				if (pep_counts[size_idx][m]>0)
				{
					cout << pep_counts[size_idx][m];
					if (pep_counts[size_idx][m]<minimal_size)
						cout << " *";
					cout << "\t" << pep_counts[size_idx][m]/(double)pep_counts[size_idx][0] << "\t";
				}
				else
					cout << "  -\t 0\t";

			cout << endl;
		}
		cout << endl;

		if (! file_path_prefix)
			continue;

		
		for (size_idx=0; size_idx<num_size_idxs; size_idx++)
		{
			int m;
			for (m=MOBILE; m<=NONMOBILE; m++)
			{
				vector<int> idxs_to_write = pep_idxs[size_idx][m];

				cout << endl;
				cout << "set " << charge << " " << size_idx << " " << m << endl;
			
				
				// try adding samples from adjacent size_idxs
				if (0 && idxs_to_write.size()<minimal_size && size_idx>0)
				{
					const vector<int>& prev_size_idxs = pep_idxs[size_idx-1][m];
					int i;

					// sort samples according to size
					vector<basicity_pair> pairs;
					pairs.clear();
					
					for (i=0; i<prev_size_idxs.size(); i++)
					{
						basicity_pair bp;
						bp.idx = prev_size_idxs[i];
						bp.basicity_score = all_tps[charge][prev_size_idxs[i]].pm_with_19;
						pairs.push_back(bp);
					}
					sort(pairs.begin(),pairs.end());

					for (i=pairs.size()-1; i>=0; i--)
					{
						if (idxs_to_write.size() == minimal_size)
							break;
						idxs_to_write.push_back(pairs[i].idx);
					//	cout << setprecision(2) << all_tps[charge][pairs[i].idx].pm_with_19 << " " <<
					//		all_tps[charge][pairs[i].idx].mobility << " ";
						
					}
					cout << endl;
					cout << " .. took " << pairs.size()-i << " samples from size " << size_idx-1 << endl;
				}

				// try adding samples from adjacent size_idxs
				if (0 && idxs_to_write.size()<minimal_size && size_idx<num_size_idxs-1)
				{
					const vector<int>& next_size_idxs = pep_idxs[size_idx+1][m];
					int i;
						// sort samples according to size
					vector<basicity_pair> pairs;
					pairs.clear();
					
					for (i=0; i<next_size_idxs.size(); i++)
					{
						basicity_pair bp;
						bp.idx = next_size_idxs[i];
						bp.basicity_score = all_tps[charge][next_size_idxs[i]].pm_with_19;
						pairs.push_back(bp);
					}
					sort(pairs.begin(),pairs.end());

					for (i=0; i<pairs.size(); i++)
					{
						if (idxs_to_write.size() == minimal_size)
							break;
						idxs_to_write.push_back(pairs[i].idx);
					//	cout << setprecision(2) << all_tps[charge][pairs[i].idx].pm_with_19 << " " <<
					//		all_tps[charge][pairs[i].idx].mobility << " ";
					}
					cout << endl;
					cout << " .. took " << i << " samples from size " << size_idx+1 << endl;
				}

				// try adding samples from another mobility state
				// choose most appropriate samepls (in terms of basicity)
				if (idxs_to_write.size()<minimal_size)
				{

					int mobility_to_take_from=0;

					if (m==MOBILE || m == NONMOBILE)
					{
						mobility_to_take_from=PARTIALLYMOBILE;
					}
					else if (pep_idxs[size_idx][MOBILE].size()>
							 pep_idxs[size_idx][NONMOBILE].size())
					{
						mobility_to_take_from = MOBILE;
					}
					else
						mobility_to_take_from = NONMOBILE;

					if (idxs_to_write.size() + pep_counts[size_idx][mobility_to_take_from] < minimal_size)
					{
						cout << "Error: insufficient number of training samples for charge " <<
							charge << " size idx " << size_idx << endl;
						continue;
					} 

					vector<basicity_pair> pairs;
					int i;
					pairs.clear();
					
					for (i=0; i<pep_idxs[size_idx][mobility_to_take_from].size(); i++)
					{
						basicity_pair bp;
						bp.idx = pep_idxs[size_idx][mobility_to_take_from][i];
						bp.basicity_score = all_tps[charge][bp.idx].get_basicity_score();
						pairs.push_back(bp);
					}

					sort(pairs.begin(),pairs.end());

					if (mobility_to_take_from>m)
					{
						// take least basic peptides
						int i;
						for (i=0; i<pairs.size() && idxs_to_write.size()<minimal_size; i++)
							idxs_to_write.push_back(pairs[i].idx);

						cout << " .. took " << i << " samples mobility " << mobility_to_take_from << endl;
					}
					else
					{
						// take most basic peptides
						int i;
						for (i=pairs.size()-1; i>=0 && idxs_to_write.size()<minimal_size; i--)
						{
							idxs_to_write.push_back(pairs[i].idx);
						//	cout << pairs[i].basicity_score << " ";
						}
						cout << endl;

						cout << " .. took " << pairs.size()-i << " samples mobility " << mobility_to_take_from << endl;
					}
				}

				sort(idxs_to_write.begin(),idxs_to_write.end());

				if (idxs_to_write.size()<minimal_size)
				{
					cout << "Error: could not add enough samples to create set for charge " << charge << " size " << size_idx << " mobility " << m << endl;
					exit(1);
				}
				
				
				if (! test_path_prefix)
				{
					char fname[256];
					sprintf(fname,"%s_%d_%d_%d.txt",file_path_prefix,charge,size_idx,m);

					ofstream ofs(fname);
					ofs << idxs_to_write.size() << endl;
					for (i=0; i<idxs_to_write.size(); i++)
						all_tps[charge][idxs_to_write[i]].write_to_stream(ofs);
					ofs.close();

					cout << "Wrote " << idxs_to_write.size() << " to " << fname << endl;
				}
				else
				{
					char fname[256], ts_name[256];
					const int num_test = int(prop_ts * (float)idxs_to_write.size());
					vector<size_t>  ts_idxs;
					vector<bool> ts_ind;
					chooseKFromN(num_test,idxs_to_write.size(),ts_idxs);
					ts_ind.resize(idxs_to_write.size(),false);
					int i;
					for (i=0; i<ts_idxs.size(); i++)
						ts_ind[ts_idxs[i]]=true;

					sprintf(fname,"%s_%d_%d_%d.txt",file_path_prefix,charge,size_idx,m);
					ofstream ofs(fname);
					ofs << idxs_to_write.size() - num_test<< endl;
					for (i=0; i<idxs_to_write.size(); i++)
						if (! ts_ind[i])
							all_tps[charge][idxs_to_write[i]].write_to_stream(ofs);
					ofs.close();

					cout << "Wrote " << idxs_to_write.size() - num_test << " to " << fname << endl;

					sprintf(ts_name,"%s_%d_%d_%d.txt",test_path_prefix,charge,size_idx,m);
					ofstream ofs_test(ts_name);
					ofs_test << num_test<< endl;
					for (i=0; i<idxs_to_write.size(); i++)
						if (ts_ind[i])
							all_tps[charge][idxs_to_write[i]].write_to_stream(ofs_test);
					ofs_test.close();


					cout << "Wrote " << num_test << " to " << ts_name << endl;

				}

				// find mix/max sizes and min/max basicity scores
				mass_t min_pm_with_19 = 1000000.0;
				mass_t max_pm_with_19 = 0.0;
				float min_basicity_score = 10000.0;
				float max_basicity_score = -1.0;

				for (i=0; i<idxs_to_write.size(); i++)
				{
					const TrainingPeptide& tp = all_tps[charge][idxs_to_write[i]];
					const float bs = tp.get_basicity_score();

					if (tp.pm_with_19>max_pm_with_19)
						max_pm_with_19 = tp.pm_with_19;

					if (tp.pm_with_19<min_pm_with_19)
						min_pm_with_19 = tp.pm_with_19;

					if (bs<min_basicity_score)
						min_basicity_score = bs;

					if (bs>max_basicity_score)
						max_basicity_score = bs;
				}
				cout << "pm ranges: " << setprecision(2) << min_pm_with_19 << " - " << max_pm_with_19 << endl;
				cout << "basicity : " << min_basicity_score << " - " << max_basicity_score << endl;

				cout << "-------------------------------------" << endl;
			}
		}
	}
}


/***************************************************************************
Reads the sample tps into the dataset and creats sample for a given 
fragment type idx
****************************************************************************/
void PeakRankModel::read_training_peptides_into_rank_boost_dataset(
										int frag_type_idx,
										int spec_charge,
										const vector<TrainingPeptide>& sample_tps,
										RankBoostDataset& rank_ds,
										vector<float>& peak_intens,
										vector<PeakStart>& peak_starts,
										vector<float>& max_annotated_intens) const
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	rank_ds.clear();
	peak_intens.clear();
	peak_starts.clear();
	max_annotated_intens.clear();

	int i;
	for (i=0; i<sample_tps.size(); i++)
	{
		const TrainingPeptide& tp = sample_tps[i];
		const int frag_pos = tp.get_frag_idx_pos(frag_type_idx);

		if (frag_pos<0)
			continue;

		PeakStart ps;

		ps.peptide_sample_idx = i;
		ps.peak_start_idx=rank_ds.get_num_samples();
		ps.num_peaks=0;

		// scan all cuts for max inten
		float max_ann_inten=0;
		int f;
		for (f=0; f<tp.intens.size(); f++)
		{
			int c;
			for (c=1; c<tp.intens[f].size()-1; c++)
				if (tp.intens[f][c]>max_ann_inten)
					max_ann_inten=tp.intens[f][c];
		}

		if (max_ann_inten<=0)
			continue;
		
		const vector<float>& intens = tp.intens[frag_pos];
		vector<int> cut_ranks;

		// gives integer values to cut idxs according to intensity
		calc_cut_ranks(intens,cut_ranks);

		mass_t c_term_mass = tp.n_mass;
		int cut_idx;
		for (cut_idx=1; cut_idx<intens.size(); cut_idx++)
			c_term_mass+=aa2mass[tp.amino_acids[cut_idx-1]];

		mass_t cut_mass=tp.n_mass;
		for (cut_idx=1; cut_idx<intens.size(); cut_idx++)
		{
			cut_mass+=aa2mass[tp.amino_acids[cut_idx-1]];
			if (intens[cut_idx]>=0)
			{
				RankBoostSample rbs;
				const FragmentType& frag = config->get_fragment(frag_type_idx);

				if (feature_set_type == 0)
				{
					fill_simple_peak_features(tp.amino_acids, cut_idx, 
						cut_mass, tp.pm_with_19, spec_charge, frag, rbs);
				}
				else if (feature_set_type == 1)
				{
					fill_advanced_peak_features(tp.amino_acids, cut_idx, 
						cut_mass, tp.pm_with_19, spec_charge, frag, rbs);
				}
				else if (feature_set_type == 2)
				{
					fill_partial_denovo_peak_features(tp.n_mass, c_term_mass, tp.amino_acids, cut_idx, 
						cut_mass, tp.pm_with_19, spec_charge, frag, tp.best_n_removed, tp.best_c_removed, rbs);
				}
				else
				{
					cout << "Error: feature set type not supported: " << feature_set_type << endl;
					exit(1);
				}

				rbs.groupIndex = i;
				rbs.rank_in_group = cut_ranks[cut_idx];

				rbs.tag1 = rank_ds.get_num_samples(); // sample idx
				rbs.tag2 = cut_idx; // cut idx
				rbs.tag3 = tp.amino_acids.size(); // peptide length, this tag is used to filter
												  // the error estimates for a specific length

				rank_ds.add_sample(rbs);
				peak_intens.push_back(intens[cut_idx]);
			}	
		}

		ps.num_peaks = rank_ds.get_num_samples() - ps.peak_start_idx;
		peak_starts.push_back(ps);
		max_annotated_intens.push_back(max_ann_inten);
	}

	rank_ds.set_num_groups(sample_tps.size());
}



/***************************************************************************
Initializes the models according to specified max charge defaults (if not 
already initialized). Assumes all training files have the same prefix
and differ on the charge/size/mobility idxs).
****************************************************************************/
void PeakRankModel::train_all_partition_models(
								int frag_fill_type,
								char *prefix_path,
								int	  sel_charge,
								int   sel_size_idx,
								int	  sel_mobility,
								int	  frag_type_idx,
								char *report_dir,
								int   num_rounds,
								char *test_set,
								int	  test_peptide_length,
								char *stop_signal_file,
								weight_t max_weight_ratio)
{
	const int max_charge = size_thresholds.size() -1;

	this->feature_set_type = frag_fill_type;

	if (partition_models.size()<max_charge+1)
		partition_models.resize(max_charge+1);

	int charge; 
	for (charge=1; charge<=max_charge; charge++)
	{
		if (sel_charge>0 && sel_charge != charge)
			continue;

		const int num_size_idxs = size_thresholds[charge].size()+1;
		if (partition_models[charge].size()<num_size_idxs)
			partition_models[charge].resize(num_size_idxs);

		cout << "FRAG FILL TYPE: " << frag_fill_type << endl;
		cout << "CHARGE " << charge << "  , # size idxs: " << num_size_idxs << endl;
		cout << "Test peptide length: " << test_peptide_length << endl;
		int size_idx;
		for (size_idx=0; size_idx<num_size_idxs; size_idx++)
		{
			if (partition_models[charge][size_idx].size()<=NONMOBILE)
				partition_models[charge][size_idx].resize(NONMOBILE+1,NULL);

			int mobility;
			for (mobility=MOBILE; mobility<=NONMOBILE; mobility++)
			{
				if ((sel_charge>=0 && sel_charge != charge) ||
					(sel_size_idx>=0 && sel_size_idx != size_idx) ||
					(sel_mobility>=0 && sel_mobility != mobility) )
					continue;

				char tr_file[256];
				sprintf(tr_file,"%s_%d_%d_%d.txt",prefix_path,charge,size_idx,mobility);

				if (! partition_models[charge][size_idx][mobility])
					partition_models[charge][size_idx][mobility] = new PartitionModel;

				cout << "Training models for charge " << charge << " size " << size_idx<< " mobility " <<
					mobility << endl;
				cout << "Max weight ratio " << max_weight_ratio << endl;

				partition_models[charge][size_idx][mobility]->set_partition_name(get_peak_rank_model_name(),
					charge, size_idx, mobility);
			
				partition_models[charge][size_idx][mobility]->train_partition_model(this,
					tr_file,charge,size_idx,mobility,frag_type_idx,report_dir,
					num_rounds,test_set, test_peptide_length, stop_signal_file, 
					max_weight_ratio);

				cout << endl;
			}
		}
	}
}


/**********************************************************************
Lists the idxs of the partition models
charge size_idx mobility (not fragment)
***********************************************************************/
void PeakRankModel::list_all_model_idxs()
{
	int charge;
	for (charge=1; charge<size_thresholds.size(); charge++)
	{
		int size_idx;
		for (size_idx = 0; size_idx<=size_thresholds[charge].size(); size_idx++)
		{
			int mobility;
			for (mobility=MOBILE; mobility<=NONMOBILE; mobility++)
			{
				cout << charge << " " << size_idx << " " << mobility << endl;
			}
		}
	}
}


void PeakRankModel::print_model_init_stats() const
{
	if (this->feature_set_type<=2)
	{
		const int num_tabs = 2+(peak_rank_model_name.length()+1)/8;
		string gap_str ="";
		int i;
		for (i=0; i<num_tabs; i++)
			gap_str += "\t";

		cout << gap_str;
		
		for (i=0; i<10; i++)
			cout << " " << i;
		cout << endl;

		int charge;
		for (charge=1; charge<size_thresholds.size(); charge++)
		{
			int size_idx;
			for (size_idx = 0; size_idx<=size_thresholds[charge].size(); size_idx++)
			{
				int mobility;
				for (mobility=MOBILE; mobility<=NONMOBILE; mobility++)
				{
					if (this->partition_models[charge][size_idx][mobility])
						partition_models[charge][size_idx][mobility]->print_model_stats();
				}
			}
		}
		cout << endl;
	}
}


// set intens so maximal intensity of the fragments is 1.0
void find_ranks(const vector<intensity_t>& intens, vector<int>& ranks)
{
	vector<inten_pair> pairs;

	pairs.resize(intens.size());

	int num_with_inten=0;
	int i;
	for (i=0; i<intens.size(); i++)
	{
		pairs[i].idx=i;
		pairs[i].inten = intens[i];
		if (intens[i]>=0)
			num_with_inten++;
	}

	sort(pairs.begin(),pairs.end());

	const int missing_peak_rank = (intens.size()+num_with_inten)/2;

	ranks.resize(intens.size());
	for (i=0; i<pairs.size(); i++)
	{
		if (pairs[i].inten <0)
		{
			ranks[pairs[i].idx]=-1;
		}
		else if (pairs[i].inten == 0)
		{
			ranks[pairs[i].idx]=missing_peak_rank;
		}
		else
			ranks[pairs[i].idx]=i+1;
	}
}


void normalize_intens(vector<intensity_t>& intens)
{
	if (intens.size()<=1)
		return;

	float max_inten = intens[1];
	int i;
	for (i=2; i<intens.size(); i++)
		if (intens[i]>max_inten)
			max_inten = intens[i];

	if (max_inten<=0.0)
		return;

	float one_over = 1.0 / max_inten;

	for (i=1; i<intens.size(); i++)
		intens[i] *= one_over;
}


void PeakRankModel::set_mass_detection_defaults()
{
	max_detected_mass = 2000.0;

	charge_min_mass_coefficients.resize(4,0);
	charge_min_mass_coefficients[1]=0.24;
	charge_min_mass_coefficients[2]=0.12;
	charge_min_mass_coefficients[3]=0.08;
}

void PeakRankModel::set_size_thresholds()
{
	size_thresholds.resize(4);

	size_thresholds[1].push_back(1150.0);
	size_thresholds[1].push_back(1400.0);

	size_thresholds[2].push_back(1100.0);
	size_thresholds[2].push_back(1300.0);
	size_thresholds[2].push_back(1600.0);
	size_thresholds[2].push_back(1900.0);
	size_thresholds[2].push_back(2400.0);

	size_thresholds[3].push_back(1950.0);
	size_thresholds[3].push_back(2450.0);
	size_thresholds[3].push_back(3000.0);
}



mass_t PeakRankModel::calc_min_detected_mass(mass_t pm_with_19, int charge) const
{
	if (charge>= charge_min_mass_coefficients.size())
		charge = charge_min_mass_coefficients.size()-1;
	return (charge_min_mass_coefficients[charge]*pm_with_19);
}



bool PeakRankModel::read_peak_rank_model(Config *_config, const char *name, bool silent_ind,
				int specific_charge, int specific_size, int specific_mobility)
{
	config = _config;
	peak_rank_model_name = string(name);

	const string model_name_prefix = config->get_resource_dir() + "/" + string(name);
	const string main_model_file = model_name_prefix + "_rank_model.txt";


	ifstream main_stream(main_model_file.c_str());

	if (! main_stream.is_open() || ! main_stream.good())
	{
		cout << "Error: couldn't open file for reading: " << main_model_file << endl;
		exit(1);
	}

	string m_name;
	int num_aa_labels=0;

	feature_set_type = -1;

	main_stream >> m_name;
	main_stream >> feature_set_type;
	main_stream >> num_aa_labels;

	if (feature_set_type<0 || feature_set_type>5 || num_aa_labels<19)
	{
		cout << "Error: bad input parameters in PeakRankModel!" << endl;
		exit(1);
	}


	model_aa_labels.resize(num_aa_labels);
	int i;
	for (i=0; i<num_aa_labels; i++)
		main_stream >> model_aa_labels[i];

	convert_session_aas_to_model_aas();

	max_detected_mass=-1;
	main_stream >> this->max_detected_mass;

	int num_min_ratios=-1;
	main_stream >> num_min_ratios;
	this->charge_min_mass_coefficients.resize(num_min_ratios,1);

	for (i=0; i<num_min_ratios; i++)
		main_stream >> charge_min_mass_coefficients[i];

	int num_charges=0;
	main_stream >> num_charges;

	size_thresholds.resize(num_charges);
	partition_models.resize(num_charges);
	for (i=0; i<num_charges; i++)
	{
		int charge=-1;
		int num_sizes=0;

		main_stream >> charge;
		main_stream >> num_sizes;
		
		if (charge != i)
		{
			cout << "Error: charge mismatch is size thresholds!" << endl;
			exit(1);
		}

		size_thresholds[i].resize(num_sizes);
		int j;
		for (j=0; j<num_sizes; j++)
			main_stream >> size_thresholds[i][j];
		
		if (num_sizes>0)
		{
			partition_models[i].resize(num_sizes+1);
			for (j=0; j<partition_models[i].size(); j++)
				partition_models[i][j].resize(4);
		}
	}

	if (! silent_ind)
		cout << "Peak rank feature set type: " << feature_set_type << endl;


	// read parition models
	int num_models_read=0;
	int charge;
	for (charge=1; charge<partition_models.size(); charge++)
	{
		if (specific_charge>=0 && charge != specific_charge)
			continue;

		int size_idx;
		for (size_idx=0; size_idx<partition_models[charge].size(); size_idx++)
		{

			if (specific_size>=0 && size_idx != specific_size)
				continue;

			int mobility;
			for (mobility=MOBILE; mobility<=NONMOBILE; mobility++)
			{

				if (specific_mobility>=0 && mobility != specific_mobility)
					continue;

				partition_models[charge][size_idx][mobility] = new PartitionModel;
				partition_models[charge][size_idx][mobility]->set_feature_set_type(feature_set_type);
	
				num_models_read+=partition_models[charge][size_idx][mobility]->read_partition_model(model_name_prefix,
					config,charge,size_idx,mobility);
				partition_models[charge][size_idx][mobility]->set_partition_name(name,charge,size_idx,mobility);
			}
		}
	}

	if (! silent_ind)
		cout << "Read " << num_models_read << " fragment rank models..." << endl;
	return true;
}


void PeakRankModel::write_peak_rank_model(char *name, char *out_dir)
{
	const string model_name_prefix = (out_dir? string(out_dir) : config->get_resource_dir())
									  + "/" + string(name);
	const string main_model_file = model_name_prefix + "_rank_model.txt";
	

	ofstream main_stream(main_model_file.c_str());

	if (! main_stream.good())
	{
		cout << "Error: couldn't open file for writing: " << main_model_file << endl;
		exit(1);
	}

	// name
	main_stream << name << " " << this->feature_set_type << endl;

	// aa labels
	main_stream << model_aa_labels.size();
	int a;
	for (a=0; a<model_aa_labels.size(); a++)
		main_stream << " " << model_aa_labels[a];
	main_stream << endl;

	// min max detected masses
	int charge;
	main_stream << setprecision(5) << max_detected_mass << endl;
	main_stream << charge_min_mass_coefficients.size();
	for (charge=0; charge<charge_min_mass_coefficients.size(); charge++)
		main_stream << " " << charge_min_mass_coefficients[charge];
	main_stream << endl;

	// size thresholds
	main_stream << size_thresholds.size() << endl;
	for (charge=0; charge<size_thresholds.size(); charge++)
	{
		main_stream << charge << " " << size_thresholds[charge].size();
		int size_idx;
		for (size_idx=0; size_idx<size_thresholds[charge].size(); size_idx++)
			main_stream << " " << setprecision(5) << size_thresholds[charge][size_idx];
		main_stream << endl;
	}

	// write the parition models
	int models_written = 0;
	for (charge=0; charge<partition_models.size(); charge++)
	{
		int size_idx;
		for (size_idx=0; size_idx<partition_models[charge].size(); size_idx++)
		{
			int mobility;
			for (mobility=MOBILE; mobility<=NONMOBILE; mobility++)
			{
				if (partition_models[charge][size_idx][mobility])
				{
					int n=partition_models[charge][size_idx][mobility]->write_partition_model(model_name_prefix);
					models_written += n;
				}
			}
		}
	}

	cout << "Wrote rank model: " << name << endl;
	cout << "A total of " << models_written << " fragment models were written..." << endl;
}


void TrainingPeptide::create_training_peptide(const PeakRankModel& rm, 
											  const AnnotatedSpectrum& as)
{
	const Config *config = as.getConfig();
	amino_acids = as.getPeptide().get_amino_acids();
	length = amino_acids.size();
	charge = as.getCharge();
	mobility = get_proton_mobility(as.getPeptide(),charge);
	pm_with_19 = as.get_true_mass_with_19();

	const mass_t min_detected_mass = rm.calc_min_detected_mass(as.get_true_mass_with_19(),charge);
	const mass_t max_detected_mass = rm.get_max_detected_mass();
	const vector<Breakage>& breakages = as.get_breakages();

	const bool verbose=false;


	// find all participating fragment types
	frag_idxs.clear();
	int b;
	for (b=1; b<breakages.size()-1; b++)
	{
		int f;
		for (f=0; f<breakages[b].fragments.size(); f++)
		{
			const int frag_idx = breakages[b].fragments[f].frag_type_idx;
			if (breakages[b].fragments[f].intensity>0.01 &&
				this->get_frag_idx_pos(frag_idx)<0)
				frag_idxs.push_back(frag_idx);
		}
	}

	intens.resize(frag_idxs.size());

//	cout << "Frags: " << ranks.size() << endl;

	if (verbose)
		cout << as.getPeptide().as_string(config) << endl;
	
	int f;
	for (f=0; f<frag_idxs.size(); f++)
	{
		const int frag_idx = frag_idxs[f];
		const FragmentType& frag = config->get_fragment(frag_idx);
		vector<mass_t>& frag_intens = intens[f];

		frag_intens.resize(length,0);
		int num_inten_peaks=0;

		int i;
		for (i=0; i<length; i++)
		{
			const int pos = breakages[i].get_position_of_frag_idx(frag_idx);
			if (pos>=0)
			{
				frag_intens[i]=breakages[i].fragments[pos].intensity;
				num_inten_peaks++;
			}
		}


		frag_intens[0]=num_inten_peaks;

		for (i=1; i<frag_intens.size(); i++)
		{
			const mass_t exp_mass = frag.calc_expected_mass(breakages[i].mass,pm_with_19);
			if (exp_mass<min_detected_mass || exp_mass>max_detected_mass)
			{
				frag_intens[i]=-999.0;
			}
	
		}

		if (verbose)
		{
			cout << frag.label << "\t";
			for (i=0; i<intens[f].size(); i++)
				cout << " " << fixed << setprecision(3) << intens[f][i];
			cout << endl;
		}
	}

	if (verbose)
		cout << endl;
}

void TrainingPeptide::get_ranks_for_frag_idx(int frag_idx, vector<int>& ranks) const
{
	const int pos = get_frag_idx_pos(frag_idx);
	if (pos<0)
	{
		cout << "Error: ranks not collect for frag " << frag_idx << endl;
		exit(1);
	}
	

	find_ranks(intens[pos],ranks);
}

// made up score, higher means more basic amino acids
float TrainingPeptide::get_basicity_score() const
{
	float basicity_score=0;
	int i;

	for (i=0; i<this->amino_acids.size(); i++)
	{
		if (amino_acids[i] == Arg)
			basicity_score+=1.0;
		if (amino_acids[i] == Lys)
			basicity_score+=0.55;
		if (amino_acids[i] == His)
			basicity_score+=0.3;
	}
	return basicity_score;
}


void TrainingPeptide::write_to_stream(ofstream& ofs) const
{
	if (frag_idxs.size() != intens.size())
	{
		cout << "Error: frag_idxs not same size ss intens!" << endl;
		exit(1);
	}

	// NP3 GOT precision from 3 to 5
	ofs << this->charge << " " << this->mobility << " " << fixed << setprecision(5) << this->pm_with_19 << endl;
	ofs << this->amino_acids.size();
	int i;
	for (i=0; i<amino_acids.size(); i++)
		ofs << " " << amino_acids[i];
	ofs << endl;
	ofs << this->frag_idxs.size();
	for (i=0; i<frag_idxs.size(); i++)
		ofs << " " << frag_idxs[i];
	ofs << endl;
	for (i=0; i<intens.size(); i++)
	{
		ofs << intens[i].size();
		int j;
		for (j=0; j<intens[i].size(); j++)
			ofs << " " << intens[i][j];
			
		ofs << endl;
	}
}

bool TrainingPeptide::read_from_stream(ifstream& ifs)
{
	char buff[2048];
	ifs.getline(buff,512);

	istringstream is(buff);

	is >> charge >> mobility >> pm_with_19;
	if (charge<=0 || charge >100 || mobility<MOBILE || mobility>NONMOBILE ||
		pm_with_19<50 || pm_with_19>1000000	)
	{
		return false;
	}
	
	length=0;
	ifs.getline(buff,1024);
	is.clear();
	is.str(buff);
	is >> length;
	if (length<=0 || length>1000)
	{
		return false;
	}
	int i;
	this->amino_acids.resize(length,-1);
	for (i=0; i<length; i++)
	{
		is >> amino_acids[i];
		if (amino_acids[i]<0 || amino_acids[i]>1000)
		{
			return false;
		}
	}

	int num_frags=0;
	ifs.getline(buff,2048);
	is.clear();
	is.str(buff);
	is >> num_frags;
	if (num_frags<1 || num_frags>100)
	{
		return false;
	}
	
	this->frag_idxs.resize(num_frags,-1);
	for (i=0; i<num_frags; i++)
	{
		is >> frag_idxs[i];
		if (frag_idxs[i]<0 || frag_idxs[i]>1000)
		{
			return false;
		}
	}
	intens.resize(num_frags);
	
	for (i=0; i<num_frags; i++)
	{
		ifs.getline(buff,2048);
		istringstream is(buff);

		int size=0;
		is >> size;
		intens[i].resize(size,0);
		int j;
		for (j=0; j<size; j++)
			is >> intens[i][j];
	}
	return true;
}

void  TrainingPeptide::print(Config *config, ostream& ofs) const
{
	const vector<string>& aa2label = config->get_aa2label();
	const vector<mass_t>& aa2mass  = config->get_aa2mass();
	int i;
	mass_t m=n_mass;

	ofs << this->best_n_removed << " " << this->best_c_removed << " " << this->n_mass << " (" << pm_with_19 << ") ";
	for (i=0; i<amino_acids.size(); i++)
	{
		ofs << aa2label[amino_acids[i]];
		m+=aa2mass[amino_acids[i]];
	}
	ofs << " " << m << endl;
}

/*
void read_data_into_training_peptides(const FileManager& fm, 
									  Config *config, 
									  PeakRankModel& rm,
									  vector<TrainingPeptide>& tps)
{
	FileSet fs;
	fs.select_all_files(fm);

	while (1)
	{
		SingleSpectrumFile *ssf;
		AnnotatedSpectrum as;
		if (! fs.get_next_spectrum(fm,config,&as,&ssf))
			break;

		as.annotate_spectrum(as.get_true_mass_with_19());

		TrainingPeptide tp;
		tp.create_training_peptide(rm,as);
		if (tp.frag_idxs.size()==0)
			continue;
		tps.push_back(tp);

	}

}
*/

void convert_list_to_trianing_peptide_file(char *list, char *tp_file,
										   char *model_name, char *ptm_line)
{
/*	PeakRankModel rm;
	AllScoreModels model;
	FileManager fm;
	vector<TrainingPeptide> all_tps;

	if (! model_name)
	{
		model.read_model("LTQ_LOW_TRYP");
	}
	else
		model.read_model(model_name);

	Config *config= model.get_config();

	if (! ptm_line)
	{
		char fixed_ptm_line[] = {"C+57:M+16:Q-17"};
		config->apply_selected_PTMs(fixed_ptm_line);
	}
	else
		config->apply_selected_PTMs(ptm_line);


	rm.set_mass_detection_defaults();
	fm.init_from_list_file(config,list);

	//read_data_into_training_peptides(fm,config,rm,all_tps);
	write_training_peptides_to_file(tp_file,all_tps);*/
}


void read_training_peptides_from_file(char *file, vector<TrainingPeptide>& all_tps,
									  int num_tp)
{

	ifstream ifs(file);
	char buff[64];

	if (! ifs.is_open())
	{
		cout << "Error: couldn't open: " << file << endl;
		exit(1);
	}

	int total_num_tp;
	ifs.getline(buff,64);
	sscanf(buff,"%d",&total_num_tp);

	if (num_tp<0 || num_tp>total_num_tp)
		num_tp = total_num_tp;

	
	const int start_idx=all_tps.size();
	all_tps.resize(start_idx+num_tp);

	int num_errors=0;
	int i;
	for (i=0; i<num_tp; i++)
		if (! all_tps[start_idx+i].read_from_stream(ifs))
			num_errors++;
	
	if (num_errors>0)
	{
		cout << "Warning: encountered " << num_errors << " errors reading " << num_tp << " training peptides..." << endl;
	}
	ifs.close();
}


void write_training_peptides_to_file(char *file, const vector<TrainingPeptide>& all_tps)
{
	ofstream ofs(file);
	ofs << all_tps.size() << endl;
	int i;
	for (i=0; i<all_tps.size(); i++)
	{
		all_tps[i].write_to_stream(ofs);
	}
	ofs.close();
}



void select_training_peptides(const vector<TrainingPeptide>& all_tps, 
							  vector<int>& selected_idxs,
							  int charge, 
							  int mobility, 
							  int min_length, 
							  int max_length, 
							  mass_t min_pm_with_19, 
							  mass_t max_pm_with_19)
{
	selected_idxs.clear();

	int i;
	for (i=0; i<all_tps.size(); i++)
	{
		const TrainingPeptide& tp = all_tps[i];
		if (charge>0 && charge != tp.charge)
			continue;
		
		if (mobility>0 && mobility != tp.mobility)
			continue;

		if (tp.length<min_length || tp.length>max_length ||
			tp.pm_with_19<min_pm_with_19 || tp.pm_with_19>max_pm_with_19)
			continue;

		selected_idxs.push_back(i);
	}
}




// for tps of a given charge
void size_mobility_stats(const vector<TrainingPeptide>& all_tps)
{
	vector< vector<int> > counts; // size / mobiliy (0 = all);
	int i;
	
	counts.resize(100);
	for (i=0; i<counts.size(); i++)
		counts[i].resize(4,0);

	for (i=0; i<all_tps.size(); i++)
	{
		int size_idx = (int)(all_tps[i].pm_with_19 * 0.01);

		if (all_tps[i].mobility<MOBILE || all_tps[i].mobility>NONMOBILE)
		{
			cout << "Error: bad mobility value!" << endl;
			exit(1);
		}

		if (size_idx<3 || size_idx>97)
		{
			cout << "Error: bad size_idx!" << endl;
			exit(1);
		}

		counts[size_idx][0]++;
		counts[size_idx][all_tps[i].mobility]++;
	}

	for (i=0; i<99; i++)
	{
		if (counts[i][0]==0)
			continue;

		cout << i*100;
		int j;
		for (j=0; j<4; j++)
		{
			cout << "\t" << counts[i][j];
			counts[99][j]+=counts[i][j];
		}
		cout << endl;	
	}
	
	cout << "---------------------------------------" << endl;
	
	for (i=0; i<4; i++)
		cout << "\t" << counts[99][i];
	cout << endl << endl;
}

// for tps of a given charge
void aa_composition_stats(const vector<TrainingPeptide>& all_tps, Config *config)
{
	const int num_aas = config->get_max_session_aa_idx()+1;

	vector< vector<int> > counts; // size / mobiliy (0 = all);
	
	int i;
	
	counts.resize(num_aas+1);
	for (i=0; i<=num_aas; i++)
		counts[i].resize(4,0);

	for (i=0; i<all_tps.size(); i++)
	{
		int size_idx = (int)(all_tps[i].pm_with_19 * 0.01);

		if (all_tps[i].mobility<MOBILE || all_tps[i].mobility>NONMOBILE)
		{
			cout << "Error: bad mobility value!" << endl;
			exit(1);
		}

		if (size_idx<3 || size_idx>97)
		{
			cout << "Error: bad size_idx!" << endl;
			exit(1);
		}

		int j;
		for (j=0; j<all_tps[i].amino_acids.size(); j++)
		{
			counts[all_tps[i].amino_acids[j]][0]++;
			counts[all_tps[i].amino_acids[j]][all_tps[i].mobility]++;
		}
	}

	for (i=0; i<num_aas; i++)
	{
		int j;
		for (j=0; j<4; j++)
			counts[num_aas][j]+=counts[i][j];
	}

	for (i=0; i<num_aas; i++)
	{
		if (counts[i][0]==0)
			continue;

		cout << config->get_aa2label()[i];
		int j;
		for (j=0; j<4; j++)
		{
			cout << "\t" << counts[i][j] << "\t" << fixed << setprecision(3) << 
				counts[i][j]/(double)counts[num_aas][j];
		}
		cout << endl;	
	}
	
	cout << "--------------------------------------------------------" << endl;
	
	for (i=0; i<4; i++)
		cout << "\t" << counts[num_aas][i] << "\t";
	cout << endl << endl;
}


void fragment_detection_stats(const vector<TrainingPeptide>& all_tps, Config *config)
{
	const int num_frags = config->get_all_fragments().size();

	vector< vector<int> > counts; // size / mobiliy (0 = all);
	vector< vector<int> > totals;
	
	int i;
	
	counts.resize(num_frags);
	for (i=0; i<=num_frags; i++)
		counts[i].resize(4,0);

	totals = counts;

	for (i=0; i<all_tps.size(); i++)
	{
		int size_idx = (int)(all_tps[i].pm_with_19 * 0.01);

		if (all_tps[i].mobility<MOBILE || all_tps[i].mobility>NONMOBILE)
		{
			cout << "Error: bad mobility value!" << endl;
			exit(1);
		}

		if (size_idx<3 || size_idx>97)
		{
			cout << "Error: bad size_idx!" << endl;
			exit(1);
		}

		const TrainingPeptide& tp = all_tps[i];
		int f;
		for (f=0; f<tp.frag_idxs.size(); f++)
		{
			const int frag_idx = tp.frag_idxs[f];
			int a;
			for (a=0; a<tp.intens[f].size(); a++)
			{
				if (tp.intens[f][a]>=0)
				{
					totals[frag_idx][0]++;
					totals[frag_idx][tp.mobility]++;
				}

				if (tp.intens[f][a]>0)
				{
					counts[frag_idx][0]++;
					counts[frag_idx][tp.mobility]++;
				}
			}
		}
	}

	for (i=0; i<num_frags; i++)
	{
	
		if (totals[i][0]==0)
			continue;

		cout << config->get_all_fragments()[i].label;
		int j;
		for (j=0; j<4; j++)
		{
			cout << "\t" << totals[i][j] << "\t" << fixed << setprecision(3) << 
				counts[i][j]/(double)totals[i][j];
		}
		cout << endl;	
	}
	
	cout << "--------------------------------------------------------" << endl;
}



void generate_size_reports()
{
/*	PeakRankModel rm;
	AllScoreModels model;
	vector<FileManager> fms;

	model.read_model("LTQ_LOW_TRYP"); 
	Config *config= model.get_config();

	char ptm_line[] = "C+57:M+16:Q-17";
	config->apply_selected_PTMs(ptm_line);

	rm.set_mass_detection_defaults();

	int charge;
	for (charge=1; charge<=3; charge++)
	{

		vector<TrainingPeptide> tps;

		tps.clear();

		char shew_file[256],hek_file[256];

		sprintf(shew_file,"C:\\Work\\msms5\\NewScore\\tps\\Shew_98_%d_unique_tps.txt",charge);
		sprintf(hek_file,"C:\\Work\\msms5\\NewScore\\tps\\HEK_98_%d_unique_tps.txt",charge);

		read_training_peptides_from_file(shew_file,tps);
		read_training_peptides_from_file(hek_file,tps);

		cout << "Cahrge " << charge << "  Read " << tps.size() << " peptides..." << endl;
		size_mobility_stats(tps);
		cout << endl;
	}
*/
}




