#include "EdgeModel.h"
#include "PrmGraph.h"
#include "AnnotatedSpectrum.h"
#include "AllScoreModels.h"
#include "ME_REG.h"

// parameters for all regional models (magic numbers)
float EdgeModel::weight_single_state_score;
float EdgeModel::weight_single_offset_score;
float EdgeModel::weight_multi_state_score;
float EdgeModel::weight_multi_offset_score;
float EdgeModel::weight_transfer_score;
float EdgeModel::weight_combo_score;
float EdgeModel::multi_aa_penalty;
float EdgeModel::bad_pair_penalty;
float EdgeModel::problematic_pair_penalty;

vector<float> EdgeModel::tol_limits;


int RegionalEdgeModel::calc_tol_bin_idx(const Breakage* n_break, const Breakage *c_break,
						 mass_t exp_edge_mass) const
{
	mass_t total_offset=0;
	int num_offsets=0;
	int i;

	for (i=0; i<strong_frag_idxs.size(); i++)
	{
		const int n_pos = n_break->get_position_of_frag_idx(strong_frag_idxs[i]);
		const int c_pos = c_break->get_position_of_frag_idx(strong_frag_idxs[i]);
		if (n_pos<0 || c_pos<0)
			continue;

		const mass_t n_mass = n_break->fragments[n_pos].mass;
		const mass_t c_mass = c_break->fragments[c_pos].mass;
		const mass_t offset = fabs(fabs((n_mass-c_mass)*strong_frag_charges[i])-exp_edge_mass);
		total_offset+=offset;
		num_offsets++;
	}

	total_offset *= one_over_tolerance;
		
	return EdgeModel::get_tol_bin_idx((float)total_offset);
}



void RegionalEdgeModel::set_params(const Config *config, int c, int s, int r)
{
	charge = c;
	region_idx = r;
	size_idx = s;

	const RegionalFragments& rf = config->get_regional_fragments(c,s,r);
	strong_frag_idxs = rf.get_strong_frag_type_idxs();
	num_states = 1;

	strong_frag_charges.resize(strong_frag_idxs.size());
	int i;
	for (i=0; i<strong_frag_idxs.size(); i++)
	{
		strong_frag_charges[i]=config->get_fragment(strong_frag_idxs[i]).charge;
		num_states *=2;
	}

	one_over_tolerance = 1.0/config->getTolerance();
}


void RegionalEdgeModel::train_regional_edge_model(void *model_ptr, const SpectraAggregator& sa, 
												  SpectraList& sl)
{
	AllScoreModels  *model = (AllScoreModels *)model_ptr;
	const Config *config = model->get_config();
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const vector<int>& org_aa = config->get_org_aa();
	const mass_t tolerance = config->getTolerance();
	const int num_limits = EdgeModel::get_num_tol_limits();

	bool verbose = true;

	vector<double> good_single_state_counts, bad_single_state_counts;
	vector< vector< double > > good_single_off_counts, bad_single_off_counts;
	vector<double> good_multi_state_counts, bad_multi_state_counts;
	vector< vector< double > > good_multi_off_counts, bad_multi_off_counts;

	vector< vector< double > > good_trans_counts, bad_trans_counts;

	vector< vector< double > > good_aa_combo_counts, bad_aa_combo_counts;

	good_single_state_counts.resize(num_states,2.0);
	bad_single_state_counts.resize(num_states,2.0);
	good_multi_state_counts.resize(num_states,2.0);
	bad_multi_state_counts.resize(num_states,2.0);

	int i;
	good_single_off_counts.resize(num_states);
	bad_single_off_counts.resize(num_states);
	good_multi_off_counts.resize(num_states);
	bad_multi_off_counts.resize(num_states);
	for (i=0; i<num_states; i++)
	{
		good_single_off_counts[i].resize(num_limits,2.0);
		bad_single_off_counts[i].resize(num_limits,2.0);
		good_multi_off_counts[i].resize(num_limits,2.0);
		bad_multi_off_counts[i].resize(num_limits,2.0);
	}

	good_trans_counts.resize(num_states);
	bad_trans_counts.resize(num_states);
	for (i=0; i<num_states; i++)
	{
		good_trans_counts[i].resize(num_states,5.0);
		bad_trans_counts[i].resize(num_states,5.0);
	}

	good_aa_combo_counts.resize(Val+1);
	bad_aa_combo_counts.resize(Val+1);
	for (i=0; i<=Val; i++)
	{
		good_aa_combo_counts[i].resize(Val+1,1.0);
		bad_aa_combo_counts[i].resize(Val+1,1.0);
	}


	if (verbose)
	{
		cout << endl;
		cout << "Computing Edge Model for " << charge << " " << size_idx << " " << region_idx << endl << endl;
	}


	double num_good_combos=0;
	double num_bad_combos=0;
	int sc;
	for (sc=0; sc<sl.getNumHeaders(); sc++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(sc);
		Peptide pep;
		pep.parseFromString(config, header->getPeptideStr());
		const vector<int> pep_aas = pep.get_amino_acids();
		const mass_t pm_with_19 = pep.get_mass_with_19();
		
		AnnotatedSpectrum as;
		as.readSpectrum(sa,header);
		as.annotate_spectrum(pm_with_19,header->getCharge(),true);

		PrmGraph prm;
		prm.create_graph_from_spectrum(model,&as,pm_with_19,charge);
		const int num_nodes = prm.get_num_nodes();
		vector<int> good_node_inds;
		vector<mass_t> correct_break_masses;
		pep.calc_expected_breakage_masses(config,correct_break_masses);
		good_node_inds.resize(num_nodes,0);
			
		int i;
		for (i=0; i<correct_break_masses.size(); i++)
		{
			int max_node_idx = prm.get_max_score_node(correct_break_masses[i],tolerance);
			if (max_node_idx>=0)
				good_node_inds[max_node_idx]=1;
	
			max_node_idx = prm.get_max_score_node(pm_with_19 - correct_break_masses[i],tolerance);
			if (max_node_idx>=0)
				good_node_inds[max_node_idx]=2;	
		}

		const vector<Node>& nodes = prm.get_nodes();
		const vector<MultiEdge>& edges = prm.get_multi_edges();

		// collect edge data
		const vector<MultiEdge>& multi_edges = prm.get_multi_edges();
		for (i=0; i<multi_edges.size(); i++)
		{
			const MultiEdge& edge = multi_edges[i];
			if (prm.get_node(edge.n_idx).breakage.region_idx != region_idx)
				continue;

			// state offset data
			bool good_edge=false;
			if (good_node_inds[edge.n_idx] == 1 && good_node_inds[edge.c_idx] == 1)
			{
				good_edge=true;
			}
			else if (good_node_inds[edge.n_idx] == 0 && good_node_inds[edge.c_idx] == 0)
			{
				good_edge=false;
			}
			else
				continue;

			const Breakage *n_break = &prm.get_node(edge.n_idx).breakage;
			const Breakage *c_break = &prm.get_node(edge.c_idx).breakage;

			if (! n_break ||  ! c_break )
				continue;

			const int n_state = calc_state(n_break);
			const int c_state = calc_state(c_break);

	//		if (n_state==0 || c_state == 0)
	//			continue;

			const int intersec_state = this->calc_intersec_state(n_state,c_state);

			if (good_edge)
			{
				if (edge.num_aa==1)
				{
					good_single_state_counts[intersec_state]++;
				}
				else
					good_multi_state_counts[intersec_state]++;
			}
			else
				if (edge.num_aa==1)
				{	
					bad_single_state_counts[intersec_state]++;
				}
				else
					bad_multi_state_counts[intersec_state]++;

			int v;
			for (v=0; v<edge.variant_ptrs.size(); v++)
			{
				int *v_ptr = edge.variant_ptrs[v];
				const int num_aa = *v_ptr++;

				mass_t exp_edge_mass=0;
				int j;
				for (j=0; j<num_aa; j++)
					exp_edge_mass+=aa2mass[v_ptr[j]];

				const int offset_bin = this->calc_tol_bin_idx(n_break,c_break,exp_edge_mass);
				if (good_edge)
				{
					int j;
					for (j=0; j<correct_break_masses.size(); j++)
						if (fabs(correct_break_masses[j]-n_break->mass)<1.0)
								break;
					if (j>=correct_break_masses.size()-1)
					{
						continue;
					}
					else 
						if (  (num_aa == 2 && (v_ptr[0] != pep_aas[j] || v_ptr[1] != pep_aas[j+1])) ||
							  (num_aa == 1 && (v_ptr[0] != pep_aas[j]) )  )
							continue;

					if (num_aa==1)
					{
						good_single_off_counts[intersec_state][offset_bin]++;
					}
					else
						good_multi_off_counts[intersec_state][offset_bin]++;

					good_trans_counts[n_state][c_state]++;
				}
				else
				{
					if (num_aa==1)
					{
						bad_single_off_counts[intersec_state][offset_bin]++;
					}
					else
						bad_multi_off_counts[intersec_state][offset_bin]++;

					bad_trans_counts[n_state][c_state]++;
				}

				// fill reg samples
				if (num_aa==2)
				{
					const int n_aa = v_ptr[0];
					const int c_aa = v_ptr[1];
					// make sure variant is good
					if (good_edge)
					{
						good_aa_combo_counts[org_aa[n_aa]][org_aa[c_aa]]++;
						num_good_combos++;
					}
					else
					{
						bad_aa_combo_counts[org_aa[n_aa]][org_aa[c_aa]]++;
						num_bad_combos++;
					}
				}
			}
		}

		if (sc>0 && sc %500 == 0)
		{
			cout << sc << "/" << sl.getNumHeaders() << " ..." << endl;
		}
	}

	// compute probs
	double n_good_single_state_counts=0, n_bad_single_state_counts=0;
	double n_good_multi_state_counts=0, n_bad_multi_state_counts=0;
	double n_good_trans_counts=0, n_bad_trans_counts=0;
	vector<double> n_good_single_off_counts,n_bad_single_off_counts;
	vector<double> n_good_multi_off_counts, n_bad_multi_off_counts;
	
	int j;
	for (i=0; i<good_single_state_counts.size(); i++)
	{
		n_good_single_state_counts+=good_single_state_counts[i];
		n_good_multi_state_counts+=good_multi_state_counts[i];
	}

	for (i=0; i<bad_single_state_counts.size(); i++)
	{
		n_bad_single_state_counts+=bad_single_state_counts[i];
		n_bad_multi_state_counts+=bad_multi_state_counts[i];
	}

	n_good_single_off_counts.resize(num_states,0);
	n_bad_single_off_counts.resize(num_states,0);
	n_good_multi_off_counts.resize(num_states,0);
	n_bad_multi_off_counts.resize(num_states,0);
	for (i=0; i<num_states; i++)
		for (j=0; j<num_limits; j++)
		{
			n_good_single_off_counts[i]+=good_single_off_counts[i][j];
			n_bad_single_off_counts[i]+=bad_single_off_counts[i][j];
			n_good_multi_off_counts[i]+=good_multi_off_counts[i][j];
			n_bad_multi_off_counts[i]+=bad_multi_off_counts[i][j];
		}
	
	n_good_trans_counts=0;
	n_bad_trans_counts=0;
	for (i=0; i<num_states; i++)
		for (j=0; j<num_states; j++)
		{
			n_good_trans_counts+=good_trans_counts[i][j];
			n_bad_trans_counts+=bad_trans_counts[i][j];
		}

	log_odds_state.resize(num_states);
	multi_log_odds_state.resize(num_states);
	for (i=0; i<num_states; i++)
	{
		log_odds_state[i]=log(good_single_state_counts[i]/n_good_single_state_counts)-
		log(bad_single_state_counts[i]/n_bad_single_state_counts);

		multi_log_odds_state[i]=log(good_multi_state_counts[i]/n_good_multi_state_counts)-
		log(bad_multi_state_counts[i]/n_bad_multi_state_counts);
	}

	log_odds_offset.resize(num_states);
	multi_log_odds_offset.resize(num_states);
	log_odds_transfer.resize(num_states);
	for (i=0; i<num_states; i++)
	{
		log_odds_offset[i].resize(num_limits);
		multi_log_odds_offset[i].resize(num_limits);
		log_odds_transfer[i].resize(num_states);

		int j;
		for (j=0; j<num_limits; j++)
		{
			log_odds_offset[i][j]=log(good_single_off_counts[i][j]/n_good_single_off_counts[i])-
				log(bad_single_off_counts[i][j]/n_bad_single_off_counts[i]);
			if (j>0 && log_odds_offset[i][j]>log_odds_offset[i][j-1])
				log_odds_offset[i][j]=log_odds_offset[i][j-1];

			multi_log_odds_offset[i][j]=log(good_multi_off_counts[i][j]/n_good_multi_off_counts[i])-
				log(bad_multi_off_counts[i][j]/n_bad_multi_off_counts[i]);
			if (j>0 && multi_log_odds_offset[i][j]>multi_log_odds_offset[i][j-1])
				multi_log_odds_offset[i][j]=multi_log_odds_offset[i][j-1];
		}


		for (j=0; j<num_states; j++)
		{
			log_odds_transfer[i][j]=log(good_trans_counts[i][j]/n_good_trans_counts) - 
									log(bad_trans_counts[i][j]/n_bad_trans_counts);
			if ((i==0 || j==0) && log_odds_transfer[i][j]>0)
				log_odds_transfer[i][j]=0;
		}
	}

	has_values = true;

	if (verbose)
	{
		cout << "Single offset:" << endl;
		for (i=0;  i<num_states; i++)
		{
			cout << setprecision(0) << fixed << i << "\t" << n_good_single_off_counts[i] << "/" << 
			n_bad_single_off_counts[i] << "\t" << setprecision(2) << log_odds_state[i] << " : ";
			for (j=0; j<num_limits; j++)
				cout << "\t" << log_odds_offset[i][j];
			cout << endl ;
		}
		cout << endl;

		cout << "Multi offset:" << endl;
		for (i=0;  i<num_states; i++)
		{
			cout << setprecision(0) << fixed << i << "\t" << n_good_multi_off_counts[i] << "/" << 
			n_bad_multi_off_counts[i] << "\t" << setprecision(2) << multi_log_odds_state[i] << " : ";
			for (j=0; j<num_limits; j++)
				cout << "\t" << multi_log_odds_offset[i][j];
			cout << endl ;
		}
		cout << endl;

		cout << "Transfer odds:" << endl;
		cout << "\t";
		for (i=0; i<num_states; i++)
			cout << "  " << i << "\t";
		cout << endl;
		for (i=0; i<num_states; i++)
		{
			cout << i << "\t";
			for (j=0; j<num_states; j++)
				cout << log_odds_transfer[i][j] << "\t";
			cout << endl;
		}

	}
}


bool RegionalEdgeModel::write_model(const char *path) const
{
	if (! has_values)
		return false;

	int i;

	ofstream ofs(path);
	if (! ofs.is_open() || ! ofs.good())
	{
		cout << "Error: couldn't open edge model file for writing: " << path << endl;
		exit(1);
	}

	ofs << charge << " " << size_idx << " " << region_idx << " " << num_states << endl;
	ofs << fixed << setprecision(3);

	ofs << strong_frag_idxs.size();
	for (i=0; i<strong_frag_idxs.size(); i++)
		ofs << " " << strong_frag_idxs[i];
	ofs << endl;

	ofs << strong_frag_charges.size();
	for (i=0; i<strong_frag_charges.size(); i++)
		ofs << " " << strong_frag_charges[i];
	ofs << endl;

	ofs << log_odds_state.size();
	for (i=0; i<log_odds_state.size(); i++)
		ofs << " " << log_odds_state[i];
	ofs << endl;
	
	ofs << log_odds_offset.size() << endl;
	for (i=0; i<log_odds_offset.size(); i++)
	{
		int j;
		ofs << log_odds_offset[i].size();
		for (j=0; j<log_odds_offset[i].size(); j++)
			ofs << " " << log_odds_offset[i][j];
		ofs << endl;
	}

	ofs << multi_log_odds_state.size();
	for (i=0; i<multi_log_odds_state.size(); i++)
		ofs << " " << multi_log_odds_state[i];
	ofs << endl;
	
	ofs << multi_log_odds_offset.size() << endl;
	for (i=0; i<multi_log_odds_offset.size(); i++)
	{
		int j;
		ofs << multi_log_odds_offset[i].size();
		for (j=0; j<multi_log_odds_offset[i].size(); j++)
			ofs << " " << multi_log_odds_offset[i][j];
		ofs << endl;
	}

	ofs << log_odds_transfer.size() << endl;
	for (i=0; i<log_odds_transfer.size(); i++)
	{
		int j;
		ofs << log_odds_transfer[i].size();
		for (j=0; j<log_odds_transfer[i].size(); j++)
			ofs << " " << log_odds_transfer[i][j];
		ofs << endl;
	}

	return true;
}


bool RegionalEdgeModel::read_model(const char *path)
{
	
	ifstream ifs(path);
	if (! ifs.is_open() || ! ifs.good())
		return false;

	ifs >> charge >> size_idx >> region_idx >> num_states;

	int i;
	int ns;
	ifs >> ns;
	strong_frag_idxs.resize(ns);
	for (i=0; i<ns; i++)
		ifs >> strong_frag_idxs[i];

	ifs >> ns;
	strong_frag_charges.resize(ns);
	for (i=0; i<ns; i++)
		ifs >> strong_frag_charges[i];

	ifs >> ns;
	log_odds_state.resize(ns);
	for (i=0; i<ns; i++)
		ifs >> log_odds_state[i];
	
	ifs >> ns;
	log_odds_offset.resize(ns);
	for (i=0; i<ns; i++)
	{
		int ls;
		int j;
		ifs >> ls;
		log_odds_offset[i].resize(ls);
		for (j=0; j<ls; j++)
			ifs >> log_odds_offset[i][j];
	}

	ifs >> ns;
	multi_log_odds_state.resize(ns);
	for (i=0; i<ns; i++)
		ifs >> multi_log_odds_state[i];
	
	ifs >> ns;
	multi_log_odds_offset.resize(ns);
	for (i=0; i<ns; i++)
	{
		int ls;
		int j;
		ifs >> ls;
		multi_log_odds_offset[i].resize(ls);
		for (j=0; j<ls; j++)
			ifs >> multi_log_odds_offset[i][j];
	}

	ifs >> ns;
	log_odds_transfer.resize(ns);

	for (i=0; i<ns; i++)
	{
		int ls;
		int j;
		ifs >> ls;
		log_odds_transfer[i].resize(ls);

		for (j=0; j<ls; j++)
			ifs >> log_odds_transfer[i][j];
	}

	has_values = true;

	return true;
}


string EdgeModel::make_regional_name(const char *model_name, int c, int s, int r) const
{
	ostringstream oss;
	oss << model_name << "_edge_" << c << "_" << s << "_" << r << ".txt";

	return oss.str();
}


void EdgeModel::write_edge_models(const char *model_name)
{
	if (! this->ind_was_initialized)
		return;
			
	const string dir_name = config_->get_resource_dir() + "/" + model_name + "_EDGE" + "/";
	const string file_path = dir_name + string(model_name) + "_edge_model.txt";

	ofstream ofs(file_path.c_str());
	if (! ofs.is_open() || ! ofs.good())
	{
		cout << "Error: couldn't open edge model for writing: " << file_path << endl;
		exit(1);
	}

	ofs << "#WEIGHT_SINGLE_STATE_SCORE " << setprecision(3) << fixed << weight_single_state_score << endl;
	ofs << "#WEIGHT_SINGLE_OFFSET_SCORE " << weight_single_offset_score << endl;
	ofs << "#WEIGHT_MULTI_STATE_SCORE " << weight_multi_state_score << endl;
	ofs << "#WEIGHT_MULTI_OFFSET_SCORE " << weight_multi_offset_score << endl;
	ofs << "#WEIGHT_TRANSFER_SCORE " << weight_transfer_score << endl;
	ofs << "#WEIGHT_COMBO_SCORE " << weight_combo_score << endl;
	ofs << "#MULTY_AA_PENALTY " << multi_aa_penalty << endl;
	ofs << "#BAD_DOUBLE_PENALTY " << bad_pair_penalty << endl;
	ofs << "#PROBLEMATIC_PAIR_PENALTY " << problematic_pair_penalty << endl;
	ofs << "#TOLS " << tol_limits.size();
	int i;
	for (i=0; i<tol_limits.size(); i++)
		ofs << " " << tol_limits[i];
	ofs << endl;

	if (double_combo_scores.size()<Val)
	{
		cout << "Error: double combo table not filled!" << endl;
		exit(1);
	}
	for (i=Ala; i<=Val; i++)
	{
		int j;
		for (j=Ala; j<=Val; j++)
			ofs << double_combo_scores[i][j] << " ";
		ofs << endl;
	}

	ofs.close();

	int num_written=0;
	for (i=0; i<this->regional_edge_models.size(); i++)
	{
		int j;
		for (j=0; j<regional_edge_models[i].size(); j++)
		{
			int k;
			for (k=0; k<regional_edge_models[i][j].size(); k++)
			{
				if (regional_edge_models[i][j][k] && regional_edge_models[i][j][k]->has_values)
				{
					string path = dir_name + make_regional_name(model_name,i,j,k);
					if (regional_edge_models[i][j][k]->write_model(path.c_str()))
						num_written++;
				}
			}
		}
	}

	cout << endl << "Wrote " << num_written << " edge models" << endl;
}

void EdgeModel::read_edge_models(const Config *config, const char *model_name, bool silent_ind)
{
	config_ = config;

	const string dir_name = config->get_resource_dir() + "/" + model_name + "_EDGE" + "/";
	const string file_path = dir_name + string(model_name) + "_edge_model.txt";

	

	ifstream ifs(file_path.c_str());
	if (! ifs.is_open() || ! ifs.good())
	{
		cout << "Error: couldn't open edge model for reading: " << file_path << endl;
		exit(1);
	}

	int i,num_tols = 0;
	string dummy;
	ifs >> dummy >> weight_single_state_score;
	ifs >> dummy >> weight_single_offset_score;
	ifs >> dummy >> weight_multi_state_score;
	ifs >> dummy >> weight_multi_offset_score;
	ifs >> dummy >> weight_transfer_score;
	ifs >> dummy >> weight_combo_score;
	ifs >> dummy >> multi_aa_penalty;
	ifs >> dummy >> bad_pair_penalty;
	ifs >> dummy >> problematic_pair_penalty;
	ifs >> dummy >> num_tols;

	tol_limits.resize(num_tols);
	for (i=0; i<num_tols; i++)
		ifs >> tol_limits[i];

	// read double combo table
	double_combo_scores.clear();
	double_combo_scores.resize(Val+1);
	for (i=0; i<=Val; i++)
		double_combo_scores[i].resize(Val+1,0);

	char buff[1024];
	for (i=Ala; i<=Val; i++)
	{
		ifs.getline(buff,1024);
		istringstream iss(buff);
		int j;
		for (j=Ala; j<=Val; j++)
			iss >> double_combo_scores[i][j];
	}

	const vector< vector< mass_t > >& size_thresholds = config->get_size_thresholds();

	// resize according to threhsolds, assume that the last mass is POS_INF
	regional_edge_models.resize(size_thresholds.size());

	ind_was_initialized = true; // ignore the regional files
	return;
	
	int num_read=0;
	int c;
	for (c=0; c<size_thresholds.size(); c++)
	{
		if (size_thresholds[c].size()>0)
			regional_edge_models[c].resize(size_thresholds[c].size());

		int s;
		for (s=0; s<regional_edge_models[c].size(); s++)
		{
			const int num_regions = config->get_num_regions(c,s);
			regional_edge_models[c][s].resize(num_regions,NULL);

			int r;
			for (r=0; r<regional_edge_models[c][s].size(); r++)
			{
				string path = dir_name + make_regional_name(model_name,c,s,r);
				
				ifstream ifs(path.c_str());
				if (ifs.is_open() && ifs.good())
				{
					ifs.close();
			
					regional_edge_models[c][s][r] = new RegionalEdgeModel;
					if (regional_edge_models[c][s][r]->read_model(path.c_str()))
						num_read++;
				}
				else
					ifs.close();
			}
		}
	}
	if (! silent_ind)
		cout << "Read " << num_read << " edge models." << endl;

	ind_was_initialized = true;
}




void EdgeModel::train_all_edge_models(const SpectraAggregator& sa, void *model_ptr, int specific_charge)
{
	AllScoreModels *model = (AllScoreModels *)model_ptr;
	config_ = model->get_config();
	
	const vector< vector< mass_t > >& size_thresholds = config_->get_size_thresholds();

	// resize according to threhsolds, assume that the last mass is POS_INF
	init_edge_model_defaults();
	regional_edge_models.resize(size_thresholds.size());

	int c;
	for (c=0; c<size_thresholds.size(); c++)
		if (size_thresholds[c].size()>0)
			regional_edge_models[c].resize(size_thresholds[c].size());


	for (c=0; c<regional_edge_models.size(); c++)
	{
		if (specific_charge>0 && c != specific_charge)
			continue;

		int s;
		for (s=0; s<regional_edge_models[c].size(); s++)
		{
			const int num_regions = config_->get_num_regions(c,s);
			regional_edge_models[c][s].resize(num_regions,NULL);

			// train all regions
			SpectraList sl(sa);
			mass_t min_mz=0;
			mass_t max_mz = 10000;
			if (s>0)
				min_mz = (size_thresholds[c][s-1]+c-1)/(mass_t)c;
			if (s< size_thresholds[c].size())
				max_mz = (size_thresholds[c][s]+c-1)/(mass_t)c;

			sl.selectHeaders(min_mz, max_mz, c, c);

			if (sl.getNumHeaders()<300)
			{
				cout << "Warning: not enough spectra for charge " << c << " size " << s << ", not training edge models!" << endl;
				continue;
			}
			sl.randomlyReduceListToSize(20000);

			// set region frag_idxs
			int r;
			for (r=0; r<num_regions; r++)
			{
				regional_edge_models[c][s][r] = new RegionalEdgeModel;
				regional_edge_models[c][s][r]->set_params(config_,c,s,r);
				regional_edge_models[c][s][r]->train_regional_edge_model(model,sa,sl);
			}
		}
	}

	// train double combo model
	if (1)
	{
		SpectraList sl(sa);
		sl.selectAllAggregatorHeaders();
		sl.randomlyReduceListToSize(250000);

		AllScoreModels  *model = (AllScoreModels *)model_ptr;
		const Config* config = model->get_config();
		const mass_t tolerance = config->getTolerance();
		const vector<mass_t>& aa2mass = config->get_aa2mass();
		const vector<int>& org_aa = config->get_org_aa();
	
		bool verbose = true;

		vector< vector< double > > good_aa_combo_counts, bad_aa_combo_counts;

		int i;
		good_aa_combo_counts.resize(Val+1);
		bad_aa_combo_counts.resize(Val+1);
		for (i=0; i<=Val; i++)
		{
			good_aa_combo_counts[i].resize(Val+1,1.0);
			bad_aa_combo_counts[i].resize(Val+1,1.0);
		}

		if (verbose)
		{
			cout << endl;
			cout << "Computing double edge model..." << endl;
		}

		double num_good_combos=0;
		double num_bad_combos=0;
		int sc;
		for (sc=0; sc<sl.getNumHeaders(); sc++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(sc);
			Peptide pep;
			pep.parseFromString(config, header->getPeptideStr());
			const vector<int> pep_aas = pep.get_amino_acids();
			const mass_t pm_with_19 = pep.get_mass_with_19();
			const int charge = header->getCharge();

			AnnotatedSpectrum as;
			as.readSpectrum(sa,header);
			as.annotate_spectrum(pm_with_19,charge,true);

			PrmGraph prm;
			prm.create_graph_from_spectrum(model,&as,pm_with_19,charge);
			const int num_nodes = prm.get_num_nodes();
			vector<int> good_node_inds;
			vector<mass_t> correct_break_masses;
			pep.calc_expected_breakage_masses(config,correct_break_masses);
			good_node_inds.resize(num_nodes,0);
				
			int i;
			for (i=0; i<correct_break_masses.size(); i++)
			{
				int max_node_idx = prm.get_max_score_node(correct_break_masses[i],tolerance);
				if (max_node_idx>=0)
					good_node_inds[max_node_idx]=1;
		
				max_node_idx = prm.get_max_score_node(pm_with_19 - correct_break_masses[i],tolerance);
				if (max_node_idx>=0)
					good_node_inds[max_node_idx]=2;	
			}

			const vector<Node>& nodes = prm.get_nodes();
			const vector<MultiEdge>& edges = prm.get_multi_edges();

			// collect edge data
			const vector<MultiEdge>& multi_edges = prm.get_multi_edges();
			for (i=0; i<multi_edges.size(); i++)
			{
				const MultiEdge& edge = multi_edges[i];
				if (edge.num_aa !=2)
					continue;

				// state offset data
				bool good_edge=false;
				if (good_node_inds[edge.n_idx] == 1 && good_node_inds[edge.c_idx] == 1)
				{
					good_edge=true;
				}
				else if (good_node_inds[edge.n_idx] == 0 && good_node_inds[edge.c_idx] == 0)
				{
					good_edge=false;
				}
				else
					continue;

				const Breakage *n_break = &prm.get_node(edge.n_idx).breakage;
				const Breakage *c_break = &prm.get_node(edge.c_idx).breakage;

				if (! n_break ||  ! c_break )
					continue;

				int v;
				for (v=0; v<edge.variant_ptrs.size(); v++)
				{
					int *v_ptr = edge.variant_ptrs[v];
					const int num_aa = *v_ptr++;

					if (good_edge)
					{
						int j;
						for (j=0; j<correct_break_masses.size(); j++)
							if (fabs(correct_break_masses[j]-n_break->mass)<1.0)
									break;
						if (j>=correct_break_masses.size()-1)
						{
							continue;
						}
						else 
							if ( num_aa == 2 && (v_ptr[0] != pep_aas[j] || v_ptr[1] != pep_aas[j+1])) 
								continue;
					}

					const int n_aa = v_ptr[0];
					const int c_aa = v_ptr[1];

					// make sure variant is good
					if (good_edge)
					{
						good_aa_combo_counts[org_aa[n_aa]][org_aa[c_aa]]++;
						num_good_combos++;
					}
					else
					{
						bad_aa_combo_counts[org_aa[n_aa]][org_aa[c_aa]]++;
						num_bad_combos++;
					}
				}
			}

			if (sc>0 && sc %1000 == 0)
			{
				cout << sc << "/" << sl.getNumHeaders() << " ..." << endl;
			}
		}
	

		cout << "Combo table:" << endl;
		double_combo_scores.resize(Val+1);


		for (i=Ala; i<=Val; i++)
		{
			double_combo_scores[i].resize(Val+1,0);
			if (i==Ile)
				continue;
			int j;
			for (j=Ala; j<=Val; j++)
				if (j!=Ile)
					double_combo_scores[i][j]=log(good_aa_combo_counts[i][j]/num_good_combos)-
											  log(bad_aa_combo_counts[i][j]/num_bad_combos);

		}

		const vector<string>& aa2label = config->get_aa2label();
		cout << fixed << setprecision(1) << setw(5);
		cout << "      ";
		for (i=Ala; i<=Val; i++)
			cout << aa2label[i] << "     ";
		cout << endl;

		for (i=Ala; i<=Val; i++)
		{
			cout << aa2label[i] << "     ";
			int j;
			for (j=Ala; j<=Val; j++)
			{
				printf("%.1f  ",double_combo_scores[i][j]);
				if (double_combo_scores[i][j]>=0)
					printf(" ");
			}
			cout << endl;
		}	
	}

	ind_was_initialized = true;
}



void EdgeModel::init_edge_model_defaults()
{
	weight_single_state_score  = 1.0;
	weight_single_offset_score = 1.0;
	weight_multi_state_score   = 1.0;
	weight_multi_offset_score  = 1.0;
	weight_transfer_score      = 1.0;
	weight_combo_score         = 3.0;
	multi_aa_penalty           = -8.0;
	bad_pair_penalty		   = -5.0;
	problematic_pair_penalty   = -2.0;

	tol_limits.clear();
	tol_limits.push_back(0.1);
	tol_limits.push_back(0.25);
	tol_limits.push_back(0.5);
	tol_limits.push_back(0.9);
	tol_limits.push_back(POS_INF);

	ind_was_initialized = false;
}

// assigns probs and scores to graph edges
void EdgeModel::score_graph_edges(PrmGraph& prm) const
{
	if (! this->ind_was_initialized)
		return;

	const vector<int>& org_aa = config_->get_org_aa();
	const vector<mass_t>& aa2mass = config_->get_aa2mass();
	const vector< vector< bool > >& allowed_double_edge = config_->get_allowed_double_edge();
	const vector< vector< bool > >& double_edge_with_same_mass_as_single = config_->get_double_edge_with_same_mass_as_single();
	const int charge = prm.charge;
	const int size_idx = prm.size_idx;

	const vector<Node>& nodes = prm.nodes;
	const int num_edges = prm.multi_edges.size();
	int i;


	for (i=0; i<num_edges; i++)
	{
		MultiEdge& edge = prm.multi_edges[i];
		const Breakage *n_break = &prm.get_node(edge.n_idx).breakage;
		const Breakage *c_break = &prm.get_node(edge.c_idx).breakage;

		const bool add_only_positive_state_scores = (edge.n_idx == 0 || edge.c_idx == nodes.size()-1 || 
												     nodes[edge.c_idx].type == NODE_DIGEST);
		
		edge.max_variant_score = NEG_INF;

		int v;
		for (v=0; v<edge.variant_ptrs.size(); v++)
		{
			int *v_ptr = edge.variant_ptrs[v];
			const int num_aa = *v_ptr++;

			float edge_score = 0;

			// this is in the EdgeModel, common to all charge/size/regions
			if (num_aa>1)
			{
				edge_score += multi_aa_penalty;
				if (num_aa>2)
					edge_score += (num_aa-2)*0.5*multi_aa_penalty;

				const int n_aa = v_ptr[0];
				const int c_aa = v_ptr[num_aa-1];
				const int org_n_aa = org_aa[n_aa];
				const int org_c_aa = org_aa[c_aa];
				edge_score += weight_combo_score * double_combo_scores[org_n_aa][org_c_aa];

				if (! allowed_double_edge[org_n_aa][org_c_aa])
					edge_score += bad_pair_penalty;

				if (double_edge_with_same_mass_as_single[org_n_aa][org_c_aa])
					edge_score += problematic_pair_penalty;
			}

			edge.variant_scores[v] += edge_score;

			if (edge.variant_scores[v]>edge.max_variant_score)
				edge.max_variant_score = edge.variant_scores[v];
		}
	}

/*	for (i=0; i<num_edges; i++)
	{
		MultiEdge& edge = prm.multi_edges[i];
		const Breakage *n_break = &prm.get_node(edge.n_idx).breakage;
		const Breakage *c_break = &prm.get_node(edge.c_idx).breakage;
		const int region_idx = (n_break ? n_break->region_idx : 0);
		const RegionalEdgeModel* rem = regional_edge_models[charge][size_idx][region_idx];
		if (! rem)
		{
			cout << "Error: no edge model for " << charge << " " << size_idx << " " << region_idx << "!" << endl;
			exit(1);
		}

		const int n_state = (n_break ? rem->calc_state(n_break) : 0);
		const int c_state = (c_break ? rem->calc_state(c_break) : 0);		
		const int intersec_state = rem->calc_intersec_state(n_state,c_state);

		const bool add_only_positive_state_scores = (edge.n_idx == 0 || edge.c_idx == nodes.size()-1 || 
												     nodes[edge.c_idx].type == NODE_DIGEST);
		
		edge.max_variant_score = NEG_INF;

		int v;
		for (v=0; v<edge.variant_ptrs.size(); v++)
		{
			int *v_ptr = edge.variant_ptrs[v];
			const int num_aa = *v_ptr++;

			float edge_score = 0;

			// this portion needs to be trained for each charge/size/region
			if (rem->has_values)
			{
				mass_t exp_edge_mass=0;
				int j;
				for (j=0; j<num_aa; j++)
					exp_edge_mass+=aa2mass[v_ptr[j]];

				const int offset_bin = rem->calc_tol_bin_idx(n_break,c_break,exp_edge_mass);
			
				if (num_aa == 1)
				{
					const float add_score = weight_single_state_score * rem->log_odds_state[intersec_state] +
											weight_single_offset_score * rem->log_odds_offset[intersec_state][offset_bin];

					if (! add_only_positive_state_scores || add_score>0)
						edge_score += add_score;
					
				}
				else
				{
					edge_score += multi_aa_penalty;
					if (num_aa>2)
						edge_score += (num_aa-2)*0.5*multi_aa_penalty;

					const float add_score = weight_multi_state_score * rem->multi_log_odds_state[intersec_state] +
											weight_multi_offset_score * rem->multi_log_odds_offset[intersec_state][offset_bin];

					if (! add_only_positive_state_scores || add_score>0)
						edge_score += add_score;
				}
				edge_score += weight_transfer_score * rem->log_odds_transfer[n_state][c_state];
			}

			// this is in the EdgeModel, common to all charge/size/regions
			if (num_aa>1)
			{
				const int n_aa = v_ptr[0];
				const int c_aa = v_ptr[num_aa-1];
				const int org_n_aa = org_aa[n_aa];
				const int org_c_aa = org_aa[c_aa];
				edge_score += weight_combo_score * double_combo_scores[org_n_aa][org_c_aa];

				if (! allowed_double_edge[org_n_aa][org_c_aa])
					edge_score += bad_pair_penalty;

				if (double_edge_with_same_mass_as_single[org_n_aa][org_c_aa])
					edge_score += problematic_pair_penalty;
			}

			edge.variant_scores[v] += edge_score;

			if (edge.variant_scores[v]>edge.max_variant_score)
				edge.max_variant_score = edge.variant_scores[v];
		}
	}	*/
}


