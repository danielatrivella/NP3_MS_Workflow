#include "PrmNodeScoreModel.h"
#include "PrmGraph.h"
#include "AllScoreModels.h"

void PrmNodeScoreModel::read_score_model(Config *config, const char *name, bool silent_ind)
{
	config_ = config;

	// resize regional breakage score models according to regional fragment sets
	const vector< vector< vector< RegionalFragments > > >& all_rfs = config_->get_regional_fragment_sets();

	int c;
	RegionalPrmNodeScoreModels_.resize(all_rfs.size());
	for (c=0; c<all_rfs.size(); c++)
	{
		RegionalPrmNodeScoreModels_[c].resize(all_rfs[c].size());
		int s;
		for (s=0; s<all_rfs[c].size(); s++)
		{
			RegionalPrmNodeScoreModels_[c][s].resize(all_rfs[c][s].size());
			int r;
			for (r=0; r<RegionalPrmNodeScoreModels_[c][s].size(); r++)
				if (! RegionalPrmNodeScoreModels_[c][s][r].get_was_initialized())
					RegionalPrmNodeScoreModels_[c][s][r].init(config_, c, s, r);
		}
	}
	
	int num_models_read=0;
	for (c=0; c<RegionalPrmNodeScoreModels_.size(); c++)
	{
		int s;
		for (s=0; s<RegionalPrmNodeScoreModels_[c].size(); s++)
		{
			int r;
			for (r=0; r<RegionalPrmNodeScoreModels_[c][s].size(); r++)
			{
				if (RegionalPrmNodeScoreModels_[c][s][r].read_regional_score_model(name,silent_ind))
				{
					if (! silent_ind)
						cout << "Read breakage model: " << c << " " << s << " " << r << endl;
					num_models_read++;
				}
			}
		}
	}

	if (! silent_ind)
		cout << "Read " << num_models_read << " regional breakage models.." << endl;

	// reaf prm_normalize model file if it exists
	read_prm_normalizer_values();

	ind_was_initialized=true;

}
	
void PrmNodeScoreModel::write_score_model(const char *name) const
{
	int num_models_read=0;
	int c;
	for (c=0; c<RegionalPrmNodeScoreModels_.size(); c++)
	{
		int s;
		for (s=0; s<RegionalPrmNodeScoreModels_[c].size(); s++)
		{
			int r;
			for (r=0; r<RegionalPrmNodeScoreModels_[c][s].size(); r++)
			{
				if (RegionalPrmNodeScoreModels_[c][s][r].get_was_initialized())
					RegionalPrmNodeScoreModels_[c][s][r].write_regional_score_model(name);
				
			}
		}
	}
}


struct NodeType {
	int   type;
	int   region;
	float org_score;
	float mod_score;
};

/*
void PrmNodeScoreModel::learn_prm_normalizer_values(const FileManager& fm)
{
	const float step = 0.5;
	const float min_delta = -1.0;
	const float max_delta = 7.0;
	const float target_mid_ratio = 0.96;
	const float target_side_ratio = 0.94;


	config_->set_use_spectrum_charge(1);

	regional_prm_normalizers.resize(RegionalPrmNodeScoreModels_.size());
	int c;
	for (c=0; c<RegionalPrmNodeScoreModels_.size(); c++)
	{
		regional_prm_normalizers[c].resize(RegionalPrmNodeScoreModels_[c].size());
		int s;
		for (s=0; s<RegionalPrmNodeScoreModels_[c].size(); s++)
			regional_prm_normalizers[c][s].resize(RegionalPrmNodeScoreModels_[c][s].size(),0);
	}
	

	const vector< vector<mass_t> >& mass_threshes = config_->get_size_thresholds();
	for (c=1; c<regional_prm_normalizers.size(); c++)
	{
		int s;
		for (s=0; s<regional_prm_normalizers[c].size(); s++)
		{
			const mass_t min_mass = (s == 0 ? 0 : mass_threshes[c][s-1]);
			const mass_t max_mass =  mass_threshes[c][s];
			const int num_regions = regional_prm_normalizers[c][s].size();
			
			cout << "Finding normalizers for charge " << c << " size " << s << "  (masses " << min_mass << " - " <<
				max_mass << ")" << endl;

			FileSet fs;
			BasicSpecReader bsr;

			fs.select_files_in_mz_range(fm,min_mass/c,max_mass/c,c);
			fs.randomly_reduce_ssfs(2000);

			vector< vector< NodeType > > all_prms;
			const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

			if (fs.get_total_spectra()<50)
			{
				cout << "Insufficient number of spectra... skipping" << endl;
				continue;
			}

			int sc;
			for (sc=0; sc<all_ssf.size(); sc++)
			{
				PrmGraph prm;
				static vector<QCPeak> peaks;
				SingleSpectrumFile *ssf = all_ssf[sc];
				if (peaks.size()<ssf->num_peaks)
				{
					int new_size = ssf->num_peaks*2;
					if (new_size<2500)
						new_size=2500;
					peaks.resize(new_size); 
				}

//				const int num_peaks = bsr.read_basic_spec(config_,fm,ssf,&peaks[0]);	
//				if (num_peaks<5)
//					continue;

				// convert peak list ot a spectrum with charge (if original charge ==0)
				// the spectrum gets charge 2, but the true charge is computed from the data
			
				// TODO fix the de novo stuff, need to somehow send a model, or create prmgraph
				// without the extra scores
				Spectrum s;
				s.init_from_QCPeaks(config_,&peaks[0],num_peaks,ssf);

				vector<mass_t> pms_with_19;
				vector<int>    charges;

				pms_with_19.clear();
				charges.clear();		
				
				BasicSpectrum bs;
				bs.ssf = ssf;
				bs.peaks = &peaks[0];
				bs.num_peaks = num_peaks;

				// output m/z and prob values for the different charge states
				
				select_pms_and_charges(&config,bs,pms_with_19,charges);
				if (pms_with_19.size()<=0)
					continue;
			
				s.setCharge(charges[0]);
				init_model_for_scoring_spectrum(&s);
		
				// TO DO fix this issue, it needs the AllScoreModels class to create a spectrum
				// prm.create_graph_from_spectrum(this,&s,pms_with_19[0]);

				vector<NodeType> spec_prms;
				vector<mass_t>   exp_masses;
				const mass_t true_mass_with_19 = s.get_true_mass_with_19();
				s.getPeptide().calc_expected_breakage_masses(&config,exp_masses);

				int i;
				for (i=1; i<prm.get_num_nodes()-1; i++)
				{
					const Node& node = prm.get_node(i);
					if (node.score == 0)
						continue;
					
					NodeType nt;

					nt.type = 0;
					int j;
					for (j=0; j<exp_masses.size(); j++)
						if (fabs(exp_masses[j]-node.mass)<config.getTolerance())
						{
							nt.type=1;
							break;
						}
					
					if (nt.type<=0)
					{
						int j;
						for (j=0; j<exp_masses.size(); j++)
							if (fabs(true_mass_with_19 - exp_masses[j] -node.mass-MASS_PROTON)<config.getTolerance())
							{
								nt.type=2;
								break;
							}
					}

					nt.org_score = node.score;
					nt.mod_score = node.score;
					nt.region = node.breakage.region_idx;
					spec_prms.push_back(nt);
				}
				all_prms.push_back(spec_prms); 
			}
		
	
			vector< vector< double > > per_pre, per_suf, per_covered;
			vector<float> deltas;

			per_pre.resize(num_regions);
			per_suf.resize(num_regions);
			per_covered.resize(num_regions);

			float delta;
			for (delta = min_delta; delta<=max_delta; delta+=step )
			{
				// perform mods
				int a;
				for (a=0; a<all_prms.size(); a++)
				{
					int b;
					for (b=0; b<all_prms[a].size(); b++)
					{
						NodeType& nt = all_prms[a][b];
						if (nt.org_score< -delta)
						{
							nt.mod_score = NEG_INF;
							continue;
						}
						
				
						nt.mod_score = nt.org_score + delta;
					}
				}

				// compute stats (if score is negative treat as 0)
				vector<double> num_pre,num_suf;
				vector<double> num_pre_wpos, num_suf_wpos;
				vector<double> score_pre, score_suf, total_score;
			

				num_pre.resize(num_regions,0);
				num_suf.resize(num_regions,0);
				num_pre_wpos.resize(num_regions,0);
				num_suf_wpos.resize(num_regions,0);
				score_pre.resize(num_regions,0);
				score_suf.resize(num_regions,0);
				total_score.resize(num_regions,0);
				
				for (a=0; a<all_prms.size(); a++)
				{
					int b;
					for (b=0; b<all_prms[a].size(); b++)
					{
						const int   type =    all_prms[a][b].type;
						const float score =   all_prms[a][b].mod_score;
						const int   region =  all_prms[a][b].region;

						if (type == 1)
						{
							num_pre[region]++;
							if (score>0)
							{
								num_pre_wpos[region]++;
								score_pre[region]+= score;
							}
						}

						if (type == 2)
						{
							num_suf[region]++;
							if (score>0)
							{
								num_suf_wpos[region]++;
								score_suf[region]+=score;
							}
						}

						if (score>0)
							total_score[region]+=score;
					}
				}

				
				deltas.push_back(delta);
				int r;
				for (r=0; r<num_regions; r++)
				{
					per_pre[r].push_back(num_pre_wpos[r]/num_pre[r]);
					per_suf[r].push_back(num_suf_wpos[r]/num_suf[r]);
					per_covered[r].push_back((score_pre[r]+score_suf[r])/total_score[r]);
				}
			}

			// report
			int r;
			for (r=0; r<num_regions; r++)
			{
				cout << endl << "Region " << r << endl;
				int d;
				for (d=0; d<deltas.size(); d++)
					cout << "\t" << deltas[d];
				cout << endl << "% Pre";
				for (d=0; d<per_pre[r].size(); d++)
					cout << "\t" << per_pre[r][d];
				cout << endl << "% Suf";
				for (d=0; d<per_suf[r].size(); d++)
					cout << "\t" << per_suf[r][d];
				cout << endl << "% Cov";
				for (d=0; d<per_covered[r].size(); d++)
					cout << "\t" << per_covered[r][d];
				cout << endl;

				// select
				float target_val = target_mid_ratio;
				if (r==0 || r == num_regions-1)
					target_val = target_side_ratio;

				float best_val=POS_INF;
				float best_delta=0;

				for (d=0; d<deltas.size(); d++)
					if (fabs(per_pre[r][d]-target_val)<best_val)
					{
						best_val = fabs(per_pre[r][d]-target_val);
						best_delta = deltas[d];
					}
				
				cout << "Chose delta = " << best_delta << endl << endl;
				regional_prm_normalizers[c][s][r]=best_delta;
			}	
		}
	}
}
*/


bool PrmNodeScoreModel::read_prm_normalizer_values()
{

	string file_path = config_->get_resource_dir() + "/" + config_->getModelName() + "_prm_norm.txt";

	regional_prm_normalizers.resize(RegionalPrmNodeScoreModels_.size());
	int c;
	for (c=0; c<RegionalPrmNodeScoreModels_.size(); c++)
	{
		regional_prm_normalizers[c].resize(RegionalPrmNodeScoreModels_[c].size());
		int s;
		for (s=0; s<RegionalPrmNodeScoreModels_[c].size(); s++)
			regional_prm_normalizers[c][s].resize(RegionalPrmNodeScoreModels_[c][s].size(),0);
	}


	ifstream in_stream(file_path.c_str());
	if (! in_stream.good())
	{
	//	cout << "Error: couldn't open prm_norm file for reading: " << file_path << endl;
	//	exit(1);
		return false;
	}


	while (! in_stream.eof())
	{
		char line_buff[256];
		in_stream.getline(line_buff,256);

		istringstream iss(line_buff);
		int c=0,s=-1,n=0;
		iss >> c >> s >> n;

		if (c<=0 || s<0 || n<=0)
			continue;

		if (c > regional_prm_normalizers.size() ||
			s > regional_prm_normalizers[c].size() ||
			n > regional_prm_normalizers[c][s].size() )
		{
			cout << "Error: mismatch between regional models and prm normalizers: " << c << " " << s << " " << n << endl;
			exit(1);
		}

		int r;
		for (r=0; r<n; r++)
		{
			iss >> regional_prm_normalizers[c][s][r];
		}
	}

	in_stream.close();
	return true;
}

void PrmNodeScoreModel::write_prm_normalizer_values() const
{
	string file_path = config_->get_resource_dir() + "/" + config_->getModelName() + "_prm_norm.txt";

	ofstream out_stream(file_path.c_str());
	if (! out_stream.good())
	{
		cout << "Error: couldn't open prm_norm file for writing: " << file_path << endl;
		exit(1);
	}

	out_stream << setprecision(3);
	int c;
	for (c=1; c<regional_prm_normalizers.size(); c++)
	{
		if (regional_prm_normalizers[c].size()==0)
			continue;
		int s;
		for (s=0; s<regional_prm_normalizers[c].size(); s++)
		{
			out_stream << c << " " << s << " " << regional_prm_normalizers[c][s].size();
			int r;
			for (r=0; r<regional_prm_normalizers[c][s].size(); r++)
				out_stream << " " << regional_prm_normalizers[c][s][r];
			out_stream << endl;
		}
	}
	out_stream.close();
}


/****************************************************************************
This function should be used for adjusting the scores of the PrmGraph
when printing of the node scores only
*****************************************************************************/
void PrmNodeScoreModel::normalize_prm_scores(PrmGraph &prm) const
{
	const int charge   = prm.get_charge();
	const int size_idx = prm.get_size_idx();
	const vector<score_t>& region_norms = this->regional_prm_normalizers[charge][size_idx];

	const int num_nodes = prm.get_num_nodes();
	int i;
	for (i=1; i<num_nodes-1; i++)
	{
		Node& node = prm.get_non_const_node(i);

		if (node.type == NODE_DIGEST)
			continue;

		const int region = node.breakage.region_idx;
		const score_t norm_val = region_norms[region];
		if (node.score<-norm_val)
		{
			node.score=NEG_INF;
			continue;
		}

		if  (node.score>norm_val)
			continue;

		node.score += (norm_val - node.score)*0.5;
	}
	
}

// simple Dancik-style score
void PrmNodeScoreModel::score_breakage(Spectrum *spec, Breakage *breakage, bool verbose) const
{
	const RegionalPrmNodeScoreModel& breakage_model = RegionalPrmNodeScoreModels_[breakage->parent_charge]
														  [breakage->parent_size_idx][breakage->region_idx];

	if (0 || breakage_model.get_has_all_breakage_models())
	{
	}
	else
	{
		const vector<int>& frag_type_idxs = breakage_model.get_frag_type_idxs();
		const vector<score_t>& frag_inten_scores = breakage_model.get_frag_inten_scores();
		const vector<score_t>& frag_no_inten_scores = breakage_model.get_frag_no_inten_scores();
		score_t breakage_score=0;
		int i;
		for (i=0; i<frag_type_idxs.size(); i++)
		{
			const int frag_idx = frag_type_idxs[i];
			if (breakage->is_frag_type_visible(frag_idx))
			{
				if (breakage->get_position_of_frag_idx(frag_idx)>=0)
				{
					breakage_score += frag_inten_scores[frag_idx];
				}
				else
					breakage_score += frag_no_inten_scores[frag_idx];
			}
		}
		breakage->score = breakage_score;
	}
}

/*
void PrmNodeScoreModel::train_score_model(const char *name, const FileManager& fm, 
						int charge, int size_idx, int region_idx)
{
	// resize regional breakage score models according to regional fragment sets
	const vector< vector< vector< RegionalFragments > > >& all_rfs = config_->get_regional_fragment_sets();

	int c;
	RegionalPrmNodeScoreModels_.resize(all_rfs.size());
	for (c=0; c<all_rfs.size(); c++)
	{
		RegionalPrmNodeScoreModels_[c].resize(all_rfs[c].size());
		int s;
		for (s=0; s<all_rfs[c].size(); s++)
		{
			RegionalPrmNodeScoreModels_[c][s].resize(all_rfs[c][s].size());
			int r;
			for (r=0; r<RegionalPrmNodeScoreModels_[c][s].size(); r++)
				if (! RegionalPrmNodeScoreModels_[c][s][r].get_was_initialized())
					RegionalPrmNodeScoreModels_[c][s][r].init(config_,c,s,r);
		}
	}


	// train models
	for (c=1; c<RegionalPrmNodeScoreModels_.size(); c++)
	{
		if (RegionalPrmNodeScoreModels_.size() == 0 || (charge>0 && charge != c))
			continue;

		if (fm.get_num_spectra(c)<200)
		{
			cout << "WARNING: insufficient number of spectra to train breakage model for charge " << c << endl;
			cout <<	"		  only " << fm.get_num_spectra(c) << " spectra were found so this charge is being skipped!" << endl << endl;
			continue;
		}

		int s;
		for (s=0; s<RegionalPrmNodeScoreModels_[c].size(); s++)
		{
			if (size_idx>=0 && s != size_idx)
				continue;

			int r;
			for (r=0; r<RegionalPrmNodeScoreModels_[c][s].size(); r++)
			{
				if (region_idx>=0 && r != region_idx)
					continue;
				
				RegionalPrmNodeScoreModels_[c][s][r].train_regional_score_model(this,name,fm);
			}
		}
	}

	// train PRM normalizer values
	cout << endl << "Training PRM normalizer vlaues..." << endl;

// TODO fix this issue, it needs to use the AllScoreModels class
//	learn_prm_normalizer_values(fm);

	ind_was_initialized=true;
}
*/


void PrmNodeScoreModel::trainNodeScoreModels(void* allScoreModelsVoidPointer,
											 const char *name, 
											 const SpectraAggregator& sa,
											 int specificCharge, 
											 int specificSize, 
											 int specificRegion)
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsVoidPointer);
	config_ = allScoreModels->get_config();
	// resize regional breakage score models according to regional fragment sets
	const vector< vector< vector< RegionalFragments > > >& all_rfs = config_->get_regional_fragment_sets();

	int c;
	RegionalPrmNodeScoreModels_.resize(all_rfs.size());
	for (c=0; c<all_rfs.size(); c++)
	{
		RegionalPrmNodeScoreModels_[c].resize(all_rfs[c].size());
		int s;
		for (s=0; s<all_rfs[c].size(); s++)
		{
			RegionalPrmNodeScoreModels_[c][s].resize(all_rfs[c][s].size());
			int r;
			for (r=0; r<RegionalPrmNodeScoreModels_[c][s].size(); r++)
				if (! RegionalPrmNodeScoreModels_[c][s][r].get_was_initialized())
					RegionalPrmNodeScoreModels_[c][s][r].init(config_,c,s,r);
		}
	}


	// train models
	for (c=1; c<RegionalPrmNodeScoreModels_.size(); c++)
	{
		if (RegionalPrmNodeScoreModels_.size() == 0 || (specificCharge>0 && specificCharge != c))
			continue;

		if (sa.getNumSpectraWithCharge(c)<200)
		{
			cout << "WARNING: insufficient number of spectra to train breakage model for charge " << c << endl;
			cout <<	"		  only " << sa.getNumSpectraWithCharge(c) << " spectra were found so this charge is being skipped!" << endl << endl;
			continue;
		}

		int s;
		for (s=0; s<RegionalPrmNodeScoreModels_[c].size(); s++)
		{
			if (specificSize>=0 && s != specificSize)
				continue;

			int r;
			for (r=0; r<RegionalPrmNodeScoreModels_[c][s].size(); r++)
			{
				if (specificRegion>=0 && r != specificRegion)
					continue;
				
				RegionalPrmNodeScoreModels_[c][s][r].trainRegionalScoreModel(allScoreModelsVoidPointer, name, sa);
			}
		}
	}

	// train PRM normalizer values
	cout << endl << "Training PRM normalizer vlaues..." << endl;

// TODO fix this issue, it needs to use the AllScoreModels class
//	learn_prm_normalizer_values(fm);

	ind_was_initialized=true;
}



void PrmNodeScoreModel::score_peptide_node_combos(void* model_ptr, 
												  PrmGraph *prm, const Peptide& peptide ) const
{

	const vector<int>& org_aas		= config_->get_org_aa();
	const vector<mass_t>& aa2mass	= config_->get_aa2mass();
	const vector<MultiEdge>& multi_edges = prm->get_multi_edges();
	const int num_nodes		   = prm->get_num_nodes();
	const vector<int>& pep_aas = peptide.get_amino_acids();
	const int num_pep_aas = pep_aas.size();

	mass_t p_mass=0;
	int aa_idx=0;

	int i;
	for (i=0; i<num_nodes; i++)
	{
		Node& node = prm->get_non_const_node(i);
		const RegionalPrmNodeScoreModel& score_model = 
			RegionalPrmNodeScoreModels_[prm->get_charge()][prm->get_size_idx()][node.breakage.region_idx];

		int in_edge_idx=NEG_INF,  in_edge_variant=NEG_INF;
		int out_edge_idx=NEG_INF, out_edge_variant=NEG_INF;

	//	cout << "N: " << node.mass << endl;
		while (aa_idx<pep_aas.size() && fabs(p_mass-node.mass)>0.1)
		{
			p_mass += aa2mass[pep_aas[aa_idx]];
			aa_idx++;
	//		cout << aa_idx << "\t" << p_mass << endl;
		}
		
		if (aa_idx == pep_aas.size() && i != num_nodes-1)
		{
			int j;
			for (j=0; j<num_nodes; j++)
				cout << j << "\t" << prm->get_node(j).mass << endl;

			cout << endl << "PEP:" << endl;
			vector<mass_t> exp_masses;
			peptide.calc_expected_breakage_masses(config_, exp_masses);
			for (j=0; j<exp_masses.size(); j++)
				cout << j << "\t" << exp_masses[j] << endl;

			cout << "Error: mismatch between nodes and peptide!" << endl;
			exit(1);
		}

		if (node.in_edge_idxs.size()>0)
		{
			int j;
			for (j=0; j<node.in_edge_idxs.size(); j++)
			{
				const int edge_idx = node.in_edge_idxs[j];
				const MultiEdge& in_edge = multi_edges[edge_idx];
				const int num_aa = in_edge.num_aa;

				if (num_aa>aa_idx)
					continue;

				const int var_idx = in_edge.get_variant_idx(num_aa,&pep_aas[aa_idx-num_aa]);
				if (var_idx<0)
					continue;

				in_edge_idx = edge_idx;
				in_edge_variant = var_idx;
				break;
			}
		}

		if (node.out_edge_idxs.size()>0)
		{
			int j;
			for (j=0; j<node.out_edge_idxs.size(); j++)
			{
				const int edge_idx = node.out_edge_idxs[j];
				const MultiEdge& out_edge = multi_edges[edge_idx];
				const int num_aa = out_edge.num_aa;

				if (num_aa + aa_idx >num_pep_aas)
					continue;

				const int var_idx = out_edge.get_variant_idx(num_aa,&pep_aas[aa_idx]);
				if (var_idx<0)
					continue;

				out_edge_idx = edge_idx;
				out_edge_variant = var_idx;
				break;
			}
		}

		const AllScoreModels* model = static_cast<AllScoreModels*>(model_ptr);
		BreakageInfo info;
		prm->fill_breakage_info(model,&info,i,in_edge_idx,in_edge_variant,out_edge_idx,out_edge_variant);
	
		node.score_combos.clear();

	//	cout << in_edge_idx << " " << in_edge_variant << " " << out_edge_idx << " " << out_edge_variant << "\t";

		info.score = score_model.score_a_single_breakage_combo(prm, node, &node.breakage, info);
	
		node.score_combos[ScoreComboLoc(info)]=info.score;
		node.score = info.score; 
		node.breakage.score = node.score;

	//	cout << node.score << endl;
	}
	prm->set_has_node_combo_scores(true);
}


/***************************************************************************
Scores the Gap-Gap combination for all nodes and sets it as the node score
****************************************************************************/
void PrmNodeScoreModel::initial_combos_score(void *model_ptr, PrmGraph *prm) const
{
	
	const vector<int>& org_aas = config_->get_org_aa();
	const vector<MultiEdge>& multi_edges = prm->get_multi_edges();
	const int num_nodes = prm->get_num_nodes();
	const mass_t mid_mass = prm->get_pm_with_19()*0.5;
	
	const AllScoreModels* model = static_cast<AllScoreModels*>(model_ptr);

	int i;
	for (i=0; i<num_nodes; i++)
	{
		Node& node = prm->get_non_const_node(i);
		const RegionalPrmNodeScoreModel& score_model = 
			RegionalPrmNodeScoreModels_[prm->get_charge()][prm->get_size_idx()][node.breakage.region_idx];

		vector<BreakageInfo> infos;

		BreakageInfo double_gap_info;
		prm->fill_breakage_info(model,&double_gap_info,i,NEG_INF,NEG_INF,NEG_INF,NEG_INF);
		infos.push_back(double_gap_info);

		// fill all combos for a digest node
		if (node.type == NODE_DIGEST)
		{
			int j;
			for (j=0; j<node.in_edge_idxs.size(); j++)
			{
				const int in_edge_idx = node.in_edge_idxs[j];
				const MultiEdge& in_edge = multi_edges[in_edge_idx];
				const int num_in_variants = in_edge.get_num_variants();

				int k;
				for (k=0; k<num_in_variants; k++)
				{
					int a;
					for (a=0; a<node.out_edge_idxs.size(); a++)
					{
						const int out_edge_idx = node.out_edge_idxs[a];
						const MultiEdge& out_edge = multi_edges[out_edge_idx];
						const int num_out_variants = out_edge.get_num_variants();

						int b;
						for (b=0; b<num_out_variants; b++)
						{
							BreakageInfo info;
							prm->fill_breakage_info(model,&info,i,in_edge_idx,k,out_edge_idx,b);
							if (info.connects_to_C_term || info.connects_to_N_term)
								infos.push_back(info);
						}
					}
				}
			}
		}
		else
		{
			int j;
			for (j=0; j<node.in_edge_idxs.size(); j++)
			{
				const int in_edge_idx = node.in_edge_idxs[j];
				const MultiEdge& in_edge = multi_edges[in_edge_idx];
				const int num_in_variants = in_edge.get_num_variants();

				int k;
				for (k=0; k<num_in_variants; k++)
				{
					BreakageInfo c_gap_info;
					prm->fill_breakage_info(model,&c_gap_info,i,in_edge_idx,k,NEG_INF,NEG_INF);
					infos.push_back(c_gap_info);
				}
			}
		
			for (j=0; j<node.out_edge_idxs.size(); j++)
			{
				const int out_edge_idx = node.out_edge_idxs[j];
				const MultiEdge& out_edge = multi_edges[out_edge_idx];
				const int num_out_variants = out_edge.get_num_variants();

				int k;
				for (k=0; k<num_out_variants; k++)
				{
					BreakageInfo n_gap_info;
					prm->fill_breakage_info(model,&n_gap_info,i,NEG_INF,NEG_INF,out_edge_idx,k);
					infos.push_back(n_gap_info);
				}
			}
		}

		node.score_combos.clear();
		int j;
		for (j=0; j<infos.size(); j++)
			infos[j].score = score_model.score_a_single_breakage_combo(prm, node, &node.breakage, infos[j]);

		for (j=0; j<infos.size(); j++)
			node.score_combos[ScoreComboLoc(infos[j])]=infos[j].score;
		
		score_t max_score=NEG_INF;
		for (j=1; j<infos.size(); j++)
			if (infos[j].score>max_score)
				max_score=infos[j].score;

		node.score = max_score; 
		node.breakage.score = node.score;
	}

	prm->set_has_node_combo_scores(true);
}


// performs scoring on demand (if the combo was not previously scored, calculates
// values, otherwise returns hashed value
score_t PrmNodeScoreModel::get_node_combo_score(PrmGraph *prm, int node_idx, 
										 int in_edge_idx,  int in_var_idx, 
										 int out_edge_idx, int out_var_idx) const
{
	const Node& node = prm->get_node(node_idx);
	ScoreComboLoc loc(in_edge_idx,out_edge_idx,in_var_idx,out_var_idx);

	const map<ScoreComboLoc,score_t>::const_iterator it = node.score_combos.find(loc);
	if (it != node.score_combos.end())
		return it->second;

	if (node.type == NODE_N_TERM || node.type == NODE_C_TERM)
	{
		const score_t terminal_score=config_->get_terminal_score();
		Node& non_const_node = prm->get_non_const_node(node_idx);
		non_const_node.score_combos[loc]=terminal_score;
		return terminal_score;
	}

	const RegionalPrmNodeScoreModel& score_model = 
			RegionalPrmNodeScoreModels_[prm->get_charge()][prm->get_size_idx()][node.breakage.region_idx];
	Node& non_const_node = prm->get_non_const_node(node_idx);

	BreakageInfo info;
	prm->fill_breakage_info(prm->get_model(),&info,node_idx,in_edge_idx,in_var_idx,out_edge_idx,out_var_idx);

	// score not here, calculate the combo score, store it, and return it
	const int charge = prm->get_charge();
	const int size_idx = prm->get_size_idx();
	Spectrum *spec = prm->get_source_spectrum();
	mass_t pm_with_19 = prm->get_pm_with_19();
	const RegionalPrmNodeScoreModel& rsm = RegionalPrmNodeScoreModels_[charge][size_idx][info.breakage->region_idx];

	score_t combo_score = rsm.score_a_single_breakage_combo(prm, non_const_node, info.breakage,info);

	non_const_node.score_combos[loc]=combo_score;
	non_const_node.breakage.score = combo_score;

	return combo_score;

}

