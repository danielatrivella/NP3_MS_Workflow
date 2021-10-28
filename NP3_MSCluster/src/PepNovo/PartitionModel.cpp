#include "PeakRankModel.h"





/*****************************************************************************
Creates the phi domain for the rank boost training.
phi contains pairs of idx0,idx1 and weight. 
idx0 and idx1 are peak sample indexes while weight is the relative intensity
between them. The weight is normalized according to the maximal intensity of a
peak from the same frag type in the peptide (e.g., the strongest y is given 
intensity 1 while a missing y is given intensity 0. The weight of x,y is set to 
intensity[y]-intensity[x] (it is assumed that for all pairs, y is stronger than x).
******************************************************************************/
void create_phi_list_from_samples(const vector<float>&     peak_intens,
								  const vector<PeakStart>& peak_starts,
								  const vector<float>&     max_annotated_intens,
								  vector<SamplePairWeight>& phi)
{
	const weight_t min_norm_weight = 0.1;

	int i;
	phi.clear();
	for (i=0; i<peak_starts.size(); i++)
	{
		const int start_idx = peak_starts[i].peak_start_idx;
		const int num_peaks = peak_starts[i].num_peaks;
	
		vector<float> intens;

		intens.resize(num_peaks,-1);
		float max_inten = -1;

		int j;
		for (j=0; j<num_peaks; j++)
		{
			intens[j]=peak_intens[start_idx+j];
			if (intens[j]>max_inten)
				max_inten=intens[j];
		}

		if (max_inten<=0)
		{
//			cout << "Error: could not find a max inten for peaks from sample " << i << endl;
//			exit(1);
			continue;
		}

		const float one_over_max = 1.0/max_inten;
		for (j=0; j<num_peaks; j++)
			intens[j] *= one_over_max;

		const float max_annotated_inten = max_annotated_intens[i]*one_over_max;

		// look at all pairs, only consider pairs for which inten[idx1]>inten[idx0]
		for (j=0; j<num_peaks; j++)
		{
			int k;
			for (k=0; k<num_peaks; k++)
				if (intens[k]>intens[j])
				{
					// ignore a pair of peaks if the difference between them is less than
					// 0.15 of the maximum intensity (there might be too many misordered 
					// pairs in that range
					if ( (intens[j]/intens[k]) > 0.85)
						continue;

					// calc weight of sample
					// hlaf is relative weight of strong peak (compared to all annotated peaks
					// in the spectrum). the other half is relative to the difference between
					// the two peaks being compared).

				
					weight_t pair_weight = 0.6*(intens[k]-intens[j]) + 0.4 *intens[k]/max_annotated_inten;
					if (pair_weight>1.0)
						pair_weight = 1.0;

					if (pair_weight<min_norm_weight)
						continue;

					// don't let very weak peaks participate in a missing peak pair
					if (intens[j]==0 && pair_weight<0.25)
						continue;

					phi.push_back(SamplePairWeight(start_idx+j,start_idx+k,pair_weight));
				}
		}
	}

	cout << "Created " << phi.size() << " pairs for " << peak_starts.size() << " peptide samples" << endl;
}


/******************************************************************************
The main function at the partition model level for calculating the peaks' rank
scores. Returns results in the PeptidePeakPrediction structure, in there there
is a 2D table of score (rows are the fragment type, columns are the cut idxs).
number of rows equals the number of fragment idxs in the ppp
*******************************************************************************/
void PartitionModel::calc_peptides_peaks_rank_scores(
									  const PeakRankModel *prank,
									  const PeptideSolution& sol,
									  mass_t min_detected_mass, 
									  mass_t max_detected_mass,
									  PeptidePeakPrediction& ppp,
									  int   feature_set_type,
									  const vector<int>* ptr_frag_type_idxs) const
{
	Config *config = prank->get_config();
	const int num_frags =		 fragment_type_idxs.size();
	const mass_t mass_with_19  = sol.pm_with_19;
	const Peptide& pep = sol.pep;
	const vector<int>& amino_acids = pep.get_amino_acids();
	vector<mass_t> exp_cuts;

	pep.calc_expected_breakage_masses(config,exp_cuts);

	ppp.amino_acids = pep.get_amino_acids();
	ppp.num_frags = num_frags;
	ppp.frag_idxs = fragment_type_idxs;
	ppp.rank_scores.resize(fragment_type_idxs.size());

	int f;
	const mass_t n_mass = pep.get_n_gap();

	// calculate a single set of ranks across the combined set of fragments
	if (this->feature_set_type >= 3)
	{

		const int start_cut_idx = (sol.reaches_n_terminal ? 1 : 0);
		const int last_cut_idx  = (sol.reaches_c_terminal ? exp_cuts.size()-1 : exp_cuts.size());
		const mass_t n_mass = exp_cuts[0];
		const mass_t c_mass = exp_cuts[exp_cuts.size()-1];

		ppp.combined_peak_scores.clear();
		ppp.frag_idxs = fragment_type_idxs;
		for (f=0; f<num_frags; f++)
		{
			const int frag_idx = fragment_type_idxs[f];
			const FragmentType& fragment = config->get_fragment(frag_idx);

			ppp.rank_scores[f].resize(exp_cuts.size(),NEG_INF);
			
			int cut_idx;
			for (cut_idx=start_cut_idx; cut_idx<last_cut_idx; cut_idx++)
			{
				const mass_t cut_mass = exp_cuts[cut_idx];
				const mass_t peak_mass = fragment.calc_expected_mass(cut_mass,mass_with_19);
				if (peak_mass<min_detected_mass || peak_mass>max_detected_mass)
					continue;

				
				RankBoostSample rbs;

				if (feature_set_type == 3)
				{
					fill_combined_peak_features(prank, amino_acids, cut_idx, cut_mass, sol, fragment, f, rbs);
				}
				else if (feature_set_type == 4)
				{
					fill_combined_dnv_peak_features(prank, n_mass, c_mass, amino_acids, cut_idx, 
						cut_mass, sol, fragment, f, rbs);
				}
				else if (feature_set_type == 5)
				{
					fill_combined_simple_peak_features(prank, amino_acids, cut_idx, cut_mass, sol, fragment, f, rbs);
				}
				else
				{
					cout << "Error: feature set type not supported: " << feature_set_type << endl;
					exit(1);
				}
				ppp.rank_scores[f][cut_idx]= this->combined_frag_boost_model.calc_rank_score(rbs);
			}
		}
		return;
	}


	// calculate separate sets of ranks for each fragment type
	if (this->feature_set_type<=2)
	{
		for (f=0; f<num_frags; f++)
		{
			const int frag_idx = fragment_type_idxs[f];

			// only compute results for frags in the list
			if (ptr_frag_type_idxs)
			{
				int i;
				for (i=0; i<(*ptr_frag_type_idxs).size(); i++)
					if ( (*ptr_frag_type_idxs)[i] == frag_idx)
						break;
				if (i==(*ptr_frag_type_idxs).size())
					continue;
			}

			ppp.rank_scores[frag_idx].resize(exp_cuts.size(),NEG_INF);

			const FragmentType& fragment = config->get_fragment(frag_idx);
			int cut_idx;
			for (cut_idx=1; cut_idx<exp_cuts.size()-1; cut_idx++)
			{
				const mass_t cut_mass = exp_cuts[cut_idx];
				const mass_t peak_mass = fragment.calc_expected_mass(cut_mass,mass_with_19);
				if (peak_mass<min_detected_mass || peak_mass>max_detected_mass)
					continue;

				RankBoostSample rbs;
		
				if (feature_set_type == 0)
				{
					prank->fill_simple_peak_features(amino_acids, cut_idx, 
						cut_mass, mass_with_19, charge, fragment, rbs);
				}
				else if (feature_set_type == 1)
				{
					prank->fill_advanced_peak_features(amino_acids, cut_idx, 
						cut_mass, mass_with_19, charge, fragment, rbs);
				}
				else if (feature_set_type == 2)
				{
					prank->fill_partial_denovo_peak_features(pep.get_n_gap(), exp_cuts[exp_cuts.size()-1],
						amino_acids,cut_idx, cut_mass, 
						mass_with_19, charge, fragment, ppp.most_basic_missing_on_n, 
						ppp.most_basic_missing_on_c , rbs);
				}
				else
				{
					cout << "Error: feature set type not supported: " << feature_set_type << endl;
					exit(1);
				}

				ppp.rank_scores[frag_idx][cut_idx]=frag_models[f].calc_rank_score(rbs);
			}
		}
	}
}


int PartitionModel::read_partition_model(const string& path, Config *config,
										 int _charge, int _size_idx, int _mobility)
{
	const int max_global_frag_idx = config->get_all_fragments().size();

	charge = _charge;
	size_idx = _size_idx;
	mobility   = _mobility;
	max_frag_idx = -1;

	fragment_type_idxs.clear();
	frag_models.clear();

	if (this->feature_set_type <=2)
	{

		int f;
		for (f=0; f<max_global_frag_idx; f++)
		{
			ostringstream oss;
			oss << path << "_" << charge << "_" << this->size_idx << "_" << mobility <<
				"_" << f << ".txt";

	//		cout << "checking: " << oss.str();
			
			ifstream boost_file(oss.str().c_str());

			if (! boost_file.is_open()  ||! boost_file.good())
			{
		//		cout << " -" << endl;
				continue;
			}

			
			fragment_type_idxs.push_back(f);

			max_frag_idx=f;
			
			RankBoostModel rbm;

			frag_models.push_back(rbm);
			frag_models[frag_models.size()-1].read_rankboost_model(boost_file);
			boost_file.close();
		//	cout << " +" << endl;
		}
		return fragment_type_idxs.size();
	}
	else
	{
		if (read_combined_partition_model(path,config,_charge,_size_idx,_mobility)>0)
			return 1;
		return 0;
	}
	
}


/*********************************************************************

**********************************************************************/
int PartitionModel::write_partition_model(const string& path)
{
	int f;
	int num_written=0;
	if (this->feature_set_type <= 2)
	{
		for (f=0; f<this->fragment_type_idxs.size(); f++)
		{
			const int frag_idx = fragment_type_idxs[f];

			if (! frag_models[f].get_ind_was_initialized())
				continue;

			ostringstream oss;
			oss << path << "_" << charge << "_" << this->size_idx << "_" << mobility <<
				"_" << frag_idx << ".txt";

			cout << "Writing: " << oss.str() << endl;
			
			fstream boost_file(oss.str().c_str(),ios::out);

			if (! boost_file.is_open() || ! boost_file.good())
			{
				cout << "Error: couldn't open file for wiriting!" << endl;
				exit(1);
			}
			
			frag_models[f].write_rankboost_model(boost_file,true);
			boost_file.close();
			num_written++;
		}
		return num_written;
	}
	else
		return write_combined_partition_model(path);

	return num_written;
}




/*****************************************************************************
Takes as an input a training file of peptide samples.
For each fragment type with sufficient count, it creates a dataset and
trains the rank model.
******************************************************************************/
void PartitionModel::train_partition_model(
									   PeakRankModel *prank, 
									   char *sample_file_path,
									   int	_charge,
									   int  _size_idx,
									   int  _mobility,
									   int frag_idx_to_train,
									   char *report_dir,
									   int max_num_rounds,
									   char *test_set_file,
									   int   test_peptide_length,
									   char *stop_signal_file,
									   weight_t max_weight_ratio)
{
	const float min_ratio_detected = 0.15;
	Config *config = prank->get_config();
	const vector<FragmentType>& all_fragments = config->get_all_fragments();
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const int num_frags = all_fragments.size();

	if (max_num_rounds<0)
		max_num_rounds = 1000;

	charge = _charge;
	size_idx = _size_idx;
	mobility   = _mobility;
	
	vector<TrainingPeptide> sample_tps, test_tps;
	vector<int> frag_counts;
	vector<int> frag_detected;
	vector<int> length_counts;

	read_training_peptides_from_file(sample_file_path,sample_tps);

	cout << "Read " << sample_tps.size() << " training tps...";

	int num_tps_to_add = 0;

	if (prank->get_feature_set_type() == 2)
	{
		if (sample_tps.size()<25000)
			num_tps_to_add = 1;
		if (sample_tps.size()<15000)
			num_tps_to_add = 2;
		if (sample_tps.size()<10000)
			num_tps_to_add = 3;
		if (sample_tps.size()<5000)
			num_tps_to_add = 4;

		cout << "Adding at most " << num_tps_to_add << " per tp." << endl;
	}

	if (prank->get_feature_set_type() == 2)
		convert_tps_to_partial_denovo(config,sample_tps,num_tps_to_add);

	if (test_set_file)
	{
		read_training_peptides_from_file(test_set_file, test_tps);
		cout << "Read " << test_tps.size() << " test tps...";
		if (prank->get_feature_set_type() == 2)
			convert_tps_to_partial_denovo(config,test_tps,num_tps_to_add);
	}

	// Create intial report on dataset
	frag_counts.resize(num_frags,0);
	frag_detected.resize(num_frags,0);
	length_counts.resize(200,0);

	int numH=0,numK=0, numR=0;
	int i;
	for (i=0; i<sample_tps.size(); i++)
	{
		const TrainingPeptide& tp = sample_tps[i];
		int f;

		for (f=0; f<tp.intens.size(); f++)
		{
			int j;
			for (j=1; j<tp.intens[f].size(); j++)
			{
				if (tp.intens[f][j]>=0)
					frag_counts[tp.frag_idxs[f]]++;
				if (tp.intens[f][j]>0)
					frag_detected[tp.frag_idxs[f]]++;
			}
		}

		int j;
		for (j=0; j<tp.amino_acids.size(); j++)
		{
			if (tp.amino_acids[j]==His)
				numH++;
			if (tp.amino_acids[j]==Lys)
				numK++;
			if (tp.amino_acids[j]==Arg)
				numR++;
		}
		length_counts[tp.amino_acids.size()]++;
	}

	// report and select fragments for training
	cout << "# training peptides: " << sample_tps.size() << endl;
	cout << "Avg #R: " << numR/(double)sample_tps.size() << endl;
	cout << "Avg #K: " << numK/(double)sample_tps.size() << endl;
	cout << "Avg #H: " << numH/(double)sample_tps.size() << endl;

	cout << endl << "Sample lengths:" << endl;
	for (i=0; i<length_counts.size(); i++)
		if (length_counts[i]>0)
			cout << i << "\t" << length_counts[i] << endl;
	cout << endl;
		
	for (i=0; i<all_fragments.size(); i++)
	{
		float ratio = (float)frag_detected[i]/frag_counts[i];
		cout << all_fragments[i].label << "\t" << frag_detected[i] << " / " << 
			frag_counts[i] << "\t = " << setprecision(3) << ratio << endl;

		if (ratio>=min_ratio_detected)
		{
			if (frag_idx_to_train<0 || frag_idx_to_train == i)
				fragment_type_idxs.push_back(i);
		}
	}
	cout << endl;

	cout << "Max weight ratio: " << max_weight_ratio << endl;

	if (fragment_type_idxs.size() == 0)
	{
		cout << "No models to train!" << endl;
		return;
	}


	// Train each selected model

	frag_models.resize(fragment_type_idxs.size());
	int f;
	for (f=0; f<fragment_type_idxs.size(); f++)
	{
		const int frag_idx = fragment_type_idxs[f];

		if (frag_idx>0 && frag_idx != frag_idx_to_train)
			continue;

		const int frag_charge = all_fragments[frag_idx].charge;

		cout << "Training frag " << frag_idx << " (" << config->get_fragment(frag_idx).label <<")" << endl;
		
		// fill RankBoostSamples and create rank ds
		RankBoostDataset         rank_ds,test_ds;
		
		vector<float>			 peak_intens;
		vector<PeakStart>        peak_starts;
		vector<float>			 max_annotated_intens;

		// Train RankBoostModel

		cout << "TRAINING..." << endl;

		// initialize and read the test set if it exists
		RankBoostDataset *test_set_ptr=NULL;
		if (test_set_file)
		{
			vector<float>			 test_peak_intens;
			vector<PeakStart>        test_peak_starts;
			vector<float>			 test_max_annotated_intens;

			cout << "Reading test tps..." << endl; 

			prank->read_training_peptides_into_rank_boost_dataset(frag_idx, charge, 
				test_tps, test_ds, test_peak_intens, test_peak_starts, 
				test_max_annotated_intens);

			cout << "Creating test phi list..." << endl;

			create_phi_list_from_samples(test_peak_intens, test_peak_starts, 
				test_max_annotated_intens, test_ds.get_non_const_phi_support());

			test_ds.compute_total_phi_weight();
			test_set_ptr = &test_ds;

			// choose length (try shorte peptide if not eonough samples, go for the max)
			if (test_peptide_length == 0)
			{
				vector<int> test_length_counts;
				test_length_counts.resize(200,0);
				const vector<RankBoostSample>& samples = test_ds.get_samples();
				vector<int> sizes;
				sizes.resize(test_ds.get_num_groups(),0);
				int i;
				for (i=0; i<samples.size(); i++)
					sizes[samples[i].groupIndex]=samples[i].tag3;

				for (i=1; i<sizes.size(); i++)
					test_length_counts[sizes[i]]++;
				
				int max=0;
				for (i=0; i<200; i++)
				{
					if (test_length_counts[i]>=500)
						break;
					if (test_length_counts[i]>test_length_counts[max])
						max=i;
				}

				if (i<200)
				{
					test_peptide_length = i;
				}
				else
					test_peptide_length = max;
			}
			cout << "test length " << test_peptide_length << endl;
		}
	
		
		cout << "read training tps..." << endl;

		prank->read_training_peptides_into_rank_boost_dataset(frag_idx, charge,  
			sample_tps, rank_ds, peak_intens, peak_starts, max_annotated_intens);

		RankBoostModel& boost = frag_models[f];
		boost.init_rankboost_model_feature_names(prank->get_binary_names(),prank->get_real_names());

		cout << "create training phi list..." << endl;

		create_phi_list_from_samples(peak_intens,peak_starts, max_annotated_intens, 
			rank_ds.get_non_const_phi_support());

		cout << "initializing boost..." << endl;

		rank_ds.compute_total_phi_weight();
		rank_ds.initialize_potenital_lists();
		rank_ds.initialize_binary_one_lists(prank->get_binary_names().size());
		rank_ds.initialize_binary_ordered_phi_lists(boost.get_ptr_to_binary_feature_names());
		rank_ds.initialzie_real_feature_table(prank->get_real_names().size());

		rank_ds.set_max_ratio_for_regular_update(max_weight_ratio);
				
		boost.init_rankboost_model_for_training(rank_ds,40,100);
		
		rank_ds.initialize_real_vote_lists(boost);

	//	boost.summarize_features(rank_ds.get_samples());

		char report_prefix[512];
		if (report_dir)
			sprintf(report_prefix,"%s/%s_%d",report_dir, partition_name.c_str(),frag_idx);

		vector<idx_weight_pair> miss_pairs;

		boost.train_rankboost_model(rank_ds,
									max_num_rounds, 
									&miss_pairs, 
									test_set_ptr, 
									test_peptide_length, 
									report_prefix, 
									stop_signal_file);

		

		// final report
		if (report_dir)
		{
			char name_buff[512];

		
			sprintf(name_buff,"%s_train_miss_pairs.txt",report_prefix);
			ofstream report_stream(name_buff);
			if (! report_stream.is_open() || ! report_stream.good())
			{
				cout << "Error: couldn't open pairs report file for writing:" << name_buff << endl;
				exit(1);
			}

			simple_print_peak_pairs(miss_pairs, sample_tps, rank_ds, prank, frag_idx, 250, report_stream);
			report_stream.close();


			sprintf(name_buff,"%s_feature_list.txt",report_prefix);
			ofstream feature_stream(name_buff);
			if (! feature_stream.is_open() || ! feature_stream.good())
			{
				cout << "Error: couldn't feature_stream file for writing:" << name_buff << endl;
				exit(1);
			}
			cout << "[...";
			boost.ouput_importance_ranked_feature_list(rank_ds,feature_stream);
			cout << " ...]" << endl;
			feature_stream.close();


			// write model (also compresses features and deletes the default values)
			sprintf(name_buff,"%s_model.txt",report_prefix);
			ofstream model_stream(name_buff);
			boost.write_rankboost_model(model_stream,true);
			model_stream.close();
		}	
		else // send to cout
			simple_print_peak_pairs(miss_pairs, sample_tps, rank_ds, prank, frag_idx, 100);
	}
}



void  PartitionModel::set_partition_name(const string& peak_rank_model_name, 
									     int charge, int size_idx, int mobility)
{
	char buff[256];
	
	sprintf(buff,"%s_%d_%d_%d",peak_rank_model_name.c_str(),charge,size_idx,mobility);

	partition_name  = string(buff);
}



/*************************************************************************
**************************************************************************/
void PartitionModel::simple_print_peak_pairs(
						  const vector<idx_weight_pair>& pair_idxs, 
						  const vector<TrainingPeptide>& tps,
						  const RankBoostDataset& ds,
						  PeakRankModel *prank,
						  int frag,
						  int max_examples,
						  ostream& os) const
{
	const int max_examples_for_tp = 2;
	const vector<SamplePairWeight>& phi_support= ds.get_phi_support();
	vector<int> tp_idx_counts;
	tp_idx_counts.clear();

	if (this->feature_set_type > 2)
	{
		cout << "Error: this funcion is not designed for feature set " << feature_set_type << endl;
		exit(1);
	}

	int f;
	for (f=0; f<this->fragment_type_idxs.size(); f++)
		if (fragment_type_idxs[f]==frag)
			break;
	if (f== fragment_type_idxs.size())
	{
		cout << "Error: partition does not have an initialized model for frag " << frag << endl;
		exit(1);
	}
	const int frag_pos=f;

	int p_idx;
	int counter = 0;
	for (p_idx=0; p_idx<pair_idxs.size(); p_idx++)
	{
		const int pair_idx = pair_idxs[p_idx].idx;
		const int x0_idx = phi_support[pair_idx].idx1;
		const int x1_idx = phi_support[pair_idx].idx2;
		const RankBoostSample sam_x0 = ds.get_sample(x0_idx);
		const RankBoostSample sam_x1 = ds.get_sample(x1_idx);
		const int tp_idx = sam_x0.groupIndex;
		const int x0_cut_idx = sam_x0.tag2;
		const int x1_cut_idx = sam_x1.tag2;
		const TrainingPeptide& tp = tps[tp_idx];
		const RankBoostModel& boost_model = frag_models[frag_pos];

		const int num_binary_features = boost_model.get_num_binary_features();
		const int num_real_features   = boost_model.get_num_real_features();
		const vector<string>& binary_feature_names = boost_model.get_binary_feature_names();
		const vector<string>& real_feature_names   = boost_model.get_real_feature_names();
		const vector<string>& model_aa_labels = prank->get_model_aa_labels();

		const vector<int>& bin_updates  = boost_model.get_binary_update_counts();
		const vector<int>& real_updates = boost_model.get_real_update_counts();

		int j;
		for (j=0; j<tp_idx_counts.size(); j++)
			if (tp_idx_counts[j]==tp_idx)
				break;

		if (j==tp_idx_counts.size())
		{
			tp_idx_counts.push_back(tp_idx);
		}
		else
		{
			if (tp_idx_counts[j]>=max_examples_for_tp)
				continue;
			tp_idx_counts[j]++;
		}

		counter++;

		vector<int> amino_acids;
		prank->convert_aas_to_model_aas(tp.amino_acids,amino_acids);

		if (sam_x0.groupIndex <0 || sam_x0.tag2<0 || sam_x1.groupIndex <0 || sam_x1.tag2<0)
		{
			cout << "Error: bad tags attached to samples:" << endl;
			cout << sam_x0.groupIndex << " " << sam_x0.tag1 << " " << sam_x0.tag2 << endl;
			cout << sam_x1.groupIndex << " " << sam_x1.tag1 << " " << sam_x1.tag2 << endl;
			exit(1);
		}
		
		if (sam_x0.groupIndex != sam_x1.groupIndex)
		{
			cout << "Error: groupIndex mismathces!" << endl;
			exit(1);
		}

		const int pos = tp.get_frag_idx_pos(frag);
		if (pos<0)
		{
			cout << "Error: frag not in tp intens!" << endl;
			exit(1);
		}
		const vector<float>& intens = tp.intens[pos];

	
		os << counter << "\t" << phi_support[pair_idx].weight << "\t" << pair_idxs[p_idx].weight << "\t";

		int a;
		for (a=0; a<amino_acids.size(); a++)
		{
			if (a == x0_cut_idx)
				os << " 0 ";

			if (a == x1_cut_idx)
				os << " 1 ";

			os << model_aa_labels[amino_acids[a]];
		}
		os << endl;



		os << "x0: cut " << x0_cut_idx << "  inten " << intens[x0_cut_idx] << "\t score " <<
			boost_model.calc_rank_score(sam_x0) << endl;

		os << "x1: cut " << x1_cut_idx << "  inten " << intens[x1_cut_idx] << "\t score " <<
			boost_model.calc_rank_score(sam_x1) << endl;

		os << endl << endl;

		if (max_examples>0 && counter>=max_examples)
			break;
	}
}


void   PartitionModel::print_model_stats() const
{
	if (this->feature_set_type <=2)
	{
		cout << this->partition_name << " \t";
		int i;
		for (i=0; i<10; i++)
		{
			int j;
			for (j=0; j<this->fragment_type_idxs.size(); j++)
				if (fragment_type_idxs[j]==i && frag_models[j].get_ind_was_initialized())
				{
					cout << " *";
					break;
				}
			if (j==fragment_type_idxs.size())
				cout << "  ";
		}
		cout << endl;
	}
	else
	{
		cout << this->partition_name << " \t";
		int i;
		for (i=0; i<10; i++)
		{
			int j;
			for (j=0; j<this->fragment_type_idxs.size(); j++)
				if (fragment_type_idxs[j]==i && this->combined_frag_boost_model.get_ind_was_initialized())
				{
					cout << " *";
					break;
				}
			if (j==fragment_type_idxs.size())
				cout << "  ";
		}
		cout << endl;
	}
}

/***********************************************************************
makes tables listing features and final scores
Only makes table if the predictions match
************************************************************************/
bool PeakRankModel::make_peak_prediction_table(
			const PeptideSolution& sol,
			const vector< vector<intensity_t> >& intens,
			int num_peaks) const
{
	PeptidePeakPrediction ppp;
	calc_peptide_predicted_scores(sol, ppp);

	// the ppp includes a table of rank scores (rows are actual frag idxs, not relative
	// position in the frag_type_idxs).

	// reduce intensities to the same dimensionality
	const int num_frags = ppp.frag_idxs.size();
	vector< vector< float> > observed_intens;
	observed_intens.resize(num_frags);

	int i,f;
	for (f=0; f<num_frags; f++)
	{
		const int frag_idx = ppp.frag_idxs[f];
		observed_intens[f]=intens[frag_idx]; 
	}

	// calculate the ranks and mapping between predicted and observed
	vector< vector<int> > observed_ranks, predicted_ranks;
	calc_combined_peak_ranks(observed_intens, observed_ranks);
	calc_combined_peak_ranks(ppp.rank_scores, predicted_ranks);

	vector<int> sel_frags, sel_idxs;
	vector< float > intensities;
	
	
	int rank;
	for (rank=0; rank<num_peaks; rank++)
	{
		bool good_pred=false;
		for (f=0; f<num_frags; f++)
		{
			int i;
			for (i=0; i<predicted_ranks[f].size(); i++)
			{
				if (predicted_ranks[f][i] == rank &&
					observed_ranks[f][i]  == rank)
				{
					good_pred=true;
					sel_frags.push_back(f);
					sel_idxs.push_back(i);
					intensities.push_back(intens[f][i]);
					break;
				}
			}
		}
		if (! good_pred)
			return false;
	}

//	cout << "#sel_frags: " << sel_frags.size() << endl;
	

	// calc specific peak vectors and collect data
	vector< vector< string> > feature_names;
	vector< vector< float > > feature_values;
	vector< vector< float > > feature_scores;
	vector< float > total_scores;


	feature_names.resize(num_peaks);
	feature_values.resize(num_peaks);
	feature_scores.resize(num_peaks);
	total_scores.resize(num_peaks,0);
	

	const Peptide& pep = sol.pep;
	const mass_t pm_with_19 = sol.pm_with_19;
	const int spec_charge = sol.charge;
	const int mobility = get_proton_mobility(pep,spec_charge);
	const int size_idx =  get_size_group(spec_charge,pm_with_19);
	
	if (! partition_models[spec_charge][size_idx][mobility])
	{
		cout << "Error: no rank partition model for " <<
			spec_charge << " " << size_idx << " " << mobility << endl;
		exit(1);
	}

	if (mobility != 2)
		return false;

	const mass_t min_detected_mass = calc_min_detected_mass(pm_with_19, spec_charge);
	const mass_t max_detected_mass = get_max_detected_mass();


	const vector<int>& amino_acids = pep.get_amino_acids();
	vector<mass_t> exp_cuts;

	pep.calc_expected_breakage_masses(config,exp_cuts);

	const mass_t n_mass = pep.get_n_gap();

	// calculate a single set of ranks across the combined set of fragments
	const int start_cut_idx = (sol.reaches_n_terminal ? 1 : 0);
	const int last_cut_idx  = (sol.reaches_c_terminal ? exp_cuts.size()-1 : exp_cuts.size());
	const mass_t c_mass = exp_cuts[exp_cuts.size()-1];


	int max_l=0;
	for (i=0; i<sel_frags.size(); i++)
	{
		const int frag_idx=sel_frags[i];
		const int cut_idx = sel_idxs[i];
	
		const FragmentType& fragment = config->get_fragment(frag_idx);

		const mass_t cut_mass = exp_cuts[cut_idx];
		const mass_t peak_mass = fragment.calc_expected_mass(cut_mass,pm_with_19);
		
		RankBoostSample rbs;

		for (f=0; f<num_frags; f++)
			if (ppp.frag_idxs[f] == frag_idx)
				break;

	//	cout << "Frag: " << fragment.label << " fi:" << frag_idx << " f:" << f << endl;

		if (f==num_frags)
		{
			cout << "Error: bad frag!!!!" << endl;
			exit(1);
		}
		

		partition_models[spec_charge][size_idx][mobility]->fill_combined_simple_peak_features(
			this, amino_acids, cut_idx, cut_mass, sol, fragment, f, rbs);
				
//		partition_models[spec_charge][size_idx][mobility]->fill_combined_peak_features(	
//			this, amino_acids, cut_idx, cut_mass, sol, fragment, f, rbs);
			
		total_scores[i] = partition_models[spec_charge][size_idx][mobility]->combined_frag_boost_model.calc_rank_score_with_details(
									rbs,feature_names[i],feature_values[i],feature_scores[i]);
							
			
		if (feature_names[i].size()>max_l)
			max_l = feature_names[i].size();
	}


	cout << "Size: " << size_idx << " Mobility: " << mobility << endl;


	// print results
	for (i=0; i<num_peaks; i++)
	{
		cout << config->get_fragment(sel_frags[i]).label << " " <<
			sel_idxs[i];
		
		if (i<num_peaks-1)
		{
			cout << " & ";
		}
		else
			cout << "\\\\" << endl;
	}

	cout << setprecision(2) << fixed;
	for (i=0; i<num_peaks; i++)
	{
		cout << total_scores[i];
		if (i<num_peaks-1)
		{
			cout << " & ";
		}
		else
			cout << "\\\\" << endl;
	}

	for (i=0; i<num_peaks; i++)
	{
		cout << intensities[i];
		if (i<num_peaks-1)
		{
			cout << " & ";
		}
		else
			cout << "\\\\" << endl;
	}

	for (i=0; i<max_l; i++)
	{
		int j;
		for (j=0; j<num_peaks; j++)
		{
			if (feature_names[j].size()<=i)
			{
				cout << "           &  ";  
			}
			else
			{
				cout << feature_names[j][i] << " " << feature_values[j][i] << " & ";
				if (feature_scores[j][i]>0)
				{
					cout << "+";
				}
				cout << feature_scores[j][i];
			}

			if (j<num_peaks-1)
			{
				cout << " & ";
			}
			else
				cout << "\\\\" << endl;
		}
	}



	return true;
}


