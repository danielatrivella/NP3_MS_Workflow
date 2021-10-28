#include "PeakRankModel.h"
#include "PepNovo_auxfun.h"

extern const int num_RKH_combos;
extern const int num_RKH_pairs;
extern const string combos[];
extern const float RKH_pair_matrix[6][6];


/************************************************************************

  A note about the different structure of functions/classes when the 
  combined features are used. Since each partition has different fragments,
  which induce different fragment names, most of the model informaiton was
  pushed into the partition level, rather than the PeakRankModel level, with
  the other models that looked at each fragment independently.
*************************************************************************/




int PartitionModel::read_combined_partition_model(const string& path, Config *config,
							 int _charge, int _size_idx, int _mobility)
{
	const int max_global_frag_idx = config->get_all_fragments().size();

	charge = _charge;
	size_idx = _size_idx;
	mobility   = _mobility;
	max_frag_idx = -1;

	fragment_type_idxs.clear();
		
	ostringstream oss;
	oss << path << "_" << charge << "_" << this->size_idx << "_" << mobility <<	"_model.txt";
	
	ifstream boost_file(oss.str().c_str());

	if (! boost_file.is_open()  ||! boost_file.good())
	{
		cout << "Warnind: Couldn\'t find " << oss.str() << endl;
		return 0;
	}

	// read feature type and frags
	char buff[256];
	boost_file.getline(buff,256);
	this->feature_set_type = NEG_INF;
	int num_frags = NEG_INF;
	istringstream iss(buff);

	iss >> feature_set_type >> num_frags;
	if (this->feature_set_type <=2)
	{
		cout << "Error: read_combined_partition_model works only on frag_set 3 and up!" << endl;
		exit(1);
	}

	int f;
	for (f=0; f<num_frags; f++)
	{
		int frag=NEG_INF;
		iss >> frag;
		if (frag<0)
		{
			cout << "Error: bad frags in combined model!" << endl;
			exit(1);
		}
		fragment_type_idxs.push_back(frag);
		if (frag>max_frag_idx)
			max_frag_idx = frag;
	}
			
	combined_frag_boost_model.read_rankboost_model(boost_file);
	boost_file.close();

	num_features_per_frag= (combined_frag_boost_model.get_num_real_features()-1) / 
						   fragment_type_idxs.size();
	
	return fragment_type_idxs.size();
}
	

int PartitionModel::write_combined_partition_model(const string& path)
{
	if (this->feature_set_type <=2)
	{
		cout << "Error: read_combined_partition_model works only on frag_set 3 and up!" << endl;
		exit(1);
	}

	ostringstream oss;
	fstream boost_file(path.c_str(),ios::out);

	if (! boost_file.is_open() || ! boost_file.good())
	{
		cout << "Error: couldn't open file for wiriting!" << endl;
		exit(1);
	}

	
	boost_file << feature_set_type << " " << fragment_type_idxs.size();
	int f;
	for (f=0; f<fragment_type_idxs.size(); f++)
		boost_file << " " << fragment_type_idxs[f];
	boost_file << endl;
		
	combined_frag_boost_model.write_rankboost_model(boost_file,true);

	return 1;
}


/*****************************************************************************
Creates the phi domain for the rank boost training.
phi contains pairs of idx0,idx1 and weight. 
idx0 and idx1 are peak sample indexes while weight is the relative intensity
between them. The weight is normalized according to the maximal intensity of a
peak from the same frag type in the peptide (e.g., the strongest y is given 
intensity 1 while a missing y is given intensity 0. The weight of x,y is set to 
intensity[y]-intensity[x] (it is assumed that for all pairs, y is stronger than x).
******************************************************************************/
void create_phi_list_from_combined_peak_samples(const vector<float>& peak_intens,
												const vector<int>  & peak_frag_types,    
											    const vector<PeakStart>& peak_starts,
												const int max_frag_idx,
												vector<SamplePairWeight>& phi)
{
	const weight_t min_norm_weight = 0.2;
	int i;
	phi.clear();
	for (i=0; i<peak_starts.size(); i++)
	{
		const int start_idx = peak_starts[i].peak_start_idx;
		const int num_peaks = peak_starts[i].num_peaks;
	
		vector<float> max_intens_per_type;
		max_intens_per_type.resize(max_frag_idx+1,0); // 

		vector<float> intens;
		intens.resize(num_peaks,-1);
		float max_inten = -1;
		int j;
		for (j=0; j<num_peaks; j++)
		{
			const int   peak_idx = start_idx+j;
			const float peak_inten = peak_intens[peak_idx];
			const int   peak_frag_type = peak_frag_types[peak_idx];
				
			intens[j] = peak_inten;
			if (peak_inten>max_inten)
				max_inten=peak_inten;
			if (peak_inten>max_intens_per_type[peak_frag_type])
				max_intens_per_type[peak_frag_type]=peak_inten;
		}

		if (max_inten<=0)
		{
			continue;
		}

		const float one_over_max = 1.0/max_inten;
		for (j=0; j<num_peaks; j++)
			intens[j] *= one_over_max;

		for (j=0; j<=max_frag_idx; j++)
			max_intens_per_type[j]*= one_over_max;


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
					if ( (intens[j]/intens[k]) > 0.7)
						continue;
				
					weight_t pair_weight = 0.4 + 0.6*intens[k] - 0.4*(intens[j]/intens[k]);
				
					if (pair_weight<min_norm_weight)
						continue;

					// don't let very weak peaks participate in a missing peak pair
					if (intens[j]==0 && intens[k]<0.15)
						continue;

					phi.push_back(SamplePairWeight(start_idx+j,start_idx+k,pair_weight));
				}
		}
	}

	cout << "Created " << phi.size() << " pairs for " << peak_starts.size() << " peptide samples" << endl;
}


void PartitionModel::train_combined_partition_model( 
								const PeakRankModel *prank, 
								char *sample_file_path,
								int	_charge,
								int  _size_idx,
								int  _mobility,
								int  num_frags_to_choose, 
								char *report_dir,
								int  max_num_rounds,
								char *test_set_file,
								int	  test_peptide_length,
								char *stop_signal_file,
								weight_t max_weight_ratio)
{
	const float min_ratio_detected = 0.15;
	Config *config = prank->get_config();
	const vector<FragmentType>& all_fragments = config->get_all_fragments();
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const int num_all_frags = all_fragments.size();

	if (max_num_rounds<0)
		max_num_rounds = 1000;

	charge = _charge;
	size_idx = _size_idx;
	mobility   = _mobility;
	
	vector<TrainingPeptide> sample_tps, test_tps;
	vector<int> frag_counts;
	vector<int> frag_detected;
	vector<int> length_counts;

	
	if (this->feature_set_type <=2)
	{
		cout << "Error: train_combined_partition_model works only on frag_set 3 and up!" << endl;
		exit(1);
	}

	read_training_peptides_from_file(sample_file_path,sample_tps);

	cout << "Read " << sample_tps.size() << " training tps...";

	int num_tps_to_add = 0; // for de novo models
	if (prank->get_feature_set_type() == 4)
	{
		if (sample_tps.size()<15000)
			num_tps_to_add = 1;
		if (sample_tps.size()<10000)
			num_tps_to_add = 2;
		if (sample_tps.size()<5000)
			num_tps_to_add = 3;

		if (charge>=3 || size_idx>=3)
			num_tps_to_add /= 2;
	
		cout << "Adding at most " << num_tps_to_add << " per tp." << endl;
	}

	if (prank->get_feature_set_type() == 4)
	{
		cout << "Converting train: " << sample_tps.size() << " --> ";
		convert_tps_to_partial_denovo(config,sample_tps,num_tps_to_add,false);
		cout << sample_tps.size() << endl;
	}

	if (test_set_file)
	{
		read_training_peptides_from_file(test_set_file, test_tps);
		cout << "Read " << test_tps.size() << " test tps...";
		if (prank->get_feature_set_type() == 4)
		{
			cout << "Converting train: " << test_tps.size() << " --> ";
			convert_tps_to_partial_denovo(config,test_tps,num_tps_to_add);
			cout << test_tps.size() << endl;
		}
	}

	// Create intial report on dataset
	frag_counts.resize(num_all_frags,0);
	frag_detected.resize(num_all_frags,0);
	length_counts.resize(200,0);

	int numH=0,numK=0, numR=0;
	int i,f;
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
		
	vector<float> frag_ratios;
	frag_ratios.resize(all_fragments.size(),0);

	cout << "Detected fragments: " << endl;
	for (i=0; i<all_fragments.size(); i++)
	{
		const float ratio = (frag_counts[i]>0 ?  (float)frag_detected[i]/frag_counts[i] : 0);
		cout << all_fragments[i].label << "\t" << frag_detected[i] << " / " << 
			frag_counts[i] << "\t = " << setprecision(3) << ratio << endl;
		if (ratio>=min_ratio_detected)
			frag_ratios[i]=ratio;
	}
	cout << endl;

	// choose the most abundant fragments
	cout << "Selected: ";
	this->fragment_type_idxs.clear();
	int max_frag_idx=-1;
	for (f=0; f<num_frags_to_choose; f++)
	{
		float max_ratio=0;
		int max_idx=-1;
		int j;
		for (j=0; j<frag_ratios.size(); j++)
			if (frag_ratios[j]>max_ratio)
			{
				max_ratio=frag_ratios[j];
				max_idx=j;
			}
		if (max_idx<0)
			break;

		fragment_type_idxs.push_back(max_idx);
		frag_ratios[max_idx]=-1;
		cout << " " << all_fragments[max_idx].label;
		if (max_idx>max_frag_idx)
			max_frag_idx=max_idx;
	}
	cout << endl;

	if (fragment_type_idxs.size() == 0)
	{
		cout << "No fragments to train!" << endl;
		return;
	}


	if (prank->get_feature_set_type() == 3)
	{
		set_combined_feature_names_in_rankboost_model(prank);
	}
	else if (prank->get_feature_set_type() == 4)
	{
		set_combined_dnv_feature_names_in_rankboost_model(prank);
	}
	else if (prank->get_feature_set_type() == 5)
	{
		set_combined_simple_feature_names_in_rankboost_model(prank);
	}
	else
	{
		cout << "Error: no feature names supported for type " << prank->get_feature_set_type() << endl;
		exit(1);
	}
			
	// fill RankBoostSamples and create rank ds
	RankBoostDataset         rank_ds,test_ds;
	cout << "TRAINING..." << endl;
	cout << "Max weight ratio: " << max_weight_ratio << endl;

	// initialize and read the test set if it exists
	RankBoostDataset *test_set_ptr=NULL;
	if (test_set_file)
	{
		vector<float>			 test_peak_intens;
		vector<PeakStart>        test_peak_starts;
		vector<int>				 test_peak_frag_types;

		cout << endl << "Reading test tps..." << endl; 

		prank->read_training_peptides_into_combined_rank_boost_dataset(charge, size_idx, mobility,
			test_tps, test_ds, test_peak_intens, test_peak_starts, test_peak_frag_types);

		cout << "Test peak starts: " << test_peak_starts.size() << endl;
		cout << "Test peak intens: " << test_peak_intens.size() << endl;
		cout << "Test peak frags : " << test_peak_frag_types.size() << endl;

		cout << "Creating test phi list..." << endl;

		create_phi_list_from_combined_peak_samples(test_peak_intens, test_peak_frag_types,
			test_peak_starts,  max_frag_idx, test_ds.get_non_const_phi_support());

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
		cout << "test length " << test_peptide_length << endl << endl;
	}
		
	vector<float>		 train_peak_intens;
	vector<PeakStart>    train_peak_starts;
	vector<int>			 train_peak_frag_types;
	cout << endl << "reading training tps..." << endl;

	prank->read_training_peptides_into_combined_rank_boost_dataset(charge, size_idx, mobility,
		sample_tps, rank_ds, train_peak_intens, train_peak_starts, train_peak_frag_types);

	cout << "Train peak starts: " << train_peak_starts.size() << endl;
	cout << "Train peak intens: " << train_peak_intens.size() << endl;
	cout << "Train peak frags : " << train_peak_frag_types.size() << endl;

	cout << "create training phi list..." << endl;

	create_phi_list_from_combined_peak_samples(train_peak_intens, train_peak_frag_types,
			train_peak_starts,  max_frag_idx, rank_ds.get_non_const_phi_support());

	cout << "initializing boost..." << endl;
	rank_ds.compute_total_phi_weight();
	rank_ds.initialize_potenital_lists();
	rank_ds.initialize_binary_one_lists(combined_frag_boost_model.get_num_binary_features());
	rank_ds.initialize_binary_ordered_phi_lists(combined_frag_boost_model.get_ptr_to_binary_feature_names());
	rank_ds.initialzie_real_feature_table(combined_frag_boost_model.get_num_real_features());

	rank_ds.set_max_ratio_for_regular_update(max_weight_ratio);
				
	combined_frag_boost_model.init_rankboost_model_for_training(rank_ds,40,100);
		
	rank_ds.initialize_real_vote_lists(combined_frag_boost_model);

	char report_prefix[512];
	if (report_dir)
		sprintf(report_prefix,"%s/%s",report_dir, partition_name.c_str());
	
	if (report_prefix)
	{
		char name_buff[512];
		sprintf(name_buff,"%s_feature_summary.txt",report_prefix);
		ofstream sum_stream(name_buff);
		combined_frag_boost_model.summarize_features(rank_ds.get_samples(),sum_stream);
		sum_stream.close();
	}

	// write same info as the write_combined_model function performs
	vector<string> model_header_strings;
	ostringstream oss;
	oss << feature_set_type << " " << this->fragment_type_idxs.size();
	for (f=0; f<fragment_type_idxs.size(); f++)
		oss << " " << fragment_type_idxs[f];

	model_header_strings.push_back(oss.str());

	vector<idx_weight_pair> miss_pairs;
	combined_frag_boost_model.train_rankboost_model(rank_ds,
								max_num_rounds, 
								&miss_pairs, 
								test_set_ptr, 
								test_peptide_length, 
								report_prefix, 
								stop_signal_file,
								&model_header_strings);
		
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

		
		cout << "[...";
		combined_frag_boost_model.ouput_importance_ranked_feature_list(rank_ds,feature_stream);
		cout << " ...]" << endl;
		feature_stream.close();

		// write model (also compresses features and deletes the default values)
		sprintf(name_buff,"%s_model.txt",report_prefix);
		ofstream model_stream(name_buff);
	
		
		vector<peak_rank_stat> dummy_stats;

		cout << "Before reduction error: " << setprecision(7) << 
			combined_frag_boost_model.calc_prediction_error(test_ds, dummy_stats, test_peptide_length) << endl;
		write_combined_partition_model(name_buff);
		cout << "After  reduction error: " << setprecision(7) << 
			combined_frag_boost_model.calc_prediction_error(test_ds, dummy_stats, test_peptide_length) << endl;
	}	
}



void PartitionModel::fill_combined_peak_features(
									const PeakRankModel *prank,
									const  vector<int>& org_amino_acids,
									const int    cut_idx,
									const mass_t cut_mass,
									const PeptideSolution& sol,
									const FragmentType& frag,
									const int   position_idx_in_model_fragment_type_idxs,
									RankBoostSample& sample) const
{
	const vector<string>& model_aa_labels = prank->get_model_aa_labels();
	const vector<int>& session_aas_to_model_aas = prank->get_session_aas_to_model_aas();
	const int length = org_amino_acids.size();
	const int num_aas = model_aa_labels.size();
	const mass_t pm_with_19 = sol.pm_with_19;
	const int spec_charge = sol.charge;

	int rankFeatureIdx=0;
	int i;

	vector<int> amino_acids;
	prank->convert_aas_to_model_aas(org_amino_acids, amino_acids);

	if (amino_acids.size() != org_amino_acids.size())
	{
		cout << "Error: aa size mismatch!" << endl;
		exit(1);
	}

	if (cut_idx<=0 || cut_idx>=amino_acids.size())
	{
		cout << "Error: cut_idx is bad!" << endl;
		exit(1);
	}

	// need to use the special Idx variables and not the regular enumerations
	const int HisIdx = session_aas_to_model_aas[His];
	const int LysIdx = session_aas_to_model_aas[Lys];
	const int ArgIdx = session_aas_to_model_aas[Arg];
	const int SerIdx = session_aas_to_model_aas[Ser];
	const int ThrIdx = session_aas_to_model_aas[Thr];
	const int ProIdx = session_aas_to_model_aas[Pro];
	const int GlyIdx = session_aas_to_model_aas[Gly];
	const int AlaIdx = session_aas_to_model_aas[Ala];
	const int LeuIdx = session_aas_to_model_aas[Leu];
	const int AsnIdx = session_aas_to_model_aas[Asn];
	const int AspIdx = session_aas_to_model_aas[Asp];
	const int GluIdx = session_aas_to_model_aas[Glu];

	// special N C side aa indicators
	int num_nH=0, num_cH=0;
	int num_nK=0, num_cK=0;
	int num_nR=0, num_cR=0;
	
	for (i=0; i<cut_idx; i++)
	{
		if (amino_acids[i] == HisIdx)
			num_nH++;

		if (amino_acids[i] == LysIdx)
			num_nK++;

		if (amino_acids[i] == ArgIdx)
			num_nR++;
	}

	for (i=cut_idx; i<length; i++)
	{
		if (amino_acids[i] == HisIdx)
			num_cH++;

		if (amino_acids[i] == LysIdx)
			num_cK++;

		if (amino_acids[i] == ArgIdx)
			num_cR++;
	}

	// MASS / LOCATION FEATURES (REAL + BINARY)
	const mass_t max_detected_mass = prank->get_max_detected_mass();
	const mass_t exp_peak_mass = frag.calc_expected_mass(cut_mass,pm_with_19);
	const mass_t min_obs_mass = prank->calc_min_detected_mass(pm_with_19,spec_charge);
	const mass_t max_obs_mass = (pm_with_19>max_detected_mass ? max_detected_mass : pm_with_19);
	
	const float peak_mass_prop = ((exp_peak_mass - min_obs_mass)/(max_obs_mass - min_obs_mass));
	const float rounded_peak_prop = 0.02*floor(peak_mass_prop * 50.0);

	const mass_t dis_from_min = 20.0*floor((exp_peak_mass - min_obs_mass)*0.05);
	const mass_t dis_from_max = 20.0*floor((max_obs_mass  - exp_peak_mass)*0.05);
	const mass_t dis_from_minmax = (dis_from_min<dis_from_max ? dis_from_min : prank->get_max_detected_mass() - dis_from_max);

	const int RKH_n_combo_idx = calc_RKH_combo_idx(num_nR,num_nK,num_nH);
	const int RKH_c_combo_idx = calc_RKH_combo_idx(num_cR,num_cK,num_cH);
	const int RKH_pair_idx = (RKH_n_combo_idx * num_RKH_combos) + RKH_c_combo_idx;
	const float RKH_liniar_pair_idx = RKH_pair_matrix[RKH_n_combo_idx][RKH_c_combo_idx];

	const float cut_prop = 0.02 * floor((cut_mass / pm_with_19)*50);
	const float cut_n_mass = 20 * floor(cut_mass*0.05);
	const float cut_c_mass = 20 * floor((pm_with_19-cut_mass)*0.05);
	const int cut_dis_from_n = cut_idx;
	const int cut_dis_from_c = length-cut_idx;

	const int third=length/3;
	const int two_thirds = length - third;
	float cut_pos=NEG_INF;
	if (cut_idx<=4 || (cut_idx<=third && cut_idx<8))
	{
		cut_pos = (float)cut_idx;
	}
	else if (cut_idx>=length-4 || (cut_idx>=two_thirds && cut_idx>length-8))
	{
		cut_pos = (float)(20-length+cut_idx);
	}
	else
		cut_pos = 10.0 + cut_prop;

	const int n_aa = amino_acids[cut_idx-1];
	const int c_aa = amino_acids[cut_idx];
	const int n_term_aa = amino_acids[0];
	const int c_term_aa = amino_acids[length-1];


	sample.clear();

	// first feature is always 0, tells what fragment this is
	sample.add_real_feature(0,(float)position_idx_in_model_fragment_type_idxs);

	// this is the new starting index for features in this sample. The offset
	// is computaed according to the fragment type
	int f_idx =1 + position_idx_in_model_fragment_type_idxs * num_features_per_frag; // num_features_per_frag_type

	if (num_features_per_frag==0)
	{
		cout << "Error: num_features_per_frag is 0!" << endl;
		exit(1);
	}

	// add general position features
	sample.add_real_feature(f_idx++,dis_from_minmax);
	sample.add_real_feature(f_idx++,rounded_peak_prop);
	sample.add_real_feature(f_idx++,RKH_pair_idx);
	sample.add_real_feature(f_idx++,RKH_liniar_pair_idx);
	sample.add_real_feature(f_idx++,cut_prop);
	sample.add_real_feature(f_idx++,cut_pos);
	sample.add_real_feature(f_idx++,cut_n_mass);
	sample.add_real_feature(f_idx++,cut_c_mass);
	sample.add_real_feature(f_idx++,(float)cut_dis_from_n);
	sample.add_real_feature(f_idx++,(float)cut_dis_from_c);
	
	const float position_var = rounded_peak_prop; // this will be used for all amino acid features

	// fill aa flanking features N side
	if (cut_idx>0)
		sample.add_real_feature(f_idx+amino_acids[cut_idx-1],position_var);
	f_idx+=num_aas;
	if (cut_idx>1)
		sample.add_real_feature(f_idx+amino_acids[cut_idx-2],position_var);
	f_idx+=num_aas;
	if (cut_idx>2)
		sample.add_real_feature(f_idx+amino_acids[cut_idx-3],position_var);
	f_idx+=num_aas;

	// fill aa flanking features C side
	if (cut_idx<length)
		sample.add_real_feature(f_idx+amino_acids[cut_idx],position_var);
	f_idx+=num_aas;
	if (cut_idx<length-1)
		sample.add_real_feature(f_idx+amino_acids[cut_idx+1],position_var);
	f_idx+=num_aas;
	if (cut_idx<length-2)
		sample.add_real_feature(f_idx+amino_acids[cut_idx+2],position_var);
	f_idx+=num_aas;

	// fill cut pair features X-Y
	sample.add_real_feature(f_idx+(n_aa*num_aas+c_aa),position_var);
	f_idx+=(num_aas*num_aas);

	// fill aa count features (up top 3 aa's away from cut)
	vector<int> n_aa_counts, c_aa_counts;
	n_aa_counts.resize(num_aas+1,0);
	c_aa_counts.resize(num_aas+1,0);

	for (i=0; i<cut_idx-3; i++)
		n_aa_counts[amino_acids[i]]++;

	for (i=cut_idx+3; i<length; i++)
		c_aa_counts[amino_acids[i]]++;

	int a;
	for (a=0; a<num_aas; a++)
		if (n_aa_counts[a]>0)
			sample.add_real_feature(f_idx+a,n_aa_counts[a]);
	f_idx+=num_aas;

	for (a=0; a<num_aas; a++)
		if (c_aa_counts[a]>0)
			sample.add_real_feature(f_idx+a,c_aa_counts[a]);
	f_idx+=num_aas;

	// add Cut\terminal features
	sample.add_real_feature(f_idx+n_term_aa,cut_pos);
	f_idx+=num_aas;
	sample.add_real_feature(f_idx+c_term_aa,cut_pos);
	f_idx+=num_aas;

	// X|
	if (cut_idx>0)
	{
		sample.add_real_feature(f_idx+amino_acids[cut_idx-1],cut_pos);
		f_idx+=num_aas;
		if (c_term_aa == LysIdx)
			sample.add_real_feature(f_idx+amino_acids[cut_idx-1],cut_pos);
		f_idx+=num_aas;
		if (c_term_aa == ArgIdx)
			sample.add_real_feature(f_idx+amino_acids[cut_idx-1],cut_pos);
		f_idx+=num_aas;
	}
	else
		f_idx+=3*num_aas;

	// X_|
	if (cut_idx>1)
	{
		sample.add_real_feature(f_idx+amino_acids[cut_idx-2],cut_pos);
		f_idx+=num_aas;
		if (c_term_aa == LysIdx)
			sample.add_real_feature(f_idx+amino_acids[cut_idx-2],cut_pos);
		f_idx+=num_aas;
		if (c_term_aa == ArgIdx)
			sample.add_real_feature(f_idx+amino_acids[cut_idx-2],cut_pos);
		f_idx+=num_aas;
	}
	else
		f_idx+=3*num_aas;

	// |X
	if (cut_idx<length)
	{
		sample.add_real_feature(f_idx+amino_acids[cut_idx],cut_pos);
		f_idx+=num_aas;
		if (c_term_aa == LysIdx)
			sample.add_real_feature(f_idx+amino_acids[cut_idx],cut_pos);
		f_idx+=num_aas;
		if (c_term_aa == ArgIdx)
			sample.add_real_feature(f_idx+amino_acids[cut_idx],cut_pos);
		f_idx+=num_aas;
	}
	else
		f_idx+=3*num_aas;

	// |_X
	if (cut_idx<length-1)
	{
		sample.add_real_feature(f_idx+amino_acids[cut_idx+1],cut_pos);
		f_idx+=num_aas;
		if (c_term_aa == LysIdx)
			sample.add_real_feature(f_idx+amino_acids[cut_idx+1],cut_pos);
		f_idx+=num_aas;
		if (c_term_aa == ArgIdx)
			sample.add_real_feature(f_idx+amino_acids[cut_idx+1],cut_pos);
		f_idx+=num_aas;
	}
	else
		f_idx+=3*num_aas;

}



void PartitionModel::fill_combined_simple_peak_features(
								 const PeakRankModel *prank,
								 const  vector<int>& org_amino_acids,
								 const int    cut_idx,
								 const mass_t cut_mass,
								 const PeptideSolution& sol,
								 const FragmentType& frag,
								 const int   position_idx_in_model_fragment_type_idxs,
								 RankBoostSample& sample,
								 bool  verbose) const
{
	const vector<string>& model_aa_labels = prank->get_model_aa_labels();
	const vector<int>& session_aas_to_model_aas = prank->get_session_aas_to_model_aas();
	const int length = org_amino_acids.size();
	const int num_aas = model_aa_labels.size();
	const mass_t pm_with_19 = sol.pm_with_19;
	const int spec_charge = sol.charge;

	if (verbose)
	{
		cout << "Filling: " << frag.label<< endl;
		cout << "Pos: " << position_idx_in_model_fragment_type_idxs << endl;
	}

	int rankFeatureIdx=0;
	int i;

	vector<int> amino_acids;
	prank->convert_aas_to_model_aas(org_amino_acids, amino_acids);

	if (amino_acids.size() != org_amino_acids.size())
	{
		cout << "Error: aa size mismatch!" << endl;
		exit(1);
	}

	if (cut_idx<=0 || cut_idx>=amino_acids.size())
	{
		cout << "Error: cut_idx is bad!" << endl;
		exit(1);
	}

	// need to use the special Idx variables and not the regular enumerations
	const int HisIdx = session_aas_to_model_aas[His];
	const int LysIdx = session_aas_to_model_aas[Lys];
	const int ArgIdx = session_aas_to_model_aas[Arg];
	const int SerIdx = session_aas_to_model_aas[Ser];
	const int ThrIdx = session_aas_to_model_aas[Thr];
	const int ProIdx = session_aas_to_model_aas[Pro];
	const int GlyIdx = session_aas_to_model_aas[Gly];
	const int AlaIdx = session_aas_to_model_aas[Ala];
	const int LeuIdx = session_aas_to_model_aas[Leu];
	const int AsnIdx = session_aas_to_model_aas[Asn];
	const int AspIdx = session_aas_to_model_aas[Asp];
	const int GluIdx = session_aas_to_model_aas[Glu];

	// special N C side aa indicators
	int num_nH=0, num_cH=0;
	int num_nK=0, num_cK=0;
	int num_nR=0, num_cR=0;
	
	for (i=0; i<cut_idx; i++)
	{
		if (amino_acids[i] == HisIdx)
			num_nH++;

		if (amino_acids[i] == LysIdx)
			num_nK++;

		if (amino_acids[i] == ArgIdx)
			num_nR++;
	}

	for (i=cut_idx; i<length; i++)
	{
		if (amino_acids[i] == HisIdx)
			num_cH++;

		if (amino_acids[i] == LysIdx)
			num_cK++;

		if (amino_acids[i] == ArgIdx)
			num_cR++;
	}

	// MASS / LOCATION FEATURES (REAL + BINARY)
	const mass_t max_detected_mass = prank->get_max_detected_mass();
	const mass_t exp_peak_mass = frag.calc_expected_mass(cut_mass,pm_with_19);
	const mass_t min_obs_mass = prank->calc_min_detected_mass(pm_with_19,spec_charge);
	const mass_t max_obs_mass = (pm_with_19>max_detected_mass ? max_detected_mass : pm_with_19);
	
	const float peak_mass_prop = ((exp_peak_mass - min_obs_mass)/(max_obs_mass - min_obs_mass));
	const float rounded_peak_prop = 0.05*floor(peak_mass_prop * 20.0);

	const mass_t dis_from_min = 20.0*floor((exp_peak_mass - min_obs_mass)*0.05);
	const mass_t dis_from_max = 20.0*floor((max_obs_mass  - exp_peak_mass)*0.05);

	const int cut_dis_from_n = cut_idx;
	const int cut_dis_from_c = length-cut_idx;


	const int n_aa = amino_acids[cut_idx-1];
	const int c_aa = amino_acids[cut_idx];
	const int n_term_aa = amino_acids[0];
	const int c_term_aa = amino_acids[length-1];

	const int half = amino_acids.size() /2 ;
	float cut_pos=NEG_INF;
	if (cut_idx<half)
	{
		cut_pos = (float)cut_idx;
		if (cut_idx>9)
			cut_pos=10;
	}
	else
	{
		cut_pos = 20 - (amino_acids.size()-cut_idx);
		if (cut_pos<10)
			cut_pos=10;
	}
	

	sample.clear();

	// first feature is always 0, tells what fragment this is
	sample.add_real_feature(0,(float)position_idx_in_model_fragment_type_idxs);

	// this is the new starting index for features in this sample. The offset
	// is computaed according to the fragment type
	int f_idx =1 + position_idx_in_model_fragment_type_idxs * num_features_per_frag; // num_features_per_frag_type
	
	if (num_features_per_frag==0)
	{
		cout << "Error: num_features_per_frag is 0!" << endl;
		exit(1);
	}

//	cout << "Starting with :" << f_idx << " (" << num_features_per_frag << ")" << endl;

	// add general position features
	if (dis_from_min<dis_from_max)
	{
		sample.add_real_feature(f_idx,dis_from_min);
	}
	else
	{
		sample.add_real_feature(f_idx+1,dis_from_max);
	}
	f_idx+=2;

	sample.add_real_feature(f_idx++,rounded_peak_prop);

//	sample.add_real_feature(f_idx++,(float)cut_pos);
	
	const float position_var = 1.0;//rounded_peak_prop; // this will be used for all amino acid features

	// fill aa flanking features N side
	if (cut_idx>0)
		sample.add_real_feature(f_idx+amino_acids[cut_idx-1],position_var);
	f_idx+=num_aas;
	if (cut_idx>1)
		sample.add_real_feature(f_idx+amino_acids[cut_idx-2],position_var);
	f_idx+=num_aas;
	if (cut_idx>2)
		sample.add_real_feature(f_idx+amino_acids[cut_idx-3],position_var);
	f_idx+=num_aas;

	// fill aa flanking features C side
	if (cut_idx<length)
		sample.add_real_feature(f_idx+amino_acids[cut_idx],position_var);
	f_idx+=num_aas;
	if (cut_idx<length-1)
		sample.add_real_feature(f_idx+amino_acids[cut_idx+1],position_var);
	f_idx+=num_aas;
	if (cut_idx<length-2)
		sample.add_real_feature(f_idx+amino_acids[cut_idx+2],position_var);
	f_idx+=num_aas;

//	// fill cut pair features X-Y
//	sample.add_real_feature(f_idx+(n_aa*num_aas+c_aa),1);
//	f_idx+=(num_aas*num_aas);

	// fill aa count features (up top 3 aa's away from cut)
	vector<int> n_aa_counts, c_aa_counts;
	n_aa_counts.resize(num_aas+1,0);
	c_aa_counts.resize(num_aas+1,0);

	for (i=0; i<cut_idx-3; i++)
		n_aa_counts[amino_acids[i]]++;

	for (i=cut_idx+3; i<length; i++)
		c_aa_counts[amino_acids[i]]++;

	int a;
	for (a=0; a<num_aas; a++)
		if (n_aa_counts[a]>0)
			sample.add_real_feature(f_idx+a,n_aa_counts[a]);
	f_idx+=num_aas;

	for (a=0; a<num_aas; a++)
		if (c_aa_counts[a]>0)
			sample.add_real_feature(f_idx+a,c_aa_counts[a]);
	f_idx+=num_aas;

	// add Cut\terminal features
	sample.add_real_feature(f_idx+n_term_aa,position_var);
	f_idx+=num_aas;
	sample.add_real_feature(f_idx+c_term_aa,position_var);
	f_idx+=num_aas;
}


void PartitionModel::set_combined_simple_feature_names_in_rankboost_model(const PeakRankModel *prank)
{
	const vector<string>& model_aa_labels = prank->get_model_aa_labels();
	const int num_aas = model_aa_labels.size();
	Config *config = prank->get_config();
	vector<string> frag_feature_names;
	vector<string> all_feature_names;

	frag_feature_names.clear();

	frag_feature_names.push_back("Dis Min");
	frag_feature_names.push_back("Dis Max");
	frag_feature_names.push_back("Peak prop [Min-Max]");
//	frag_feature_names.push_back("Cut Pos");

	int i;
	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|_");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "_|_");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "__|_");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is _|" + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is _|_" + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is _|__" + model_aa_labels[i] );

//	for (i=0; i<num_aas; i++)
//	{
//		int j;
//		for (j=0; j<num_aas; j++)
//			frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|" + model_aa_labels[j]);
//	}

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("# N-side " + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("# C-side " + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("N-term aa is " + model_aa_labels[i] + " cut pos");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("C-term aa is " + model_aa_labels[i] + " cut pos");



	// create all feature names
	all_feature_names.clear();
	all_feature_names.push_back("FRAG TYPE");
	int f;
	for (f=0; f<this->fragment_type_idxs.size(); f++)
	{
		const string frag_label = config->get_fragment(fragment_type_idxs[f]).label;
		int i;
		for (i=0; i<frag_feature_names.size(); i++)
		{
			string feature_name = frag_label + ": " + frag_feature_names[i];
			all_feature_names.push_back(feature_name);
		}
	}

	num_features_per_frag = frag_feature_names.size();
	vector<string> empty;
	empty.clear();
	combined_frag_boost_model.init_rankboost_model_feature_names(empty, all_feature_names);
}


void PartitionModel::set_combined_feature_names_in_rankboost_model(const PeakRankModel *prank)
{
	const vector<string>& model_aa_labels = prank->get_model_aa_labels();
	const int num_aas = model_aa_labels.size();
	Config *config = prank->get_config();
	vector<string> frag_feature_names;
	vector<string> all_feature_names;

	frag_feature_names.clear();

	frag_feature_names.push_back("Dis Min/Max");
	frag_feature_names.push_back("Peak prop [Min-Max]");
	frag_feature_names.push_back("RHK pair idx");
	frag_feature_names.push_back("RHK liniar pair idx");
	frag_feature_names.push_back("Cut prop [0-M+19]");
	frag_feature_names.push_back("Cut pos");
	frag_feature_names.push_back("Cut N mass");
	frag_feature_names.push_back("Cut C mass");
	frag_feature_names.push_back("Cut idx from N");
	frag_feature_names.push_back("Cut idx from C");
//	push_back_all_RHK_pairs(frag_feature_names,"Cut prop");

	int i;
	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|_");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "_|_");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "__|_");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is _|" + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is _|_" + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is _|__" + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
	{
		int j;
		for (j=0; j<num_aas; j++)
			frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|" + model_aa_labels[j]);
	}

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("# N-side " + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("# C-side " + model_aa_labels[i] );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("N-term aa is " + model_aa_labels[i] + ", cut pos" );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("C-term aa is " + model_aa_labels[i] + ", cut pos" );

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|, cut pos");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|, cut pos, C-term is K");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "|, cut pos, C-term is R");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "_|, cut pos");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "_|, cut pos, C-term is K");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is " + model_aa_labels[i] + "_|, cut pos, C-term is R");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |" + model_aa_labels[i] + ", cut pos");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |" + model_aa_labels[i] + ", cut pos, C-term is K");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |" + model_aa_labels[i] + ", cut pos, C-term is R");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |_" + model_aa_labels[i] + ", cut pos");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |_" + model_aa_labels[i] + ", cut pos, C-term is K");

	for (i=0; i<num_aas; i++)
		frag_feature_names.push_back("Cut is |_" + model_aa_labels[i] + ", cut pos, C-term is R");

	// create all feature names
	all_feature_names.clear();
	all_feature_names.push_back("FRAG TYPE");
	int f;
	for (f=0; f<this->fragment_type_idxs.size(); f++)
	{
		const string frag_label = config->get_fragment(fragment_type_idxs[f]).label;
		int i;
		for (i=0; i<frag_feature_names.size(); i++)
		{
			string feature_name = frag_label + ": " + frag_feature_names[i];
			all_feature_names.push_back(feature_name);
		}
	}

	num_features_per_frag = frag_feature_names.size();
	vector<string> empty;
	empty.clear();
	combined_frag_boost_model.init_rankboost_model_feature_names(empty, all_feature_names);
}



struct inten_trip {
	inten_trip() : frag_type_idx(int(NEG_INF)), idx(int(NEG_INF)), inten(int(NEG_INF)) {};
	inten_trip(int _f, int _i, float _n) : frag_type_idx(_f), idx(_i), inten(_n) {};
	bool operator< (const inten_trip& other) const
	{
		return inten>other.inten;
	}

	int frag_type_idx;
	int idx;
	float inten;
};


/***********************************************************************
calculates the relative rank of the cuts. If intensity is 0, the rank 
999 is given (rank[0] also is 999). The calculation is done across the
different fragment types given (not only one kind).
************************************************************************/
void calc_combined_peak_ranks(const vector< vector<float> >& intens, 
							  vector< vector<int> >& peak_ranks)
{
	vector<inten_trip> trips;
	int i;
	for (i=0; i<intens.size(); i++)
	{
		int j;
		for (j=1; j<intens[i].size(); j++)
			if (intens[i][j]>0)
			trips.push_back(inten_trip(i,j,intens[i][j]));
	}
	sort(trips.begin(),trips.end());
	
	peak_ranks.clear();
	peak_ranks.resize(intens.size());
	for (i=0; i<intens.size(); i++)
		peak_ranks[i].resize(intens[i].size(),999);

	for (i=0; i<trips.size(); i++)
		peak_ranks[trips[i].frag_type_idx][trips[i].idx]=i;
}


void PeakRankModel::read_training_peptides_into_combined_rank_boost_dataset(
										int spec_charge,
										int size_idx,
										int mobility,
										const vector<TrainingPeptide>& sample_tps,
										RankBoostDataset& rank_ds,
										vector<float>& peak_intens,
										vector<PeakStart>& peak_starts,
										vector<int>& peak_frag_types) const
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const PartitionModel* part_model = partition_models[spec_charge][size_idx][mobility];
	const vector<int>& frag_type_idxs = part_model->fragment_type_idxs;

	double random_thresh = 0.8;
	if (sample_tps.size()<12000)
		random_thresh=1.0;

	if (sample_tps.size()>15000)
		random_thresh=0.7;

	if (spec_charge>=3 && size_idx >=2 && mobility > MOBILE && sample_tps.size()>12000 )
		random_thresh = 0.7;

	if (sample_tps.size()>30000)
		random_thresh = 0.6;

	rank_ds.clear();
	peak_intens.clear();
	peak_frag_types.clear();
	peak_starts.clear();
	
	if (this->feature_set_type <=2)
	{
		cout << "Error: read_training_peptides_into_combined_rank_boost_dataset works only on frag set 3 and up!" << endl;
		exit(1);
	}

	int i;
	for (i=0; i<sample_tps.size(); i++)
	{
		const TrainingPeptide& org_tp = sample_tps[i];

		if (random_thresh<1.0 && myRandom() > random_thresh)
			continue;

		// create a new tp only with the frags we are interested in, and in the order
		// of frag_type_idxs
		// the first intensity and last intensity are -999 (so the intens size is length+1)
		vector< vector<float> > all_intens;
		all_intens.resize(frag_type_idxs.size());
		int f;
		for (f=0; f<frag_type_idxs.size(); f++)
		{
			const int frag_type_idx = frag_type_idxs[f];
			const int frag_pos = org_tp.get_frag_idx_pos(frag_type_idx);

			if (frag_pos<0)
			{
				all_intens[f].clear();
				continue;
			}
			all_intens[f]=org_tp.intens[frag_pos];
			all_intens[f][0]=-999.0;

			if (all_intens[f].size() == org_tp.length)
				all_intens[f].push_back(-999.0);
		}
		

		// scan all cuts for max inten
		float max_ann_inten=0;
		for (f=0; f<all_intens.size(); f++)
		{
			if (all_intens[f].size()==0)
				continue;

			int c;
			for (c=0; c<all_intens[f].size()-1; c++)
				if (all_intens[f][c]>max_ann_inten)
					max_ann_inten=all_intens[f][c];
		}

		if (max_ann_inten<=0)
			continue;

		// gives integer values to cut idxs/frags according to intensity
		vector< vector<int> > peak_ranks;
		calc_combined_peak_ranks(all_intens, peak_ranks);

		mass_t c_term_mass = org_tp.n_mass;
		int aa_idx;
		for (aa_idx=0; aa_idx<org_tp.amino_acids.size(); aa_idx++)
			c_term_mass+=aa2mass[org_tp.amino_acids[aa_idx]];

		PeakStart ps;
		ps.peptide_sample_idx = i;
		ps.peak_start_idx=rank_ds.get_num_samples();
		ps.num_peaks=0;

		PeptideSolution sol;
		sol.pep.set_peptide_aas(org_tp.amino_acids);
		sol.pep.calc_mass(config);
		sol.charge = org_tp.charge;
		sol.most_basic_aa_removed_from_n = org_tp.best_n_removed;
		sol.most_basic_aa_removed_from_c = org_tp.best_c_removed;
		sol.reaches_n_terminal = (org_tp.n_mass<1.0);
		sol.reaches_c_terminal = (c_term_mass + 25 > org_tp.pm_with_19);
		sol.pm_with_19 = org_tp.pm_with_19;

		for (f=0; f<frag_type_idxs.size(); f++)
		{
			if (all_intens[f].size()==0)
				continue;

			const int frag_type_idx = frag_type_idxs[f];
			const vector<float> & frag_intens = all_intens[f];
			const FragmentType& frag = config->get_fragment(frag_type_idx);

			// fill samples
			mass_t cut_mass=org_tp.n_mass;
			int cut_idx;
			for (cut_idx=0; cut_idx<frag_intens.size(); cut_idx++)
			{
				if (cut_idx>0)
					cut_mass+=aa2mass[org_tp.amino_acids[cut_idx-1]];

				if (frag_intens[cut_idx]>=0)
				{
					RankBoostSample rbs;
					
					if (feature_set_type == 3)
					{
						part_model->fill_combined_peak_features(this,
								 org_tp.amino_acids, cut_idx, cut_mass, sol, frag, f, rbs);
					}
					else if (feature_set_type == 4)
					{
						part_model->fill_combined_dnv_peak_features(this, org_tp.n_mass, c_term_mass, 
							org_tp.amino_acids, cut_idx, cut_mass, sol, frag, f, rbs);
					}
					else if (feature_set_type == 5)
					{
						part_model->fill_combined_simple_peak_features(this,
								 org_tp.amino_acids, cut_idx, cut_mass, sol, frag, f, rbs);
					}
					else
					{
						cout << "Error: feature set type not supported: " << feature_set_type << endl;
						exit(1);
					}

					rbs.groupIndex = i;
					rbs.rank_in_group = peak_ranks[f][cut_idx];

					rbs.tag1 = rank_ds.get_num_samples(); // sample idx
					rbs.tag2 = cut_idx;                   // cut idx
					rbs.tag3 = org_tp.amino_acids.size(); // peptide length, this tag is used to filter
													      // the error estimates for a specific length
					rank_ds.add_sample(rbs);
					peak_intens.push_back(frag_intens[cut_idx]);
					peak_frag_types.push_back(frag_type_idx);
				}	
			}
		}

		ps.num_peaks = rank_ds.get_num_samples() - ps.peak_start_idx;
		peak_starts.push_back(ps);
	}

	rank_ds.set_num_groups(sample_tps.size());
}



void PeakRankModel::train_all_combined_partition_models(
								int frag_fill_type,
								char *prefix_path,
								int	  sel_charge,
								int   sel_size_idx,
								int	  sel_mobility,
								int	  num_frags,
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

				cout << "Training combined models for charge " << charge << " size " << size_idx<< " mobility " <<
					mobility << endl;
				cout << "Max weight ratio " << max_weight_ratio << endl;

				partition_models[charge][size_idx][mobility]->set_partition_name(get_peak_rank_model_name(),
					charge, size_idx, mobility);

				partition_models[charge][size_idx][mobility]->set_feature_set_type(feature_set_type);

			
				partition_models[charge][size_idx][mobility]->train_combined_partition_model(this,
					tr_file,charge,size_idx,mobility, num_frags, report_dir,
					num_rounds,test_set, test_peptide_length, stop_signal_file, 
					max_weight_ratio);

				cout << endl;
			}
		}
	}
}



