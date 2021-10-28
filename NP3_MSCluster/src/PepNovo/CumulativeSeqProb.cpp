#include "CumulativeSeqProb.h"
#include "AllScoreModels.h"
#include "PeptideRankScorer.h"
#include "DeNovoSolutions.h"
#include "PepNovo_auxfun.h"

const int num_special_idxs = sizeof(special_rank_feature_idxs)/sizeof(int);
const int num_rank_level_idxs = sizeof(rank_level_idxs)/sizeof(int);

const char * Cumulative_Seq_Prob_CSP_names[]={
"CSP_CONST",	"CSP_IND_QUAL_10",	"CSP_IND_QUAL_25",	"CSP_IND_QUAL_50",	"CSP_IND_QUAL_75",	"CSP_IND_QUAL_90",	"CSP_IND_QUAL_95",	"CSP_IND_QUAL_98",	"CSP_IND_QUAL_100",	"CSP_IND_RANK_SCORE_TOP_LESS_MINUS_5",	"CSP_IND_RANK_SCORE_TOP_LESS_MINUS_4",	"CSP_IND_RANK_SCORE_TOP_LESS_MINUS_3",	"CSP_IND_RANK_SCORE_TOP_LESS_MINUS_2",	"CSP_IND_RANK_SCORE_TOP_LESS_MINUS_1",	"CSP_IND_RANK_SCORE_TOP_LESS_0",	
"CSP_IND_RANK_SCORE_TOP_LESS_1",	"CSP_IND_RANK_SCORE_TOP_LESS_2",	"CSP_IND_RANK_SCORE_TOP_LESS_3",	"CSP_IND_RANK_SCORE_TOP_LESS_4",	"CSP_IND_RANK_SCORE_TOP_LESS_5",	"CSP_IND_RANK_SCORE_TOP_LESS_6",	"CSP_IND_RANK_SCORE_TOP_LESS_7",	"CSP_IND_RANK_SCORE_MORE_7",	"CSP_IND_PNV_SCORE_TOP_LESS_0",	"CSP_IND_PNV_SCORE_TOP_LESS_10",	"CSP_IND_PNV_SCORE_TOP_LESS_20",	"CSP_IND_PNV_SCORE_TOP_LESS_30",	
"CSP_IND_PNV_SCORE_TOP_LESS_40",	"CSP_IND_PNV_SCORE_TOP_LESS_50",	"CSP_IND_PNV_SCORE_TOP_LESS_60",	"CSP_IND_PNV_SCORE_TOP_LESS_70",	"CSP_IND_PNV_SCORE_TOP_LESS_80",	"CSP_IND_PNV_SCORE_TOP_LESS_90",	"CSP_IND_PNV_SCORE_TOP_LESS_100",	"CSP_IND_PNV_SCORE_TOP_LESS_110",	"CSP_IND_PNV_SCORE_TOP_LESS_120",	"CSP_IND_PNV_SCORE_TOP_LESS_130",	"CSP_IND_PNV_SCORE_TOP_LESS_140",	"CSP_IND_PNV_SCORE_TOP_MORE_140",	
"CSP_IND_AVG_PROB_TOP_LESS_10",	"CSP_IND_AVG_PROB_TOP_LESS_25",	"CSP_IND_AVG_PROB_TOP_LESS_50",	"CSP_IND_AVG_PROB_TOP_LESS_75",	"CSP_IND_AVG_PROB_TOP_LESS_90",	"CSP_IND_AVG_PROB_TOP_LESS_95",	"CSP_IND_AVG_PROB_TOP_MORE_95",	"CSP_0_CONST",	"CSP_0_REL_RANK",	"CSP_0_RANK_SCR",	"CSP_0_DEL_RANK",	"CSP_0_AVG_PROB",	"CSP_0_MIN_PROB_1",	"CSP_0_MIN_PROB_2",	"CSP_0_NUM_AA",	"CSP_1_CONST",	
"CSP_1_REL_RANK",	"CSP_1_RANK_SCR",	"CSP_1_DEL_RANK",	"CSP_1_AVG_PROB",	"CSP_1_MIN_PROB_1",	"CSP_1_MIN_PROB_2",	"CSP_1_NUM_AA",	"CSP_2_CONST",	"CSP_2_REL_RANK",	"CSP_2_RANK_SCR",	"CSP_2_DEL_RANK",	"CSP_2_AVG_PROB",	"CSP_2_MIN_PROB_1",	"CSP_2_MIN_PROB_2",	"CSP_2_NUM_AA",	"CSP_3_CONST",	"CSP_3_REL_RANK",	"CSP_3_RANK_SCR",	"CSP_3_DEL_RANK",	"CSP_3_AVG_PROB",	"CSP_3_MIN_PROB_1",	
"CSP_3_MIN_PROB_2",	"CSP_3_NUM_AA",	"CSP_4_CONST",	"CSP_4_REL_RANK",	"CSP_4_RANK_SCR",	"CSP_4_DEL_RANK",	"CSP_4_AVG_PROB",	"CSP_4_MIN_PROB_1",	"CSP_4_MIN_PROB_2",	"CSP_4_NUM_AA",	"CSP_IND_RANK_SCORE_4_LESS_MINUS_5",	"CSP_IND_RANK_SCORE_4_LESS_MINUS_4",	"CSP_IND_RANK_SCORE_4_LESS_MINUS_3",	"CSP_IND_RANK_SCORE_4_LESS_MINUS_2",	"CSP_IND_RANK_SCORE_4_LESS_MINUS_1",	"CSP_IND_RANK_SCORE_4_LESS_0",	
"CSP_IND_RANK_SCORE_4_LESS_1",	"CSP_IND_RANK_SCORE_4_LESS_2",	"CSP_IND_RANK_SCORE_4_LESS_3",	"CSP_IND_RANK_SCORE_4_LESS_4",	"CSP_IND_RANK_SCORE_4_LESS_5",	"CSP_IND_RANK_SCORE_4_LESS_6",	"CSP_IND_RANK_SCORE_4_LESS_7",	"CSP_IND_RANK_SCORE_4_MORE_7",	"CSP_IND_PNV_SCORE_4_LESS_0",	"CSP_IND_PNV_SCORE_4_LESS_10",	"CSP_IND_PNV_SCORE_4_LESS_20",	"CSP_IND_PNV_SCORE_4_LESS_30",	"CSP_IND_PNV_SCORE_4_LESS_40",	
"CSP_IND_PNV_SCORE_4_LESS_50",	"CSP_IND_PNV_SCORE_4_LESS_60",	"CSP_IND_PNV_SCORE_4_LESS_70",	"CSP_IND_PNV_SCORE_4_LESS_80",	"CSP_IND_PNV_SCORE_4_LESS_90",	"CSP_IND_PNV_SCORE_4_LESS_100",	"CSP_IND_PNV_SCORE_4_LESS_110",	"CSP_IND_PNV_SCORE_4_LESS_120",	"CSP_IND_PNV_SCORE_4_LESS_130",	"CSP_IND_PNV_SCORE_4_LESS_140",	"CSP_IND_PNV_SCORE_4_MORE_140",	"CSP_IND_AVG_PROB_4_LESS_10",	"CSP_IND_AVG_PROB_4_LESS_25",	
"CSP_IND_AVG_PROB_4_LESS_50",	"CSP_IND_AVG_PROB_4_LESS_75",	"CSP_IND_AVG_PROB_4_LESS_90",	"CSP_IND_AVG_PROB_4_LESS_95",	"CSP_IND_AVG_PROB_4_MORE_95",	"CSP_5_CONST",	"CSP_5_REL_RANK",	"CSP_5_RANK_SCR",	"CSP_5_DEL_RANK",	"CSP_5_AVG_PROB",	"CSP_5_MIN_PROB_1",	"CSP_5_MIN_PROB_2",	"CSP_5_NUM_AA",	"CSP_6_CONST",	"CSP_6_REL_RANK",	"CSP_6_RANK_SCR",	"CSP_6_DEL_RANK",	"CSP_6_AVG_PROB",	
"CSP_6_MIN_PROB_1",	"CSP_6_MIN_PROB_2",	"CSP_6_NUM_AA",	"CSP_7_CONST",	"CSP_7_REL_RANK",	"CSP_7_RANK_SCR",	"CSP_7_DEL_RANK",	"CSP_7_AVG_PROB",	"CSP_7_MIN_PROB_1",	"CSP_7_MIN_PROB_2",	"CSP_7_NUM_AA",	"CSP_8_CONST",	"CSP_8_REL_RANK",	"CSP_8_RANK_SCR",	"CSP_8_DEL_RANK",	"CSP_8_AVG_PROB",	"CSP_8_MIN_PROB_1",	"CSP_8_MIN_PROB_2",	"CSP_8_NUM_AA",	"CSP_IND_RANK_SCORE_8_LESS_MINUS_5",	
"CSP_IND_RANK_SCORE_8_LESS_MINUS_4",	"CSP_IND_RANK_SCORE_8_LESS_MINUS_3",	"CSP_IND_RANK_SCORE_8_LESS_MINUS_2",	"CSP_IND_RANK_SCORE_8_LESS_MINUS_1",	"CSP_IND_RANK_SCORE_8_LESS_0",	"CSP_IND_RANK_SCORE_8_LESS_1",	"CSP_IND_RANK_SCORE_8_LESS_2",	"CSP_IND_RANK_SCORE_8_LESS_3",	"CSP_IND_RANK_SCORE_8_LESS_4",	"CSP_IND_RANK_SCORE_8_LESS_5",	"CSP_IND_RANK_SCORE_8_LESS_6",	"CSP_IND_RANK_SCORE_8_LESS_7",	
"CSP_IND_RANK_SCORE_8_MORE_7",	"CSP_IND_PNV_SCORE_8_LESS_0",	"CSP_IND_PNV_SCORE_8_LESS_10",	"CSP_IND_PNV_SCORE_8_LESS_20",	"CSP_IND_PNV_SCORE_8_LESS_30",	"CSP_IND_PNV_SCORE_8_LESS_40",	"CSP_IND_PNV_SCORE_8_LESS_50",	"CSP_IND_PNV_SCORE_8_LESS_60",	"CSP_IND_PNV_SCORE_8_LESS_70",	"CSP_IND_PNV_SCORE_8_LESS_80",	"CSP_IND_PNV_SCORE_8_LESS_90",	"CSP_IND_PNV_SCORE_8_LESS_100",	"CSP_IND_PNV_SCORE_8_LESS_110",	
"CSP_IND_PNV_SCORE_8_LESS_120",	"CSP_IND_PNV_SCORE_8_LESS_130",	"CSP_IND_PNV_SCORE_8_LESS_140",	"CSP_IND_PNV_SCORE_8_MORE_140",	"CSP_IND_AVG_PROB_8_LESS_10",	"CSP_IND_AVG_PROB_8_LESS_25",	"CSP_IND_AVG_PROB_8_LESS_50",	"CSP_IND_AVG_PROB_8_LESS_75",	"CSP_IND_AVG_PROB_8_LESS_90",	"CSP_IND_AVG_PROB_8_LESS_95",	"CSP_IND_AVG_PROB_8_MORE_95",	"CSP_9_CONST",	"CSP_9_REL_RANK",	"CSP_9_RANK_SCR",	
"CSP_9_DEL_RANK",	"CSP_9_AVG_PROB",	"CSP_9_MIN_PROB_1",	"CSP_9_MIN_PROB_2",	"CSP_9_NUM_AA",	"CSP_10_CONST",	"CSP_10_REL_RANK",	"CSP_10_RANK_SCR",	"CSP_10_DEL_RANK",	"CSP_10_AVG_PROB",	"CSP_10_MIN_PROB_1",	"CSP_10_MIN_PROB_2",	"CSP_10_NUM_AA"};



CumulativeSeqProbModel::CumulativeSeqProbModel()
{
	ind_model_was_initialized = false;

	if (num_special_idxs != num_special_ranks ||
		num_rank_level_idxs != num_rank_levels)
	{
		cout << "Error: mismatch in rank levels of CumulativeSeqProbModel!" << endl;
		exit(1);
	}
}

CumulativeSeqProbModel::~CumulativeSeqProbModel()
{
	int i,j,k;

	for (i=0; i<me_models.size(); i++)
		for (j=0; j<me_models[i].size(); j++)
			for (k=0; k<me_models[i][j].size(); k++)
				if (me_models[i][j][k])
					delete me_models[i][j][k];
}



bool CumulativeSeqProbModel::read_model(const Config *config, ifstream& ifs)
{
	const vector< vector<mass_t> >& size_thresholds = config->get_size_thresholds();
	this->me_models.resize(size_thresholds.size());
	int c;


	for (c=1; c<size_thresholds.size(); c++)
	{
		me_models[c].resize(size_thresholds[c].size());
		int s;
		for (s=0; s<size_thresholds[c].size(); s++)
			me_models[c][s].resize(num_special_ranks,NULL);

	}

	char buff[128];
	while (! ifs.eof())
	{
		ifs.getline(buff,128);


		if (strncmp("#MODELS",buff,7))
			continue;	

		int charge=-1,size_idx=-1;
		if (sscanf(buff,"#MODELS\t%d\t%d",&charge,&size_idx) != 2)
		{
			cout << "Error: bad line in CSP model:" << buff << endl;
			exit(1);
		}


		if (charge<0 || size_idx<0)
			continue;

		if (charge>=me_models.size() || size_idx>=me_models[charge].size())
		{
			cout << "Error: illegal model charge and size for cumulative seq prob: " << charge << " " << size_idx << endl;
			exit(1);
		}

		int n;
		for (n=0; n<num_special_ranks; n++)
		{
			me_models[charge][size_idx][n]=new ME_Regression_Model;
			me_models[charge][size_idx][n]->read_regression_model(ifs);
		}
	}

	return true;
}

void compute_avg_min1_min2(const vector<PathPos>& positions, float& avg, float& min1, float& min2)
{
	vector<float> probs;
	probs.clear();

	avg=0;
	min1=0;
	min2=0;

	int i;
	for (i=0; i<positions.size(); i++)
	{
		const float p =positions[i].edge_variant_prob;
		if (p>=0)
		{
			probs.push_back(p);
			avg+=p;
		}
	}


	if (probs.size()<=0)
		return;

	avg /= (float)probs.size();
	sort(probs.begin(),probs.end());
	min1 = probs[0];

	if (probs.size()==1)
		return;

	min2=probs[1];
}





// fills full me samples only for rank levels that equal exactly the special levels (0,1,2,4,8,...)
// other ranks get only a partial set of features (that correspond to the addition after the last full sample
void fill_sparse_me_samples(
							const vector<SeqPath>& solutions,
							const vector<score_pair>& pair_order,
							float spectrum_quality,
							vector<ME_Regression_Sample>& me_samples)
{
	int i;
	
	if (solutions.size()<=0)
		return;

	me_samples.clear();
	me_samples.resize(solutions.size());

	int last_full_sample_idx    = -1;
	int last_special_pos		= -1;
	int last_full_rank_level    = -1;

	const float first_rank_score = solutions[pair_order[0].idx].rerank_score;


	for (i=0; i<pair_order.size(); i++)
	{
		const int sol_idx = pair_order[i].idx;
		const SeqPath& curr_sol = solutions[sol_idx];
		ME_Regression_Sample& curr_me_sample = me_samples[i];

		curr_me_sample.clear();
		if (i==0)
		{
			curr_me_sample.add_feature(CSP_CONST,1.0);
			int j;
			for (j=0; j<num_qual_levels; j++)
				if (spectrum_quality<qual_levels[j])
					break;
			curr_me_sample.add_feature(CSP_IND_QUAL_10+j,1.0);
		}

		const float rank_score = curr_sol.rerank_score;
		const float pnv_score  = curr_sol.path_score;
		const int   num_aa	   = curr_sol.get_num_aa()-6;

		float avg,min1,min2;
		compute_avg_min1_min2(curr_sol.positions,avg,min1,min2);

		// fill a full sample
 		if (last_full_rank_level< num_rank_levels && i==rank_levels[last_full_rank_level+1])
		{
			last_full_rank_level++;

			// copy previous full sample features
			if (last_full_sample_idx>=0)
				curr_me_sample = me_samples[last_full_sample_idx];

			// check if we need to fill special score features
			if (last_special_pos<num_special_ranks &&
				i == special_ranks[last_special_pos+1])
			{
				last_special_pos++;
				const int base_idx =  special_rank_feature_idxs[last_special_pos];

				int j;
				for (j=0; j<num_rank_score_levels; j++)
					if (rank_score<rank_score_levels[j])
						break;

				curr_me_sample.add_feature(base_idx + j,1.0);

				for (j=0; j<num_pnv_score_levels; j++)
					if (pnv_score<pnv_score_levels[j])
						break;

				curr_me_sample.add_feature(base_idx + num_rank_score_levels + j, 1.0);
				
				for (j=0; j<num_avg_prob_levels; j++)
					if (avg<avg_prob_levels[j])
						break;

				curr_me_sample.add_feature(base_idx + num_rank_score_levels + num_pnv_score_levels + j, 1.0);
			}
			
			// add the rank level features
			int feat_idx = rank_level_idxs[last_full_rank_level];
			int del_rank = 0;
			if (i>0)
				del_rank = i-last_full_sample_idx;

			// CSP_0_CONST,	CSP_0_REL_RANK,	CSP_0_RANK_SCR, CSP_0_DEL_RANK,	CSP_0_AVG_PROB,	CSP_0_MIN_PROB_1,	CSP_0_MIN_PROB_2,	CSP_0_NUM_AA,

			curr_me_sample.add_feature(feat_idx++,1.0); 
			curr_me_sample.add_feature(feat_idx++,del_rank);
			curr_me_sample.add_feature(feat_idx++,rank_score);
			curr_me_sample.add_feature(feat_idx++,first_rank_score - rank_score);
			curr_me_sample.add_feature(feat_idx++,avg);
			curr_me_sample.add_feature(feat_idx++,min1);
			curr_me_sample.add_feature(feat_idx++,min2);
			curr_me_sample.add_feature(feat_idx++,num_aa);

			last_full_sample_idx = i;
		}
		else // fill a sparse (delta sample)
		{
			int feat_idx = rank_level_idxs[last_full_rank_level+1];
			
			curr_me_sample.add_feature(feat_idx++,1.0); 
			curr_me_sample.add_feature(feat_idx++,i-last_full_sample_idx);
			curr_me_sample.add_feature(feat_idx++,rank_score);
			curr_me_sample.add_feature(feat_idx++,first_rank_score - rank_score);
			curr_me_sample.add_feature(feat_idx++,avg);
			curr_me_sample.add_feature(feat_idx++,min1);
			curr_me_sample.add_feature(feat_idx++,min2);
			curr_me_sample.add_feature(feat_idx++,num_aa);
		}
	}
}



void make_complete_me_sample(
							const vector<ME_Regression_Sample>& sparse_me_samples,
							int rank_idx,
							ME_Regression_Sample& full_me_sample)
{

	int  level_idx;

	for (level_idx=0; level_idx<num_rank_levels; level_idx++)
	{
		if (rank_levels[level_idx] > rank_idx)
			break;
	}

	level_idx--;
	full_me_sample = sparse_me_samples[rank_levels[level_idx]];

	
	
	if (rank_levels[level_idx] == rank_idx)
		return;

	int j;
	for (j=0; j<sparse_me_samples[rank_idx].f_vals.size(); j++)
		full_me_sample.f_vals.push_back(sparse_me_samples[rank_idx].f_vals[j]);
}


// select what samples to use out of the total triaing set
void select_sample_idxs(int model_idx, int num_solutions, int max_num_per_spectrum, vector<int>& sam_idxs)
{
	
	const int min_rank = special_ranks[model_idx];
	const int max_rank = (model_idx == num_special_ranks -1? POS_INF : special_ranks[model_idx+1]);

	sam_idxs.clear();
	if (max_rank - min_rank-1<=max_num_per_spectrum)
	{
		int i;
		for (i=min_rank; i<max_rank; i++)
			sam_idxs.push_back(i);
		return;
	}

	// if not all samples are to be selected, find how many bins we are looking at

	int first_bin = -1;
	int last_bin  = -1;

	int i;
	for (i=0; i<num_rank_levels; i++)
		if (rank_levels[i] == min_rank)
			break;

	first_bin = i;

	for ( ; i<num_rank_levels; i++)
		if (rank_levels[i] == max_rank)
			break;
	
	last_bin = i;

	if (first_bin == num_rank_levels || last_bin == num_rank_levels)
	{
		cout << "Error: bad bins (mismatch between special ranks and rank levels!" << endl;
		exit(1);
	}

	sam_idxs.push_back(min_rank);
	if (max_rank<num_solutions-1)
	{
		sam_idxs.push_back(max_rank-1);
	}
	else
		sam_idxs.push_back(num_solutions-1);

	const int num_sams_per_bin = (int)( ((float)max_num_per_spectrum/(last_bin - first_bin + 1)+0.5) );
	int bin_idx;
	for (bin_idx = first_bin; bin_idx<=last_bin; bin_idx++)
	{
		int min_sam_idx=(bin_idx>0 ? rank_levels[bin_idx-1] : 0);
		int max_sam_idx=(bin_idx<num_rank_levels-1 ? rank_levels[bin_idx] : num_solutions-1);
		if (min_sam_idx>=num_solutions)
			break;

		int j;
		for (j=0; j<num_sams_per_bin; j++)
		{
			const int rand_sam_idx = (int)(myRandom()*(max_sam_idx-min_sam_idx)+min_sam_idx);
			if (rand_sam_idx<num_solutions-1)
				sam_idxs.push_back(rand_sam_idx);
		}
	}

	sort(sam_idxs.begin(),sam_idxs.end());
}


/*
void CumulativeSeqProbModel::train_seq_prob_models(const FileManager& fm, 
							   void *model_ptr,
							   int tag_length, // 0 - equals full denovo
							   int specific_charge, 
							   int specific_size_idx)
{
	const int num_samples_per_spec = 10;
	const int num_rerank_sols = 2000;
	const int num_models = num_special_ranks;
	vector<int> num_features;
	num_features.resize(num_special_ranks);

	int j;
	for (j=1; j<num_special_ranks; j++)
		num_features[j-1]=special_rank_feature_idxs[j];
	num_features[num_special_ranks-1]=CSP_NUM_FIELDS;

	vector<int> num_tags;
	num_tags.resize(10);
	if (tag_length>0)
		num_tags[tag_length]=500+(tag_length-3)*300;

	AllScoreModels *model = (AllScoreModels *)model_ptr;
	Config *config = model->get_config();

	const vector< vector<mass_t> >& size_threshes =  config->get_size_thresholds();
	
	PeptideRankScorer *rank_model = (PeptideRankScorer *)model->get_rank_model_ptr(1);
	static vector<PrmGraph *> prm_ptrs;
	prm_ptrs.resize(8,NULL);

	// resize model arrays
	
	me_models.resize(size_threshes.size());
	int charge;
	for (charge=1; charge<size_threshes.size(); charge++)
	{
		me_models[charge].resize(size_threshes[charge].size());
		int size_idx;
		for (size_idx=0; size_idx<size_threshes[charge].size(); size_idx++)
			me_models[charge][size_idx].resize(num_models,NULL);
	}
	
	config->set_use_spectrum_charge(1);

	// train models
	for (charge=1; charge<size_threshes.size(); charge++)
	{
		if (specific_charge>0 && charge != specific_charge)
			continue;

		int size_idx;
		for (size_idx=0; size_idx<size_threshes[charge].size(); size_idx++)
		{

			if (specific_size_idx>=0 && size_idx != specific_size_idx)
				continue;

			const mass_t min_mass = (size_idx == 0 ? 0 : size_threshes[charge][size_idx-1]);
			const mass_t max_mass = size_threshes[charge][size_idx];
			FileSet fs;

			fs.select_files_in_mz_range(fm,min_mass/charge,max_mass/charge,charge);
			fs.randomly_reduce_ssfs(10000);

			cout << "Training cumulative sequence probability models for charge " << charge << " size " << size_idx << ", found " << 
				fs.get_total_spectra() << " spectra.." << endl;

			if (! rank_model->get_ind_part_model_was_initialized(charge,size_idx))
			{
				cout << "Error: no re rank model for de novo sequences of these models..." << endl;
				exit(1);
			}
	
			vector<ME_Regression_DataSet> sample_sets;
			sample_sets.resize(num_models);

			BasicSpecReader bsr;
			const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
			int sc;
			for (sc=0; sc<all_ssf.size(); sc++)
			{
				static vector<QCPeak> peaks;
				SingleSpectrumFile *ssf = all_ssf[sc];
				if (peaks.size()<ssf->num_peaks)
				{
					int new_size = ssf->num_peaks*2;
					if (new_size<2500)
						new_size=2500;
					peaks.resize(new_size); 
				}


				cout << sc << "\t" << getRandomSeed() << "\t";
				ssf->print_ssf_stats(config);

	//			cout << sc << "\t(" << all_ssf.size() << ")" << endl;

				if (sc>0 && sc % 100 == 0)
				{
					cout << "Done " << sc << "/" << all_ssf.size() << endl;
				}

				const int num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0]);
				Spectrum s;
				s.init_from_QCPeaks(config,&peaks[0],num_peaks,ssf);

				vector<SeqPath> solutions;
				solutions.clear();

				vector<mass_t> pms_with_19;
				vector<int>    charges;
				pms_with_19.clear();
				charges.clear();		
				BasicSpectrum bs;
				bs.ssf = ssf;
				bs.peaks = &peaks[0];
				bs.num_peaks = num_peaks;

				const string pep_str = bs.ssf->peptide.as_string(config);

				// output m/z and prob values for the different charge states
				vector<PmcSqsChargeRes> pmcsqs_res;
				//model->select_pms_and_charges(config,bs,pms_with_19,charges,&pmcsqs_res);

				const float spectrum_quality = pmcsqs_res[charges[0]].sqs_prob;
				
				if (pms_with_19[0]<100)
				{
					continue;
				}
				
				vector<score_pair> pair_order;
				if (tag_length<1) // full de novo
				{
					generate_denovo_solutions_from_several_pms(
						prm_ptrs,
						model,
						&s,
						true, 
						300,
						7,
						16,
						pms_with_19,
						charges,
						solutions,
						false);
				
					//rank_model->score_denovo_sequences(solutions,ssf,&peaks[0],num_peaks,pair_order,size_idx);
					int j;
					for (j=0; j<solutions.size(); j++)
						solutions[j].rerank_score = pair_order[j].score;

					sort(pair_order.begin(),pair_order.end());
				}
				else // tags
				{
					//generate_tags(prm_ptrs, model, bs, &s, num_tags, tag_length,
					//		  pms_with_19, charges, solutions);

					pair_order.resize(solutions.size());
					int j;
					for (j=0; j<solutions.size(); j++)
						pair_order[j].idx=j;

				}

				int j;
				for (j=0; j<solutions.size(); j++)
				{
					const int sol_idx = (tag_length<1 ? pair_order[j].idx : j);
					if (solutions[sol_idx].check_if_correct(pep_str,config))
						break;
				}

				const int rank_corr = j;

				if (solutions.size()<100)
					continue;

				for (j=0; j<solutions.size(); j++)
				{
					const int sol_idx = (tag_length<1 ? pair_order[j].idx : j);
					solutions[sol_idx].prm_ptr->calc_amino_acid_probs(solutions[sol_idx],j);
				}

		//		cout << sc << "\t" << rank_corr;

				// add training samples to training sets
				vector<ME_Regression_Sample> sparse_samples;
				fill_sparse_me_samples(solutions, pair_order, spectrum_quality, sparse_samples);

	

				int model_idx;
				for (model_idx=0; model_idx<num_models; model_idx++)
				{
					vector<int> sam_idxs;
					select_sample_idxs(model_idx, solutions.size(), 20, sam_idxs);

					int j;
					for (j=0; j<sam_idxs.size(); j++)
					{
						const int sam_idx = sam_idxs[j];
						const int sol_idx = (tag_length<1 ? pair_order[sam_idx].idx : sam_idx);
					
						ME_Regression_Sample full_sam;

						make_complete_me_sample(sparse_samples, sam_idx, full_sam);

						full_sam.label = (sam_idx>=rank_corr ? 0 : 1);
						full_sam.weight=1.0;

						sample_sets[model_idx].add_sample(full_sam);

						int k;
						for (k=0; k<full_sam.f_vals.size(); k++)
						{
							if (full_sam.f_vals[k].f_idx >= num_features[model_idx])
							{
								cout << "MODEL TYPE: " << model_idx << "  sam_idx: " << sam_idx << 
									" num solutions: " << solutions.size() << endl;
								cout << "All sam idxs: ";
								int q;
								for (q=0; q<sam_idxs.size(); q++)
									cout << " " << q << ":" << sam_idxs[q];
								cout << endl << " bad sample:" << endl;
								full_sam.print(Cumulative_Seq_Prob_CSP_names);
								cout << "-------------------";
								cout << "All sparse samples:";
								cout << "--------------------";

								int a;
								for (a=0; a<sparse_samples.size() && a<150; a++)
								{
									sparse_samples[a].print(Cumulative_Seq_Prob_CSP_names);
								}

								exit(1);
							}
						}

					//	cout << "RANK " << j << endl;
					//	full_sam.print(Cumulative_Seq_Prob_CSP_names);
					//	cout << endl;
					}

				//	cout << "\t" << sample_sets[model_idx].num_samples;
				}
			//	cout << endl;
			}

			int type;
			for (type=0; type<sample_sets.size(); type++)
			{
				sample_sets[type].num_classes=2;
				sample_sets[type].num_features = num_features[type];
				if (! me_models[charge][size_idx][type])
					me_models[charge][size_idx][type] = new ME_Regression_Model;

				cout << endl << "Training ME model for type " << type << ", #features = " << num_features[type] << endl;
				sample_sets[type].tally_samples();
				sample_sets[type].print_summary();
				sample_sets[type].print_feature_summary(cout,Cumulative_Seq_Prob_CSP_names);
				me_models[charge][size_idx][type]->train_cg(sample_sets[type],400,1E-5);
				me_models[charge][size_idx][type]->print_ds_probs(sample_sets[type]);

		//	int j;
		//		for (j=0; j<sample_sets[type].num_features; j++)
		//		{
		//			sample_sets[type].report_feature_statistics(j,Cumulative_Seq_Prob_CSP_names[j]);
		//			cout << endl;
		//		}
			}

			// write models
			char fname[128];
			sprintf(fname,"%s_CSP_%d_%d_%d.txt",config->getModelName().c_str(),tag_length, charge,size_idx);
			
			cout << "Writing model: " << fname << endl;

			ofstream ofs(fname);
			ofs << "#MODELS\t" << charge << "\t" << size_idx << endl;
			for (type=0; type<sample_sets.size(); type++)
				me_models[charge][size_idx][type]->write_regression_model(ofs);
			ofs.close();
		}
	}

	ind_model_was_initialized = true; 
}
*/


void CumulativeSeqProbModel::calc_cumulative_seq_probs(
											int first_sol_charge,
											int first_sol_size_idx,
											float spectrum_quality,
											const vector<score_pair>& score_idx_order,
											vector<SeqPath>& solutions) const
{
	// check that ME models are avavilable for this sequence
	bool bad_models=false;
	if (me_models.size()<=first_sol_charge ||
		me_models[first_sol_charge].size()<=first_sol_size_idx)
	{
		bad_models=true;
	}
	else
	{
		int n;
		for (n=0; n<num_special_ranks; n++)
			if (! me_models[first_sol_charge][first_sol_size_idx][n] ||
				! me_models[first_sol_charge][first_sol_size_idx][n]->get_has_weights())
				break;

		if (n<num_special_ranks)
			bad_models=true;
	}

	if (bad_models)
	{
		cout << "Error: not all cumulative seq prob models intialized for charge " << 
			first_sol_charge << " size " << first_sol_size_idx << endl;
			exit(1);
	}


	vector<ME_Regression_Sample> sparse_samples;
	fill_sparse_me_samples(solutions, score_idx_order, spectrum_quality,sparse_samples);

	const int num_solutions = solutions.size();
	int rank=0;
	while (rank<num_solutions)
	{
		const int sam_idx = score_idx_order[rank].idx;
		float prob=0;
		ME_Regression_Sample full_sam;
		make_complete_me_sample(sparse_samples,rank,full_sam);

	//	full_sam.print(Cumulative_Seq_Prob_CSP_names);

		const int sol_charge = solutions[sam_idx].charge;
		
		// only calculate the probability if we are using a sequence with the same charge as
		// the top of the list, otherwise, using the sequence won't improve the cumulative probability
		if (sol_charge == first_sol_charge)
		{
			const int model_idx = get_model_idx(rank);
			prob = me_models[first_sol_charge][first_sol_size_idx][model_idx]->p_y_given_x(0,full_sam);
		}
		
	
		// there might be irregularities when switching between ranks (different models 
		// might be used for different ranks). This fix makes sure the probability is 
		// monotonically increasing
		if (rank>0)
		{
			const float prev_prob = solutions[score_idx_order[rank-1].idx].cumulative_seq_prob;
			if (prev_prob>prob)
				prob=prev_prob;
		}

		solutions[score_idx_order[rank].idx].cumulative_seq_prob=prob;

		// update idx and fill in the skipped probs (not every rank is comupted)
		int prev_rank=rank;
		if (rank<20)
		{
			rank++;
		} 
		else if (rank<50)
		{
			rank+=2;
		}
		else if (rank<100)
		{
			rank+=5;
		}
		else
			rank+=10;

		for (int r=prev_rank+1; r<rank && r<num_solutions; r++)
			solutions[score_idx_order[r].idx].cumulative_seq_prob=prob;
	}
}

