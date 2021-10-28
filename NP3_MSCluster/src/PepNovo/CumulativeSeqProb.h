#ifndef __CUMULATIVESEQPROB_H__
#define __CUMULATIVESEQPROB_H__

#include "ME_REG.h"
#include "Config.h"


struct SeqPath;

const float qual_levels[]		= {0.10,0.25,0.5,0.75,0.9,0.95,0.98,static_cast<float>(POS_INF)};
const int num_qual_levels		= sizeof(qual_levels)/sizeof(float);
const int rank_levels[]			= {0,1,2,4,8,16,32,64,128,256,POS_INF}; // at these ranks we fill the 
const int num_rank_levels		= sizeof(rank_levels)/sizeof(int);
const float rank_score_levels[] = {-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,static_cast<float>(POS_INF)};
const int num_rank_score_levels = sizeof(rank_score_levels)/sizeof(float);
const float pnv_score_levels[]	= {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,static_cast<float>(POS_INF)};
const int num_pnv_score_levels	= sizeof(pnv_score_levels)/sizeof(float);
const float avg_prob_levels[]	= {0.1,0.25,0.5,0.75,0.9,0.95,static_cast<float>(POS_INF)};
const int num_avg_prob_levels	= sizeof(avg_prob_levels)/sizeof(float);
const int special_ranks[]		= {0,8,128};	// at these ranks we fill the big set of rank score and pnv score features
const int num_special_ranks		= sizeof(special_ranks)/sizeof(int);

typedef enum Cumulative_Seq_Prob_CSP {
	CSP_CONST,

	// special rank features (rank 0)
	CSP_IND_QUAL_10,	CSP_IND_QUAL_25,	CSP_IND_QUAL_50,	CSP_IND_QUAL_75,	
	CSP_IND_QUAL_90,	CSP_IND_QUAL_95,	CSP_IND_QUAL_98,	CSP_IND_QUAL_100,

	CSP_IND_RANK_SCORE_TOP_LESS_MINUS_5, CSP_IND_RANK_SCORE_TOP_LESS_MINUS_4, CSP_IND_RANK_SCORE_TOP_LESS_MINUS_3,
	CSP_IND_RANK_SCORE_TOP_LESS_MINUS_2, CSP_IND_RANK_SCORE_TOP_LESS_MINUS_1, CSP_IND_RANK_SCORE_TOP_LESS_0,
	CSP_IND_RANK_SCORE_TOP_LESS_1, CSP_IND_RANK_SCORE_TOP_LESS_2, CSP_IND_RANK_SCORE_TOP_LESS_3, CSP_IND_RANK_SCORE_TOP_LESS_4,
	CSP_IND_RANK_SCORE_TOP_LESS_5, CSP_IND_RANK_SCORE_TOP_LESS_6, CSP_IND_RANK_SCORE_TOP_LESS_7, CSP_IND_RANK_SCORE_MORE_7,

	CSP_IND_PNV_SCORE_TOP_LESS_0,	CSP_IND_PNV_SCORE_TOP_LESS_10,  CSP_IND_PNV_SCORE_TOP_LESS_20,
	CSP_IND_PNV_SCORE_TOP_LESS_30,	CSP_IND_PNV_SCORE_TOP_LESS_40,  CSP_IND_PNV_SCORE_TOP_LESS_50,  CSP_IND_PNV_SCORE_TOP_LESS_60, 
	CSP_IND_PNV_SCORE_TOP_LESS_70,	CSP_IND_PNV_SCORE_TOP_LESS_80,  CSP_IND_PNV_SCORE_TOP_LESS_90,  CSP_IND_PNV_SCORE_TOP_LESS_100, 
	CSP_IND_PNV_SCORE_TOP_LESS_110,	CSP_IND_PNV_SCORE_TOP_LESS_120, CSP_IND_PNV_SCORE_TOP_LESS_130, CSP_IND_PNV_SCORE_TOP_LESS_140, CSP_IND_PNV_SCORE_TOP_MORE_140,

	CSP_IND_AVG_PROB_TOP_LESS_10,	CSP_IND_AVG_PROB_TOP_LESS_25,	CSP_IND_AVG_PROB_TOP_LESS_50,	CSP_IND_AVG_PROB_TOP_LESS_75,	
	CSP_IND_AVG_PROB_TOP_LESS_90,	CSP_IND_AVG_PROB_TOP_LESS_95,	CSP_IND_AVG_PROB_TOP_MORE_95,

	// rank level features
	CSP_0_CONST,	CSP_0_REL_RANK,	CSP_0_RANK_SCR, CSP_0_DEL_RANK,	CSP_0_AVG_PROB,	CSP_0_MIN_PROB_1,	CSP_0_MIN_PROB_2,	CSP_0_NUM_AA,
	CSP_1_CONST,	CSP_1_REL_RANK,	CSP_1_RANK_SCR, CSP_1_DEL_RANK,	CSP_1_AVG_PROB,	CSP_1_MIN_PROB_1,	CSP_1_MIN_PROB_2,	CSP_1_NUM_AA,
	CSP_2_CONST,	CSP_2_REL_RANK,	CSP_2_RANK_SCR, CSP_2_DEL_RANK,	CSP_2_AVG_PROB,	CSP_2_MIN_PROB_1,	CSP_2_MIN_PROB_2,	CSP_2_NUM_AA,
	CSP_3_CONST,	CSP_3_REL_RANK,	CSP_3_RANK_SCR, CSP_3_DEL_RANK,	CSP_3_AVG_PROB,	CSP_3_MIN_PROB_1,	CSP_3_MIN_PROB_2,	CSP_3_NUM_AA,
	CSP_4_CONST,	CSP_4_REL_RANK,	CSP_4_RANK_SCR, CSP_4_DEL_RANK,	CSP_4_AVG_PROB,	CSP_4_MIN_PROB_1,	CSP_4_MIN_PROB_2,	CSP_4_NUM_AA,

	// special rank features (rank 8 = level 4)
	CSP_IND_RANK_SCORE_4_LESS_MINUS_5, CSP_IND_RANK_SCORE_4_LESS_MINUS_4, CSP_IND_RANK_SCORE_4_LESS_MINUS_3,
	CSP_IND_RANK_SCORE_4_LESS_MINUS_2, CSP_IND_RANK_SCORE_4_LESS_MINUS_1, CSP_IND_RANK_SCORE_4_LESS_0,
	CSP_IND_RANK_SCORE_4_LESS_1, CSP_IND_RANK_SCORE_4_LESS_2, CSP_IND_RANK_SCORE_4_LESS_3, CSP_IND_RANK_SCORE_4_LESS_4,
	CSP_IND_RANK_SCORE_4_LESS_5, CSP_IND_RANK_SCORE_4_LESS_6, CSP_IND_RANK_SCORE_4_LESS_7, CSP_IND_RANK_SCORE_4_MORE_7,

	CSP_IND_PNV_SCORE_4_LESS_0,		CSP_IND_PNV_SCORE_4_LESS_10,  CSP_IND_PNV_SCORE_4_LESS_20,
	CSP_IND_PNV_SCORE_4_LESS_30,	CSP_IND_PNV_SCORE_4_LESS_40,  CSP_IND_PNV_SCORE_4_LESS_50,  CSP_IND_PNV_SCORE_4_LESS_60, 
	CSP_IND_PNV_SCORE_4_LESS_70,	CSP_IND_PNV_SCORE_4_LESS_80,  CSP_IND_PNV_SCORE_4_LESS_90,  CSP_IND_PNV_SCORE_4_LESS_100, 
	CSP_IND_PNV_SCORE_4_LESS_110,	CSP_IND_PNV_SCORE_4_LESS_120, CSP_IND_PNV_SCORE_4_LESS_130, CSP_IND_PNV_SCORE_4_LESS_140, CSP_IND_PNV_SCORE_4_MORE_140,

	CSP_IND_AVG_PROB_4_LESS_10,	CSP_IND_AVG_PROB_4_LESS_25,	CSP_IND_AVG_PROB_4_LESS_50,	CSP_IND_AVG_PROB_4_LESS_75,	
	CSP_IND_AVG_PROB_4_LESS_90,	CSP_IND_AVG_PROB_4_LESS_95,	CSP_IND_AVG_PROB_4_MORE_95,

	CSP_5_CONST,	CSP_5_REL_RANK,	CSP_5_RANK_SCR, CSP_5_DEL_RANK,	CSP_5_AVG_PROB,	CSP_5_MIN_PROB_1,	CSP_5_MIN_PROB_2,	CSP_5_NUM_AA,
	CSP_6_CONST,	CSP_6_REL_RANK,	CSP_6_RANK_SCR, CSP_6_DEL_RANK,	CSP_6_AVG_PROB,	CSP_6_MIN_PROB_1,	CSP_6_MIN_PROB_2,	CSP_6_NUM_AA,
	CSP_7_CONST,	CSP_7_REL_RANK,	CSP_7_RANK_SCR, CSP_7_DEL_RANK,	CSP_7_AVG_PROB,	CSP_7_MIN_PROB_1,	CSP_7_MIN_PROB_2,	CSP_7_NUM_AA,
	CSP_8_CONST,	CSP_8_REL_RANK,	CSP_8_RANK_SCR, CSP_8_DEL_RANK,	CSP_8_AVG_PROB,	CSP_8_MIN_PROB_1,	CSP_8_MIN_PROB_2,	CSP_8_NUM_AA,

	// special rank features (rank 128 = level 8)
	CSP_IND_RANK_SCORE_8_LESS_MINUS_5, CSP_IND_RANK_SCORE_8_LESS_MINUS_4, CSP_IND_RANK_SCORE_8_LESS_MINUS_3,
	CSP_IND_RANK_SCORE_8_LESS_MINUS_2, CSP_IND_RANK_SCORE_8_LESS_MINUS_1, CSP_IND_RANK_SCORE_8_LESS_0,
	CSP_IND_RANK_SCORE_8_LESS_1, CSP_IND_RANK_SCORE_8_LESS_2, CSP_IND_RANK_SCORE_8_LESS_3, CSP_IND_RANK_SCORE_8_LESS_4,
	CSP_IND_RANK_SCORE_8_LESS_5, CSP_IND_RANK_SCORE_8_LESS_6, CSP_IND_RANK_SCORE_8_LESS_7, CSP_IND_RANK_SCORE_8_MORE_7,

	CSP_IND_PNV_SCORE_8_LESS_0,		CSP_IND_PNV_SCORE_8_LESS_10,  CSP_IND_PNV_SCORE_8_LESS_20,
	CSP_IND_PNV_SCORE_8_LESS_30,	CSP_IND_PNV_SCORE_8_LESS_40,  CSP_IND_PNV_SCORE_8_LESS_50,  CSP_IND_PNV_SCORE_8_LESS_60, 
	CSP_IND_PNV_SCORE_8_LESS_70,	CSP_IND_PNV_SCORE_8_LESS_80,  CSP_IND_PNV_SCORE_8_LESS_90,  CSP_IND_PNV_SCORE_8_LESS_100, 
	CSP_IND_PNV_SCORE_8_LESS_110,	CSP_IND_PNV_SCORE_8_LESS_120, CSP_IND_PNV_SCORE_8_LESS_130, CSP_IND_PNV_SCORE_8_LESS_140, CSP_IND_PNV_SCORE_8_MORE_140,

	CSP_IND_AVG_PROB_8_LESS_10,	CSP_IND_AVG_PROB_8_LESS_25,	CSP_IND_AVG_PROB_8_LESS_50,	CSP_IND_AVG_PROB_8_LESS_75,	
	CSP_IND_AVG_PROB_8_LESS_90,	CSP_IND_AVG_PROB_8_LESS_95,	CSP_IND_AVG_PROB_8_MORE_95,


	CSP_9_CONST,	CSP_9_REL_RANK,	CSP_9_RANK_SCR, CSP_9_DEL_RANK,	CSP_9_AVG_PROB,	CSP_9_MIN_PROB_1,	CSP_9_MIN_PROB_2,	CSP_9_NUM_AA,
	CSP_10_CONST,	CSP_10_REL_RANK, CSP_10_RANK_SCR, CSP_10_DEL_RANK, CSP_10_AVG_PROB, CSP_10_MIN_PROB_1,	CSP_10_MIN_PROB_2,	CSP_10_NUM_AA,

	CSP_NUM_FIELDS
} Cumulative_Seq_Prob_CSP;


const int special_rank_feature_idxs[]={CSP_IND_RANK_SCORE_TOP_LESS_MINUS_5,
									   CSP_IND_RANK_SCORE_4_LESS_MINUS_5,
									   CSP_IND_RANK_SCORE_8_LESS_MINUS_5};



const int rank_level_idxs[]={CSP_0_CONST, CSP_1_CONST, CSP_2_CONST, CSP_3_CONST,
							 CSP_4_CONST, CSP_5_CONST, CSP_6_CONST, CSP_7_CONST,
							 CSP_8_CONST, CSP_9_CONST, CSP_10_CONST};





class CumulativeSeqProbModel {
public:

	CumulativeSeqProbModel();
	~CumulativeSeqProbModel();

	bool read_model(const Config *config, ifstream& ifs);

	bool get_ind_initialized() const { return ind_model_was_initialized; }

	int get_model_idx(int rank) const
	{
		int j;
		for (j=1; j<num_special_ranks; j++)
			if (special_ranks[j]>rank)
				break;
		return j-1;
	}

/*	void train_seq_prob_models(const FileManager& fm, 
							   void *model_ptr,
							   int tag_length, // 0 equals full denovo
							   int specific_charge=-1, 
							   int specific_size_idx=-1);*/

	void read_seq_prob_models(Config *config, char *file_name);

	void write_seq_prob_models(const char *path) const;

	void calc_cumulative_seq_probs(
						int first_sol_charge,
						int first_sol_size_idx,
						float spectrum_quality,
						const vector<score_pair>& score_idx_order,
						vector<SeqPath>& solutions) const;


private:
	bool ind_model_was_initialized;
	vector< vector< vector< ME_Regression_Model *> > > me_models;

};

void fill_sparse_me_samples(const vector<SeqPath>& solutions,
							const vector<score_pair>& pair_order,
							float spectrum_quality,
							vector<ME_Regression_Sample>& me_samples);

void make_complete_me_sample(const vector<ME_Regression_Sample>& sparse_me_samples,
							 int sample_idx,
							 ME_Regression_Sample& full_me_sample);



#endif


