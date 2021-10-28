#include "PeptideRankScorer.h"
#include "PrmGraph.h"
#include "AllScoreModels.h"
#include "PepNovo_auxfun.h"

extern const float RKH_pair_matrix[6][6];

void find_prediction_ranks(const vector<float>& scores, vector<int>& ranks)
{
	vector<score_pair> pairs;
	pairs.resize(scores.size());
	ranks.clear();

	if (scores.size()==0)
		return;
	
	pairs[0].idx=0;
	pairs[0].score=NEG_INF;
	int i;
	for (i=1; i<scores.size(); i++)
	{
		pairs[i].idx=i;
		pairs[i].score = scores[i];
	}
	sort(pairs.begin(),pairs.end());
	ranks.resize(scores.size());
	for (i=0; i<pairs.size(); i++)
		if (pairs[i].score == NEG_INF)
		{
			ranks[pairs[i].idx]=POS_INF;
		}
		else
			ranks[pairs[i].idx]=i;

/*	cout << setprecision(1) << fixed;
	for (i=0; i<scores.size(); i++)
		cout << scores[i] << "\t";
	cout << endl;
	for (i=0; i<scores.size(); i++)
		cout << ranks[i] << "\t";
	cout << endl << endl;*/
}

struct inten_pair {
	inten_pair() : idx(int(NEG_INF)), inten(NEG_INF) {};
	inten_pair(int _i, intensity_t _n) : idx(_i), inten(_n) {};
	bool operator< (const inten_pair& other) const
	{
		return inten>other.inten;
	}
	int idx;
	intensity_t inten;
};

void find_inten_ranks(const vector<intensity_t>& intens, vector<int>& ranks)
{
	vector<inten_pair> pairs;
	pairs.resize(intens.size());
	ranks.clear();
	if (intens.size()==0)
		return;

	pairs[0].idx=0;
	pairs[0].inten=NEG_INF;

	int i;
	for (i=1; i<intens.size(); i++)
	{
		pairs[i].idx=i;
		pairs[i].inten = intens[i];
	}
	sort(pairs.begin(),pairs.end());
	ranks.resize(intens.size());
	for (i=0; i<intens.size(); i++)
		if (pairs[i].inten <= 0)
		{
			ranks[pairs[i].idx]=POS_INF;
		}
		else
			ranks[pairs[i].idx]=i;
}


struct rank_pair {
	bool operator< ( const rank_pair& other) const
	{
		return org_rank<other.org_rank;
	}

	int idx;
	int org_rank;
};

void convert_to_relative_ranks(const vector<int>& org_ranks, vector<int>& rel_ranks)
{
	const int num_ranks = org_ranks.size();
	vector<rank_pair> pairs;
	pairs.resize(num_ranks);
	if (num_ranks==0)
		return;

	int i;
	for (i=0; i<num_ranks; i++)
	{
		pairs[i].idx=i;
		pairs[i].org_rank = org_ranks[i];
	}
	sort(pairs.begin(),pairs.end());

	rel_ranks.resize(num_ranks);
	for (i=0; i<num_ranks; i++)
		if (pairs[i].org_rank >= 999)
		{
			rel_ranks[pairs[i].idx]=POS_INF;
		}
		else
			rel_ranks[pairs[i].idx]=i;
}

/***********************************************************************
dimensions of intens are (#all frags X #breakages (== peptide length +1)
************************************************************************/
void DeNovoPartitionModel::fill_combined_peak_prediction_features(
			const PeptideSolution& sol,
			const vector< vector<intensity_t> >& intens,
			const PeakRankModel *peak_model,
			RankBoostSample& rbs,
			int specific_size) const
{
	const int num_ranks_to_consider = 46;
	int f_idx = combined_ppp_start_idx;

	PeptidePeakPrediction ppp;
	peak_model->calc_peptide_predicted_scores(sol, ppp, specific_size);

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

	vector<int> pred2obs, obs2pred;
	vector<int> num_obs_for_frag, num_pred_for_frag;
	vector<float> ordered_scores,  // scores sorted according to their value
				  obs_ordered_scores;  // scores sorted according to the observed intensity rank
	pred2obs.resize(num_ranks_to_consider,999); // look at top 50 peaks
	obs2pred.resize(num_ranks_to_consider,999);
	ordered_scores.resize(num_ranks_to_consider,NEG_INF);
	obs_ordered_scores.resize(num_ranks_to_consider,NEG_INF);
	num_obs_for_frag.resize(num_frags,0);
	num_pred_for_frag.resize(num_frags,0);


	for (f=0; f<num_frags; f++)
	{
		if (observed_ranks[f].size() != predicted_ranks[f].size())
		{
			cout << "#obs  frags: " << observed_ranks.size() << endl;
			cout << "#pred frags: " << predicted_ranks.size() << endl;
			cout << "Error: mismatch in rank dimensionalities!" << endl;
			cout << f << "\tobs : " << observed_ranks[f].size() << "   pred " << predicted_ranks[f].size() << endl;
			exit(1);
		}
		const int num_ranks = predicted_ranks[f].size();
		const vector<float>& frag_rank_scores = ppp.rank_scores[f];
		const vector<float>& frag_intens = observed_intens[f];
		int i;
		for (i=0; i<num_ranks; i++)
		{
			const int obs_rank  = observed_ranks[f][i];
			const int pred_rank = predicted_ranks[f][i];
			const float pred_score = frag_rank_scores[i];

			if (pred_rank<num_ranks_to_consider)
			{
				pred2obs[pred_rank]=obs_rank;
				ordered_scores[pred_rank]=pred_score;
			}

			if (obs_rank<num_ranks_to_consider)
			{
				obs2pred[obs_rank]=pred_rank;
				obs_ordered_scores[obs_rank]=pred_score;
			}

			if (frag_intens[i]>0)
				num_obs_for_frag[f]++;

			if (frag_rank_scores[i]>NEG_INF)
				num_pred_for_frag[f]++;
		}
	}

	vector<int> missed_ranks;
	for (i=0; i<num_ranks_to_consider; i++)
		if (pred2obs[i]>900)
			missed_ranks.push_back(i);
	
	const int mobility = sol.calc_mobility();
	rbs.add_real_feature(f_idx++,(float)mobility);
	
	// fill per frag features
	for (f=0; f<num_frags && f<3; f++)
	{
		vector<int> rel_obs_ranks, rel_pred_ranks;
		convert_to_relative_ranks(observed_ranks[f], rel_obs_ranks);
		convert_to_relative_ranks(predicted_ranks[f],rel_pred_ranks);

		if (num_pred_for_frag[f]>0)
		{
			rbs.add_real_feature(f_idx+f,num_obs_for_frag[f]/(float)num_pred_for_frag[f]);
		}
	}
	f_idx+=3;

	//sums of missed ranks
	rbs.add_real_feature(f_idx++,(float)missed_ranks.size());

	int sum_missed=NEG_INF;
	if (missed_ranks.size()>4)
	{
		sum_missed = missed_ranks[0]+missed_ranks[1]+missed_ranks[2]+missed_ranks[3]+missed_ranks[4];
		rbs.add_real_feature(f_idx+(mobility-1),sum_missed);
		rbs.add_real_feature(f_idx+3,sum_missed);
	}
	f_idx+=4;

	if (missed_ranks.size()>9)
	{
		sum_missed=missed_ranks[5]+missed_ranks[6]+missed_ranks[7]+missed_ranks[8]+missed_ranks[9];
		rbs.add_real_feature(f_idx+(mobility-1),sum_missed);
		rbs.add_real_feature(f_idx+3,sum_missed);
	}
	f_idx+=4;

	if (missed_ranks.size()>14)
	{
		sum_missed=missed_ranks[10]+missed_ranks[11]+missed_ranks[12]+missed_ranks[13]+missed_ranks[14];
		rbs.add_real_feature(f_idx+(mobility-1),sum_missed);
		rbs.add_real_feature(f_idx+3,sum_missed);
	}
	f_idx+=4;

	for (i=0; i<7; i++)
		rbs.add_real_feature(f_idx++,pred2obs[i]);

	for (i=0; i<7; i++)
		rbs.add_real_feature(f_idx++,obs2pred[i]);

	for (i=0; i<9&& i*2<missed_ranks.size(); i++)
		rbs.add_real_feature(f_idx+i,missed_ranks[2*i]);
	f_idx+=9;

	for (i=0; i<7; i++)
		if (ordered_scores[i]>NEG_INF && obs_ordered_scores[i]>NEG_INF)
			rbs.add_real_feature(f_idx+i,ordered_scores[i]-obs_ordered_scores[i]);
	f_idx+=7;
		

	// calc dot prod score feature
	// normalize ranks according to f(x)=1/(1+x)
	static vector<float> one_over_rank, one_over_rank_sqr;
	if (one_over_rank.size()<1000)
	{
		one_over_rank.resize(1000);
		one_over_rank_sqr.resize(1000);
		int i;
		for (i=0; i<num_ranks_to_consider; i++)
		{
			one_over_rank[i]=(i/(1.0+(float)i));
			one_over_rank_sqr[i]=one_over_rank[i]*one_over_rank[i];
		}
	}

	float top_a=0, top_b=0;
	float bottom_a1=1E-8,bottom_a2=1E-8, bottom_b1=1E-8, bottom_b2=1E-8;
	int round_idx=0;
	for (i=0; i<num_ranks_to_consider && i<45; i++)
	{
		const int obs_rank=(pred2obs[i]>999 ? 999 : pred2obs[i]);
		const int pred_rank = (obs2pred[i]>999 ? 999 : obs2pred[i]);

		top_a += (one_over_rank[i]*one_over_rank[obs_rank]);
		bottom_a1 += one_over_rank_sqr[i];
		bottom_a2 += one_over_rank_sqr[obs_rank];

		top_b += (one_over_rank[i]*one_over_rank[pred_rank]);
		bottom_b1 += one_over_rank_sqr[i];
		bottom_b2 += one_over_rank_sqr[pred_rank];

		if (i>0 && (i+1) % 15 == 0)
		{
			rbs.add_real_feature(f_idx+2*round_idx,top_a/sqrt(bottom_a1*bottom_a2));
			rbs.add_real_feature(f_idx+2*round_idx+1,top_b/sqrt(bottom_b1*bottom_b2));
			round_idx++;
			if (round_idx==3)
				break;
		}
	}
	f_idx+=6;	
}




void DeNovoPartitionModel::fill_peak_prediction_features(
							const PeptideSolution& sol,
							const vector< vector<intensity_t> >& intens,
							const PeakRankModel *peak_model,
							RankBoostSample& rbs,
							int specific_size) const
{
	int f_idx = ppp_start_idx;

	vector< vector<int> > ppp_prediction_ranks;
	if (num_ppp_frags>1)
		ppp_prediction_ranks.resize(num_ppp_frags);
	
	PeptidePeakPrediction ppp;
	peak_model->calc_peptide_predicted_scores(sol, ppp, specific_size, &ppp_frag_type_idxs);

	int f;
	for (f=0; f<num_ppp_frags; f++)
	{
		const int frag_idx = ppp_frag_type_idxs[f];
		const vector<intensity_t>& frag_intens = intens[frag_idx];
		const vector<float>& predicted_scores  = ppp.rank_scores[frag_idx];
		int i;

		int num_obs_frags=0;
		for (i=1; i<frag_intens.size(); i++)
			if (frag_intens[i]>0)
				num_obs_frags++;

		int num_predicted=0;
		for (i=1; i<predicted_scores.size(); i++)
			if (predicted_scores[i]>NEG_INF)
				num_predicted++;
		///
		rbs.add_real_feature(f_idx++,(float)num_obs_frags);
		rbs.add_real_feature(f_idx++,(float)num_predicted);
		if (num_predicted>0)
			rbs.add_real_feature(f_idx,(float)num_obs_frags/(float)num_predicted);
		f_idx++;

	//	for (i=0; i<frag_intens.size(); i++)
	//		cout << fixed << setprecision(1) << frag_intens[i] << "\t";
	//	cout << endl;
		
		vector<int> prediction_ranks, inten_ranks;
		find_prediction_ranks(predicted_scores,prediction_ranks);
		find_inten_ranks(frag_intens, inten_ranks);

		if (num_ppp_frags>1)
			ppp_prediction_ranks[f]=prediction_ranks;
		
		vector<int> max_ranks;
		max_ranks.push_back(1);
		max_ranks.push_back(3);
		max_ranks.push_back(5);
		max_ranks.push_back(7);
		max_ranks.push_back((int)(0.1666*num_predicted));
		max_ranks.push_back((int)(0.3333*num_predicted));
		max_ranks.push_back((int)(0.5000*num_predicted));
		max_ranks.push_back((int)(0.6667*num_predicted));

		vector<int> counts;
		counts.resize(max_ranks.size(),0);

		for (i=0; i<predicted_scores.size(); i++)
			if (predicted_scores[i]>NEG_INF)
			{
				int j;
				for (j=0; j<counts.size(); j++)
					if (prediction_ranks[i]>=0 && 
						prediction_ranks[i]<max_ranks[j] && 
						inten_ranks[i]<max_ranks[j])
					counts[j]++;
			}

		for (i=0; i<4; i++)
		{
			if (max_ranks[i]>0 && num_obs_frags>=max_ranks[i])
				rbs.add_real_feature(f_idx,(float)counts[i]);
			f_idx++;
		}

		for (   ; i<counts.size(); i++)
		{
			if (max_ranks[i]>0 && num_obs_frags>=max_ranks[i])
				rbs.add_real_feature(f_idx,(float)counts[i]/(float)max_ranks[i]);
			f_idx++;
		}
	

		// find the prediction ranks for the hi
		vector<int> missing_ranks;
		missing_ranks.resize(3,prediction_ranks.size());
		int num_miss=0;
		int pred_rank;
		for (pred_rank=0; pred_rank<prediction_ranks.size(); pred_rank++)
		{
			int cut_idx;
			for (cut_idx=1; cut_idx<prediction_ranks.size(); cut_idx++)
				if (prediction_ranks[cut_idx]==pred_rank)
					break;
			
			if (cut_idx==prediction_ranks.size())
				break;

			if (frag_intens[cut_idx]<=0)
				missing_ranks[num_miss++]=pred_rank;
			if (num_miss>=3)
				break;
		}

		rbs.add_real_feature(f_idx++,(float)missing_ranks[0]);
		rbs.add_real_feature(f_idx++,(float)missing_ranks[1]);
		rbs.add_real_feature(f_idx++,(float)missing_ranks[2]);
		rbs.add_real_feature(f_idx++,(float)(missing_ranks[0]+missing_ranks[1]));
		rbs.add_real_feature(f_idx++,(float)(missing_ranks[0]+missing_ranks[1]+missing_ranks[2]));

		// calculate the prediction score offsets
		// (equals ths sum of scores for peaks that are between the peak's true
		// position and the peaks predicted position)
		vector<float> sorted_predicted_scores = predicted_scores;
		sort(sorted_predicted_scores.begin(),sorted_predicted_scores.end());

		const int max_idx = sorted_predicted_scores.size()-1;
		const int mid = sorted_predicted_scores.size()/2;
		for (i=0; i<mid; i++)
		{
			float t=sorted_predicted_scores[i];
			sorted_predicted_scores[i]=sorted_predicted_scores[max_idx-i];
			sorted_predicted_scores[max_idx-i]=t;
		}

		vector<float> score_offsets;
		score_offsets.resize(10,NEG_INF);
		for (pred_rank=0; pred_rank<10; pred_rank++)
		{
			int cut_idx;
			for (cut_idx=1; cut_idx<prediction_ranks.size(); cut_idx++)
				if (prediction_ranks[cut_idx]==pred_rank)
					break;
			
			if (cut_idx==prediction_ranks.size())
				break;

			const int obs_rank = inten_ranks[cut_idx];
			if (frag_intens[cut_idx]>0)
				score_offsets[pred_rank]=predicted_scores[cut_idx]-sorted_predicted_scores[obs_rank];

		//	cout << pred_rank << "\t" << cut_idx << "\t" << score_offsets[pred_rank] << endl;
		}

		for (i=0; i<10 && i<num_obs_frags; i++)
			rbs.add_real_feature(f_idx+i,score_offsets[i]);
	
		f_idx+=10;



		
	}


//	ppp.print_ranks_vs_intens(intens);

	vector<score_pair> pairs;
	if (num_ppp_frags>1)
	{
		const int frag1_idx = ppp_frag_type_idxs[0];
		const int frag2_idx = ppp_frag_type_idxs[1];
		int i;
		int num_comp_pairs=0;
		const int min_size = (intens[frag1_idx].size()>intens[frag2_idx].size() ? intens[frag2_idx].size() : 
									intens[frag1_idx].size());
		for (i=1; i<min_size; i++)
		{
			const float inten1 = intens[frag1_idx][i];
			const float inten2 = intens[frag2_idx][i];
			const int pred1 = ppp_prediction_ranks[0][i];
			const int pred2 = ppp_prediction_ranks[1][i];

			if (inten1>0 && inten2>0)
				num_comp_pairs++;

			if (inten1>=0 && inten2>=0 && pred1<POS_INF && pred2<POS_INF)
				pairs.push_back(score_pair(i,pred1+pred2));
		}

		rbs.add_real_feature(f_idx, (float)num_comp_pairs);

		sort(pairs.begin(),pairs.end());
		
		int counter=0;
		for (i=pairs.size()-1; i>=0; i--)
		{
			if (counter++==7)
				break;

			const int cut_idx = pairs[i].idx;
			int stat=0;
			if (intens[frag1_idx][cut_idx]>0)
				stat++;
			if (intens[frag2_idx][cut_idx]>0)
				stat++;

			rbs.add_real_feature(f_idx+counter, (float)stat);
		}
	}
	f_idx += 8;
}



void DeNovoPartitionModel::fill_prm_features(const PeptideSolution& sol, 
											 const SeqPath& path, 
											 int model_type, 
											 RankBoostSample& rbs) const
{
	int f_idx = prm_start_idx;

	const PrmGraph *prm = path.prm_ptr;
	const Config *config = prm->get_config();
	const vector<PathPos>& positions = path.positions;
	const int num_aas = path.get_num_aa();
	vector<int> amino_acids;
	path.get_amino_acids(amino_acids);
	if (! prm)
	{
		cout << "Error: SeqPath has not prm ptr!" << endl;
		exit(1);
	}

	int i;
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	mass_t pep_mass = 0;
	for (i=0; i<num_aas; i++)
		pep_mass += aa2mass[amino_acids[i]];

	const mass_t delta = (positions[positions.size()-1].mass - positions[0].mass) - pep_mass;

	if (fabs(delta)>20)
	{
		int i;
		for (i=0; i<positions.size(); i++)
		{
			cout << "Error: mismatch between peptide and Seqpath: " << i << "\t" << 
				positions[i].node_idx << " " << positions[i].mass << endl;
		}
		exit(0);
	}

	vector<float> breakage_scores;
	breakage_scores.clear();
	for (i=0; i<positions.size(); i++)
	{
		const int node_idx = positions[i].node_idx;
		if (node_idx>=0)
		{
			const Node& node = prm->get_node(node_idx);
			if (node.type != NODE_N_TERM && node.type != NODE_C_TERM)
				breakage_scores.push_back(positions[i].breakage->score);
		}
	}

	float min_consec_two_breaks = POS_INF;
	float min_consec_three_breaks = POS_INF;
	float max_consec_two_breaks = NEG_INF;
	float max_consec_three_breaks = NEG_INF;

	float total_breakage_score = NEG_INF;
	if (breakage_scores.size()>1)
	{
		total_breakage_score=0;
		for (i=0; i<breakage_scores.size()-1; i++)
		{
			const float sum_score = breakage_scores[i]+breakage_scores[i+1];
			if (sum_score<min_consec_two_breaks)
				min_consec_two_breaks = sum_score;
			if (sum_score>max_consec_two_breaks)
				max_consec_two_breaks = sum_score;

			if (i>0)
			{
				const float sum_three= sum_score + breakage_scores[i-1];
				if (sum_three<min_consec_three_breaks)
					min_consec_three_breaks = sum_three;
				if (sum_three>max_consec_three_breaks)
					max_consec_three_breaks = sum_three;
			}
			total_breakage_score+=breakage_scores[i];
		}
		total_breakage_score+=breakage_scores[i];
	}


		


	int aa_start=0;
	if (positions[0].node_idx==0 && positions[0].edge_idx<0)
		while (aa_start<positions.size() && positions[aa_start].node_idx<0)
			aa_start++;
	int aa_end = positions.size()-2;
	while (aa_end>aa_start && positions[aa_end].node_idx<0)
		aa_end--;

	int eff_num_aas = (aa_end-aa_start);
	if (sol.reaches_c_terminal)
		eff_num_aas++;

	if (eff_num_aas<6)
		eff_num_aas=6;

	const int num_breakage_scores = breakage_scores.size();
	const int delta_num_breakages = eff_num_aas-breakage_scores.size();

//	int num_missing_edges=0;
//	for (i=0; i<positions.size()-1; i++)
//		if (positions[i].node_idx>=0 && positions[i].edge_idx<0)
//			num_missing_edges++;

	sort(breakage_scores.begin(),breakage_scores.end());
	
	// add features

	int combo_idx=0;
	if (sol.reaches_n_terminal && ! sol.reaches_c_terminal)
	{
		combo_idx=1;
	}
	else if (! sol.reaches_n_terminal && sol.reaches_c_terminal)
	{
		combo_idx=2;
	}
	else if (! sol.reaches_n_terminal && ! sol.reaches_c_terminal)
		combo_idx = 3;

	const int idx_shift = ( (model_type == 1 || model_type == 3) ? combo_idx * 6 : 0);
	int sfidx = f_idx + idx_shift;

	rbs.add_real_feature(sfidx,   delta);

	if (total_breakage_score>NEG_INF)
		rbs.add_real_feature(sfidx+1, total_breakage_score);

	if (num_breakage_scores>0 && total_breakage_score>NEG_INF)
		rbs.add_real_feature(sfidx+2, total_breakage_score/num_breakage_scores);

	if (eff_num_aas>0 && total_breakage_score>NEG_INF)
		rbs.add_real_feature(sfidx+3, total_breakage_score/eff_num_aas);
	
	rbs.add_real_feature(sfidx+4, path.path_score);

	if (eff_num_aas>0)
		rbs.add_real_feature(sfidx+5, path.path_score/eff_num_aas);
	
	f_idx+= ((model_type == 1 || model_type == 3) ? 24 : 6);

	if (model_type == 1 || model_type == 3)
	{
		int path_rank = path.org_rank;
		if (model_type == 1 && path_rank > 200)
			path_rank = 200;
		if (model_type == 3 && path_rank> 75)
			path_rank = 75;

		rbs.add_real_feature(f_idx++, total_breakage_score);
		rbs.add_real_feature(f_idx++, path.path_score);
		rbs.add_real_feature(f_idx++, path_rank);
		rbs.add_real_feature(f_idx++, path.multi_path_score);

		if (path.delta_score>=0)
		{
			rbs.add_real_feature(f_idx,path.delta_score);
			if (path.delta_score<=1.5)
			{
				rbs.add_real_feature(f_idx+1,path_rank);
			}
			else if (path.delta_score<=7.5)
			{
				rbs.add_real_feature(f_idx+2,path_rank);
			}
			else if (path.delta_score<=15.0)
			{
				rbs.add_real_feature(f_idx+3,path_rank);
			}
			else
				rbs.add_real_feature(f_idx+4,path_rank);
				
		}
		f_idx+=5;	


		if (model_type == 3)
		{
			
			if (path.tag_percent_top_5>0)
				rbs.add_real_feature(f_idx,path.tag_percent_top_5);

			if (path.tag_percent_top_20>0)
				rbs.add_real_feature(f_idx+1,path.tag_percent_top_20);

			if (path.tag_percent_all>0)
				rbs.add_real_feature(f_idx+2,path.tag_percent_all);

			if (path.tag_percent_top_5>0)
			{
				rbs.add_real_feature(f_idx+3,path.org_rank);
			}
			else if (path.tag_percent_top_20>0)
			{
				rbs.add_real_feature(f_idx+4,path.org_rank);
			}
			else if (path.tag_percent_all>0)
			{
				rbs.add_real_feature(f_idx+5,path.org_rank);
			}

			if (path.multi_path_rank<POS_INF)
				rbs.add_real_feature(f_idx+5,path.multi_path_rank);
		
			f_idx+=7;
		}
	}
	
	rbs.add_real_feature(f_idx++,delta_num_breakages);
	rbs.add_real_feature(f_idx++,path.num_forbidden_nodes);
	rbs.add_real_feature(f_idx++,num_breakage_scores);
	if (num_breakage_scores>0)
		rbs.add_real_feature(f_idx,breakage_scores[0]);
	f_idx++;

	if (num_breakage_scores>2)
		rbs.add_real_feature(f_idx,breakage_scores[1]);
	f_idx++;

	if (num_breakage_scores>4)
		rbs.add_real_feature(f_idx,breakage_scores[2]);
	f_idx++;

	if (min_consec_three_breaks<POS_INF)
		rbs.add_real_feature(f_idx,min_consec_three_breaks);
	f_idx++;
	if (max_consec_three_breaks>NEG_INF)
		rbs.add_real_feature(f_idx,max_consec_three_breaks);
	f_idx++;

	if (min_consec_two_breaks<POS_INF)
		rbs.add_real_feature(f_idx,min_consec_two_breaks);
	f_idx++;
	if (max_consec_two_breaks>NEG_INF)
		rbs.add_real_feature(f_idx,max_consec_two_breaks);
	f_idx++;


	int num1=0,num2=0,num3=0,num4=0,num5=0;

	for (i=0; i<breakage_scores.size(); i++)
	{
		const float& score = breakage_scores[i];
		if (score<-10)
		{
			num1++;
		}
		else if (score<0)
		{
			num2++;
		}
		else if (score<8)
		{
			num3++;
		}
		else if (score<15)
		{
			num4++;
		}
		else
			num5++;
	}	

	rbs.add_real_feature(f_idx++,num1);
	rbs.add_real_feature(f_idx++,num2);
	rbs.add_real_feature(f_idx++,num3);
	rbs.add_real_feature(f_idx++,num4);
	rbs.add_real_feature(f_idx++,num5);
	if (num_breakage_scores>0)
	{
		rbs.add_real_feature(f_idx++,(float)num1/(float)num_breakage_scores);
		rbs.add_real_feature(f_idx++,(float)(num1+num2)/(float)num_breakage_scores);
		rbs.add_real_feature(f_idx++,(float)(num3+num4+num5)/(float)num_breakage_scores);
		rbs.add_real_feature(f_idx++,(float)(num4+num5)/(float)num_breakage_scores);
	}
	else
		f_idx+=4;


	if (sol.reaches_n_terminal && positions[0].edge_idx>=0 && num_breakage_scores>2)
	{
		int i=1;
		while (positions[i].node_idx<0)
			i++;
		rbs.add_real_feature(f_idx,positions[i].node_score);
	}
	f_idx++;

	if (sol.reaches_c_terminal && num_breakage_scores>2)
	{
		int i=positions.size()-2;
		while (positions[i].node_idx<0)
			i++;
		if (positions[i].edge_idx>=0)
		{
			rbs.add_real_feature(f_idx,positions[i].node_score);
		}
	}
	f_idx++;

	static vector<int> frag_idx_oris;
	if (frag_idx_oris.size() == 0)
	{
		const vector<FragmentType>& all_frags = config->get_all_fragments();
		frag_idx_oris.resize(all_frags.size(),NEG_INF);
		int i;
		for (i=0; i<all_frags.size(); i++)
			frag_idx_oris[i]=all_frags[i].orientation;
	}

	int num_dual_ori=0;
	int num_with1=0, num_with2=0, num_with_alot=0;
	vector<int> oris;
	for (i=0; i<positions.size(); i++)
	{
		const Breakage *breakage = positions[i].breakage;
		if (! breakage)
			continue;
		const int num_frags = breakage->fragments.size();
		if (num_frags ==1)
		{
			num_with1++;
		}
		else if (num_frags == 2)
		{
			num_with2++;
		}
		else if (num_frags > 5)
		{
			num_with_alot++;
		}

		int num_pre=0,num_suf=0;
		int j;
		for (j=0; j<breakage->fragments.size(); j++)
		{
			const int ori = frag_idx_oris[breakage->fragments[j].frag_type_idx];
			if (ori==PREFIX)
			{
				num_pre++;
			}
			else
				num_suf++;
		}
		if (num_pre == 0 || num_suf - num_pre > 5)
		{
			oris.push_back(SUFFIX);
		}
		else if (num_suf == 0 || num_pre - num_suf>5)
		{
			oris.push_back(PREFIX);
		}
		else
			oris.push_back(99);

		if (num_pre>0 && num_suf>0)
			num_dual_ori++;
	}

	int prev=99;
	int switches=0;
	for (i=0; i<oris.size(); i++)
	{
		if (oris[i] != 99)
		{
			if (prev != 99 && prev != oris[i])
				switches++;
			prev=oris[i];
		}
	}

	if (oris.size()>0)
	{
		const float one_over = 1.0/(float)oris.size();
		rbs.add_real_feature(f_idx++,(float)num_with1*one_over);
		rbs.add_real_feature(f_idx++,(float)num_with2*one_over);
		rbs.add_real_feature(f_idx++,(float)num_with_alot*one_over);
		rbs.add_real_feature(f_idx++,num_dual_ori*one_over);
	}
	else
		f_idx+=4;
	rbs.add_real_feature(f_idx++,(float)switches);

/*	feature_names.push_back("PRM #breakages with 1 frag detected");
	feature_names.push_back("PRM #breakages with 2 frag detected");
	feature_names.push_back("PRM #breakages with > 5 frags detected");

	feature_names.push_back("PRM #breakages with dual orientation frags");
	feature_names.push_back("PRM #orientation switches");*/
}


void DeNovoPartitionModel::fill_pmcsqs_features(
		const PeptideSolution& sol,
		const vector<PmcSqsChargeRes>& res,
		const PMCSQS_Scorer *pmc_model,
		RankBoostSample& rbs) const
{
	int f_idx = pmc_start_idx;

	const mass_t& mz1 = res[sol.charge].mz1;
	const mass_t& mz2 = res[sol.charge].mz2;
	const mass_t pm1_with_19 = (mz1>0 ? mz1 * sol.charge - (sol.charge-1)*MASS_PROTON : NEG_INF);
	const mass_t pm2_with_19 = (mz2>0 ? mz2 * sol.charge - (sol.charge-1)*MASS_PROTON : NEG_INF); 
	const float prob_charge = res[sol.charge].min_comp_prob;
	
	if (mz1>0)
	{
		rbs.add_real_feature(f_idx,res[sol.charge].sqs_prob);
		rbs.add_real_feature(f_idx+1,prob_charge);

		if (prob_charge>0.95)
		{
			rbs.add_real_feature(f_idx+2,pm1_with_19-sol.pep.get_mass_with_19());
		}
		else
			rbs.add_real_feature(f_idx+3,pm1_with_19-sol.pep.get_mass_with_19());

		rbs.add_real_feature(f_idx+4,res[sol.charge].score1);
	}
	
	f_idx+=5;

	if (mz2>0)
	{
		rbs.add_real_feature(f_idx++,res[sol.charge].score2);
		rbs.add_real_feature(f_idx++,pm2_with_19-sol.pep.get_mass_with_19());
	}
	else
		f_idx+=2;

	float max_prob=-1;
	int c;
	for (c=1; c<res.size(); c++)
	{
		if (c==sol.charge)
			continue;
		if (res[c].min_comp_prob>max_prob)
			max_prob=res[c].min_comp_prob;
	}

	if (max_prob>=0)
		rbs.add_real_feature(f_idx,max_prob);
	f_idx++;

	const vector< vector< PMCRankStats > >& curr_stats = pmc_model->get_curr_spec_rank_pmc_tables();
	const vector< PMCRankStats >& charge_stats = curr_stats[sol.charge];

	int i;
	int mz_idx =0;
	while (mz_idx<charge_stats.size() && charge_stats[mz_idx].m_over_z>mz1)
		mz_idx++;

	if (mz_idx>0 && mz_idx<charge_stats.size() && 
		charge_stats[mz_idx].m_over_z-mz1>mz1-charge_stats[mz_idx-1].m_over_z)
		mz_idx--;


	int best_idx=0;
	for (i=0; i<charge_stats.size(); i++)
		if (charge_stats[i].rank_score>charge_stats[best_idx].rank_score)
			best_idx = i;
	
	const float score_diff = charge_stats[best_idx].rank_score - charge_stats[mz_idx].rank_score;
	if (prob_charge>0.95)
	{
		rbs.add_real_feature(f_idx,score_diff);
	}
	else if (prob_charge>0.7)
	{
		rbs.add_real_feature(f_idx+1,score_diff);
	}
	else
	{
		rbs.add_real_feature(f_idx+2,score_diff);
	}
	f_idx+=3;	
}

void DeNovoPartitionModel::fill_composition_features(
							const PeptideSolution& sol,
							const Config *config,
							PeptideCompAssigner *comp_assigner,
							const SeqPath& path,
							RankBoostSample& rbs) const
{
	PeptideCompStats comp_stats;

	comp_assigner->fill_peptide_stats(sol.pep,comp_stats);

	const int num_aas = sol.pep.get_num_aas();
	const vector<int> sol_aas = sol.pep.get_amino_acids();

	int f_idx = comp_start_idx;

	if (sol.reaches_n_terminal)
		rbs.add_real_feature(f_idx,(float)comp_stats.start_comp[3]);
	f_idx++;

	if (sol.reaches_c_terminal)
		rbs.add_real_feature(f_idx,(float)comp_stats.end_comp[3]);
	f_idx++;

	const int *len3_counts = comp_stats.cat_counts[3];

	if (1)
	{
		rbs.add_real_feature(f_idx++,(float)(len3_counts[19]+len3_counts[20]));
		rbs.add_real_feature(f_idx++,(float)(len3_counts[15]+len3_counts[16]+len3_counts[17]+len3_counts[18]));
		rbs.add_real_feature(f_idx++,(float)(len3_counts[7]+len3_counts[8]+len3_counts[9]+len3_counts[10]+
											 len3_counts[11]+len3_counts[12]+len3_counts[13]+len3_counts[14]));
		rbs.add_real_feature(f_idx++,(float)(len3_counts[3]+len3_counts[4]+len3_counts[5]+len3_counts[6]));
		rbs.add_real_feature(f_idx++,(float)(len3_counts[1]+len3_counts[2]));
	}
	else
		f_idx+=5;
	
	int avg = 0;
	int cat;
	for (cat=1; cat<=MAX_COMP_CAT; cat++)
	{
		avg+=cat*len3_counts[cat];
		if (len3_counts[cat]>0)
			break;
	}
	int min_cat = cat++;
	for ( ; cat<=MAX_COMP_CAT; cat++)
		avg+=cat*len3_counts[cat];

	if (1)
		rbs.add_real_feature(f_idx,min_cat);
	f_idx++;

	if (num_aas>0)
		rbs.add_real_feature(f_idx,avg/(float)num_aas);
	f_idx++;
	

	vector<score_pair> pairs;
	if (path.positions.size()>3)
	{
		const int start_idx = 1;
		const int end_idx   = path.positions.size()-2;

		int i;
		for (i=start_idx; i<end_idx; i++)
		{
			const PathPos& pos = path.positions[i];
			if (pos.node_idx>0)
			{
				score_pair p;
				p.idx = i;
				p.score = pos.node_score;
				pairs.push_back(p);
			}
		}
		sort(pairs.begin(),pairs.end());

		for (i=0; i<pairs.size() && i<4; i++)
		{
			const int idx = pairs[i].idx;
			int prev_idx;
			for (prev_idx = idx-1; prev_idx>=0; prev_idx--)
				if (path.positions[prev_idx].edge_idx>=0)
					break;

			int prev_cat = -1;
			if (prev_idx>=0 && idx- prev_idx<=3)
			{
				prev_cat = comp_assigner->get_aa_category(idx-prev_idx,&sol_aas[prev_idx],
					(prev_idx == 0 && sol.reaches_n_terminal), false);
			}

			int next_cat = -1;
			if (path.positions[idx].edge_idx>=0)
			{
				int next_idx;
				for (next_idx = idx+1; next_idx<=num_aas; next_idx++)
					if (path.positions[next_idx].node_idx>0)
						break;
				if (next_idx<=num_aas && next_idx-idx<=3)
				{
					next_cat = comp_assigner->get_aa_category(next_idx-idx,&sol_aas[idx],
					false, (next_idx == num_aas && sol.reaches_c_terminal));
				}
			}

			int span_aas[2]={path.positions[idx-1].aa,path.positions[idx].aa};
			int span_cat = comp_assigner->get_aa_category(2,span_aas, false, false);

			rbs.add_real_feature(f_idx+i*3,(float)prev_cat);
			rbs.add_real_feature(f_idx+i*3+1,(float)next_cat);
			rbs.add_real_feature(f_idx+i*3+2,(float)span_cat);
		}
	}
	f_idx+=12;

	const vector<int>& org_aas = config->get_org_aa();
	vector<int> aa_counts;
	int max_val=2;
	if (sol_aas.size()>14)
		max_val=3;
	if (sol_aas.size()>22)
		max_val=4;
	aa_counts.resize(Val+1,0);
	int i;
	for (i=0; i<sol_aas.size(); i++)
		aa_counts[org_aas[sol_aas[i]]]++;
	aa_counts[Leu]+=aa_counts[Ile];

	int c=0;
	for (i=Ala; i<=Val; i++)
	{
		if (i==Ile)
			continue;
		
		if (aa_counts[i]>0)
		{
			rbs.add_real_feature(f_idx+c,(aa_counts[i]>max_val ? max_val : aa_counts[i]));
		}
		c++;
	}
	f_idx+=c;

	const int max_idx = sol_aas.size()-1;
	int num_W=0;
	int num_Q=0;
	int num_N=0;
	int num_XG=0;
	for (i=0; i<max_idx; i++)
	{
		if (path.positions[i].edge_idx>=0)
			continue;

		const int aa1=sol_aas[i];
		const int aa2=sol_aas[i+1];

		if (aa2 == Gly && (aa1 == Glu || aa1 == Gly || aa1 == Ala))
			num_XG++;

		if (aa1 == Gly && aa2 == Gly)
		{
			num_N++;
			continue;
		}

		if ((aa1 == Ala && aa2 == Gly) ||
			(aa1 == Gly && aa2 == Ala))
		{
			num_Q++;
			continue;
		}

		if ((aa1 == Glu && aa2 == Gly) ||
			(aa1 == Gly && aa2 == Glu) ||
			(aa1 == Ala && aa2 == Asp) ||
			(aa1 == Asp && aa2 == Ala) ||
			(aa1 == Val && aa2 == Ser) ||
			(aa1 == Ser && aa2 == Val))
		{
			num_W++;
			continue;
		}
	}

	const int num_problematic = (num_W+num_Q+num_N);
	if (num_problematic>0)
		rbs.add_real_feature(f_idx,(float)num_problematic);
	f_idx++;

	if (num_W>0)
		rbs.add_real_feature(f_idx,(float)num_W);
	f_idx++;
	if (num_Q>0)
		rbs.add_real_feature(f_idx,(float)num_Q);
	f_idx++;
	if (num_N>0)
		rbs.add_real_feature(f_idx,(float)num_N);
	f_idx++;
	if (num_XG>0)
		rbs.add_real_feature(f_idx,(float)num_XG);
	f_idx++;

/*	feature_names.push_back("PEP COMP #double EG,GE,AD,DA,VS,SV");
	feature_names.push_back("PEP COMP #double ");
	feature_names.push_back("PEP COMP #double GG");
	feature_names.push_back("PEP COMP #double GA");
	feature_names.push_back("PEP COMP #double AG");
	feature_names.push_back("PEP COMP #double SL");*/

}

void DeNovoPartitionModel::fill_peak_offset_features(
								   const Config *config,
								   const PeptideSolution& sol,
								   const vector< vector<mass_t> >& masses,
								   const vector< vector<intensity_t> >& intens,
								   RankBoostSample& rbs) const
{
	int f_idx = peak_offset_start_idx;

	vector<mass_t> break_masses;
	sol.pep.calc_expected_breakage_masses(config,break_masses);
	const int max_break = break_masses.size()-1;
	const mass_t pm_with_19 = sol.pm_with_19;

	int f;
	for (f=0; f<ppp_frag_type_idxs.size() && f<2; f++)
	{
		const int frag_idx = ppp_frag_type_idxs[f];
		const vector<intensity_t>& frag_intens = intens[frag_idx];
		const vector<mass_t>&	   frag_masses = masses[frag_idx];
		const FragmentType& fragment = config->get_fragment(frag_idx);
		vector<mass_t> exp_peak_masses;
		vector<int>	   inten_ranks, rank_positions;

		exp_peak_masses.resize(break_masses.size(),NEG_INF);
		find_inten_ranks(frag_intens, inten_ranks);
		rank_positions.resize(inten_ranks.size(),NEG_INF);
		int c;
		for (c=1; c<inten_ranks.size(); c++)
			if (frag_intens[c]>0)
				rank_positions[inten_ranks[c]]=c;

		// self offset features

		mass_t max_self_off=0;
		mass_t avg_self_off=0;
		int	   num_frags_detected=0;

		int b;
		for (b=1; b<=max_break; b++)
		{
			exp_peak_masses[b]=fragment.calc_expected_mass(break_masses[b],pm_with_19);
			if (frag_intens[b]<=0)
				continue;

			const float& frag_mass = frag_masses[b];
			const float& peak_mass = exp_peak_masses[b];
			mass_t offset=fabs(frag_masses[b]-exp_peak_masses[b]);

			if (offset>3.0)
			{
				cout << "Error: bad peak offset calculations!" << endl;

				exp_peak_masses[b]=fragment.calc_expected_mass(break_masses[b],pm_with_19);
				
				const float& frag_mass = frag_masses[b];
				const float& peak_mass = exp_peak_masses[b];
				mass_t offset=fabs(frag_masses[b]-exp_peak_masses[b]);

				exit(1);
			}
			avg_self_off+=offset;
			if (offset>max_self_off)
				max_self_off=offset;
			num_frags_detected++;
		}

		rbs.add_real_feature(f_idx++,num_frags_detected);

		if (num_frags_detected==0)
		{
			f_idx+=7;
			continue;
		}
		
		rbs.add_real_feature(f_idx++,max_self_off);
		if (num_frags_detected>0)
			rbs.add_real_feature(f_idx,(avg_self_off/num_frags_detected));
		f_idx++;
		
		if (num_frags_detected<5)
		{
			f_idx+=4;
			continue;
		}

		// between peak offsets
		int num_gaps=0;
		mass_t max_consec_offset=0;
		mass_t avg_consec_offset=0;

		for (b=1; b<exp_peak_masses.size(); b++)
			if (frag_intens[b]>0)
				break;
		int prev=b++;

		bool in_gap=false;
		int num_gap_offsets=0;
		for ( ; b<exp_peak_masses.size(); b++)
			if (frag_intens[b]<=0)
			{
				if (! in_gap)
				{
					num_gaps++;
					in_gap=true;
				}
			}
			else
			{
				in_gap=false;
				const mass_t offset=fabs(exp_peak_masses[b] -
										 exp_peak_masses[prev] -
										 frag_masses[b] + 
										 frag_masses[prev]);
				if (fabs(offset)>3.0)
				{
					
					cout << "Exp " << b << " " << exp_peak_masses[b] << endl;
					cout << "Exp " << prev << " " << exp_peak_masses[prev] << endl;
					cout << "Frag " << b << " " << frag_masses[b] << endl;
					cout << "Frag " << prev << " " << frag_masses[prev] << endl;
					exit(0);
				}



				if (offset>max_consec_offset)
					max_consec_offset=offset;
				avg_consec_offset+=offset;
				num_gap_offsets++;
			}

		rbs.add_real_feature(f_idx++,max_consec_offset);
		if (num_gap_offsets>0)
			rbs.add_real_feature(f_idx,avg_consec_offset/num_gap_offsets);
		f_idx++;
		
		// Peak grab features
		vector<mass_t> offset_diffs;
		offset_diffs.clear();
		for (b=1; b<max_break-2; b++)
		{
			if (frag_intens[b]<=0)
				continue;

			if (frag_intens[b+1]>0 && frag_intens[b+2]>0)
			{
				const mass_t offset1 = fabs(exp_peak_masses[b+2]-exp_peak_masses[b]-frag_masses[b+2]+frag_masses[b]);
				const mass_t offset2 = fabs(exp_peak_masses[b+1]-exp_peak_masses[b]-frag_masses[b+1]+frag_masses[b]);
				const mass_t offset3 = fabs(exp_peak_masses[b+2]-exp_peak_masses[b+1]-frag_masses[b+2]+frag_masses[b+1]);
				offset_diffs.push_back(fabs(offset1-offset2-offset3));
			}

			if (frag_intens[b+1]>0 && frag_intens[b+2]<=0 && frag_intens[b+3]>0)
			{
				const mass_t offset1 = fabs(exp_peak_masses[b+3]-exp_peak_masses[b]-frag_masses[b+3]+frag_masses[b]);
				const mass_t offset2 = fabs(exp_peak_masses[b+1]-exp_peak_masses[b]-frag_masses[b+1]+frag_masses[b]);
				const mass_t offset3 = fabs(exp_peak_masses[b+3]-exp_peak_masses[b+1]-frag_masses[b+3]+frag_masses[b+1]);
				offset_diffs.push_back(fabs(offset1-offset2-offset3));
			}

			if (frag_intens[b+1]<=0 && frag_intens[b+2]>0 && frag_intens[b+3]>0)
			{
				const mass_t offset1 = fabs(exp_peak_masses[b+3]-exp_peak_masses[b]-frag_masses[b+3]+frag_masses[b]);
				const mass_t offset2 = fabs(exp_peak_masses[b+2]-exp_peak_masses[b]-frag_masses[b+2]+frag_masses[b]);
				const mass_t offset3 = fabs(exp_peak_masses[b+3]-exp_peak_masses[b+2]-frag_masses[b+3]+frag_masses[b+2]);
				offset_diffs.push_back(fabs(offset1-offset2-offset3));
			}
		}
		sort(offset_diffs.begin(),offset_diffs.end());

		int i,counter=0;
		for (i=offset_diffs.size()-1; i>=0; i--)
		{
			if (offset_diffs[i]==0)
				break;

			rbs.add_real_feature(f_idx+counter,offset_diffs[i]);

			if (++counter==3)
				break;
		}
		f_idx+=3; 
	}
}


void DeNovoPartitionModel::fill_ann_peak_features(const PeptideSolution& sol,
								const vector< vector<mass_t> >& masses,
								const vector< vector<intensity_t> >& intens,
								const AnnotatedSpectrum& as,
								RankBoostSample& rbs) const
{
	int f_idx = ann_peak_start_idx;

	const intensity_t total_inten = as.getTotalPeakIntensity();
	intensity_t ann_inten=0;
	int num_ann_inten=0;
	vector<int> frag_counts;
	frag_counts.resize(intens.size(),0);
	int f;
	for (f=0; f<intens.size(); f++)
	{
		const int num_intens = intens[f].size();
		int c;
		for (c=1; c<num_intens; c++)
		{	
			const intensity_t& frag_inten = intens[f][c];
			if (frag_inten>0)
			{
				ann_inten += frag_inten;
				frag_counts[f]++;
			}
		}
		num_ann_inten+=frag_counts[f];
	}

	int length = sol.pep.get_num_aas();
	if (length<7)
		length=7;

	rbs.add_real_feature(f_idx++,sol.pm_with_19-as.get_org_pm_with_19());
	rbs.add_real_feature(f_idx++,length);
	if (total_inten>0)
		rbs.add_real_feature(f_idx,(float)ann_inten/total_inten);
	f_idx++;
	if (as.getNumPeaks()>0)
		rbs.add_real_feature(f_idx,(float)num_ann_inten/as.getNumPeaks());
	f_idx++;

	const Peak* const peaks = as.getPeaks();
	const vector< vector<PeakAnnotation> >& peak_anns = as.get_peak_annotations();

	int count25=0;
	int count_half=0;
	int count_top_third=0;
	int count_mid_third=0;
	int count_last_third=0;
	const int half_num_peaks = (as.getNumPeaks()/2);
	const int third_num_peaks = (as.getNumPeaks()/3);
	const int two_third_num_peaks = third_num_peaks*2;

	int i;
	for (i=0; i<as.getNumPeaks(); i++)
	{
		if (peak_anns[i].size()==0)
			continue;

		const int rank = as.get_peak_rank(i);
		if (rank>=two_third_num_peaks)
		{
			count_last_third++;
			continue;
		}

		if (rank>=third_num_peaks)
		{
			count_mid_third++;
		}
		else
			count_top_third++;


		if (rank<half_num_peaks)
			count_half++;
		if (rank<25)
			count25++;
	}
	rbs.add_real_feature(f_idx++,(float)count25);
	rbs.add_real_feature(f_idx++,(float)count_half);
	rbs.add_real_feature(f_idx++,(float)(count_top_third-count_mid_third));
	rbs.add_real_feature(f_idx++,(float)(count_top_third-count_last_third));
	rbs.add_real_feature(f_idx++,(float)(count_mid_third-count_last_third));

	const int num_frags = as.getConfig()->get_all_fragments().size();

	int max_f = 7;
	if (num_frags<max_f)
		max_f = num_frags;
	if (intens.size()<max_f)
		max_f=intens.size();
	for (f=0; f<max_f; f++)
	{
		rbs.add_real_feature(f_idx+f,(float)frag_counts[f]);
	}
	f_idx+=7;
}


void DeNovoPartitionModel::fill_inten_balance_features(const Config *config,
													   const PeptideSolution& sol, 
													   const SeqPath& path,
													   RankBoostSample& rbs) const
{
	int f_idx = inten_balance_start_idx;

	const vector<int>& amino_acids = sol.pep.get_amino_acids();
	const vector<PathPos>& positions = path.positions;
	int n_idx = 0;

	while (n_idx<positions.size() && 
		( ! positions[n_idx].breakage || positions[n_idx].breakage->fragments.size() == 0))
		n_idx++;

	int c_idx = positions.size()-1;
	while (c_idx>0 && 
		( ! positions[c_idx].breakage || positions[c_idx].breakage->fragments.size() == 0))
		c_idx--;

	const int nc_diff = c_idx-n_idx;

	rbs.add_real_feature(f_idx++,(float)nc_diff);
	if (nc_diff<4)
		return;

	vector<int> pos_idxs;
	const int mid_idx = (n_idx + c_idx)/2;
	int i;
	for (i=mid_idx-2; i<=mid_idx+2; i++)
		if (positions[i].node_idx>0)
			pos_idxs.push_back(i);


	intensity_t pre_inten=0;
	intensity_t suf_inten=0;

	for (i=0; i<pos_idxs.size(); i++)
	{
		Breakage *breakage = positions[pos_idxs[i]].breakage;

		int j;
		for (j=0; j<breakage->fragments.size();  j++)
		{
			const BreakageFragment& bf = breakage->fragments[j];
			const FragmentType& frag = config->get_fragment(bf.frag_type_idx);
			if (frag.orientation == PREFIX)
			{
				pre_inten += bf.intensity;
			}
			else
				suf_inten += bf.intensity;
		}
	}
	const intensity_t sum_inten = pre_inten + suf_inten;
	if (sum_inten<=0)
		return;

	const float pre_ratio = pre_inten/sum_inten;

	intensity_t all_pre_inten=0;
	intensity_t all_suf_inten=0;

	for (i=0; i<positions.size(); i++)
	{
		Breakage *breakage = positions[i].breakage;
		if (! breakage)
			continue;

		int j;
		for (j=0; j<breakage->fragments.size();  j++)
		{
			const BreakageFragment& bf = breakage->fragments[j];
			const FragmentType& frag = config->get_fragment(bf.frag_type_idx);
			if (frag.orientation == PREFIX)
			{
				all_pre_inten += bf.intensity;
			}
			else
				all_suf_inten += bf.intensity;
		}
	}
	const intensity_t all_sum_inten = all_pre_inten + all_suf_inten;
	if (all_sum_inten<=0)
		return;

	const float all_pre_ratio = all_pre_inten/all_sum_inten;


	
	// special N C side aa indicators
	int num_nH=0, num_cH=0;
	int num_nK=0, num_cK=0;
	int num_nR=0, num_cR=0;
	
	for (i=0; i<=mid_idx; i++)
	{
		if (amino_acids[i] == His)
			num_nH++;

		if (amino_acids[i] == Lys)
			num_nK++;

		if (amino_acids[i] == Arg)
			num_nR++;
	}

	// uses regular amino acid codes
	if (sol.most_basic_aa_removed_from_n>0)
	{
		if (sol.most_basic_aa_removed_from_n == His)
		{
			num_nH++;
		}
		else if (sol.most_basic_aa_removed_from_n == Lys) 
		{
			num_nK++;
		}
		else if (sol.most_basic_aa_removed_from_n == Arg) 
			num_nR++;
	}

	for (i=mid_idx; i<positions.size()-1; i++)
	{
		if (amino_acids[i] == His)
			num_cH++;

		if (amino_acids[i] == Lys)
			num_cK++;

		if (amino_acids[i] == Arg)
			num_cR++;
	}

	// uses regular amino acid codes
	if (sol.most_basic_aa_removed_from_c>0)
	{
		if (sol.most_basic_aa_removed_from_c == His)
		{
			num_cH++;
		} else if (sol.most_basic_aa_removed_from_c == Lys) 
		{
			num_cK++;
		}
		else if (sol.most_basic_aa_removed_from_c == Arg) 
			num_cR++;
	}

	const int RKH_n_combo_idx = calc_RKH_combo_idx(num_nR,num_nK,num_nH);
	const int RKH_c_combo_idx = calc_RKH_combo_idx(num_cR,num_cK,num_cH);
	const float RKH_liniar_pair_idx = RKH_pair_matrix[RKH_n_combo_idx][RKH_c_combo_idx];

	rbs.add_real_feature(f_idx++,(float)RKH_n_combo_idx);
	rbs.add_real_feature(f_idx++,(float)RKH_c_combo_idx);
	rbs.add_real_feature(f_idx++,RKH_liniar_pair_idx);

	if (RKH_liniar_pair_idx<=-4)
	{
		rbs.add_real_feature(f_idx,pre_ratio);
		rbs.add_real_feature(f_idx+5,all_pre_ratio);
	} 
	else if (RKH_liniar_pair_idx<=-2)
	{
		rbs.add_real_feature(f_idx+1,pre_ratio);
		rbs.add_real_feature(f_idx+6,all_pre_ratio);
	}
	else if (RKH_liniar_pair_idx<=1)
	{
		rbs.add_real_feature(f_idx+2,pre_ratio);
		rbs.add_real_feature(f_idx+7,all_pre_ratio);
	}
	else if (RKH_liniar_pair_idx<=3)
	{
		rbs.add_real_feature(f_idx+3,pre_ratio);
		rbs.add_real_feature(f_idx+8,all_pre_ratio);
	}
	else
	{
		rbs.add_real_feature(f_idx+4,pre_ratio);
		rbs.add_real_feature(f_idx+9,all_pre_ratio);
	}

	f_idx+=10;

/*	feature_names.push_back("INTEN BAL c_idx - n_idx");
	feature_names.push_back("INTEN BAL RHK N");
	feature_names.push_back("INTEN BAL RHK C");
	feature_names.push_back("INTEN BAL RHK pair");
	feature_names.push_back("INTEN BAL prefix prop, pair -4,-5");
	feature_names.push_back("INTEN BAL prefix prop, pair -2,-3");
	feature_names.push_back("INTEN BAL prefix prop, pair -1,0,+1");
	feature_names.push_back("INTEN BAL prefix prop, pair +2,+3");
	feature_names.push_back("INTEN BAL prefix prop, pair +4,+5");*/
}

void DeNovoPartitionModel::fill_tryp_terminal_features(const PeptideSolution& sol, 
													   const SeqPath& path,
													   RankBoostSample& rbs) const
{
	int f_idx = tryp_terminal_start_idx;

	const Peptide& pep = sol.pep;
	const int num_aas = pep.get_num_aas();
	const vector<int>& amino_acids = pep.get_amino_acids();

	int first_aa = amino_acids[0];
	int last_aa  = amino_acids[num_aas-1];
	int aa_before = pep.get_aa_before();
	int aa_after  = pep.get_aa_after();

	int num_bad_tryp=0;
	int num_tryp=0;
	if (aa_before<Ala || aa_before == Lys || aa_before == Arg || ! sol.reaches_n_terminal)
		num_tryp++;

	if (sol.type>=0)
	{
		if (last_aa == Lys || last_aa == Arg || aa_after <Ala)
			num_tryp++;
	}
	else
	{
		if (last_aa == Lys || last_aa == Arg || ! sol.reaches_c_terminal)
			num_tryp++;
	}

	// count missed tryptic (not including Pro)
	int num_missed_tryp=0;
	int i;
	for (i=0; i<amino_acids.size()-1; i++)
		if ((amino_acids[i] == Arg || amino_acids[i] == Lys) && amino_acids[i+1] != Pro)
			num_missed_tryp++;

	if (! sol.reaches_c_terminal && (last_aa == Lys || last_aa == Arg))
		num_missed_tryp++;

	rbs.add_real_feature(f_idx++,num_tryp);
	rbs.add_real_feature(f_idx++,num_missed_tryp);
	if (sol.reaches_c_terminal)
		rbs.add_real_feature(f_idx,last_aa);
	f_idx++;

	if (path.positions.size()>2)
	{
		PathPos digest_pos = path.positions[path.positions.size()-2];
		if (digest_pos.node_idx>0)
		{
			const int num_frags = digest_pos.breakage->fragments.size();
			if (last_aa == Arg)
			{
				rbs.add_real_feature(f_idx,num_frags);
			}
			else if (last_aa == Lys)
			{
				rbs.add_real_feature(f_idx+1,num_frags);
			}
			else
				rbs.add_real_feature(f_idx+2,num_frags);
		}
	}
	f_idx += 3;

	if (sol.reaches_n_terminal)
	{
		if (last_aa == Arg)
		{
			rbs.add_real_feature(f_idx,first_aa);
		}
		else if (last_aa == Lys)
		{
			rbs.add_real_feature(f_idx+1,first_aa);
		}
		else
			rbs.add_real_feature(f_idx+2,first_aa);
	}
	f_idx+=3;
}



void DeNovoPartitionModel::fill_PTM_peak_features(
								const Config *config,
								const PeptideSolution& sol,
								const vector< vector<mass_t> >& masses,
								const vector< vector<intensity_t> >& intens,
								const AnnotatedSpectrum& as,
								RankBoostSample& rbs) const
{
	static const int m16_aa_idx = config->get_aa_from_label("M+16");
	static const int b_frag_idx = config->get_frag_idx_from_label("b");
	static const int y_frag_idx = config->get_frag_idx_from_label("y");
	static const mass_t tolerance = (config->getTolerance()>0.2 ? 
									 config->getTolerance() * 0.5 : config->getTolerance());
	int f_idx = PTM_peak_start_idx;

	// find the most prominant occerences of M+16
	if (m16_aa_idx>0)
	{	
		const vector<int>& amino_acids = sol.pep.get_amino_acids();
		vector<score_pair> pairs;
		pairs.clear();
		int i;
		for (i=0; i<amino_acids.size(); i++)
			if (amino_acids[i] == m16_aa_idx)
			{
				float total_inten=0;
				int f;
				for (f=0; f<ppp_fragments.size(); f++)
				{
					if (ppp_fragments[f].orientation == PREFIX)
					{
						const float& inten = intens[ppp_frag_type_idxs[f]][i+1];
						if (inten>0)
							total_inten+= inten;
					}
					else
					{
						const float& inten = intens[ppp_frag_type_idxs[f]][i];
						if (inten>0)
							total_inten+= inten;
					}
				}
				pairs.push_back(score_pair(i,total_inten));
			}
		sort(pairs.begin(),pairs.end());
		vector<mass_t> break_masses;
		sol.pep.calc_expected_breakage_masses(config,break_masses);
		if (pairs.size()>0)
		{
			const int num_aa = amino_acids.size();
			for (i=0; i<pairs.size() && i<MAX_NUM_M16s; i++)
			{
				const int aa_idx=pairs[i].idx;

				if (intens[b_frag_idx].size()>aa_idx+1)
				{
					const float b_inten = intens[b_frag_idx][aa_idx+1];
					const int pre_minus_63_idx = as.findPeakWithMaxIntensity(break_masses[aa_idx+1]-63.0,tolerance);
					
					float ratio_pre_minus_63=-10.0;
					if (b_inten == 0)
					{
						if (pre_minus_63_idx>=0)
							ratio_pre_minus_63 = 11.0;
					}
					else
					{
						ratio_pre_minus_63 = ( pre_minus_63_idx<0 ? 0 : as.getPeakIntensity(pre_minus_63_idx)/b_inten);
						if (ratio_pre_minus_63>10.0)
							ratio_pre_minus_63=10.0;
					}
					rbs.add_real_feature(f_idx+i*2,ratio_pre_minus_63);
				}

				if (intens[y_frag_idx].size()>aa_idx)
				{
					const mass_t pep_mass = sol.pep.get_mass();
					const float y_inten = intens[y_frag_idx][aa_idx];
					const int suf_minus_45_idx = as.findPeakWithMaxIntensity((pep_mass - break_masses[aa_idx])-45.0,tolerance);
				
					float ratio_suf_minus_45=-10;
				
					if (y_inten == 0)
					{
						if (suf_minus_45_idx>=0)
							ratio_suf_minus_45 = 11.0;
					}
					else
					{
						ratio_suf_minus_45 = ( suf_minus_45_idx<0 ? 0 : as.getPeakIntensity(suf_minus_45_idx)/y_inten);
						if (ratio_suf_minus_45>10.0)
							ratio_suf_minus_45=10.0;
	
					}
					rbs.add_real_feature(f_idx+i*2+1,ratio_suf_minus_45);
				}
			}
		}
	}

	f_idx += 4;
}



void PeptideRankScorer::init_tables(bool silent_ind)
{
	AllScoreModels* zsm = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeakRankModel *& peak_model = zsm->get_peak_prediction_model_ptr(model_type);

	if (! zsm || ! zsm->get_ind_pmcsqs_was_intialized())
	{
		cout << "Error: must first initialize the fragment model!" << endl;
		exit(1);
	}

	if (! peak_model || peak_model->get_size_thresholds().size()<=0)
	{
		cout << "Error: must first initialize the peak model!" << endl;
		exit(1);
	}

	const vector< vector< mass_t> >& size_thresholds = peak_model->get_size_thresholds();

	this->dnv_part_models.resize(size_thresholds.size());
	int i;

	if (! silent_ind)
		cout << "Init tables:" << endl;

	for (i=1; i<dnv_part_models.size(); i++)
	{
		dnv_part_models[i].resize(size_thresholds[i].size()+1,NULL);

		if (! silent_ind)
			cout << i << "\t" << size_thresholds[i].size()+1 << endl;
	}
}


void DeNovoPartitionModel::init_features(int model_type, int _charge, int _size_idx,
							const vector<int>& ppp_frags, Config *config)
{
	charge = _charge;
	size_idx = _size_idx;

	ScalingFactor def_scale;
	def_scale.max_pm_with_19 = POS_INF;
	def_scale.score_shift = 0;
	def_scale.score_scale = 1.0;
	scaling_factors.clear();
	scaling_factors.push_back(def_scale);

	feature_names.clear();

	ind_was_initialized = true;

	use_PTM_peak_features	   = true;
	use_tryp_terminal_features = true;
	use_ann_peak_features = true;
	use_inten_balance_features = (model_type != 3);
	use_peak_offset_features = true;
	use_comp_features = true;
	use_pmc_features = (model_type != 3);
	use_prm_features = true;
	use_ppp_features = false;
	use_combined_ppp_features = true;

	if (use_combined_ppp_features && use_ppp_features)
	{
		cout << "Error: must choose combine or regular ppp features, not both!" << endl;
		exit(1);
	}

	num_ppp_frags=ppp_frags.size();
	ppp_frag_type_idxs=ppp_frags;
	ppp_fragments.resize(ppp_frag_type_idxs.size());
	int f;
	for (f=0; f<ppp_frag_type_idxs.size(); f++)
		ppp_fragments[f] = config->get_fragment(ppp_frag_type_idxs[f]);


	// features for special PTM peaks (like the ones for M+16)
	PTM_peak_start_idx = feature_names.size();
	if (use_PTM_peak_features)
	{
		feature_names.push_back("M+16 first p-63 ratio");
		feature_names.push_back("M+16 first s-45 ratio");
		feature_names.push_back("M+16 second p-63 ratio");
		feature_names.push_back("M+16 second s-45 ratio");
	}

	// tryptic terminal features
	tryp_terminal_start_idx = feature_names.size();
	if (use_tryp_terminal_features)
	{
		feature_names.push_back("TRYP #num good tryp terminals");
		feature_names.push_back("TRYP #num missed tryp terminals");
		feature_names.push_back("TRYP C-term AA");
		feature_names.push_back("TRYP #frags at digest when C-term is R");
		feature_names.push_back("TRYP #frags at digest when C-term is K");
		feature_names.push_back("TRYP #frags at digest when C-term is other");
		feature_names.push_back("TRYP AA at N-terminal When C-term is R");
		feature_names.push_back("TRYP AA at N-terminal When C-term is K");
		feature_names.push_back("TRYP AA at N-terminal When C-term is other");
	}

	// Ann peak features
	ann_peak_start_idx = feature_names.size();
	if (use_ann_peak_features)
	{
		feature_names.push_back("ANN PEAK diff from org pm_with_19");
		feature_names.push_back("ANN PEAK # aas in peptide");
		feature_names.push_back("ANN PEAK %ann intensity");
		feature_names.push_back("ANN PEAK %ann peaks");
		feature_names.push_back("ANN PEAK #ann in top 25");
		feature_names.push_back("ANN PEAK #ann in top half (up to 50)");
		feature_names.push_back("ANN PEAK #ann in top third - #ann in mid third");
		feature_names.push_back("ANN PEAK #ann in top third - #ann in last third");
		feature_names.push_back("ANN PEAK #ann in mid third - #ann in last third");

		const vector<FragmentType>& all_fragments = config->get_all_fragments();
		int f;
		for (f=0; f<all_fragments.size() && f<7; f++)
		{
			const string frag_label = all_fragments[f].label;
			feature_names.push_back( "ANN PEAK #" + frag_label + " annotated");
		}
	}

	// inten balance features
	inten_balance_start_idx = feature_names.size();
	if (use_inten_balance_features)
	{
		feature_names.push_back("INTEN BAL c_idx - n_idx");
		feature_names.push_back("INTEN BAL RHK N");
		feature_names.push_back("INTEN BAL RHK C");
		feature_names.push_back("INTEN BAL RHK pair");
		feature_names.push_back("INTEN BAL prefix prop, pair -4,-5");
		feature_names.push_back("INTEN BAL prefix prop, pair -2,-3");
		feature_names.push_back("INTEN BAL prefix prop, pair -1,0,+1");
		feature_names.push_back("INTEN BAL prefix prop, pair +2,+3");
		feature_names.push_back("INTEN BAL prefix prop, pair +4,+5");

		feature_names.push_back("INTEN BAL all prefix prop, pair -4,-5");
		feature_names.push_back("INTEN BAL all prefix prop, pair -2,-3");
		feature_names.push_back("INTEN BAL all prefix prop, pair -1,0,+1");
		feature_names.push_back("INTEN BAL all prefix prop, pair +2,+3");
		feature_names.push_back("INTEN BAL all prefix prop, pair +4,+5");
	}

	// Peak offset features
	peak_offset_start_idx=feature_names.size();
	if (use_peak_offset_features)
	{
		int f;
		for (f=0; f<ppp_frag_type_idxs.size() && f<2; f++)
		{
			const int frag_idx = ppp_frag_type_idxs[f];
			const string frag_label = config->get_fragment(frag_idx).label;
			const string prefix = "PEAK OFF " + frag_label + " ";

			feature_names.push_back(prefix+"num frags detected");
			feature_names.push_back(prefix+"max self offset");
			feature_names.push_back(prefix+"avg self offset");

			feature_names.push_back(prefix+"max consecutive offset");
			feature_names.push_back(prefix+"avg consecutive offset");

			// peak grab feature
			feature_names.push_back(prefix+"grab offset #1");
			feature_names.push_back(prefix+"grab offset #2");
			feature_names.push_back(prefix+"grab offset #3"); 
		}

		if (f<2)
		{
			int i;
			for (i=0; i<8; i++)
				feature_names.push_back("PEAK OFF dummy");
		}
	}

	// Peptide composition features
	comp_start_idx=feature_names.size();
	if (use_comp_features)
	{
		feature_names.push_back("PEP COMP start cat N (len 3)");
		feature_names.push_back("PEP COMP end cat C (len 3)");
		feature_names.push_back("PEP COMP len 3 # cat 19-20");
		feature_names.push_back("PEP COMP len 3 # cat 15-18");
		feature_names.push_back("PEP COMP len 3 # cat 7-14");
		feature_names.push_back("PEP COMP len 3 # cat 3-6");
		feature_names.push_back("PEP COMP len 3 # cat 1-2");
		feature_names.push_back("PEP COMP min cat, len 3");
		feature_names.push_back("PEP COMP avg cat, len 3");	
		feature_names.push_back("PEP COMP before cat score 1");
		feature_names.push_back("PEP COMP after cat score 1");
		feature_names.push_back("PEP COMP span cat score 1");
		feature_names.push_back("PEP COMP before cat score 2");
		feature_names.push_back("PEP COMP after cat score 2");
		feature_names.push_back("PEP COMP span cat score 2");
		feature_names.push_back("PEP COMP before cat score 3");
		feature_names.push_back("PEP COMP after cat score 3");
		feature_names.push_back("PEP COMP span cat score 3");
		feature_names.push_back("PEP COMP before cat score 4");
		feature_names.push_back("PEP COMP after cat score 4");
		feature_names.push_back("PEP COMP span cat score 4");

		const vector<string>& aa2label = config->get_aa2label();
		int a;
		for (a=Ala; a<=Val; a++)
		{
			if (a==Ile)
				continue;
			feature_names.push_back("PEP COMP #aa " + aa2label[a]);	
		}

		feature_names.push_back("PEP COMP #problematic double combos");
		feature_names.push_back("PEP COMP #double combo=W");
		feature_names.push_back("PEP COMP #double combo=Q");
		feature_names.push_back("PEP COMP #double combo=N");
		feature_names.push_back("PEP COMP #double problematic combos with XG");
	}

	// PMCSQS features
	pmc_start_idx=feature_names.size();
	if (use_pmc_features)
	{
		feature_names.push_back("PMCSQS sqs prob for peptide charge");
		feature_names.push_back("PMCSQS prob for peptide charge");
		feature_names.push_back("PMCSQS mass diff from pm1, prob>0.95");
		feature_names.push_back("PMCSQS mass diff from pm1, prob<=0.95");
		feature_names.push_back("PMCSQS score1 for peptide charge");
		feature_names.push_back("PMCSQS score2 for peptide charge");
		feature_names.push_back("PMCSQS mass diff from pm2");
		feature_names.push_back("PMCSQS max  prob for other charges");
		feature_names.push_back("PMCSQS score diff from max score with this charge, prob>=0.95");
		feature_names.push_back("PMCSQS score diff from max score with this charge, 0.95>prob>=0.7");
		feature_names.push_back("PMCSQS score diff from max score with this charge, prob<0.7");
	}

	prm_start_idx = feature_names.size();
	if (use_prm_features) // use these feature only with de novo
	{
		const string term_combo[4]={"N/C","N/-C","-N/C","-N/-C"};
		int i;
		for (i=0; i<4; i++)
		{
			if (i>0 && (model_type != 1 && model_type != 3))
				break;

			feature_names.push_back("PRM " + term_combo[i] + " delta mass");
			feature_names.push_back("PRM " + term_combo[i] + " total breakage score");
			feature_names.push_back("PRM " + term_combo[i] + " average breakage score");
			feature_names.push_back("PRM " + term_combo[i] + " normalized average breakage score");
			feature_names.push_back("PRM " + term_combo[i] + " path score");
			feature_names.push_back("PRM " + term_combo[i] + " average path score");
		}

		// adjust to random prob in spectrum (so model can be used with FT)
		if (model_type == 1 || model_type == 3)
		{
			feature_names.push_back("PRM path score");
			feature_names.push_back("PRM total breakage score");
			feature_names.push_back("PRM SeqPath rank");
			feature_names.push_back("PRM multipath score");

			feature_names.push_back("PRM delta score");
			feature_names.push_back("PRM rank, delta score<=1.5");
			feature_names.push_back("PRM rank, 1.5<delta score<=7.5");
			feature_names.push_back("PRM rank, 7.5<delta score<=15");
			feature_names.push_back("PRM rank, delta score>15");

			if (model_type == 3)
			{
				feature_names.push_back("PRM tag, percent in top 5 denovo");
				feature_names.push_back("PRM tag, percent in top 20 denovo");
				feature_names.push_back("PRM tag, percent in all denovo");
				feature_names.push_back("PRM tag, rank if in top 5");
				feature_names.push_back("PRM tag, rank if in top 5-20");
				feature_names.push_back("PRM tag, rank if in top 20-all");
				feature_names.push_back("PRM tag, highest full denovo rank");
			}
		}

		feature_names.push_back("PRM delta num breakage scores (missing)");
	//	feature_names.push_back("PRM num missing edges");
		feature_names.push_back("PRM num forbidden node pairs");
		feature_names.push_back("PRM num breakage scores");
		feature_names.push_back("PRM breakage score min 1");
		feature_names.push_back("PRM breakage score min 2");
		feature_names.push_back("PRM breakage score min 3");
		feature_names.push_back("PRM breakage score min consecutive 3");
		feature_names.push_back("PRM breakage score max consecutive 3");
		feature_names.push_back("PRM breakage score min consecutive 2");
		feature_names.push_back("PRM breakage score max consecutive 2");
		feature_names.push_back("PRM #breakage scores below -10"); 
		feature_names.push_back("PRM #breakage scores 0 - -10");
		feature_names.push_back("PRM #breakage scores 0 - 8");
		feature_names.push_back("PRM #breakage scores 8 - 15");
		feature_names.push_back("PRM #breakage scores above 15");
		feature_names.push_back("PRM %breakage scores below -10"); 
		feature_names.push_back("PRM %breakage scores below 0");
		feature_names.push_back("PRM %breakage scores above 0");
		feature_names.push_back("PRM %breakage scores above 8");
		feature_names.push_back("PRM Score connected to N-terminal");
		feature_names.push_back("PRM Score connected to C-terminal");

		feature_names.push_back("PRM %breakages with 1 frag detected");
		feature_names.push_back("PRM %breakages with 2 frag detected");
		feature_names.push_back("PRM %breakages with > 5 frags detected");
		feature_names.push_back("PRM %breakages with dual orientation frags");
		feature_names.push_back("PRM #orientation switches");
	}
	
	// Peak prediction features
	ppp_start_idx=feature_names.size();
	if (use_ppp_features)
	{
		int f;
		for (f=0; f<ppp_frag_type_idxs.size(); f++)
		{
			const int frag_idx = ppp_frag_type_idxs[f];
			const string frag_label = config->get_fragment(frag_idx).label;
			const string prefix = "PPP " + frag_label + " ";
			feature_names.push_back(prefix+"# observed frags");
			feature_names.push_back(prefix+"# predicted frags");
			feature_names.push_back(prefix+"observation ratio");
			feature_names.push_back(prefix+"# observed frags in top 1 predicted");
			feature_names.push_back(prefix+"# observed frags in top 3 predicted");
			feature_names.push_back(prefix+"# observed frags in top 5 predicted");
			feature_names.push_back(prefix+"# observed frags in top 7 predicted");
			feature_names.push_back(prefix+"% observed frags in top 1/6 predicted");
			feature_names.push_back(prefix+"% observed frags in top 1/3 predicted");
			feature_names.push_back(prefix+"% observed frags in top 1/2 predicted");
			feature_names.push_back(prefix+"% observed frags in top 2/3 predicted");
			feature_names.push_back(prefix+" predicted rank of first missing peak");
			feature_names.push_back(prefix+" predicted rank of second missing peak");
			feature_names.push_back(prefix+" predicted rank of third missing peak");
			feature_names.push_back(prefix+" predicted rank of first+second missing peak");
			feature_names.push_back(prefix+" predicted rank of first+second+third missing peak");
			feature_names.push_back(prefix+"score offset of rank 1");
			feature_names.push_back(prefix+"score offset of rank 2");
			feature_names.push_back(prefix+"score offset of rank 3");
			feature_names.push_back(prefix+"score offset of rank 4");
			feature_names.push_back(prefix+"score offset of rank 5");
			feature_names.push_back(prefix+"score offset of rank 6");
			feature_names.push_back(prefix+"score offset of rank 7");
			feature_names.push_back(prefix+"score offset of rank 8");
			feature_names.push_back(prefix+"score offset of rank 9");
			feature_names.push_back(prefix+"score offset of rank 10");
		}
	
		feature_names.push_back("PPP #comp pairs");
		feature_names.push_back("PPP stat of predicted pair #1");
		feature_names.push_back("PPP stat of predicted pair #2");
		feature_names.push_back("PPP stat of predicted pair #3");
		feature_names.push_back("PPP stat of predicted pair #4");
		feature_names.push_back("PPP stat of predicted pair #5");
		feature_names.push_back("PPP stat of predicted pair #6");
		feature_names.push_back("PPP stat of predicted pair #7");
	}

	combined_ppp_start_idx = feature_names.size();
	if (use_combined_ppp_features)
	{
		feature_names.push_back("COMB PPP mobility");
		feature_names.push_back("COMP PPP frag 1 obs_ratio");
		feature_names.push_back("COMP PPP frag 2 obs_ratio");
		feature_names.push_back("COMP PPP frag 3 obs_ratio");
		feature_names.push_back("COMP PPP num missed peaks");
		
		feature_names.push_back("COMP PPP MOBILE sum ranks of missed 1-5");
		feature_names.push_back("COMP PPP PARTMOBILE sum ranks of missed 1-5");
		feature_names.push_back("COMP PPP NONMOBILE sum ranks of missed 1-5");
		feature_names.push_back("COMP PPP sum ranks of missed 1-5");

		feature_names.push_back("COMP PPP MOBILE sum ranks of missed 6-10");
		feature_names.push_back("COMP PPP PARTMOBILE sum ranks of missed 6-10");
		feature_names.push_back("COMP PPP NONMOBILE sum ranks of missed 6-10");
		feature_names.push_back("COMP PPP sum ranks of missed 6-10");

		feature_names.push_back("COMP PPP MOBILE sum ranks of missed 11-15");
		feature_names.push_back("COMP PPP PARTMOBILE sum ranks of missed 11-15");
		feature_names.push_back("COMP PPP NONMOBILE sum ranks of missed 11-15");
		feature_names.push_back("COMP PPP sum ranks of missed 11-15");

		int i;
		for (i=0; i<7; i++)
		{
			char buff[64];
			sprintf(buff,"COMB PPP observed rank of predicted rank %d",i+1);
			feature_names.push_back(buff);
		}

		for (i=0; i<7; i++)
		{
			char buff[64];
			sprintf(buff,"COMB PPP predicted rank of observed rank %d",i+1);
			feature_names.push_back(buff);
		}

		for (i=0; i<9; i++)
		{
			char buff[64];
			sprintf(buff,"COMB PPP rank of missed #%d",2*i+1);
			feature_names.push_back(buff);
		}

		for (i=0; i<7; i++)
		{
			char buff[64];
			sprintf(buff,"COMB PPP delta score #%d",i+1);
			feature_names.push_back(buff);
		}

		for (i=0; i<3; i++)
		{
			char buff[64];
			sprintf(buff,"COMB PPP dot prod pred-obs top %d",(i+1)*15);
			feature_names.push_back(buff);

			sprintf(buff,"COMB PPP dot prod obs-pred top %d",(i+1)*15);
			feature_names.push_back(buff);
		}

	}
}


struct mass_tuple {
	mass_tuple() : idx(-1), mass(NEG_INF), score(NEG_INF) {};
	mass_tuple(int i, mass_t m, float s) : idx(i), mass(m), score(s) {}
	bool operator< (const mass_tuple& other) const
	{
		return (mass<other.mass);
	}
	int idx;
	mass_t mass;
	float score;
};

/************************************************************************************
Calculates the shift and scale values for the score, so the distribution of 
samples that have rank 1 (i.e., the top scoring database miss) have a mean 0 and
variance 1. This will enable a comparison of scores accros different models.
*************************************************************************************/
void DeNovoPartitionModel::set_shifts_and_scales_for_db(Config *config,
														const RankBoostDataset& train_ds , 
														const vector<string>& peptide_strings)
{
	cout << "Setting shift and scale parameters.. " << endl;

	const vector<RankBoostSample>& samples = train_ds.get_samples();

	// sort sample indices according to mass
	vector<mass_tuple> bad_tuples,good_tuples;
	bad_tuples.clear();
	good_tuples.clear();

	
	int prev_group_idx=NEG_INF;
	int i;
	for (i=0; i<samples.size(); i++)
	{
		if (samples[i].rank_in_group>0)
		{
			cout << "Problem with sample " << i << ", expecting rank 0!" << endl;
			exit(1);
		}

		score_t score = boost_model.calc_rank_score(samples[i]);
		mass_tuple mt(i,samples[i].float_tag1,score);
		good_tuples.push_back(mt);

		vector<mass_t> good_cuts;

		Peptide good_pep;
		good_pep.parseFromString(config,peptide_strings[i]);
		good_pep.calc_expected_breakage_masses(config,good_cuts);

		samples[i];
		const int good_idx  = i;
		const int groupIndex = samples[i].groupIndex;
		mass_tuple max_mt;

		max_mt.score = NEG_INF;
		i++;
		while (i<samples.size() && samples[i].groupIndex == groupIndex)
		{
			if (samples[i].tag2 != SOL_INCORRECT_DB) // only use incorrect db samples for bad
			{
				i++;
				continue;
			}

			Peptide bad_pep;
			bad_pep.parseFromString(config,peptide_strings[i]);
			vector<mass_t> bad_cuts;
			bad_pep.calc_expected_breakage_masses(config,bad_cuts);

			if (fabs(bad_cuts[bad_cuts.size()-1]-good_cuts[good_cuts.size()-1])>10)
			{
				cout << "Mismatch with peptide strings in recalibration!" << endl;
				exit(1);
			}

			if (compare_cut_lists(config->getTolerance(),good_cuts,bad_cuts))
			{
				i++;
				continue;
			}

			score_t score = boost_model.calc_rank_score(samples[i]);
			if (score>max_mt.score)
			{
				max_mt.score = score;
				max_mt.mass  = samples[i].float_tag1;
				max_mt.idx   = i;
			}
			i++;
		}

		if (max_mt.score>NEG_INF)
			bad_tuples.push_back(max_mt);

		if (i<samples.size())
			i--;
	}


	sort(bad_tuples.begin(),bad_tuples.end());
	sort(good_tuples.begin(),good_tuples.end());

	cout << "Got: " << good_tuples.size() << " good, and " << bad_tuples.size() << " bad scores..." << endl;

	mass_t increment = 100;
	if (good_tuples[0].mass>1500)
		increment = 200;
	if (good_tuples[0].mass>2000)
		increment = 500;


	vector<mass_t> threshes;
	mass_t next_mass = good_tuples[0].mass+increment;
	int last_idx=0;
	for (i=0; i<good_tuples.size(); i++)
	{
		if (good_tuples[i].mass>next_mass && i-last_idx>300)
		{
			threshes.push_back(good_tuples[i].mass);
			next_mass = good_tuples[i].mass + increment;
			last_idx = i;
		}	
	}

	if (threshes.size()>0)
		threshes.pop_back();

	threshes.push_back(POS_INF);


	cout << "Scaling scores..." << endl;
	scaling_factors.resize(threshes.size());

	int good_idx=0;
	int badIndex=0;
	int last_bad_idx=0;
	int r;
	for (r=0; r<threshes.size(); r++)
	{
		vector<score_t> good_scores;
		vector<score_t> bad_scores;

		while (good_idx<good_tuples.size() && good_tuples[good_idx].mass<threshes[r])
		{
			good_scores.push_back(good_tuples[good_idx].score);
			good_idx++;
		}

		while (badIndex<bad_tuples.size() && bad_tuples[badIndex].mass<threshes[r])
		{
			bad_scores.push_back(bad_tuples[badIndex].score);
			badIndex++;
		}
	
		scaling_factors[r].max_pm_with_19 = threshes[r];
	
		score_t good_mean=0, good_sd=1.0;
		score_t bad_mean=0,  bad_sd=1.0;
		calc_mean_sd(good_scores,&good_mean,&good_sd);
		calc_mean_sd(bad_scores,&bad_mean,&bad_sd);
		
		cout << setprecision(3) << fixed;
		cout << endl << endl << "Scaling factor : " << r << endl;
		cout << "For max mass " << threshes[r] << " (" << badIndex - last_bad_idx << " scores)" << endl;
		cout << "Good, before, mean= " << fixed << setprecision(3) << good_mean << "  sd= " << good_sd << endl;
		cout << "Bad,  before, mean= " << fixed << setprecision(3) << bad_mean << "  sd= " << bad_sd << endl << endl;
	
		last_bad_idx=badIndex;

		const mass_t score_shift = -bad_mean;
		const mass_t score_scale = 1.0 / bad_sd;
		scaling_factors[r].score_shift = score_shift;
		scaling_factors[r].score_scale = score_scale;

		for (i=0; i<good_scores.size(); i++)
		{
			good_scores[i] += score_shift;
			good_scores[i] *= score_scale;
		}

		for (i=0; i<bad_scores.size(); i++)
		{
			bad_scores[i] += score_shift;
			bad_scores[i] *= score_scale;
		}


		calc_mean_sd(good_scores,&good_mean,&good_sd);
		calc_mean_sd(bad_scores,&bad_mean,&bad_sd);

		cout << "Good, after, mean= " << fixed << setprecision(3) << good_mean << "  sd= " << good_sd << endl;
		cout << "Bad,  after, mean= " << fixed << setprecision(3) << bad_mean << "  sd= " << bad_sd << endl;

	}
}



void DeNovoPartitionModel::write_denovo_partition_model_header_to_strings(int model_type,
														vector<string>& header_strings) const
{
	header_strings.clear();

	ostringstream oss;
	
	oss << charge << " " << size_idx << endl;
	header_strings.push_back(oss.str()); oss.str("");
	oss << scaling_factors.size() << " " << fixed << setprecision(3);
	int i;
	for (i=0; i<scaling_factors.size(); i++)
		oss << " " << scaling_factors[i].score_shift << " " << scaling_factors[i].score_scale << 
			   " " <<scaling_factors[i].max_pm_with_19;
	oss << endl;
	header_strings.push_back(oss.str()); oss.str("");

	oss << ppp_frag_type_idxs.size();
	int f;
	for (f=0; f<this->ppp_frag_type_idxs.size(); f++)
		oss << " " << ppp_frag_type_idxs[f];
	oss << endl;
	header_strings.push_back(oss.str()); oss.str("");

	oss << use_PTM_peak_features << " " << PTM_peak_start_idx << endl;
	header_strings.push_back(oss.str()); oss.str("");

	oss << use_tryp_terminal_features << " " << tryp_terminal_start_idx << endl;
	header_strings.push_back(oss.str()); oss.str("");

	oss << use_ann_peak_features << " " << ann_peak_start_idx << endl;
	header_strings.push_back(oss.str()); oss.str("");

	oss << use_inten_balance_features << " " << inten_balance_start_idx << endl;
	header_strings.push_back(oss.str()); oss.str("");

	oss << use_peak_offset_features << " " << peak_offset_start_idx << endl;
	header_strings.push_back(oss.str()); oss.str("");

	oss << use_comp_features << " " << comp_start_idx << endl;
	header_strings.push_back(oss.str()); oss.str("");

	oss << use_pmc_features << " " <<  pmc_start_idx << endl;
	header_strings.push_back(oss.str()); oss.str("");

	oss << use_prm_features << " " << prm_start_idx << endl;
	header_strings.push_back(oss.str()); oss.str("");

	oss << use_ppp_features << " " << ppp_start_idx << endl;
	header_strings.push_back(oss.str()); oss.str("");

	oss << use_combined_ppp_features << " " << combined_ppp_start_idx << endl;
	header_strings.push_back(oss.str()); oss.str("");
}


void DeNovoPartitionModel::write_denovo_partition_model(int model_type, const char *path)
{
	ofstream ofs(path);
	if (!ofs.is_open() || !ofs.good())
	{
		cout << "Error: couldn't open file for writing: " << path << endl;
		exit(1);
	}
	

	vector<string> headers;
	write_denovo_partition_model_header_to_strings(model_type, headers);

	int i;
	for (i=0; i<headers.size(); i++)
		ofs << headers[i];

	boost_model.write_rankboost_model(ofs);

	ofs.close();
}


bool DeNovoPartitionModel::read_denovo_part_model(const char *path, Config *config)
{
	ifstream ifs(path);
	if (! ifs.is_open() || ! ifs.good())
	{
		return false;
	}

	ifs >> charge >> size_idx;
	int num_scaling_factors;
	ifs >> num_scaling_factors;
	scaling_factors.resize(num_scaling_factors);
	int i;
	for (i=0; i<num_scaling_factors; i++)
	{
		ifs >> scaling_factors[i].score_shift >> scaling_factors[i].score_scale >> scaling_factors[i].max_pm_with_19;
		if (scaling_factors[i].score_scale <0 || scaling_factors[i].score_scale>500 || scaling_factors[i].max_pm_with_19 < 500)
		{
			cout << "Error: bad scaling factor line in " << path << endl;
			exit(0);
		}
	}
	ifs >> num_ppp_frags;
	int f;
	ppp_frag_type_idxs.resize(num_ppp_frags,NEG_INF);
	for (f=0; f<num_ppp_frags; f++)
		ifs >> ppp_frag_type_idxs[f];

	ppp_fragments.resize(num_ppp_frags);
	for (f=0; f<ppp_frag_type_idxs.size(); f++)
		ppp_fragments[f] = config->get_fragment(ppp_frag_type_idxs[f]);

	ifs >> use_PTM_peak_features >> PTM_peak_start_idx;
	ifs >> use_tryp_terminal_features >> tryp_terminal_start_idx;
	ifs >> use_ann_peak_features >> ann_peak_start_idx;
	ifs >> use_inten_balance_features >> inten_balance_start_idx;
	ifs >> use_peak_offset_features >> peak_offset_start_idx;
	ifs >> use_comp_features >> comp_start_idx;
	ifs >> use_pmc_features >> pmc_start_idx;
	ifs >> use_prm_features >> prm_start_idx;
	ifs >> use_ppp_features >> ppp_start_idx;
	ifs >> use_combined_ppp_features >> combined_ppp_start_idx;

	if (! boost_model.read_rankboost_model(ifs))
		return false;

	this->ind_was_initialized = true;

	return true;
}
