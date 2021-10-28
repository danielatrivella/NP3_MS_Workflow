#include "PeakRankStats.h"
#include "AllScoreModels.h"
#include "AnnotatedSpectrum.h"


extern const char* mobility_labels[];

/*
void create_training_sets()
{
	PeakRankModel rm;
	AllScoreModels model;
	vector<FileManager> fms;

	model.read_model("LTQ_LOW_TRYP");
	Config *config= model.get_config();

	config->apply_selected_PTMs("C+57:M+16:Q-17");

	rm.set_mass_detection_defaults();
	rm.set_size_thresholds();

	vector< vector<TrainingPeptide> > all_tps;

	all_tps.resize(4);

	char *names[]={"HEK","Shew","ecoli","SPut","SFri","Root","Dicty"};
//	char *names[]={"ecoli","SPut","SFri"};
//	char *names[]={"ecoli"};
//	char *names[]={"Dicty"};
	const int num_names = sizeof(names)/sizeof(char *);

	
	int f;
	for (f=0; f<num_names; f++)
	{
		int charge;
		for (charge=1; charge<=3; charge++)
		{
			char path[256];
			sprintf(path,"C:\\Work\\msms5\\NewScore\\tps\\%s_98_%d_unique_tps.txt",names[f],charge);
			int before_size = all_tps[charge].size();
			read_training_peptides_from_file(path,all_tps[charge]);
			cout << names[f] << "\tcharge " << charge << ": read " << all_tps[charge].size() - before_size << endl;
		}
	}

//	rm.partition_training_samples(all_tps);
	rm.partition_training_samples(all_tps,"C:\\Work\\msms5\\NewScore\\sams\\tr",
										  "C:\\Work\\msms5\\NewScore\\sams\\ts",
										  4750,
										  0.25);
}
*/

/*
void make_proton_mobility_summary(FileManager& fm, Config *config)
{

	vector< vector<int> > counts;
	size_t i;

	counts.resize(4);
	for (i=0; i<4; i++)
		counts[i].resize(50,0);

	FileSet fs;
	fs.select_all_files(fm);

	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

	
	for (i=0; i<all_ssf.size(); i++)
	{
		const int mobility = get_proton_mobility(all_ssf[i]->peptide,all_ssf[i]->charge);
		const int length = all_ssf[i]->peptide.get_num_aas();

		counts[0][length]++;
		counts[mobility][length]++;
		
	//	cout << all_ssf[i]->charge << ": " << 
	//		all_ssf[i]->peptide.as_string(config) << " " << mobility << endl;
	}

	vector<int> sums;

	sums.resize(4,0);
	cout << "\tMOB\tPART\tNON\tTOTAL" << endl;
	for (i=0; i<counts[0].size(); i++)
	{
		if (counts[0][i]==0)
			continue;

		cout << i << "\t" << counts[1][i] << "\t" << counts[2][i] << "\t" << 
				counts[3][i] << "\t" << counts[0][i] << endl;
		
		size_t j;
		for (j=0; j<4; j++)
			sums[j]+=counts[j][i];
	}

	cout << "--------------------------------------------" << endl;
	cout << "\t " << sums[1] << "\t" << sums[2] << "\t" << sums[3]  << "\t" << sums[0] << endl;
	cout << fixed << setprecision(3) << "\t" << sums[1]/double(sums[0]) << "\t" <<
		sums[2]/double(sums[0]) << "\t" << sums[3]/double(sums[0]) << endl;
	
}
*/

/*
void make_frag_stat_summary(FileManager& fm, Config *config, int charge)
{

	FileSet fs;
	fs.select_all_files(fm);

	const int b_idx = config->get_frag_idx_from_label("b");
	const int y_idx = config->get_frag_idx_from_label("y");

	int num_b=0,num_y=0,total=0;

	while (1)
	{
		SingleSpectrumFile *ssf;
		AnnotatedSpectrum as;
		if (! fs.get_next_spectrum(fm,config,&as,&ssf))
			break;

		if (ssf->charge != charge)
			continue;
		
		as.annotate_spectrum(as.get_true_mass_with_19());

		const int length = as.getPeptide().get_num_aas();
		const vector<Breakage>& breakages = as.get_breakages();

		cout << "b: ";
		int i;
		for (i=1; i<breakages.size()-1; i++)
		{
			const int b_pos = breakages[i].get_position_of_frag_idx(b_idx);
			if (b_pos>=0)
			{
				cout << i << " ";
				num_b++;
			}
			else
				cout << "  ";
		}
		cout << " y: ";
		for (i=1; i<breakages.size()-1; i++)
		{
			const int y_pos = breakages[i].get_position_of_frag_idx(y_idx);
			if (y_pos>=0)
			{
				cout << i << " ";
				num_y++;
			}
			else
				cout << "  ";
		}

		total += length-1;
		cout << endl;
	}

	cout << "Total:" << endl;
	cout << "b: " << num_b << "/" << total << " = " << setprecision(3) << 
		double(num_b)/total << endl;
	cout << "y: " << num_y << "/" << total << " = " << setprecision(3) << 
		double(num_y)/total << endl;

}
*/

/*

void make_min_max_peak_mass_reports(FileManager& fm, Config *config,
									int charge, char *name)
{
	FileSet fs;
	fs.select_all_files(fm);

	const string out_dir = "C:/Work/msms5/NewScore/reports";
	char out_name[256];

	sprintf(out_name,"%s/%s_minmax_rep_%d.txt",out_dir.c_str(),name,charge);
	ofstream out(out_name);

	while (1)
	{
		SingleSpectrumFile *ssf;
		AnnotatedSpectrum as;
		if (! fs.get_next_spectrum(fm,config,&as,&ssf))
			break;

		as.annotate_spectrum(as.get_true_mass_with_19());
		const int num_peaks = as.getNumPeaks();

		const int length = as.getPeptide().get_num_aas();
		const vector<Breakage>& breakages = as.get_breakages();

		mass_t min_peak_mass = as.getPeakMass(0);
		mass_t max_peak_mass = as.getPeakMass(num_peaks-1);
		mass_t pm_with_19 = as.get_true_mass_with_19();

		if (max_peak_mass>pm_with_19)
			max_peak_mass = pm_with_19;

		out << pm_with_19 << "\t" << min_peak_mass << "\t" << max_peak_mass << endl;


	}
	out.close();

}

*/



struct aa_wr {
	bool operator< (const aa_wr& other) const
	{
		return (wr<other.wr);
	}

	int aa;
	double wr;
};


void report_wr_group_mobile(Config *config, 
			   const vector< vector< vector< vector<double> > > >& wr,
			   const vector< vector< vector< vector<int> > > >& counts_wr,
			   const vector< vector< vector< vector<int> > > >& lens_wr,
			   const vector<string>& pos_labels,
			   ostream& os = cout)
{
	const int min_count = 7;
	const int num_groups = wr.size();
	const vector<string>& aa2labels = config->get_aa2label();
	vector< vector< vector< vector<string> > > > table_data; // group / mobility / distance / ordered aas

	table_data.resize(wr.size());
	int g;
	for (g=0; g<num_groups; g++)
	{
		table_data[g].resize(NONMOBILE+1);

		int m;
		for (m=MOBILE; m<=NONMOBILE; m++)
		{
			table_data[g][m].resize(wr[g][m].size());
			int d;
			for (d=0; d<wr[g][m].size(); d++)
			{
				const vector<double>& wr_gmd = wr[g][m][d];
				const vector<int>& counts_wr_gmd = counts_wr[g][m][d];
				const vector<int>& lens_wr_gmd = lens_wr[g][m][d];
				vector<string>& table = table_data[g][m][d];
				vector<aa_wr> order;

				table.clear();
				order.clear();
				int aa;
				for (aa=0; aa<wr_gmd.size(); aa++)
				{
					const int aa_count = counts_wr_gmd[aa];
					if (aa_count>=min_count)
					{
						aa_wr x;
						x.aa=aa;
						x.wr = (wr_gmd[aa] * (double)lens_wr_gmd[aa]) /
							((double)aa_count*aa_count);
						order.push_back(x);
					}
				}
				sort(order.begin(),order.end());

				table.push_back("AA\tAvg.R\tCount ");
				double avg_rank = 0;
				int total_count = 0;
				int j;
				for (j=0; j<order.size(); j++)
				{
					int aa = order[j].aa;
					char line[64];
					sprintf(line,"%s\t%.3f\t%5d ",aa2labels[aa].c_str(),(float)order[j].wr,
						counts_wr_gmd[aa]);
					table.push_back(string(line));
					avg_rank+=order[j].wr * counts_wr_gmd[aa];
					total_count += counts_wr_gmd[aa];
				}
				if (total_count>0)
					avg_rank /= total_count;
				char line[64];
				sprintf(line,"All\t%.3f\t%6d",avg_rank,total_count);
				table.push_back(string(line));
			}
		}
	}

	// output results
	int m;
	for (m=MOBILE; m<=NONMOBILE; m++)
	{
		const int max_d = wr[0][m].size();
		int d;
		for (d=0; d<max_d; d++)
		{
			int max_table_length = table_data[0][m][d].size();
			int g;
			for (g=1; g<num_groups; g++)
				if (max_table_length<table_data[g][m][d].size())
					max_table_length=table_data[g][m][d].size();

			os << endl << mobility_labels[m] << "\t" << pos_labels[d] << endl;
			for (g=0; g<num_groups; g++)
			{
				os << "\t  G" << g << "\t       ";
				if (g<num_groups-1)
					os << "| ";
			}
			os << endl;

			int t;
			for (t=0; t<max_table_length-1; t++)
			{
				int g;
				for (g=0; g<num_groups; g++)
				{
					if (t>=table_data[g][m][d].size()-1)
					{
						os << "\t\t       ";
						if (g<num_groups-1)
							os << "| ";
					}
					else
					{
						os << table_data[g][m][d][t];
						if (g<num_groups-1)
							os << " | ";
					}
				}
				os << endl;
			}
			for (g=0; g<num_groups; g++)
			{
				os << table_data[g][m][d][table_data[g][m][d].size()-1];
				if (g<num_groups-1)
					os << " | ";
			}
			os << endl;
		}
	}
}



struct nc_wr {
	bool operator< (const nc_wr& other) const
	{
		return (wr<other.wr);
	}

	int n_aa,c_aa;
	double wr;
};


void report_wrnc_group_mobile(Config *config, 
			   const vector< vector< vector< vector<double> > > >& wrnc,
			   const vector< vector< vector< vector<int> > > >& counts_wrnc,
			   const vector< vector< vector< vector<int> > > >& lens_wrnc,
			   const vector<string>& pos_labels,
			   ostream& os = cout)
{
	const int min_count = 5;
	const int num_groups = wrnc.size();
	const vector<string>& aa2labels = config->get_aa2label();
	vector< vector< vector<string> > > table_data; // group / mobility / distance / ordered aas

	table_data.resize(wrnc.size());
	int g;
	for (g=0; g<num_groups; g++)
	{
		table_data[g].resize(NONMOBILE+1);

		int m;
		for (m=MOBILE; m<=NONMOBILE; m++)
		{
			const vector< vector<double> >& wrnc_gm = wrnc[g][m];
			const vector< vector<int> >&    counts_wrnc_gm = counts_wrnc[g][m];
			const vector< vector<int> >&    lens_wrnc_gm = lens_wrnc[g][m];
			vector<string>& table = table_data[g][m];
			vector<nc_wr> order;

			table.clear();
			order.clear();

			int naa;
			for (naa=Ala; naa<wrnc_gm.size(); naa++)
			{
				int caa;
				for (caa=Ala; caa<wrnc_gm[naa].size(); caa++)
				{
					
					const int aa_count = counts_wrnc_gm[naa][caa];
					if (aa_count>=min_count)
					{
						nc_wr x;
						x.n_aa=naa;
						x.c_aa=caa;
						x.wr = (wrnc_gm[naa][caa] * (double)lens_wrnc_gm[naa][caa]) /
							((double)aa_count*aa_count);
						order.push_back(x);
					}
				}
			}

			sort(order.begin(),order.end());

			table.push_back("N-C\tAvg.R\tCount ");
			double avg_rank = 0;
			int total_count = 0;
			int j;
			for (j=0; j<order.size(); j++)
			{
				int naa = order[j].n_aa;
				int caa = order[j].c_aa;
				char line[64];
				sprintf(line,"%s-%s\t%.3f\t%5d ",aa2labels[naa].c_str(),aa2labels[caa].c_str(),
					(float)order[j].wr, counts_wrnc_gm[naa][caa]);
				table.push_back(string(line));
				avg_rank+=order[j].wr * counts_wrnc_gm[naa][caa];
				total_count += counts_wrnc_gm[naa][caa];
			}
			if (total_count>0)
				avg_rank /= total_count;
			char line[64];
			sprintf(line,"All\t%.3f\t%6d",avg_rank,total_count);
			table.push_back(string(line));
		}
	}

	// output results
	int m;
	for (m=MOBILE; m<=NONMOBILE; m++)
	{
		int max_table_length = table_data[0][m].size();
		int g;
		for (g=1; g<num_groups; g++)
			if (max_table_length<table_data[g][m].size())
				max_table_length=table_data[g][m].size();

		os << endl << mobility_labels[m] << endl;
		for (g=0; g<num_groups; g++)
		{
			os << "\t  G" << g << "\t       ";
			if (g<num_groups-1)
				os << "| ";
		}
		os << endl;

		int t;
		for (t=0; t<max_table_length-1; t++)
		{
			int g;
			for (g=0; g<num_groups; g++)
			{
				if (t>=table_data[g][m].size()-1)
				{
					os << "\t\t       ";
					if (g<num_groups-1)
						os << "| ";
				}
				else
				{
					os << table_data[g][m][t];
					if (g<num_groups-1)
						os << " | ";
				}
			}
			os << endl;
		}
		for (g=0; g<num_groups; g++)
		{
			os << table_data[g][m][table_data[g][m].size()-1];
			if (g<num_groups-1)
				os << " | ";
		}
		os << endl;
	}

}

/*
void center_cleavage_reports(string frag_label, int charge)
{
	PeakRankModel rm;
	AllScoreModels model;
	vector<FileManager> fms;

	model.read_model("LTQ_LOW_TRYP");
	Config *config= model.get_config();

	config->apply_selected_PTMs("C+57:M+16:Q-17");

	rm.set_mass_detection_defaults();

	const int frag_idx = config->get_frag_idx_from_label(frag_label);

	vector<TrainingPeptide> tps;

	tps.clear();

//	char *names[]={"HEK","Shew","ecoli","SPut","SFri","Root"};
//	char *names[]={"ecoli","SPut","SFri"};
	char *names[]={"Root"};
	const int num_names = sizeof(names)/sizeof(char *);

	cout << "CHARGE : " << charge << endl << endl;

	int f;
	for (f=0; f<num_names; f++)
	{
		char path[256];
		sprintf(path,"C:\\Work\\msms5\\NewScore\\tps\\%s_98_%d_unique_tps.txt",names[f],charge);
		int before_size = tps.size();
		read_training_peptides_from_file(path,tps);
		cout << names[f] << " : read " << tps.size() - before_size << endl;
	}
	
	cout << "Total Read " << tps.size() << " peptides..." << endl << endl;

	cout << "AA Composition: " << endl;
	aa_composition_stats(tps,config);

	cout << "Fragment detection stats: " << endl;
	fragment_detection_stats(tps,config);
	cout << endl;

	const int num_aas = config->get_max_session_aa_idx()+1;

	cout << "Stats for " << frag_label << " charge " << charge << endl;

	const int num_groups = 4;
	vector< vector< vector<int> > > idxs;

	idxs.resize(num_groups);
	
	int g;
	for (g=0; g<num_groups; g++)
	{
		idxs[g].resize(4);

		const int min_length = 8 + 4*g;
		const int max_length = 11 + 4*g;

		int m;
		for (m=MOBILE; m<=NONMOBILE; m++)
		{
			select_training_peptides(tps,idxs[g][m],charge,m,min_length,max_length);
			cout << g << " " << m << " " << idxs[g][m].size() << endl;
		}
	}



	vector< vector< vector< vector<double> > > > all_wr; // n3 n2 n1 c1 c2 c3 / aa
	vector< vector< vector< vector<int> > > >    all_lens_wr;
	vector< vector< vector< vector<int> > >  >   all_counts_wr;
	vector< vector< vector< vector<double> > > > all_wrnc;
	vector< vector< vector< vector<int> > > >    all_counts_wrnc;
	vector< vector< vector< vector<int> > > >    all_lens_wrnc;

	all_wr.resize(num_groups);
	all_lens_wr.resize(num_groups);
	all_counts_wr.resize(num_groups);
	all_wrnc.resize(num_groups);
	all_counts_wrnc.resize(num_groups);
	all_lens_wrnc.resize(num_groups);

	for (g=0; g<num_groups; g++)
	{
		all_wr[g].resize(NONMOBILE+1);
		all_lens_wr[g].resize(NONMOBILE+1);
		all_counts_wr[g].resize(NONMOBILE+1);
		all_wrnc[g].resize(NONMOBILE+1);
		all_counts_wrnc[g].resize(NONMOBILE+1);
		all_lens_wrnc[g].resize(NONMOBILE+1);

		int m;
		for (m=MOBILE; m<=NONMOBILE; m++)
		{

			vector< vector<double> >& wr = all_wr[g][m];
			vector< vector<int> >&    lens_wr = all_lens_wr[g][m];
			vector< vector<int> >&    counts_wr = all_counts_wr[g][m];
			vector< vector<double> >& wrnc      = all_wrnc[g][m];
			vector< vector<int> >&    counts_wrnc = all_counts_wrnc[g][m];
			vector< vector<int> >&    lens_wrnc   = all_lens_wrnc[g][m];

			// init
			wr.resize(6);
			counts_wr.resize(6);
			lens_wr.resize(6);
			int i;
			for (i=0; i<6; i++)
			{
				wr[i].resize(num_aas,0);
				counts_wr[i].resize(num_aas,0);
				lens_wr[i].resize(num_aas,0);
			}

			wrnc.resize(num_aas);
			counts_wrnc.resize(num_aas);
			lens_wrnc.resize(num_aas);
			for (i=0; i<num_aas; i++)
			{
				wrnc[i].resize(num_aas,0);
				counts_wrnc[i].resize(num_aas,0);
				lens_wrnc[i].resize(num_aas,0);
			}
		
			// fill
			for (i=0; i<idxs[g][m].size(); i++)
			{
				const TrainingPeptide& tp = tps[idxs[g][m][i]];

				if (tp.get_frag_idx_pos(frag_idx)<0)
					continue;

				int min_cut=-1,max_cut=-1;
				int half_len = tp.length/2;

				if (tp.length % 2 == 0)
				{
					min_cut = half_len-1;
					max_cut = half_len+1;
				}
				else
				{
					min_cut = half_len;
					max_cut = half_len+1;
				}
				min_cut -= g/2;
				max_cut += g/2;

				vector<int> ranks;
				tp.get_ranks_for_frag_idx(frag_idx,ranks);
				
				int num_in_range = 0;
				int r;
				for (r=1; r<ranks.size(); r++)
					if (ranks[r]>=0)
						num_in_range++;

				int cut_idx;
				for (cut_idx=min_cut; cut_idx<=max_cut; cut_idx++)
				{
					vector<int> tp_aas;
					tp_aas.resize(6);
					int a;
					for (a=0; a<6; a++)
						tp_aas[a] = tp.amino_acids[cut_idx-3+a];

					
					// fill in the weighted ranks for single aa

					const double weighted_cut_rank = ranks[cut_idx]/(double)num_in_range;
					if (weighted_cut_rank<0)
						continue;

					int i;
					for (i=0; i<6; i++)
					{
						const int aa = tp_aas[i];
						const double& wriaa=wr[i][aa];
						const int&    counts_wriaa = counts_wr[i][aa];
						
						wr[i][aa] += weighted_cut_rank;
						counts_wr[i][aa]++;
						lens_wr[i][aa] += num_in_range;
					}

					const int n_aa = tp_aas[2];
					const int c_aa = tp_aas[3];
					wrnc[n_aa][c_aa]+=weighted_cut_rank;
					counts_wrnc[n_aa][c_aa]++;
					lens_wrnc[n_aa][c_aa] += num_in_range;
				}
			}
		}
	}

	vector<string> pos_labels;

	pos_labels.push_back("N-3");
	pos_labels.push_back("N-2");
	pos_labels.push_back("N-1");
	pos_labels.push_back("C+1");
	pos_labels.push_back("C+2");
	pos_labels.push_back("C+3");

	// report results
	report_wr_group_mobile(config,  all_wr,   all_counts_wr,   all_lens_wr,   pos_labels);

	report_wrnc_group_mobile(config,all_wrnc, all_counts_wrnc, all_lens_wrnc, pos_labels);

}
*/

/*
void n_terminal_cleavage_reports(string frag_label, int cut_idx)
{
	PeakRankModel rm;
	AllScoreModels model;
	vector<FileManager> fms;

	model.read_model("LTQ_LOW_TRYP");
	Config *config= model.get_config();

	config->apply_selected_PTMs("C+57:M+16:Q-17");

	rm.set_mass_detection_defaults();

//	fms.resize(4);
//	fms[1].init_from_list_file(config,"C:\\Work\\msms5\\NewScore\\Shew_98_1_unique_mgf_list.txt");
//	fms[2].init_from_list_file(config,"C:\\Work\\msms5\\NewScore\\Shew_98_2_unique_mgf_list.txt");
//	fms[3].init_from_list_file(config,"C:\\Work\\msms5\\NewScore\\Shew_98_3_unique_mgf_list.txt");


	const int frag_idx = config->get_frag_idx_from_label(frag_label);

	vector< vector<TrainingPeptide> > tps;
	tps.resize(4);
	int c;
	for (c=2; c<=2; c++)
	{
//		read_data_into_training_peptides(fms[c],config,rm,tps[c]);
		read_training_peptides_from_file("C:\\Work\\msms5\\NewScore\\tps\\HEK_98_2_unique_tps.txt",tps[c]);
	}


	const int num_aas = config->get_aa2mass().size();
	const int num_groups = 4;
	const int num_pos = cut_idx+2;

	vector< vector< vector<int> > > idxs;

	idxs.resize(num_groups);
	
	int g;
	for (g=0; g<num_groups; g++)
	{
		idxs[g].resize(4);

		const int min_length = 7 + 2*g;
		const int max_length = 8 + 2*g;

		int m;
		for (m=MOBILE; m<=NONMOBILE; m++)
		{
			select_training_peptides(tps[2],idxs[g][m],2,m,min_length,max_length);
			cout << g << " " << m << " " << idxs[g][m].size() << endl;
		}
	}



	vector< vector< vector< vector<double> > > > all_wr;  // g / m / aa_pos [0..cut+2] / aa
	vector< vector< vector< vector<int> > > >    all_lens_wr;
	vector< vector< vector< vector<int> > >  >   all_counts_wr;
	vector< vector< vector< vector<double> > > > all_wrnc;
	vector< vector< vector< vector<int> > > >    all_counts_wrnc;
	vector< vector< vector< vector<int> > > >    all_lens_wrnc;

	all_wr.resize(num_groups);
	all_lens_wr.resize(num_groups);
	all_counts_wr.resize(num_groups);
	all_wrnc.resize(num_groups);
	all_counts_wrnc.resize(num_groups);
	all_lens_wrnc.resize(num_groups);

	for (g=0; g<num_groups; g++)
	{
		all_wr[g].resize(NONMOBILE+1);
		all_lens_wr[g].resize(NONMOBILE+1);
		all_counts_wr[g].resize(NONMOBILE+1);
		all_wrnc[g].resize(NONMOBILE+1);
		all_counts_wrnc[g].resize(NONMOBILE+1);
		all_lens_wrnc[g].resize(NONMOBILE+1);

		int m;
		for (m=MOBILE; m<=NONMOBILE; m++)
		{

			vector< vector<double> >& wr = all_wr[g][m];
			vector< vector<int> >&    lens_wr = all_lens_wr[g][m];
			vector< vector<int> >&    counts_wr = all_counts_wr[g][m];
			vector< vector<double> >& wrnc      = all_wrnc[g][m];
			vector< vector<int> >&    counts_wrnc = all_counts_wrnc[g][m];
			vector< vector<int> >&    lens_wrnc   = all_lens_wrnc[g][m];

			
			// init
			wr.resize(num_pos);
			counts_wr.resize(num_pos);
			lens_wr.resize(num_pos);
			int i;
			for (i=0; i<num_pos; i++)
			{
				wr[i].resize(num_aas,0);
				counts_wr[i].resize(num_aas,0);
				lens_wr[i].resize(num_aas,0);
			}

			wrnc.resize(num_aas);
			counts_wrnc.resize(num_aas);
			lens_wrnc.resize(num_aas);
			for (i=0; i<num_aas; i++)
			{
				wrnc[i].resize(num_aas,0);
				counts_wrnc[i].resize(num_aas,0);
				lens_wrnc[i].resize(num_aas,0);
			}
		
			// fill
			for (i=0; i<idxs[g][m].size(); i++)
			{
				const TrainingPeptide& tp = tps[2][idxs[g][m][i]];

				if (tp.get_frag_idx_pos(frag_idx)<0)
					continue;
		
				vector<int> ranks;
				tp.get_ranks_for_frag_idx(frag_idx,ranks);
				if (ranks[cut_idx]<0) // not in range
					continue;

				// fill in the weighted ranks for single aa
				int num_in_range = 0;
				int r;
				for (r=1; r<ranks.size(); r++)
					if (ranks[r]>=0)
						num_in_range++;

				const double weighted_cut_rank = ranks[cut_idx]/(double)num_in_range;
				if (weighted_cut_rank<0)
					continue;

				int i;
				for (i=0; i<num_pos; i++)
				{
					const int aa = tp.amino_acids[i];
					const double& wriaa=wr[i][aa];
					const int&    counts_wriaa = counts_wr[i][aa];
						
					wr[i][aa] += weighted_cut_rank;
					counts_wr[i][aa]++;
					lens_wr[i][aa] += num_in_range;
				}

				const int n_aa = tp.amino_acids[0];
				const int c_aa = tp.amino_acids[1];
				wrnc[n_aa][c_aa]+=weighted_cut_rank;
				counts_wrnc[n_aa][c_aa]++;
				lens_wrnc[n_aa][c_aa] += num_in_range;
			}
		}
	}

	vector<string> pos_labels;

	int i;
	for (i=0; i<num_pos; i++)
	{
		char label[64];
		if (i>=cut_idx)
		{
			sprintf(label,"Cut %d, AA +%d",cut_idx, i-cut_idx+1);
		}
		else
			sprintf(label,"Cut %d, AA %d",cut_idx, i-cut_idx);
		pos_labels.push_back(label);
	}
	

	// report results
	report_wr_group_mobile(config,  all_wr,   all_counts_wr,   all_lens_wr,   pos_labels);

	report_wrnc_group_mobile(config,all_wrnc, all_counts_wrnc, all_lens_wrnc, pos_labels);

}

*/

/*

void c_terminal_cleavage_reports(string frag_label, int cut_idx_offset)
{
	PeakRankModel rm;
	AllScoreModels model;
	vector<FileManager> fms;

	model.read_model("LTQ_LOW_TRYP");
	Config *config= model.get_config();

	config->apply_selected_PTMs("C+57:M+16:Q-17");

	rm.set_mass_detection_defaults();

	if (cut_idx_offset>-1)
	{
		cout << "Error: cut_idx_offset must be -1 or less (it is from the c-terminal)" << endl;
		exit(1);
	}

//	fms.resize(4);
//	fms[1].init_from_list_file(config,"C:\\Work\\msms5\\NewScore\\Shew_98_1_unique_mgf_list.txt");
//	fms[2].init_from_list_file(config,"C:\\Work\\msms5\\NewScore\\Shew_98_2_unique_mgf_list.txt");
//	fms[3].init_from_list_file(config,"C:\\Work\\msms5\\NewScore\\Shew_98_3_unique_mgf_list.txt");


	const int frag_idx = config->get_frag_idx_from_label(frag_label);

	vector< vector<TrainingPeptide> > tps;
	tps.resize(4);
	int c;
	for (c=2; c<=2; c++)
	{
//		read_data_into_training_peptides(fms[c],config,rm,tps[c]);
		read_training_peptides_from_file("C:\\Work\\msms5\\NewScore\\tps\\HEK_98_2_unique_tps.txt",tps[c]);
	}

	const int num_aas = config->get_aa2mass().size();
	const int num_groups = 4;
	const int num_pos = 2-cut_idx_offset;

	vector< vector< vector<int> > > idxs;

	idxs.resize(num_groups);
	
	int g;
	for (g=0; g<num_groups; g++)
	{
		idxs[g].resize(4);

		const int min_length = 7 + 2*g;
		const int max_length = 8 + 2*g;

		int m;
		for (m=MOBILE; m<=NONMOBILE; m++)
		{
			select_training_peptides(tps[2],idxs[g][m],2,m,min_length,max_length);
			cout << g << " " << m << " " << idxs[g][m].size() << endl;
		}
	}



	vector< vector< vector< vector<double> > > > all_wr;  // g / m / aa_pos [0..cut+2] / aa
	vector< vector< vector< vector<int> > > >    all_lens_wr;
	vector< vector< vector< vector<int> > >  >   all_counts_wr;
	vector< vector< vector< vector<double> > > > all_wrnc;
	vector< vector< vector< vector<int> > > >    all_counts_wrnc;
	vector< vector< vector< vector<int> > > >    all_lens_wrnc;

	all_wr.resize(num_groups);
	all_lens_wr.resize(num_groups);
	all_counts_wr.resize(num_groups);
	all_wrnc.resize(num_groups);
	all_counts_wrnc.resize(num_groups);
	all_lens_wrnc.resize(num_groups);

	for (g=0; g<num_groups; g++)
	{
		all_wr[g].resize(NONMOBILE+1);
		all_lens_wr[g].resize(NONMOBILE+1);
		all_counts_wr[g].resize(NONMOBILE+1);
		all_wrnc[g].resize(NONMOBILE+1);
		all_counts_wrnc[g].resize(NONMOBILE+1);
		all_lens_wrnc[g].resize(NONMOBILE+1);

		int m;
		for (m=MOBILE; m<=NONMOBILE; m++)
		{

			vector< vector<double> >& wr = all_wr[g][m];
			vector< vector<int> >&    lens_wr = all_lens_wr[g][m];
			vector< vector<int> >&    counts_wr = all_counts_wr[g][m];
			vector< vector<double> >& wrnc      = all_wrnc[g][m];
			vector< vector<int> >&    counts_wrnc = all_counts_wrnc[g][m];
			vector< vector<int> >&    lens_wrnc   = all_lens_wrnc[g][m];

			
			// init
			wr.resize(num_pos);
			counts_wr.resize(num_pos);
			lens_wr.resize(num_pos);
			int i;
			for (i=0; i<num_pos; i++)
			{
				wr[i].resize(num_aas,0);
				counts_wr[i].resize(num_aas,0);
				lens_wr[i].resize(num_aas,0);
			}

			wrnc.resize(num_aas);
			counts_wrnc.resize(num_aas);
			lens_wrnc.resize(num_aas);
			for (i=0; i<num_aas; i++)
			{
				wrnc[i].resize(num_aas,0);
				counts_wrnc[i].resize(num_aas,0);
				lens_wrnc[i].resize(num_aas,0);
			}
		
			// fill
			for (i=0; i<idxs[g][m].size(); i++)
			{
				const TrainingPeptide& tp = tps[2][idxs[g][m][i]];

				if (tp.get_frag_idx_pos(frag_idx)<0)
					continue;

				const int cut_idx = tp.amino_acids.size() + cut_idx_offset;
				vector<int> ranks;
				tp.get_ranks_for_frag_idx(frag_idx,ranks);
				if (ranks[cut_idx]<0) // not in range
					continue;

				// fill in the weighted ranks for single aa
				int num_in_range = 0;
				int r;
				for (r=1; r<ranks.size(); r++)
					if (ranks[r]>=0)
						num_in_range++;

				const double weighted_cut_rank = ranks[cut_idx]/(double)num_in_range;
				if (weighted_cut_rank<0)
					continue;
				
				const int num_aa = tp.amino_acids.size();

				int d;
				for (d=0; d<num_pos; d++)
				{
					const int aa_idx = num_aa - num_pos + d;
					const int aa = tp.amino_acids[aa_idx];
					const double& wriaa=wr[d][aa];
					const int&    counts_wriaa = counts_wr[d][aa];
						
					wr[d][aa] += weighted_cut_rank;
					counts_wr[d][aa]++;
					lens_wr[d][aa] += num_in_range;
				}

				const int n_aa = tp.amino_acids[num_aa-2];
				const int c_aa = tp.amino_acids[num_aa-1];
				wrnc[n_aa][c_aa]+=weighted_cut_rank;
				counts_wrnc[n_aa][c_aa]++;
				lens_wrnc[n_aa][c_aa] += num_in_range;
			}
		}
	}

	vector<string> pos_labels;

	int i;
	for (i=0; i<num_pos; i++)
	{
		int p = i-num_pos;
		char label[64];
		if (p>=0)
		{
			sprintf(label,"Cut %d, AA +%d",cut_idx_offset, p+1);
		}
		else
			sprintf(label,"Cut %d, AA %d",cut_idx_offset, p);
		pos_labels.push_back(label);
	}
	

	// report results
	report_wr_group_mobile(config,  all_wr,   all_counts_wr,   all_lens_wr,   pos_labels);

	report_wrnc_group_mobile(config,all_wrnc, all_counts_wrnc, all_lens_wrnc, pos_labels);
}
*/

/*
void proline_cleavage_reports(string frag_label, int charge)
{
	PeakRankModel rm;
	AllScoreModels model;
	vector<FileManager> fms;

	model.read_model("LTQ_LOW_TRYP");
	Config *config= model.get_config();
	config->apply_selected_PTMs("C+57:M+16:Q-17");
	const vector<string>& aa2label = config->get_aa2label();

	rm.set_mass_detection_defaults();

	const int frag_idx = config->get_frag_idx_from_label(frag_label);

	vector<TrainingPeptide> tps;

	tps.clear();

	char *names[]={"HEK","Shew","ecoli","SPut","SFri","Root"};
//	char *names[]={"ecoli","SPut","SFri"};
//	char *names[]={"Root"};
	const int num_names = sizeof(names)/sizeof(char *);

	cout << "CHARGE : " << charge << endl << endl;

	int f;
	for (f=0; f<num_names; f++)
	{
		char path[256];
		sprintf(path,"C:\\Work\\msms5\\NewScore\\tps\\%s_98_%d_unique_tps.txt",names[f],charge);
		int before_size = tps.size();
		read_training_peptides_from_file(path,tps);
		cout << names[f] << " : read " << tps.size() - before_size << endl;
	}
	
	cout << "Total Read " << tps.size() << " peptides..." << endl << endl;

	cout << "AA Composition: " << endl;
	aa_composition_stats(tps,config);

	cout << "Fragment detection stats: " << endl;
	fragment_detection_stats(tps,config);
	cout << endl;

	const int num_aas = config->get_max_session_aa_idx()+1;

	cout << "Stats for " << frag_label << " charge " << charge << endl;

	const int min_length = 8;
	const int max_length = 14;

	int m;
	for (m=MOBILE; m<=NONMOBILE; m++)
	{
		vector<int> idxs;
		select_training_peptides(tps,idxs,charge,m,min_length,max_length);
		cout << "Examining mobility " << m << " " << idxs.size() << endl;
	
		vector<int> good_sam_idxs;  
		vector<int> good_aa_idxs;
		vector<int> bad_sam_idxs;
		vector<int> bad_aa_idxs;
		
		int i;
		for (i=0; i<idxs.size(); i++)
		{
			const TrainingPeptide& tp = tps[idxs[i]];

			if (tp.get_frag_idx_pos(frag_idx)<0)
					continue;

			int aa;
			for (aa=3; aa<tp.length-3; aa++)
			{
				if (tp.amino_acids[aa] != Pro)
					continue;

				vector<int> ranks;
				tp.get_ranks_for_frag_idx(frag_idx,ranks);
				
				if (ranks[aa]<0 || ranks[aa+1]<0)
					continue;

				if (ranks[aa]<ranks[aa+1] && ranks[aa+1]==tp.length)
				{
					good_sam_idxs.push_back(idxs[i]);
					good_aa_idxs.push_back(aa);
				}
				else if (ranks[aa]>ranks[aa+1] && ranks[aa]==tp.length)
				{
					bad_sam_idxs.push_back(idxs[i]);
					bad_aa_idxs.push_back(aa);

					if (1)
					{
						cout << bad_aa_idxs.size() << "\t";
						int j;
						for (j=0; j<aa; j++)
							cout << aa2label[tp.amino_acids[j]];
						cout << " ";
						cout << aa2label[tp.amino_acids[j]];
						for (j=aa+1; j<tp.amino_acids.size(); j++)
							cout << aa2label[tp.amino_acids[j]];

						cout << "  ";
						cout << ranks[aa] << " " << ranks[aa+1] << endl;
					}
				}
			}
		}

		
		cout << "Good samples: " << good_sam_idxs.size() << endl;
		cout << "Bad  samples: " << bad_sam_idxs.size() << endl;
	}
}
*/




