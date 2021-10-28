#include "PeptideRankScorer.h"
#include "AllScoreModels.h"
#include "PepNovo_auxfun.h"
#include <iostream>

struct InspectResultsLine {
	bool parse_from_fields(Config *config,
						   const vector<string>& fields);

	string	SpectrumFile;
	int		scan;
	string	Annotation;
	string	Protein;
	int		Charge;
	float	MQScore;
	float	Score;
	int		Length;
	float	TotalPRMScore;
	float	MedianPRMScore;
	float	FractionY;
	float	FractionB;
	float	Intensity;
	int		NTT;
	float	p_value;
	float	F_Score;
	float	DeltaScore;
	float	DeltaScoreOther;
	int		RecordNumber;
	int		DBFilePos;
	int		SpecFilePos;

	int		orgRank;
	int		aaBefore;
	int		aaAfter;
	Peptide	pep;
};

struct ScanCandidateSet {

	bool add_new_line(const InspectResultsLine& res); // adds only if has the same file_idx, scan as others

	void recalbirate_scores(Config *config); // resorts and fills in the delta score fields 
											// according tho the value of the Score field

	void output_to_stream(ostream& os, int num_lines=-1	) const;

	int scan;
	vector<InspectResultsLine> results;
};



bool ScanCandidateSet::add_new_line(const InspectResultsLine& res)
{
	if (results.size() == 0)
	{
		scan = res.scan;
		results.push_back(res);
		results[0].orgRank = 0;
		return true;
	}

	if ( scan != res.scan)
		return false;

	results.push_back(res);
	results[results.size()-1].orgRank = results.size()-1;
	return true;
}


/***********************************************************************
Returns true if one of the mass lists is contained in the other (up to 
the tolerance level);
************************************************************************/
bool compare_cut_lists(const mass_t tolerance,
					   const vector<mass_t>& a_masses,
					   const vector<mass_t>& b_masses)
{
	const int num_a_cuts = a_masses.size();
	const int num_b_cuts = b_masses.size();
	int num_same=0;
	int a=0, b=0;
	
	const int min_needed = (int)(0.8 * (num_a_cuts<num_b_cuts ? num_a_cuts : num_b_cuts));

	while (a<num_a_cuts && b<num_b_cuts)
	{
		const mass_t a_mass = a_masses[a];
		const mass_t b_mass = b_masses[b];
		if (fabs(a_mass-b_mass)<tolerance)
		{
			num_same++;
			a++;
			b++;
			continue;
		}
		if (a_mass<b_mass)
		{
			a++;
			continue;
		}
		b++;
	}
	if (num_same >= min_needed)
		return true;

	// check shifted
	const mass_t shift = a_masses[num_a_cuts-1] - b_masses[num_b_cuts-1];
	a=0;
	b=0;
	num_same=0;
	while (a<num_a_cuts && b<num_b_cuts)
	{
		const mass_t a_mass = a_masses[a];
		const mass_t b_mass = b_masses[b]+shift;
		if (fabs(a_mass-b_mass)<tolerance)
		{
			num_same++;
			a++;
			b++;
			continue;
		}
		if (a_mass<b_mass)
		{
			a++;
			continue;
		}
		b++;
	}
	if (num_same >= min_needed)
		return true;

	
	return false;
}


/************************************************************************
Reorders the set's solutions according to the scores.
Updates the Delta score and Delta score other according to the new score.
*************************************************************************/
void ScanCandidateSet::recalbirate_scores(Config *config)
{
	vector<score_pair> pairs; 
	
	pairs.resize(results.size());
	int i;
	for (i=0; i<pairs.size(); i++)
	{
		pairs[i].idx=i;
		pairs[i].score = results[i].Score;
	}

	sort(pairs.begin(),pairs.end());

	vector<InspectResultsLine> sorted_results;
	sorted_results.resize(results.size());
	for (i=0; i<pairs.size(); i++)
		sorted_results[i]=results[pairs[i].idx];

	results=sorted_results;
	const int num_results = results.size();

	if (num_results==1)
	{
		results[0].DeltaScore=0;
		results[0].DeltaScoreOther=0;
		return;
	}

	vector< vector<mass_t> > exp_cut_masses;
	exp_cut_masses.resize(num_results);
	for (i=0; i<num_results; i++)
		results[i].pep.calc_expected_breakage_masses(config,exp_cut_masses[i]);
	
	results[0].DeltaScore = results[0].Score - results[1].Score;
	results[0].DeltaScoreOther = results[0].DeltaScore;
	for (i=1; i<num_results; i++)
	{
		results[i].DeltaScore = results[i].Score - results[0].Score;
		results[i].DeltaScoreOther = results[i].DeltaScore; // default value
	}

	const mass_t tolerance = config->getTolerance();
	vector<bool> similar_cuts;
	similar_cuts.resize(num_results,true);
	for (i=1; i<num_results; i++)
		similar_cuts[i]=compare_cut_lists(tolerance,exp_cut_masses[0],exp_cut_masses[i]);
	
	// check if we need to correct the delta other
	if (similar_cuts[1] && num_results>2)
	{
		int j;
		for (j=2; j<num_results-1; j++)
			if (! similar_cuts[j])
				break;
		results[0].DeltaScoreOther = results[0].Score - results[j].Score;
	}
	
	// don't change the delta other of the lower ranks even if they are similar to the top scoring one
	// only the top score is what should be considerd

/*	for (i=1; i<num_results; i++)
	{
		if (similar_cuts[i])
		{
			int j;
			for (j=1; j<num_results-1; j++)
			{
				if (j==i)
					continue;
				if (! compare_cut_lists(tolerance,exp_cut_masses[i],exp_cut_masses[j]))
					break;
			}
			results[i].DeltaScoreOther = results[i].Score - results[j].Score;
		}
	}*/
}



bool InspectResultsLine::parse_from_fields(Config *config,
										   const vector<string>& fields)
{
	if (fields.size() != 20)
	{
		cout<< "Error: inspect results line has " << fields.size() << ", expecting 20" << endl;
		exit(1);
	}

	SpectrumFile = fields[0];

	if (sscanf(fields[1].c_str(),"%d",&scan) != 1 ||
		scan<0 || scan>100000000)
		error("scan");

	Annotation = fields[2];
	Protein	   = fields[3];

	if (sscanf(fields[4].c_str(),"%d",&Charge) != 1 ||
		Charge<0 || Charge>20)
		error("Charge");

	if (sscanf(fields[5].c_str(),"%f",&MQScore) != 1 ||
		Score<NEG_INF || Score>POS_INF)
		error("MQScore");
	
	if (sscanf(fields[6].c_str(),"%d",&Length) != 1 ||
		Length<1 || Length>POS_INF)
		error("Length");
	
	if (sscanf(fields[7].c_str(),"%f",&TotalPRMScore) != 1 ||
		TotalPRMScore<NEG_INF || TotalPRMScore>POS_INF)
		error("TotalPRMScore");

	if (sscanf(fields[8].c_str(),"%f",&MedianPRMScore) != 1 ||
		MedianPRMScore<NEG_INF || MedianPRMScore>POS_INF)
		error("MedianPRMScore");

	if (sscanf(fields[9].c_str(),"%f",&FractionY) != 1 ||
		FractionY<0 || FractionY>1000)
		error("FractionY");

	if (sscanf(fields[10].c_str(),"%f",&FractionB) != 1 ||
		FractionB<0 || FractionB>1000)
		error("FractionB");

	if (sscanf(fields[11].c_str(),"%f",&Intensity) != 1 ||
		Intensity<0)
		error("Intensity");

	if (sscanf(fields[12].c_str(),"%d",&NTT) != 1 ||
		NTT<0 || NTT>3)
		error("NTT");

	if (sscanf(fields[13].c_str(),"%f",&p_value) != 1)
		error("p_value");

	if (sscanf(fields[14].c_str(),"%f",&F_Score) != 1)
		error("F_Score");

	if (sscanf(fields[15].c_str(),"%f",&DeltaScore) != 1)
		error("DeltaScore");

	if (sscanf(fields[16].c_str(),"%f",&DeltaScoreOther) != 1)
		error("DeltaScoreOther");

	if (sscanf(fields[17].c_str(),"%d",&RecordNumber) != 1)
		error("RecordNumber");

	if (sscanf(fields[18].c_str(),"%d",&DBFilePos) != 1)
		error("DBFilePos");

	if (sscanf(fields[19].c_str(),"%d",&SpecFilePos) != 1)
		error("SpecFilePos");

	Score = MQScore;

	const vector<int>& char2aa = config->get_char2aa();
	const int ann_length = Annotation.length();

	if ((Annotation[1] != '.') || (Annotation[ann_length-2] != '.'))
	{
		cout << "Error: bad annotation format: " << Annotation << endl;
		cout << "Expecting X.XXXXXXXXX.X" << endl;
		cout << "Ann1   : " << Annotation[1] << endl;
		cout << "Ann n-2: " << Annotation[ann_length-2] << endl;
		exit(1);
	}

//	cout << "|" << Annotation << "|" << endl;
	aaBefore = char2aa[Annotation[0]];
	aaAfter	 = char2aa[Annotation[ann_length-1]];

	pep.parseFromString(config,Annotation.substr(2,ann_length-4));
	
	return true;
}

void ScanCandidateSet::output_to_stream(ostream& os, int num_lines) const
{
	int i;
	for (i=0; i<this->results.size(); i++)
	{
		if (i==num_lines)
			break;

		os << results[i].SpectrumFile << "\t" << results[i].scan << "\t" << results[i].Annotation << "\t" << results[i].Protein << "\t";
		os << results[i].Charge << "\t" << results[i].Score << "\t" << results[i].Length << "\t";
		os << results[i].TotalPRMScore <<"\t" << results[i].MedianPRMScore << "\t" << results[i].FractionY << "\t";
		os << results[i].FractionB << "\t" << results[i].Intensity << "\t" << results[i].NTT << "\t";
		os << results[i].p_value << "\t" << results[i].F_Score << "\t" << results[i].DeltaScore << "\t";
		os << results[i].DeltaScoreOther << "\t" << results[i].RecordNumber << "\t" << results[i].DBFilePos << "\t";
		os << results[i].SpecFilePos << endl;

	}
}


/***************************************************************************************
This function touches up inspect search results by rescoring the sequences returned by
inspect. The function produces a new inspect results file with the scores (and delta scores)
replaced.
****************************************************************************************/
void PeptideRankScorer::rescore_inspect_results(char *spectra_file, 
											   char *inspect_res, 
											   char *new_res_file) const
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(this->allScoreModelsPtr_);
	Config *config = allScoreModels->get_config();

	ifstream org_res(inspect_res);

	if (!  org_res.is_open() || ! org_res.good())
	{
		cout << "Error: couldn't open original inspect results file for reading:" << inspect_res << endl;
		exit(1);
	}

	ofstream new_res(new_res_file);
	if (! new_res.is_open() || ! new_res.good())
	{
			cout << "Error: couldn't open new inspect results file for writing:" << new_res_file << endl;
		exit(1);
	}

	char line_buff[1024];
	org_res.getline(line_buff,1024);

	bool read_line  = true;
	vector<string> field_names;
	if (line_buff[0] != '#')
	{
		read_line = false;
	}
	else
	{
		string header = string(line_buff);
		split_string(header,field_names);

		int i;
		for (i=0; i<field_names.size(); i++)
			cout << i << "\t" << field_names[i] << endl;
	}


	vector<ScanCandidateSet> cand_sets;
	vector<int> scan_mapping;
	cand_sets.clear();
	scan_mapping.resize(100000,-1);
	
	while (! org_res.eof())
	{
		vector<string> fields;

		if (read_line)
		{
			org_res.getline(line_buff,1024);
			if (org_res.gcount() < 5)
				continue;
		}
		else
		{
			read_line = true;
		}

		split_string(line_buff,fields);
		InspectResultsLine res;

		res.parse_from_fields(config,fields);

		if (cand_sets.size()==0 || ! cand_sets[cand_sets.size()-1].add_new_line(res))
		{
			ScanCandidateSet new_set;
			new_set.add_new_line(res);
			
			if (new_set.scan>=scan_mapping.size())
				scan_mapping.resize(2*scan_mapping.size(),-1);

			scan_mapping[new_set.scan]=cand_sets.size();
			cand_sets.push_back(new_set);
		}
	}
	org_res.close();

	cout << "Read results for " << cand_sets.size() << " scans..." << endl;

	SpectraAggregator sa;
	sa.initializeFromSpectraFilePath(spectra_file, config);
	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();

	cout << "Read " <<  sl.getNumHeaders() << " spectra headers..." << endl;

	vector<bool> spectrum_indicators(cand_sets.size(),false);
	size_t num_found =0;
	for (size_t i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
		const int scan_number = header->getScanNumber();
		if (scan_mapping[scan_number]<0)
			continue;

		AnnotatedSpectrum as;
		as.readSpectrum(sa, header);

		spectrum_indicators[scan_mapping[scan_number]]=true;
		num_found++;

		ScanCandidateSet& cand_set = cand_sets[scan_mapping[scan_number]];
		vector<PeptideSolution> peptide_sols(cand_set.results.size());;
		for (size_t j=0; j<cand_set.results.size(); j++)
		{
			InspectResultsLine& inspect_res = cand_set.results[j];
			PeptideSolution& sol = peptide_sols[j];

			sol.pep = inspect_res.pep;
			sol.pm_with_19 = sol.pep.get_mass_with_19();
			sol.charge = inspect_res.Charge;
			sol.reaches_n_terminal = true;
			sol.reaches_c_terminal = true;
		}

		vector<score_pair> scores;
		scoreCompleteSequences(peptide_sols, as, scores);
	//	score_complete_sequences(peptide_sols,ssf,peaks,num_peaks,scores);

		for (size_t j=0; j<scores.size(); j++)
			cand_set.results[j].Score = scores[j].score;

		cand_set.recalbirate_scores(config);

		vector<string> pep_strings;
		pep_strings.resize(scores.size());
		int max_len =0;
		for (size_t j=0; j<cand_set.results.size(); j++)
		{
			pep_strings[j]=cand_set.results[j].pep.as_string(config);
			if (pep_strings[j].length()>max_len)
				max_len = pep_strings[j].length();
		}

		if (1)
		{
			cand_set.output_to_stream(new_res,10);
		}
		else
		{
			for (size_t j=0; j<cand_set.results.size(); j++)
			{
				cout << cand_set.scan << " " << cand_set.results[j].Charge << "\t";

				cout << cand_set.results[j].Protein.substr(0,3) << " " << pep_strings[j];
				if (pep_strings[j].length()<max_len)
				{
					for (size_t k=pep_strings[j].length(); k<max_len; k++)
						cout << " ";
				}
				cout << "\t" << cand_set.results[j].MQScore << "\t" << cand_set.results[j].Score << "\t" <<
				cand_set.results[j].DeltaScore << "\t" << cand_set.results[j].DeltaScoreOther << endl;
			}
			cout << endl;
		}
	}

	if (num_found<cand_sets.size())
	{
		cout << "Warning: found only " << num_found << "/" << cand_sets.size() << " of the scans scored by InsPecT!" << endl;
	}
	else
	{
		cout << "All scored scans found in spectrum file." << endl;
	}
}



/***************************************************************************************
This function touches up inspect search results by rescoring the sequences returned by
inspect. The function produces a new inspect results file with the scores (and delta scores)
replaced.
****************************************************************************************/
void PeptideRankScorer::recalibrate_inspect_delta_scores(char *spectra_file, 
											   char *inspect_res, 
											   char *new_res_file) const
{
/*	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(this->allScoreModelsPtr_);
	Config *config = allScoreModels->get_config();

	ifstream org_res(inspect_res);

	if (!  org_res.is_open() || ! org_res.good())
	{
		cout << "Error: couldn't open original inspect results file for reading:" << inspect_res << endl;
		exit(1);
	}

	ofstream new_res(new_res_file);
	if (! new_res.is_open() || ! new_res.good())
	{
			cout << "Error: couldn't open new inspect results file for writing:" << new_res << endl;
		exit(1);
	}

	char line_buff[1024];
	org_res.getline(line_buff,1024);

	bool read_line  = true;
	vector<string> field_names;
	if (line_buff[0] != '#')
	{
		read_line = false;
	}
	else
	{
		string header = string(line_buff);
		split_string(header,field_names);

		int i;
		for (i=0; i<field_names.size(); i++)
			cout << i << "\t" << field_names[i] << endl;
	}


	vector<ScanCandidateSet> cand_sets;
	vector<int> scan_mapping;
	cand_sets.clear();
	scan_mapping.resize(100000,-1);
	
	while (! org_res.eof())
	{
		vector<string> fields;

		if (read_line)
		{
			org_res.getline(line_buff,1024);
			if (org_res.gcount() < 5)
				continue;
		}
		else
		{
			read_line = true;
		}

		split_string(line_buff,fields);
		InspectResultsLine res;

		res.parse_from_fields(config,fields);

		if (cand_sets.size()==0 || ! cand_sets[cand_sets.size()-1].add_new_line(res))
		{
			ScanCandidateSet new_set;
			new_set.add_new_line(res);
			
			if (new_set.scan>=scan_mapping.size())
				scan_mapping.resize(2*scan_mapping.size(),-1);

			scan_mapping[new_set.scan]=cand_sets.size();
			cand_sets.push_back(new_set);
		}
	}
	org_res.close();

	cout << "Read results for " << cand_sets.size() << " scans..." << endl;

	FileManager fm;
	FileSet     fs;
	fm.init_from_file(config,spectra_file);
	fs.select_all_files(fm);
	const vector<SingleSpectrumFile *>& all_ssfs = fs.get_ssf_pointers();

	cout << "Read " <<  all_ssfs.size() << " spectra headers..." << endl;

	BasicSpecReader bsr;
	QCPeak *peaks = new QCPeak[5000];

	vector<bool> spectrum_indicators;
	spectrum_indicators.resize(cand_sets.size(),false);

	int num_found =0;
	int i;
	for (i=0; i<all_ssfs.size(); i++)
	{
		SingleSpectrumFile *ssf = all_ssfs[i];
		
		const int scan_number = ssf->get_scan();
		if (scan_mapping[scan_number]<0)
			continue;

		const int num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);

		spectrum_indicators[scan_mapping[scan_number]]=true;
		num_found++;

		ScanCandidateSet& cand_set = cand_sets[scan_mapping[scan_number]];
		
		cand_set.recalbirate_scores(config);

		cand_set.output_to_stream(new_res,10);
	}

	if (num_found<cand_sets.size())
	{
		cout << "Warning: found only " << num_found << "/" << cand_sets.size() << " of the scans scored by InsPecT!" << endl;
	}
	else
	{
		cout << "All scored scans found in spectrum file." << endl;
	}


	delete [] peaks; */
}



