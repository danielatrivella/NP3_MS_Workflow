    #include "Fragmentation.h"
#include "Config.h"


const char *neutral_loss_labels[]={"","-H2O","-NH3","-H2OH2O","-H2ONH3","-NH3NH3"};
const mass_t neutral_loss_offsets[]={0,-MASS_H2O,-MASS_NH3,-MASS_H2OH2O,-MASS_H2ONH3,-MASS_NH3NH3};
const int num_neutral_losses = sizeof(neutral_loss_labels)/sizeof(char *);

const char  *prefix_base_labels[]={"a","b","c"};
const mass_t prefix_base_offsets[]={-26.9871,MASS_PROTON,18.0343}; 
const int num_prefix_bases = sizeof(prefix_base_labels)/sizeof(char *);

const char  *suffix_base_labels[]={"x","y","z"};
const mass_t suffix_base_offsets[]={45.9976,MASS_OHHH,2.99};
const int num_suffix_bases = sizeof(suffix_base_labels)/sizeof(char *);


void FragmentType::read_fragment(istream& is)
{
	char line[128];
	is.getline(line,128);
	istringstream iss(line);
	char direction;
	iss >> direction;
	if (direction == 'p')
	{
		orientation = PREFIX;
	}
	else if (direction == 's')
	{
		orientation = SUFFIX;
	}
	else
	{
		cout << "Error reading fragment: " << line << endl;
		exit(1);
	}
	iss >> charge >> offset >> label;
}

void FragmentType::write_fragment(ostream& os) const
{
	if (orientation == PREFIX)
	{
		os << "p ";
	}
	else
		os << "s ";

	os << charge << " ";
	os << fixed << setprecision(5) << offset << " ";
	os << label << endl;
}

/********************************************************
Creates a string label from the fragment information.
*********************************************************/
void FragmentType::make_frag_label(mass_t tolerance)
{
	

	label = (orientation == PREFIX ? "p" : "s");

	if (charge>1)
	{
		ostringstream os;
		os << charge;
		label+= os.str();
	}

	if (offset != 0)
	{
		ostringstream os;
		os << fixed << setprecision(1) << offset;
		if (offset>0)
			label += "+";
		label += os.str();
	}

	// use standard labels for known fragments


	char **base_labels;
	mass_t *base_offsets;
	int num_bases;

	if (orientation == PREFIX)
	{
		base_labels  = (char **)prefix_base_labels;
		base_offsets = (mass_t *)prefix_base_offsets;
		num_bases    = num_prefix_bases;
	}
	else
	{
		base_labels  = (char **)suffix_base_labels;
		base_offsets = (mass_t *)suffix_base_offsets;
		num_bases    = num_suffix_bases;
	}

	// check all possible fragments
	mass_t tol = tolerance;
	
	if (tolerance>=0.5)
	{
		tol = (0.8 * tolerance)/charge;
	}
	else if (tolerance>0.1)
		tol = tolerance / charge;
	
	int b;
	for (b=0; b<num_bases; b++)
	{
		int n;
		for (n=0; n<num_neutral_losses; n++)
		{
			mass_t calc_offset = (base_offsets[b] + (charge - 1)*MASS_PROTON + neutral_loss_offsets[n])/(mass_t)charge;
			if (fabs(calc_offset-offset)<tol)
			{
				label = base_labels[b];
				if (charge>1)
					label+=char('0'+charge);
				
				if (n>0)
					label+=neutral_loss_labels[n];

				offset = (base_offsets[b] + neutral_loss_offsets[n] + (charge - 1)*MASS_PROTON) / charge;
				break;
			}
		}
		if (n<num_neutral_losses)
			break;
	}
}	
	



// removes fragments that appear to be isotopic peaks of previously
// selected fragments such as b+1, y+2 etc.
void FragmentTypeSet::remove_isotopic_fragments(mass_t tolerance, bool verbose)
{
	if (verbose)
		cout << "Removing isotopic fragments..." << endl;

	mass_t smallTolerance = (tolerance < 0.1 ? tolerance : tolerance * 0.5);

	sort_fragments_according_to_probs();

	int f;
	for (f=0; f<this->fragments.size(); f++)
	{
		FragmentType& frag = fragments[f];

		if (frag.label.length()<1)
			frag.make_frag_label(tolerance);

		// don't check known labels
		if (frag.label[0] != 's' &&  frag.label[0] != 'p')
			continue;

		int j;
		for (j=0; j<f; j++)
		{
			FragmentType& previous = fragments[j];

			if (previous.orientation == frag.orientation &&
				previous.charge      == frag.charge)
			{
				mass_t offset_diff = frag.offset - previous.offset;

				offset_diff *= frag.charge;

				int k;
				for (k=1; k<=3; k++)
					if (fabs(offset_diff - k*MASS_ISO)<smallTolerance)
						break;

				if (k<=3) 
				{
					if (verbose)
						cout << " --- removing isotopic frag " << frag.label << endl;
					frag.prob = -1;
					frag.spec_count = 0;
					break;
				}
			}
		}
	}

	sort_fragments_according_to_probs();

	while (fragments.size()>0 && fragments[fragments.size()-1].prob<0)
		fragments.pop_back();
}


void FragmentTypeSet::sort_fragments_according_to_probs()
{
	sort(fragments.begin(),fragments.end());
}

/*************************************************
Outputs labels in single line
**************************************************/
void FragmentTypeSet::output_fragment_labels(ostream& os) const
{
	int i;
	for (i=0; i<fragments.size()-1; i++)
		os << fragments[i].label << " ";
	os << fragments[i].label << endl;
}


void FragmentTypeSet::print() const
{
	int i;
	
	cout << setw(4) << left << " "
		     << setw(10) << left << "Label " 
		//	 << setw(5) << "Pidx"
		//	 <<  "offset"
			 << endl;
	for (i=0; i<fragments.size(); i++)
		cout << setw(4) << left <<i 
		     << setw(10) << left << fragments[i].label
		//	 << setw(5) << fragments[i].parent_frag_idx
		//	 << fragments[i].offset_from_parent_frag
			 << endl;
}






/***************************************************************
For each fragment idx its parent fragment is defined as the
fragment with the highest probability that has the same charge
and oreientation as the fragment in question.
****************************************************************/
void FragmentTypeSet::set_parent_frag_idxs()
{
	int i;

	for (i=0; i<fragments.size(); i++)
	{
		int j;
		for (j=0; j<i; j++)
		{
			if (fragments[j].charge == fragments[i].charge &&
				fragments[j].orientation == fragments[i].orientation)
				break;
		}

		// found parent_frag_idx
		if (j<i)
		{
			fragments[i].parent_frag_idx = j;
			fragments[i].offset_from_parent_frag = fragments[i].offset - fragments[j].offset;
		}
		else
		{
			fragments[i].parent_frag_idx=-1;
			fragments[i].offset_from_parent_frag =0;
		}
	}
}






void RegionalFragments::select_fragments_with_minimum_prob(score_t min_prob, int max_num_frags)
{
//	cout << "select " << frag_type_idxs.size() << " -> ";
	while (frag_type_idxs.size()>0 && frag_probs[frag_type_idxs.size()-1]<min_prob)
	{
		frag_type_idxs.pop_back();
		frag_probs.pop_back();
	}

	if (max_num_frags>0)
	{
		while (frag_type_idxs.size()>max_num_frags)
		{
			frag_type_idxs.pop_back();
			frag_probs.pop_back();
		}
	}
//	cout << frag_type_idxs.size() << endl;
}


void FragmentCombo::print_combo(const Config *config, ostream & os) const
{
	int i;
	for (i=0; i<this->frag_inten_idxs.size(); i++)
		os << config->get_fragment(frag_inten_idxs[i]).label << " ";
	os << "|";
	for (i=0; i<this->frag_no_inten_idxs.size(); i++)
		os << " " << config->get_fragment(frag_no_inten_idxs[i]).label;
	os << endl;
}


void FragmentCombo::read_combo(Config *config, istream& is)
{
	frag_inten_idxs.clear();
	frag_no_inten_idxs.clear();

	vector<string> tokens;
	char buff[128];
	is.getline(buff,128);

	istringstream iss(buff);
	string s;
	while (iss >> s)
		tokens.push_back(s);

	// find pos of separator
	int i;
	for (i=0; i<tokens.size(); i++)
		if (! strcmp(tokens[i].c_str(),"|"))
			break;
	if (i== tokens.size())
	{
		cout << "Error: couldn't find '|' for FragmentCombo : " << buff << endl;
		exit(1);
	}
	const int sep_pos = i;
	
	// get inten idxs
	for (i=0; i<sep_pos; i++)
	{
		int f_idx=config->get_frag_idx_from_label(tokens[i]);
		if (f_idx>=0)
		{
			frag_inten_idxs.push_back(f_idx);
		}
		else
			break;
	}
	
	// get no inten idxs
	for (i=sep_pos+1; i<tokens.size(); i++)
	{
		int f_idx=config->get_frag_idx_from_label(tokens[i]);
		if (f_idx>=0)
		{
			frag_no_inten_idxs.push_back(f_idx);
		}
		else
			break;
	}
}




struct idx_prob {
	bool operator< (const idx_prob& other) const
	{
		return (prob>other.prob);
	}
	int idx;
	score_t prob;
};



void RegionalFragments::sort_by_prob()
{
	vector<idx_prob> ip;

	if (frag_type_idxs.size() != frag_probs.size())
	{
		cout << "Error: probs and frag_idxs mismatch!" << endl;
		exit(1);
	}
	
	ip.resize(frag_type_idxs.size());
	int i;
	for (i=0; i<frag_type_idxs.size(); i++)
	{
		ip[i].idx=frag_type_idxs[i];
		ip[i].prob= frag_probs[i];
	}
	sort(ip.begin(),ip.end());
	for (i=0; i<ip.size(); i++)
	{
		frag_type_idxs[i]=ip[i].idx;
		frag_probs[i]=ip[i].prob;
	}

	strong_frag_type_idxs.clear();
	frag_type_combos.clear();
}


/***************************************************
Sets default fragments for PepNovo charge 2
****************************************************/
void RegionalFragments::init_pepnovo_types(int charge, Config *config)
{
	const char* pep_labels[]={"y","b","a","y-H2O","b-H2O","a-H2O","y-NH3","b-NH3",
                           "a-NH3","y-H2OH2O","b-H2OH2O","y-NH3H2O","b-NH3H2O",
						   "y2","b2","y3","b3","y2-H2O","b2-H2O","y2-NH3","b2-NH3",
						   "y3-H2O","b3-H2O","y3-NH3","b3-NH3"};
	frag_type_idxs.resize(13);
	int i;

	for (i=0; i<13; i++)
	{
		string label(pep_labels[i]);
		frag_type_idxs[i]=config->get_frag_idx_from_label(label);
	}
		
	if (charge>1)
	{
		frag_type_idxs.resize(15);
		for (i=13; i<15; i++)
		{
			string label(pep_labels[i]);
			frag_type_idxs[i]=config->get_frag_idx_from_label(label);
		}
	}

	
	if (charge>2)
	{
		frag_type_idxs.resize(21);
		for (i=15; i<21; i++)
		{
			string label(pep_labels[i]);
			frag_type_idxs[i]=config->get_frag_idx_from_label(label);
		}
	}

//	cout << ">> ";
//	for (i=0; i<frag_type_idxs.size(); i++)
//		cout << frag_type_idxs[i] << " ";
//	cout << endl;

	frag_probs.resize(config->get_all_fragments().size(),0);
}




void RegionalFragments::init_with_all_types(int charge, Config *config)
{
	int i;
	frag_type_idxs.clear();
	this->frag_type_combos.clear();
	this->strong_frag_type_idxs.clear();

	const vector<FragmentType>& all_fragments = config->get_all_fragments();
	for (i=0; i<all_fragments.size(); i++)
		if (all_fragments[i].charge<=charge)
			frag_type_idxs.push_back(i);

	frag_probs.resize(all_fragments.size(),0);
}






struct idx_pair {
	bool operator< (const idx_pair& other) const
	{
		return (prob>other.prob);
	}
	int f_idx;
	score_t prob;
};

/******************************************************************
// makes the order of the fragments in the frag_idxs and frag_probs 
// vectors be in descending probability order
*******************************************************************/
void RegionalFragments::sort_according_to_frag_probs()
{
	int i;
	vector<idx_pair> pairs;
	if (frag_type_idxs.size() != frag_probs.size())
	{
		cout << "Error: mismtach in frag_idxs and frag_probs sizes!" << 
			frag_type_idxs.size() << " <=> " << frag_probs.size() << endl;
		exit(1);
	}

	pairs.resize(frag_type_idxs.size());
	for (i=0; i<frag_type_idxs.size(); i++)
	{
		pairs[i].f_idx=frag_type_idxs[i];
		pairs[i].prob =frag_probs[i];
	}
	sort(pairs.begin(),pairs.end());
	for (i=0; i<pairs.size(); i++)
	{
		frag_type_idxs[i]=pairs[i].f_idx;
		frag_probs[i]=pairs[i].prob;
	}
}





/*********************************************************************
// chooses all fragments that have high enough probability to be strong
**********************************************************************/
void RegionalFragments::select_strong_fragments(Config *config,
												score_t min_prob, 
												int max_num_strong)
{
	const vector<FragmentType>& all_fragments = config->get_all_fragments();
	int i;

	bool got_prefix = false;
	bool got_suffix = false;

	strong_frag_type_idxs.clear();

	for (i=0; i<frag_type_idxs.size() && i <max_num_strong; i++)
	{
		if (i>=2 && frag_probs[i]<min_prob)
			continue;

		strong_frag_type_idxs.push_back(frag_type_idxs[i]);

		if (all_fragments[frag_type_idxs[i]].orientation == PREFIX)
			got_prefix=true;
		if (all_fragments[frag_type_idxs[i]].orientation == SUFFIX)
			got_suffix=true;
	}

	// try and add an additional frag type...
	if (! got_prefix)
	{
		for (i=0; i<frag_type_idxs.size(); i++)
		{
			if (all_fragments[frag_type_idxs[i]].orientation == PREFIX &&
				frag_probs[i]>min_prob*0.5)
			{
				strong_frag_type_idxs.push_back(frag_type_idxs[i]);
				break;
			}
		}
	}

	if (! got_suffix)
	{
		for (i=0; i<frag_type_idxs.size(); i++)
		{
			if (all_fragments[frag_type_idxs[i]].orientation == SUFFIX &&
				frag_probs[i]>min_prob*0.5)
			{
				strong_frag_type_idxs.push_back(frag_type_idxs[i]);
				break;
			}
		}
	}
}


void Config::set_all_regional_fragment_relationships()
{
	int charge;
	for (charge=0; charge<regionalFragmentSets_.size(); charge++)
	{
		int size_idx;
		for (size_idx=0; size_idx<regionalFragmentSets_[charge].size(); size_idx++)
		{
			int region_idx;
			for (region_idx=0; region_idx<regionalFragmentSets_[charge][size_idx].size(); region_idx++)
				if (regionalFragmentSets_[charge][size_idx][region_idx].get_num_fragments()>0)
				{
					regionalFragmentSets_[charge][size_idx][region_idx].set_fragment_relationships(this);

				//	if (charge==2 && size_idx==1 && region_idx ==0)
				//		regional_fragment_sets[charge][size_idx][region_idx].print_fragment_relationships(this);
				}
			
		}
	}
}

void Config::print_all_regional_fragment_relationships() const
{
	int charge;
	for (charge=0; charge<regionalFragmentSets_.size(); charge++)
	{
		int size_idx;
		for (size_idx=0; size_idx<regionalFragmentSets_[charge].size(); size_idx++)
		{
			int region_idx;
			for (region_idx=0; region_idx<regionalFragmentSets_[charge][size_idx].size(); region_idx++)
				if (regionalFragmentSets_[charge][size_idx][region_idx].get_num_fragments()>0)
				{
					cout << "MODEL " << charge << " " << size_idx << " " << region_idx << endl;
					regionalFragmentSets_[charge][size_idx][region_idx].print_fragment_relationships(this);
				}
			
		}
	}
}

/****************************************************************************
This function sets the idxs of parents of the fragments (i.e., they appear in
a higher rank in the table. We distinguish between two cases for the parents:
same charge orientation, and other.
The function also chooses for each strong fragment a 
*****************************************************************************/
void RegionalFragments::set_fragment_relationships(Config *config)
{
	const vector<FragmentType>& all_fragments = config->get_all_fragments();
	const int num_frags = frag_type_idxs.size();

	mirrorFragmentIndexes_.clear();
	parentFragmentIndexes_.clear();
	parentFragmentsWithSameChargeAndOrientation_.clear();
	
	mirrorFragmentIndexes_.resize(num_frags);
	parentFragmentIndexes_.resize(num_frags);
	parentFragmentsWithSameChargeAndOrientation_.resize(num_frags);
	
	int i;
	for (i=0; i<num_frags; i++)
	{
		const int curr_frag_idx = frag_type_idxs[i];
		const FragmentType& curr_frag = all_fragments[curr_frag_idx];
		int j;
		for (j=0; j<i; j++)
		{
			const int other_frag_idx = frag_type_idxs[j];
			const FragmentType& other_frag = all_fragments[other_frag_idx];

			parentFragmentIndexes_[i].push_back(other_frag_idx);

			if (curr_frag.charge == other_frag.charge &&
				curr_frag.orientation == other_frag.orientation)
				parentFragmentsWithSameChargeAndOrientation_[i].push_back(other_frag_idx);
			
			if (curr_frag.orientation != other_frag.orientation)
			{
				const mass_t sum_offset = (curr_frag.offset * curr_frag.charge) + (other_frag.offset * other_frag.charge);
				const int sum_charge = (curr_frag.charge + other_frag.charge);
				if (fabs(sum_offset - MASS_PROTON*(sum_charge-1)-MASS_OHHH)<0.1) 
				{
					mirrorFragmentIndexes_[i].push_back(other_frag_idx);
					mirrorFragmentIndexes_[j].push_back(curr_frag_idx);
				}
			}
		}	
	}
}

void RegionalFragments::print_fragment_relationships(const Config *config) const
{
	const vector<FragmentType>& all_fragments = config->get_all_fragments();
	const int num_frags = frag_type_idxs.size();
	int i;
	for (i=0; i<num_frags; i++)
	{
		const int curr_frag_idx = frag_type_idxs[i];
		const FragmentType& curr_frag = all_fragments[curr_frag_idx];
		int j;

		cout << i <<"\t" << curr_frag.label <<"\tparents  " << parentFragmentIndexes_[i].size() << " ";
		for (j=0; j<parentFragmentIndexes_[i].size(); j++)
			cout << "\t" << all_fragments[parentFragmentIndexes_[i][j]].label;
		cout << endl;
		
		cout <<  "\t\tparents with same " << parentFragmentsWithSameChargeAndOrientation_[i].size() << " ";
		for (j=0; j<parentFragmentsWithSameChargeAndOrientation_[i].size(); j++)
			cout << "\t" << all_fragments[parentFragmentsWithSameChargeAndOrientation_[i][j]].label;
		cout << endl;
		
		cout << "\t\tmirror frags " << mirrorFragmentIndexes_[i].size() << " ";
		for (j=0; j<mirrorFragmentIndexes_[i].size(); j++)
			cout << "\t" << all_fragments[mirrorFragmentIndexes_[i][j]].label;
		cout << endl;

		cout << endl << endl;
	}
}




