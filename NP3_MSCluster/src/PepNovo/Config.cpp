#include "Config.h"
#include "SpectraList.h"
#include "FragmentSelection.h"
#include "PepNovo_auxfun.h"


void Config::computeSizeThresholds(const SpectraAggregator& sa)
{
	vector< vector<mass_t> > pmsWith19;

	pmsWith19.resize(sa.getMaxSepctrumCharge()+1);
	int badPmCounts=0;

	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();
	

	int i;
	for (i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* const ssh = sl.getSpectrumHeader(i);
		if (ssh->getPmWith19()>0)
		{
			pmsWith19[ssh->getCharge()].push_back(ssh->getPmWith19());
		}
		else
			badPmCounts++;
	}

	const int minCharge = (sa.getMinSpectrumCharge() > 0 ? sa.getMinSpectrumCharge()  : 1);
	int charge;
	for (charge=minCharge; charge<=sa.getMaxSepctrumCharge(); charge++)
		sort(pmsWith19[charge].begin(),pmsWith19[charge].end());

	// report
	cout << "Examining " << sl.getNumHeaders() << " spectra headers:" << endl;
	
	for (charge=minCharge; charge<=sa.getMaxSepctrumCharge(); charge++)
		cout << "Charge " << charge << ", " << pmsWith19[charge].size() << endl;
	if (badPmCounts>0)
		cout << "Warning: " << badPmCounts << " spectra lacked charge/precursor mass information!" << endl;

	// print a histogram of masses vs. charges
	cout << endl << "Histogram of precursor masses for each charge:" << endl;
	cout		 << "----------------------------------------------" << endl;
	cout << "  %";
	for (charge=minCharge; charge<=sa.getMaxSepctrumCharge(); charge++)
		cout << "\t   " << charge;
	
	cout << endl;
	for (i=0; i<20; i++)
	{
		cout << setprecision(2) << fixed << (i+1.0)/20.0;
		for (charge=minCharge; charge<=sa.getMaxSepctrumCharge(); charge++)
		{
			int idx = (i == 19 ? 
					   pmsWith19[charge].size()-1 :
					   (pmsWith19[charge].size() / 20) * (i+1) );
			cout << "\t" << pmsWith19[charge][idx];
		}
		cout << endl;
	}
	cout << endl;
	
	// partition
	massThresholdsForSizes_.resize(sa.getMaxSepctrumCharge()+1);

	cout << endl << "Partitioning:" << endl;
	for (charge = minCharge; charge<=sa.getMaxSepctrumCharge(); charge++)
	{
		const int numSpectra = pmsWith19[charge].size();
		vector<mass_t>& spectraMasses = pmsWith19[charge];

		massThresholdsForSizes_[charge].clear();

		if (numSpectra<5000)
		{
			massThresholdsForSizes_[charge].push_back(POS_INF);
			continue;
		}

		
		if (numSpectra<20000)
		{
			int idx = numSpectra/2;
			massThresholdsForSizes_[charge].push_back(spectraMasses[idx]);
			massThresholdsForSizes_[charge].push_back(POS_INF);
			continue;
		}

		if (numSpectra<50000)
		{
			// use 3 sizes
			
			int idx1 = numSpectra/3;
			int idx2 = idx1 * 2;

			massThresholdsForSizes_[charge].push_back(spectraMasses[idx1]);
			massThresholdsForSizes_[charge].push_back(spectraMasses[idx2]);
			massThresholdsForSizes_[charge].push_back(POS_INF);

			continue;
		}

		// for large datasets of spectra, use the following rule
		const mass_t minPmDiff = 200.0;
		const int    minNumSpectraPerModel = 16000;
		const int	 maxNumSizes = (numSpectra / minNumSpectraPerModel > 5 ?  
									5 : numSpectra / minNumSpectraPerModel);

		const int	expNumSpectraPerModel = numSpectra / maxNumSizes;

		if (maxNumSizes <=1)
		{
			massThresholdsForSizes_[charge].clear();
			massThresholdsForSizes_[charge].push_back(POS_INF);
			return;
		}

		int i=0;
		int prevIndex=0;
		int nextIndex = expNumSpectraPerModel;
		mass_t lastMass =0;

		massThresholdsForSizes_[charge].clear();

		for (i=0; i<numSpectra; i++)
		{
			if (spectraMasses[i]-lastMass > minPmDiff && i>=nextIndex)
			{
				massThresholdsForSizes_[charge].push_back(spectraMasses[i]);
				lastMass = spectraMasses[i];
				prevIndex = i;

				const int numSizes = massThresholdsForSizes_[charge].size();
				if ((numSizes == maxNumSizes - 1) || (numSpectra - i < 1.5*minNumSpectraPerModel))
					break;
				
				nextIndex = i + ((numSpectra - i)/(maxNumSizes - numSizes));
				i = nextIndex;
			}
		}
		massThresholdsForSizes_[charge].push_back(POS_INF);

		const int numSizes = massThresholdsForSizes_[charge].size();
		// round masses to nearest 10 Da
		for (i=0; i<numSizes-1; i++)
			massThresholdsForSizes_[charge][i] = 10.0 * 
				(static_cast<int>(massThresholdsForSizes_[charge][i])/10);

		// make final spectra counts
		vector<int> counts;
		counts.resize(numSizes, 0);
		for (i=0; i<numSpectra; i++)
		{
			int j;
			for (j=0; j<numSizes; j++)
				if (spectraMasses[i]<massThresholdsForSizes_[charge][j])
					break;
			counts[j]++;
		}

		// report:
		cout << "Charge " << charge << ", " << massThresholdsForSizes_[charge].size() <<
			" sizes: " << setprecision(2) << fixed << endl;

		for (i=0; i<numSizes; i++)
			cout <<  i << "\t" << (i>0 ? massThresholdsForSizes_[charge][i-1] : 0) << " to "
				 << massThresholdsForSizes_[charge][i] << "\t( " << counts[i] << " spectra)" << endl;
		cout << endl;
	}
	cout << "Done partitioning according to charge/size." << endl;
}



bool Config::selectFragmentIonTypes(const SpectraAggregator& sa,
							int	  maxNumberFragmentsPerRegion, 
							float minFragmentProbability )
{

	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();


	cout << "Using " << sl.getNumHeaders() << " spectra headers:" << endl;
	const int minCharge = (sa.getMinSpectrumCharge() > 0 ? sa.getMinSpectrumCharge()  : 1);

	int charge;
	for (charge=minCharge; charge<=sa.getMaxSepctrumCharge(); charge++)
		cout << "Charge " << charge << ", " << sa.getNumSpectraWithCharge(charge) << endl;
	cout << endl;

	
	// select potential fragment type using the fragment offset test
	all_fragments.clear_set();
	FragmentTypeSet potentialFragmentTypes;
	createFragmentsAccordingToOffsetCounts(sa, this, potentialFragmentTypes, minFragmentProbability);

	// add these fragments to the existing set
	addFragmentTypes(potentialFragmentTypes);

	// For each charge, find the most abundant fragments 
	for (charge=minCharge; charge<=sa.getMaxSepctrumCharge(); charge++)
	{
		cout << "Generating fragment sets for charge " << charge << endl;

		// check that there is a minimal number of files...
		const int numSpectraWithCharge = sa.getNumSpectraWithCharge(charge);

		if (numSpectraWithCharge < MINIMAL_NUMBER_SPECTRA_FOR_FRAGMENT_SELECTION)
			continue;

		init_regional_fragment_set_defaults(0,charge);

		selectRegionalFragmentSets(sa, this , charge, true);

		applyCuttoffsToRegionalSets(1.1, maxNumberFragmentsPerRegion);

		// select strong, combos...
		int max_num_combos = maxNumberFragmentsPerRegion > 0 ? maxNumberFragmentsPerRegion : 2;
		if (max_num_combos>3)
			max_num_combos = 3;
		
		select_strong_fragments(charge, 0.5, 3);
		selectFragmentCombinations(sa, this, charge, 3);
	}

	string fragments_file = get_resource_dir() + "/" + getModelName() + "_fragments.txt";
	ofstream os(fragments_file.c_str(),ios::out);
	print_fragments(os);
	set_fragments_file(fragments_file);
	os.close();

	string regional_fragment_sets_file = get_resource_dir() + "/" + getModelName() + "_fragment_sets.txt";
	os.open(regional_fragment_sets_file.c_str(),ios::out);
	print_regional_fragment_sets(os);
	set_regional_fragment_sets_file(regional_fragment_sets_file);
	os.close();

	return true;
}




void Config::learnTolerancesFromData(const SpectraAggregator& sa, mass_t initialToleranceEstimate)
{	
	set_all_regional_fragment_relationships();
	selectTwoOverallStrongFragments();

	cout << "Calculating precursor mass tolerance..." << endl;
	mass_t pm_tol = calculatePrecursorMassTolerance(sa, 0.95);

	cout << "Calculating fragment mass tolerance..." << endl;
	mass_t tol    = calculateFragmentMassTolerance(sa , initialToleranceEstimate*1.2, 0.96);

	setPrecursorMassTolerance(pm_tol);

	if (pm_tol<0.075)
	{
		set_need_to_estimate_pm(0);
	}

	if (pm_tol <0.000001)
	{
		pm_tol = tol;
	}

	if (pm_tol<tol)
	{
		set_tolerance(tol+pm_tol);
	}
	else
		set_tolerance(tol);
}


// determines the parent mass tolerance for which *cuttoff_prob* of the abundant fragments
// are caught
mass_t Config::calculatePrecursorMassTolerance(const SpectraAggregator& sa, 
										       float cutoffProbability) const
{
	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();

	vector<float> offsets;
	int i;
	for (i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);

		if (header->getPeptideStr().length() < 3)
			continue;

		if (header->getCharge()<1)
			continue;

		Peptide pep;
		pep.parseFromString(this, header->getPeptideStr());
		mass_t offset =  header->getOriginalPmWith19()  - pep.get_mass_with_19();
		offsets.push_back(offset);
	}

	// find the offset that keeps the desired proportion of fragments
	sort(offsets.begin(),offsets.end());
	int count=0;
	int targetCount = static_cast<int>(((1.0 - cutoffProbability)*offsets.size()));
	int leftSideIndex=0;
	int rightSideIndex=offsets.size()-1;
	mass_t cutoffOffset=-1;
	while (count<targetCount)
	{
		if (fabs(offsets[leftSideIndex])>offsets[rightSideIndex])
		{
			leftSideIndex++;
		}
		else
			rightSideIndex--;

		if (++count >= targetCount)
		{
			if (fabs(offsets[leftSideIndex]) > fabs(offsets[rightSideIndex]))
			{
				cutoffOffset = fabs(offsets[leftSideIndex]); 
			}
			else
				cutoffOffset = fabs(offsets[rightSideIndex]);

			break;
		}
	}

	cout << "Precursor mass offset for " << setprecision(4) << cutoffProbability << " is " << cutoffOffset << endl;
	return cutoffOffset;
}


// determines the tolerance for which *cuttoff_prob* of the abundant fragments
// are caught
mass_t Config::calculateFragmentMassTolerance(const SpectraAggregator& sa, 
											  mass_t maxTolerance,
											  float cutoffProbability) const
{
	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();

	vector<int> fragmentTypeIndexes;
	fragmentTypeIndexes.clear();
	
	if (get_strong_type1_idx()>=0)
		fragmentTypeIndexes.push_back(get_strong_type1_idx());

	if (get_strong_type2_idx()>=0)
		fragmentTypeIndexes.push_back(get_strong_type2_idx());

	if (fragmentTypeIndexes.size()==0)
	{
		cout << endl <<"Warning: no strong fragments selected, using maximal tolerance!!!" << endl << endl;
		return maxTolerance;
	}

	cout << "Using the following fragments to calculate fragment ion tolerance: ";
	int f;
	for (f=0; f<fragmentTypeIndexes.size(); f++)
		cout << "\t" << this->get_fragment(fragmentTypeIndexes[f]).label;
	cout << endl;

	
	vector<float> offsetCountBins;
	vector<float> offsets;

	offsetCountBins.resize(41,0);

	int totalFragmentCount=0;
	
	int i;
	for (i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);

		Spectrum s;
		if (! s.readSpectrum(sa,header))
			continue;

		vector<mass_t> breakageMasses;
		mass_t trueMassWith19;

		s.getPeptide().calc_expected_breakage_masses(this,breakageMasses);
		trueMassWith19 =  s.getPeptide().get_mass_with_19();

		if (breakageMasses.size()<3)
			continue;		
	
		// loop on fragments first, so high count fragments get precedence over
		// low count fragments that are actually due to b/y ions of previous or
		// next amino acids
		int f;
		for (f=0; f<fragmentTypeIndexes.size(); f++)
		{
			const FragmentType& frag = get_fragment(fragmentTypeIndexes[f]);
			int b;

			
			for (b=1; b<breakageMasses.size()-1; b++)
			{
				mass_t breakageMass = breakageMasses[b];

				const mass_t expectedFragmentMass = frag.calc_expected_mass(breakageMass, trueMassWith19);
				const int peakIndex = s.findPeakWithMaxIntensity(expectedFragmentMass, maxTolerance);

				if (peakIndex>=0)
				{
					
					totalFragmentCount++;
					mass_t offset =  s.getPeakMass(peakIndex) - expectedFragmentMass;

					int binIndex = 20 + (int)((offset / maxTolerance)*20);
					if (binIndex<0)
						binIndex=0;
					if (binIndex>40)
						binIndex=40;

					offsetCountBins[binIndex]++;
					offsets.push_back(offset);
				}
			}
		}
	}

	mass_t toleranceIncrement = maxTolerance * 0.05;
	cout << "Tolerance bin histogram: " << endl;
	for (i=0; i<=40; i++)
		cout << setprecision(4) << (i-20)*toleranceIncrement << "\t" << 
			    offsetCountBins[i]/totalFragmentCount << endl;

	// find the offset that keeps the desired proportion of fragments
	sort(offsets.begin(),offsets.end());

	const int targetCount = (int)((1.0 - cutoffProbability)*totalFragmentCount);
	int count=0;
	int leftIndex=0;
	int rightIndex=offsets.size()-1;
	mass_t cutoffOffset=-1;
	while (count<targetCount)
	{
		if (fabs(offsets[leftIndex])>offsets[rightIndex])
		{
			leftIndex++;
		}
		else
			rightIndex--;

		if (++count == targetCount)
		{
			if (fabs(offsets[leftIndex])>fabs(offsets[rightIndex]))
			{
				cutoffOffset = fabs(offsets[leftIndex]); 
			}
			else
				cutoffOffset = fabs(offsets[rightIndex]);

			break;
		}
	}

	cout << "offset for " << cutoffProbability << " of the fragments is " << cutoffOffset << endl;
	return cutoffOffset;
}


void Config::init_with_defaults()
{
	int i;

	mass_spec_type=ESI_MASS_SPEC; // default type
	config_file = "";
	model_name = "";

	ind_read_PTM_file = false;

	max_n_term_mod = 0;
	max_c_term_mod = 0;
	min_n_term_mod = 0;
	min_c_term_mod = 0; 

	tolerance = 0.5;
	pm_tolerance = 2.5;

	local_window_size=200;
	max_number_peaks_per_local_window=15;
	number_of_strong_peaks_per_local_window=10;
//	random_prob = 0.1;

	max_combo_mass = 0;

	max_edge_length = 2;

	max_charge_for_size=10;       // the size  of size_thresholds

	set_digest_type(0);

	need_to_estimate_pm = 1;

	use_spectrum_charge = 0;

	use_spectrum_mz = 0;

	filter_flag = 1;

	need_to_normalize = 1;

	itraq_mode = 0;

	terminal_score = 10;

	digest_score = 10;

	forbidden_pair_penalty = 25;
	strong_type1_idx =-1;
	strong_type2_idx =-1;

	
	min_ranges.clear();
	max_ranges.clear();
	min_exclude_range=999999;
	max_exclude_range=NEG_INF;

	init_standard_aas();

	session_aas = standard_aas;

	init_model_size_and_region_thresholds();

	int c;
	for (c=max_charge_for_size; c>=0; c--)
		init_regional_fragment_set_defaults(0,c);

	// these all operate on the original aas
	session_tables.init_for_standard_aas();

	// insert labels of original aas
	label2aa.clear();
	const vector<string>& aa2label = get_aa2label();
	for (i=0; i<aa2label.size(); i++)
		label2aa.insert(STRING2INT_MAP::value_type(aa2label[i],i));

	calc_aa_combo_masses();
	set_aa_variants();
	fill_allowed_double_edges();

	init_allowed_node_masses();

	char PTM_file[256];
	strcpy(PTM_file,get_resource_dir().c_str());
	strcat(PTM_file,"/");
	strcat(PTM_file,"PepNovo_PTMs.txt"); 

	read_PTM_file(PTM_file);

	indWasInitialized_ = true;
//	all_fragments.print();
}





void Config::init_standard_aas()
{
	int i;
	standard_aas.clear();
	for (i=Ala; i<=Val; i++)
		standard_aas.push_back(i);
}


void Config::init_model_size_and_region_thresholds()
{


	massThresholdsForSizes_.clear();
	massThresholdsForSizes_.resize(max_charge_for_size+1);

	// region thresholds
	region_thresholds.resize(2);
	region_thresholds[0]=0.225;
	region_thresholds[1]=0.775;
}



void Config::add_exclude_range(mass_t min_range, mass_t max_range)
{
	if (min_range>=max_range)
		return;

	this->min_ranges.push_back(min_range);
	this->max_ranges.push_back(max_range);
	if (this->min_exclude_range>min_range)
		min_exclude_range = min_range;

	if (this->max_exclude_range<max_range)
		max_exclude_range = max_range;
}


/********************************************************
Sets the diegest type.
Assumes that each digest type has a set of N or C terminal
preferred amino acids
*********************************************************/
void Config::set_digest_type(int type)
{
	digest_type = type;
	n_term_digest_aas.clear();
	c_term_digest_aas.clear();

	if (digest_type == NON_SPECIFIC_DIGEST)
		return;

	if (digest_type == TRYPSIN_DIGEST)
	{
		c_term_digest_aas.push_back(Lys);
		c_term_digest_aas.push_back(Arg);
		return;
	}

	cout << "Error: digest type not supported: " << type << endl;
	exit(1);
}

/*********************************************************
Inits the regional fragment sets (resizes according to the
charges, size_idxs, and region_idxs). If set_type = 0, 
each region is given all fragments that are relevant.
If set type = 0, uses all
**********************************************************/
void Config::init_regional_fragment_set_defaults(int set_type, int charge)
{
	const int num_regions = region_thresholds.size()+1;

	if (max_charge_for_size <= charge)
	{
		max_charge_for_size = charge;
		regionalFragmentSets_.resize(max_charge_for_size+1);
	}

	regionalFragmentSets_[charge].resize(massThresholdsForSizes_[charge].size());
	int j;
	for (j=0; j<massThresholdsForSizes_[charge].size(); j++)
	{
		int k;
		regionalFragmentSets_[charge][j].resize(num_regions);
		for (k=0; k<num_regions; k++)
		{
			if (set_type == 0)
			{
				regionalFragmentSets_[charge][j][k].init_with_all_types(charge,this);
			}
			else
			{
				cout << "Error: bad option for init!"<< endl;
				exit(1);
			}
		}
	}
}



// selects the fragments to be used, uses a cutoff that is X times random prob)
void Config::applyCuttoffsToRegionalSets(score_t min_prob_coef, int maxNumberFragmentsPerRegion)
{
	int c;
	for (c=1; c<regionalFragmentSets_.size(); c++)
	{
		int s;
		for (s=0; s<regionalFragmentSets_[c].size(); s++)
		{
			int r;
			
			for (r=0; r<regionalFragmentSets_[c][s].size(); r++)
			{
				float min_prob = min_prob_coef * get_regional_random_probability(c,s,r);
				regionalFragmentSets_[c][s][r].select_fragments_with_minimum_prob(min_prob,
					maxNumberFragmentsPerRegion);
			}
		}
	}
}


// selects the fragments to be used, uses a cutoff that is X times random prob)
void Config::select_fragments_in_sets(score_t min_prob_coef, int maxNumberFragmentsPerRegion)
{
	int c;
	for (c=1; c<regionalFragmentSets_.size(); c++)
	{
		int s;
		for (s=0; s<regionalFragmentSets_[c].size(); s++)
		{
			int r;
			
			for (r=0; r<regionalFragmentSets_[c][s].size(); r++)
			{
				float min_prob = min_prob_coef * get_regional_random_probability(c,s,r);
				regionalFragmentSets_[c][s][r].select_fragments_with_minimum_prob(min_prob,
					maxNumberFragmentsPerRegion);
			}
		}
	}
}





// For each regional fragments selects all fragments that have a minimal probability
// to be strong.
// also look for the two strongest fragment types and store them in the
// strong_type1_idx and strong_type2_idx variables
void Config::select_strong_fragments(int charge,
									 score_t min_prob, 
									 int max_num_strong, 
									 bool verbose)
{
	const vector<FragmentType>& allFragmentTypes = get_all_fragments();
	vector<float> probCounts;

	
	this->all_strong_fragment_type_idxs.clear();

	int s;
	for (s=0; s<regionalFragmentSets_[charge].size(); s++)
	{
		int r;
		for (r=0; r<regionalFragmentSets_[charge][s].size(); r++)
		{
			regionalFragmentSets_[charge][s][r].select_strong_fragments(this,min_prob,max_num_strong);
			const vector<int>& strongFragmentIndexes = regionalFragmentSets_[charge][s][r].get_strong_frag_type_idxs();

			// see if these frags should be put in the strongest two fragments
			int i;
			for (i=0; i<strongFragmentIndexes.size(); i++)
			{
				int fragmentIndex = strongFragmentIndexes[i];

				if (strong_type1_idx == -1)
				{
					strong_type1_idx = fragmentIndex;
					continue;
				}

				if (strong_type2_idx == -1 &&
					fragmentIndex != strong_type1_idx)
				{
					strong_type2_idx = fragmentIndex;
					continue;
				}

					// make sure that type2 has the lower probability
				if (allFragmentTypes[strong_type1_idx].prob<
					allFragmentTypes[strong_type2_idx].prob)
				{
					int tmp;
					// NP3 GOT bug - it wasnt flipping type2 with type1 when the last had a lower prob
					tmp=strong_type2_idx;
					strong_type2_idx=strong_type1_idx;
					strong_type1_idx=tmp;
				}

				if (allFragmentTypes[fragmentIndex].prob>
					allFragmentTypes[strong_type2_idx].prob)
				{
					strong_type2_idx = fragmentIndex;
			

					// make sure that type2 has the lower probability
					if (allFragmentTypes[strong_type1_idx].prob<
						allFragmentTypes[strong_type2_idx].prob)
					{
						int tmp;
						tmp=strong_type2_idx;
						strong_type2_idx=strong_type1_idx;
						strong_type1_idx=tmp;
					}
				}
			}
		}
	}	
}


void Config::selectTwoOverallStrongFragments()
{
	const vector<FragmentType>& allFragmentTypes = get_all_fragments();
	vector<float> probCounts;

	probCounts.resize(allFragmentTypes.size(),0);

	int charge;
	for (charge=1; charge<regionalFragmentSets_.size(); charge++)
	{
		int s;
		for (s=0; s<regionalFragmentSets_[charge].size(); s++)
		{
			int r;
			for (r=0; r<regionalFragmentSets_[charge][s].size(); r++)
			{
				const vector<int>& strongFragmentIndexes = regionalFragmentSets_[charge][s][r].get_strong_frag_type_idxs();

				// see if these frags should be put in the strong two fragments
				int i;
				for (i=0; i<strongFragmentIndexes.size(); i++)
				{
					const int fragmentIndex = strongFragmentIndexes[i];
					probCounts[fragmentIndex] += regionalFragmentSets_[charge][s][r].get_frag_prob(fragmentIndex);
				}
			}
		}
	}

	int r;
	for (r=1; r<=2; r++)
	{
		float maxFragProb=0.0;
		int	  maxFragIndex = -1;

		int i;
		for (i=0; i<probCounts.size(); i++)
		{
			if (probCounts[i]>maxFragProb)
			{
				maxFragProb = probCounts[i];
				maxFragIndex  = i;
			}
		}

		probCounts[maxFragIndex]=-1.0;
		if (r == 1)
		{
			strong_type1_idx = maxFragIndex;
		}
		else
			strong_type2_idx = maxFragIndex;
	}

}

// returns the region idx for a (breakge) mass
int	 Config::calc_region_idx(mass_t breakageMass, mass_t pm_with_19, int charge,
					 mass_t min_peak_mass, mass_t max_peak_mass) const
{
//	return 0;

	const int n_threshes = region_thresholds.size();
	int i;
	vector<mass_t> threshes;
	threshes.resize(region_thresholds.size());
	for (i=0; i<region_thresholds.size(); i++)
		threshes[i] = pm_with_19 * region_thresholds[i];

	if (charge<=2)
	{
		mass_t alt_first_thresh = (min_peak_mass>0) ? min_peak_mass + 150 : threshes[0];
		mass_t alt_last_thresh =  (max_peak_mass>0) ? max_peak_mass - 150 : threshes[n_threshes-1];
		
		if (breakageMass<threshes[0] || breakageMass<alt_first_thresh)
			return 0;

		if (breakageMass>threshes[n_threshes-1] || breakageMass>alt_last_thresh)
			return n_threshes;


		for (i=1; i<n_threshes; i++)
			if (breakageMass<threshes[i])
				break;
		return i;
	}
	else
	{
		if (breakageMass<threshes[0])
			return 0;

		if (breakageMass>threshes[n_threshes-1])
			return n_threshes;


		for (i=1; i<n_threshes; i++)
			if (breakageMass<threshes[i])
				break;
		return i;
	}
}


void Config::set_size_thresholds_according_to_set_of_masses(int charge,
											vector<mass_t>& spectra_masses)
{
	const int num_spectra = spectra_masses.size();

	if (massThresholdsForSizes_.size()<=charge)
		massThresholdsForSizes_.resize(charge+1);

	massThresholdsForSizes_[charge].clear();

	if (num_spectra<2000)
	{
		massThresholdsForSizes_[charge].push_back(POS_INF);
		return;
	}

	sort(spectra_masses.begin(),spectra_masses.end());

	if (num_spectra<10000)
	{
		int idx = num_spectra/2;
		massThresholdsForSizes_[charge].push_back(spectra_masses[idx]);
		massThresholdsForSizes_[charge].push_back(POS_INF);
		return;
	}

	if (num_spectra<20000)
	{
		// use 3 sizes
		
		int idx1 = num_spectra/3;
		int idx2 = idx1+idx1;

		massThresholdsForSizes_[charge].push_back(spectra_masses[idx1]);
		massThresholdsForSizes_[charge].push_back(spectra_masses[idx2]);
		massThresholdsForSizes_[charge].push_back(POS_INF);

		return;
	}

	// for large datasets of spectra, use the following rule
	const mass_t min_pm_diff = 200.0;
	const int    min_num_spectra_per_model = 20000;
	const int	 num_models = spectra_masses.size()/min_num_spectra_per_model;

	if (num_models <=1)
	{
		massThresholdsForSizes_[charge].clear();
		massThresholdsForSizes_[charge].push_back(POS_INF);
		return;
	}

	int i=0;
	int prev_i=0;
	int next_idx = min_num_spectra_per_model;
	mass_t last_mass =0;

	massThresholdsForSizes_[charge].clear();

	while (i<spectra_masses.size())
	{
		if (spectra_masses[i]-last_mass>min_pm_diff && 
			i>= next_idx)
		{
			massThresholdsForSizes_[charge].push_back(spectra_masses[i]);
			cout << "charge " << charge << "\t size " << massThresholdsForSizes_[charge].size() << " \tmasss " << 
				setprecision(1) << fixed << massThresholdsForSizes_[charge][massThresholdsForSizes_[charge].size()-1] << 
				"\t  (spectra " << i-prev_i << ")" << endl;

			last_mass = spectra_masses[i];
			prev_i = i;

			int md_count = massThresholdsForSizes_[charge].size();
			if (md_count == num_models -1 || spectra_masses.size() - i < 1.5*min_num_spectra_per_model)
				break;
			
			next_idx = i + (spectra_masses.size() - i)/(num_models-md_count);
	
		}
		i++;
	}

	massThresholdsForSizes_[charge].push_back(POS_INF);
	
	cout << "charge " << charge << "\t size " << massThresholdsForSizes_[charge].size() << 
		" \tmasss " << setprecision(1) << fixed << massThresholdsForSizes_[charge][massThresholdsForSizes_[charge].size()-1] <<
		"\t  (spectra " << spectra_masses.size()-prev_i << ")" << endl;
	
}

void Config::print_size_thresholds() const
{
	int charge;
	for (charge=0; charge<massThresholdsForSizes_.size(); charge++)
	{
		cout << charge << "\t" << massThresholdsForSizes_[charge].size();
		int t;
		for (t=0; t<massThresholdsForSizes_[charge].size(); t++)
			cout << "\t" << massThresholdsForSizes_[charge][t];
		cout << endl;
	}
}


int Config::determine_charge(Peptide& pep, mass_t m_over_z) const
{
	int c;
	pep.calc_mass(this);
	mass_t pep_mass_with_19 = pep.get_mass()+MASS_OHHH;

	for (c=1; c<10; c++)
	{
		mass_t calc_mass_with_19 = m_over_z * c - (1.0023 * (c-1));
		if (fabs(calc_mass_with_19- pep_mass_with_19)<10.0)
			return c;
	}
	return -1;
}


int	 Config::get_max_session_aa_idx() const
{
	const vector<int>& aas = get_session_aas();
	int max_aa=-1;
	int i;
	for (i=0; i<aas.size(); i++)
		if (aas[i]>max_aa)
			max_aa=aas[i];

	return max_aa;
}



/***********************************************************
// returns the idx of an aa from its label
// -1 if label is not found
***********************************************************/
int Config::get_aa_from_label(const string& label) const
{
	STRING2INT_MAP::const_iterator iter = label2aa.find(label);

	if (iter == label2aa.end())
		return -1;

	return (*iter).second;
}

int Config::get_aa_with_position_from_label(const string& label, int position) const
{
	STRING2INT_MAP::const_iterator iter;
	for (iter = label2aa.begin();iter != label2aa.end(); iter++)
	{
		int aa_idx = (*iter).second;
		if (session_tables.get_aa_position(aa_idx) != position ||
			label[0] != (*iter).first[0])
			continue;

		if (label == session_tables.get_aa2label(aa_idx))
			return aa_idx;
	}


	return -1;
}


void Config::print_fragments(ostream &os) const
{
	int i;
	os << "#ALL FRAGMENTS " << all_fragments.fragments.size() << endl;
	for (i=0; i<all_fragments.fragments.size(); i++)
		all_fragments.fragments[i].write_fragment(os);
}


void Config::print_regional_fragment_sets(ostream &os) const
{
	os << "#FRAGMENT SETS" << endl;

	int c;
	for (c=regionalFragmentSets_.size()-1; c>=0; c--)
	{
		int s;
		for (s=regionalFragmentSets_[c].size()-1; s>=0; s--)
		{
			int r;
			for (r=regionalFragmentSets_[c][s].size()-1; r>=0; r--)
			{
				int f;
				const vector<int>& fragmentIndexs = regionalFragmentSets_[c][s][r].get_frag_type_idxs();

				if (fragmentIndexs.size()>0)
				{
					os << c << " " << s << " " << r << " " << fragmentIndexs.size() << " " <<
						setprecision(5) << regionalFragmentSets_[c][s][r].get_rand_prob() << endl;
					int i;

					// write fragments
					for (f=0; f<fragmentIndexs.size(); f++)
					{
						os << left << setw(9) << all_fragments.get_fragment(fragmentIndexs[f]).label 
							<< " " << setprecision(4) << regionalFragmentSets_[c][s][r].get_frag_prob(f)<< "  " << endl;
						//	setprecision(5) << all_fragments.get_fragment(fragmentIndexs[f]).offset << endl;
					}

					// write strong
					const vector<int>& strong = regionalFragmentSets_[c][s][r].get_strong_frag_type_idxs();
					os << "strong " << strong.size();
					for (i=0; i<strong.size(); i++)
						os << " " << get_fragment(strong[i]).label;
					os << endl;

					// write combo
					const vector<FragmentCombo>& combos = regionalFragmentSets_[c][s][r].get_frag_type_combos();
					os << "combos " << combos.size() << endl;
					for (i=0; i<combos.size(); i++)
						combos[i].print_combo(this,os);

					os << endl;
				}
			}
		}
	}
}





/************************************************************
Reads the fragments from the stream
*************************************************************/
void Config::read_fragments(istream& is)
{
	int i,num_frags=-1;
	char buff[128];
	is.getline(buff,128);

	if (sscanf(buff,"#ALL FRAGMENTS %d",&num_frags) != 1)
	{
		cout << "Error: bad line in fragments file: " << buff << endl;

		if (strlen(buff)<4)
		{
			cout << "This error can possibly be due to Wnix/Windows issues.";
			cout << "Try running dos2unix on all \".txt\" files in the Models directory (run \"dos2unix Models/*.*\" and \"dos2unix Models/*/*.*\")." << endl;
		}
		exit(1);
	}

	all_fragments.clear_set();
	all_strong_fragment_type_idxs.clear();
	
	for (i=0; i<num_frags; i++)
	{
		FragmentType ft;
		ft.read_fragment(is);
		all_fragments.add_fragment_type(ft);
	}
}

/************************************************************
reads the fragment sets from a stream
*************************************************************/
void Config::read_regional_fragment_sets(istream& is)
{
	int i,num_frags=-1;
	char buff[128];
	is.getline(buff,128);

	
	if (strncmp(buff,"#FRAGMENT SETS",12))
	{
		cout << "Error: bad line in fragments file: " << buff << endl;
		if (strlen(buff)<4)
		{
			cout << "This error can possibly be due to Wnix/Windows issues.";
			cout << "Try running dos2unix on all \".txt\" files in the Models directory (run \"dos2unix Models/*.*\" and \"dos2unix Models/*/*.*\")." << endl;
		}

		exit(1);
	}

	regionalFragmentSets_.clear();

	while ( ! is.eof())
	{
		int charge,size_idx,region_idx,num_fragments;
		is.getline(buff,128);
		if (is.gcount()<2)
			continue;

		// some other parameter set is in this file, put back line and return
		if (buff[0]=='#')
		{
			int i;
			for (i=strlen(buff)-1; i>=0; i--)
				is.putback(buff[i]);
			return;
		}

		float rand_prob=-1;
		if (sscanf(buff,"%d %d %d %d %f",&charge,&size_idx,&region_idx,&num_fragments,
			&rand_prob) != 5)
		{
			cout << "Error: bad line in fragments file: " << buff << endl;
			if (strlen(buff)<4)
			{
				cout << "This error can possibly be due to Wnix/Windows issues.";
				cout << "Try running dos2unix on all \".txt\" files in the Models directory (run \"dos2unix Models/*.*\" and \"dos2unix Models/*/*.*\")." << endl;
			}

			exit(1);
		}
		
		// make sure fragment_sets vectors are large enough
		if (regionalFragmentSets_.size()<charge+1)
			regionalFragmentSets_.resize(charge+1);
		if (regionalFragmentSets_[charge].size()<size_idx+1)
			regionalFragmentSets_[charge].resize(size_idx+1);
		if (regionalFragmentSets_[charge][size_idx].size()<region_idx+1)
			regionalFragmentSets_[charge][size_idx].resize(region_idx+1);
		
		RegionalFragments& rf = regionalFragmentSets_[charge][size_idx][region_idx];

		rf.set_rand_prob(rand_prob);

		rf.frag_type_idxs.clear();

		// read fragments
		for (i=0; i<num_fragments; i++)
		{
			string  label;
			score_t prob = 0;
			is.getline(buff,128);
			istringstream iss(buff);

			iss >> label >> prob;
			
			const int fragmentIndex = this->get_frag_idx_from_label(label);
			if (fragmentIndex<0)
			{
				cout << "Error: unrecognized fragment: " << label << endl;
				exit(1);
			}

			all_fragments.fragments[fragmentIndex].prob += prob;

			rf.frag_type_idxs.push_back(fragmentIndex);
			rf.frag_probs.push_back(prob);
		}
		
		// read strong fragments
		int i;
		is.getline(buff,128);
		if (strncmp(buff,"strong",6))
		{
			cout << "Error: expected to see list of strong frags! : " << buff << endl;
			exit(1);
		}
		int num_strong=0;
		istringstream iss(buff+7);
		iss >> num_strong;
		rf.strong_frag_type_idxs.clear();
		for (i=0; i<num_strong; i++)
		{
			string label;
			iss >> label;
			int f_idx=get_frag_idx_from_label(label);
			if (f_idx>=0)
			{
				rf.strong_frag_type_idxs.push_back(f_idx);
				int j;
				for (j=0; j<all_strong_fragment_type_idxs.size(); j++)
					if (all_strong_fragment_type_idxs[j] == f_idx)
						break;
				if (j==all_strong_fragment_type_idxs.size())
					all_strong_fragment_type_idxs.push_back(f_idx);
			}
			else
				break;	
		}

		// read combos
		is.getline(buff,128);
		if (strncmp(buff,"combos",6))
		{
			cout << "Error: expected to see list of combos! : " << buff << endl;
			exit(1);
		}
		int num_combos;
		istringstream iss2(buff+7);
		iss2 >> num_combos;
		rf.frag_type_combos.clear();
		for (i=0; i<num_combos; i++)
		{
			FragmentCombo fc;
			fc.read_combo(this,is);
			int j;
			int  strong_frag_pos=-1;
			for (j=0; j<fc.frag_inten_idxs.size(); j++)
			{
				int k;
				for (k=0; k<rf.strong_frag_type_idxs.size(); k++)
				{
					if (fc.frag_inten_idxs[j] == rf.strong_frag_type_idxs[k])
					{
						strong_frag_pos = j;
						break;
					}
				}
			}

			if (strong_frag_pos<0)
			{
				cout << "Combo: ";
				fc.print_combo(this);
				cout << "Error: combo must have at least one strong fragment type with intensity!" << endl;
				exit(1);
			}

			// swap to make the strong frag first
			if (strong_frag_pos>0)
			{
				int tmp;
				tmp=fc.frag_inten_idxs[strong_frag_pos];
				fc.frag_inten_idxs[strong_frag_pos] = fc.frag_inten_idxs[0];
				fc.frag_inten_idxs[0] = tmp;
			}
			rf.frag_type_combos.push_back(fc);
		}
	}

	// clone regional fragment sets if needed
	if (regionalFragmentSets_[2].size()>0)
	{
		if (regionalFragmentSets_[1].size() == 0)
			clone_regional_fragment_sets(2,1);

		if (regionalFragmentSets_.size()<=4)
			regionalFragmentSets_.resize(5);

		if (regionalFragmentSets_[3].size() == 0)
			clone_regional_fragment_sets(2,3);

		if (regionalFragmentSets_[4].size() == 0)
			clone_regional_fragment_sets(3,4);
	}



}


void Config::clone_regional_fragment_sets(int source_charge, int target_charge)
{
	if (regionalFragmentSets_.size()<=target_charge)
		regionalFragmentSets_.resize(target_charge+1);

	regionalFragmentSets_[target_charge] = regionalFragmentSets_[source_charge];
}



bool read_mass_type_line(const char* prefix, char *line, mass_t& val)
{
	int len = strlen(prefix);
	if (strncmp(prefix,line,len))
		return false;
	istringstream is(line+len);
	is >> val;
	return true;
}

bool read_score_type_line(const char* prefix, char *line, score_t& val)
{
	int len = strlen(prefix);
	if (strncmp(prefix,line,len))
		return false;
	istringstream is(line+len);
	is >> val;
	return true;
}

/**********************************************************************
// parses a line that is assumed to be from a config file
// all parameters are assumed to start with
// #CONF <PARAMETER_NAME> VALUES
***********************************************************************/
void Config::parse_config_parameter(char *current_line)
{
	char buff[256];
	int number;

	// fragments file
	if ( sscanf(current_line,"#CONF FRAGMENTS_FILE %s",buff) == 1)
	{
		string path = resource_dir + "/" + string(buff);

		fstream fs(path.c_str(),ios::in);
		if (! fs.good())
		{
			cout << "Error openening fragment sets: " << path << endl;
			exit(1);
		}
		this->read_fragments(fs);
		fragments_file=buff;
		return;
	}

	// regional fragment sets
	if ( sscanf(current_line,"#CONF REGIONAL_FRAGMENT_SETS_FILE %s",buff) == 1)
	{
		string path = resource_dir + "/" + string(buff);

		fstream fs(path.c_str(),ios::in);
		if (! fs.good())
		{
			cout << "Error openening fragment sets: " << buff << endl;
			exit(1);
		}
		this->read_regional_fragment_sets(fs);
		regional_fragment_sets_file=buff;
		set_all_regional_fragment_relationships();
		return;
	}


	// mass spec type
	if ( sscanf(current_line,"#CONF MASS_SPEC_TYPE %d",&number) == 1)
	{
		mass_spec_type = number;
		return;
	}

	// digest
	if ( sscanf(current_line,"#CONF DIGEST_TYPE %d",&number) == 1)
	{
		set_digest_type(number);
		return;
	}


	// resource dir
	if ( sscanf(current_line,"#CONF RESOURCE_DIR %s",buff) == 1)
	{
		resource_dir = string(buff);
		return;
	}

	if ( sscanf(current_line,"#CONF MAX_NUMBER_OF_PEAKS_PER_LOCAL_WINDOW %d",&number) == 1)
	{
		max_number_peaks_per_local_window = number;
		return;
	}

	if ( sscanf(current_line,"#CONF NUMBER_OF_STRONG_PEAKS_PER_LOCAL_WINDOW %d",&number) == 1)
	{
		number_of_strong_peaks_per_local_window = number;
		return;
	}

	if (read_mass_type_line("#CONF LOCAL_WINDOW_SIZE",current_line,local_window_size))
		return;

	if (read_mass_type_line("#CONF TOLERANCE",current_line,tolerance))
		return;

	if (read_mass_type_line("#CONF PM_TOLERANCE",current_line,pm_tolerance))
		return;

	if ( sscanf(current_line,"#CONF MAX_EDGE_LENGTH %d",&number) == 1)
	{
		max_edge_length = number;
		return;
	}

	
	if (read_score_type_line("#CONF TERMINAL_SCORE",current_line,terminal_score))
		return;

	if (read_score_type_line("#CONF DIGEST_SCORE",current_line,digest_score))
		return;

	if (read_score_type_line("#CONF FORBIDDEN_PAIR_PENALTY",current_line,forbidden_pair_penalty))
		return;

	if (sscanf(current_line,"#CONF NEED_TO_CORRECT_PM %d",&number) == 1)
	{
		need_to_estimate_pm = number;
		return;
	}


	if (!strncmp("#CONF SIZE_THRESHOLDS",current_line,21))
	{
		istringstream is(current_line+22);
		int i,num_sizes;

		is >> num_sizes;
		massThresholdsForSizes_.resize(num_sizes);
		for (i=0; i<num_sizes; i++)
		{
			int j,size;
			is >> size;
			massThresholdsForSizes_[i].resize(size);
			for (j=0; j<size; j++)
				is >> massThresholdsForSizes_[i][j];
		}
		return;
	}

	if (!strncmp("#CONF REGION_THRESHOLDS",current_line,22))
	{
		istringstream is(current_line+23);
		int i,num_sizes;

		is >> num_sizes;
		region_thresholds.resize(num_sizes); 
		for (i=0; i<num_sizes; i++)
			is >> region_thresholds[i];

		return;
	}

	if ( !strncmp(current_line,"#CONF MASS_EXCLUDE_RANGE",23))
	{
		istringstream is(current_line+24);
		mass_t min_inp_range=-1.0;
		mass_t max_inp_range=-2.0;

		is >> min_inp_range >> max_inp_range;
		if (max_inp_range<min_inp_range)
		{
			cout << "Error: paramter \"#CONF MASS_EXCLUDE_RANGE\" should be followed a pair of numbers: minimal mass range and maximal mass range to exclude!" << endl;
			cout << "BAD Line: " << current_line << endl;
		}

		add_exclude_range(min_inp_range,max_inp_range);
	
		return;
	}

	if (current_line[0] == '#' || strlen(current_line)<3)
		return;

	cout << "Warning: bad line in config file: " << current_line << endl;
}




void Config::read_config(const char* file_name)
{
	config_file = file_name;
	
	init_with_defaults();
	
	string path = resource_dir + "/" + string(file_name);

	fstream fs(path.c_str(),ios::in);

	while (! fs.eof() && fs.good())
	{
		char buff[1024];

		fs.getline(buff,1024);
		if (fs.gcount()<4)
			continue;

		parse_config_parameter(buff);
	}

	indWasInitialized_ = true;
}



void Config::write_config()
{

	if (fragments_file.length()>0)
	{
		fragments_file = model_name + "_fragments.txt";
		string path = resource_dir + "/" + fragments_file;
		ofstream fs(path.c_str(),ios::out);
		print_fragments(fs);
	}

	if (regional_fragment_sets_file.length()>0)
	{
		regional_fragment_sets_file = model_name + "_fragment_sets.txt";
		string path = resource_dir + "/" + regional_fragment_sets_file;
		ofstream rfs(path.c_str(),ios::out);
		print_regional_fragment_sets(rfs);
	}

	ofstream os(config_file.c_str(), ios::out);
	print_config_parameters(os);
}



/****************************************************************
	calclates the masses of the different aa_combos
	combos are sorted aa lists.
	lazy version, hardcoded, more elegant to do recursive way 
	(without real recursion).

  And also finds the maximum combo mass

  Fills in the aa_edge_combos vector which holds all lengths together
*****************************************************************/
void Config::calc_aa_combo_masses()
{
	const vector<int>& session_aas = get_session_aas();
	const vector<mass_t>& aa2mass = get_aa2mass();
	const int num_aas = session_aas.size();
	const int last_aa = session_aas[num_aas-1];


	vector< vector<AA_combo> > aa_combo_by_length; // first dim is edge length
	aa_combo_by_length.resize(max_edge_length+1);

	int i;

	if (max_edge_length > MAX_EDGE_SIZE)
	{
		cout << "Error: code doesn't support edges larger than "<< MAX_EDGE_SIZE << endl;
		exit(0);
	}

	if (max_edge_length>5)
	{
		cout << "Error: code doesn't support such large edges!"<< endl;
		exit(0);
	}

	aa_combo_by_length[1].clear();
	for (i=0; i<num_aas; i++)
	{
		const int aa1 = session_aas[i];
		if (aa1 == Ile)
			continue;

		AA_combo ac;
		ac.amino_acids[0]=aa1;
		ac.num_aa=1;
		ac.total_mass=aa2mass[aa1];
		
		aa_combo_by_length[1].push_back(ac);
	}
	sort(aa_combo_by_length[1].begin(),aa_combo_by_length[1].end());
	

	if (max_edge_length>=2)
	{
		aa_combo_by_length[2].clear();
		for (i=0; i<num_aas; i++)
		{
			int j;
			const int aa1 = session_aas[i];
			const mass_t mass1 = aa2mass[aa1];
			if (aa1 == Ile)
				continue;

			for (j=i; j<num_aas; j++)
			{
				const int aa2 = session_aas[j];
				if (aa2 == Ile)
					continue;
		
				AA_combo ac;
				ac.amino_acids[0]=aa1;
				ac.amino_acids[1]=aa2;
				ac.num_aa=2;
				ac.total_mass+=mass1 + aa2mass[aa2];
				
				aa_combo_by_length[2].push_back(ac);
			}
		}
		sort(aa_combo_by_length[2].begin(),aa_combo_by_length[2].end());
	}

	if (max_edge_length>=3)
	{
		aa_combo_by_length[3].clear();
		for (i=0; i<num_aas; i++)
		{
			int j;
			const int aa1 = session_aas[i];
			const mass_t mass1 = aa2mass[aa1];
			if (aa1 == Ile)
				continue;

			for (j=i; j<num_aas; j++)
			{
				const int aa2 = session_aas[j];
				if (aa2 == Ile)
					continue;

				const mass_t mass2 = mass1 + aa2mass[aa2];
				int k;

				for (k=j; k<num_aas; k++)
				{
					const int aa3 = session_aas[k];
					if (aa3 == Ile)
						continue;

					AA_combo ac;
					ac.amino_acids[0]=aa1;
					ac.amino_acids[1]=aa2;
					ac.amino_acids[2]=aa3;
					ac.num_aa=3;
					ac.total_mass+=mass2 + aa2mass[aa3];
					
					aa_combo_by_length[3].push_back(ac);
				}
			}
		}
		sort(aa_combo_by_length[3].begin(),aa_combo_by_length[3].end());
	}


	if (max_edge_length>=4)
	{
		aa_combo_by_length[4].clear();
		for (i=0; i<num_aas; i++)
		{
			int j;
			const int aa1 = session_aas[i];
			const mass_t mass1 = aa2mass[aa1];
			if (aa1 == Ile)
				continue;

			for (j=i; j<num_aas; j++)
			{
				const int aa2 = session_aas[j];
				if (aa2 == Ile)
					continue;

				const mass_t mass2 = mass1 + aa2mass[aa2];
				int k;

				for (k=j; k<num_aas; k++)
				{
					const int aa3 = session_aas[k];
					if (aa3 == Ile)
						continue;

					const mass_t mass3 = mass2 + aa2mass[aa3];
					int l;

					for (l=k; l<num_aas; l++)
					{
						const int aa4 = session_aas[l];
						if (aa4 == Ile)
							continue;

						AA_combo ac;
						ac.amino_acids[0]=aa1;
						ac.amino_acids[1]=aa2;
						ac.amino_acids[2]=aa3;
						ac.amino_acids[3]=aa4;
						ac.num_aa=4;
						ac.total_mass+=mass3 + aa2mass[aa4];
						
						aa_combo_by_length[4].push_back(ac);
					}
				}
			}
		}
		sort(aa_combo_by_length[4].begin(),aa_combo_by_length[4].end());
	}

	if (max_edge_length>=5)
	{
		aa_combo_by_length[5].clear();
		for (i=0; i<num_aas; i++)
		{
			int j;
			const int aa1 = session_aas[i];
			const mass_t mass1 = aa2mass[aa1];
			if (aa1 == Ile)
				continue;

			for (j=i; j<num_aas; j++)
			{
				const int aa2 = session_aas[j];
				if (aa2 == Ile)
					continue;

				const mass_t mass2 = mass1 + aa2mass[aa2];
				int k;

				for (k=j; k<num_aas; k++)
				{
					const int aa3 = session_aas[k];
					if (aa3 == Ile)
						continue;

					const mass_t mass3 = mass2 + aa2mass[aa3];
					int l;

					for (l=k; l<num_aas; l++)
					{
						const int aa4 = session_aas[l];
						if (aa4 == Ile)
							continue;

						const mass_t mass4 = mass3 + aa2mass[aa4];
						int p;

						for (p=l; p<num_aas; p++)
						{
							const int aa5 = session_aas[p];
							if (aa5 == Ile)
								continue;

							AA_combo ac;
							ac.amino_acids[0]=aa1;
							ac.amino_acids[1]=aa2;
							ac.amino_acids[2]=aa3;
							ac.amino_acids[3]=aa4;
							ac.amino_acids[4]=aa5;
							ac.num_aa=5;
							ac.total_mass+=mass4 + aa2mass[aa5];
						
							aa_combo_by_length[5].push_back(ac);
						}
					}
				}
			}
		}
		sort(aa_combo_by_length[5].begin(),aa_combo_by_length[5].end());
	}
	
	max_combo_mass = aa_combo_by_length[max_edge_length][aa_combo_by_length[max_edge_length].size()-1].total_mass + 1.0;

	// fill in aa_edge_combos
	aa_edge_combos.clear();
	aa_edge_combos.reserve((int)(1.5 * aa_combo_by_length[max_edge_length].size()));
	for (i=1; i<aa_combo_by_length.size(); i++)
	{
		int j;
		for (j=0; j<aa_combo_by_length[i].size(); j++)
			aa_edge_combos.push_back(aa_combo_by_length[i][j]);	
	}

	sort(aa_edge_combos.begin(),aa_edge_combos.end());

	// make edge variant vector and fill in the variant pointers
	// and make combo_idxs_by_length vectors
	combo_idxs_by_length.resize(max_edge_length+1);
	for (i=0; i<=max_edge_length; i++)
		combo_idxs_by_length[i].clear();

	variant_vector.clear();
	variant_vector.reserve(aa_edge_combos.size()*5);

	for (i=0; i<aa_edge_combos.size(); i++)
	{
		AA_combo& combo = aa_edge_combos[i];
		combo_idxs_by_length[combo.num_aa].push_back(i);

		if (combo.num_aa == 1)
		{
			combo.num_variants = 1;
			combo.variant_start_idx = variant_vector.size();

			variant_vector.push_back(1);
			variant_vector.push_back(combo.amino_acids[0]);
			continue;
		}


		if (combo.num_aa == 2)
		{
			if (combo.amino_acids[0]==combo.amino_acids[1])
			{
				combo.num_variants = 1;
				combo.variant_start_idx = variant_vector.size();

				variant_vector.push_back(2);
				variant_vector.push_back(combo.amino_acids[0]);
				variant_vector.push_back(combo.amino_acids[1]);
				continue;
			}
			else
			{
				combo.num_variants = 2;
				combo.variant_start_idx = variant_vector.size();

				variant_vector.push_back(2);
				variant_vector.push_back(combo.amino_acids[0]);
				variant_vector.push_back(combo.amino_acids[1]);
				variant_vector.push_back(2);
				variant_vector.push_back(combo.amino_acids[1]);
				variant_vector.push_back(combo.amino_acids[0]);
				continue;
			}
		}

		// generate variant using permutation generating function
		vector<int> org_perm;
		vector< vector<int> > all_perms;

		int i;
		org_perm.clear();
		for (i=0; i<combo.num_aa; i++)
			org_perm.push_back(combo.amino_acids[i]);

		generate_all_permutations(org_perm,all_perms);
		combo.num_variants = all_perms.size();
		combo.variant_start_idx = variant_vector.size();

		for (i=0; i<all_perms.size(); i++)
		{
			int j;
			variant_vector.push_back(combo.num_aa);
			for (j=0; j<combo.num_aa; j++)
				variant_vector.push_back(all_perms[i][j]);
		}
		continue;
	}

	// fill the combo_start_idxs
	int last_combo_idx = aa_edge_combos.size()-1;
	mass_t largest_combo_mass = aa_edge_combos[last_combo_idx].total_mass;
	int last_mass_idx = (int)(largest_combo_mass + 1.0);

	combo_start_idxs.resize(last_mass_idx+5,-1);

	for (i=0; i<aa_edge_combos.size(); i++)
	{
		int mass_idx = (int)(aa_edge_combos[i].total_mass);
		if (combo_start_idxs[mass_idx]<0)
			combo_start_idxs[mass_idx]=i;
	}

	// fill in combo idxs for entries with -1
	int previous=-1;
	for (i=0; i<combo_start_idxs.size(); i++)
	{
		if (combo_start_idxs[i]<0)
		{
			combo_start_idxs[i]=previous;
		}
		else
			previous = combo_start_idxs[i];
	}

	
}


const int * Config::get_first_variant_ptr(int combo_idx) const
{ 
	return &variant_vector[aa_edge_combos[combo_idx].variant_start_idx]; 
}


// returns a pointer and number of variants which have masses in the given mass ranges
int Config::get_ptrs_for_combos_in_mass_range(mass_t min_mass, mass_t max_mass, 
												   int& num_combos) const
{
	num_combos = 0;
	int min_idx = (int)(min_mass);

	if (min_idx>=combo_start_idxs.size())
		return -1;

	int combo_idx = combo_start_idxs[min_idx];
	while (combo_idx<aa_edge_combos.size())
	{
		if (aa_edge_combos[combo_idx].total_mass>max_mass)
			return -1;

		if (aa_edge_combos[combo_idx].total_mass>=min_mass)
			break;

		combo_idx++;
	}

	if (combo_idx == aa_edge_combos.size())
		return -1;

	int combos_start_idx = combo_idx;
	while (combo_idx<aa_edge_combos.size())
	{
		if (aa_edge_combos[combo_idx].total_mass>max_mass)
			break;

		combo_idx++;
	}

	num_combos = combo_idx-combos_start_idx;
	
	return combos_start_idx;
}

// returns true if there is a combo that contains the ordered variant of the given aas
bool Config::combos_have_variant(const vector<int>& combos, int num_aa, int *var_aas) const
{
	int c;
	for (c=0; c<combos.size(); c++)
	{
		const AA_combo& combo =  aa_edge_combos[combos[c]];
		int pos = combo.variant_start_idx;
		int v;

		for (v=0; v<combo.num_variants; v++)
			if (variant_vector[pos++]==num_aa)
			{
				int a;
				for (a=0; a<num_aa; a++)
					if (variant_vector[pos+a] != var_aas[a])
						break;
				
				if (a==num_aa)
					return true;

				pos += num_aa;
			}
	}
	return false;
}






/**********************************************************
	calculates the aa_variants vectors (for terminalss and
	amino acids Ala-Val
***********************************************************/
void Config::set_aa_variants()
{
	const vector<int>& org_aas = this->get_org_aa();
	int a;
	aa_variants.clear();
	aa_variants.resize(Val+1);

	aa_variants[N_TERM].push_back(N_TERM);
	aa_variants[C_TERM].push_back(C_TERM);
	aa_variants[Xle].push_back(Xle);

	for (a=0; a<session_aas.size(); a++)
	{
		int aa = session_aas[a];
		int org_aa=org_aas[aa];
		aa_variants[org_aa].push_back(aa);
	}
}


/**********************************************************
Fills the array of amino acid combinations that can be
double edeges. Based on Yingying Huang et al  2003
***********************************************************/
void Config::fill_allowed_double_edges(bool allow_all)
{
	int i;
	int max_aa = Val;
	for (i=0; i<session_aas.size(); i++)
		if (session_aas[i]>max_aa)
			max_aa= session_aas[i];

	int num_aa = max_aa+1;

	allowed_double_edge.clear();
	double_edge_with_same_mass_as_single.clear();

	allowed_double_edge.resize(num_aa);
	double_edge_with_same_mass_as_single.resize(num_aa);

	for (i=0; i<num_aa; i++)
	{
		allowed_double_edge[i].resize(num_aa,allow_all);
		double_edge_with_same_mass_as_single[i].resize(num_aa,false);
	}

	if (tolerance<0.05)
		allow_all=true;


	if (! allow_all)
	{
		for (i=Ala; i<=Val; i++)
		{
			allowed_double_edge[Pro][i]=true;
			allowed_double_edge[Gly][i]=true;
			allowed_double_edge[Ser][i]=true;
		}

	//	allowed_double_edge[Gly][Ala]=false;
		allowed_double_edge[Ser][Ser]=false;
		allowed_double_edge[Ser][Tyr]=false;
		allowed_double_edge[Ser][Pro]=false;
		allowed_double_edge[Ser][His]=false;
		allowed_double_edge[Ser][Gly]=false;
		allowed_double_edge[Ser][Val]=false;



		allowed_double_edge[Thr][Val]=true;
		allowed_double_edge[Thr][Gln]=true;
		allowed_double_edge[Thr][His]=true;
		allowed_double_edge[Thr][Glu]=true;
		allowed_double_edge[Thr][Asp]=true;

		allowed_double_edge[Asn][Asp]=true;
		allowed_double_edge[Asn][Glu]=true;
		allowed_double_edge[Asn][His]=true;
		allowed_double_edge[Asn][Ile]=true;
		allowed_double_edge[Asn][Leu]=true;
		allowed_double_edge[Asn][Met]=true;
		allowed_double_edge[Asn][Asn]=true;
		allowed_double_edge[Asn][Gln]=true;
		allowed_double_edge[Asn][Val]=true;
		allowed_double_edge[Asn][Tyr]=true;

		allowed_double_edge[Tyr][Asp]=true;
		allowed_double_edge[Tyr][Glu]=true;
		allowed_double_edge[Tyr][Thr]=true;

		for (i=Ala; i<=Val; i++)
		{
			allowed_double_edge[i][Lys]=true;
			allowed_double_edge[i][Arg]=true;
		}

		allowed_double_edge[Phe][Glu]=true;

		allowed_double_edge[Ala][Leu]=true;
		allowed_double_edge[Leu][Ala]=true;
		allowed_double_edge[Ala][Ala]=true;
		allowed_double_edge[Trp][Glu]=true;

		// add these double edges, give penalty if they are used as a double edge
		allowed_double_edge[Ala][Gly]=true;
		allowed_double_edge[Gly][Ala]=true;
		allowed_double_edge[Gly][Gly]=true;
		allowed_double_edge[Ala][Asp]=true;
		allowed_double_edge[Asp][Ala]=true;
		allowed_double_edge[Val][Ser]=true;
		allowed_double_edge[Ser][Val]=true;
		allowed_double_edge[Gly][Glu]=true;
		allowed_double_edge[Glu][Gly]=true;

	}	

	double_edge_with_same_mass_as_single[Ala][Gly]=true;
	double_edge_with_same_mass_as_single[Gly][Ala]=true;
	double_edge_with_same_mass_as_single[Gly][Gly]=true;
	double_edge_with_same_mass_as_single[Ala][Asp]=true;
	double_edge_with_same_mass_as_single[Asp][Ala]=true;
	double_edge_with_same_mass_as_single[Val][Ser]=true;
	double_edge_with_same_mass_as_single[Ser][Val]=true;
	double_edge_with_same_mass_as_single[Gly][Glu]=true;
	double_edge_with_same_mass_as_single[Glu][Gly]=true;

	// add fix for PTMs
	const vector<int>& org_aas = session_tables.get_org_aa();
	for (i=0; i<session_aas.size(); i++)
	{
		int j;
		for (j=0; j<session_aas.size(); j++)
		{
			int aa1 = session_aas[i];
			int aa2 = session_aas[j];

			if (org_aas[aa1]>=Ala && org_aas[aa2]>=Ala &&
				allowed_double_edge[org_aas[aa1]][org_aas[aa2]])
				allowed_double_edge[aa1][aa2]=true;
		}
	}
}


/************************************************************
Add new fragments if they do not appear in list
*************************************************************/
void Config::addFragmentTypes(const FragmentTypeSet& fts)
{
	int i;
	for (i=0; i<fts.get_fragments().size(); i++)
	{
		all_fragments.add_fragment_type(fts.get_fragment(i));
	}
}



void Config::print_config_parameters(ostream& os) const
{
	int i;

	if (fragments_file.length()>0)
		os << "#CONF FRAGMENTS_FILE " << fragments_file << endl;

	if (regional_fragment_sets_file.length()>0)
		os << "#CONF REGIONAL_FRAGMENT_SETS_FILE " << regional_fragment_sets_file << endl;

	if (aa_combo_file.length()>0)
		os << "#CONF AA_COMBO_FILE " << aa_combo_file << endl;

//	if (resource_dir.length()>0)
//		os << "#CONF RESOURCE_DIR " << resource_dir << endl;


	os << "#CONF MASS_SPEC_TYPE " << mass_spec_type << "(currently only option 1) " << endl;

	os << "#CONF DIGEST_TYPE " << digest_type << " (" << NON_SPECIFIC_DIGEST << " - non specific, " <<
		 TRYPSIN_DIGEST << " - trypsin)" << endl;

	os << "#CONF MAX_NUMBER_OF_PEAKS_PER_LOCAL_WINDOW " << 
		max_number_peaks_per_local_window << endl;

	os << "#CONF NUMBER_OF_STRONG_PEAKS_PER_LOCAL_WINDOW " <<
		number_of_strong_peaks_per_local_window << endl;


	os << "#CONF LOCAL_WINDOW_SIZE " << setprecision(4) << local_window_size << endl;

	os << "#CONF TOLERANCE " << setprecision(4) << tolerance << endl;

	os << "#CONF PM_TOLERANCE " << setprecision(4) << pm_tolerance << endl;

	os << "#CONF MAX_EDGE_LENGTH " << max_edge_length << endl;

	os << "#CONF TERMINAL_SCORE " << terminal_score << endl;

	os << "#CONF DIGEST_SCORE " << digest_score << endl;

	os << "#CONF FORBIDDEN_PAIR_PENALTY " << forbidden_pair_penalty << endl;

	os << "#CONF NEED_TO_CORRECT_PM " << need_to_estimate_pm << " (0 - no, 1 - yes)" << endl;

	os << "#CONF SIZE_THRESHOLDS " << massThresholdsForSizes_.size() << " ";
	for (i=0; i<massThresholdsForSizes_.size(); i++)
	{
		int j;
		os << massThresholdsForSizes_[i].size() << " ";
		for (j=0; j<massThresholdsForSizes_[i].size(); j++)
			os << fixed << setprecision(4) << massThresholdsForSizes_[i][j] << " ";
	}
	os << endl;

	os << "#CONF REGION_THRESHOLDS " << region_thresholds.size();
	for (i=0; i<region_thresholds.size(); i++)
		os << " " << setprecision(4) << region_thresholds[i];

	os << endl;

	for (i=0; i<min_ranges.size(); i++)
	{
		os << "#CONF MASS_EXCLUDE_RANGE " << setprecision(3) << min_ranges[i] << " " << max_ranges[i] << endl;
	}
}


/************************************************************
outputs the selected aas for the different regions
*************************************************************/
void Config::print_table_aas(const ConversionTables& table, 
							 const vector<int>& aas) const
{
	int i;
	cout << aas.size() << " amino acids:" << endl;

	for (i=0; i<aas.size(); i++)
		cout << setw(6) << left << aas[i] << setprecision(4) << right << fixed << setw(8) 
			  << table.get_aa2mass(aas[i]) << "   " << left << 
			  table.get_aa2label(aas[i]) << endl;
}

void Config::print_session_aas() const
{
	cout << endl << "AMINO ACIDS" << endl;
	print_table_aas(session_tables,session_aas);
	cout << "N_TERM " << session_tables.get_aa2mass(N_TERM) << endl;
	cout << "C_TERM " << session_tables.get_aa2mass(C_TERM) << endl;
}


void Config::print_aa_variants() const
{
	int i;
	const vector<string>& aa2label = get_aa2label();

	for (i=0; i<aa_edge_combos.size(); i++)
	{
		cout << left << i << " " << aa_edge_combos[i].num_variants << " , " <<
			aa_edge_combos[i].variant_start_idx << "  ";
		int j;

		int pos = aa_edge_combos[i].variant_start_idx;
		cout << "(" << pos << ") ";
		for (j=0; j<aa_edge_combos[i].num_variants; j++)
		{
			int k;
			int num_aa = variant_vector[pos++];
			cout << " ";
			for (k=0; k<num_aa; k++)
				cout << aa2label[variant_vector[pos++]];
		}
		cout << endl;
	}
}



