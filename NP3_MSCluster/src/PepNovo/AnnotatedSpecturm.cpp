#include "AnnotatedSpectrum.h"




/***************************************************************
Finds all the peaks in the given spectrum that support a breakage
at a certain location.
****************************************************************/
void annotate_breakage(Spectrum *spec,
					   mass_t pm_with_19, 
					   int peptide_size_idx,
					   Breakage& breakage)
{
	if (breakage.region_idx <0 || breakage.mass <0)
	{
		cout << "Must first set breakage mass and region!" << endl;
		exit(1);
	}

	const Config *config = spec->getConfig();
	const vector<FragmentType>& all_fragments = config->get_all_fragments();

	const RegionalFragments & rf = config->get_regional_fragments(spec->getCharge(),
													peptide_size_idx, breakage.region_idx);

	const vector<int>& frag_type_idxs = config->get_regional_fragment_type_idxs(spec->getCharge(),
													peptide_size_idx, breakage.region_idx);

	const mass_t min_peak_mass = spec->get_min_peak_mass();
	const mass_t max_peak_mass = spec->get_max_peak_mass();
	const mass_t tolerance = config->getTolerance();
	const mass_t small_tolerance = tolerance * 0.5;

	const int spec_charge = spec->getCharge();


	// if the fragments vector has a size>0, we can assume
	// that it contains previous fragments

	breakage.frag_type_idxs_not_visible.clear();

	bool has_previous_frags = false;
	if (breakage.fragments.size() >0)
		has_previous_frags=true;

	breakage.parent_charge = spec->getCharge();
	breakage.parent_size_idx = peptide_size_idx;

	int f;
	for (f=0; f<frag_type_idxs.size(); f++)
	{
		const int frag_type_idx = frag_type_idxs[f];
		const FragmentType& ft = all_fragments[frag_type_idx];
		
		if (ft.charge> spec_charge)
			continue;

		// if somehow already annotated, skip
		if (has_previous_frags && breakage.get_position_of_frag_idx(frag_type_idxs[f])>= 0)
			continue;
			
		// try use parent frag for expected position of this fragment
		mass_t exp_frag_mass;

		const int& parent_frag_pos = breakage.get_position_of_frag_idx(ft.parent_frag_idx);
		bool used_parent_frag_for_exp_mass = false;
		if (parent_frag_pos>=0 && 
			breakage.fragments[parent_frag_pos].mass > 0)
		{
			exp_frag_mass = breakage.fragments[parent_frag_pos].mass +
							ft.offset_from_parent_frag;

			used_parent_frag_for_exp_mass= true;
		}
		else
			exp_frag_mass = ft.calc_expected_mass(breakage.mass,pm_with_19);


		// check if this frag is visible, if not, add it to no-visible list
		if (exp_frag_mass<min_peak_mass || exp_frag_mass>max_peak_mass)
		{
			breakage.frag_type_idxs_not_visible.push_back(frag_type_idx);
		//	cout << "NOZ VIZ: " << frag_type_idx << "  " << exp_frag_mass << "(" <<
		//		min_peak_mass << " , " << max_peak_mass << ")" << endl;
			continue;
		}

		// look for top peak. if we already found the parent fragment,
		// then we can use a smaller tolerance, since this peak should 
		// align closely with the more common parent peak
		int peak_idx = spec->findPeakWithMaxIntensity(exp_frag_mass, 
			(used_parent_frag_for_exp_mass ? small_tolerance : tolerance ) );
			
		if (peak_idx<0)
			continue;

		// check if this peak was already used for a fragment with a similar charge 
		// and orientaition (for instance -H2O and -NH3)
		bool add_frag = true;
		if (ft.offset>0)
		{
			int j;
			for (j=0; j<breakage.fragments.size(); j++)
			{
				const int j_frag_type_idx = breakage.fragments[j].frag_type_idx;
				if (breakage.fragments[j].peak_idx == peak_idx &&
					all_fragments[j_frag_type_idx].charge == ft.charge &&
					all_fragments[j_frag_type_idx].orientation == ft.orientation)
				{
					mass_t dis1 = fabs(breakage.fragments[j].mass - breakage.fragments[j].expected_mass);
					mass_t dis2 = fabs(breakage.fragments[j].mass - exp_frag_mass);
					if (dis1<dis2)
					{
						add_frag = false;
						break;
					}
					else // remove the frag annotation for frag j because the peak 
						 // is better for the frag f
					{
						breakage.remove_fragment(j_frag_type_idx);
						break;
					}
				}
			}
		}

		if (! add_frag)
			continue;
	
		BreakageFragment brf;
		brf.peak_idx=peak_idx;
		brf.expected_mass=exp_frag_mass;
		brf.frag_type_idx = frag_type_idx;
		brf.intensity = spec->getPeakIntensity(peak_idx);
		brf.mass = spec->getPeakMass(peak_idx);
		brf.is_strong_fragment = rf.is_a_strong_frag_type(frag_type_idx);

		breakage.add_fragment(brf);
	}			
}


void AnnotatedSpectrum::annotate_spectrum(mass_t pm_with_19, int spec_charge, bool reset_annotations)
{
	int i;

	if (! config_ || numPeaks_<1)
	{
		cout << "Error: spectrum must be read before annotation!" << endl;
		exit(1);
	}

	if (breakages.size() == 0 || reset_annotations)
	{
		breakages.clear();
		peak_annotations.clear();
		peak_annotations.resize(numPeaks_);
	}

	if (pm_with_19 < 0)
	{
		cout << "Error: negative mass given to annotate_spectrum!" << endl;
		exit(1);
	}

	int org_charge= charge_;
	if (spec_charge>0)
		charge_=spec_charge;

	if (this->peptide_.get_length() == 0)
		return;

	vector<mass_t> exp_breakage_masses;
	peptide_.calc_expected_breakage_masses(config_,exp_breakage_masses);
	
	if (breakages.size() != exp_breakage_masses.size())
		breakages.resize(exp_breakage_masses.size());

	// loop on all breakages and annotate them
	for (i=0; i<exp_breakage_masses.size(); i++)
	{
		const int region_idx = config_->calc_region_idx(exp_breakage_masses[i],
							pm_with_19, charge_, minimalPeakMass_, maximalPeakMass_);

		breakages[i].mass = exp_breakage_masses[i];
		breakages[i].region_idx = region_idx;	
		annotate_breakage((Spectrum *)this,pm_with_19, sizeIndex_, breakages[i]);
	}

	// assign annotations to peaks
	for (i=0; i<breakages.size(); i++)
	{
		int f;
		for (f=0 ; f<breakages[i].fragments.size(); f++)
		{
			if (breakages[i].fragments[f].mass>0)
			{
				PeakAnnotation p;
				p.frag_type_idx = breakages[i].fragments[f].frag_type_idx;
				p.breakage_idx = i;

				const FragmentType& ft = config_->get_fragment(p.frag_type_idx);
				
				int idx = (ft.orientation == PREFIX ) ? i : breakages.size() - i -1;
				ostringstream os;
				os << idx << ft.label;
				p.label = os.str();
				peak_annotations[breakages[i].fragments[f].peak_idx].push_back(p);
			}
		}
	}
	charge_ = org_charge;
}


// chooses the charge that gives a good_pm_with_19
void AnnotatedSpectrum::set_charge_from_seq_and_m_over_z()
{
	mass_t seq_mass = get_true_mass_with_19();
	int c;
	for (c=1; c<10; c++)
	{
		mass_t p19 = mOverZ_ * c - (c-1)*MASS_PROTON;
		if (fabs(p19-seq_mass)<10.0)
		{
			charge_ = c;
			break;
		}
	}

	if (c == 10)
	{
		cout << "Error: couldn't find mataching charge for sequence!" << endl;
		cout << "Seq : " << this->peptide_.as_string(config_) << " " << seq_mass << endl;
		cout << "m/z : m_over_z" << endl;
		exit(1);
	}
}


// how many of the expected fragment peaks were observed
void AnnotatedSpectrum::get_number_observed_frags(const vector<int>& frag_types, 
												  int& num_obs, int &num_exp) const
{
	vector<mass_t> break_masses;

	peptide_.calc_expected_breakage_masses(config_,break_masses);
	mass_t true_mass = peptide_.get_mass() + MASS_OHHH;

	num_obs=0;
	num_exp=0;
	int f;
	for (f=0; f<frag_types.size(); f++)
	{
		const FragmentType& ft = config_->get_fragment(frag_types[f]);

		int b;
		for (b=1; b<break_masses.size(); b++)
		{
			mass_t exp_mass = ft.calc_expected_mass(break_masses[b],true_mass);
			if (exp_mass > minimalPeakMass_ && exp_mass<maximalPeakMass_)
			{	
				num_exp++;

				if (findPeakWithMaxIntensity(exp_mass,config_->getTolerance())>=0)
					num_obs++;
			}
		}
	}
}

int AnnotatedSpectrum::get_num_observed_frags(int frag_idx) const
{
	vector<mass_t> break_masses;

	peptide_.calc_expected_breakage_masses(config_,break_masses);
	mass_t true_mass = peptide_.get_mass() + MASS_OHHH;

	int num_obs=0;
	const FragmentType& ft = config_->get_fragment(frag_idx);

	int b;
	for (b=0; b<break_masses.size(); b++)
	{
		mass_t exp_mass = ft.calc_expected_mass(break_masses[b],true_mass);
	
		if (findPeakWithMaxIntensity(exp_mass,config_->getTolerance())>=0)
			num_obs++;
	
	}
	return num_obs;
}

int AnnotatedSpectrum::get_num_annotated_peaks() const
{
	int i,n=0;

	for (i=0; i<this->peak_annotations.size(); i++)
		if (peak_annotations[i].size()>0)
			n++;

	return n;
}


float AnnotatedSpectrum::get_explianed_intensity() const
{
	if (peak_annotations.size()==0 || (numPeaks_ != peak_annotations.size()) )
	{
		cout << "Error: mismatch in peaks annotations!" << endl;
		exit(1);
	}

	float tot_inten=0;
	float ann_inten=0;
	int i;
	for (i=0; i<numPeaks_; i++)
	{
		tot_inten += peaks_[i].intensity;
		if (peak_annotations[i].size()>0)
			ann_inten += peaks_[i].intensity;
	}
	
	return (ann_inten/tot_inten);
}

/**************************************************************************
Checks if the spectrum has a suffcient stretch of b or y ladders
***************************************************************************/
bool AnnotatedSpectrum::has_stretch_of_b_or_y(int min_stretch_length, int max_skip)
{
	vector<mass_t> b_ladder,y_ladder;
	int b_frag_idx = config_->get_frag_idx_from_label("b");
	int y_frag_idx = config_->get_frag_idx_from_label("y");
	const vector<mass_t>& aa2mass = config_->get_aa2mass();
	mass_t tolerance = config_->getTolerance();

	if (breakages.size() ==0)
	{
		cout << "Error: need to first annotate spectrum!" << endl;
		exit(1);
	}

	b_ladder.resize(breakages.size(),0);
	b_ladder[0]=1;
	b_ladder[b_ladder.size()-1]=get_true_mass()+1;
	y_ladder.resize(breakages.size(),0);
	y_ladder[0]=get_true_mass_with_19();
	y_ladder[y_ladder.size()-1]=19;

	int i;
	vector<mass_t> exp_b_masses,exp_y_masses;
	peptide_.calc_expected_breakage_masses(config_,exp_b_masses);
	exp_y_masses = exp_b_masses;

	for (i=0; i<exp_y_masses.size(); i++)
		exp_y_masses[i] = get_true_mass_with_19() - exp_y_masses[i];
	for (i=0; i<exp_b_masses.size(); i++)
		exp_b_masses[i]+=0.95;


	for (i=0; i<breakages.size(); i++)
	{
		int b_pos= breakages[i].get_position_of_frag_idx(b_frag_idx);
		if (b_pos>=0)
			b_ladder[i]=breakages[i].fragments[b_pos].mass;

		int y_pos = breakages[i].get_position_of_frag_idx(y_frag_idx);
		if (y_pos>=0)
			y_ladder[i]=breakages[i].fragments[y_pos].mass;
	}

	
	const vector<int>& amino_acids = peptide_.get_amino_acids();

//	for (i=0; i<b_ladder.size(); i++)
//		cout << y_ladder[i] << " " << exp_y_masses[i] << endl;


	// check for stretch
	for (i=0; i<= b_ladder.size()-min_stretch_length; i++)
	{
		int gaps=0;
		int j;

		mass_t last_b = exp_b_masses[i];
		mass_t exp_mass =0;
		for (j=0; j<min_stretch_length; j++)
		{
			if (j>0)
				exp_mass += aa2mass[amino_acids[i+j-1]];

			if (b_ladder[i+j]>0 && (j==0 || 
									fabs(b_ladder[i+j]-last_b - exp_mass)<=0.5 * tolerance ) )
			{
				gaps=0;
				exp_mass =0;
				last_b = exp_b_masses[i+j];
			}
			else
				gaps++;

			if (gaps>max_skip)
				break;
		}

		if (j==min_stretch_length)
			return true;
	} 

	// check for stretch
	int num_aa = amino_acids.size();
	for (i=0; i<= y_ladder.size()-min_stretch_length; i++)
	{
		int gaps=0;
		int j;
		mass_t last_y = exp_y_masses[i];
		mass_t exp_mass =0;
		for (j=0; j<min_stretch_length; j++)
		{
			if (j>0)
				exp_mass += aa2mass[amino_acids[i+j-1]];

			if (y_ladder[i+j]>0 && (j==0 || 
									fabs(last_y - y_ladder[i+j]- exp_mass)<=0.5*tolerance ) )
			{
				gaps=0;
				exp_mass =0;
				last_y = exp_y_masses[i+j];
			}
			else
				gaps++;

			if (gaps>max_skip)
				break;
		}
		if (j==min_stretch_length)
			return true;
	}

//	exit(0);
	return false;
}


/********************************************************************************
Extracts tables of peak intensities and peak masses.
#rows = all fragments
# column = # num breakages (peptide length + 1)
*********************************************************************************/
void AnnotatedSpectrum::extract_annotated_intens_and_masses(
											 vector< vector<intensity_t> >& intens,
											 vector< vector<mass_t> >&      masses) const
{
	const int num_frags = config_->get_all_fragments().size();
	const int num_breakages = breakages.size();
	int i;
	intens.resize(num_frags);
	masses.resize(num_frags);
	for (i=0; i<num_frags; i++)
	{
		intens[i].resize(num_breakages,0);
		masses[i].resize(num_breakages,0);
	}

	for (i=0; i<peak_annotations.size(); i++)
	{
		int j;
		for (j=0; j<peak_annotations[i].size(); j++)
		{
			const PeakAnnotation& ann = peak_annotations[i][j];
			if (ann.breakage_idx>0 && ann.frag_type_idx<num_frags)
			{
				intens[ann.frag_type_idx][ann.breakage_idx]=peaks_[i].intensity;
				masses[ann.frag_type_idx][ann.breakage_idx]=peaks_[i].mass;
			}
		}
	}
}







void AnnotatedSpectrum::print_annotated(ostream& os) const
{
		int i;
	if (title_.length() > 0)
		os << "#TITLE " << title_ << endl;
	if (peptide_.get_num_aas()>0)
		os << "#SEQ " << peptide_.as_string(config_) << endl;

	os << fixed << setprecision(2) << originalPmWith19_ << " " << charge_ << endl;

	cout << "Peak  Mass        Inten   Rank   Annotation" << endl;
	int strong_idx=0;
	for (i=0; i<numPeaks_; i++)
	{
		os << left << setw(5) << i << setw(8) << setprecision(NUM_SIG_DIGITS) 
			<< fixed << right  << peaks_[i].mass;
		os << setw(12) << right  << setprecision(2) << peaks_[i].intensity;
		if (ranks_.size()>i)
			os << setw(4) << right << ranks_[i];

		if (strong_idx<strongPeakIndexes_.size())
		{
			if (strongPeakIndexes_[strong_idx]==i)
			{
				strong_idx++;
				os << " *";
			}
			else
				os << "  ";
		}
		else
			os << "  ";

		int j;
		for (j=0; j<peak_annotations[i].size(); j++)
			os << "  " << peak_annotations[i][j].label;
			
		os << endl;
	}
}


void print_dataset_spectra_by_stats(Config *config, char *mgf_file)
{
/*	FileManager fm;
	FileSet     fs;

	fm.init_from_mgf(config,mgf_file);
	fs.select_all_files(fm);

	int b_frag_idx = config->get_frag_idx_from_label("b");
	int y_frag_idx = config->get_frag_idx_from_label("y");

	int counter=0;
	while (1)
	{
		AnnotatedSpectrum as;

		if (! fs.get_next_spectrum(fm,config,&as))
			break;

		as.annotate_spectrum(as.get_true_mass_with_19());

		if (counter==0)
			as.print_expected_by();

		int num_b = as.get_num_observed_frags(b_frag_idx);
		int num_y = as.get_num_observed_frags(y_frag_idx);
		float exp_int = as.get_explianed_intensity();

		
		cout << counter++ << " " << num_b << " " << num_y << " " << exp_int << endl;
	}*/
}




