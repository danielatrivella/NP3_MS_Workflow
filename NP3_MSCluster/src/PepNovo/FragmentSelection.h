#ifndef __FRAGMENTSELECTION_H__
#define __FRAGMENTSELECTION_H__

#include "SpectraAggregator.h"
#include "AnnotatedSpectrum.h"
#include "Config.h"


void createFragmentsAccordingToOffsetCounts(const SpectraAggregator& sa, 
										    const Config* const config,
										    FragmentTypeSet& fts, 
										    float mininmalFragmentProbability = 0.05);







// returns the probablility of observing a the different types of
// fragments for different charges/sizes/regions
void collectProbabilitiesOfFragments(const SpectraAggregator& sa, 
									 Config *config,
									 vector< vector< vector< vector<double> > > >& frag_probs, // charge , size, region, frag_idx
									 vector< vector< vector< vector<double> > > >& in_range_counts,
									 double &avg_rand, 
									 int &num_files_used, 
									 int charge = 0,
									 bool verbose = false);




// Selects for a charge state's regional models the set of fragments
// that should be used.
void selectRegionalFragmentSets(const SpectraAggregator& sa, Config *config, 
								int charge, bool verbose = true);


void selectFragmentCombinations(const SpectraAggregator& sa, Config* config,
								int charge, int maxNumCombos = 3);


/*********************************************************************
Adds the counts for peaks around a breakage to their respective bins

**********************************************************************/
void add_offset_counts_arround_mass(vector<int>& counts, Spectrum *spec,
									mass_t min_mass, mass_t max_mass, 
									mass_t bin_coef, mass_t break_mass,
									int charge);


void select_fragments_from_bins(vector<int>& counts, FragmentTypeSet& fts, int max_num_frags, 
								int charge, int orientation, mass_t min_offset_mass, 
								mass_t bin_coef, mass_t tolerance);


void calculateTrueFragmentProbabilities(const SpectraAggregator& sa,
							 FragmentTypeSet& fts,  float minProbability);


#endif

