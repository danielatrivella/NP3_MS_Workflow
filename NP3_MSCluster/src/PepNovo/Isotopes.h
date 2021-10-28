#ifndef __ISOTOPES_H__
#define __ISOTOPES_H__

#include "PepNovo_includes.h"


// returns the ratios of the expected isotopic peaks
// bases numbers on the averagine statistics
void calc_expected_iso_ratios(mass_t peak_mass, vector<float>& ratios,
						int num_ratios=5);


#endif


