#ifndef __MSCLUSTERAUXFUNS_H__
#define __MSCLUSTERAUXFUNS_H__

#include "../Common/includes.h"
#include "../PepNovo/PepNovo_includes.h"
#include "../PepNovo/AnnotationFile.h"

struct MsParameterStruct; // fwd dclr
class  AllScoreModels; // fwd dclr

inline
size_t computeMzIndex(mass_t m, mass_t mzIncrement, float tol)
{
    // NP3 deal with the mz limit problem
	size_t idx = (static_cast<size_t>(m/mzIncrement));
	if (m - mzIncrement * idx < tol)
		idx--;

	return idx;
}

inline
size_t computeDatRoundedMz(mass_t m, mass_t mzIncrement)
{
	return (static_cast<size_t>(computeMzIndex(m,mzIncrement, 0.0) * mzIncrement));
}


void splitDatListIntoSameMzGroups(const vector<string>& datList, 
								  vector< vector<string> >& datGroups);

void sortDatPathsAccordingToMz(vector<string>& datPaths);

/*! \fn convertSpectra
	\brief converts a file or list to another ouput.

	@param params The parameter struct for the run; it can use the following arguments:
	- inputFile or list
	- fileConversionType
	- outDir (optional)
	- peakDensity (optional)
	- sqsThreshold (optional)
	@param model a pointer the model class
*/
void convertSpectra(const MsParameterStruct* params, AllScoreModels* model);

void convertArchive(const MsParameterStruct* params, AllScoreModels* model);


void readIdsTitleFromIdFile(const MsParameterStruct* params, map<string,int>& idTitles);


#endif



