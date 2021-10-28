#ifndef __SPECTRALSIMILARITY_H__
#define __SPECTRALSIMILARITY_H__

#include <cmath>
#include "MsClusterIncludes.h"

const float MASS_TO_INT_RATIO = 1000000.0;

inline
int convertMassToInt(mass_t mass)
{
	return (static_cast<int>(MASS_TO_INT_RATIO * mass));
}

inline
mass_t convertIntToMass(int m)
{
	return (static_cast<float>(m)/MASS_TO_INT_RATIO);
}



/*! if SIMPLE_DISTANCE is defined the similarity distances between spectra are computed
    using only the mass (intensity of each peak = 1). this method works suprisngly well
    with the spectra that were tested. It also saves alot of memory space and processing time.*/
//#define SIMPLE_DISTANCE

/*!
	@struct DistancePeak
	\brief A peak used to compute distances between spectra.

	There are two modes for this structure:
	\arg If a simple distance is used (when \c SIMPLE_DISTANCE is defined), then only the mass is used
	  (each peak is assumed to have intensity 1)
    \arg If simple distance is not used, then there is a field called adjusted intensity which holds
	  the transformed intensity of the peak (this intensity only goes into the dot-product computation)
*/
struct DistancePeak {
	DistancePeak() : massAsInt(0), intensity(0.0), adjustedIntensity(0.0) {}
	bool operator< ( const DistancePeak& rhs) const
	{
		return (massAsInt<rhs.massAsInt);
	}

	int	  massAsInt;	/// The mass of the peaks (in Da converted to int with )
	float intensity;	/// The intensity of the peak (as observed in the spectrum)					
	float adjustedIntensity; /// The transformed intensity (used only for the dot-product computation)
};

/*!
	@struct DistancePeakList
	\breif Holds a set of peaks that represnt the spectrum when similarity distance is computed
*/
struct DistancePeakList {
	DistancePeakList() : numPeaks(0), sumSqrAdjustedIntensity(0.0) {}

	void print() const {
		cout << fixed << setprecision(3) << "Num distance: " << numPeaks << endl;
		for (int i=0; i<numPeaks; i++)
			cout << i << "\t" << convertIntToMass(peaks[i].massAsInt) << "\t" << peaks[i].intensity << "\t" << peaks[i].adjustedIntensity << endl;
		cout << endl;
	}
	int			 numPeaks;					/// number of distance peaks (is <= MAX_NUM_PEAKS_FOR_DISTANCE)
	float		 sumSqrAdjustedIntensity;	/// (for simple distance == numPeaks)
	DistancePeak peaks[MAX_NUM_PEAKS_FOR_DISTANCE]; /// holds the distance peaks
};


/********************************************************************************/
/*! @brief The function used to compute the similarity (distance) between two spectra

Performs a dot-product computation on the two DistancePeakLists. Since most distance
computations are between spectra of differnet spectra, it is faster to first ignore
intensities, and only if the similarity exceeds a certain threshold, do we compute
its accruate value.

@param listA A pointer to the DistancePeakList of spectrum A.
@param listB A pointer to the DistancePeakList of spectrum B.
@param tolerance The maximal mass distance between two peaks for which we still consider them as overlapping.
@return the similarity value between 0-1 (higher means more similar).*/
/*******************************************************************************/


inline 
float computeSimilarity(const DistancePeakList* listA,
				 const DistancePeakList* listB,
				 int   toleranceAsInt)
{
	assert( toleranceAsInt>0 );

	const DistancePeak* peaksA = listA->peaks;
	const DistancePeak* peaksB = listB->peaks;
	int idxA=0, idxB=0;

	float topSum = 0.0;
	while (idxA < listA->numPeaks && idxB < listB->numPeaks)
	{
		const int diff = peaksA[idxA].massAsInt - peaksB[idxB].massAsInt;
		if (diff<=0)
		{
			if (diff + toleranceAsInt >=0)
			{
				topSum += (peaksA[idxA++].adjustedIntensity*peaksB[idxB++].adjustedIntensity);
			}
			else
				++idxA;
		}
		else
		{
			if (diff <= toleranceAsInt)
			{
				topSum += (peaksA[idxA++].adjustedIntensity*peaksB[idxB++].adjustedIntensity);
			}
			else
				++idxB;
		}
	}
	assert(topSum>=0.0 && listA->sumSqrAdjustedIntensity>0.0 && listB->sumSqrAdjustedIntensity>0.0);
	const float similarity=(topSum / sqrt(listA->sumSqrAdjustedIntensity * listB->sumSqrAdjustedIntensity));

	return (similarity>=1.0 ? 1.0 : similarity);
	
}



inline
float computeSimilarityMixed(const DistancePeakList* listA,
						const DistancePeakList* listB,
						int   toleranceAsInt)
{
	assert( toleranceAsInt>0 );

	const DistancePeak* peaksA = listA->peaks;
	const DistancePeak* peaksB = listB->peaks;

	int matchCount = 0;
	int idxA=0, idxB=0;
	while (idxA<listA->numPeaks && idxB<listB->numPeaks)
	{
		const int diff = peaksA[idxA].massAsInt - peaksB[idxB].massAsInt;
		if (diff<=0)
		{
			if (diff + toleranceAsInt >=0)
			{
				++matchCount;
				++idxB;
			}
			++idxA;
		}
		else
		{
			if (diff<=toleranceAsInt)
			{
				++matchCount;
				++idxA;
			}
			++idxB;
		}
	}

	const float similarityEstimate = static_cast<float>(matchCount + matchCount) / static_cast<float>(listA->numPeaks + listB->numPeaks);
	if (similarityEstimate < MIN_SIMILARITY_FOR_ACCURATE_COMPUTATION)
		return similarityEstimate;

	// perform accurate computation that takes into account peak intensities
	idxA=0, idxB=0;
	float topSum = 0.0;
	while (idxA<listA->numPeaks && idxB<listB->numPeaks)
	{
		const int diff = peaksA[idxA].massAsInt - peaksB[idxB].massAsInt;
		if (diff<=0)
		{
			if (diff + toleranceAsInt >=0)
			{
				topSum += (peaksA[idxA++].adjustedIntensity*peaksB[idxB++].adjustedIntensity);
			}
			else
				++idxA;
		}
		else
		{
			if (diff <= toleranceAsInt)
			{
				topSum += (peaksA[idxA++].adjustedIntensity*peaksB[idxB++].adjustedIntensity);
			}
			else
				++idxB;
		}
	}

	assert(listA->sumSqrAdjustedIntensity>0.0 && listB->sumSqrAdjustedIntensity>0.0);
	const float similarity=(topSum / sqrt(listA->sumSqrAdjustedIntensity * listB->sumSqrAdjustedIntensity));
	return (similarity>1.0 ? 1.0 : similarity);
}



inline
float sharedPeakProportion(const DistancePeakList* listA,
						   const DistancePeakList* listB,
						   int   toleranceAsInt)
{
	assert( toleranceAsInt>0 );

	const DistancePeak* peaksA = listA->peaks;
	const DistancePeak* peaksB = listB->peaks;

	int matchCount = 0;
	int idxA=0, idxB=0;
	while (idxA<listA->numPeaks && idxB<listB->numPeaks)
	{
		const int diff = peaksA[idxA].massAsInt - peaksB[idxB].massAsInt;
		if (diff<=0)
		{
			if (diff + toleranceAsInt >=0)
			{
				++matchCount;
				++idxB;
			}
			++idxA;
		}
		else
		{
			if (diff<=toleranceAsInt)
			{
				++matchCount;
				++idxA;
			}
			++idxB;
		}
	}
	const float peakProportion = (static_cast<float>(matchCount + matchCount) / static_cast<float>(listA->numPeaks + listB->numPeaks));
	return (peakProportion > 1.0 ? 1.0 : peakProportion);
}


















#endif


