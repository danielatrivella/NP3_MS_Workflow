#ifndef __CLUSTER_H__
#define __CLUSTER_H__

/*! @file Cluster.h
	\brief Holds the class Cluster and a few supporting struct.
*/

#include "MsClusterIncludes.h"
#include "../PepNovo/AnnotationFile.h"
#include "../PepNovo/PeakList.h"
#include "../Common/VectorAllocator.h"

#include "SpectralSimilarity.h"



/*! if SIMPLE_DISTANCE is defined the similarity distances between spectra are computed
    using only the mass (intensity of each peak = 1). this method works suprisngly well
    with the spectra that were tested. It also saves alot of memory space and processing time.*/
//#define SIMPLE_DISTANCE

/*
	@struct DistancePeak
	\brief A peak used to compute distances between spectra.

	There are two modes for this structure:
	\arg If a simple distance is used (when \c SIMPLE_DISTANCE is defined), then only the mass is used
	  (each peak is assumed to have intensity 1)
    \arg If simple distance is not used, then there is a field called adjusted intensity which holds
	  the transformed intensity of the peak (this intensity only goes into the dot-product computation)
*
struct DistancePeak {
	bool operator< ( const DistancePeak& rhs) const
	{
		return (mass<rhs.mass);
	}

	float mass;	/// The mass of the peaks (in Da)
	float intensity; /// The intensity of the peak (as observed in the spectrum)	
						// this field is kept even with the simple distance because it is used
						// by the function that selects the consensus peask (selectDistancePeaksFromMultipleSigneltons)
						// if memeory becomes an issue, it can be removed from here and only
						// kept in that function
#ifndef SIMPLE_DISTANCE
	float adjustedIntensity; /// The transformed intensity (used only for the dot-product computation)
#endif	//SIMPLE_DISTANCE
};

/*!
	@struct DistancePeakList
	\breif Holds a set of peaks that represnt the spectrum when similarity distance is computed
*
struct DistancePeakList {
	DistancePeakList() : numPeaks(0), sumSqrAdjustedIntensity(0.0) {}

	void print() const {
		cout << fixed << setprecision(3) << "Num distance: " << numPeaks << endl;
		for (int i=0; i<numPeaks; i++)
		{
#ifdef SIMPLE_DISTANCE
			cout <<  i << "\t" << peaks[i].mass << endl;
#else
			cout << i << "\t" << peaks[i].mass << "\t" << peaks[i].intensity << "\t" << peaks[i].adjustedIntensity << endl;
#endif // SIMPLE_DISTANCE
		}
		cout << endl;
	}
	int			 numPeaks;					/// number of distance peaks (is <= MAX_NUM_PEAKS_FOR_DISTANCE)
	float		 sumSqrAdjustedIntensity;	/// (for simple distance == numPeaks)
	DistancePeak peaks[MAX_NUM_PEAKS_FOR_DISTANCE]; /// holds the distance peaks
};


/********************************************************************************
/*! The function used to compute the similarity (distance) between two spectra

Performs a dot-product computation on the two DistancePeakLists. If SIMPLE_DISTANCE
is used, this amounts to a simple normalized shared peak count. Otherwise this is
a dot-procut of two vectors of adjusted intensities

@param listA A pointer to the DistancePeakList of spectrum A.
@param listB A pointer to the DistancePeakList of spectrum B.
@param tolerance The maximal mass distance between two peaks for which we still consider them as overlapping.
@return the similarity value between 0-1 (higher means more similar).*
/*******************************************************************************
inline
float computeSimilarity(const DistancePeakList* listA,
						const DistancePeakList* listB,
						mass_t tolerance)
{
	assert( tolerance>0.0 );
	const DistancePeak* peaksA = listA->peaks;
	const DistancePeak* peaksB = listB->peaks;
	float topSum=0;
	int idxA=0, idxB=0;
	while (idxA<listA->numPeaks && idxB<listB->numPeaks)
	{
		if (fabs(peaksA[idxA].mass - peaksB[idxB].mass)<= tolerance)
		{
#ifdef SIMPLE_DISTANCE
			++topSum;
			++idxA;
			++idxB;
#else
			topSum+= (peaksA[idxA++].adjustedIntensity*peaksB[idxB++].adjustedIntensity);
#endif // SIMPLE_DISTANCE
			continue;
		}

		if (peaksA[idxA].mass < peaksB[idxB].mass)
		{
			idxA++;
		}
		else
			idxB++;
	}

	assert(listA->sumSqrAdjustedIntensity>0.0 && listB->sumSqrAdjustedIntensity>0.0);
	const float similarity=(topSum / sqrt(listA->sumSqrAdjustedIntensity * listB->sumSqrAdjustedIntensity));
	assert( similarity>=0.0 && similarity<=1.0);
	return similarity;
}

*/


/*!	@struct MassCount
	\brief A simple struct used when computing consensus spectra.
	
	This struct is used to remember for each peak mass how many spectra
	have peaks that can "cover" this mass (i.e., this mass falls in their [min,max] range.
*/
struct MassCount {
	MassCount() : mass(0.0), numSpectra(0) {}
	MassCount(mass_t m, clusterIdx_t n) : mass(m), numSpectra(n) {}
	const bool operator< (const MassCount& rhs) const
	{
		return (mass<rhs.mass);
	}
	mass_t        mass;
	clusterIdx_t  numSpectra;
};



/*! @class Cluster
	\brief This is a \c PeakList with additional supporitng information that is used in the clustering process.
*/
class Cluster : public PeakList {
	friend class MsClusterDataStorage;
public:

	Cluster() : indInPlay_(1), clusterTotalRT_(0.0), clusterTotalPrecursorIntensity_(0.0), clusterRT_(0.0),
				clusterRTMin_(0.0), clusterRTMax_(1000000.0), clusterTotalRTMin_(0.0), clusterTotalRTMax_(1000000.0),
				clusterCharge_(0), clusterMOverZ_(0.0), totalNonPrecursorIntensity_(0.0),
				clusterIdx_(MAX_CLUSTER_IDX),  assignedClusterIdx_(MAX_CLUSTER_IDX),
				clusterSize_(1), singletonDistancePeaks_(0), singletonIdxVector_(0), lowestSimilarityPos_(0),
				minSimilarityToJoin_(2.0), numSimilarityPairs_(0)  {}

	~Cluster();

	bool operator< (const Cluster& rhs) const
	{
		return (clusterMOverZ_<rhs.clusterMOverZ_);
	}

	static mass_t getPeakIndexTolerance() { return peakIndexTolerance_; }
	static mass_t getIsoTolerance()		  { return isoTolerance_; }
	static mass_t getPeakTolerance()	  { return peakTolerance_; }
	static int	  getPeakIndexToleranceAsInt() { return peakIndexToleranceAsInt_; }


	bool		 getIndInPlay()			 const { return indInPlay_; }
	void		 setIndInPlay(bool b)		   { indInPlay_ = b; }
	clusterIdx_t getClusterIdx()         const { return clusterIdx_; }
	void		 setClusterIdx(clusterIdx_t idx) { clusterIdx_ = idx; }
	clusterIdx_t getAssignedClusterIdx() const { return assignedClusterIdx_; }
	unsigned int getClusterSize()        const { return clusterSize_; }
	mass_t		 getClusterMOverZ()		 const { return clusterMOverZ_; }
	short		 getClusterCharge()		 const { return clusterCharge_; }
	void		 setClusterCharge(int c)	   { clusterCharge_ = static_cast<short>(c); }
	const AllocatedVector<clusterIdx_t>* const getSingletonVector() const { return singletonIdxVector_; }
	const DistancePeakList* getDistancePeaks() const { return &distancePeaks_; }
	const DistancePeakList* getSingletonPeaks() const { return (singletonDistancePeaks_ ? singletonDistancePeaks_ : &distancePeaks_); }
	const unsigned int* const getTopPeakIdxs() const { return topPeakIdxs_; }
	const float* const			getBestSimilarityValues()	   const { return bestSimilarityValues_; }
	const clusterIdx_t* const	getBestSimilarityClusterIdxs() const { return bestSimilarityClusterIdxs_; }
	float		getMinSimilarityToJoin() const { return minSimilarityToJoin_; }
	int			getNumSimilarityPairs()  const { return numSimilarityPairs_; }

	// NP3 GOT RT getter
	float getClusterRT() const { return clusterRT_; }
	static float getRtTolerance() { return rtTolerance_; }
	// NP3 GOT RT width getters
	float getClusterRTMin() const { return clusterRTMin_; }
	float getClusterRTMax() const { return clusterRTMax_; }
	intensity_t getClusterTotalPrecursorIntensity() const {return clusterTotalPrecursorIntensity_;}
	// NP3 scale factor
	static float getScaleFactor() { return scale_factor_; };


	/*! \fn getNumSingletonsIncludingSelf
	\brief Computes the number of actual \c PeakList that are contained in the cluster.

	This number should not be confused with the cluster size, since when clustering multiple datasets,
	a single spectrum might represent several singletons, however we do not have access to their peak lists any more.
	*/
	clusterIdx_t getNumSingletonsIncludingSelf() const {
		if (! singletonIdxVector_)
			return 1;
		return (singletonIdxVector_->getSize() + 1);
	}


	/*! \fn updateSimilarityValueInPosition
	Updates the similarity value in a specific position.
	@param i The position in the chached vectors.
	@param similarity The similarity value between this cluster and the other cluster.
	@param clusterIdx The index of the other cluster.
	*/
	void updateSimilarityValueInPosition(size_t i, float similarity, clusterIdx_t clusterIdx)
	{
		bestSimilarityValues_[i]=similarity;
		bestSimilarityClusterIdxs_[i]=clusterIdx;
	}

	/*!
	Updates the cached similarity value with another cluster. It only updates if
	the similarity is higher than one of the other cached pairs (which is then replaced by
	this new cluster
	@param similarity The similarity value between this cluster and the other cluster.
	@param clusterIdx The index of the other cluster.
	*/
	void updateSimilarityValues(float similarity, clusterIdx_t clusterIdx)
	{
		if (similarity > bestSimilarityValues_[lowestSimilarityPos_])
		{
			bestSimilarityValues_[lowestSimilarityPos_] = similarity;
			bestSimilarityClusterIdxs_[lowestSimilarityPos_] = clusterIdx;
			lowestSimilarityPos_ = 0;
			for (int i=1; i<NUM_TOP_SIMILARITIES_TO_SAVE; i++)
				if (bestSimilarityValues_[i]<bestSimilarityValues_[lowestSimilarityPos_])
					lowestSimilarityPos_ = i;
		}
	}

	/*!
	Updates the similarity value for clusterIdx. Checks if the clsuterIdx is already
	there, if it is, the function replaces the value; if not, the function updates the position
	with the lowest similarity.
	@param similarity The similarity value between this cluster and the other cluster.
	@param clusterIdx The index of the other cluster.
	*/
	void updateSimilarityValuesCheckIfClusterThere(float similarity, clusterIdx_t clusterIdx)
	{
		int i;
		for (i=0; i<NUM_TOP_SIMILARITIES_TO_SAVE; i++)
			if (bestSimilarityClusterIdxs_[i]==clusterIdx)
			{
				bestSimilarityValues_[i]=0.0;
				lowestSimilarityPos_=i;
				break;
			}

		if (similarity > bestSimilarityValues_[lowestSimilarityPos_])
		{
			bestSimilarityValues_[lowestSimilarityPos_] = similarity;
			bestSimilarityClusterIdxs_[lowestSimilarityPos_] = clusterIdx;
			lowestSimilarityPos_ = 0;
			for (int i=1; i<NUM_TOP_SIMILARITIES_TO_SAVE; i++)
				if (bestSimilarityValues_[i]<bestSimilarityValues_[lowestSimilarityPos_])
					lowestSimilarityPos_ = i;
		}
	}

	/*!
	Initialzes a new single spectrum cluster (singleton) with the peaks and other information.

	@param idx The unique cluster index (this number is unqiue for this process). It must be
			between 0 and MAX_CLUSTER_IDX-1.
	@param header A pointer the the header of the peak list for this cluster.
	@param peaks A pointer to an array of \c Peak instances (the peak list of the singleton)
	@param numPeaks Number of peaks in \c peaks
	*/
	bool createNewCluster(clusterIdx_t idx,
						  const SingleSpectrumHeader* header,
						  const Peak* peaks,
						  int   numPeaks);

	// copies the distance peaks to the singleton distance peaks (so this information is not lost when p-values are computed)
	void backupDistancePeaks();

	static void setTolerances(mass_t fragmentTolerance, float percursorPPMs=MAX_FLOAT);
	// NP3 GOT RT tolerance
	static void setRtTolerance(float rtTol);
	// NP3 scale factor
	static void setScaleFactor(float factor);

	/*!
	Copies a cluster but maintains an empty singletonIdxVector_ which needs to be created manually.
	*/
	void copyWithoutSingletonIdxVector(const Cluster& rhs);


	void detailedPrint() const;

	bool checkIfTopPeaksHaveMatch(const Cluster& other) const
	{
		for (size_t i=0; i<NUM_PEAKS_FOR_HEURISTIC; i++)
			for (size_t j=0; j<NUM_PEAKS_FOR_HEURISTIC; j++)
			{
				const unsigned int diff = (topPeakIdxs_[i] >= other.topPeakIdxs_[j] ? topPeakIdxs_[i] - other.topPeakIdxs_[j] :  other.topPeakIdxs_[j] - topPeakIdxs_[i]);
				if (diff<=2)
					return true;
			}
		return false;
	}

protected:

	static mass_t peakIndexTolerance_; // The tolerance for computing peak indexes (should be isoTolerance/2)
	static mass_t peakTolerance_;	/// The tolerance in Da used for fragment peaks
	static mass_t isoTolerance_;	/// The tolerance in Da used for isotopic peaks (usually tighter than \c peakTolerance_ ).
	static mass_t precursorPPM_;	/// The tolerance in Da used for precusor mass (used with high-res data)
	static int	  peakIndexToleranceAsInt_;

	bool		 indInPlay_; /// states if this cluster can be joined ( equivalent to asking clusterIdx_ == assignedClusterIdx_)

	short		 clusterCharge_; /// The charge assigned to the cluster (0 if unknown)
	mass_t       clusterMOverZ_; /// The precursor m/z of the cluster (weighted average of the cluster members' precursor m/z)
	intensity_t  totalNonPrecursorIntensity_; /// Sum of the peaks' intensities

	// NP3 GOT RT mean and tolerance
	float clusterTotalRT_; // stores the weighted sum of the cluster members
	intensity_t clusterTotalPrecursorIntensity_; // stores the sum of the clsuters members precursor intensity
	float clusterRT_; // stores the cluster current rt weighted mean
	static float rtTolerance_; // the tolerance in seconds for clusters retention time - clusters with rt abs diff greater than this wont be joined
	// NP3 Got RT with in cluster
	float clusterTotalRTMin_;
	float clusterTotalRTMax_;
	float clusterRTMin_;
	float clusterRTMax_;
	// NP3 scale factor
	static float scale_factor_;

	clusterIdx_t clusterIdx_;		   /// The unique index of this cluster (must be between 0 and \c MAX_CLUSTER_IDX-1)
	clusterIdx_t assignedClusterIdx_; /// The index to which this cluster was added (initally the same as \c clusterIdx_)

	unsigned int clusterSize_;      /*! Number of spectra in this cluster (this does not equal number of singletons in the cluster!) */


	DistancePeakList* singletonDistancePeaks_; /// these are allocated only if the distance peaks get overwritten because this is the main cluster
	DistancePeakList  distancePeaks_; /// The peaks that represent the cluster in the similarity computations


	AllocatedVector<clusterIdx_t>* singletonIdxVector_; // the indexs of the spectra participating in the cluster

	unsigned int topPeakIdxs_[NUM_PEAKS_FOR_HEURISTIC]; /*! These are the bin indexes for the top \c NUM_PEAKS_FOR_HEURISTIC peaks
													       (the indexes designate lists in \c clusterIdxsLists_ )
													        The indexes are computed using the function \c computeMzIndex().*/

	int			  lowestSimilarityPos_;	/*! Points to which of the \c NUM_TOP_SIMILARITIES_TO_SAVE entries has the lowest similarity
										   this refers to a position in the arrays \c bestSimilarityValues_ and \c bestSimilarityClusterIdxs_ .*/

	float		  minSimilarityToJoin_; /*! States the minimal simialrity needed to join this cluster to another. This number
											is computed according to the number of distance peaks and the number of pairs being compared. */

	int			  numSimilarityPairs_; /*! The number of pairs of spectra that were considered when computing the similarity (affects the threshold and p-value). */

	float         bestSimilarityValues_[NUM_TOP_SIMILARITIES_TO_SAVE];  /*! The similarity values with the clusters listed in \c bestSimilarityClusterIdxs_ */

	clusterIdx_t  bestSimilarityClusterIdxs_[NUM_TOP_SIMILARITIES_TO_SAVE]; /*! The indexes of the clusters with which this cluster has maximal similarites */


	bool checkDistancePeaksOk() const;
	bool selectDistancePeaksAndTopIdxs();
	void setAdjustedIntensities();
	void adjustMaxPossiblePeakCounts(const vector<MassCount>& minMassCounts,
									 const vector<MassCount>& maxMassCounts);
};


/*! \fn readAnnotatedSpectraIntoClusters
	Reads spectra annoated spectra files into clusters.

	@param list List of paths (or a single spectra path).
	@param config The models config.
	@param sa The SpectraAggregator class which will hold all the read headers.
	@param clusters The clusters which will hold the spectra.
	@param annotationMap a mapping of all peptide annotations to the different clusters which hold the peptide/charge.
*/
void readAnnotatedSpectraIntoClusters(const string& list,
									  const Config* config,
									  SpectraAggregator& sa,
									  vector<Cluster>& clusters,
									  map<Annotation, vector<size_t> >& annotationMap);


#endif
