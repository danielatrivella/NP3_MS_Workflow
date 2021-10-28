#ifndef __MSCLUSTERDATASTORAGE_H__
#define __MSCLUSTERDATASTORAGE_H__

#include "Cluster.h"
#include "DatFileManager.h"
#include "ClusterOutputter.h"
#include "MsArchive.h"
#include "ClustFile.h"
#include "PeakWeightTable.h"

#include "../PepNovo/Config.h"

struct MsParameterStruct; // fwd dclr
class  MsSimilarityModel;

struct ClusterPositionPair {
	ClusterPositionPair() : lowerPos(0), higherPos(0) {}
	ClusterPositionPair(clusterIdx_t low, clusterIdx_t high) : lowerPos(low), higherPos(high) {}
	bool operator< (const ClusterPositionPair& rhs) const
	{
		return (lowerPos<rhs.lowerPos ||
				lowerPos == rhs.lowerPos && higherPos<rhs.higherPos);
	}

	bool operator != (const ClusterPositionPair& rhs) const
	{
		return (lowerPos != rhs.lowerPos || higherPos != rhs.higherPos);
	}

	bool operator == (const ClusterPositionPair& rhs) const
	{
		return (lowerPos == rhs.lowerPos && higherPos == rhs.higherPos);
	}

	clusterIdx_t lowerPos;
	clusterIdx_t higherPos;
};

class MsClusterDataStorage {
	friend class MsClusterAlgorithm;
public:
	MsClusterDataStorage() : config_(0), params_(0), simModel_(0), outDir_(std::string()), 
		nameStr_(std::string()), datasetIdx_(0), batch_(0), nextClusterPos_(0),  
		nextPeakPos_(0), runningClusterIdx_(0), runningOutputIdx_(1), totalSpectraWritten_(0),
		totalClustersWritten_(0), verboseLevel_(0), minSqsForSingleton_(0.0) {}
	// NP3 GOT change runningOutputIdx to 1
	~MsClusterDataStorage();

	size_t getNextClusterPos()      const { return nextClusterPos_; }

	size_t getNextPeakPos() const { return nextPeakPos_; }

	void allocateMemory( float availableGB, bool verbose = true);

	/*! \brief Adds new batches of spectra for clustering.

		First removes some spectra to make room for the new ones (tries to keep all spectra that
		are in the window that can still be compared). After that adds in new spectra in the next
		empty locations (added spectra are also entered as clusters).

		@param Window the m/z window size (how much margin needs to be kept)
		@param numSpectraFreed Holds the number of spectra that were actually freed (a return value)

		\return The number of spectra added. 
	*/
	clusterIdx_t addNewSpectra( mass_t window, clusterIdx_t* numSpectraFreed = 0);

	/*! \brief Initializes the various data structures invovled in storing and writing the spectra data.

	The function examines all the input dat files and takes not of statistics such as total
	numebr of spectra/peaks needed and the min/mz precursor m/z.
	*/
	void initialize(MsParameterStruct* params, const Config* config, const MsSimilarityModel* simModel);

	clusterIdx_t getTotalSpectraToCluster() const { return datFileManager_.getTotalSpectraToCluster(); }

	size_t findPosOfClusterWithMzAbove(mass_t mz, size_t startPos=0) const;

	// makes a list of all pairs of psoitions that need to have their similarity computed
	// (these pairs should have a cluster m/z that is at most *window* apart from each other)
	// the indexes are taken from the lists in newClusterIdxsLists_, and it is assumed that
	// every one of the higherPosition clusters is in play
	// each pair has one element from the range [firstPosLower,  lastPosLower)
	// and the second element form the range    [firstPosHigher, lastPosHigher)
	// the pairs are stored and sorted in the vector similarityPairPositions_
	size_t makeListOfSimilarities(size_t firstPosLower,  size_t lastPosLower,
								  size_t firstPosHigher, size_t lastPosHigher,
								  mass_t window);

	// Computes for each cluster the minimal similarity needed to join the cluster to another
	// this is based on the number of pairs involving this cluster in the list (n), the number of distand peaks (d)
	// and the probability of a mixed cluster (p). The function uses the cdf values in the similarity model to find
	// a threshold s, s.t. 1-[cdf(s,d)]^n > p
	void computeMinSimilarityForJoiningFromLists(double mixedClusterProb, size_t lowerPos, size_t higherPos);

	Cluster* getCluster(clusterIdx_t clusterIdx) {
		if (clusterIdx<clusterAlloc_[0].getClusterIdx())
			return 0;
		const size_t offset = clusterIdx - clusterAlloc_[0].getClusterIdx();
		if (offset>=nextClusterPos_)
			return 0;
		return &clusterAlloc_[offset];
	}

	const Cluster* getCluster(clusterIdx_t clusterIdx) const {
		if (clusterIdx<clusterAlloc_[0].getClusterIdx())
			return 0;
		const size_t offset = clusterIdx - clusterAlloc_[0].getClusterIdx();
		if (offset>=nextClusterPos_)
			return 0;
		return &clusterAlloc_[offset];
	}

	size_t getCluserPosition(clusterIdx_t clusterIdx) const
	{
		if (clusterIdx<clusterAlloc_[0].getClusterIdx())
			return MAX_SIZE_T;
		const size_t offset = clusterIdx - clusterAlloc_[0].getClusterIdx();
		if (offset>=nextClusterPos_)
			return MAX_SIZE_T;
		return offset;
	}

	float getMaxSimilarityForCluster(Cluster* cluster, clusterIdx_t& otherIdx) const
	{
		const clusterIdx_t* bestSimilarityClusterIdxs = cluster->getBestSimilarityClusterIdxs();
		const float*		bestSimilarityValues	  = cluster->getBestSimilarityValues();
		otherIdx=MAX_CLUSTER_IDX;
		float maxSimilarity=0.0;
		for (size_t i=0; i<NUM_TOP_SIMILARITIES_TO_SAVE; i++)
		{
			if (bestSimilarityValues[i]>maxSimilarity)
			{
				const clusterIdx_t idx = bestSimilarityClusterIdxs[i];
				const Cluster* otherCluster = getCluster(idx);
				if (! otherCluster || ! otherCluster->getIndInPlay())
					continue;

				maxSimilarity = bestSimilarityValues[i];
				otherIdx      = bestSimilarityClusterIdxs[i];
			}
		}
		return maxSimilarity;
	}

private:

	const Config* config_;
	MsParameterStruct* params_;
	const MsSimilarityModel* simModel_;

	string		 outDir_;
	string		 nameStr_;

	int datasetIdx_;
	int batch_;

	size_t nextClusterPos_;	  // the position of the next cluster that can be written
	size_t nextPeakPos_;

	clusterIdx_t runningClusterIdx_; // a unique index given to every spectrum that is read for clustering
	
	clusterIdx_t runningOutputIdx_;  // a unique index given to every new cluster that is outputted (this
									 // cluster was first formed in this run, not a cluster that got added onto)
									
	clusterIdx_t totalSpectraWritten_; // equals the sum of singletons in totalClustersWritten_

	clusterIdx_t totalClustersWritten_;

	int verboseLevel_;

	MsArchive   archive1_, archive2_;

	PeakWeightTable peakWeightTable_; // for creating consensuses

	DatFileManager datFileManager_;

	ClusterOutputter clusterOutputter_;
	float minSqsForSingleton_; // don't output a singleton cluster if it has an sqs below this value

	// storage vectors
	vector<SingleSpectrumHeader>  headerAlloc_;      // original singletion headers are sttored here
	vector<Peak>                  peakAlloc_;        // original singletion peaks are stored here
	vector<Cluster>               clusterAlloc_;     // the cluester structures


	// These are used for maintaining lists for each peak, which spectra have it as one of the top K
	MyVectorAllocator<clusterIdx_t>          clusterIdxAllocation_;

	vector< AllocatedVector<clusterIdx_t>* > newClusterIdxsLists_;  /*! these are the lists that hold for each mass
																     all non-archived spectra that have a peak with that mass
															         in the top k = NUM_PEAKS_FOR_HEURISTIC
															         the first dim is peak mass in Cluster::peakIndexTolerance increments */


	vector<ClusterPositionPair> similarityPairPositions_; // used to store which similarity computations should be performed

	map<string, ClusterClustEntry> clustEntries_; /// used for merging two archives
	ClustFile					   clustFile_;   /// used for mergin two clusters

	// private functions

	clusterIdx_t getFirstClusterIdx() const { return clusterAlloc_[0].getClusterIdx(); }


	// shifts all relevant memory left (pushing away the clusters on the lower m/z range)
	// the sift is done to all storage vectors (peaks, headers, clusters, lists, etc.)
	// this (expensive) action is to be taken before a new m/z slice of spectra is added for 
	// clustering. Every cluster that is removed is outputed (to DAT and/or MGF)
	void shiftStorage(size_t posToMoveToZero, mass_t windowSize);

	// Sets the SQS value needed for an unpaired singleton in order for it to be written in the output.
	// This threshold is set according to the spectrum density
	void setMinSqsForSingleton(float minSimilarity);


	/*! Writes the clusters to the output. This function uses the functions \c writeRegularClustersToOutput
		and \c writeMergedArchiveClustersToOutput (if there are archives involved). */
	void writeClustersToOutput(size_t n);

	// writes the first n clusters in clusterAlloc_ to the appropriate files
	size_t writeRegularClustersToOutput(size_t n, size_t& numSkipped, size_t& numLowSqs);

	// writes the first n clusters (where the clusters can be merged from two archives)
	size_t writeMergedArchiveClustersToOutput(size_t n, size_t& numSkipped, size_t& numLowSqs);



	// obtains the cluster record from the relevant clust files
	bool getClusterEntry(const Cluster& cluster, ClusterClustEntry& clustEntry);

	void UpdateClusterIdxsLists(clusterIdx_t minClusterIdx); // removes all cluster idxs that have a lower

	// make cluster member string. Makes the text entry that should appear in the clust file (without the first line)
	size_t makeClusterMembersStringAndSelectTitleAndPeptide(clusterIdx_t idx, ostringstream& oss, 
									const char*& existingTitle, string& maximalPeptide) const;

	// NP3 make cluster member string. Makes the text entry that should appear in the clust file (without the first line)
	vector<string> makeClusterMembersStringAndSelectTitleAndPeptide_CSV(clusterIdx_t idx,
															const char*& existingTitle, string& maximalPeptide);

	size_t makeClustMemberStringForMergedClusters(const Cluster& cluster, 
								ostringstream& oss, string& existingTitle, string& maximalPeptide, bool verb=false);

	// adds addedIdx to the clsuter mainIdx.
	// performs all the needed actions including updating the info in the added clusters singletons.
	void joinClusters(Cluster* mainCluster, Cluster* addedCluster, int context=0);

	void selectDistancePeaksFromMultipleSigneltons(Cluster* clsuter);

	void makeConsensus(Cluster* cluster, SingleSpectrumHeader* clusterHeader=0);
	void makeConsensusForSmallCluster(Cluster* cluster, SingleSpectrumHeader* clusterHeader, bool indSetConsensusParameters= true);
	void makeConsensusForLargeCluster(Cluster* cluster, SingleSpectrumHeader* clusterHeader);

	void setConsensusParameters(Cluster* cluster, SingleSpectrumHeader* clusterHeader = NULL);

	/*! \brief A debug function that tests that all the m/zs are in a reasonable range */
	bool testSingletonsWithGoodMzs(const Cluster* cluster, mass_t windowSize) const;
};



#endif

