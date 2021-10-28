#ifndef __DATFILEMANAGER_H__
#define __DATFILEMANAGER_H__

/*! @file DatFileManager.h
	Holds the class DatFileManager.
*/

#include "../PepNovo/DatFile.h"

class MsArchive; // fwd dclr
class Cluster; // fwd dclr

/*! \struct DatBatch
A DatBatch is a set of dat files that all fall inside a specific min/max m/z range
(all other Dat files are completely out of this range).
*/
struct DatBatch {
	DatBatch() : startIdx(MAX_UINT), endIdx(0), numSpectra(0), numPeaks(0), minMz(-1.0), maxMz(-1.0) {} 

	size_t startIdx;
	size_t endIdx;
	clusterIdx_t numSpectra;
	longInt8_t	 numPeaks;
	mass_t minMz;
	mass_t maxMz;
};

/*! \class DatFileManager
	Holds and managers all DatFiles that need to be read for a clustering job.

	This structure is used for the acutal clustering job, not the dat creation or splitting.
	The dat creation is done with the class DatFileWriter.
*/
class DatFileManager {
public:

	DatFileManager() : config_(0), listOfDatFiles_(""), nextBatchToOpen_(0), 
		minMOverZ_(MAX_FLOAT), maxMOverZ_(0.0), maxPeakMass_(0.0), totalSpectraToCluster_(0),
		maxNumSpectraInBatch_(0), maxNumPeaksInBatch_(0) {}

	void init(const string& datListPath, const Config* config);
	void init(const vector<string>& datPaths, const Config* config);

	void initFromTwoArchives(const vector<string>& datPaths1, const vector<string>& datPaths2, const Config* config);

	/*! \brief returns the next batch of files.

		@param nextBatch holds the entry for the next batch of dat files

		@return true if this is a good batch, false otherwise.
	*/
	 bool getNextBatch(DatBatch& nextBatch);

	 void printBatches() const;


	 const DatFile& getDatFile(size_t datFileIdx) const { return datFiles_[datFileIdx]; }
	 const vector<DatFile>& getDatFiles() const { return datFiles_; }

	// returns how many peaks and spectra are needed for the next batch (or to finish up the current
	// one if it is already open)
	void getRequirementsForNextBatch(clusterIdx_t& numSpectra, longInt8_t& numPeaks) const;

	// gets a set of dat spectra which do not exceed the avilable memory constraints
	// returns false if there are no more spectra availble.
	bool getNewDatSpectraStats(size_t availableSpectraPlaces, 
							   size_t availablePeaks,
							   vector<DatSpectrumStats>& datStats);


	/*! @breif collects all the needed datstats from the current open one
	Makes sure the mz ranges do not exceed maxMz (the minimum of the batch is not above that)
	and the total number of datstats does not exceed numPositionsAvailable
	*/
	bool fillDatStats(vector<DatSpectrumStats>& datStats, mass_t minMz, mass_t maxMz, size_t numPositionsAvailable);

	// determines the offsets of where to write the cluster and peak information
	// the spectra are sorted in datStats according to precursor m/z
	void setDatStatsRelativePositions(vector<DatSpectrumStats>& datStats) const;

	void readSpectraToStorage(SingleSpectrumHeader* headers,
							  Peak*				    peaks,
							  vector<DatSpectrumStats>& datStats,
							  size_t& numSpectraRead,
							  size_t& numPeaksRead);

	/*! @brief reads spectra to clusters according to order defined by datStats.

	Does not store peak lists (only top peaks for similarity computation)
	@return number of clusters read.
	*/
	size_t readSpectraToClusterStorage(SingleSpectrumHeader* headers,
									   Cluster*			   clusters,
									   vector<DatSpectrumStats>& datStats);


	/*! @brief finds that maximal number of spectra that can belong to batches spanning [-window,+window]
	*/
	size_t getMaxSpaceNeededToCover(mass_t window) const;

	/*! \brief Generates list of clust files that corresponds to dat files that were marked as opened
	*/
	void generateListOfClustPaths(vector<string>& clustPaths) const;


	/*! \bried Writes specified spectra to the output (using their dat index and positions)
	*/
	size_t writeSpectraToOutput(const string& outDir, const string& outName, vector<DatSpectrumStats>& datStats);



	mass_t getMaxMOverZ() const { return maxMOverZ_; }
	mass_t getMinMOverZ() const { return minMOverZ_; }
	clusterIdx_t getTotalSpectraToCluster() const { return totalSpectraToCluster_; }

	clusterIdx_t getMaxNumSpectraInBatch() const { return maxNumSpectraInBatch_; }
	longInt8_t   getMaxNumPeaksInBatch() const { return maxNumPeaksInBatch_; }

private:

	const Config* config_;

	string listOfDatFiles_;
	size_t nextBatchToOpen_;
	mass_t minMOverZ_;
	mass_t maxMOverZ_;
	mass_t maxPeakMass_;
	clusterIdx_t   totalSpectraToCluster_;

	clusterIdx_t maxNumSpectraInBatch_;
	longInt8_t	 maxNumPeaksInBatch_;

	vector<DatFile>			 datFiles_;
	vector<DatBatch>		 datBatches_;
	vector<DatSpectrumStats> unreadSpectra_;

	vector<bool>			 openedDatFileInds_;


	void  splitIntoBatches();
	
};



#endif


