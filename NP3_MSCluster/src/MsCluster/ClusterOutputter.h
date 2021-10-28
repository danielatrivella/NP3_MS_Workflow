#ifndef __CLUSTEROUTPUTTER_H__
#define __CLUSTEROUTPUTTER_H__

#include "Cluster.h"
#include "MsModClusterWriter.h"
#include "../PepNovo/DatFile.h"
#include "../Common/BufferedOutputFile.h"

struct MsParameterStruct; // fwd dclr

const clusterIdx_t sizeThresholdsForReport[]={1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,100,200,500,1000000000};
const clusterIdx_t numSizeThresholdsForReport=sizeof(sizeThresholdsForReport)/sizeof(clusterIdx_t);


struct OutputFilesForMzSlice {
public:
	OutputFilesForMzSlice() : datPath_(std::string()), mgfPath_(std::string()), clustPath_(std::string()),
		verboseLevel_(0), maxNumClustersPerFile_(20000), 
		runningFileIdx_(MAX_SIZE_T), numClustersInFile_(0), indWriteMgf_(0), indWriteDat_(0), 
		indIsInitialized_(0), indIsOpen_(0) {}

	~OutputFilesForMzSlice() { close(); }

	void init(const MsParameterStruct* params, const string& nameWithMassLabel);

	void open();

	void writeClusterToOutputFiles(const Cluster& cluster, 
								   const string& clustFileEntry,
								   const SingleSpectrumHeader* header,
								   bool  outputThisClusterToMod);
	void close();

	void setMaxNumClustersPerFile(size_t n) { maxNumClustersPerFile_ = n; }

	bool getIndIsInitialized() const { return indIsInitialized_; }
	bool getIndIsOpen() const { return indIsOpen_; }
	size_t getRunningFileIdx() const { return runningFileIdx_; }
	const vector<string>& getDatPaths() const { return datPaths_; }
	const vector<string>& getMgfPaths() const { return mgfPaths_; }
	const vector<string>& getClustPaths() const { return clustPaths_; }

	static ModClusterWriter modClusterWriter; 

private:

	string datPath_, mgfPath_, clustPath_;

	int verboseLevel_;

	size_t maxNumClustersPerFile_;
	size_t runningFileIdx_;
	size_t numClustersInFile_;
	
	bool indWriteMgf_;
	bool indWriteDat_;
	bool indIsInitialized_;
	bool indIsOpen_;
	bool indOutputModClusters_;

	BufferedOutputFile mgfOut_;
	BufferedOutputFile clustOut_;
	DatFile datOut_;

	vector<string> datPaths_, mgfPaths_, clustPaths_;

};







class ClusterOutputter {
public:

	ClusterOutputter() : params_(0), mzIncrement_(MAJOR_MZ_INCREMENT_FOR_DAT), numClustersWritten_(0), totalSpectraCount_(0),
		indWasInitialized_(0), lastIndexClosed_(0) {}

	~ClusterOutputter() { closeAll(); }

	void init(const MsParameterStruct* params, mass_t maxMOverZ);
				
	void outputCluster(const Cluster& cluster, 
					   const string&  clustFileEntry,
					   bool  indCloseOutputFiles,
					   bool	 outputToMod);

	// closes files and collects paths
	void closeAll();

	/*! \breif Closes all open ouptut files below a certain mass (gives a couple indexes room for safety)
		
		@param mz The precursor m/z value of the spectra outputted this round.
	*/
	void closeAllWithMassBelow(mass_t mz);

	void reportClusterStats() const;

	const vector<string>& getAllClustPaths() const { return allClustPaths_; }
	const vector<string>& getAllDatPaths() const { return allDatPaths_; }
	const vector<string>& getAllMgfPaths() const { return allMgfPaths_; }

	clusterIdx_t getNumClustersWritten() const { return numClustersWritten_; }
	longInt8_t   getTotalSpectraCount() const { return totalSpectraCount_; }

private:
	
	const MsParameterStruct* params_;
	int				verboseLevel_;
	int			    numToAddToDsIdx_; // When merging two archives the dataset idxs in the clust listings needs to be updated
									 // to conform with the merged dataset idxs (this number is added only to the idxs of the
									// spectra in the second archive, and equals the number of datasets in the first archive.
	mass_t			mzIncrement_;
	clusterIdx_t	numClustersWritten_;
	longInt8_t		totalSpectraCount_;

	vector<clusterIdx_t> clusterCounts_;
	vector<clusterIdx_t> singletonCounts_;
	vector<longInt8_t>   totalSizeCounts_; // includes sizes accumlated from multiple datasets

	
	bool indWasInitialized_;

	size_t lastIndexClosed_; // the last slice index to be closed

	vector<OutputFilesForMzSlice> outputFiles_;

	vector<string> allDatPaths_, allMgfPaths_, allClustPaths_; // paths of files that were created
};




#endif

