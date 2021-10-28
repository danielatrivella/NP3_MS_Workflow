#ifndef __CLUSTFILE_H__
#define __CLUSTFILE_H__

/*! \file ClustFile.h
	\brief holds the information of the cluster membership for an output file (dat or mgf)
*/

#include "MsClusterIncludes.h"

/*! @struct SingletonClustEntry
	Holds membership information for a singleton in a cluster, does not hold
	peptide sequence (if one is written)
*/
struct SingletonClustEntry {
	SingletonClustEntry() : datasetIdx(0), fileIdx(0), scanNumber(MIN_INT), mz(0.0),
		similarityToConsensus(0.0), pValue(1.0), charge(0), peptide(std::string()) {}

	bool operator< (const SingletonClustEntry& rhs) const
	{
		return (datasetIdx<rhs.datasetIdx ||
			    (datasetIdx==rhs.datasetIdx && fileIdx<rhs.fileIdx) ||
				(datasetIdx==rhs.datasetIdx && fileIdx==rhs.fileIdx && scanNumber<rhs.scanNumber));
	}

	int datasetIdx;
	int fileIdx;
	int scanNumber;
	mass_t mz;
	float  similarityToConsensus;
	float  pValue;
	int    charge;
	string peptide;
};


/*! @struct ClusterClustEntry
	Holds the information for creating  a cluster's clust memeber entry
*/
struct ClusterClustEntry {
	ClusterClustEntry() : clusterIdx(MAX_CLUSTER_IDX), clusterSize(0), mz(0.0), charge(0), title(std::string()), peptide(std::string()), 
		entries(0) {}
	clusterIdx_t clusterIdx;
	clusterIdx_t clusterSize;
	mass_t		 mz;
	int			 charge;
	string		 title;
	string		 peptide;
	vector<SingletonClustEntry> entries;
};


class ClustFile {
public:
	ClustFile() : datasetIdx_(MIN_INT), fileIdx_(MIN_INT), path_(std::string()), fileSize_(0),
		 buffer_(0), bufferSize_(0), currentPos_(0), bufferEnd_(0) {}
	
	~ClustFile();

	bool read(const char* path, int datasetIdx, int fileIdx);

	char* getNextEntryPointer(string* titleString=NULL);

	bool parseEntry(char* pointer, ClusterClustEntry& entry) const;

	int getDatasetIdx() const { return datasetIdx_; }
	int getFileIdx() const    { return fileIdx_; }

private:
	int datasetIdx_; /// genertion index for the cluster file
	int fileIdx_;	 /// file index in the list
	string path_;	 /// actual full path to file

	size_t fileSize_;
	char*  buffer_;
	size_t bufferSize_;
	char*  currentPos_;
	char* bufferEnd_;
};




#endif

