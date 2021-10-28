#ifndef __MSCARCHIVE_H__
#define __MSCARCHIVE_H__
#include <unistd.h>

/*! @file MsArchive.h
	Holds the MsArchive class
*/

#include "../Common/includes.h"
#include "MsClusterIncludes.h"

class ClusterOutputter;   // fwd dclr
struct MsParameterStruct; // fwd dclr
class Config;
class AllScoreModels;


/*! @struct DataSetRecord
	\brief Container that holds statistics for every dataset that got added to the archive.
*/
struct DataSetRecord {
	DataSetRecord() : index(0), dateAdded(0), numClusters(0), numSingletons(0), 
		fragTolerance(0.34), precursorTolerance(2.5), precursorPPM(1000), fileList(std::string()), comment(std::string()) {} 

	bool readRecordFromBuffer(const char* buffer);
	void writeRecordToBuffer(ostream& oss, size_t index) const;
	void computeTotalsFromDatasets();

	unsigned int index;
	unsigned int dateAdded;
	longInt8_t	 numClusters;
	longInt8_t	 numSingletons;

	mass_t		 fragTolerance;
	mass_t		 precursorTolerance;
	float		 precursorPPM;
	string		 fileList;
	string		 comment;
};

/*! @class MsArchive
	\brief Container holds archive information:
	-# <tt>Dataset statistics<\tt>. There is a line for each dataset that ever got added to the archive (holds information
	like date, number of clusters, number of singletons, and comments.
	-# archive statistics: number of clusters, singletons, identified clusters.
	-# <tt>Archive file paths<\tt>. All archive files are assumed to reside in a single directory.
	paths to a set of dat files that constitue an archive. Also, provides the
	capability to read the cluster assignment information for the archived clusterss.

	This class manages the information involved in maintaining the spectral archive.
	The archive directory is assumed to conain all the archive files:
	-# <tt>A main text file</tt>.
	-# <tt>Archive ".dat" and ".clust" files</tt>. All dat files are assumed to reside in the same directory as the main file.
	In addition, that same direcotry should contain the ".clust" files (one for each dat file) that
	hold the cluster membership information.
	
	In addition to storing the data, the class also provides the functionallity of reading the clusrter membership
	information and providing for writing the cluster membership information.
*/
class MsArchive {
public:

	MsArchive() : archiveName_("Archive"), archivePath_("."), firstCreated_(std::string()), archiveGeneration_(0),
				  totalClusters_(0), totalSingletons_(0), totalIdentifiedClusters_(0){}


	const string& getArchiveName() const { return archiveName_; }
	const string& getArchivePath() const { return archivePath_; }
	longInt8_t	 getTotalSingletons() const { return totalSingletons_; }
	size_t getNumDatasets() const { return archiveDatasets_.size(); }
	size_t getGeneration() const { return archiveGeneration_; }
	
	void clear();
	

	/*! @brief Reads an archive.

	It is assumed that the archive files are all in the same directory as the main file.
	@param mainArchiveFilePath Full path to the main archive file (the text file that has some parameters
	and a list of the archive file names).
	*/
	bool readArchive(const string& mainArchiveFilePath);

	
	/*! \brief Writes an archive.

	This function writes the information about the archive. The assumption made is that all files are
	written to the same directory (and are later stored in the same directory when the archive is read).

	@param clusterOutputter The ClusterOutputter currently being used to write the cluster data.
	@param params The structure with the run parameters (name, generation, batch, outDir).
	*/
	bool writeArchive(const ClusterOutputter& clusterOutputter, const MsParameterStruct* params, const Config* config);

	void addArchive(const MsParameterStruct* params, MsArchive& other);

	void updateGenerationIfNeeded(string& outputName){ if (outputName == archiveName_) archiveGeneration_++; }

	/*! Exapands the list of names to include the archive path.*/
	void generateArchiveDatPaths(vector<string>& paths) const;

	void printStatistics(const Config* config);


	/*! @brief searches spectra agianst the archive
		
		the query spectra are given via --list. For each spectrum returns upto --max-results clusters
		each having a match p-value of at least --max-pvalue
	*/
	void searchArchive(const MsParameterStruct* params, const AllScoreModels* model) const;


	/*! \brief creates an archive from mgf files and an id file (keeps only spectra with
	ids) */
	void createArchiveFromMgf(MsParameterStruct* params, AllScoreModels* model);

private:
	string archiveName_;
	string archivePath_; // all archive files are assumed to reside in the same directory
	string firstCreated_;
	
	// These attributes only get written when the archive file is written
	// if archive is merged/added to in batches, this needs to be done
	// after all batches are done.
	int			 archiveGeneration_;
	longInt8_t	 totalClusters_;
	longInt8_t	 totalSingletons_;
	longInt8_t	 totalIdentifiedClusters_;

	vector<DataSetRecord> archiveDatasets_;
	vector<string>        datNames_; // file names without directory path prefix and without ".dat" or ".clust" suffixes
};


void correctDatError(const char* datList, const Config* config);




#endif

