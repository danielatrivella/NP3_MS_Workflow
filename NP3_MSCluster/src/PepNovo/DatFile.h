#ifndef __DATFILE_H__
#define __DATFILE_H__

#include "../PepNovo/PeakList.h"

class SingleSpectrumHeader; // fwd dclr
class Cluster;				// fwd dclr
struct Peak;


const size_t DAT_HEADER_SIZE = sizeof(clusterIdx_t) + sizeof(longInt8_t) + 2*sizeof(mass_t);
const size_t DAT_BUFFER_SIZE = 1179648;

// offset to the position of the sqs field in the header of a DAT spectrum
// this location can be computed by looking in SingleSpectrumHeader::scanDatSpectrumHeaderFromBuffer()
const size_t DAT_OFFSET_OF_SQS = 2*sizeof(unsigned int) + 4*sizeof(mass_t) + 3*sizeof(short) 
								 + 4*sizeof(int) + 2*sizeof(float);



struct DatSpectrumStats {
public:

	// two spectra can have the same m/z (since there is only limited precision)
	// this can mess things up when they are written (violate some of the asserts)
	// so to avoid this problem we can use a more sophisticated sort that will not
	// assign a random to spectra with the same m/z
	bool operator< (const DatSpectrumStats& rhs) const
	{
		return (mOverZ <  rhs.mOverZ ||
			    mOverZ == rhs.mOverZ && fileIdx <rhs.fileIdx ||
				mOverZ == rhs.mOverZ && fileIdx==rhs.fileIdx && filePosition <rhs.filePosition);
	}


	unsigned int fileIdx;
	long         filePosition;
	mass_t       mOverZ;
	unsigned int numPeaks;
	size_t clusterWritePos; // filled in after sorting the stats
	size_t peakWritePos;   // filled in after sorting the stats 
};




class DatFile {
public:

	DatFile() : indOpen_(0), fileIdx_(MAX_UINT), datasetIndex_(MIN_INT), typeIndex_(0),
		path_(""), numSpectra_(0),  numPeaks_(0), minMOverZ_(MAX_FLOAT), maxMOverZ_(0.0), 
		bufferEndStreamPosition_(0), bufferSize_(0), bufferPosition_(0), bufferEnd_(0), 
		numSpectraRead_(0), stream_(0), buffer_(0) {}

	~DatFile();

	bool operator< (const DatFile& rhs) const
	{
		return (minMOverZ_ < rhs.minMOverZ_);
	}

	// reads all the header information from a dat file and adds ito to the vector
	bool readAndAddDatSpectrumStats(vector<DatSpectrumStats>& headerStats, unsigned int fileIdx);

	bool readDatFile( const Config* config,
						   PeakList*			 peakListStorage,
						   SingleSpectrumHeader* headersStorage, // 
						   Peak*	   		     peaksStroge);

	void readDatSpectraToStorage( const Config* config,
								  const vector<DatSpectrumStats>& datStats,
								  const int datFileIdx,
								  SingleSpectrumHeader* headersStorage, // 
								  Peak*				    peaksStroge,
								  size_t& numSpectraReadFromDat,
								  size_t& numPeaksReadFromDat, 
								  size_t startIdx,
								  size_t endIdx );

	/*! @brief reads the dat file into clusters

	Peaks are not retained, only used to init clusters and then discarded
	*/
	void readDatSpectraToClusterStorage(  const Config* config,
										  const vector<DatSpectrumStats>& datStats,
										  const int datFileIdx,
										  SingleSpectrumHeader* headersStorage, // 
										  Cluster*				clusterStroge,
										  size_t& numSpectraReadFromDat,
										  size_t startIdx,
										  size_t endIdx );

	void writeDatSpectraToOutput(  const Config* config,
								   const vector<DatSpectrumStats>& datStats,
								   const int datFileIdx,		   
								   size_t startIdx,
								   size_t endIdx,
								   size_t& numSpectraReadFromDat,
								   void* modClusterWriter);

	// reads spectrum to buffer directly without processing, reports the number of bytes
	// If this is the last spectrum in the file or there is problem returns false, otherwise true.
	bool getNextDatSpectrum(char*& spec, size_t& numBytes, float* sqs=NULL, 
							long* peakPosition=NULL, char** peaksBufferPosition=NULL);

	size_t writePeakList(const PeakList& pl, const SingleSpectrumHeader* newHeader=0);
	size_t writeSpectrumFromBuffer(const char* buffer, size_t numBytes);

	// opens reads the statistics on number of spectra, peaks, and m/z's and closes
	// doesn't allocate buffer
	void peekAtStats(const char* path);

	bool openForReading(const char* path = 0);
	bool closeAfterReading();

	bool openForWriting(const char* path, int verboseLevel = 0);
	bool closeAfterWriting();

	bool getIndOpen() const { return indOpen_; }

	void addSpectrumToStats(mass_t mz, unsigned int numPeaks);
	
	const string& getPath() const { return path_; }
	void  setPath(const string& p) { path_ = p; }

	size_t getFileIdx() const { return fileIdx_; }
	void   setFileIdx(size_t i) { fileIdx_ = i; }
	int    getDatasetIndex() const { return datasetIndex_; }
	void   setDatasetIndex(int i)   { datasetIndex_ = i; }
	longInt8_t		 getNumPeaks() const { return numPeaks_; }
	clusterIdx_t getNumSpectra() const { return numSpectra_; }
	mass_t	getMinMOverZ() const { return minMOverZ_; }
	mass_t  getMaxMOverZ() const { return maxMOverZ_; }

	short	getTypeIndex() const { return typeIndex_; }
	void	setTypeIndex(short s)  { typeIndex_ = s; }
	
private:
	bool indOpen_;

	string path_;
	size_t fileIdx_;
	int	   datasetIndex_;
	short  typeIndex_;
	
	clusterIdx_t numSpectra_;
	longInt8_t	 numPeaks_;
	mass_t minMOverZ_;
	mass_t maxMOverZ_;

	// buffer for reading/writing spectra
	long   bufferEndStreamPosition_; // holds current location of the input file stream
	size_t bufferSize_;
	size_t bufferPosition_;
	size_t bufferEnd_;
	size_t numSpectraRead_;
	FILE* stream_;
	char*  buffer_;
	

	void flushBuffer();
	void init(bool indClearPath = true);
};





#endif



