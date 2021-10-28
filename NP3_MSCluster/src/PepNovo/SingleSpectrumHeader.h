#ifndef __SINGLESPECTRUMHEADER_H__
#define __SINGLESPECTRUMHEADER_H__

#include "PepNovo_includes.h"

// forward declarations
class Peptide; 
class Config;
struct Peak;

///////////////////////////////////////////////////////////////////////////////////
/// Base clasee for single spectrum headers
/// Contains all information required to pinpoint a specific spectrum (file index 
/// and file position). In addition there is a bunch of other info that might be
/// useful. This is class should not be instantiated, but only used as a pointer
/// to other classes


class SingleSpectrumHeader {
public:

	SingleSpectrumHeader() : mOverZ_(NEG_INF), originalNumPeaks_(0), originalPmWith19_(NEG_INF), pmWith19_(MIN_FLOAT),
						charge_(0), fileType_(0), msLevel_(2), archiveSource_(0), scanNumber_(MIN_INT), indexInFile_(-1), positionInFile_(0), 
						clusterSize_(1), spectraFileIndexInList_(MIN_INT), retentionTimeMin_(-1.0),retentionTimeMax_(-1.0),
						peakArea_(0.0), peakId_(std::string()), datasetIndex_(MIN_INT), retentionTime_(0), precursorIntensity_(0),
						firstPeakMass_(-1.0), sqs_(NEG_INF), title_(std::string()), peptideStr_(std::string()) {}

	mass_t	getMOverZ()		const { return mOverZ_; }
	void	setMOverZ(mass_t mz)	  { mOverZ_ = mz; }

	unsigned int getOriginalNumPeaks()	const { return originalNumPeaks_; }
	void	setOriginalNumPeaks(int n) { originalNumPeaks_ = static_cast<unsigned int>(n); }
	
	mass_t	getOriginalPmWith19()	const { return originalPmWith19_; }
	void	setOriginalPmWith19(mass_t m) { originalPmWith19_ = m; }

	mass_t	getPmWith19()	const { return pmWith19_; }
	void	setPmWith19(mass_t m) { pmWith19_ = m; }

	short	getCharge()		const { return charge_; }
	void	setCharge(int c)	{ charge_ = static_cast<short>(c); }
	
	short	getFileType()		const { return fileType_; }
	void	setFileType(int t) { fileType_ = static_cast<short>(t); }
	
	int		getScanNumber() const { return scanNumber_; }
	void	setScanNumber(int n) { scanNumber_ = n; }

	int		getIndexInFile() const { return indexInFile_; }
	void	setIndexInFile(int i) { indexInFile_ = i; }

	short   getArchiveSource() const { return archiveSource_; }
	void    setArchiveSource(short s) { archiveSource_ = s; }
	
	int		getClusterSize() const { return clusterSize_; }
	void	setClusterSize(int i) { clusterSize_ = i; }
	
	int		getSpectraFileIndexInList()  const { return spectraFileIndexInList_; }
	void	setSpectraFileIndexInList(int i)	{ spectraFileIndexInList_ = i; }


	int		getDatasetIndex() const { return datasetIndex_; }
	void    setDatasetIndex(int i)  { datasetIndex_ = i; }
	void    setDatasetIndex(int i) const  { datasetIndex_ = i; }
	
	long	getPositionInFile()	const	{ return positionInFile_; }
	void	setPositionInFile(long i) { positionInFile_ = i; }
	
	float	getRetentionTime()	const { return retentionTime_; }
	void	setRetentionTime(float f) { retentionTime_ = f; }

	// NP3 GOT rt min and max
    float	getRetentionTimeMin()	const { return retentionTimeMin_; }
    void	setRetentionTimeMin(float f) { retentionTimeMin_ = f; }
    float	getRetentionTimeMax()	const { return retentionTimeMax_; }
    void	setRetentionTimeMax(float f) { retentionTimeMax_ = f; }
	float	getPeakArea()	const { return peakArea_; }
	void	setPeakArea(float f) { peakArea_ = f; }

    const string&	getPeakId()	const { return peakId_; }
    void  setPeakId_(const string& s) { peakId_ = s; }

    float	getPrecursorIntensity()	const { return precursorIntensity_; }
	void	setPrecursorIntensity(float f) { precursorIntensity_ = f; }

	short	getMsLevel() const { return msLevel_; }
	void	setMsLevel(int n) { msLevel_ = static_cast<short>(n); }

	mass_t	getFirstPeakMass() const { return firstPeakMass_; }
	void	setFirstPeakMass(mass_t m) { firstPeakMass_ =m; }
	
	float	getSqs()	const { return sqs_; }
	void	setSqs(float f) { sqs_ = f; }
	void	setSqs(float f) const { sqs_ = f; }
	
	const string&	getTitle()	const { return title_; }
	void  setTitle(const string& s) { title_ = s; }		
	
	const string&	getPeptideStr()	const { return peptideStr_; }
	void  setPeptideStr(const string& p) { peptideStr_ = p; }


	void printStats(const Config *config, ostream& os = cout, bool printEndl = true) const;
	void printStats( ostream& os = cout, bool printEndl = true) const;

	bool operator< (const SingleSpectrumHeader& rhs) const
	{
		return ( spectraFileIndexInList_ <  rhs.spectraFileIndexInList_ ||
				(spectraFileIndexInList_ == rhs.spectraFileIndexInList_ &&
				 positionInFile_ < rhs.positionInFile_) ) ;
	}

	/// Scans the spectrum header and extracts all information about the spectrum
	/// (except for the actual peak list). This is function relies on  
	/// separate implementations for each file type that can be read.
	bool	scanSpectrumHeader(FILE *stream, const Config *config);

	bool	scanSpectrumHeaderFromBuffer(const char* buffer, const Config *config);

	size_t  writeHeaderToDatBuffer(char* buffer) const;
	size_t  writeHeaderToMgfBuffer(char* buffer) const;

	//SingleSpectrumHeader& operator= (const SingleSpectrumHeader& rhs);

private:

	bool	scanDtaSpectrumHeader(FILE *stream, const Config *config);
	bool	scanMgfSpectrumHeader(FILE *stream, const Config *config);
	bool	scanMzxmlSpectrumHeader(FILE *stream, const Config *config);
	bool	scanPklSpectrumHeader(FILE *stream, const Config *config);
	bool	scanMs2SpectrumHeader(FILE *stream, const Config *config);
	bool	scanDatSpectrumHeaderFromBuffer(const char* buffer, const Config *config);

	mass_t	     mOverZ_;
	unsigned int originalNumPeaks_; /// value before filtering

	mass_t	originalPmWith19_;
	mass_t	pmWith19_;
	
	short	charge_;
	short   fileType_;       /// IFT_DTA, IFT_MGF, IFT_MZXML, IFT_DAT, ...
	short	msLevel_;
	short	archiveSource_;	 /// used when merging two arcives
	int		scanNumber_;	 /// the value of the scan field in the mzXML, or SCANS= in mgf
	int		indexInFile_;	 /// e.g., = number of BEGIN IONS before this one
	long	positionInFile_; /// position pointer in stream
	int		clusterSize_;
	int		spectraFileIndexInList_; /*! The index of the original spectrum file (MGF or mzXML)
									     For DAT files this value needs to be set manually if the Dat index is needed
							             (for instance if running de novo on a DAT file as input) */

	mutable int datasetIndex_; /*! For clustering jobs, this can be set to different values so we can tell 
							    what list/generation the spectrum came from. */

	float	retentionTime_;
	float	precursorIntensity_;

	// NP3 GOT rt min and max
    float	retentionTimeMin_;
    float	retentionTimeMax_;
    float peakArea_;
    string peakId_;
	
	mutable mass_t	firstPeakMass_;    /// for debug purposes
	mutable float	sqs_;
	string  title_;
	string  peptideStr_; /// holds the peptide structure if it exists
};



#endif


