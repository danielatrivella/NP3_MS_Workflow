#ifndef __SPECTRAFILE_H__
#define __SPECTRAFILE_H__

#include "SingleSpectrumHeader.h"


const size_t XML_BUFFER_SIZE = 524288;
const size_t XML_BUFFER_HALF_SIZE = 262144;

/*!
	\class SpectraFile
	A class for holding the spectra information (headers) for a single input file.
*/
class SpectraFile {
public:

	SpectraFile() : fileType_(-1),		
					minSpectrumMz_(POS_INF),	maxSpectrumMz_(NEG_INF),
					minSpectrumCharge_(POS_INF),maxSpectrumCharge_(NEG_INF), 
					filePath_(std::string()),   unzippedFilePath_(std::string())
	{
					spectraCountsPerCharge_.clear();
					headers_.clear();
	}

	~SpectraFile();

	int		getFileType()	const { return fileType_; }
	int		getNumHeaders()	   const { return headers_.size(); }
	mass_t  getMinSpectrumMz() const { return minSpectrumMz_; }
	mass_t  getMaxSepctrumMz() const { return maxSpectrumMz_; }
	int		getMinSpectrumCharge() const { return minSpectrumCharge_; }
	int		getMaxSepctrumCharge() const { return maxSpectrumCharge_; }

	const char*	getFilePath() const { return (unzippedFilePath_.length()>0 ? unzippedFilePath_.c_str() : filePath_.c_str()); }

	const vector<int>& getSpectraCountsPerCharge() const { return spectraCountsPerCharge_; }
	
	const vector<SingleSpectrumHeader>& getHeaders() const { return headers_; }

	const SingleSpectrumHeader* getHeader(int i) const { return &headers_[i]; }

	/// This function is used to do the initial read ("once over") and collect
	/// the spectra information (and put it in headers)
	int	scanFile(const char *filePath, int datasetIdx, int fileIndexInList, 
		const Config* config, bool removeDuplicates = true, bool overwriteExisitngLocations = false);


	/// This function reads the peaks from the spectrum file (list of <mass,intensity>
	/// It is supplied with a header with all needed information.
	/// The peak buffer *peaks* is assumed to have enough memory to read the numebr
	/// of peaks specified in *header*.
	/// This function relies on a separate implementation for each file type
	int	readPeakList(FILE* stream, const SingleSpectrumHeader* header, 
					 Peak* peaks, const Config* config) const;

	void setAlldatasetIdxs(int idx);

private:


	int	scanDtaFile(const char *filePath, int datasetIdx, int fileIndexInList, const Config* config, bool removeDuplicates = true);
	int	scanMgfFile(const char *filePath, int datasetIdx, int fileIndexInList, const Config* config, bool removeDuplicates = true);
	int	scanMzxmlFile(const char *filePath, int datasetIdx, int fileIndexInList, const Config* config, bool removeDuplicates = true);
	int	scanDatFile(const char *filePath, int datasetIdx, int fileIndexInList, const Config* config, bool indOverwriteDatLocation = false);
	int	scanPklFile(const char *filePath, int fileIndexInList, Config* config);
	int	scanMs2File(const char *filePath, int datasetIdx, int fileIndexInList, const Config* config, bool removeDuplicates = true);

	int	readDtaPeakList(FILE* stream, const SingleSpectrumHeader* header, Peak* peaks, const Config* config) const;
	int	readMgfPeakList(FILE* stream, const SingleSpectrumHeader* header, Peak* peaks, const Config* config) const;
	int	readMzxmlPeakList(FILE* stream, const SingleSpectrumHeader* header, Peak* peaks, const Config* config) const;
	int	readDatPeakList(FILE* stream, const SingleSpectrumHeader* header, Peak* peaks, const Config* config) const;
	int	readPklPeakList(FILE* stream, const SingleSpectrumHeader* header, Peak* peaks, Config* config);
	int	readMs2PeakList(FILE* stream, const SingleSpectrumHeader* header, Peak* peaks, const Config* config) const;

	void	tallySpectraStats();
	
	int		fileType_;
	mass_t	minSpectrumMz_,		maxSpectrumMz_;
	int		minSpectrumCharge_, maxSpectrumCharge_;

	string	filePath_;

	string  unzippedFilePath_; // used in cases where there is a zip file for example

	vector<int> spectraCountsPerCharge_;
	vector<SingleSpectrumHeader> headers_;
};




#endif


