#ifndef __SPECTRAAGGREGATOR_H__
#define __SPECTRAAGGREGATOR_H__

#include "SpectraFile.h"


//////////////////////////////////////////////////////////////////////////////
/// Holds all information concerning the input spectra
class SpectraAggregator {
public:

	SpectraAggregator() : config_(NULL), inputListPath_(""), firstFileIdx_(0),
		numTotalSpectra_(0),
		minSpectrumMz_(POS_INF),     maxSpectrumMz_(NEG_INF), 
		minSpectrumCharge_(POS_INF), maxSpectrumCharge_(NEG_INF), 
		currentOpenStream_(NULL),	 currentOpenFileIndex_(-1) {};

	~SpectraAggregator() { if (currentOpenStream_) fclose(currentOpenStream_); }

	mass_t  getMinSpectrumMz() const { return minSpectrumMz_; }
	mass_t  getMaxSepctrumMz() const { return maxSpectrumMz_; }
	int		getMinSpectrumCharge() const { return minSpectrumCharge_; }
	int		getMaxSepctrumCharge() const { return maxSpectrumCharge_; }

	const vector<int>& getSpectraCountsPerCharge() const { return spectraCountsPerCharge_; }
	int		getNumSpectraWithCharge(int c) const { return spectraCountsPerCharge_[c]; }

	const Config* getConfig() const { return config_; }

	int getNumSpectraFiles() const { return spectraFiles_.size(); }
	size_t getNumTotalSpectra() const { return numTotalSpectra_; }

	size_t getFirstFileIdx() const { return firstFileIdx_; }
	void   setFirstFileIdx(size_t i) { firstFileIdx_ = i; }

	const SpectraFile& getSpectraFile(int i) const { return spectraFiles_[i]; }

	int	initializeFromTextFile(const char* inputListPath, const Config* config, int fileIdx=0, bool overwriteExisitngLocations=false);

	int	initializeFromSpectraFilePath(const char* filePath, const Config* config, int datasetIdx=0, int fileIdx=0, bool overwriteExisitngLocations=false);

	int	initializeFromListOfPaths(const vector<string>& paths, const Config* config, int datasetIdx=0, 
		int startFileIdx = 0, bool overwriteExisitngLocations = false);

	int	readPeakList(const SingleSpectrumHeader* header, Peak* peaks) const;

	/// File nmae must be in format ORGXXX
	void setDatasetIdxAccordingToFileName();

private:

	const Config*	config_;

	string	inputListPath_;

	size_t  firstFileIdx_; // in case files are a batch from a larger job

	int		numTotalSpectra_;

	mass_t	minSpectrumMz_;
	mass_t	maxSpectrumMz_;

	int		minSpectrumCharge_;
	int		maxSpectrumCharge_;

	mutable	FILE*	currentOpenStream_;
	mutable int		currentOpenFileIndex_;

	vector<int> spectraCountsPerCharge_;
	vector<int> spectraCountsPerType_;
	vector<SpectraFile>   spectraFiles_;
};


#endif

