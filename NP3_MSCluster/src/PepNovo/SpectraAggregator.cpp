#include "SpectraAggregator.h"
#include "PepNovo_auxfun.h"

int SpectraAggregator::initializeFromTextFile(const char* inputListPath, const Config *config, int fileIdx,
											  bool overwriteExisitngLocations)
{
	config_ = config;
	vector<string> paths;

	firstFileIdx_ = readListOfPaths(inputListPath, paths);

	return (initializeFromListOfPaths(paths, config, fileIdx, overwriteExisitngLocations));
}



int	SpectraAggregator::initializeFromSpectraFilePath(const char* filePath, 
													 const Config *config,
													 int datasetIdx,
													 int fileIdx,
													 bool overwriteExisitngLocations)
{
	config_ = config;
	vector<string> paths;
	paths.push_back(filePath);

	return (initializeFromListOfPaths(paths, config, datasetIdx, fileIdx, overwriteExisitngLocations));
}



int	SpectraAggregator::initializeFromListOfPaths(const vector<string>& paths,
												 const Config* config,
												 int datasetIdx,
												 int startFileIdx,
												 bool overwriteExisitngLocations)
{
	config_ = config;

	int numFilesWithSpectra=0;

	spectraCountsPerType_.clear();
	spectraCountsPerType_.resize(IFT_NUM_FILE_TYPES,0);
	spectraCountsPerCharge_.clear();


	spectraFiles_.clear();
	spectraFiles_.resize(paths.size());

	for (int i=0; i<paths.size(); i++)
	{
		SpectraFile& sf = spectraFiles_[i];
		const int numHeadersRead = sf.scanFile(paths[i].c_str(), datasetIdx, i+startFileIdx, config_,
											   true, overwriteExisitngLocations);
		if (numHeadersRead>0)
		{
			++numFilesWithSpectra;
		}
		else
			continue;

		numTotalSpectra_ += numHeadersRead;
		spectraCountsPerType_[sf.getFileType()] +=numHeadersRead;

		if (sf.getMaxSepctrumCharge() > maxSpectrumCharge_)
		{
			maxSpectrumCharge_ = sf.getMaxSepctrumCharge();
			spectraCountsPerCharge_.resize(maxSpectrumCharge_+1,0);
		}

		if (sf.getMinSpectrumCharge() < minSpectrumCharge_)
			minSpectrumCharge_ = sf.getMinSpectrumCharge();

		if (sf.getMaxSepctrumMz() > maxSpectrumMz_)
			maxSpectrumMz_ = sf.getMaxSepctrumMz();

		if (sf.getMinSpectrumMz() < minSpectrumMz_)
			minSpectrumMz_ = sf.getMinSpectrumMz();

		const vector<int>& fileSpectraPerCharge = sf.getSpectraCountsPerCharge();
		int charge;
		for (charge=0; charge<fileSpectraPerCharge.size(); charge++)
			if (fileSpectraPerCharge[charge]>0)
				spectraCountsPerCharge_[charge] += fileSpectraPerCharge[charge];
	}


	return numFilesWithSpectra;
}



int	SpectraAggregator::readPeakList(const SingleSpectrumHeader* header, Peak* peaks) const
{
	const int fileIndex = header->getSpectraFileIndexInList();
	if (fileIndex <0)
		error("Illegal file index read from header: ",fileIndex);

	assert(fileIndex<spectraFiles_.size());

	if (fileIndex != currentOpenFileIndex_)
	{
		if (currentOpenStream_)
			fclose(currentOpenStream_);

		if (header->getFileType() == IFT_MZXML || header->getFileType() == IFT_DAT)
		{
			currentOpenStream_= fopen(spectraFiles_[fileIndex].getFilePath(),"rb");
		}
		else
			currentOpenStream_= fopen(spectraFiles_[fileIndex].getFilePath(),"rb");

		if (! currentOpenStream_)
		{
			cout << "Errorz: couldn't open spectrum file for reading: " << 
				spectraFiles_[fileIndex].getFilePath() << endl;
			cout << "File index = " << fileIndex << endl;
			exit(1);
		}
		currentOpenFileIndex_ = fileIndex;
	}
	if ( fseek(currentOpenStream_, header->getPositionInFile(),0) )
	{
		cout << "Error: could not skip to position " << header->getPositionInFile() << endl;
		cout << "Spectra file: " << spectraFiles_[fileIndex].getFilePath() << endl;
		exit(1);
	}

	return (spectraFiles_[fileIndex].readPeakList(currentOpenStream_, header, peaks, config_));
}



/// File nmae must be in format ORGXXX
void SpectraAggregator::setDatasetIdxAccordingToFileName()
{
	for (size_t i=0; i<spectraFiles_.size(); i++)
	{
		int idx =0;
		string name;
		getFileNameWithoutExtension(spectraFiles_[i].getFilePath(),name);

		if (sscanf(name.c_str(),"ORG%d",&idx) == 1)
			spectraFiles_[i].setAlldatasetIdxs(idx);
	}
}



