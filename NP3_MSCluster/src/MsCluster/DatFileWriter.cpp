#include "DatFileWriter.h"
#include "MsClusterAuxfuns.h"
#include "MsParameterStruct.h"
#include "../PepNovo/SpectraList.h"
#include "../PepNovo/PeakList.h"
#include "../PepNovo/ScanList.h"
#include "../PepNovo/MetaList.h"

void DatFileWriter::init( mass_t mzIncrement)
{
	mzIncrement_ = mzIncrement;
	config_ = model_->get_config();
	maxBinIdx_ = static_cast<size_t>(maxMz_ / mzIncrement_)+2;
	datFiles_.clear();
	datFiles_.resize(maxBinIdx_);
	fileCounters_.clear();
	fileCounters_.resize(maxBinIdx_,0);
	fileSizes_.clear();
	fileSizes_.resize(maxBinIdx_,0);
	numSpectraWrittenFirstPass_ = 0;
	numSpectraWrittenSecondPass_ =0;
	numSpectraReadFromOriginalFiles_ = 0; 
	numOriginalPaths_  = 0;

	datPaths_.clear();
}

void DatFileWriter::closeAllOpenDats()
{
	for (size_t i=0; i<datFiles_.size(); i++)
		if (datFiles_[i].closeAfterWriting())
			datPaths_.push_back(datFiles_[i].getPath());
}

DatFileWriter::~DatFileWriter()
{
	closeAllOpenDats();
}



string DatFileWriter::writeDatPaths(bool indSort)
{
	string outPath = datDir_ + "/" + datName_ + "_dat_list.txt";
	FILE* stream = fopen(outPath.c_str(),"w");
	if (! stream)
		error("Couldn't open file for writing: ",outPath.c_str());

	if (indSort)
		sortDatPathsAccordingToMz(datPaths_);

	for (size_t i=0; i<datPaths_.size(); i++)
		fprintf(stream,"%s\n",datPaths_[i].c_str());

	fclose(stream);
	return outPath;
}

void DatFileWriter::addPeakListToDat(const PeakList& pl, bool ignoreMzIdx)
{
	const mass_t mz = pl.getHeader()->getMOverZ();
	const size_t idx = (ignoreMzIdx ? 0 :computeMzIndex(mz, mzIncrement_, 0.0));
	if (! ignoreMzIdx && (idx<1 || idx>= datFiles_.size()))
	{
		cout << endl << "Error with the following header:" << endl;
		pl.getHeader()->printStats();
		cout << "M/Z = " << mz << "   idx = " << idx << " (dat files " << datFiles_.size() << ")" << endl;
		error("m/z value out of range!");
	}

	if (! datFiles_[idx].getIndOpen())
	{
		ostringstream oss;
		if (ignoreMzIdx)
		{
			oss << datDir_ << "/" << datName_ << "_" << fileCounters_[idx] << ".dat";
		}
		else
		{
			size_t roundedMz = computeDatRoundedMz(mz, mzIncrement_);
			oss << datDir_ << "/" << datName_ << "_" << roundedMz << "_" << fileCounters_[idx] << ".dat";
		}

		datFiles_[idx].openForWriting(oss.str().c_str(), verboseLevel_);
	}

	// This is an ugly workaround that plants the new number of peaks in the peak lists header
	const size_t orgNumPeaks = pl.getHeader()->getOriginalNumPeaks();
	SingleSpectrumHeader* header = const_cast<SingleSpectrumHeader*>(pl.getHeader());
	
	assert(pl.getNumPeaks()>0 && pl.getNumPeaks()<100000);
	assert( header->getFileType() != IFT_MZXML || header->getScanNumber() >= 0);
	header->setOriginalNumPeaks(pl.getNumPeaks());
	fileSizes_[idx] += datFiles_[idx].writePeakList(pl);
	header->setOriginalNumPeaks(orgNumPeaks);

	if (fileSizes_[idx]>maxDatFileSize_)
	{
		datPaths_.push_back(datFiles_[idx].getPath());
		datFiles_[idx].closeAfterWriting();
		fileCounters_[idx]++;
		fileSizes_[idx]=0;
	}
}

// writes the spectrum directly from a buffer
void DatFileWriter::writeDatSpectrumDirectly(const char* spec, size_t numBytes)
{
	const unsigned int numPeaks =  *(reinterpret_cast<const unsigned int*>(spec + sizeof(unsigned int)));
	const mass_t mz = *(reinterpret_cast<const mass_t*>(spec + 2*sizeof(unsigned int))); // skip spectrum size and number of peaks
	const size_t idx = computeMzIndex(mz, mzIncrement_, 0.0);

	assert(numPeaks>0);
	
	if (idx<1 || idx>= datFiles_.size())
	{
		cout << "M/Z = " << mz << endl;
		error("Precursor m/z value out of range!");
	}
	
	if (! datFiles_[idx].getIndOpen())
	{
		ostringstream oss;
		size_t roundedMz = computeDatRoundedMz(mz, mzIncrement_);
		oss << datDir_ << "/" << datName_ << "_" << roundedMz << "_" << fileCounters_[idx] << ".dat";
		datFiles_[idx].openForWriting(oss.str().c_str(), verboseLevel_);
	}

	fileSizes_[idx] += datFiles_[idx].writeSpectrumFromBuffer(spec, numBytes);
	datFiles_[idx].addSpectrumToStats(mz, numPeaks);
	
	if (fileSizes_[idx]>maxDatFileSize_)
	{
		datPaths_.push_back(datFiles_[idx].getPath());
		datFiles_[idx].closeAfterWriting();
		fileCounters_[idx]++;
		fileSizes_[idx]=0;
	}
}





/********************************************************************************
Due to the limit on the number of different open file descriptors, the DAT
creation is done in two stages. First we convert the data into DAT files
using a large mz increment (e.g. 25 Da). No qaulity filtration is peformed at 
this stage.
*********************************************************************************/
string DatFileWriter::convertDataToDatFirstPass(const MsParameterStruct* params)
{
	const string& orgList  = (params->spectraListToLoad.length()>0 ? params->spectraListToLoad : params->list);
	const string& metaList = params->metaList;
	const string& datDir   = params->tmpDir;
	const string& datName  = params->outputName;
	float sqsThreshold     = params->sqsThreshold;
	size_t  fileStartIdx   = params->startFileIdx;
	int     verboseLevel   = params->verboseLevel;

	map<string,int> idTitles;
	if (params->gotCreateArchiveFromMgfs)
		readIdsTitleFromIdFile(params, idTitles);

	datDir_  = datDir;
	datName_ = datName + "_R1";
	verboseLevel_ = verboseLevel;

	init(MAJOR_MZ_INCREMENT_FOR_DAT);

	cout << endl << "Pass 1: reading spectra files and writing to dat with " 
			<< MAJOR_MZ_INCREMENT_FOR_DAT << " Da increments." << endl;
	cout         << "----------------------------------------------------------------------" << endl << endl;

	PMCSQS_Scorer* pmcsqsModel = const_cast<PMCSQS_Scorer*>(model_->get_pmcsqs_ptr());
	if (sqsThreshold>0.0 && ! pmcsqsModel->getIndInitializedSqs())
		error("Sqs model not initialized!, need a valid sqs model if using a filtering threshold!");

	if (sqsThreshold>0.0)
		cout << "Filtering spectra with SQS threshold of " << sqsThreshold << endl;

	vector<SinglePath> paths;
	size_t firstFileIdxInList = 0;
	if (orgList.length())
	{
		vector<string> regularPaths;
		firstFileIdxInList  = readListOfPaths(orgList.c_str(), regularPaths);
		if (fileStartIdx == 0 && firstFileIdxInList>0)
			fileStartIdx = firstFileIdxInList;

		numOriginalPaths_ = regularPaths.size();
		if (verboseLevel_>0)
		{
			cout << "Read " << regularPaths.size() << " paths to spectra files." << endl;
			cout << "Converting data to DAT, using m/z increment of " << fixed << setprecision(2) << mzIncrement_ << endl;
		}

		paths.resize(regularPaths.size());
		for (size_t i=0; i<regularPaths.size(); i++)
		{
			paths[i].path = regularPaths[i];
			paths[i].datasetIdx = (params->datasetIdx == MAX_INT ? 0 : params->datasetIdx);
			paths[i].idxInList = fileStartIdx + i;
		}

		createDirIfDoesNotExist(params->outDir.c_str());
		ostringstream oss;
		oss << params->outputStub << "_" << (params->datasetIdx == MAX_INT ? 0 : params->datasetIdx) << "_spec_list.txt";
		writeListOfPaths(oss.str().c_str(), regularPaths);
	}
	else
	{
		assert( metaList.length()>0);
		MetaList ml;
		ml.readMetaList(metaList.c_str());
		ml.writeLists(params->outputName.c_str(), "_spec_list.txt");
		paths = ml.getSinglePaths();
		
	}


	ScanListManager sem;
	if (params->exclusionList.length()>0)
	{
		const size_t numExclusions = sem.initialize(params->exclusionList.c_str(), 
				    								params->minMz - 5.0, 
													params->maxMz + 5.0);
		if (verboseLevel_>0)
			cout << "Read " << numExclusions << " from " << params->exclusionList << endl;
	}

	if (paths.size() == 0)
		return (std::string(""));

	int numFilesWithoutSpectra = 0;
	size_t peakBufferSize = 10000;
	Peak*  peakBuffer = new Peak[peakBufferSize];
	
	numSpectraWrittenFirstPass_      = 0;
	numSpectraReadFromOriginalFiles_ = 0;
	
	map<string,int> numTimes;

	for (size_t i=0; i<paths.size(); i++)
	{
		const double fileStartTime = time(NULL);

		if (verboseLevel_>0)
			cout << i << "\tExtracting from: " << paths[i].path << " ["
				 << paths[i].datasetIdx << " : " << paths[i].idxInList << "]" << endl;

		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(paths[i].path.c_str(), config_, 
			paths[i].datasetIdx ,
			paths[i].idxInList, 
			params->gotOverwriteLocations);

		SpectraList sl(sa);
        //cout << "test2" << endl;
		sl.selectAllAggregatorHeaders();
        //cout << "test3" << endl;
		sl.removeExcludedScans(sem);
        //cout << "test4" << endl;
		if (verboseLevel_>0)
			cout << "\tFound " << sl.getNumHeaders() << " spectra...";

		if (sl.getNumHeaders() == 0)
		{
			numFilesWithoutSpectra++;
			cout << endl << endl;
			cout.flush();
			continue;
		}
		size_t numExtracted =0;
		for (size_t j=0; j<sl.getNumHeaders(); j++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(j);
			if (header->getOriginalNumPeaks()>1e6)
				continue;

			if (header->getOriginalNumPeaks()>= peakBufferSize)
			{
				delete [] peakBuffer;
				peakBufferSize = header->getOriginalNumPeaks()*2;
				peakBuffer = new Peak[peakBufferSize];
			}

			if (params->gotCreateArchiveFromMgfs && idTitles.size()>0)
			{
				if (idTitles.find(header->getTitle()) == idTitles.end())
					continue; // don't write spectra if the id file was supplied and the title is not there
			}
			numTimes[header->getTitle()]++;
		//	if (numTimes[header->getTitle()]>1 && header->getTitle().length()>0)
		//		cout << endl << "Warning: header appears multiple times: " << header->getTitle() << endl;

			// HACK (bad design)
			// save original generation idx and index in list
			// The problem I am trying to solve is how to keep the indexes written
			// in dat files (generation and index in file), yet still be able to read 
			// the file in the current list of paths (which has a different index)
			// solution (for next version)
			// have separate attributes (originaldatasetIdx, originalFileIndex)
			// these never change no matter where the spectrum gets moved!
			int originalDatasetIdx = header->getDatasetIndex();
			int originalIndexInList= header->getSpectraFileIndexInList();

			SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
			nonConstHeader->setSpectraFileIndexInList(0);
			PeakList pl;
			pl.setPeaksPtr( peakBuffer );

            // NP3 GOT change few peaks from 7 to 1 @@@
			if (pl.readPeaksToBuffer(sa, header, peakBuffer) < 1) // if not enough peaks read, skip this spectrum
				continue;

			if (! params->gotOverwriteLocations)
			{
				nonConstHeader->setDatasetIndex(originalDatasetIdx);
				nonConstHeader->setSpectraFileIndexInList(originalIndexInList);
			}
			else
			{
				nonConstHeader->setDatasetIndex(paths[i].datasetIdx);
				nonConstHeader->setSpectraFileIndexInList(paths[i].idxInList);
			}
			pl.initializePeakList(config_, true);
			numSpectraReadFromOriginalFiles_++;

            // NP3 GOT change few peaks from 7 to 1 @@@
			if (pl.getNumPeaks()<1) // don't bother with spectra with too few peaks
				continue;

			if (pl.getNumPeaks()>100000)
			{
				header->printStats();
				cout << "num peaks: " << pl.getNumPeaks() << endl;
				error("Too many peaks in spectrum, something went wrong!");
			}

			if (pmcsqsModel && (sqsThreshold>0.0 || params->gotCorrectPM ))
			{
				size_t maxCharge=0;
				const float sqs = pmcsqsModel->calculateSqsScore(config_, pl, &maxCharge);

				if (sqs<sqsThreshold || maxCharge == 0)
					continue;

				header->setSqs(sqs);
				if (params->gotCorrectPM)
				{
					PmcSqsChargeRes res;
					pmcsqsModel->computeBestMzValuesForCharge(pl, maxCharge, config_->get_pm_tolerance(), res);

					//cout << header->getMOverZ() << " : " << maxCharge << "\t" << res.mz1 << "\t" << res.score1 << "\t" << res.mz2 << "\t" << res.score2 << endl;
					SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
					nonConstHeader->setOriginalPmWith19(header->getMOverZ());
					
					// this is a wrong charge assignment, use original m/z
					if (fabs(res.mz1-header->getMOverZ())>8.0)
					{	
						nonConstHeader->setMOverZ(header->getMOverZ());
						nonConstHeader->setCharge(header->getCharge());
					}
					else
					{
						nonConstHeader->setMOverZ(res.mz1);
						nonConstHeader->setCharge(maxCharge);
					}
				}
			}
			addPeakListToDat(pl);

			numExtracted++;
		}
		numSpectraWrittenFirstPass_ += numExtracted;

		if (verboseLevel_>0)
		{
			const double fileEndTime = time(NULL);
			cout << " Wrote " << numExtracted << " to dat files (this took " 
				<< fileEndTime-fileStartTime << " sec.)" << endl << endl;
			cout.flush();
		}
	}

	closeAllOpenDats();

	if (peakBuffer)
		delete [] peakBuffer;

	// summary
	if (verboseLevel_>0)
	{
		cout << endl << "SUMMARY (first pass):" << endl;
		cout         << "---------------------" << endl;
		cout << "Wrote " << datPaths_.size() << " dat files to " << datDir_ << endl;
		cout << "These files contain " << numSpectraWrittenFirstPass_ 
			 << " spectra (from a total of " << numSpectraReadFromOriginalFiles_ << " that were read)" << endl;
	}

	if (numFilesWithoutSpectra>0)
	{
		cout << endl << "Warning: encountered " << numFilesWithoutSpectra 
			 << " spectra files for which no spectra were read." << endl << endl;
	}
	
	if (numSpectraWrittenFirstPass_ == 0)
		error("Did not write any spectra in first pass! Exiting.");

	// returns the path to the list of created dat files
	return (writeDatPaths());
}


// reads directly from the dat file and writes it in the DatFileWriter
// without performing any processing (used when splitting dat files)
size_t DatFileWriter::writeDatFileToOtherDats(const string& datPath, 
											  const string& datDir, 
											  const string& datName,
											  float sqsThreshold,
											  size_t& totalSpectraRead)
{
	DatFile dat;
	dat.openForReading(datPath.c_str());

	totalSpectraRead = 0;
	size_t numSpectraWritten=0;
	size_t numBytes = 0;
	float  sqs;
	char* spec = 0;
	bool continueReading=true;
	while ( continueReading )
	{
		continueReading = dat.getNextDatSpectrum(spec, numBytes, &sqs);
		totalSpectraRead++;

		// check if the sqs threshold is high enough
		if (sqsThreshold>0.0 && sqs < sqsThreshold)
		{
			continue;
		}

		writeDatSpectrumDirectly(spec, numBytes);
		numSpectraWritten++;
	}

	dat.closeAfterReading();
	return numSpectraWritten;
}


/*******************************************************************************
The second pass of dat creation has two purposes. First, it splits the dat
files into smaller files with 1 Da m/z increments. Second, it performs the quality
filtration and removes spectra with low quality. The filtration can be done in two
ways: either give a fixed qualtiy threhsold like 0.1, and throw away all spectra
that have an SQS score below that value; or another option is to supply the function
with a maximal fraction of data that should be removed e.g., 50%, and thus ensure
that good spectra are not liekly to be thrown away.
********************************************************************************/
string DatFileWriter::convertDataToDatSecondPass(const string& firstPassDatList,
												 const string& datDir, 
												 const string& datName, 
												 float sqsThreshold,
												 int verboseLevel)
{
	datDir_ = datDir;
	datName_ = datName;
	verboseLevel_ = verboseLevel;

	readListOfPaths(firstPassDatList.c_str(), datPaths_);

	vector< vector<string> > datGroups;
	splitDatListIntoSameMzGroups(datPaths_, datGroups);

	cout << endl << "Pass 2: Splitting dat files into smaller files with 1 Da increments." << endl;
	cout		 << "--------------------------------------------------------------------" << endl << endl;

	numSpectraWrittenFirstPass_  = 0;
	numSpectraWrittenSecondPass_ = 0;
	vector<string> allGroupPaths; // collects paths from all gorups
	for (size_t i=0; i<datGroups.size(); i++)
	{
		DatFileWriter groupWriter(model_);
		groupWriter.init(1.0);
		groupWriter.datName_ = datName;
		groupWriter.datDir_  = datDir;

		for (size_t j=0; j<datGroups[i].size(); j++)
		{
			if (verboseLevel_>0)
				cout << i << "\t" << j+1 << "\t" << datGroups[i][j] << endl;

			size_t numRead = 0;
			numSpectraWrittenSecondPass_ += groupWriter.writeDatFileToOtherDats(datGroups[i][j], 
															datDir, datName, sqsThreshold, numRead);
			numSpectraWrittenFirstPass_ += numRead;
		}

		groupWriter.closeAllOpenDats();

		// extract file name from datPaths_;
		for (size_t j=0; j<groupWriter.datPaths_.size(); j++)
			allGroupPaths.push_back(groupWriter.datPaths_[j]);

		if (verboseLevel_>1)
			cout << "All group paths: " << allGroupPaths.size() << endl;
	}
	cout << endl << endl;
	datPaths_ = allGroupPaths;

	string secondPassDatList = writeDatPaths();

	cout << endl << "SUMMARY (second pass):" << endl;
	cout         << "---------------------" << endl;

	if (numSpectraReadFromOriginalFiles_ >0)
		cout << "Read  " << numSpectraReadFromOriginalFiles_ << " spectra from original input files." << endl;
	if (numSpectraWrittenFirstPass_ > 0)
		cout << "Wrote " << numSpectraWrittenFirstPass_ << " spectra into dat files in the first pass." << endl;

	cout << "Wrote " << numSpectraWrittenSecondPass_ << " spectra to dat files in the second pass." << endl;
	if (sqsThreshold>0.0 && numSpectraReadFromOriginalFiles_>0)
	{
		double diff = numSpectraReadFromOriginalFiles_ - numSpectraWrittenSecondPass_;
		assert(diff>=0.0);
		cout << "Removed " << diff << " spectra due to low quality,"
			<< " which is " << setprecision(3) << fixed
			<< diff/static_cast<double>(numSpectraReadFromOriginalFiles_) << " of the data." << endl;
	}
	cout << "The second pass spectra were written to " << datPaths_.size() << " dat files." << endl;
	cout << "(list of dat files written to " << secondPassDatList << ")" << endl;

	// delete temporary files
#ifdef WIN32
	string cmd = "del " + datDir + "\\" + datName + "_R1*.*";
#else
	string cmd = "rm -rf " + datDir + "/" + datName + "_R1*";
#endif
	int returnValue = system(cmd.c_str());

	return secondPassDatList;
}



