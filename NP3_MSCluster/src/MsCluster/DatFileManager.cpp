#include "DatFileManager.h"
#include "MsClusterAuxfuns.h"
#include "MsArchive.h"
#include "Cluster.h"
#include "MsModClusterWriter.h"
#include "../Common/auxfun.h"


void DatFileManager::init(const string& datListPath, const Config* config)
{
	vector<string> datPathsFromList;

	// get dat files set up and sorted according to minMOverZ_
	readListOfPaths(datListPath.c_str(), datPathsFromList);

	init(datPathsFromList, config);
}

	
void DatFileManager::init(const vector<string>& datPathsFromList, const Config* config)
{
	config_ = config;
	maxMOverZ_ = 0;
	totalSpectraToCluster_ = 0;


	datFiles_.clear();
	datFiles_.resize(datPathsFromList.size());

	openedDatFileInds_.clear();
	openedDatFileInds_.resize(datPathsFromList.size(),0);

	clusterIdx_t numSpectra=0;
	mass_t maxMargin = 0.0;
	for (size_t i=0; i<datPathsFromList.size(); i++)
	{
		datFiles_[i].peekAtStats(datPathsFromList[i].c_str());
		numSpectra += datFiles_[i].getNumSpectra();
		if (datFiles_[i].getMaxMOverZ() > maxMOverZ_)
			maxMOverZ_ = datFiles_[i].getMaxMOverZ();
		if (datFiles_[i].getMinMOverZ() < minMOverZ_)
			minMOverZ_ = datFiles_[i].getMinMOverZ();

		if (datFiles_[i].getMaxMOverZ() - datFiles_[i].getMinMOverZ() > 1.0)
			cout << datFiles_[i].getPath() << "\t" << datFiles_[i].getMaxMOverZ() << " - " << datFiles_[i].getMinMOverZ() << " = "
					<< datFiles_[i].getMaxMOverZ() - datFiles_[i].getMinMOverZ() << endl;

		if (datFiles_[i].getMaxMOverZ() - datFiles_[i].getMinMOverZ()>maxMargin)
			maxMargin = datFiles_[i].getMaxMOverZ() - datFiles_[i].getMinMOverZ();
		
	}

	cout << "Max margin in dat files: " << maxMargin << endl;

	const clusterIdx_t maxAllowed = numeric_limits<clusterIdx_t>::max() - 10;
	if (numSpectra > maxAllowed )
	{
		cout << "You are trying to cluster more than " << numSpectra << " spectra in a single batch!" << endl;
		cout << "That is too many, and it exceeds the limit of " << maxAllowed << " per batch." << endl;
		cout << "Please split this clustering job into several batches..." << endl;
		error("Terminating because input size is too large!");
	}

	totalSpectraToCluster_ = numSpectra;

	sort(datFiles_.begin(),datFiles_.end());

	int numBad=0;
	for (size_t i=0; i<datFiles_.size(); i++)
	{
		datFiles_[i].setFileIdx(i);
	/*	if (datFiles_[i].getMinMOverZ()>datFiles_[i].getMaxMOverZ() || datFiles_[i].getMinMOverZ()>20000.0)
		{
			cout << "Bad dat file: " << datFiles_[i].getPath() << endl;
			numBad++;
		}*/
	}
	assert(numBad==0);

	splitIntoBatches();

	nextBatchToOpen_  = 0;
	
	size_t maxBatchSize=0;
	for (size_t i=0; i<datBatches_.size(); i++)
		if (datBatches_[i].numSpectra>maxBatchSize)
		{
			maxBatchSize = datBatches_[i].numSpectra;
			assert(maxBatchSize<999999999);
		}
	unreadSpectra_.reserve(maxBatchSize);
}


void DatFileManager::initFromTwoArchives(const vector<string>& datPaths1, const vector<string>& datPaths2, const Config* config)
{
	config_ = config;
	maxMOverZ_ = 0;
	totalSpectraToCluster_ = 0;

	
	datFiles_.resize(datPaths1.size() + datPaths2.size());
	openedDatFileInds_.clear();
	openedDatFileInds_.resize(datFiles_.size(),0);

	clusterIdx_t numSpectra=0;
	for (size_t i=0; i<datPaths1.size(); i++)
	{
		datFiles_[i].peekAtStats(datPaths1[i].c_str());
		datFiles_[i].setTypeIndex(1);
		numSpectra += datFiles_[i].getNumSpectra();
		if (datFiles_[i].getMaxMOverZ() > maxMOverZ_)
			maxMOverZ_ = datFiles_[i].getMaxMOverZ();
		if (datFiles_[i].getMinMOverZ() < minMOverZ_)
			minMOverZ_ = datFiles_[i].getMinMOverZ();
	}

	for (size_t i=0; i<datPaths2.size(); i++)
	{
		size_t idx = datPaths1.size() + i;
		datFiles_[idx].peekAtStats(datPaths2[i].c_str());
		datFiles_[idx].setTypeIndex(2);
		numSpectra += datFiles_[idx].getNumSpectra();
		if (datFiles_[i].getMaxMOverZ() > maxMOverZ_)
			maxMOverZ_ = datFiles_[idx].getMaxMOverZ();
		if (datFiles_[i].getMinMOverZ() < minMOverZ_)
			minMOverZ_ = datFiles_[idx].getMinMOverZ();
	}

	const clusterIdx_t maxAllowed = numeric_limits<clusterIdx_t>::max() - 10;
	if (numSpectra > maxAllowed )
	{
		cout << "You are trying to cluster more than " << numSpectra << " spectra in a single batch!" << endl;
		cout << "That is too many, and it exceeds the limit of " << maxAllowed << " per batch." << endl;
		cout << "Please split this clustering job into several batches..." << endl;
		error("Terminating because input size is too large!");
	}

	totalSpectraToCluster_ = numSpectra;

	sort(datFiles_.begin(),datFiles_.end());

	int numBad=0;
	for (size_t i=0; i<datFiles_.size(); i++)
	{
		datFiles_[i].setFileIdx(i);
		if (datFiles_[i].getMinMOverZ()>datFiles_[i].getMaxMOverZ() || datFiles_[i].getMinMOverZ()>20000.0)
		{
			cout << "Bad dat file: " << datFiles_[i].getPath() << endl;
			numBad++;
		}
	}
	assert(numBad==0);

	splitIntoBatches();
	nextBatchToOpen_  = 0;
	
	size_t maxBatchSize=0;
	for (size_t i=0; i<datBatches_.size(); i++)
		if (datBatches_[i].numSpectra>maxBatchSize)
		{
			maxBatchSize = datBatches_[i].numSpectra;
			assert(maxBatchSize<999999999);
		}
	unreadSpectra_.reserve(maxBatchSize);
}

// assumes datFiles_ is sorted according to increasing minMOverZ_
void DatFileManager::splitIntoBatches()
{
	datBatches_.clear();
	maxNumSpectraInBatch_ = 0;
	maxNumPeaksInBatch_   = 0;

	for (size_t i=0; i<datFiles_.size(); i++)
	{
		DatBatch batch;
		batch.startIdx = i;
		batch.endIdx   = i;
		batch.numPeaks   = datFiles_[i].getNumPeaks();
		batch.numSpectra = datFiles_[i].getNumSpectra();
		batch.minMz	     = datFiles_[i].getMinMOverZ();
		batch.maxMz		 = datFiles_[i].getMaxMOverZ();

	/*	if (datFiles_[i].getMinMOverZ()>datFiles_[i].getMaxMOverZ())
		{
			cout << datFiles_[i].getMinMOverZ() << "\t" << datFiles_[i].getMaxMOverZ() << endl;
			cout << "Bad: " << datFiles_[i].getPath() << endl;
		}
		assert(batch.minMz<=batch.maxMz);*/


		size_t j;
		for ( j=i+1; j<datFiles_.size(); j++)
		{
			if (datFiles_[j].getMinMOverZ() >= batch.maxMz)
			{	
				datBatches_.push_back(batch);
				break;
			}

			batch.endIdx = j;
			batch.numPeaks += datFiles_[j].getNumPeaks();
			batch.numSpectra += datFiles_[j].getNumSpectra();
			if (datFiles_[j].getMaxMOverZ() > batch.maxMz)
				batch.maxMz = datFiles_[j].getMaxMOverZ();
		}
		if (j == datFiles_.size())
			datBatches_.push_back(batch);

		i=j-1;

	}

	for (size_t i=0; i<datBatches_.size(); i++)
	{
		if (datBatches_[i].numSpectra > maxNumSpectraInBatch_)
			maxNumSpectraInBatch_ = datBatches_[i].numSpectra;

		if (datBatches_[i].numPeaks > maxNumPeaksInBatch_)
			maxNumPeaksInBatch_ = datBatches_[i].numPeaks;
	}
}

/*******************************************************************************
// gets a set of dat spectra which do not exceed the avilable memory constraints
// returns false if there are no more spectra availble
********************************************************************************/
bool DatFileManager::getNewDatSpectraStats(size_t availableSpectraPlaces, 
										   size_t availablePeaks,
										   vector<DatSpectrumStats>& datStats)
{
	size_t assignedSpectra = 0;
	size_t assignedPeaks   = 0;
	size_t datStartIdx = datStats.size();
	
	// first take spectra from the unread ones, if there are any
	size_t numUnreadToUSe=0;
	if (unreadSpectra_.size()>0)
	{
		for (numUnreadToUSe=0; numUnreadToUSe<unreadSpectra_.size(); numUnreadToUSe++)
		{
			if (assignedSpectra == availableSpectraPlaces ||
				assignedPeaks+unreadSpectra_[numUnreadToUSe].numPeaks > availablePeaks)
				break;

			assignedSpectra++;
			assignedPeaks += unreadSpectra_[numUnreadToUSe].numPeaks; 
		}

		assert(numUnreadToUSe>0);
		datStats.resize(datStartIdx + numUnreadToUSe);
		memcpy(&datStats[datStartIdx], &unreadSpectra_[0], sizeof(DatSpectrumStats)*numUnreadToUSe);

		if (numUnreadToUSe==unreadSpectra_.size())
		{
			unreadSpectra_.clear();
		}
		else
		{
			const size_t numLeft = unreadSpectra_.size() - numUnreadToUSe;
			memmove(&unreadSpectra_[0], &unreadSpectra_[numUnreadToUSe], sizeof(DatSpectrumStats)*numLeft);
			unreadSpectra_.resize(numLeft);
		}
	}

	assert(assignedSpectra <= availableSpectraPlaces);
	assert(assignedPeaks   <= availablePeaks );

	// Try adding complete batches
	mass_t minMzAdded = 0.0;
	bool indAddedABatch = false;
	while (nextBatchToOpen_ < datBatches_.size())
	{
		if (! indAddedABatch)
			minMzAdded = datBatches_[nextBatchToOpen_].minMz;

		// We do not want to load more than 40Da slices because this
		// will cause too many ouput files to be open simultaneously
		// and the program might crash becuase there aren't enoug file descriptors
		if ( indAddedABatch && datBatches_[nextBatchToOpen_].minMz - minMzAdded > 39.666)
			break;

		// if there is enough room, copy the whole batch to the datStats
		if (datBatches_[nextBatchToOpen_].numPeaks + assignedPeaks <= availablePeaks &&
			datBatches_[nextBatchToOpen_].numSpectra + assignedSpectra <= availableSpectraPlaces)
		{
			datStats.reserve(datStats.size()+datBatches_[nextBatchToOpen_].numSpectra);
			for (size_t fileIdx =  datBatches_[nextBatchToOpen_].startIdx; 
				        fileIdx <= datBatches_[nextBatchToOpen_].endIdx;
						fileIdx++)
				datFiles_[fileIdx].readAndAddDatSpectrumStats(datStats, fileIdx);

			assignedPeaks   += datBatches_[nextBatchToOpen_].numPeaks;
			assignedSpectra += datBatches_[nextBatchToOpen_].numSpectra;
			nextBatchToOpen_++;
			indAddedABatch = true;
			continue;
		}

		assert(assignedSpectra <= availableSpectraPlaces);
		assert(assignedPeaks   <= availablePeaks );

		// we've added some spectra already, but cannot fit a whole new batch so we'll stop here
		if (assignedSpectra>0)
			break;

		// this batch is too large, but we need to give some spectra anyway so we'll copy them
		// to the unread spectra area and take as many as possible
		assert(unreadSpectra_.size() == 0);
		unreadSpectra_.reserve(datBatches_[nextBatchToOpen_].numSpectra);
		for (size_t fileIdx =  datBatches_[nextBatchToOpen_].startIdx; 
				    fileIdx <= datBatches_[nextBatchToOpen_].endIdx;
				    fileIdx++)
				datFiles_[fileIdx].readAndAddDatSpectrumStats(unreadSpectra_, fileIdx);

		sort(unreadSpectra_.begin(), unreadSpectra_.end());

		size_t numUnreadToUSe=0;
		for (numUnreadToUSe=0; numUnreadToUSe<unreadSpectra_.size(); numUnreadToUSe++)
		{
			if (assignedSpectra == availableSpectraPlaces ||
				assignedPeaks + unreadSpectra_[numUnreadToUSe].numPeaks >=  availablePeaks)
				break;

			assignedSpectra++;
			assignedPeaks += unreadSpectra_[numUnreadToUSe].numPeaks;
		}

		assert(assignedSpectra <= availableSpectraPlaces);
		assert(assignedPeaks  <= availablePeaks );
		assert(numUnreadToUSe < unreadSpectra_.size());
		
		datStats.resize(numUnreadToUSe);
		memcpy(&datStats[datStartIdx], &unreadSpectra_[0], sizeof(DatSpectrumStats)*numUnreadToUSe);
		size_t numLeft = unreadSpectra_.size() - numUnreadToUSe;
		memmove(&unreadSpectra_[0], &unreadSpectra_[numUnreadToUSe], sizeof(DatSpectrumStats)*numLeft);
		unreadSpectra_.resize(numLeft);
		nextBatchToOpen_++; // this batch is all used up (the spectra are either in datStats or in
							// unreadSpectra_. Next time we open a new batch
		break;
	}
	
	assert(assignedSpectra <= availableSpectraPlaces);
	assert(assignedPeaks <= availablePeaks );

	// set the clusterPos and peakPos fields
	setDatStatsRelativePositions(datStats);

	return (nextBatchToOpen_ < datBatches_.size());
}




bool DatFileManager::fillDatStats(vector<DatSpectrumStats>& datStats, mass_t minMz, mass_t maxMz, size_t numPositionsAvailable)
{
	datStats.clear();
	while (nextBatchToOpen_ < datBatches_.size())
	{
		if (datBatches_[nextBatchToOpen_].maxMz < minMz)
		{
			nextBatchToOpen_++;
			continue;
		}

		if (datBatches_[nextBatchToOpen_].minMz > maxMz)
			break;
		
		if (datBatches_[nextBatchToOpen_].numSpectra + datStats.size() > numPositionsAvailable)
			break;

		for (size_t fileIdx =  datBatches_[nextBatchToOpen_].startIdx; fileIdx <= datBatches_[nextBatchToOpen_].endIdx;fileIdx++)
			datFiles_[fileIdx].readAndAddDatSpectrumStats(datStats, fileIdx);

		nextBatchToOpen_++;
	}

	setDatStatsRelativePositions(datStats);
	return (nextBatchToOpen_<datBatches_.size());
}


// returns how many peaks and spectra are needed for the next batch (or to finish up the current
// one if it is already open)
void DatFileManager::getRequirementsForNextBatch(clusterIdx_t& numSpectra, longInt8_t& numPeaks) const
{
	numSpectra=0;
	numPeaks=0;

	if (unreadSpectra_.size()>0)
	{
		numSpectra = unreadSpectra_.size();
		for (size_t i=0; i<unreadSpectra_.size(); i++)
			numPeaks   += unreadSpectra_[i].numPeaks;
		return;
	}
	
	if (nextBatchToOpen_ == datBatches_.size())
		return;

	numSpectra = datBatches_[nextBatchToOpen_].numSpectra;
	numPeaks   = datBatches_[nextBatchToOpen_].numPeaks;
}


// determines the offsets of where to write the cluster and peak information
// the spectra are sorted in datStats according to precursor m/z and assigned positions
void DatFileManager::setDatStatsRelativePositions(vector<DatSpectrumStats>& datStats) const
{
	sort(datStats.begin(),datStats.end());

	size_t nextPeakPos    = 0;
	for (size_t i=0; i<datStats.size(); i++)
	{
		datStats[i].clusterWritePos = i;
		datStats[i].peakWritePos    = nextPeakPos;
		nextPeakPos += datStats[i].numPeaks;
	}
}

bool compDatStasFilePos(const DatSpectrumStats& lhs, const DatSpectrumStats& rhs)
{
	return ( lhs.fileIdx < rhs.fileIdx ||
		    (lhs.fileIdx == rhs.fileIdx && lhs.filePosition <  rhs.filePosition) );
}


void DatFileManager::readSpectraToStorage(SingleSpectrumHeader* headers,
										  Peak*				  peaks,
										  vector<DatSpectrumStats>& datStats,
										  size_t& numSpectraRead,
										  size_t& numPeaksRead)
{
	numSpectraRead=0;
	numPeaksRead=0;

	sort(datStats.begin(), datStats.end(), compDatStasFilePos);
	size_t startIdx=0;

	while (startIdx<datStats.size())
	{
		const size_t fileIdx = datStats[startIdx].fileIdx;
		size_t endIdx = startIdx;
		while (endIdx<datStats.size() && datStats[endIdx].fileIdx == fileIdx)
		{
			endIdx++;
		}
		endIdx--; // always need to step one index back

		assert(datStats[startIdx].fileIdx == datStats[endIdx].fileIdx && datStats[endIdx].fileIdx == fileIdx);

		if (! datFiles_[fileIdx].openForReading())
			error("Couldn't open dat file for reading: ", fileIdx);
		
		size_t numPeaksReadFromDat = 0;
		size_t numSpectraReadFromDat = 0;
		datFiles_[fileIdx].readDatSpectraToStorage(config_,  datStats, fileIdx, headers, peaks, 
			numSpectraReadFromDat, numPeaksReadFromDat,  startIdx, endIdx );

		openedDatFileInds_[fileIdx] = true;

		assert( numSpectraReadFromDat == endIdx - startIdx +1 );

		datFiles_[fileIdx].closeAfterReading();
		startIdx = endIdx+1;
		numSpectraRead += numSpectraReadFromDat;
		numPeaksRead += numPeaksReadFromDat;
	}
}


void DatFileManager::generateListOfClustPaths(vector<string>& clustPaths) const
{
	clustPaths.clear();
	for (size_t i=0; i<openedDatFileInds_.size(); i++)
	{
		if (! openedDatFileInds_[i])
			continue;

		const string& datPath = datFiles_[i].getPath();

		assert(datPath.substr(datPath.length()-3,3) == "dat");
		size_t pos = datPath.find_last_of("/\\");
		assert(pos != string::npos);
		string clustPath = std::string();
		if (datPath.substr(pos-3,3) == "dat")
		{
			clustPath = datPath.substr(0,pos-3);
		}
		else if (datPath.substr(pos-4,3) == "dat")
		{
			clustPath = datPath.substr(0,pos-4);
		}
		else
			error("Could not parse dat path correctly: ",datPath.c_str());

		clustPath += "clust";
		clustPath += datPath.substr(pos,datPath.length() - pos -3);
		clustPath += "clust";

	//	cout << "CLUST PATH: " << clustPath << endl;
		clustPaths.push_back(clustPath);
	}
}



size_t DatFileManager::readSpectraToClusterStorage(SingleSpectrumHeader* headers,
								   				   Cluster*			   clusters,
									               vector<DatSpectrumStats>& datStats)
{
	size_t numSpectraRead=0;
	
	sort(datStats.begin(), datStats.end(), compDatStasFilePos);
	size_t startIdx=0;

	while (startIdx<datStats.size())
	{
		const size_t fileIdx = datStats[startIdx].fileIdx;
		size_t endIdx = startIdx;
		while (endIdx<datStats.size() && datStats[endIdx].fileIdx == fileIdx)
		{
			endIdx++;
		}
		endIdx--; // always need to step one index back

		assert(datStats[startIdx].fileIdx == datStats[endIdx].fileIdx);

		if (! datFiles_[fileIdx].openForReading())
			error("Couldn't open dat file for reading: ", fileIdx);
		
		size_t numSpectraReadFromDat = 0;

		// read dat file
		datFiles_[fileIdx].readDatSpectraToClusterStorage(config_,  datStats, fileIdx,
			headers, clusters, numSpectraReadFromDat, startIdx, endIdx );

		openedDatFileInds_[fileIdx] = true;

		assert( numSpectraReadFromDat == endIdx - startIdx +1 );

		datFiles_[fileIdx].closeAfterReading();
		startIdx = endIdx+1;
		numSpectraRead += numSpectraReadFromDat;
		
	}
	return numSpectraRead;
}


/*! \bried Writes specified spectra to the output (using their dat index and positions)
*/
size_t DatFileManager::writeSpectraToOutput(const string& outDir, const string& outName, vector<DatSpectrumStats>& datStats)
{
	size_t numSpectraWritten=0;
	
	sort(datStats.begin(), datStats.end(), compDatStasFilePos);

	// remove duplicates!
	size_t pos=0;
	for (size_t i=1; i<datStats.size(); i++)
		if (datStats[i].fileIdx != datStats[pos].fileIdx ||
			datStats[i].filePosition != datStats[pos].filePosition)
			datStats[++pos] = datStats[i];
	datStats.resize(pos+1);

	size_t startIdx=0;

	ModClusterWriter spectraWriter;
	string path = outDir + "/" + outName + "_conensus_spectra";
	spectraWriter.init(path);
	void* spectraWriterVoidPtr = reinterpret_cast<void*>(&spectraWriter);

	while (startIdx<datStats.size())
	{
		const size_t fileIdx = datStats[startIdx].fileIdx;
		size_t endIdx = startIdx;
		while (endIdx<datStats.size() && datStats[endIdx].fileIdx == fileIdx)
		{
			endIdx++;
		}
		endIdx--; // always need to step one index back

		assert(datStats[startIdx].fileIdx == datStats[endIdx].fileIdx);

		if (! datFiles_[fileIdx].openForReading())
			error("Couldn't open dat file for reading: ", fileIdx);
		
		size_t numSpectraWrittenFromDat = 0;

		// read dat file
		datFiles_[fileIdx].writeDatSpectraToOutput(config_,  datStats, fileIdx,
			startIdx, endIdx, numSpectraWrittenFromDat, spectraWriterVoidPtr );


		assert( numSpectraWrittenFromDat == endIdx - startIdx +1 );

		datFiles_[fileIdx].closeAfterReading();
		startIdx = endIdx+1;
		numSpectraWritten += numSpectraWrittenFromDat;
		
	}
	return numSpectraWritten;
}



bool DatFileManager::getNextBatch(DatBatch& nextBatch)
{
	if (nextBatchToOpen_ >= datBatches_.size())
		return false;

	nextBatch =  datBatches_[nextBatchToOpen_++];
	return true;
}


size_t DatFileManager::getMaxSpaceNeededToCover(mass_t window) const
{
	mass_t spanSize = 2.0 * window;
	size_t maxSpectra = 0;
	for (size_t i=0; i<datBatches_.size(); i++)
	{
		mass_t mzToCover = datBatches_[i].minMz + spanSize;
		size_t total = 0;

		size_t j=i;
		while (j<datBatches_.size())
		{
			if (datBatches_[j].minMz > mzToCover)
				break;
			total += datBatches_[j].numSpectra;
			j++;
		}

		if (total > maxSpectra)
			maxSpectra = total;
	}
	return maxSpectra;
}

void DatFileManager::printBatches() const
{
	for (size_t i=0; i<datBatches_.size(); i++)
	{
		cout << i << "\t" << datBatches_[i].minMz << "\t" << datBatches_[i].maxMz << "\t" << datBatches_[i].numSpectra << endl;
	}
}


void DatFile::readDatSpectraToClusterStorage(  const Config* config,
									  const vector<DatSpectrumStats>& datStats,
									  const int datFileIdx,
									  SingleSpectrumHeader* headersStorage, // 
									  Cluster*				clusterStorage,
									  size_t& numSpectraReadFromDat,
									  size_t startDatStatIdx,
									  size_t endDatStatIdx )
{
	numSpectraReadFromDat = 0;
	fileIdx_ = datFileIdx;

	assert(datStats.size());
	if (startDatStatIdx == MAX_SIZE_T)
		startDatStatIdx = 0;
	if (endDatStatIdx == MAX_SIZE_T)
		endDatStatIdx = datStats.size()-1;

	assert(startDatStatIdx<=endDatStatIdx && endDatStatIdx<datStats.size());
	assert( fileIdx_ == MAX_SIZE_T || 
		   (datStats[startDatStatIdx].fileIdx == datStats[endDatStatIdx].fileIdx && datStats[startDatStatIdx].fileIdx == fileIdx_) );

	if (! indOpen_)
		error("Trying to read from closed dat file!");

	if (datStats.size() == 0)
		return;

	const size_t numSpectraToRead = endDatStatIdx - startDatStatIdx +1;
	const longInt8_t startFilePos = datStats[startDatStatIdx].filePosition;
	longInt8_t bufferStartOffset  = startFilePos;

	fseek(stream_, startFilePos, SEEK_SET); // go to first spectrum

	Peak* peaksTarget = new Peak[100000];
	
	while (numSpectraReadFromDat<numSpectraToRead)
	{
		const size_t bufferTail = bufferEnd_ - bufferPosition_;
		if (bufferEnd_ == 0 ||
			bufferEnd_ == DAT_BUFFER_SIZE &&  bufferTail < 131072)
		{
			// shunt end of buffer
			if (bufferPosition_>0)
			{
				memmove(buffer_, buffer_ + bufferPosition_, bufferEnd_ - bufferPosition_);
				bufferEnd_ = bufferEnd_ - bufferPosition_;
				bufferStartOffset += bufferPosition_;
				bufferPosition_ = 0;
			}
			bufferEnd_ += fread(buffer_ + bufferEnd_, 1, DAT_BUFFER_SIZE - bufferEnd_, stream_);
		}

		if (bufferEnd_ == 0)
			return;

		const size_t datStatsIdx   = numSpectraReadFromDat + startDatStatIdx;
		const DatSpectrumStats& ds = datStats[datStatsIdx];

		assert( numSpectraReadFromDat <= numSpectraToRead );
		assert( datStats[datStatsIdx].fileIdx == MAX_SIZE_T || datStats[datStatsIdx].fileIdx == fileIdx_ );

		// goto start of spectrum
		const char *specStart = buffer_ + bufferPosition_;
		const unsigned int spectrumSize = *(reinterpret_cast<const unsigned int*>(specStart));

		// check if this spectrum should actually be read (it might need to be skipped if the
		// whole DAT file 1 m/z slice cannot be fit in memory, so only parts of the dat files
		// will get read.
		if (bufferStartOffset + bufferPosition_ < datStats[datStatsIdx].filePosition)
		{
			bufferPosition_ += spectrumSize;
			continue;
		}

		assert( bufferStartOffset + bufferPosition_ == datStats[datStatsIdx].filePosition );
		assert( datStats[datStatsIdx].filePosition - bufferStartOffset < DAT_BUFFER_SIZE );

		// read spectrum header
		SingleSpectrumHeader& header = headersStorage[datStats[datStatsIdx].clusterWritePos];
		header.setFileType(IFT_DAT);
		header.scanSpectrumHeaderFromBuffer(specStart, config);
		if (! config->getKeepOriginalDatasetIdx())
			header.setDatasetIndex(datasetIndex_); 
		header.setIndexInFile(static_cast<int>(numSpectraReadFromDat));
		header.setArchiveSource(typeIndex_);

		// these values overwrite the values read in the files since the file values are not used
		// storing the position and the dat file idx in the header makes it easy to 
		// write the spectra files that are matched
		header.setPositionInFile(datStats[datStatsIdx].filePosition);
		header.setSpectraFileIndexInList(datFileIdx);

		if (header.getScanNumber()<0 && header.getScanNumber() != MIN_INT)
		{
			cout << "Problem with spectrum: " << numSpectraReadFromDat << endl;
			header.printStats();
		}
		assert( header.getScanNumber() == MIN_INT || header.getScanNumber()>= 0);

		// if no file index was written (for example this might be a cluster) then use
		// the datFile index
		if (header.getSpectraFileIndexInList()<0)
			header.setSpectraFileIndexInList(datFileIdx);

		// copy peaks
		const unsigned int numPeaks = header.getOriginalNumPeaks();
		const size_t copySize = numPeaks*sizeof(Peak);
		const char* peaksSource = specStart + spectrumSize - copySize;

		assert(copySize < spectrumSize);
		memcpy(reinterpret_cast<char*>(peaksTarget), peaksSource , copySize);

		// make sure that all peak counts are ok
		if (header.getClusterSize()<=1)
		{
			for (size_t i=0; i<numPeaks; i++)
			{
				peaksTarget[i].count = 1;
				peaksTarget[i].maxPossible = 1;
			}
		}

		// sanity chaeck
		// TODO: place in _DEBUG mode
		bool foundError = (numPeaks == 0 || peaksTarget[0].mass<=0.0 || peaksTarget[0].intensity < 0.0);
		size_t i=0;
		if (! foundError)
			for (i=1; i<numPeaks; i++)
				if (peaksTarget[i].mass<peaksTarget[i-1].mass ||peaksTarget[i].intensity < 0.0)
				{
					foundError=true;
					break;
				}
		if (foundError)
		{
			cout << "Error reading spectrum from Dat buffer!" << endl;
			cout << "Error at peak " << i << endl;
			cout << "Num peaks: " << numPeaks << endl;
			cout << "Dataset index: " << header.getDatasetIndex() << endl;
			cout << "Cluster size : " << header.getClusterSize() << endl;
			header.printStats();
			
			for (size_t i=1; i<numPeaks; i++)
			{
				cout << i << "\t";
				peaksTarget[i].print();
			}
			exit(1);
		}
		// create cluster
		Cluster& cluster = clusterStorage[datStats[datStatsIdx].clusterWritePos];
		cluster.createNewCluster(numSpectraReadFromDat, &header, peaksTarget, numPeaks);

		bufferPosition_ += spectrumSize;
		numSpectraReadFromDat++;
	}

	delete [] peaksTarget;
}


void DatFile::writeDatSpectraToOutput(  const Config* config,
								   const vector<DatSpectrumStats>& datStats,
								   const int datFileIdx,		   
								   size_t startDatStatIdx,
								   size_t endDatStatIdx,
								   size_t& numSpectraReadFromDat,
								   void* modClusterWriterVoidPtr)
{
	ModClusterWriter* clusterWriter = reinterpret_cast<ModClusterWriter*>(modClusterWriterVoidPtr);

	numSpectraReadFromDat = 0;
	fileIdx_ = datFileIdx;

	assert(datStats.size());
	if (startDatStatIdx == MAX_SIZE_T)
		startDatStatIdx = 0;
	if (endDatStatIdx == MAX_SIZE_T)
		endDatStatIdx = datStats.size()-1;

	assert(startDatStatIdx<=endDatStatIdx && endDatStatIdx<datStats.size());
	assert( fileIdx_ == MAX_SIZE_T || 
		   (datStats[startDatStatIdx].fileIdx == datStats[endDatStatIdx].fileIdx && datStats[startDatStatIdx].fileIdx == fileIdx_) );

	if (! indOpen_)
		error("Trying to read from closed dat file!");

	if (datStats.size() == 0)
		return;

	const size_t numSpectraToRead = endDatStatIdx - startDatStatIdx +1;
	const longInt8_t startFilePos = datStats[startDatStatIdx].filePosition;
	longInt8_t bufferStartOffset  = startFilePos;

	fseek(stream_, startFilePos, SEEK_SET); // go to first spectrum

	Peak* peaksTarget = new Peak[100000];
	char* mgfBuffer = new char[1000000]; // should be enough for any MGF
	
	while (numSpectraReadFromDat<numSpectraToRead)
	{
		const size_t bufferTail = bufferEnd_ - bufferPosition_;
		if (bufferEnd_ == 0 ||
			bufferEnd_ == DAT_BUFFER_SIZE &&  bufferTail < 131072)
		{
			// shunt end of buffer
			if (bufferPosition_>0)
			{
				memmove(buffer_, buffer_ + bufferPosition_, bufferEnd_ - bufferPosition_);
				bufferEnd_ = bufferEnd_ - bufferPosition_;
				bufferStartOffset += bufferPosition_;
				bufferPosition_ = 0;
			}
			bufferEnd_ += fread(buffer_ + bufferEnd_, 1, DAT_BUFFER_SIZE - bufferEnd_, stream_);
		}

		if (bufferEnd_ == 0)
			return;

		const size_t datStatsIdx   = numSpectraReadFromDat + startDatStatIdx;
		const DatSpectrumStats& ds = datStats[datStatsIdx];

		assert( numSpectraReadFromDat <= numSpectraToRead );
		assert( datStats[datStatsIdx].fileIdx == MAX_SIZE_T || datStats[datStatsIdx].fileIdx == fileIdx_ );

		// goto start of spectrum
		const char *specStart = buffer_ + bufferPosition_;
		const unsigned int spectrumSize = *(reinterpret_cast<const unsigned int*>(specStart));

		// check if this spectrum should actually be read (it might need to be skipped if the
		// whole DAT file 1 m/z slice cannot be fit in memory, so only parts of the dat files
		// will get read.
		if (bufferStartOffset + bufferPosition_ < datStats[datStatsIdx].filePosition)
		{
			bufferPosition_ += spectrumSize;
			continue;
		}

		if (bufferStartOffset + bufferPosition_ != datStats[datStatsIdx].filePosition )
		{
			cout << "datStatsIdx       = " << datStatsIdx << endl;
			cout << "bufferStartOffset = " << bufferStartOffset << endl;
			cout << "bufferPosition_   = " << bufferPosition_ << endl;
			cout << "filePosition      = " << datStats[datStatsIdx].filePosition << " ( " << bufferStartOffset + bufferPosition_ << ")" << endl;
			for (size_t i=0; i<10 && i<datStats.size(); i++)
				cout << "\t" << i << "\t" << datStats[i].filePosition << endl;
			exit(1);
		}

		assert( bufferStartOffset + bufferPosition_ == datStats[datStatsIdx].filePosition );
		assert( datStats[datStatsIdx].filePosition - bufferStartOffset < DAT_BUFFER_SIZE );

		// read spectrum header
		SingleSpectrumHeader header;
		header.setFileType(IFT_DAT);
		header.scanSpectrumHeaderFromBuffer(specStart, config);	
		header.setDatasetIndex(datasetIndex_); // this is the file number in the dat list
		header.setIndexInFile(static_cast<int>(numSpectraReadFromDat));
		header.setArchiveSource(typeIndex_);

		if (header.getScanNumber()<0 && header.getScanNumber() != MIN_INT)
		{
			cout << "Problem with spectrum: " << numSpectraReadFromDat << endl;
			header.printStats();
		}
		assert( header.getScanNumber() == MIN_INT || header.getScanNumber()>= 0);

		// if no file index was written (for example this might be a cluster) then use
		// the datFile index
		if (header.getSpectraFileIndexInList()<0)
			header.setSpectraFileIndexInList(datFileIdx);

		// copy peaks
		const unsigned int numPeaks = header.getOriginalNumPeaks();
		const size_t copySize = numPeaks*sizeof(Peak);
		const char* peaksSource = specStart + spectrumSize - copySize;

		assert(copySize < spectrumSize);
		memcpy(reinterpret_cast<char*>(peaksTarget), peaksSource , copySize);

		// make sure that all peak counts are ok
		if (header.getClusterSize()<=1)
		{
			for (size_t i=0; i<numPeaks; i++)
			{
				peaksTarget[i].count = 1;
				peaksTarget[i].maxPossible = 1;
			}
		}

		// sanity chaeck
		// TODO: place in _DEBUG mode
		bool foundError = (numPeaks == 0 || peaksTarget[0].mass<=0.0 || peaksTarget[0].intensity < 0.0);
		size_t i=0;
		if (! foundError)
			for (i=1; i<numPeaks; i++)
				if (peaksTarget[i].mass<peaksTarget[i-1].mass ||peaksTarget[i].intensity < 0.0)
				{
					foundError=true;
					break;
				}
		if (foundError)
		{
			cout << "Error reading spectrum from Dat buffer!" << endl;
			cout << "Error at peak " << i << endl;
			cout << "Num peaks: " << numPeaks << endl;
			cout << "Dataset index: " << header.getDatasetIndex() << endl;
			cout << "Cluster size : " << header.getClusterSize() << endl;
			header.printStats();
			
			for (size_t i=1; i<numPeaks; i++)
			{
				cout << i << "\t";
				peaksTarget[i].print();
			}
			exit(1);
		}
		// create cluster
		
		PeakList pl;
		pl.setHeader(&header);
		pl.setPeaksPtr(peaksTarget);
		pl.setNumPeaks(numPeaks);
		size_t n=pl.outputToMgfBuffer(mgfBuffer, &header);
		if (n>=1000000)
			error("MGF buffer in cluster out too small, increase size!");

		clusterWriter->writeToBuffer(mgfBuffer, n);

		bufferPosition_ += spectrumSize;
		numSpectraReadFromDat++;
	}

	delete [] peaksTarget;
	delete [] mgfBuffer;
}

