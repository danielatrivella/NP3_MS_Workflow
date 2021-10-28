#include "DatFile.h"
#include "../Common/auxfun.h"


DatFile::~DatFile()
{
	if (indOpen_)
		fclose(stream_);

	if (buffer_)
		delete [] buffer_;

}


void DatFile::init( bool indClearPath)
{
	if (indOpen_)
		fclose(stream_);

	if (! buffer_ || bufferSize_<DAT_BUFFER_SIZE)
	{
		if (buffer_)
			delete [] buffer_;

		buffer_ = new char[DAT_BUFFER_SIZE];
		bufferSize_ = DAT_BUFFER_SIZE;
	}

	if (indClearPath)
	{
		path_.clear();
		fileIdx_ = MAX_UINT;
	}
	numSpectra_ = 0;
	numPeaks_ = 0;
	minMOverZ_ = MAX_FLOAT;
	maxMOverZ_ = 0.0;
	bufferPosition_ = 0;
	bufferEnd_		= 0;
	numSpectraRead_ = 0;

	indOpen_ = false;
}


// opens reads the statistics on number of spectra, peaks, and m/z's and closes
// doesn't allocate buffer
void DatFile::peekAtStats(const char* path)
{
	numSpectra_ = 0;
	numPeaks_ = 0;
	minMOverZ_ = MAX_FLOAT;
	maxMOverZ_ = 0.0;
	bufferPosition_ = 0;
	bufferEnd_		= 0;
	numSpectraRead_ = 0;
	path_ = path;
	
	stream_ = fopen(path_.c_str(), "rb");
	if (! stream_)
		error("Couldn't open dat file for reading: ",path);

	// read first numbers in the file
	char headerBuffer[DAT_HEADER_SIZE+1];
	if (fread(headerBuffer, 1, DAT_HEADER_SIZE, stream_) != DAT_HEADER_SIZE)
	{
		fclose(stream_);
		error("Could not read DAT file correctly: ",path);
	}

	char* p = headerBuffer;
	numSpectra_ = *(reinterpret_cast<clusterIdx_t*>(p));
	p+= sizeof(clusterIdx_t);
	numPeaks_   = *(reinterpret_cast<longInt8_t*>(p));
	p+= sizeof(longInt8_t);
	minMOverZ_  = *(reinterpret_cast<mass_t*>(p));
	p+= sizeof(mass_t);
	maxMOverZ_  = *(reinterpret_cast<mass_t*>(p));
	p+= sizeof(mass_t);

	if (minMOverZ_<0 || maxMOverZ_>20000.0)
		error("Invalid m/z values in dat file: ",path);

	fclose(stream_);

	indOpen_ = false;
}




bool DatFile::openForReading(const char* path)
{
	if (path)
	{
		init();
		path_ = path;
	}
	else
		init(false); // keep existing path
	
	stream_ = fopen(path_.c_str(), "rb");
	if (! stream_)
		error("Couldn't open DAT file for reading: ",path);
	bufferEndStreamPosition_ = 0;

	// read first numbers in the file
	if (fread(buffer_,1,DAT_HEADER_SIZE,stream_) != DAT_HEADER_SIZE)
	{
		fclose(stream_);
		error("Could not read DAT file correctly: ",path);
	}
	bufferEndStreamPosition_+=DAT_HEADER_SIZE;

	char* p = buffer_;
	numSpectra_ = *(reinterpret_cast<clusterIdx_t*>(p));
	p+= sizeof(clusterIdx_t);
	numPeaks_   = *(reinterpret_cast<longInt8_t*>(p));
	p+= sizeof(longInt8_t);
	minMOverZ_  = *(reinterpret_cast<mass_t*>(p));
	p+= sizeof(mass_t);
	maxMOverZ_  = *(reinterpret_cast<mass_t*>(p));
	p+= sizeof(mass_t);

	indOpen_ = true;

	return true;
}

bool DatFile::closeAfterReading()
{
	if (! indOpen_)
		return false;

	fclose(stream_);

	numSpectra_ =0;
	numPeaks_ =0;
	minMOverZ_ = MAX_FLOAT;
	maxMOverZ_ = 0.0;
	bufferPosition_ = 0;

	if (buffer_)
		delete [] buffer_;

	bufferEndStreamPosition_ = 0;
	bufferEnd_ = 0;
	bufferPosition_ = 0;
	buffer_ = 0;

	indOpen_ = false;
	return true;
}



void DatFile::flushBuffer()
{
	if (bufferPosition_ == 0)
		return;

	if (fwrite(buffer_, 1, bufferPosition_, stream_) != bufferPosition_)
		error("Error writing buffer to file ",path_.c_str());

	bufferPosition_ = 0;
}


bool DatFile::openForWriting(const char* path, int verboseLevel)
{
	init();
	path_ = path;
	static int counter =0;

	stream_ = fopen(path,"wb");
	counter++;
	if (! stream_)
	{
		cout << "This is file #" << counter << endl;
		error("Couldn't open file for writing: ",path);
	}
	indOpen_ = true;

	if (verboseLevel>1)
		cout << "\tOPENED " << counter << ": " << path << endl;
	static char bufferMessage[]={"This message should be erased...This message should be erased...This message should be erased..."};
	
	// write space for the 4 numbers (numSpectra_, numPeaks_, minMz_, maxMz_
	fwrite(bufferMessage, 1, DAT_HEADER_SIZE, stream_);

	return true;
}


bool DatFile::closeAfterWriting()
{
	if (! indOpen_)
		return false;

	if (bufferPosition_>0)
		flushBuffer();

	// go back to start of file and write the final stats
	fseek(stream_, 0, SEEK_SET);
	clusterIdx_t* cit=reinterpret_cast<clusterIdx_t*>(buffer_);
	*cit++ = numSpectra_;
	longInt8_t* i8t = reinterpret_cast<longInt8_t*>(cit);
	*i8t++ = numPeaks_;
	mass_t *mt = reinterpret_cast<mass_t*>(i8t);
	*mt++ = minMOverZ_;
	*mt++ = maxMOverZ_;

	fwrite(buffer_, 1, DAT_HEADER_SIZE, stream_);
	fclose(stream_);

	numSpectra_ =0;
	numPeaks_ =0;
	minMOverZ_ = MAX_FLOAT;
	maxMOverZ_ = 0.0;
	bufferPosition_ = 0;
	bufferEndStreamPosition_ =0;

	if (buffer_)
		delete [] buffer_;

	buffer_ = NULL;

	indOpen_ = false;
	return true;
}

size_t DatFile::writePeakList(const PeakList& pl, const SingleSpectrumHeader* newHeader)
{
	assert(indOpen_);
	size_t requiredSpace = sizeof(unsigned int); // length of spectrum in bytes
	requiredSpace += pl.getNumPeaks() * sizeof(Peak) + sizeof(SingleSpectrumHeader);
	requiredSpace += pl.getHeader()->getTitle().length();
	requiredSpace += pl.getHeader()->getPeptideStr().length();
	requiredSpace += 2*sizeof(unsigned short int); // length of title + length of peptide str

	if (bufferSize_ - bufferPosition_ >= requiredSpace + 100)
		flushBuffer();

	if (! newHeader)
		newHeader = pl.getHeader();

	size_t bytesWritten = pl.writeToDatBuffer(buffer_, newHeader);
	bufferPosition_ += bytesWritten;
	addSpectrumToStats(pl.getHeader()->getMOverZ(), pl.getNumPeaks());

	return bytesWritten;
}

void DatFile::addSpectrumToStats(mass_t mz, unsigned int numPeaksInSpectrum)
{
	numSpectra_++;
	numPeaks_ += numPeaksInSpectrum;
	if (mz > maxMOverZ_)
		maxMOverZ_ = mz;
	if (mz < minMOverZ_)
		minMOverZ_ = mz;
}

size_t DatFile::writeSpectrumFromBuffer(const char* buffer, size_t numBytes)
{
	assert(indOpen_);
	return (fwrite(buffer, 1, numBytes, stream_)); 
}


/****************************************************************
This function is used when reading the whole Dat file
*****************************************************************/
bool DatFile::getNextDatSpectrum(char*& specStart, size_t& numBytes, float* sqs, 
								 long* peaksFilePosition, char** peaksBufferPosition)
{
	if (bufferEnd_ == 0 ||
		bufferEnd_ == DAT_BUFFER_SIZE && bufferEnd_ - bufferPosition_ < 131072)
	{
		// shunt end of buffer
		if (bufferPosition_>0)
		{
			memmove(buffer_, buffer_ + bufferPosition_, bufferEnd_ - bufferPosition_);
			bufferEnd_ = bufferEnd_ - bufferPosition_;
			bufferPosition_ = 0;
		}
		const size_t numBytesRead = fread(buffer_ + bufferEnd_, 1, DAT_BUFFER_SIZE - bufferEnd_, stream_);
		bufferEnd_ += numBytesRead;
		bufferEndStreamPosition_ += static_cast<long>(numBytesRead);
	}

	if (bufferEnd_ == 0)
		return false;

	specStart = buffer_ + bufferPosition_;
	unsigned int* ui = reinterpret_cast<unsigned int*>(specStart);
	numBytes = *ui++;
	unsigned int numPeaks = *ui;

	float* fl = reinterpret_cast<float*>(specStart + DAT_OFFSET_OF_SQS);
	if (sqs)
		*sqs = *fl;

	// if needed give exact pointer for seeking later on
	if (peaksFilePosition)
		*peaksFilePosition = bufferEndStreamPosition_ - bufferEnd_ + bufferPosition_
							 + numBytes - numPeaks*sizeof(Peak);
	if (peaksBufferPosition)
		*peaksBufferPosition = buffer_ + bufferPosition_ + numBytes - numPeaks*sizeof(Peak);

	// move to next
	bufferPosition_ += numBytes;

	if (numBytes<50 || numBytes>=131072)
		error("Dat size not in range of expected size (50 to 131072 bytes): ",numBytes);

	numSpectraRead_++;
	return (numSpectraRead_ < numSpectra_);
}

// reads all the header information from a dat file and adds ito to the vector
bool DatFile::readAndAddDatSpectrumStats(vector<DatSpectrumStats>& headerStats, unsigned int fileIdx)
{
	if (! indOpen_)
		openForReading();

	if (numSpectra_ == 0)
		return false;

	fseek(stream_, DAT_HEADER_SIZE, SEEK_SET);

	const size_t statStartIdx = headerStats.size();
	headerStats.resize(headerStats.size()+numSpectra_);

	size_t numSpectraExamined = 0;
	longInt8_t bufferStartOffset    = DAT_HEADER_SIZE; // the dat header
	long	   nextSpectrumPosition = bufferStartOffset;

	while (numSpectraExamined < numSpectra_)
	{
		// use large margin 131072 to avoid possibility that the spectrum
		// overreaces the buffer boundry, as this causes problems

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
			return (numSpectraExamined == numSpectra_);

		assert(bufferPosition_ < bufferEnd_);

		DatSpectrumStats& stats = headerStats[statStartIdx + numSpectraExamined];
		stats.fileIdx      = fileIdx;
		stats.filePosition = bufferStartOffset + bufferPosition_;

		char* specStart = buffer_ + bufferPosition_;
		char* p = specStart;
		const unsigned int spectrumSize = *(reinterpret_cast<unsigned int*>(p)); // add number of bytes to skip
		bufferPosition_ += spectrumSize;

		p+= sizeof(unsigned int);
		stats.numPeaks     = *(reinterpret_cast<unsigned int*>(p));
		p+= sizeof(unsigned int);
		stats.mOverZ	   = *(reinterpret_cast<mass_t*>(p));

		if (spectrumSize>=131072)
		{
			cout << endl << "ERROR: Spectrum size too large, how did this happen?" << endl;
			cout << this->path_ << endl;
			cout << "SIZE     : " << spectrumSize << endl;
			cout << "File idx : " << stats.fileIdx << endl;
			cout << "File pos : " << stats.filePosition << endl;
			cout << "m/z      : " << stats.mOverZ << endl;
			cout << "num peaks: " << stats.numPeaks << endl;
			cout.flush();
		}
		assert(spectrumSize < 131072);
		assert(stats.filePosition == nextSpectrumPosition);

		nextSpectrumPosition += spectrumSize;
		numSpectraExamined++;
	}
	closeAfterReading();

	return (numSpectraExamined == numSpectra_);
}


bool DatFile::readDatFile( const Config* config,
						   PeakList*			 peakListStorage,
						   SingleSpectrumHeader* headersStorage, // 
						   Peak*	   		     peaksStroge)
{
	size_t numSpectraReadFromDat = 0;
	size_t numPeaksReadFromDat = 0;
	longInt8_t bufferStartOffset  = DAT_HEADER_SIZE;

	if (! indOpen_)
		error("Trying to read from closed dat file!");


	fseek(stream_, DAT_HEADER_SIZE, SEEK_SET); // go to first spectrum
	
	while (numSpectraReadFromDat < numSpectra_)
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
			return true;

		assert( numSpectraReadFromDat < numSpectra_ );
		
		// goto start of spectrum
		const char *specStart = buffer_ + bufferPosition_;
		const unsigned int spectrumSize = *(reinterpret_cast<const unsigned int*>(specStart));


		// read spectrum header
		SingleSpectrumHeader& header = headersStorage[numSpectraReadFromDat];
		header.setFileType(IFT_DAT);
		Peak*			peaksTarget  = &peaksStroge[numPeaksReadFromDat];
		header.scanSpectrumHeaderFromBuffer(specStart, config);
		if (! config->getKeepOriginalDatasetIdx())
			header.setDatasetIndex(datasetIndex_); // this is the file number in the dat list
		header.setIndexInFile(static_cast<int>(numSpectraReadFromDat));

		if (header.getScanNumber()<0 && header.getScanNumber() != MIN_INT)
		{
			cout << "Problem with spectrum: " << numSpectraReadFromDat << endl;
			header.printStats();
		}
		assert( header.getScanNumber() == MIN_INT || header.getScanNumber()>= 0);

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

		peakListStorage[numSpectraReadFromDat].setHeader(&headersStorage[numSpectraReadFromDat]);
		peakListStorage[numSpectraReadFromDat].setPeaksPtr(peaksTarget);

		numSpectraReadFromDat++;
		numPeaksReadFromDat += numPeaks;

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

		bufferPosition_ += spectrumSize;
	}

	return (numSpectraReadFromDat == numSpectra_);
}


void DatFile::readDatSpectraToStorage(	const Config* config,
										const vector<DatSpectrumStats>& datStats,
										const int datFileIdx,
										SingleSpectrumHeader*           headersStorage, // 
										Peak*				            peaksStroge,
										size_t& numSpectraReadFromDat,
										size_t& numPeaksReadFromDat, 
										size_t  startDatStatIdx,
										size_t  endDatStatIdx )
{
	numSpectraReadFromDat = 0;
	numPeaksReadFromDat = 0;

	assert(datStats.size());
	if (startDatStatIdx == MAX_SIZE_T)
		startDatStatIdx = 0;
	if (endDatStatIdx == MAX_SIZE_T)
		endDatStatIdx = datStats.size()-1;

	assert(startDatStatIdx<=endDatStatIdx && endDatStatIdx<datStats.size());
	assert( fileIdx_ == MAX_SIZE_T || (datStats[startDatStatIdx].fileIdx == datStats[endDatStatIdx].fileIdx && datStats[startDatStatIdx].fileIdx == fileIdx_) );

	if (! indOpen_)
		error("Trying to read from closed dat file!");

	if (datStats.size() == 0)
		return;

	const size_t numSpectraToRead = endDatStatIdx - startDatStatIdx + 1;
	const longInt8_t startFilePos = datStats[startDatStatIdx].filePosition;
	longInt8_t bufferStartOffset  = startFilePos;

	fseek(stream_, startFilePos, SEEK_SET); // go to first spectrum
	
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

		const size_t datStatsIdx = startDatStatIdx + numSpectraReadFromDat;
		const DatSpectrumStats& ds = datStats[datStatsIdx];

		assert( numSpectraReadFromDat < numSpectraToRead );
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
		Peak*			peaksTarget  = &peaksStroge[datStats[datStatsIdx].peakWritePos];
		header.scanSpectrumHeaderFromBuffer(specStart, config);
		if (! config->getKeepOriginalDatasetIdx())
			header.setDatasetIndex(datasetIndex_); 
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

		numSpectraReadFromDat++;
		numPeaksReadFromDat += numPeaks;

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

		bufferPosition_ += spectrumSize;
	}
}



