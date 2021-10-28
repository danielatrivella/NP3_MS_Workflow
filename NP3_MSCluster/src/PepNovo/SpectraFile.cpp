#include "SpectraFile.h"
#include "DatFile.h"
#include "BasicDataStructs.h"
#include "PepNovo_auxfun.h"
#include "base64.h"

#ifndef WIN32
#include <unistd.h>
#endif

SpectraFile::~SpectraFile()
{
	if (unzippedFilePath_.length()>0)
	{
		unlink(unzippedFilePath_.c_str());
	}
}

void SpectraFile::tallySpectraStats()
{
	spectraCountsPerCharge_.resize(6,0);
	minSpectrumCharge_ = POS_INF;
	maxSpectrumCharge_ = NEG_INF;
	minSpectrumMz_ = POS_INF;
	maxSpectrumMz_ = NEG_INF;
	
	int i;
	for (i=0; i<headers_.size(); i++)
	{
		const SingleSpectrumHeader& ssh = headers_[i];
		if (ssh.getCharge()>maxSpectrumCharge_)
		{
			maxSpectrumCharge_ = ssh.getCharge();
			if (maxSpectrumCharge_ >= spectraCountsPerCharge_.size())
				spectraCountsPerCharge_.resize(maxSpectrumCharge_+1,0);
		}
		else if (ssh.getCharge() < minSpectrumCharge_)
			minSpectrumCharge_=  ssh.getCharge();
		
		spectraCountsPerCharge_[ssh.getCharge()]++;

		if (ssh.getMOverZ() < minSpectrumMz_)
			minSpectrumMz_ = ssh.getMOverZ();

		if (ssh.getMOverZ() > maxSpectrumMz_)
			maxSpectrumMz_ = ssh.getMOverZ();
	}
}


int SpectraFile::scanFile(const char *filePath, int datasetIdx, int fileIndexInList, 
						  const Config* config, bool removeDuplicates, bool overwriteExisitngLocations)
{
	fileType_ = getFileExtensionType(filePath);
	filePath_ = std::string(filePath);

	// if this is a zip file, try unzipping it and using the unzipped file instead
	unzippedFilePath_ = std::string();
	if (fileType_ == IFT_ZIP)
	{	
		vector<string> unzippedPaths;
		size_t numUnzips = unzipSingleFile(filePath_, unzippedPaths);
		if (numUnzips == 0)
		{
			cout << "Warning: Could not unzip file correctly: " << filePath << endl;
			return 0;
		}
		if (numUnzips > 1)
			error("Zip file should contain only one spectrum file: ", filePath);

		unzippedFilePath_ = unzippedPaths[0];
		fileType_ = getFileExtensionType(unzippedFilePath_.c_str());
	}

	if (fileType_ == IFT_MGF)
	{
		scanMgfFile(filePath, datasetIdx, fileIndexInList, config, removeDuplicates);
	}
	else if (fileType_ == IFT_MZXML)
	{
		scanMzxmlFile(filePath, datasetIdx, fileIndexInList, config);
	}
	else if (fileType_ == IFT_DAT)
	{
		scanDatFile(filePath, datasetIdx, fileIndexInList, config, overwriteExisitngLocations);
	}
	else if (fileType_ == IFT_DTA)
	{
		scanDtaFile((unzippedFilePath_.length() ? unzippedFilePath_.c_str() : filePath), 
					 datasetIdx, fileIndexInList, config, removeDuplicates);
	}
	else if (fileType_ == IFT_MS2)
	{
		scanMs2File(filePath, datasetIdx, fileIndexInList, config, removeDuplicates);
	}
	else
		error("File type not supported for: ", filePath);

	tallySpectraStats();
	return (headers_.size());
}


int	SpectraFile::scanDtaFile(const char *filePath, int datasetIdx, int fileIndexInList, const Config* config,
							 bool removeDuplicates)
{
	FILE* dtaStream=fopen(filePath,"r");
	if (! dtaStream)
	{
		cout << "Error: couldn't open dta file for reading: |" << filePath << "|" << endl;
		exit(1);
	}

	headers_.clear();
	int counter=0;
	int numSkipped =0;
	while (true)
	{
		SingleSpectrumHeader ssh;
		ssh.setFileType(IFT_DTA);
		
		if (! ssh.scanSpectrumHeader(dtaStream, config))
		{
			char tmpBuff[256];
			if (! fgets(tmpBuff,256,dtaStream)) // move forwards, there was a problem with that position
				break;							// reached eof
				
			continue;
		}

		ssh.setFileType(IFT_DTA);
		ssh.setSpectraFileIndexInList(fileIndexInList);
		ssh.setDatasetIndex(datasetIdx);
		ssh.setIndexInFile(counter++);

		if (ssh.getScanNumber() == MIN_INT)
			ssh.setScanNumber(ssh.getIndexInFile());

		if (removeDuplicates && 
			headers_.size()>0 && 
			headers_.back().getScanNumber()>=0 &&
			headers_.back().getScanNumber() == ssh.getScanNumber())
		{
			numSkipped++;
			continue;
		}

		headers_.push_back(ssh);
	}
	fclose(dtaStream);

	if (numSkipped>0)
		cout << " [ Skipped " << numSkipped << " doubles ] ";
	
	return counter;
}

int	SpectraFile::scanMs2File(const char *filePath, int datasetIdx, int fileIndexInList, const Config* config,
							 bool removeDuplicates)
{
	FILE* ms2Stream=fopen(filePath,"r");
	if (! ms2Stream)
	{
		cout << "Error: couldn't open dta file for reading: |" << filePath << "|" << endl;
		exit(1);
	}

	headers_.clear();
	int counter=0;
	int numSkipped =0;
	while (true)
	{
		SingleSpectrumHeader ssh;
		ssh.setFileType(IFT_MS2);
		
		if (! ssh.scanSpectrumHeader(ms2Stream, config))
		{
			char tmpBuff[256];
			if (! fgets(tmpBuff,256,ms2Stream)) // move forwards, there was a problem with that position
				break;							// reached eof
				
			continue;
		}

		ssh.setFileType(IFT_MS2);
		ssh.setSpectraFileIndexInList(fileIndexInList);
		ssh.setDatasetIndex(datasetIdx);
		ssh.setIndexInFile(counter++);
		if (ssh.getScanNumber() == MIN_INT)
			ssh.setScanNumber(ssh.getIndexInFile());

		if (removeDuplicates && 
			headers_.size()>0 && 
			headers_.back().getScanNumber()>=0 &&
			headers_.back().getScanNumber() == ssh.getScanNumber())
		{
			numSkipped++;
			continue;
		}

		headers_.push_back(ssh);
	}
	fclose(ms2Stream);

	if (numSkipped>0)
		cout << " [ Skipped " << numSkipped << " doubles ] ";
	
	return counter;
}

int SpectraFile::scanMgfFile(const char* filePath, int datasetIdx, int fileIndexInList, const Config* config,
							 bool removeDuplicates)
{
	FILE* mgfStream=fopen(filePath,"rb");
	if (! mgfStream)
	{
		cout << "Error: couldn't open mgf file for reading: |" << filePath << "|" << endl;
		exit(1);
	}

	headers_.clear();
	int counter=0;
	while (true)
	{
		SingleSpectrumHeader ssh;
		ssh.setFileType(IFT_MGF);
		
		if (! ssh.scanSpectrumHeader(mgfStream, config))
		{
			char tmpBuff[256];
			if (! fgets(tmpBuff,256,mgfStream)) // move forwards, there was a problem with that position
				break;							// reached eof
				
			continue;
		}

		if (ssh.getOriginalNumPeaks() == 0)
		{
			counter++;
			continue;
		}

		ssh.setFileType(IFT_MGF);
		ssh.setSpectraFileIndexInList(fileIndexInList);
		ssh.setDatasetIndex(datasetIdx);
		ssh.setIndexInFile(counter++);

		if (ssh.getScanNumber() == MIN_INT)
			ssh.setScanNumber(ssh.getIndexInFile());

		if (removeDuplicates && 
			headers_.size()>0 && 
			headers_.back().getScanNumber()>=0 &&
			headers_.back().getScanNumber() == ssh.getScanNumber())
			continue;

		headers_.push_back(ssh);
	}
	fclose(mgfStream);
	
	return counter;
}


/***********************************************************************************
This is a special function designed to overcome parsing problems I have
with mzXML in Linux enviornments. The function serially extracts spectra from
an mzXML file and stores the peak lists pairs (mass,intensity). This means
that only a single mzXML should be processed at once.
************************************************************************************/
int SpectraFile::scanMzxmlFile(const char* filePath, int datasetIdx, int fileIndexInList, const Config* config,
							   bool removeDuplicates)
{
    static char* buffer = 0;
	if (! buffer)
		buffer = (char*)calloc(XML_BUFFER_SIZE + 1, sizeof(char));
 
	FILE* mzxmlStream=fopen(filePath,"rb");
	if (! mzxmlStream)
	{
		cout << "Error: couldn't open mzxml file for reading: " << filePath << endl;
		exit(1);
	}
	char* scanStartPtr = 0;
	char* bufferEnd = buffer;
	size_t bufferStartOffset =0;
	size_t currentFilePos=0;

	bool doneReading=false;
    while (true)
    {
        // Read more data, to fill up the buffer:
     	if ( currentFilePos == 0)
		{
			scanStartPtr = buffer;
			size_t bytesRead = fread(bufferEnd, sizeof(char), XML_BUFFER_SIZE, mzxmlStream);
			doneReading = (bytesRead<XML_BUFFER_SIZE);

			bufferEnd = buffer + bytesRead;
			*bufferEnd = '\0';
			bufferStartOffset = 0;
			currentFilePos += bytesRead;

			if (bytesRead<50)
				break;
		}
		else
		{
			// try shunt half of the buffer
			if (bufferEnd-buffer > XML_BUFFER_HALF_SIZE &&
				scanStartPtr - buffer > XML_BUFFER_HALF_SIZE)
			{
				const size_t bytesLeft = bufferEnd - scanStartPtr;
				memmove(buffer, scanStartPtr, bytesLeft);
				bufferEnd = buffer + bytesLeft;

				scanStartPtr = buffer;
				const size_t bytesToRead = XML_BUFFER_SIZE - (bufferEnd - buffer);
				const size_t bytesRead = fread(bufferEnd, sizeof(char), bytesToRead, mzxmlStream);

				doneReading = (bytesRead<bytesToRead);
				if (bytesRead>0)
				{
					bufferEnd += bytesRead;
					*bufferEnd = '\0';
					currentFilePos += bytesRead;
					bufferStartOffset = currentFilePos - (bufferEnd - buffer);
				}
			}
		}

		if (doneReading && scanStartPtr > bufferEnd - 50)
			break;

        // Look for a new <scan tag opening:
        char *lastPos = bufferEnd - 5;
		char *pos = scanStartPtr-1;
		while (++pos<lastPos)
		{
			if (*pos != '<')
				continue;
			if (*(pos+1)=='s' && *(pos+2)=='c' && *(pos+3)=='a' && *(pos+4)=='n')
				break;
		}
		char* scanStr =  (pos<lastPos) ? pos : 0;
        if (!scanStr )
        {
			scanStartPtr = lastPos+1;
            continue;
        }

        char* scanNumberStr = strstr(scanStr, "num=");
		const int scanNumber = (scanNumberStr ? parseIntFromXml(scanNumberStr) : MIN_INT);

		char * retentionTimeStr = strstr(scanStr,"retentionTime=\"PT");
		const mass_t retentionTime = (retentionTimeStr ? parseMassFromXml(retentionTimeStr) : -1.0);
		
		char* peakCountStr = strstr(scanStr, "peaksCount=\"");
		if (!peakCountStr)
		{
			cout << "Warninig: bad parsing of peaks from mzxml! " << endl;
			cout << "Scan: " << scanNumber << "  Pos: " << currentFilePos << endl;
			scanStartPtr += 50;
			continue;
			
		}
		const int peakCount = parseIntFromXml(peakCountStr);
	
        char* msLevelStr = strstr(scanStr, "msLevel=");
        if (! msLevelStr)
        {
            cout << "Warning: mzXML parser encountered a scan with no MS level!" << endl; 
			cout << "Scan: " << scanNumber << "  Pos: " << currentFilePos << endl;
            scanStartPtr += 50;
			continue;
        }
        const int msLevel = parseIntFromXml(msLevelStr);
		if (msLevel <= 1)
		{
			scanStartPtr = msLevelStr + 20;
			continue;
		}
        
		char* precursorIntensityStr = strstr(scanStr,"precursorIntensity=");
		if (! precursorIntensityStr)
		{
			scanStartPtr += 50;
			continue;
		}
		const intensity_t precursorIntensity = parseMassFromXml(precursorIntensityStr);
		
		char* precursorMzStr = strstr(scanStr, "<precursorMz");
		if (!precursorMzStr && msLevel > 1)
		{
		
			cout << "Warning: mzXML parser encountered a scan with no m/z" << endl;
			cout << "Scan: " << scanNumber << "  Pos: " << currentFilePos << endl;
            scanStartPtr += 50;
			continue;
		}

		precursorMzStr = strstr(precursorMzStr, ">");
		const mass_t precursorMZ = parseMassFromXml(precursorMzStr);

		// NP3 GOT change numPeaks from 7 to 1
        if (msLevel <= 1 || scanNumber < 0 || peakCount < 1)
		{
			scanStartPtr += 50;
			continue;
		}
        
		char* peakStr = strstr(precursorMzStr, "<peaks");
		if (! peakStr)
		{
			scanStartPtr += 50;
			continue;
		}

		SingleSpectrumHeader ssh;

		ssh.setFileType(IFT_MZXML);
		ssh.setSpectraFileIndexInList(fileIndexInList);
		ssh.setDatasetIndex(datasetIdx);
		ssh.setScanNumber(scanNumber);
		ssh.setMOverZ(precursorMZ);
		ssh.setRetentionTime(retentionTime);

		// NP3 set rt min and max
        ssh.setRetentionTimeMin(retentionTime - config->get_rt_tolerance());
        ssh.setRetentionTimeMax(retentionTime + config->get_rt_tolerance());

		ssh.setPrecursorIntensity(precursorIntensity);
		ssh.setOriginalNumPeaks(peakCount);
		ssh.setMsLevel(msLevel);
		ssh.setPositionInFile(static_cast<long>(bufferStartOffset + (peakStr - buffer)));

		bool skip =(removeDuplicates && 
					headers_.size()>0 && 
					headers_.back().getScanNumber()>=0 &&
					headers_.back().getScanNumber() == ssh.getScanNumber());
		if (! skip)
			headers_.push_back(ssh);
		
		scanStartPtr = peakStr + 8 * peakCount;

		if (scanStartPtr>=buffer+XML_BUFFER_SIZE)
			scanStartPtr = buffer+XML_BUFFER_SIZE-1;

		assert(fileIndexInList>=0 && scanNumber>=0);
    }
	fclose(mzxmlStream);

	return (headers_.size());
}


/*******************************************************************************
This file scans a DAT file and stores its headers.
Should be used when DAT file is read as input for non-clustering applications
(de novo, tags, etc.). It should manually replace the spectraFileIndexInList_ field
of the header to the one given as an argument to the function.
Since Dat files are a special case, we use the DatFile class to do the processing
********************************************************************************/
int	SpectraFile::scanDatFile(const char *filePath, int datasetIdx, int fileIndexInList, 
							 const Config* config, bool indOverwriteDatLocation)
{
	DatFile dat;

	dat.openForReading(filePath);
	
	size_t totalSpectraRead = 0;
	size_t numBytes = 0;
	float  sqs;
	char* spec = 0;
	bool continueReading = (dat.getNumSpectra()>0);
	while ( continueReading )
	{
		long peaksPosition=0;
		continueReading = dat.getNextDatSpectrum(spec, numBytes, &sqs, &peaksPosition);

		SingleSpectrumHeader ssh;
		ssh.setFileType(IFT_DAT);
		ssh.scanSpectrumHeaderFromBuffer(spec, config);
		ssh.setIndexInFile(totalSpectraRead++);
		ssh.setPositionInFile(peaksPosition);

		// overwrites the fileIndexInList that was written in the dat file
		// uses the current index (the index of the dat file in the dat file list).
		if (indOverwriteDatLocation)
		{
			ssh.setSpectraFileIndexInList(fileIndexInList);
			ssh.setDatasetIndex(datasetIdx);
		}
		
		headers_.push_back(ssh);
	}

	dat.closeAfterReading();
	return (headers_.size());
}



int	SpectraFile::readPeakList(FILE* stream, const SingleSpectrumHeader* header, 
							 Peak* peaks, const Config* config) const
{

	int numPeaksRead =0;
	if (fileType_ == IFT_MGF)
	{

		numPeaksRead = readMgfPeakList(stream, header, peaks, config);

	}
	else if (fileType_ == IFT_MZXML)
	{
		numPeaksRead = readMzxmlPeakList(stream, header, peaks, config);
	}
	else if (fileType_== IFT_DAT)
	{
		numPeaksRead = readDatPeakList(stream, header, peaks, config);
	}
	else if (fileType_== IFT_MS2)
	{
		numPeaksRead = readMs2PeakList(stream, header, peaks, config);
	}
	else if (fileType_== IFT_DTA)
	{

		numPeaksRead = readDtaPeakList(stream, header, peaks, config);
	}
	else
		error("Could not read peak list for file tpye: ",fileType_);

	intensity_t totalIntensity=0;
	int i;
	for (i=0; i<numPeaksRead; i++)
	{
		totalIntensity += peaks[i].intensity;
		peaks[i].count = 1;
		peaks[i].maxPossible = 1;
	}

	return numPeaksRead;
}


int SpectraFile::readDtaPeakList(FILE* stream, const SingleSpectrumHeader* header, 
								Peak* peaks, const Config* config) const
{
	assert(header->getFileType() == IFT_DTA);
	

	mass_t		mass=-1.0;
	intensity_t intensity=-1.0;
	char buffer[256];

	if( ! fgets(buffer, 256, stream) )
			return 0;

	istringstream is(buffer);
	is >> mass >> intensity;

	peaks[0].mass = mass;
	peaks[0].intensity = intensity;

	assert(mass == header->getFirstPeakMass());

	int peakCount = 1;
	while (fgets(buffer, 256, stream))
	{
		if (buffer[0] < '0' || buffer[0] > '9')
				break;

		istringstream is(buffer);
				
		mass_t mass = -1.0;
		intensity_t intensity = -1.0;

		is >> mass >> intensity;
		// NP3 GOT change ms2 baseline intensity <= 0
		if (mass <=0.0 || intensity<=0.0)   // the peak probably got rounded down
			continue;

		// a fix to avoid problems with "0" intensity peaks 
		if (intensity == 0.0) 
			intensity = 1E-9;

		peaks[peakCount].mass		= mass;
		peaks[peakCount].intensity  = intensity;

		++peakCount;
	}

	return peakCount;
}

int SpectraFile::readMs2PeakList(FILE* stream, const SingleSpectrumHeader* header, 
								Peak* peaks, const Config* config) const
{
	assert(header->getFileType() == IFT_MS2);
	
	mass_t		mass=-1.0;
	intensity_t intensity=-1.0;
	char buffer[256];

	if( ! fgets(buffer, 256, stream) )
			return 0;

	istringstream is(buffer);
	is >> mass >> intensity;

	peaks[0].mass = mass;
	peaks[0].intensity = intensity;

	assert(mass == header->getFirstPeakMass());

	int peakCount = 1;
	while (fgets(buffer, 256, stream))
	{
		if (buffer[0] < '0' || buffer[0] > '9')
				break;

		istringstream is(buffer);
				
		mass_t mass = -1.0;
		intensity_t intensity = -1.0;

		is >> mass >> intensity;

		// NP3 GOT  ms2 baseline intensity <= 0
		if (mass <=0.0 || intensity<=0.0)   // the peak probably got rounded down
			continue;

		// a fix to avoid problems with "0" intensity peaks 
		if (intensity == 0.0) 
			intensity = 1E-9;

		peaks[peakCount].mass		= mass;
		peaks[peakCount].intensity  = intensity;

		++peakCount;
	}

	return peakCount;
}


int SpectraFile::readMgfPeakList(FILE* stream, const SingleSpectrumHeader* header, 
								 Peak* peaks, const Config* config) const
{
	assert (header->getFileType() == IFT_MGF);

	mass_t		mass=-1.0;
	intensity_t intensity=-1.0;
	char buffer[256];

	while (true)
	{

		if( ! fgets(buffer, 256, stream) )
			return 0;

		if (! strncmp("END IONS", buffer, 7))
			return 0;

		mass=-1.0;
		intensity=-1.0;
		istringstream is(buffer);
		is >> mass >> intensity;

		// NP3 GOT min mS2 intensity baseline
		if (mass>0.0 &&  intensity>0.0)
			break;
	}
	
	peaks[0].mass = mass;
	peaks[0].intensity = intensity;

	if (mass != header->getFirstPeakMass())
	{
		cout << "Mass:  " << mass << endl;
		cout << "Header:" << header->getFirstPeakMass() << endl;
		cout << "Title: " << header->getTitle() << endl;
	}
	assert(mass == header->getFirstPeakMass());

	int peakCount = 1;
	while (fgets(buffer, 256, stream))
	{
		if (! strncmp("END IONS",buffer,7))
			break;

		istringstream is(buffer);
				
		mass_t mass = -1.0;
		intensity_t intensity = -1.0;

		is >> mass >> intensity;

		if (mass <=0.0 || intensity<=0.0)   // the peak probably got rounded down
			continue;

		// a fix to avoid problems with "0" intensity peaks 
		if (intensity == 0.0) 
			intensity = 1E-9;

		peaks[peakCount].mass		= mass;
		peaks[peakCount].intensity  = intensity;
		++peakCount;
	}

	return peakCount;
}



int SpectraFile::readMzxmlPeakList(FILE* stream, 
								   const SingleSpectrumHeader* header, 
								   Peak* peaks, 
								   const Config* config) const
{
	static char*  buffer = 0;
    static char*  decodedPeakBuffer = 0;
	static float* floatBuffer = 0;
    static size_t maxNumPeaks = 0;
    //cout << "T3.72 peakCount " <<  header->getOriginalNumPeaks() << "maxPeaks " << maxNumPeaks <<  endl;
	const size_t peakCount = header->getOriginalNumPeaks();

	if (maxNumPeaks < peakCount)
	{
		maxNumPeaks = 200 + 2 * peakCount;
        // NP3 GOT ERROR when mzxml has more peaks than max - deleting not allocated pointers
		if (buffer)
			delete [] buffer;
        //cout << "T3.75" << endl;
		if (decodedPeakBuffer)
			delete [] decodedPeakBuffer;
        //cout << "T3.76" << endl;
		if (floatBuffer)
			delete [] floatBuffer;
        //cout << "T3.77" << endl;
		//buffer = new char[256 + peakCount * 12];
		//decodedPeakBuffer = new char[peakCount * 8];
		//floatBuffer = new float[peakCount *2];
	}
	buffer = new char[256 + maxNumPeaks * 12];
	decodedPeakBuffer = new char[maxNumPeaks * 8];
	floatBuffer = new float[maxNumPeaks *2];


	const size_t bytesToRead = 256 + peakCount* 12;
	const size_t bytesRead = fread(buffer,sizeof(char),bytesToRead,stream);
	buffer[bytesRead]='\0';

	if (strncmp(buffer,"<peaks",6))
		return 0;

	char* peakStr = buffer;
	int byteOrderLittle = 0;
	const char* byteOrderStr = strstr(peakStr, "byteOrder=\"");
	if (byteOrderStr)
	{
		byteOrderStr += 11;
		if (!strncmp(byteOrderStr, "network", 7))
		{
			byteOrderLittle = 0;
		}
		if (!strncmp(byteOrderStr, "big", 3))
		{
			byteOrderLittle = 0;
		}
		if (!strncmp(byteOrderStr, "little", 6))
		{
			byteOrderLittle = 1;
		}
	}
	peakStr = strstr(peakStr, ">");

	if (!peakStr)
		return 0;

	char *peakBuffer = peakStr+1;
	int trail = (peakCount % 3);
	if (!(peakCount % 3))
	{
		peakBuffer[peakCount * 32/3] = '\0';
	}
	else
		peakBuffer[(peakCount * 32/3) + trail + 1] = '\0';

	b64_decode_mio( decodedPeakBuffer, peakBuffer);
	const size_t maxIndex = 2 * peakCount;
	for (size_t floatIndex = 0; floatIndex < maxIndex; floatIndex++)
	{
	#ifdef BYTEORDER_LITTLE_ENDIAN
		if (!byteOrderLittle)
		{
			char byteSwap = decodedPeakBuffer[floatIndex*4];
			decodedPeakBuffer[floatIndex*4] = decodedPeakBuffer[floatIndex*4 + 3];
			decodedPeakBuffer[floatIndex*4 + 3] = byteSwap;
			byteSwap = decodedPeakBuffer[floatIndex*4 + 1];
			decodedPeakBuffer[floatIndex*4 + 1] = decodedPeakBuffer[floatIndex*4 + 2];
			decodedPeakBuffer[floatIndex*4 + 2] = byteSwap;
		}
	//	memcpy(floatBuffer + floatIndex, decodedPeakBuffer + floatIndex * 4, 4);
	#else
		if (byteOrderLittle)
		{
			char byteSwap = decodedPeakBuffer[floatIndex*4];
			decodedPeakBuffer[floatIndex*4] = decodedPeakBuffer[floatIndex*4 + 3];
			decodedPeakBuffer[floatIndex*4 + 3] = byteSwap;
			byteSwap = decodedPeakBuffer[floatIndex*4 + 1];
			decodedPeakBuffer[floatIndex*4 + 1] = decodedPeakBuffer[floatIndex*4 + 2];
			decodedPeakBuffer[floatIndex*4 + 2] = byteSwap;
		}
	//	memcpy(peakBuffer + floatIndex, decodedPeakBuffer + floatIndex * 4, 4);
	#endif
	}

	float *p = reinterpret_cast<float*>(decodedPeakBuffer);
	size_t peakIdx=0;
	for (size_t i=0; i<peakCount; i++)
	{
		peaks[peakIdx].mass = *p++;
		peaks[peakIdx].intensity = *p++;
		// NP3 GOT added mass check
		if (peaks[peakIdx].intensity> 0.0 && peaks[peakIdx].mass > 0.0)
			peakIdx++;
	}

	// if the precusor intensity is 0 (ms3?) add up the peak intensity
	if (header->getPrecursorIntensity()<=0.0)
	{
		intensity_t totalIntensity=0.0;
		for (size_t i=0; i<peakIdx; i++)
			totalIntensity += peaks[i].intensity;
		
		if (totalIntensity <= 0.0) // do not use this spectrum
			return 0;

		// replace the precursor intensity, becasue 0 intensity can cause problems later on
		SingleSpectrumHeader* varHeader = const_cast<SingleSpectrumHeader*>(header);
		varHeader->setPrecursorIntensity(totalIntensity);
	}

	return peakIdx; // this is the actual number of non-zero peaks read
}




int	SpectraFile::readDatPeakList(FILE* stream, 
								 const SingleSpectrumHeader* header, 
								 Peak* peaks, 
								 const Config* config) const
{
	if (fread(peaks, sizeof(Peak), header->getOriginalNumPeaks(), stream) != header->getOriginalNumPeaks())
		error("problem reading peak list in dat file!");

	return (header->getOriginalNumPeaks());
}



void SpectraFile::setAlldatasetIdxs(int idx)
{
	for (size_t i=0; i<headers_.size(); i++)
		headers_[i].setDatasetIndex(idx);
}

