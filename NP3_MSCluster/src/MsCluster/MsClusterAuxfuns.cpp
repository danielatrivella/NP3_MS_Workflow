#include "MsClusterAuxfuns.h"
#include "MsClusterIncludes.h"
#include "MsParameterStruct.h"
#include "MsModClusterWriter.h"
#include "DatFileWriter.h"
#include "../PepNovo/AllScoreModels.h"
#include "../PepNovo/DatFile.h"
#include "../PepNovo/ScanList.h"
#include "../Common/BufferedOutputFile.h"



/***********************************************************************
Assumes each dat file ends with the pattern _[MZ]_[\d+].dat .
This function splits the files into groups according to the MZ value
************************************************************************/
void splitDatListIntoSameMzGroups(const vector<string>& datList, 
								  vector< vector<string> >& datGroups)
{
	vector<string> mzLabels;

	mzLabels.clear();
	datGroups.clear();

	for (size_t i=0; i<datList.size(); i++)
	{
		const string& path = datList[i];
		size_t endPos = path.find_last_of('_');
		if (endPos == string::npos)
			continue;

		size_t startPos = path.find_last_of('_',endPos-1);
		if (startPos == string::npos)
			continue;

		string mzLabel = path.substr(startPos+1, endPos-startPos);

		size_t j;
		for (j=0; j<mzLabels.size(); j++)
			if (mzLabel == mzLabels[j])
				break;

		if (j<mzLabels.size())
		{
			datGroups[j].push_back(path);
		}
		else
		{
			datGroups.push_back(vector<string>(1,path));
			mzLabels.push_back(mzLabel);
		}
	}
}


struct FileOrder {
	bool operator< (const FileOrder& rhs) const
	{
		return (mainMz < rhs.mainMz ||
				mainMz == rhs.mainMz && secondaryIdx < rhs.secondaryIdx);
	}

	double mainMz;       // first number
	int	   secondaryIdx; // second number
	size_t fileIdx; // position in vector
};

void sortDatPathsAccordingToMz(vector<string>& datPaths)
{
	vector<FileOrder> fileOrder(datPaths.size());
	for (size_t i=0; i<datPaths.size(); i++)
	{
		const string& path = datPaths[i];
		size_t endPos = path.find_last_of('_');
		if (endPos == string::npos)
			error("coudln't parse dat file name:",path.c_str());
		size_t startPos = path.find_last_of('_',endPos-1);
		if (startPos == string::npos)
			error("coudln't parse dat file name:",path.c_str());

		string idxLabel = path.substr(endPos+1,path.length() - 4 - endPos);
		string mzLabel = path.substr(startPos+1, endPos-startPos);
		fileOrder[i].mainMz = atof(mzLabel.c_str());
		fileOrder[i].secondaryIdx = atoi(idxLabel.c_str());
		fileOrder[i].fileIdx = i;
	}
	sort(fileOrder.begin(),fileOrder.end());
	
	vector<string> newPaths;
	newPaths.resize(datPaths.size());
	for (size_t i=0; i<fileOrder.size(); i++)
		newPaths[i]=datPaths[fileOrder[i].fileIdx];

	datPaths = newPaths;
}


void convertSpectra(const MsParameterStruct* params, AllScoreModels* model)
{
	const Config* config = model->get_config();
	PMCSQS_Scorer* pmcsqsModel = const_cast<PMCSQS_Scorer*>(model->get_pmcsqs_ptr());
	if (params->sqsThreshold>0.0 && ! pmcsqsModel->getIndInitializedSqs())
		error("Sqs model not initialized!, need a valid sqs model if using a filtering threshold!");

	if (params->sqsThreshold>0.0)
		cout << "Filtering spectra with SQS threshold of " << params->sqsThreshold << endl;

	assert(params->fileConversionType==IFT_DAT || params->fileConversionType == IFT_MGF);
	string newType = "dat";
	if (params->fileConversionType == IFT_MGF)
		newType = "mgf";

	size_t fileStartIdx = 0;
	vector<string> paths;
	if (params->inputFile.length()>0)
	{
		paths.push_back(params->inputFile);
	}
	else if (params->list.length()>0)
	{
		fileStartIdx = readListOfPaths(params->list.c_str(), paths);
	}
	else if (params->datList.length()>0)
	{
		fileStartIdx = readListOfPaths(params->datList.c_str(), paths);
	}
	else
		error("No input files given for conversion (must use --file, --list, or --dat-list)");

	if (params->verboseLevel>0)
	{
		cout << "Read " << paths.size() << " paths to spectra files." << endl;
		cout << "Converting data to " << newType << endl;
	}


	int numFilesWithoutSpectra = 0;
	size_t peakBufferSize = 10000;
	Peak*  peakBuffer = new Peak[peakBufferSize];

	char*  mgfBuffer = NULL;
	if (params->fileConversionType == IFT_MGF)
		mgfBuffer = new char[131072]; // should be enough for any MGF

	createDirIfDoesNotExist(params->outDir.c_str());

	ScanListManager inclusionScans, exclusionScans;
	if (params->inclusionList.length()>0)
		inclusionScans.initialize(params->inclusionList.c_str());
	if (params->exclusionList.length()>0)
		exclusionScans.initialize(params->exclusionList.c_str());

	// loop on files
	for (size_t i=0; i<paths.size(); i++)
	{
		DatFile			   datOut;
		BufferedOutputFile mgfOut;

		const double fileStartTime = time(NULL);

		if (params->verboseLevel>0)
			cout << i << "\tExtracting from: " << paths[i] << endl;

		// read headers
		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(paths[i].c_str(), config);
		SpectraList sl(sa);
		sl.selectAllAggregatorHeaders();
		
		if (params->verboseLevel>0)
			cout << "\tFound " << sl.getNumHeaders() << " spectra...";

		if (sl.getNumHeaders() == 0)
		{
			numFilesWithoutSpectra++;
			continue;
		}

		// open converted file
		string fname;
		if (params->outputName.length()>0)
		{
			ostringstream oss;
			oss << params->outputName;
			if (paths.size()>0)
				oss << "_" << i;
			fname = oss.str();
		}
		else
		{
			getFileNameWithoutExtension(paths[i].c_str(), fname);
		}
		string outputName = params->outDir + "/" + fname + "." + newType;
		if (params->fileConversionType == IFT_MGF)
		{
			mgfOut.open(outputName.c_str(), 1<<20);
		}
		else if (params->fileConversionType == IFT_DAT)
		{
			datOut.openForWriting(outputName.c_str(), params->verboseLevel);
		}

		// read and write spectra
		size_t numExtracted =0;
		for (size_t j=0; j<sl.getNumHeaders(); j++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(j);
			SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
			nonConstHeader->setSpectraFileIndexInList(0);

			if (header->getOriginalNumPeaks()>1e6)
				continue;

			if (header->getOriginalNumPeaks()>= peakBufferSize)
			{
				delete [] peakBuffer;
				peakBufferSize = header->getOriginalNumPeaks()*2;
				peakBuffer = new Peak[peakBufferSize];
			}

			PeakList pl;
			pl.setPeaksPtr( peakBuffer );
			pl.readPeaksToBuffer(sa, header, peakBuffer);
			pl.initializePeakList(config, true);

			// NP3 GOT change few peaks from 15 to 1 @@@
			if (pl.getNumPeaks()<1) // don't bother with spectra with too few peaks
				continue;

			nonConstHeader->setSpectraFileIndexInList(fileStartIdx + i);

			if (pmcsqsModel && params->sqsThreshold>0.0)
			{
				size_t maxCharge=0;
				const float sqs = pmcsqsModel->calculateSqsScore(config, pl, &maxCharge);
				if (sqs<params->sqsThreshold || maxCharge == 0)
					continue;
				header->setSqs(sqs);
			}

			if (params->fileConversionType == IFT_DAT)
			{
				const size_t orgNumPeaks = pl.getHeader()->getOriginalNumPeaks();
				SingleSpectrumHeader* header = const_cast<SingleSpectrumHeader*>(pl.getHeader());

				if (params->datasetIdx >= 0)
					header->setDatasetIndex(params->datasetIdx);
				
				assert(pl.getNumPeaks()>0);
				assert( header->getFileType() != IFT_MZXML || header->getScanNumber() >= 0);
				header->setOriginalNumPeaks(pl.getNumPeaks());
				datOut.writePeakList(pl);
			}
			else if (params->fileConversionType == IFT_MGF)
			{
				size_t n=pl.outputToMgfBuffer(mgfBuffer, header);
				if (n>=131072)
					error("MGF buffer in cluster out too small, increase size!");

				mgfOut.writeToBuffer(mgfBuffer, n);
			}
			numExtracted++;
		}

		if (params->verboseLevel>0)
		{
			const double fileEndTime = time(NULL);
			cout << " Wrote " << numExtracted << " to " << newType << " files (this took " 
				 << fileEndTime-fileStartTime << " sec.)" << endl << endl;
			cout.flush();
		}

		// close files
		if (params->fileConversionType == IFT_MGF)
		{
			mgfOut.close();
		}
		else if (params->fileConversionType == IFT_DAT)
		{
			datOut.closeAfterWriting();
		}
	}

	// free mem
	if (peakBuffer)
		delete [] peakBuffer;
	if (mgfBuffer)
		delete [] mgfBuffer;
}


void convertArchive(const MsParameterStruct* params, AllScoreModels* model)
{
	const Config* config = model->get_config();
	if (params->datList.length() == 0)
		error("Must supply --dat-list when converting archive");
	if (params->outputName.length() == 0)
		error("Must supply --ouput-name when converting archive");

	string newType = "mgf";
	size_t fileStartIdx = 0;
	vector<string> paths;
	fileStartIdx = readListOfPaths(params->datList.c_str(), paths);

	ModClusterWriter modClusterWriter;
	string path = params->outDir + "/" + params->outputName;
	modClusterWriter.init(path, params->outputFileSize);

	if (params->verboseLevel>0)
	{
		cout << "Read " << paths.size() << " paths to spectra files." << endl;
		cout << "Converting data to " << newType << endl;
	}


	int numFilesWithoutSpectra = 0;
	size_t peakBufferSize = 10000;
	Peak*  peakBuffer = new Peak[peakBufferSize];

	char*  mgfBuffer = new char[131072]; // should be enough for any MGF

	createDirIfDoesNotExist(params->outDir.c_str());

	// loop on files
	for (size_t i=0; i<paths.size(); i++)
	{
		DatFile			   datOut;
		BufferedOutputFile mgfOut;

		const double fileStartTime = time(NULL);

		if (params->verboseLevel>0)
			cout << i << "\tExtracting from: " << paths[i] << endl;

		// read headers
		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(paths[i].c_str(), config);
		SpectraList sl(sa);
		sl.selectAllAggregatorHeaders();
		
		if (params->verboseLevel>0)
			cout << "\tFound " << sl.getNumHeaders() << " spectra...";

		if (sl.getNumHeaders() == 0)
		{
			numFilesWithoutSpectra++;
			continue;
		}

		// read and write spectra
		size_t numExtracted =0;
		for (size_t j=0; j<sl.getNumHeaders(); j++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(j);
			SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
			nonConstHeader->setSpectraFileIndexInList(0);

			if (header->getOriginalNumPeaks()>1e6)
				continue;

			if (header->getOriginalNumPeaks()>= peakBufferSize)
			{
				delete [] peakBuffer;
				peakBufferSize = header->getOriginalNumPeaks()*2;
				peakBuffer = new Peak[peakBufferSize];
			}

			PeakList pl;
			pl.setPeaksPtr( peakBuffer );
			pl.readPeaksToBuffer(sa, header, peakBuffer);
			pl.initializePeakList(config, true);

			// NP3 GOT change few peaks from 15 to 1 @@@
			if (pl.getNumPeaks()<1) // don't bother with spectra with too few peaks
				continue;

			nonConstHeader->setSpectraFileIndexInList(fileStartIdx + i);


			size_t n=pl.outputToMgfBuffer(mgfBuffer, header);
			if (n>=131072)
				error("MGF buffer in cluster out too small, increase size!");

			modClusterWriter.writeToBuffer(mgfBuffer, n);
			numExtracted++;
		}

		if (params->verboseLevel>0)
		{
			const double fileEndTime = time(NULL);
			cout << " Wrote " << numExtracted << " to " << newType << " files (this took " 
				 << fileEndTime-fileStartTime << " sec.)" << endl << endl;
			cout.flush();
		}
	}

	modClusterWriter.close();

	// free mem
	if (peakBuffer)
		delete [] peakBuffer;
	if (mgfBuffer)
		delete [] mgfBuffer;
}


void readIdsTitleFromIdFile(const MsParameterStruct* params, map<string,int>& idTitles)
{
	// read ids
	vector<string> idPaths;
	readListOfPaths(params->idListPath.c_str(), idPaths);
	string titleLine;
	char buffer[256];
	size_t numLines = 0;
	for (size_t i=0; i<idPaths.size(); i++)
	{
		ifstream ifs(idPaths[i].c_str());
		if (! ifs.good())
		{
			cout << "Warning: could not open id file: " << idPaths[i] << endl;
			continue;
		}
		ifs.getline(buffer,256);
		string firstToken = getFirstToken(buffer,ifs.gcount());
		if (firstToken.length()>0)
		{
			if (firstToken == "CLUSTER_ID" || firstToken == "CLUSTER" || firstToken[0] == '#')
				titleLine = std::string(buffer + firstToken.length());
		}

		while (ifs.good())
		{
			ifs.getline(buffer,256);
			string firstToken = getFirstToken(buffer,ifs.gcount());
			idTitles[firstToken]++;
		}
		ifs.close();
	}
	cout << "Found ids for " << idTitles.size() << " spectra" << endl;
}






