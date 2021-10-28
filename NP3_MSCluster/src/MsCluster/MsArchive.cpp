#include "MsArchive.h"
#include "ClusterOutputter.h"
#include "MsParameterStruct.h"
#include "MsSimilarityModel.h"
#include "DatFileManager.h"
#include "DatFileWriter.h"
#include "MsClusterAuxfuns.h" 
#include "../Common/auxfun.h"
#include "../PepNovo/Config.h"
#include "../PepNovo/SpectraList.h"
#include "../PepNovo/PMCSQS.h"
#include "../PepNovo/AllScoreModels.h"


bool DataSetRecord::readRecordFromBuffer(const char* buffer)
{
	istringstream iss(buffer);
	iss >> index;
	iss >> dateAdded >> numClusters >> numSingletons >> fragTolerance >> precursorTolerance >>  precursorPPM >> 
		fileList;
	
	ostringstream oss;
	string s;
	while (iss >> s)
	{
		oss << s << " ";
	}

	comment = oss.str();
//	cout << "buffer size: " << strlen(buffer) << endl;
//	cout << "gcount:      " << iss.tellg() << endl;
	return (index>=0 && dateAdded>0 && numClusters>0 && numSingletons>0);
}


void DataSetRecord::writeRecordToBuffer(ostream& oss, size_t newIndex) const
{
	oss << newIndex << "\t" << dateAdded << "\t" << numClusters << "\t" << numSingletons << "\t";
	oss << setprecision(3) << fragTolerance <<  "\t" << fixed << precursorTolerance <<  "\t" << precursorPPM << "\t";
	oss << fileList << "\t" << comment << endl;
}



void MsArchive::clear()
{
	archiveName_ = "Archive";
	archivePath_ = "."; 

	archiveDatasets_.clear();
	datNames_.clear();
}





bool MsArchive::writeArchive(const ClusterOutputter& clusterOutputter, 
							 const MsParameterStruct* params,
							 const Config* config)
{
	const vector<string>& allClustPaths = clusterOutputter.getAllClustPaths();

	// open main archive file
	if (! params->gotMergeArchives && ! params->gotCreateArchiveFromMgfs)
	{
		string listFileName = params->outputNameWithVersion + "_ds_0.txt";

		// write file list to file list area
	#if defined(WIN32) || defined(WIN64)
		string spectrumListsDir = params->archiveOutDir + "\\spectrum_lists";
		string listPath = spectrumListsDir + "\\" + listFileName;
	#else
		string spectrumListsDir = params->archiveOutDir + "/spectrum_lists";
		string listPath = spectrumListsDir + "/" + listFileName;
	#endif
		
		createDirIfDoesNotExist(spectrumListsDir.c_str());

		copyFile(params->list.c_str(), listPath.c_str());

		DataSetRecord ds;
		ds.fileList            = listFileName;
		ds.index               = 0;
		ds.fragTolerance       = config->getTolerance();
		ds.precursorTolerance  = config->get_pm_tolerance();
		ds.dateAdded		   = getDateAsInt();
		ds.numClusters		   = clusterOutputter.getNumClustersWritten();
		ds.numSingletons       = clusterOutputter.getTotalSpectraCount();

		archiveDatasets_.push_back(ds);
	}
	else // update the list file names
	{
		for (size_t i=0; i<archiveDatasets_.size(); i++)
		{
			ostringstream oss;
			oss << params->outputNameWithVersion + "_ds_" << i << ".txt";
			archiveDatasets_[i].fileList = oss.str();
		}
	}
	

	ostringstream oss;
	oss << params->outDir <<  "/" << params->outputName << "_" << archiveGeneration_ << "_"
		<< params->batchIdx << "_archive.txt";
	ofstream ofs(oss.str().c_str());
	if (! ofs.good())
		error("Couldn't open archive file for writing: ",oss.str().c_str());
	ofs << "Name          " << params->outputName << endl;
	ofs << "Generation    " << archiveGeneration_ << endl;
	ofs << "First created " << (firstCreated_.length() == 0 ? getTimeString() : firstCreated_) << endl;
	ofs << "Last updated  " << getTimeString() << endl;

	if (params->gotMergeArchives || params->gotCreateArchiveFromMgfs)
	{
		longInt8_t numSingletons=0;
		for (size_t i=0; i<archiveDatasets_.size(); i++)
			numSingletons += archiveDatasets_[i].numSingletons;
		ofs << "Number of Spectra  " << numSingletons << endl;
	}
	else
	{
		ofs << "Number of Spectra  " << clusterOutputter.getTotalSpectraCount() << endl;
		
	}
	if (! params->gotCreateArchiveFromMgfs)
	{
		ofs << "Number of Clsuters " << clusterOutputter.getNumClustersWritten() << endl;
		ofs << "Number of Files    " << allClustPaths.size() << endl;
	}
	
	ofs << endl << "Datasets: " << archiveDatasets_.size() << endl;
	for (size_t i=0; i<archiveDatasets_.size(); i++)
		archiveDatasets_[i].writeRecordToBuffer(ofs,i);
	

	ofs << endl << "List of archive file names:" << endl;

	// write files (clust or dat, without suffix?)
	if (! params->gotCreateArchiveFromMgfs)
	{
		for (size_t i=0; i<allClustPaths.size(); i++)
		{
			string fileName;
			getFileNameWithoutExtension(allClustPaths[i].c_str(), fileName);
			ofs << fileName << endl;
		}
	}
	else
	{
		for (size_t i=0; i<datNames_.size(); i++)
			ofs << datNames_[i] << endl;
	}
	
	ofs.close();



	return true;
}


bool MsArchive::readArchive(const std::string& mainArchiveFilePath)
{
	clear();
	ifstream ifs(mainArchiveFilePath.c_str());
	if (! ifs.good())
	{
		cout << "Warning: Couldn't open archive file: " << mainArchiveFilePath.c_str() << endl;
		return false;
	}
	
	archiveName_ = "no name";
	archiveDatasets_.clear();
	datNames_.clear();

	bool ok=false;
	char line[256];
	char buff[256];
	while (ifs.getline(line,256))
	{
		if (! strcmp(line,"List of archive file names:"))
		{
			ok=true;
			break;
		}

		if (! strncmp(line,"Name",4) )
		{
			sscanf(line+4,"%s",buff);
			archiveName_ = std::string(buff);
			continue;
		}
		if (! strncmp(line,"Generation",10) )
		{
			sscanf(line+11,"%d",&archiveGeneration_);
			continue;
		}
		if (! strncmp(line,"First created",14) )
		{
			firstCreated_ = std::string(line+15);
			continue;
		}

		if (! strncmp(line,"Datasets",8) )
		{
			int numDatasets=0;
			sscanf(line+9,"%d",&numDatasets);
			if (numDatasets>0)
			{
				archiveDatasets_.resize(numDatasets);
				// read dataset records
				for (int i=0; i<numDatasets; i++)
				{
					ifs.getline(line,256);
					if (! ifs.good() || ! archiveDatasets_[i].readRecordFromBuffer(line))
						error("Bad line in archive file: ",line);
				}
			}
			else
				error("Archive file contains 0 datasets!");

			continue;
		}
	}

	if (! ok || archiveName_ == "no name" || archiveDatasets_.size() == 0)
	{
		cout << "Warning: Couldn't read archive correctly: " << mainArchiveFilePath.c_str() << endl;
		return false;
	}

	char nameBuff[256];
	while (ifs.getline(line,256))
	{
		sscanf(line,"%s",nameBuff);
		if (strlen(nameBuff)>2)
			datNames_.push_back(std::string(nameBuff));
	}

	cout << "Archive read " << datNames_.size() << " file names from archive " << mainArchiveFilePath << endl;

	string dirPath;
	getDirFromFullPath(mainArchiveFilePath.c_str(), dirPath);
	ostringstream oss;
	
	const string suffix = "_archive.txt";
	archivePath_ = mainArchiveFilePath.substr(0, mainArchiveFilePath.length() - suffix.length());
	cout << "Archive files are assumed to be in: " << archivePath_ << endl;

	return true;
}


void MsArchive::addArchive(const MsParameterStruct* params, MsArchive& other)
{
#if defined(WIN32) || defined(WIN64)
	const string thisListDir = params->archiveOutDir + "\\spectrum_lists";
#else
	const string thisListDir = params->archiveOutDir + "/spectrum_lists";
#endif

	createDirIfDoesNotExist(thisListDir.c_str());

	if (params->archiveOutDir != params->archivePath)
	{
	// copy archive file list files
		for (size_t i=0; i<archiveDatasets_.size(); i++)
		{
		#if defined(WIN32) || defined(WIN64)
			string listSource = archivePath_ + "\\spectrum_lists\\" + archiveDatasets_[i].fileList;
			ostringstream oss;
			oss << thisListDir << "\\" << params->outputNameWithVersion << "_ds_" << i << ".txt";
		#else
			string listSource = archivePath_ + "/spectrum_lists/" + archiveDatasets_[i].fileList;
			ostringstream oss;
			oss << thisListDir << "/" << params->outputNameWithVersion << "_ds_" << i << ".txt";
		#endif
			string listTarget =  oss.str();
			copyFile(listSource.c_str(), listTarget.c_str());
		}
	}

	const size_t orgNumDatasets = archiveDatasets_.size();
	for (size_t i=0; i<other.archiveDatasets_.size(); i++)
	{
		archiveDatasets_.push_back(other.archiveDatasets_[i]);
		totalSingletons_ += other.archiveDatasets_[i].numSingletons;
		totalClusters_	 += other.archiveDatasets_[i].numClusters;
		
		ostringstream oss;
	#if defined(WIN32) || defined(WIN64)
		string listSource = other.archivePath_ + "\\spectrum_lists\\" + other.archiveDatasets_[i].fileList;
		oss << thisListDir << "\\" << params->outputNameWithVersion << "_ds_" << orgNumDatasets + i << ".txt";
	#else
		string listSource = other.archivePath_ + "/spectrum_lists/" + other.archiveDatasets_[i].fileList;
		oss << thisListDir << "/" << params->outputNameWithVersion << "_ds_" << orgNumDatasets + i << ".txt";
	#endif
		string listTarget =  oss.str();

		// copy other's file lists to this archives file lists
		copyFile(listSource.c_str(), listTarget.c_str());
	}

	
}


void MsArchive::generateArchiveDatPaths(vector<string>& paths) const
{
	paths.resize(datNames_.size());

#if defined(WIN32) || defined(WIN64)
	for (size_t i=0; i<datNames_.size(); i++)
		paths[i]= archivePath_ + "\\dat\\" + datNames_[i] + ".dat";
#else
	for (size_t i=0; i<datNames_.size(); i++)
		paths[i]= archivePath_ + "/dat/" + datNames_[i] + ".dat";
#endif
}	





void MsArchive::printStatistics(const Config* config)
{
	size_t numBadClusters=0;
	size_t numClusters=0, numAnnClusters=0;
	longInt8_t numSpectra=0, numAnnSpectra=0;
	map<string, size_t> peptideCounts, peptideChargeCounts;

	vector<string> archiveDatPaths;
	generateArchiveDatPaths(archiveDatPaths);

	for (size_t i=0; i<archiveDatPaths.size(); i++)
	{
		DatFile datFile;
		datFile.openForReading(archiveDatPaths[i].c_str());

		size_t numBytes = 0;
		float  sqs;
		char* specStart = 0;
		long peaksFilePosition = 0;
		char* peaksBufferPosition = 0;
		size_t numSpectraRead = 0;
		mass_t previousMz = MIN_FLOAT;
		while( datFile.getNextDatSpectrum(specStart, numBytes, &sqs, &peaksFilePosition, &peaksBufferPosition) )
		{
			SingleSpectrumHeader ssh;
			ssh.setFileType(IFT_DAT);
			ssh.scanSpectrumHeaderFromBuffer(specStart, config);
			int clusterSize = ssh.getClusterSize();
			if (clusterSize<=0 || clusterSize>10000000)
			{
				numBadClusters++;
				continue;
			}

			numClusters++;
			numSpectra+= clusterSize;
			if (ssh.getPeptideStr().length()>0)
			{
				peptideCounts[ssh.getPeptideStr()]++;
				ostringstream oss;
				oss << ssh.getPeptideStr() << ssh.getCharge();
				peptideChargeCounts[oss.str()]++;
				numAnnClusters++;
				numAnnSpectra+=clusterSize;
			}
		}
		datFile.closeAfterReading();
	}

	if (numBadClusters>0)
		cout << "Encounterd " << numBadClusters << "!!!" << endl;

	cout << "Read " << archiveDatPaths.size() << " dat files with " << numAnnClusters << "/" 
			 << numClusters << " annotated clusters." << endl;
	cout << "These clusters correspond to " << numAnnSpectra << "/" << numSpectra << " annotated spectra." << endl;
	cout << "Unique peptide/charge: " << peptideChargeCounts.size() << endl;
	cout << "Unique peptide       : " << peptideCounts.size() << endl;
}



/*! @brief reads spectra into clusters

The clusters are sorted according to the m/z.
*/
void readQuerySpectraIntoClusters(const MsParameterStruct* params, vector<Cluster>& clusters, vector<SingleSpectrumHeader>& headers, 
								  const AllScoreModels* model)
{
	const string& orgList  = (params->spectraListToLoad.length()>0 ? params->spectraListToLoad : params->list);
	float sqsThreshold     = params->sqsThreshold;
	size_t  fileStartIdx   = params->startFileIdx;
	int     verboseLevel   = params->verboseLevel;

	const Config* config = model->get_config();
	PMCSQS_Scorer* pmcsqsModel = const_cast<PMCSQS_Scorer*>(model->get_pmcsqs_ptr());
	if (sqsThreshold>0.0 && ! ( pmcsqsModel && pmcsqsModel->getIndInitializedSqs()) )
		error("Sqs model not initialized!, need a valid sqs model if using a filtering threshold!");

	if (sqsThreshold>0.0)
		cout << "Filtering spectra with SQS threshold of " << sqsThreshold << endl;

	
	vector<string> paths;
	size_t startIdx = readListOfPaths(orgList.c_str(), paths);

	size_t peakBufferSize = 100000;
	Peak*  peakBuffer = new Peak[peakBufferSize];

	SpectraAggregator sa;
	sa.initializeFromListOfPaths(paths, config, 0, startIdx);
	SpectraList sl(sa);
	sl.selectHeaders(params->minMz, params->maxMz);

	if (params->verboseLevel>0)
		cout << "\tFound " << sl.getNumHeaders() << " spectra in " << paths.size() << " files." << endl;

	clusters.clear();
	headers.clear();

	if (sl.getNumHeaders() == 0)
		return;

	// allocate space
	try
	{
		clusters.resize(sl.getNumHeaders());
		headers.resize(sl.getNumHeaders());
	}
	catch (...)
	{
		cout << "Problem allocating memory for " << sl.getNumHeaders() << " query spectra" << endl;
		cout << "Try breakin input into batches (e.g., using --miz-mz and --max-mz)" << endl;
		exit(0);
	}

	size_t numExtracted =0;
	for (size_t i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
		if (header->getOriginalNumPeaks()>100000)
			continue;
		if (header->getOriginalNumPeaks()>= peakBufferSize)
		{
			delete [] peakBuffer;
			peakBufferSize = header->getOriginalNumPeaks()*2;
			peakBuffer = new Peak[peakBufferSize];
		}

		PeakList pl;
		pl.setPeaksPtr( peakBuffer );

        // NP3 GOT change few peaks from 7 to 1 @@@
		if (pl.readPeaksToBuffer(sa, header, peakBuffer) < 1) // if not enough peaks read, skip this spectrum
			continue;

		pl.initializePeakList(config, true);
		// NP3 GOT change few peaks from 7 to 1 @@@
		if (pl.getNumPeaks()<1) // don't bother with spectra with too few peaks
			continue;

		mass_t prevMass = 0.0;
		bool goodPeaks = true;
		for (size_t i=1; i<pl.getNumPeaks(); i++)
		{
			mass_t mass = pl.getPeakMass(i);
			if (prevMass >= mass)
			{
				goodPeaks = false;
				break;
			}
			prevMass = mass;
		}
		if (! goodPeaks)
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
			const float sqs = pmcsqsModel->calculateSqsScore(config, pl, &maxCharge);
			if (sqs<sqsThreshold || maxCharge == 0)
				continue;
			header->setSqs(sqs);

			if (params->gotCorrectPM)
			{
				PmcSqsChargeRes res;
				pmcsqsModel->computeBestMzValuesForCharge(pl, maxCharge, config->get_pm_tolerance(), res);
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

		headers[numExtracted] = *header;
		if ( clusters[numExtracted].createNewCluster(numExtracted, &headers[numExtracted], pl.getPeaks(), pl.getNumPeaks()) )
			numExtracted++;
	}

	if (params->verboseLevel)
		cout << "Including " << numExtracted << " query spectra in archive search" << endl;

	clusters.resize(numExtracted);
	headers.resize(numExtracted);
}




size_t extractClusterIdxFromTitle(const string& title)
{
	const size_t pos = title.find_last_of('.');
	if (pos == string::npos)
		return MAX_SIZE_T;
	stringstream ss(title.substr(pos+1));
	size_t retVal;
	if ((ss >> retVal).fail())
		error("Could not convert to size_t: ", title.c_str(), title.substr(pos).c_str());
	return retVal;
}

struct ClusterComposition {
	ClusterComposition() : size(-1), maxDsIdx(-1), compositionStr(std::string()) {}
	int size;
	int maxDsIdx;
	string compositionStr;
};

void addClustInfoToResults(const map<string,int>& titlesToCheck, const vector<string>& clustPaths, const string& oldRes, const string& newRes)
{
	// vector that will indicate if the cluster idx number (i.e., 1772379 in ARCHPNNL_0_7.1772379) is one of a cluster that needs to be checked
	// quering the bool will make this run faster (assuming most clusters are not really reported on)
	vector<bool> clusterIdxsToCheck;

	// find max idx
	unsigned int maxIdx = 0;
	for (map<string,int>::const_iterator it=titlesToCheck.begin(); it != titlesToCheck.end(); it++)
	{
		unsigned int idx = extractClusterIdxFromTitle(it->first);
		if (idx>maxIdx)
			maxIdx = idx;
	}
	clusterIdxsToCheck.resize(maxIdx+1,0);
	for (map<string,int>::const_iterator it=titlesToCheck.begin(); it != titlesToCheck.end(); it++)
	{
		size_t idx = extractClusterIdxFromTitle(it->first);
		assert(idx != MAX_SIZE_T);
		clusterIdxsToCheck[idx]=true;
	}
	
	// create indicator vector for all idxs that are worth checking
	map<string, ClusterComposition> clustComps;
	for (size_t i=0; i<clustPaths.size(); i++)
	{
		ifstream ifs(clustPaths[i].c_str());
		if (! ifs.good())
			error("Could not open clust file for reading: ",clustPaths[i].c_str());
		char buffer[256];
		while( ifs.getline(buffer,256))
		{
			const size_t count = ifs.gcount();
			if (count<4)
				continue;
		
			const string firstToken = getFirstToken(buffer,count);
			const size_t idx = extractClusterIdxFromTitle(firstToken);
			if (idx>=clusterIdxsToCheck.size())
				continue;

			stringstream ss(buffer + firstToken.length() + 1);
			int numLines=0;
			if ( (ss>>numLines).fail())
				error("Could not parse line in clust file: ", clustPaths[i].c_str(), buffer);
			assert(numLines>0);
		
			bool needToBuildComposition = (clusterIdxsToCheck[idx] && titlesToCheck.find(firstToken) != titlesToCheck.end());
			if (needToBuildComposition)
			{
				ClusterComposition comp;
				comp.size = numLines;
				map<int,int> counts;
				for (int j=0; j<numLines; j++)
				{
					ifs.getline(buffer,256);
					int orgNum=-1;
					orgNum=atoi(buffer);
					if (orgNum<0 || orgNum>1000)
					{
						cout << "Buffer=" << buffer << endl;
						cout << "atoi(buffer)= " << atoi(buffer) << endl;
						exit(0);
					}
					counts[orgNum]++;
				}
				int maxDsIdx=-1;
				int maxDsCount = 0;

				// build composition string and find max org (the dataset with largest count)
				stringstream oss;
				bool firstPair=true;
				for (map<int,int>::const_iterator it=counts.begin(); it!=counts.end(); it++)
				{
					if (firstPair)
					{
						firstPair = false;
					}
					else
						oss << ",";
					oss << it->first << ":" << it->second;
					if (it->second > maxDsCount)
					{
						maxDsCount = it->second;
						maxDsIdx   = it->first;
					}
				}
				comp.maxDsIdx         = maxDsIdx;
				comp.compositionStr   = oss.str();
				clustComps[firstToken]= comp;
			}
			else // skip lines
			{
				for (int j=0; j<numLines; j++)
					ifs.getline(buffer,256);
			}
		}
	}

	if (titlesToCheck.size() != clustComps.size())
		cout << "WARNING: not all clust data was found " << clustComps.size() << "/" << titlesToCheck.size() << endl;

	// add composition information to results lines
	ifstream orgStream(oldRes.c_str());
	ofstream newStream(newRes.c_str());

	char buffer[256];
	orgStream.getline(buffer,256);
	newStream << buffer << endl;
	orgStream.getline(buffer,256);
	newStream << buffer << endl;
	

	while (orgStream)
	{
		orgStream.getline(buffer,256);
		if (orgStream.gcount()<3)
		{
			newStream << endl;
			continue;
		}
		const string firstToken = getFirstToken(buffer,orgStream.gcount());
		newStream << buffer;
		if (firstToken.find_last_of('.') == string::npos)
		{
			newStream << endl;
			continue;
		}
		const map<string,ClusterComposition>::const_iterator it = clustComps.find(firstToken);
		if (it == clustComps.end())
		{
			error("Could not find cluster composition for: ",buffer);
		}
		else
			newStream << "\t" << it->second.size << "\t" << it->second.maxDsIdx << "\t" << it->second.compositionStr << endl;
	}

	orgStream.close();
	newStream.close();

}


void addIdsToResults(const map<string,int>& titlesToId, const string& idListPath, const string& orgPath, string& newPath)
{
	// collect info from id files
	vector<string> idPaths;
	readListOfPaths(idListPath.c_str(), idPaths);
	map<string, string> idLines;
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
			if (titlesToId.find(firstToken) != titlesToId.end())
				idLines[firstToken] = std::string(buffer + firstToken.length() + 1);
			numLines++;
		}
		ifs.close();
	}
	cout << "Found ids for " << idLines.size() << "/" << titlesToId.size() << " unique cluster ids (id files contained total of " << numLines << " id lines)" << endl;;

	ifstream orgStream(orgPath.c_str());
	ofstream newStream(newPath.c_str());

	orgStream.getline(buffer,256);
	newStream << buffer << endl;
	orgStream.getline(buffer,256);
	newStream << buffer << "\t" << titleLine << endl;


	while (orgStream)
	{
		orgStream.getline(buffer,256);
		string firstToken = getFirstToken(buffer,orgStream.gcount());
		newStream << buffer;
		const map<string,string>::const_iterator it = idLines.find(firstToken);
		if (it != idLines.end())
			newStream << "\t" << it->second;
		newStream << endl;
	}

	orgStream.close();
	newStream.close();
}


struct SimPair {
	SimPair(float s, clusterIdx_t i) : sim(s), idx(i) {}
	bool operator< (const SimPair& rhs) const
	{
		return (sim>rhs.sim);
	}
	float sim;
	clusterIdx_t idx;
};




void MsArchive::searchArchive(const MsParameterStruct* params, const AllScoreModels* model) const
{
	const Config* config = model->get_config();
	Cluster::setTolerances(config->getTolerance());
	MsSimilarityModel simModel;
	simModel.readSimilarityModel(config);

	cout << "Start of archive search " << getTimeString() << endl;

	// read query spectra, store as clusters and sort cluster pointers according to m/z
	vector<SingleSpectrumHeader> queryHeaders;
	vector<Cluster> queryClusters;
	readQuerySpectraIntoClusters(params, queryClusters, queryHeaders, model);
	sort(queryClusters.begin(), queryClusters.end());

	for (size_t i=1; i<queryClusters.size(); i++)
		if (queryClusters[i].getClusterMOverZ()<queryClusters[i-1].getClusterMOverZ())
			error("Cluster order inconsistent!");

	// init dat (also splits into batches)
	vector<string> datPaths;
	generateArchiveDatPaths(datPaths);
	DatFileManager datManager;
	datManager.init(datPaths, config);

	size_t maxStorageSize = datManager.getMaxSpaceNeededToCover(params->mzWindow) + 10000;
	if (maxStorageSize < 25000)
		maxStorageSize = 25000;

	vector<SingleSpectrumHeader> archHeaders;
	vector<Cluster>				 archClusters;

	try {
		cout << "Allocating memory for maximal storage of " << maxStorageSize << " clusters." << endl;
		archHeaders.resize(maxStorageSize);
		archClusters.resize(maxStorageSize);
	}
	catch (...)
	{
		cout << "Problem..." << endl;
		cout << "M/Z window           = " << params->mzWindow << endl;
		cout << "Maximal storage size = " << maxStorageSize << endl;
		cout << "Header mem           = " << maxStorageSize * sizeof(SingleSpectrumHeader) << endl;
		cout << "Cluster mem		  = " << maxStorageSize * sizeof(Cluster) << endl;
		cout << "Total to cluster     = " << datManager.getTotalSpectraToCluster() << endl;
		cout << "batches:" << endl;
		datManager.printBatches();
		error("Could not allocate memory for archive clusters!");
	}


	map<string, int> matchedTitles, idTitles;
	if (params->idListPath.length()>0 && ! params->gotCompleteArchiveSearch)
		readIdsTitleFromIdFile(params, idTitles);


	// open res file
	
	vector<DatSpectrumStats> clustersToWrite;
	bool storeTitles = (params->idListPath.length()>0);

	string resPath = params->outDir + "/" + params->outputName + "_full_res.txt";
	string tmpResPath = resPath + ".tmp";
	ofstream ofs(tmpResPath.c_str());

	ofs << "MAIN_COLUMNS:\tIDX\tFILE\tSCAN\tM/Z\tNUM_RESULTS\t(NUM_COMPARED)\tTITLE" << endl;
	ofs << "RESULT_COLUMNS:\tCLUSTER_ID\tSIMILARITY\tP-VALUE\tCLUSTER_M/Z\tCLUSTER_SIZE\tMAX_DS\tFULL_COMPOSITION" << endl;
	
	const size_t numSpectra = queryClusters.size();
	size_t qIdx=0;
	size_t minStoragePos = 0; // where the minimal cluster that needs to be considered is stored
	size_t maxStoragePos = 0; // where the maximal cluster that needs to be considered is stored
	size_t lastStoragePos = 0;
	size_t totalRead=0;
	clusterIdx_t clusterCounter=0;
	bool gotMoreBatches = true;
	

	// vector lists for top idxes (each list holds the spectra idxs that have a certain top idx)
	const size_t maxIdx = computeMzIndex(2000.0, Cluster::getPeakIndexTolerance(), 0.0)+10;
	MyVectorAllocator<clusterIdx_t>          clusterIdxAllocation;
	vector< AllocatedVector<clusterIdx_t>* > clusterIdxLists(maxIdx,0); // lists for top peaks
	for (size_t i=0; i<maxIdx; i++)
		clusterIdxLists[i] = clusterIdxAllocation.allocateVector();

	for (qIdx=0; qIdx<numSpectra; qIdx++)
	{
		Cluster& queryCluster = queryClusters[qIdx];
		const mass_t minMzToCompare = queryCluster.getClusterMOverZ() - params->mzWindow;
		const mass_t maxMzToCompare = queryCluster.getClusterMOverZ() + params->mzWindow;

		// release clusters that no longer need to be compared
		if (lastStoragePos == 0 || archClusters[lastStoragePos-1].getClusterMOverZ()<maxMzToCompare)
		{
			// shunt the storage
			if (lastStoragePos>0)
			{
				size_t i;
				for (i=0; i<lastStoragePos; i++)
					if (archClusters[i].getClusterMOverZ() >= minMzToCompare)
						break;
				if (i==lastStoragePos)
				{
					lastStoragePos = 0;
					minStoragePos  = 0;
				}
				else
				{
					shuntInVector(archClusters,i,lastStoragePos-i);
					shuntInVector(archHeaders,i,lastStoragePos-i);
					lastStoragePos -= i;
					minStoragePos = 0;
					for (size_t j=0; j<lastStoragePos; j++)
						archClusters[j].setHeader(&archHeaders[j]); // correct header pointer

					for (maxStoragePos=0; maxStoragePos<lastStoragePos; maxStoragePos++)
						if (archClusters[maxStoragePos].getClusterMOverZ()>maxMzToCompare)
							break;
				}

				// remove older idxs from list
				const clusterIdx_t lowestIdx = archClusters[0].getClusterIdx();
				for (size_t i=0; i<clusterIdxLists.size(); i++)
				{
					if (! clusterIdxLists[i] || clusterIdxLists[i]->getSize() == 0)
						continue;
					const size_t size = clusterIdxLists[i]->getSize();
					const clusterIdx_t* idxs = clusterIdxLists[i]->getElementsPtr();
					const size_t pos = clusterIdxLists[i]->getFirstPositionGreaterOrEqual(lowestIdx);
					if (pos==0)
						continue;

					if (pos>=size)
					{
						clusterIdxLists[i]->relinquish();
						clusterIdxLists[i] = 0;
						continue;
					}
					clusterIdxLists[i]->shunt(pos);
				}
			}

			// read the new spectra into place and perform the relevant operations
			// collect all spectra that will be needed for next queries
			size_t q;
			for (q=qIdx+1; q<numSpectra; q++)
				if (queryClusters[q].getClusterMOverZ() - queryClusters[q-1].getClusterMOverZ() > params->mzWindow)
					break;

			mass_t maxMzForFill = (q < numSpectra ? queryClusters[q].getClusterMOverZ() + params->mzWindow : 1.0E6);
			vector<DatSpectrumStats> datStats;
			datManager.fillDatStats(datStats, minMzToCompare, maxMzForFill, maxStorageSize - lastStoragePos - 10);
			if (datStats.size()>0)
			{
				size_t numRead = datManager.readSpectraToClusterStorage(&archHeaders[lastStoragePos], &archClusters[lastStoragePos], datStats);
				for (size_t i=0; i<numRead; i++)
				{
					const clusterIdx_t clusterIdx = clusterCounter++;
					Cluster& cluster = archClusters[lastStoragePos+i];
					cluster.setClusterIdx(clusterIdx);

					// HACK if searching a partially annotated archive and --complete-search is not set
					// remove spectra that are not annotated from the search

					if (idTitles.size()>0 && ! params->gotCompleteArchiveSearch)
					{
						if (idTitles.find(cluster.getHeader()->getTitle()) == idTitles.end())
							continue;
					}
					

					// add the cluster top the relevant lists of top idxs

					// add the cluster index to the peak lists
					for (size_t j=0; j<NUM_PEAKS_FOR_HEURISTIC; j++)
					{
						const size_t idx = cluster.getTopPeakIdxs()[j];
						if (idx>0)
						{
							if (idx>=clusterIdxLists.size())
								clusterIdxLists.resize(idx + 1000);

							if (! clusterIdxLists[idx])
								clusterIdxLists[idx] = clusterIdxAllocation.allocateVector();
							clusterIdxLists[idx]->addElement(clusterIdx);
						}
					}
				}

				
				lastStoragePos += numRead;
				totalRead += numRead;
				size_t a = minStoragePos;
				size_t b = (lastStoragePos == 0 ? 0 : lastStoragePos-1);
				cout << qIdx << "\t" << totalRead << "\t" << a << " - " << b << 
					"\t( " << archClusters[a].getClusterMOverZ() << " - " << archClusters[b].getClusterMOverZ() << " )" << endl;
			}

			minStoragePos = 0;
			maxStoragePos = 0;
		}

		// update minSotragePos (the first position that needs to be considered)
		while (minStoragePos<lastStoragePos && archClusters[minStoragePos].getClusterMOverZ()<minMzToCompare)
			minStoragePos++;

		while (maxStoragePos<lastStoragePos && archClusters[maxStoragePos].getClusterMOverZ()<=maxMzToCompare)
			maxStoragePos++;
	
		if (maxStoragePos == lastStoragePos && maxStoragePos>0)
			maxStoragePos--;

		// create list of comparison cnadidates for the query spectrum
		vector<clusterIdx_t> archiveClustersToCompare;
		const clusterIdx_t lowestClusterIdx  = archClusters[minStoragePos].getClusterIdx();
		const clusterIdx_t highestClusterIdx = archClusters[maxStoragePos].getClusterIdx();
		const clusterIdx_t zeroPosClusterIdx = archClusters[0].getClusterIdx();
		const unsigned int* const topPeakIdxs = queryCluster.getTopPeakIdxs();
		assert(zeroPosClusterIdx<=lowestClusterIdx);
		assert(lowestClusterIdx<=highestClusterIdx);
		for (size_t i=0; i<NUM_PEAKS_FOR_HEURISTIC; i++)
		{
			const unsigned int peakIdx      = topPeakIdxs[i];
			const clusterIdx_t clusterIdx   = queryCluster.getClusterIdx();
			assert(peakIdx < clusterIdxLists.size());
			if (peakIdx <= 1)	// idx can be 0 if the spectrum has very few peaks
				continue;

			// always take the lists of the new clusters
			// take the clusters in idx-2, ..., idx+2
			for (size_t j=peakIdx-2; j<=peakIdx+2; j++)
			{
				if (clusterIdxLists[j] && clusterIdxLists[j]->getSize()>0)
				{
					const size_t listSize = clusterIdxLists[j]->getSize();
					const clusterIdx_t* elementList = clusterIdxLists[j]->getElementsPtr();

					// first position in the idxList that is greater than the lowest cluster idx
					size_t firstPosInList = clusterIdxLists[j]->getFirstPositionGreaterOrEqual(lowestClusterIdx);
					if (firstPosInList == MAX_SIZE_T)
						continue;
					if (lowestClusterIdx == 0)
						firstPosInList = 0;
					for (size_t k=firstPosInList; k<listSize; k++)
					{
						const clusterIdx_t clusterIdx = elementList[k];
						assert(clusterIdx >= lowestClusterIdx);
						// don't add pairs in which the lower cluster idx is equal or higher
						// to the upper cluster idx
						if (clusterIdx>highestClusterIdx)
							break;

						archiveClustersToCompare.push_back(clusterIdx);
					}
				}
			}
		}

		// compare the query cluster to the selected archive clusters
		sort(archiveClustersToCompare.begin(), archiveClustersToCompare.end());
		size_t numCompared=0;
		const DistancePeakList* queryDistancePeaks = queryCluster.getDistancePeaks();
		for (size_t i=0; i<archiveClustersToCompare.size(); i++)
		{
			if (i>0 && archiveClustersToCompare[i] == archiveClustersToCompare[i-1])
				continue;
			
			size_t pos = archiveClustersToCompare[i] - zeroPosClusterIdx;
			assert(archClusters[pos].getClusterIdx() == archiveClustersToCompare[i]);
			//float sim = computeSimilarity(queryDistancePeaks, archClusters[pos].getDistancePeaks(), Cluster::getPeakIndexToleranceAsInt());
			// NP3 GOT sim computation
			float sim = computeSimilarity(queryDistancePeaks, archClusters[pos].getDistancePeaks(), Cluster::getPeakIndexToleranceAsInt());
			queryCluster.updateSimilarityValues(sim, archiveClustersToCompare[i]);
			numCompared++;
		}

		if (numCompared == 0)
			continue;

		const float minSimNeeded = simModel.computeMinSimilarityAllowed(queryDistancePeaks->numPeaks, numCompared, params->maxPvalue);
		assert(minSimNeeded>0.0);
		vector<SimPair> simPairs;
		const float*        simValues   = queryCluster.getBestSimilarityValues();
		const clusterIdx_t* clusterIdxs = queryCluster.getBestSimilarityClusterIdxs();
		for (size_t i=0; i<NUM_TOP_SIMILARITIES_TO_SAVE; i++)
			if (simValues[i] >= minSimNeeded)
				simPairs.push_back(SimPair(simValues[i],clusterIdxs[i]));
		if (simPairs.size() ==0)
			continue;

		sort(simPairs.begin(), simPairs.end());

		const SingleSpectrumHeader* header = queryCluster.getHeader();
		ofs << qIdx << "\t" << header->getSpectraFileIndexInList() << "\t" << header->getScanNumber() << "\t";
		ofs << header->getMOverZ() << "\t" << simPairs.size() << "\t(" << numCompared << ")\t" << header->getTitle() << endl;
		for (size_t i=0; i<simPairs.size(); i++)
		{
			const Cluster& archCluster = archClusters[simPairs[i].idx - zeroPosClusterIdx];
			const string& title = archCluster.getHeader()->getTitle();
			ofs << title << "\t" 
				<< setprecision(3) << fixed << simPairs[i].sim << "\t"
				<< scientific << simModel.computePValue(queryDistancePeaks->numPeaks, numCompared, simPairs[i].sim) << "\t" 
				<< fixed << archCluster.getClusterMOverZ() << endl;

			matchedTitles[title]++;

			if (params->gotOutputMatchedSpectra)
			{
				DatSpectrumStats dss;
				dss.fileIdx      = archCluster.getHeader()->getSpectraFileIndexInList();
				dss.filePosition = archCluster.getHeader()->getPositionInFile();
				dss.mOverZ       = archCluster.getClusterMOverZ();
				clustersToWrite.push_back(dss);
			}
		}
	}
	ofs.close();

	// add clust info (cluster size, most abundant dataset, full dataset composition)
	string clustDir = archivePath_ + "/clust";

	if (checkDirPathExists(clustDir.c_str()))
	{
		cout << "Adding clust file data ... " << endl;
		vector<string> clustPaths;
		datManager.generateListOfClustPaths(clustPaths);
		string newPath = (storeTitles ? resPath + ".tmp2" : resPath);
		addClustInfoToResults(matchedTitles, clustPaths, tmpResPath, newPath); 
		unlink(tmpResPath.c_str());
		tmpResPath = newPath;
	}
	else
		cout << "No clust information found!" << endl;

	if (storeTitles)
	{
		cout << "Adding id data ... " << endl;
		// add titles to search results ....
		// pathToOpen --> resPath
		addIdsToResults(matchedTitles, params->idListPath, tmpResPath, resPath);
		unlink(tmpResPath.c_str());
	}
	

	if (params->gotOutputMatchedSpectra)
	{
		cout << "Writing matched consensus spectra ... " << endl;
		size_t numWritten = datManager.writeSpectraToOutput(params->outDir, params->outputName, clustersToWrite);
		cout << "Wrote " << numWritten << " spectra to mgf files." << endl;
	}

	cout << "End of archive search " << getTimeString() << endl;
}





/* This funciton is used to correct a bug in the clustering algorithm where
an incorrect m/z value was written into the dat files. To avoid re-clusteirng the
data, the correct m/z values are read from the clust files and written into the dat.
*/
void correctDatError(const char* datList, const Config* config)
{
	vector<string> datPaths;
	readListOfPaths(datList, datPaths);
	
	for (size_t i=0; i<datPaths.size(); i++)
	{
		string name,dir;
		getFileNameWithoutExtension(datPaths[i].c_str(), name);
		getDirFromFullPath(datPaths[i].c_str(), dir);
		string clustPath = dir.substr(0,dir.length()-3) + "clust/"+ name + ".clust";

		// read clust path
		ifstream ifs(clustPath.c_str(),ios::binary);
		if (! ifs.good())
			error("Could not open clust path: ",clustPath.c_str());

		map<string,mass_t> correctMzs;
		char buffer[256];
		while (ifs)
		{
			int n=0;
			mass_t mz=0;
			ifs.getline(buffer,256);
			if (! strncmp(buffer,"ARCHPNNL",8))
			{
				string title;
				istringstream iss(buffer);
				iss >> title >> n >> mz;
				if (mz<10.0 || mz>5000.0)
					error("line not parsed correctly: ",clustPath.c_str(), buffer);
				correctMzs[title]=mz;
			}
			for (int i=0; i<n; i++)
				ifs.getline(buffer,256);
		}
	
		string newPath = datPaths[i] + ".new";
		DatFile dat, newDat;
		dat.openForReading(datPaths[i].c_str());
		newDat.openForWriting(newPath.c_str());
		
		size_t totalSpectraRead = 0;
		size_t numBytes = 0;
		float  sqs;
		char* spec = 0;
		
			
		mass_t minMz = 1000000.0;
		mass_t maxMz = 0.0;
		bool continueReading = (dat.getNumSpectra()>0);
		while ( continueReading )
		{
			long peaksPosition=0;
			continueReading = dat.getNextDatSpectrum(spec, numBytes, &sqs, &peaksPosition);

			SingleSpectrumHeader ssh;
			ssh.setFileType(IFT_DAT);
			ssh.scanSpectrumHeaderFromBuffer(spec, config);
			string title = ssh.getTitle();
			map<string,mass_t>::const_iterator it = correctMzs.find(title);
			if (it == correctMzs.end())
				error("Could not find mz for ",title.c_str());

			size_t mzOffset = 2*sizeof(unsigned int);
			mass_t mz = it->second;
			*reinterpret_cast<mass_t*>(spec + mzOffset)=mz;
			newDat.addSpectrumToStats(mz, ssh.getOriginalNumPeaks());
			newDat.writeSpectrumFromBuffer(spec,numBytes);
		}
		newDat.closeAfterWriting();
		copyFile(newPath.c_str(),datPaths[i].c_str());
		unlink(newPath.c_str());
		cout << "processed " << datPaths[i] << endl;
	}	
}


void MsArchive::createArchiveFromMgf(MsParameterStruct* params, AllScoreModels* model)
{
	const Config* config = model->get_config();
	if (params->list.length() == 0)
		error("Must supply --list when creating archive from mgf");
	if (params->outputName.length() == 0)
		error("Must supply --ouput-name when creating archive from mgf");
	if (params->idListPath.length() == 0)
		error("Must supply --id-list when creating archive from mgf");

	// make archive
	archiveName_ = params->outputName;

	params->outputNameWithVersion = archiveName_ + "_0_0";

	// make sure dat dir is present
	createDirIfDoesNotExist(params->tmpDir.c_str(), params->verboseLevel);
	DatFileWriter dfw(model);
	const double startTime = time(NULL);

	// perform first pass
	params->datList = dfw.convertDataToDatFirstPass(params);

	string datDir  = params->outDir + "/" + params->outputNameWithVersion + "/dat";

	createDirIfDoesNotExist(params->outDir.c_str());
	string archDir = params->outDir + "/" + params->outputNameWithVersion;
	createDirIfDoesNotExist(archDir.c_str());
	createDirIfDoesNotExist(datDir.c_str(), params->verboseLevel);

	// perform second pass
	params->datList = dfw.convertDataToDatSecondPass(params->datList, datDir, params->outputName, 
		params->sqsThreshold, params->verboseLevel);

	DataSetRecord ds;
	ds.fileList            = params->list;
	ds.index               = 0;
	ds.fragTolerance       = config->getTolerance();
	ds.precursorTolerance  = config->get_pm_tolerance();
	ds.dateAdded		   = getDateAsInt();
	ds.numClusters		   = dfw.getNumSpectraWrittenSecondPass();
	ds.numSingletons       = dfw.getNumSpectraWrittenSecondPass();
	archiveDatasets_.push_back(ds);

	
	vector<string> datPaths;
	readListOfPaths(params->datList.c_str(),datPaths);
	
	for (size_t i=0; i<datPaths.size(); i++)
	{
		string fileName;
		getFileNameWithoutExtension(datPaths[i].c_str(), fileName);
		size_t pos = fileName.find_last_of("\\/");
		if (pos == string::npos)
		{
			datNames_.push_back(fileName);
		}
		else
			datNames_.push_back(fileName.substr(pos));

	}
	cout << "Generated " << datNames_.size() << " dat names" << endl;

	ClusterOutputter cs;
	writeArchive(cs, params, config);
}


