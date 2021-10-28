#include "ClusterOutputter.h"
#include "MsClusterAuxfuns.h"
#include "MsParameterStruct.h"

// static member declaration
ModClusterWriter OutputFilesForMzSlice::modClusterWriter;


void OutputFilesForMzSlice::init(const MsParameterStruct* params, const string& nameWithMassLabel)
{
	verboseLevel_ = params->verboseLevel;
	indWriteMgf_ = params->gotOutputMgf;
	indWriteDat_ = (params->gotCreateArchive || params->gotOutputDat);

	if (params->gotCreateArchive)
	{
#if defined(WIN32) || defined(WIN64)
		datPath_   = params->archiveOutDir + "\\dat\\" + nameWithMassLabel;
		clustPath_ = params->archiveOutDir + "\\clust\\" + nameWithMassLabel;
		mgfPath_   = params->archiveOutDir + "\\mgf\\" + nameWithMassLabel;
#else
		datPath_   = params->archiveOutDir + "/dat/" + nameWithMassLabel;
		clustPath_ = params->archiveOutDir + "/clust/" + nameWithMassLabel;
		mgfPath_   = params->archiveOutDir + "/mgf/" + nameWithMassLabel;
#endif
	}
	else
	{
#if defined(WIN32) || defined(WIN64)
		datPath_   = params->outDir + "\\dat\\" + nameWithMassLabel;
		clustPath_ = params->outDir + "\\clust\\" + nameWithMassLabel;
		mgfPath_   = params->outDir + "\\mgf\\" + nameWithMassLabel;
#else
		datPath_   = params->outDir + "/dat/" + nameWithMassLabel;
		clustPath_ = params->outDir + "/clust/" + nameWithMassLabel;
		mgfPath_   = params->outDir + "/mgf/" + nameWithMassLabel;
#endif
	}

	datPaths_.clear();
	mgfPaths_.clear();
	clustPaths_.clear();

	maxNumClustersPerFile_ = params->outputFileSize;
	indOutputModClusters_  = params->gotOutputModified;
	indIsInitialized_ = true;
}

void OutputFilesForMzSlice::open()
{
	if (runningFileIdx_ == MAX_SIZE_T)
	{
		runningFileIdx_=0;
	}
	else
		runningFileIdx_++;

	ostringstream nameWithNum;
	nameWithNum << clustPath_ << "_" << runningFileIdx_;

	const string clustName = nameWithNum.str() + ".clust";
	if (! clustOut_.open(clustName.c_str(), 1<<20))
		error("Couldn't open clust file:", clustName.c_str());
	clustPaths_.push_back(clustName);

	if (verboseLevel_>1)
		cout << "Opening: " << clustName << endl;

	if (indWriteMgf_)
	{
		ostringstream nameWithNum;
		nameWithNum << mgfPath_ << "_" << runningFileIdx_;
		const string mgfName = nameWithNum.str() + ".mgf";
		if (! mgfOut_.open(mgfName.c_str(), 1<<20))
			error("Couldn't open mgf file: ", mgfName.c_str());
		mgfPaths_.push_back(mgfName);

		if (verboseLevel_>1)
			cout <<"         " << mgfName << endl;
	}

	if (indWriteDat_)
	{
		ostringstream nameWithNum;
		nameWithNum << datPath_ << "_" << runningFileIdx_;
		const string datName = nameWithNum.str() + ".dat";
		if (! datOut_.openForWriting(datName.c_str()))
			error("Couldn't open dat file: ", datName.c_str());
		datPaths_.push_back(datName);
		if (verboseLevel_>1)
			cout <<"         " << datName << endl;
	}

	numClustersInFile_ = 0;
	indIsOpen_ = true;
}

void OutputFilesForMzSlice::close()
{
	if (clustOut_.isOpen())
		clustOut_.close();

	if (mgfOut_.isOpen())
		mgfOut_.close();

	if (datOut_.getIndOpen())
		datOut_.closeAfterWriting();

	indIsOpen_ = false;
}


void OutputFilesForMzSlice::writeClusterToOutputFiles(const Cluster& cluster,
													  const string& clustFileEntry,
													  const SingleSpectrumHeader* header,
													  bool  outputThisClusterToMod)
{
	static char* mgfBuffer=0;
	if (! indIsOpen_) {
        open();
        // NP3 adding CSV header to clust files
        ostringstream addHeader;
        addHeader.str(std::string());
        addHeader << "jobName, clustId, clustSize, rtMean, mzConsensus, chargeConsensus, rtMin, rtMax,precursorInt," <<
                "peptideStrConsensus,datasetIdx, fileIndex, scanNumber, peakArea, peakId, rt, mz, " <<
                "simToConsensus, pValue, charge, peptideStr\n" << clustFileEntry;
        string clustFileEntryHeader = addHeader.str();
        assert(header->getOriginalNumPeaks() == cluster.getNumPeaks());
        clustOut_.writeToBuffer(clustFileEntryHeader.c_str(), clustFileEntryHeader.length());
    } else {
        assert(header->getOriginalNumPeaks() == cluster.getNumPeaks());
        clustOut_.writeToBuffer(clustFileEntry.c_str(), clustFileEntry.length());
	}

	if (indWriteMgf_ || (indOutputModClusters_ && outputThisClusterToMod))
	{
		if (! mgfBuffer)
			mgfBuffer = new char[131072]; // should be enough for any MGF

		const size_t n=cluster.outputToMgfBuffer(mgfBuffer, header);
		if (n>=131072)
			error("MGF buffer in cluster out too small, increase size!");

		if (indOutputModClusters_ && outputThisClusterToMod)
			modClusterWriter.writeToBuffer(mgfBuffer, n);

		if (indWriteMgf_)
			mgfOut_.writeToBuffer(mgfBuffer, n);
	}

	if (indWriteDat_)
	{
		// When writing a cluster spectrum as a dat (wether it is a singleton cluster or a larger one),
		// we do not want to have information about the spectraFileIndexInList_ . This will only cause
		// confusion: there is no single spectraFileIndexInList for a cluster of more than 1 spectrum;
		// the index in the dat list is not something you write in the dat file (since we can join lists
		// of dat files), this value is given when the dat is read (and then we know what is the true index
		// of the dat file in the dat file list).
		
		const int originalSpectraFileIndexInList = header->getSpectraFileIndexInList();
		SingleSpectrumHeader*  varHeader = const_cast<SingleSpectrumHeader*>(header);
		varHeader->setMOverZ(cluster.getClusterMOverZ()); // make sure we are using cluster m/z
		varHeader->setSpectraFileIndexInList(MIN_INT);
		varHeader->setScanNumber(MIN_INT);
		varHeader->setOriginalNumPeaks(cluster.getNumPeaks()); // update the number of peaks since this is the number
															  // that gets written in the dat header
		datOut_.writePeakList(cluster, header);
		varHeader->setSpectraFileIndexInList(originalSpectraFileIndexInList);
	}

	if (++numClustersInFile_ == maxNumClustersPerFile_)
		close();
}

void ClusterOutputter::init(const MsParameterStruct* params, mass_t maxMOverZ)
{
	params_		  = params;
	verboseLevel_ = params->verboseLevel;

	size_t size = computeMzIndex(maxMOverZ, params->majorIncrement, Cluster::getPeakTolerance()) + 2;
	mzIncrement_ = params->majorIncrement;
	numClustersWritten_  = 0;
	totalSpectraCount_	 = 0;
	lastIndexClosed_	 = 0;

	clusterCounts_.clear();
	clusterCounts_.resize(numSizeThresholdsForReport, 0);
	singletonCounts_.clear();
	singletonCounts_.resize(numSizeThresholdsForReport, 0);
	totalSizeCounts_.clear();
	totalSizeCounts_.resize(numSizeThresholdsForReport, 0);

	outputFiles_.resize(size);

#if defined(WIN32) || defined(WIN64)
	string clustDir = (params->gotCreateArchive? params->archiveOutDir + "\\clust" : params->outDir + "\\clust");
#else
	string clustDir = (params->gotCreateArchive? params->archiveOutDir + "/clust " : params->outDir + "/clust");
#endif
	createDirIfDoesNotExist(clustDir.c_str());

	if (params->gotOutputMgf)
	{
#if defined(WIN32) || defined(WIN64)
	string mgfDir = (params->gotCreateArchive? params->archiveOutDir + "\\mgf" : params->outDir + "\\mgf");
#else
	string mgfDir = (params->gotCreateArchive? params->archiveOutDir + "/mgf" : params->outDir + "/mgf");
#endif
		createDirIfDoesNotExist(mgfDir.c_str());
	}

	if (params->gotCreateArchive || params->gotOutputDat)
	{
#if defined(WIN32) || defined(WIN64)
	string datDir = (params->gotCreateArchive? params->archiveOutDir + "\\dat" : params->outDir + "\\dat");
#else
	string datDir = (params->gotCreateArchive? params->archiveOutDir + "/dat" : params->outDir + "/dat");
#endif
		createDirIfDoesNotExist(datDir.c_str());
	}

	for (size_t i=0; i<size; i++)
	{
		ostringstream oss;
		oss << params->outputNameWithVersion << "_" << static_cast<size_t>(i * mzIncrement_ + 0.00001);
		outputFiles_[i].init(params_, oss.str());
	}

	if (params_->gotOutputModified)
	{
		string path = params_->outDir + std::string("/") + params->outputNameWithVersion + std::string("_mod");
		OutputFilesForMzSlice::modClusterWriter.init(path, params->outputFileSize);
	}
}

// this function changes the information in the header belonging to the
// cluster idx to reflect the updated cluster information (num peaks, m/z, etc.)
void ClusterOutputter::outputCluster(const Cluster& cluster, 
									 const string&  clustFileEntry,
									 bool  indCloseOutputFiles,
									 bool  outputToMod)
{
	const size_t idx = computeMzIndex(cluster.getClusterMOverZ(), mzIncrement_, Cluster::getPeakTolerance());
	assert(idx<outputFiles_.size());
	
	outputFiles_[idx].writeClusterToOutputFiles(cluster, clustFileEntry, cluster.getHeader(), outputToMod);

	// add cluster to statistics
	const clusterIdx_t size = cluster.getNumSingletonsIncludingSelf();
	size_t sizeIndex;
	for (sizeIndex=0; sizeIndex<numSizeThresholdsForReport-1; sizeIndex++)
		if (size<=sizeThresholdsForReport[sizeIndex])
			break;
	clusterCounts_[sizeIndex]++;
	singletonCounts_[sizeIndex]+= size;
	totalSizeCounts_[sizeIndex]+= cluster.getClusterSize();

	totalSpectraCount_ += cluster.getClusterSize();
	numClustersWritten_++;

	// close files to avoid having too many open ones
	// we assume that there is no way that there should be more than 15
	// consecutive indexes open (that represents somehow out of order clusters with
	// a difference of 10 Da in their precursor mass - this shouldn't happen in a clustering job.
	if (indCloseOutputFiles && idx>lastIndexClosed_+15)
	{
		size_t i;
		for (i=lastIndexClosed_+1; i<=idx-15; i++)
			outputFiles_[i].close();
		lastIndexClosed_ = i-1;
	}
}


// closes files and collects paths
void ClusterOutputter::closeAll()
{
	size_t n=0;
	for (size_t i=0; i<outputFiles_.size(); i++)
	{
		if (outputFiles_[i].getIndIsOpen())
			outputFiles_[i].close();

		if (outputFiles_[i].getRunningFileIdx()< MAX_SIZE_T)
			n+= outputFiles_[i].getRunningFileIdx() +1;
	}

	if (! params_)
		return;
	
	allClustPaths_.clear();
	allClustPaths_.reserve(n);
	if (params_->gotCreateArchive || params_->gotOutputDat)
	{
		allDatPaths_.clear();
		allDatPaths_.reserve(n);
	}

	if (params_->gotOutputMgf)
	{
		allMgfPaths_.clear();
		allMgfPaths_.reserve(n);
	}

	for (size_t i=0; i<outputFiles_.size(); i++)
	{
		if (outputFiles_[i].getRunningFileIdx() < MAX_SIZE_T)
		{
			for (size_t j=0; j<=outputFiles_[i].getRunningFileIdx(); j++)
			{
				allClustPaths_.push_back(outputFiles_[i].getClustPaths()[j]);
				if (params_->gotOutputMgf)
					allMgfPaths_.push_back(outputFiles_[i].getMgfPaths()[j]);
				if (params_->gotCreateArchive || params_->gotOutputDat)
					allDatPaths_.push_back(outputFiles_[i].getDatPaths()[j]);
			}
		}
	}

	if (params_->gotOutputModified)
	{
		OutputFilesForMzSlice::modClusterWriter.close();
		if (verboseLevel_>0)
			cout << "Wrote " << OutputFilesForMzSlice::modClusterWriter.getFileIdx()+1 << " files of modified clusters ("
				 << OutputFilesForMzSlice::modClusterWriter.getNumClusters() << " clusters)" << endl;
	}

	if (verboseLevel_>1)
	{
		cout << endl;
		cout << "Summary of output files:" << endl;
		cout << "------------------------" << endl;
		cout << "clust files written: " << allClustPaths_.size() << endl;
		cout << "mgf   files writeen: " << allMgfPaths_.size() << endl;
		cout << "dat   files writeen: " << allDatPaths_.size() << endl << endl;
	}

	// write lists of files
	const string outDir = (params_->gotCreateArchive? params_->archiveOutDir  : params_->outDir );
	string clustList =  outDir + "/" + params_->outputNameWithVersion + "_clust_list.txt";
	writeListOfPaths(clustList.c_str(), allClustPaths_);

	if (params_->gotOutputMgf)
	{
		string mgfList = outDir + "/" + params_->outputNameWithVersion + "_mgf_list.txt";
		writeListOfPaths(mgfList.c_str(), allMgfPaths_);
	}

	if (params_->gotCreateArchive || params_->gotOutputDat)
	{
		string datList = outDir + "/" + params_->outputNameWithVersion + "_dat_list.txt";
		writeListOfPaths(datList.c_str(), allDatPaths_);
	}
}


void ClusterOutputter::closeAllWithMassBelow(mass_t mz)
{
	const size_t idx = computeMzIndex(mz, mzIncrement_, Cluster::getPeakTolerance());
	assert (idx<outputFiles_.size());

	if (idx<3)
		return;

	for (size_t i=0; i<idx-2; i++)
	{
		if (outputFiles_[i].getIndIsOpen())
			outputFiles_[i].close();
	}
}

/****************************************************************************
Reports on the number of clusters in each size bin
*****************************************************************************/
void ClusterOutputter::reportClusterStats() const
{
	cout << endl << "Final clustering report:" << endl;
	cout		 << "------------------------" << endl << endl;
	
	longInt8_t   sumTotalSizeCounts=0;
	clusterIdx_t sumSingletons     =0;
	clusterIdx_t sumClusters	   =0;
	for (size_t i=0; i<numSizeThresholdsForReport; i++)
	{
		sumClusters+=clusterCounts_[i];
		sumSingletons+=singletonCounts_[i];
		sumTotalSizeCounts+=totalSizeCounts_[i];
	}

	if (sumTotalSizeCounts>sumSingletons)
	{
		cout << "size\t#clust\t(%)\t#sing\t#total\t(%)" << endl;
		for (size_t i=0; i<numSizeThresholdsForReport; i++)
		{
			if (i==0)
			{
				cout << sizeThresholdsForReport[i];
			}
			else if (i==numSizeThresholdsForReport-1)
			{
				cout << sizeThresholdsForReport[numSizeThresholdsForReport-2]+1 << "+";
			}
			else if (sizeThresholdsForReport[i]-sizeThresholdsForReport[i-1] == 1)
			{
				cout << sizeThresholdsForReport[i];
			}
			else
				cout << sizeThresholdsForReport[i-1]+1 << "-" << sizeThresholdsForReport[i];
		
			cout << "\t" << clusterCounts_[i] << "\t(" << static_cast<double>(clusterCounts_[i])/static_cast<double>(sumClusters)
				<< ")\t"<< singletonCounts_[i] << "\t"<< totalSizeCounts_[i] << "\t(" << setprecision(3) 
				<< fixed << static_cast<double>(totalSizeCounts_[i])/static_cast<double>(sumTotalSizeCounts) 
				<< ")" << endl;
		}

		cout << endl;
		cout << "Total: " << sumClusters << "\t1.0\t" << sumSingletons << "\t" << sumTotalSizeCounts << "\t(1.0)" << endl;
	}
	else
	{
		cout << "size\t#clust\t(%)\t#sing\t(%)" << endl;
		for (size_t i=0; i<numSizeThresholdsForReport; i++)
		{
			if (i==0)
			{
				cout << sizeThresholdsForReport[i];
			}
			else if (i==numSizeThresholdsForReport-1)
			{
				cout << sizeThresholdsForReport[numSizeThresholdsForReport-2]+1 << "+";
			}
			else if (sizeThresholdsForReport[i]-sizeThresholdsForReport[i-1] == 1)
			{
				cout << sizeThresholdsForReport[i];
			}
			else
				cout << sizeThresholdsForReport[i-1]+1 << "-" << sizeThresholdsForReport[i];
		
			cout << "\t" << clusterCounts_[i] << "\t(" << setprecision(3) << fixed 
				 << static_cast<double>(clusterCounts_[i])/static_cast<double>(sumClusters)
				 << ")\t" << singletonCounts_[i] << "\t(" 
				 << static_cast<double>(singletonCounts_[i])/static_cast<double>(sumSingletons) 
				 << ")" << endl;
		}
		cout << endl;
		cout << "Total: " << sumClusters << "\t1.0\t" << sumSingletons << "\t(1.0)" << endl;
	}
}

