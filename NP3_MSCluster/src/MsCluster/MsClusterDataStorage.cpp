#include "MsClusterDataStorage.h"
#include "MsClusterAuxfuns.h"
#include "MsParameterStruct.h"
#include "MsSimilarityModel.h"

MsClusterDataStorage::~MsClusterDataStorage()
{
	clusterAlloc_.clear();
	headerAlloc_.clear();
	peakAlloc_.clear();
}


void MsClusterDataStorage::initialize(MsParameterStruct* params, 
									  const Config* config,
									  const MsSimilarityModel* simModel)
{
	params_		  = params;
	outDir_		  = params->outDir;
	config_		  = config;
	datasetIdx_   = params->datasetIdx;
	batch_		  = params->batchIdx;
	verboseLevel_ = params->verboseLevel;

	simModel_	  = simModel;

	runningClusterIdx_ = 0;
    // NP3 GOT change runningOutputIdx_ from 0 to 1, preventing SCANS=0
	runningOutputIdx_  = 1;
	totalSpectraWritten_ = 0;
	totalClustersWritten_ = 0;

	ostringstream oss;
	oss << params->outputName << "_" << datasetIdx_ << "_" << batch_;
	nameStr_ = oss.str();

	const mass_t tolerance = (params->fragmentTolrance > 0.0 ? params->fragmentTolrance : config->getTolerance());
	Cluster::setTolerances(tolerance, params->precursorPPMs);
	// NP3 GOT RT tolerance and scale factor
	Cluster::setRtTolerance(params->rtTolerance);
	Cluster::setScaleFactor(params->scaleFactor);

	peakWeightTable_.initWeights(64, 0.15);

	if (verboseLevel_>0)
		cout << endl << "Initalizing storage, examining dat files..." << endl;
	
	oss.str("");

	if (params->gotMergeArchives)
	{
		oss << params->outputName << "_" << archive1_.getGeneration() << "_" << params->batchIdx;
	}
	else
		oss << params->outputName << "_" << params->datasetIdx << "_" << params->batchIdx;


	params_->outputNameWithVersion = oss.str();	

	if (params->gotCreateArchive)
	{
#if  defined(WIN32) || defined(WIN64)
		params_->archiveOutDir = params->outDir + "\\" + params_->outputNameWithVersion;
#else
		params_->archiveOutDir = params->outDir + "/" + params_->outputNameWithVersion;
#endif
		createDirIfDoesNotExist(params_->archiveOutDir.c_str(), params->verboseLevel);
	}
	
	if (params_->gotMergeArchives)
	{
		datFileManager_.initFromTwoArchives(params->datPaths1, params->datPaths2, config);
	}
	else
		datFileManager_.init(params->datList, config);

	if (datFileManager_.getTotalSpectraToCluster() == 0)
	{
		cout << "No spectra found in dat files... exiting!" << endl;
		exit(0);
	}

	const size_t maxIdx = computeMzIndex(datFileManager_.getMaxMOverZ()+10.0, Cluster::peakIndexTolerance_, 0.0)+2;

	// create allocaed vector lists
	newClusterIdxsLists_.resize(maxIdx,0);
	for (size_t i=0; i<newClusterIdxsLists_.size(); i++)
		newClusterIdxsLists_[i] = clusterIdxAllocation_.allocateVector();


	clusterOutputter_.init(params, datFileManager_.getMaxMOverZ());

	setMinSqsForSingleton(params->sqsThreshold);

	if (verboseLevel_>0)
	{
		cout << endl << "Found " << datFileManager_.getTotalSpectraToCluster() << 
			" spectra, in " << datFileManager_.getDatFiles().size() << " dat files, precursor m/z range " << fixed << setprecision(2) << 
			datFileManager_.getMinMOverZ() << " - " << datFileManager_.getMaxMOverZ() << endl << endl;

	}

	if (verboseLevel_>0 && params->sqsThreshold>0.0)
	{
		cout << "Minimal SQS quality score for all spectra   : " << params->sqsThreshold << endl;
		if (minSqsForSingleton_ > params->sqsThreshold)
			cout << "Unpaired singletons need an SQS of at least : " << minSqsForSingleton_ << endl;
		cout << endl;
	}
}

size_t MsClusterDataStorage::findPosOfClusterWithMzAbove(mass_t mz, size_t startPos) const
{
	static Cluster dummyCluster;
	dummyCluster.clusterMOverZ_ = mz;
	if (nextClusterPos_==0)
		return MAX_SIZE_T;

	vector<Cluster>::const_iterator it = upper_bound(clusterAlloc_.begin() + startPos, 
													 clusterAlloc_.begin() + nextClusterPos_-1,
													 dummyCluster);
	if (it == clusterAlloc_.end())
		return MAX_SIZE_T;
	return distance(clusterAlloc_.begin(), it);
}


// makes a list of all pairs of psoitions that need to have their similarity computed
// each pair has one element from the range [firstPosLower,  lastPosLower)
// and the second element form the range    [firstPosHigher, lastPosHigher)
// the pairs are stored and sorted in the vector similarityPairPositions_
// returns actual value of lastHigherPos (if the list grows too long
// and we need to terminate early)
size_t MsClusterDataStorage::makeListOfSimilarities(size_t firstPosLower,  size_t lastPosLower,
													size_t firstPosHigher, size_t lastPosHigher,
													mass_t window)
{
	similarityPairPositions_.clear();
	const clusterIdx_t clusterIdxAllocStart = clusterAlloc_[0].getClusterIdx();
	const clusterIdx_t lastIdx	     = clusterAlloc_[lastPosHigher-1].getClusterIdx() + 1;

	size_t lowerPosStart = firstPosLower;
	size_t higherPos;
	for (higherPos = firstPosHigher; higherPos<lastPosHigher; higherPos++)
	{
		if (! clusterAlloc_[higherPos].getIndInPlay())
			continue;

		const size_t similarityListSizeBefore = similarityPairPositions_.size();
		const mass_t higherMz = clusterAlloc_[higherPos].getClusterMOverZ();
		while (lowerPosStart < higherPos && higherMz - clusterAlloc_[lowerPosStart].getClusterMOverZ() > window)
		{
			lowerPosStart++;
		}

		if (lowerPosStart == higherPos)
			continue;

		assert(higherMz - clusterAlloc_[lowerPosStart].getClusterMOverZ() <= window); 
		const clusterIdx_t lowestClusterIdx = clusterAlloc_[lowerPosStart].getClusterIdx();
		for (size_t i=0; i<NUM_PEAKS_FOR_HEURISTIC; i++)
		{
			const Cluster& higherPosCluster = clusterAlloc_[higherPos];
			const unsigned int idx          = higherPosCluster.topPeakIdxs_[i];
			const clusterIdx_t higherPosIdx = higherPosCluster.getClusterIdx();
			const size_t numHigerPosDistancePeaks = higherPosCluster.getDistancePeaks()->numPeaks;
			const short	higherClusterArchiveSource = (higherPosCluster.getHeader()->getArchiveSource() != 0 ? 
													  higherPosCluster.getHeader()->getArchiveSource() : 999);
			
			assert(idx < newClusterIdxsLists_.size());

			if (idx <= 1)	// idx can be 0 if the spectrum has very few peaks
				continue;

			// always take the lists of the new clusters
			// take the clusters in idx-2, ..., idx+2
			for (size_t j=idx-2; j<=idx+2; j++)
			{
				if (newClusterIdxsLists_[j] && newClusterIdxsLists_[j]->getSize()>0)
				{
					const size_t listSize = newClusterIdxsLists_[j]->getSize();
					const clusterIdx_t* elementList = newClusterIdxsLists_[j]->getElementsPtr();

					// first position in the idxList that is greater than the lowest cluster idx
					size_t firstPosInList = newClusterIdxsLists_[j]->getFirstPositionGreaterOrEqual(lowestClusterIdx);
					if (firstPosInList == MAX_SIZE_T)
						continue;
					if (lowestClusterIdx == 0)
						firstPosInList = 0;

					// if there are no clusters with a lower idx to compare to the higherPos cluster
					// at this mass we can skip this list
					if (elementList[firstPosInList] >= higherPosIdx)
						continue;
					
					assert( getCluster(elementList[firstPosInList])->getClusterMOverZ() + window >= higherMz);
					
					const  size_t lastPosInList = newClusterIdxsLists_[j]->getSize();
					for (size_t k=firstPosInList; k<lastPosInList; k++)
					{
						// don't add pairs in which the lower cluster idx is equal or higher
						// to the upper cluster idx
						if (elementList[k]>=higherPosIdx)
							break;

						// for each position added, make sure that this cluster is still a valid
						// one to be compared to (it could have been added to another cluster)
						const size_t lowerPos = elementList[k]-clusterIdxAllocStart;
						assert(lowerPos<higherPos);

						// if two archives are not used thaen lowerPosCluster.getHeader()->getArchiveSource() will always be
						// different from higherClusterArchiveSource since higherClusterArchiveSource=999
						const Cluster& lowerPosCluster = clusterAlloc_[lowerPos];
						if (lowerPosCluster.getIndInPlay() && lowerPosCluster.getHeader()->getArchiveSource() != higherClusterArchiveSource)
						{
							const size_t numLowerPosDistancePeaks = lowerPosCluster.getDistancePeaks()->numPeaks;
							if ((numHigerPosDistancePeaks > numLowerPosDistancePeaks + ALLOWED_DIFF_IN_NUM_DISTANCE_PEAKS) ||
								(numLowerPosDistancePeaks > numHigerPosDistancePeaks + ALLOWED_DIFF_IN_NUM_DISTANCE_PEAKS) )
							{
							    //cout << "@@@@@ poucos picos diff :" << numHigerPosDistancePeaks - numLowerPosDistancePeaks << endl;
                                continue;
                            }

                            // don't compare clusters that are too far away from each other
                            //if (abs(clusterAlloc_[lowerPos].getClusterRT() - clusterAlloc_[higherPos].getClusterRT()) > Cluster::getRtTolerance())
                            //    continue;
                            // NP3 GOT rt width comparision
							if (!((clusterAlloc_[lowerPos].getClusterRT() >= clusterAlloc_[higherPos].getClusterRTMin() - Cluster::getRtTolerance() &&
								   clusterAlloc_[lowerPos].getClusterRT() <= clusterAlloc_[higherPos].getClusterRTMax() + Cluster::getRtTolerance()) ||
								  (clusterAlloc_[higherPos].getClusterRT() >= clusterAlloc_[lowerPos].getClusterRTMin() - Cluster::getRtTolerance() &&
								   clusterAlloc_[higherPos].getClusterRT() <= clusterAlloc_[lowerPos].getClusterRTMax() + Cluster::getRtTolerance())))
								continue;

							similarityPairPositions_.push_back(ClusterPositionPair(lowerPos,higherPos));
						}
					}
				}
			}
		}
		

		// check if there are redundancies in the list, after adding all
		// the pairs for the cluser, and remove them
		const size_t similarityListSizeAfter = similarityPairPositions_.size();
		if (similarityListSizeAfter - similarityListSizeBefore > 1)
		{
			sort(similarityPairPositions_.begin() + similarityListSizeBefore,
					 similarityPairPositions_.begin() + similarityListSizeAfter);
			size_t trailIdx = similarityListSizeBefore;
			size_t leadIdx  = similarityListSizeBefore;
			const size_t lastIdx = similarityListSizeAfter;
			while (++leadIdx<lastIdx)
			{
				if (similarityPairPositions_[leadIdx] != similarityPairPositions_[trailIdx])
					similarityPairPositions_[++trailIdx] = similarityPairPositions_[leadIdx];
			}
			similarityPairPositions_.resize(trailIdx+1);
		}

		if (similarityPairPositions_.size()>MAX_SIMILARITY_LIST_SIZE)
			break;
	}

	sort(similarityPairPositions_.begin(),similarityPairPositions_.end());

	return higherPos;
}

void MsClusterDataStorage::computeMinSimilarityForJoiningFromLists(double mixedClusterProb, size_t lowerPos, size_t higherPos)
{
	vector<int> counts(higherPos+100,0);
	for (size_t i=0; i<similarityPairPositions_.size(); i++)
	{
		counts[similarityPairPositions_[i].lowerPos]++;
		counts[similarityPairPositions_[i].higherPos]++;
	}

	for (size_t i=lowerPos; i<=higherPos; i++)
	{
		if (counts[i] == 0) // value should already be 2.0 any way
			continue; 

		const float simVal = simModel_->computeMinSimilarityAllowed( clusterAlloc_[i].getDistancePeaks()->numPeaks,
							    								    counts[i],
																	mixedClusterProb);
		if (simVal > clusterAlloc_[i].minSimilarityToJoin_ || clusterAlloc_[i].minSimilarityToJoin_ == 2.0)
		{
			clusterAlloc_[i].minSimilarityToJoin_ = simVal;
			clusterAlloc_[i].numSimilarityPairs_  = counts[i];
		}
	}
}


void MsClusterDataStorage::allocateMemory(float availableGB, bool verbose)
{
	if (sizeof(size_t)<8 && availableGB>3.2)
		error("Cannot allocate more than 3.2 GB with this architecture because sizeof(size_t)<8");
	
	// assume that storage will take 85% of the memory, leave memory free for other runtime allocations
	const size_t availableBytes    = static_cast<size_t>(1073741824.0 * availableGB * 0.86);
	const size_t cost = 200*sizeof(Peak) + 2*sizeof(Cluster) + 2*sizeof(SingleSpectrumHeader);
	const size_t costSimilarityList = static_cast<size_t>(MAX_SIMILARITY_LIST_SIZE*1.1);
	const size_t n = (availableBytes - costSimilarityList)/cost;

	cout << "Allocating " << setprecision(3) << availableGB << " GB." << endl;
	if (verbose)
	{
		cout << endl;
		cout << "Bytes for data storage       = " << availableBytes << endl;
        cout << "sizeof(MAX_SIMILARITY_LIST)  = " << costSimilarityList << endl;
		cout << "sizeof(Peak)                 = " << sizeof(Peak) << endl;
		cout << "sizeof(Cluster)              = " << sizeof(Cluster) << endl;
		cout << "sizeof(SingleSpectrumHeader) = " << sizeof(SingleSpectrumHeader)  << endl;
		cout << "Total cost per singleton     = " << cost << endl;
		cout << "Allocating room for N = " << n << " singletons (" << n * cost << " bytes)" << endl;
		cout << "This means that up to " << n/2 << " spectra can simultaneously be evaluated." << endl << endl;
	}

	try 
	{
		headerAlloc_.resize(2*n); // original singletion headers are sttored here
		peakAlloc_.resize(200*n);  // original singletion peaks are stored here
		clusterAlloc_.resize(2*n);
		similarityPairPositions_.resize(costSimilarityList);
	}
	catch ( ... )
	{
		error("Could not allocate memory for MsCluster data storage, try smaller allocation size!");
	}

	cout << "Allocated memory successfully..." << endl << endl;
}



// adds new batches of spectra for clustering.
// first removes some spectra to make room for the new ones.
// returns the number of spectra added. Added spectra are also entered as clusters.
clusterIdx_t MsClusterDataStorage::addNewSpectra( mass_t window, clusterIdx_t* returnedNumSpectraFreed)
{
	clusterIdx_t	numSpectraNeeded =0;
	longInt8_t		numPeaksNeeded   =0;
	datFileManager_.getRequirementsForNextBatch(numSpectraNeeded, numPeaksNeeded);
	if (numSpectraNeeded == 0)
		return 0;

	assert( nextClusterPos_ <= clusterAlloc_.size());
	assert( nextPeakPos_ <= peakAlloc_.size() );

	// compute how many spectra need to be freed
	// first count space that is not filled towards the end of the allocated space
	const clusterIdx_t clustersLeft = clusterAlloc_.size() - nextClusterPos_;
	numSpectraNeeded = (clustersLeft >= numSpectraNeeded ? 0 : numSpectraNeeded - clustersLeft);
	
	const longInt8_t peaksLeft =  peakAlloc_.size() - nextPeakPos_;
	numPeaksNeeded   = ( peaksLeft >= numPeaksNeeded ? 0 : numPeaksNeeded - peaksLeft);

	// find the number of spectra that will actually have to be evacuated (shifted out of clusterAlloc_)
	size_t numSpectraRequiredToFree =0; 
	if (numSpectraNeeded > 0 || numPeaksNeeded > 0)
	{
		size_t sumPeaks =0;
		while ( numSpectraRequiredToFree < nextClusterPos_ && // can't free more than are there
			   (sumPeaks < numPeaksNeeded || numSpectraRequiredToFree < numSpectraNeeded) )
		{
			sumPeaks += headerAlloc_[numSpectraRequiredToFree++].getOriginalNumPeaks();
		}
	}

	// we don't mind freeing this number of clusters since they are already compared to all others
	// in the m/z window
	clusterIdx_t numSpectraFullyCompared=0;
	if (runningClusterIdx_ > 0 && nextClusterPos_>0)
	{
		size_t pos=nextClusterPos_-1;
		const mass_t maxMzInAlloc = clusterAlloc_[pos].getClusterMOverZ();
		while (pos>0 &&  maxMzInAlloc - clusterAlloc_[pos].getClusterMOverZ() < window)
			pos--;

		while (pos<nextClusterPos_-1 && ! clusterAlloc_[pos].getIndInPlay())
			pos++;

		numSpectraFullyCompared = pos;
	}

	// the idx in clusterAlloc_ that will become 0 after the shifting
	// At first we will consider freeing the required number of clusters
	clusterIdx_t finalNumSpectraToFree = numSpectraFullyCompared;
	
	// check if we need to free more than we can
	if (numSpectraRequiredToFree>0 && 
		numSpectraFullyCompared <= numSpectraRequiredToFree)
	{
		// we will remove at most a 1/3 of the clusters that were not fully compared
		finalNumSpectraToFree = numSpectraFullyCompared + (nextClusterPos_ - numSpectraFullyCompared)/3;
		if (finalNumSpectraToFree >= numSpectraRequiredToFree)
			finalNumSpectraToFree = numSpectraRequiredToFree;

		size_t numPeaksNeeded = 0;
	}

	assert(finalNumSpectraToFree ==0 || finalNumSpectraToFree < nextClusterPos_ );
	
	// shifts all the clusters, headers, etc. so the cluster in position 
	// finalFreeIdx moves to position 0
	shiftStorage(finalNumSpectraToFree, window);

	if (returnedNumSpectraFreed)
		*returnedNumSpectraFreed = finalNumSpectraToFree;

	assert(clusterAlloc_.size() >= nextClusterPos_);
	assert(peakAlloc_.size() >= nextPeakPos_);

	// get spectra stats
	static vector<DatSpectrumStats> newSpectraStats;
	newSpectraStats.reserve(clusterAlloc_.size());
	const bool reachedEnd = datFileManager_.getNewDatSpectraStats(clusterAlloc_.size() - nextClusterPos_,
										  peakAlloc_.size() - nextPeakPos_,
										  newSpectraStats);

	assert(newSpectraStats.size() + nextClusterPos_ <= clusterAlloc_.size());
	assert( reachedEnd || newSpectraStats.size() + nextClusterPos_ >= 0);

	// read spectra
	size_t numSpectraReadFromDat = 0;
	size_t numPeaksReadFromDat = 0;
	datFileManager_.readSpectraToStorage(&headerAlloc_[nextClusterPos_], 
										 &peakAlloc_[nextPeakPos_], 
										 newSpectraStats,
										 numSpectraReadFromDat,
										 numPeaksReadFromDat);

	assert( numSpectraReadFromDat == newSpectraStats.size() );
	assert( nextClusterPos_ + newSpectraStats.size() <= clusterAlloc_.size());
	assert( nextPeakPos_ + numPeaksReadFromDat <= peakAlloc_.size() );
	
	// create clusters
	sort(newSpectraStats.begin(), newSpectraStats.end());
	for (size_t i=0; i<newSpectraStats.size(); i++)
	{
		const size_t clusterPos = nextClusterPos_ + i;
		const size_t peakPos	= nextPeakPos_ + newSpectraStats[i].peakWritePos;
		const clusterIdx_t clusterIdx = runningClusterIdx_++;

		assert( newSpectraStats[i].clusterWritePos == clusterPos-nextClusterPos_);
		assert( clusterPos < clusterAlloc_.size() );

		const SingleSpectrumHeader* header = &headerAlloc_[clusterPos];
		const Peak*					peaks  = &peakAlloc_[peakPos];
		const int					numPeaks = newSpectraStats[i].numPeaks;

		bool createdOk = clusterAlloc_[clusterPos].createNewCluster(clusterIdx, header, peaks, numPeaks);
		if (! createdOk)
		{
			cout << "BAD creation " << clusterIdx << " " << header->getTitle() << endl;
			clusterAlloc_[clusterPos].setIndInPlay(false);
			false;
		}
		// NP3 check if the cluster precursor intensity is bellow the baseline, if true remove it from play - GOT BASELINE (originally <= 0.0)
		if (header->getPrecursorIntensity() < INTENSITY_BASELINE)
		{
			cout << "Precursor intensity bellow baseline for " << clusterIdx << " " << header->getTitle() << endl;
			clusterAlloc_[clusterPos].setIndInPlay(false);
			continue;
		}

		// if we are not assigning charges according to the values in the spectra, 
		// all charges should be 0
		if (! params_->gotAssignCharges && header->getPeptideStr().length() == 0)
		{
			clusterAlloc_[clusterPos].clusterCharge_ = 0;
			headerAlloc_[clusterPos].setCharge(0);
		}

		assert(i==0 || clusterAlloc_[clusterPos].clusterMOverZ_ >= clusterAlloc_[clusterPos-1].clusterMOverZ_);

		// add the cluster index to the peak lists
		for (size_t j=0; j<NUM_PEAKS_FOR_HEURISTIC; j++)
		{
			const size_t idx = clusterAlloc_[clusterPos].topPeakIdxs_[j];
			if (idx>0)
			{
				if (idx>=newClusterIdxsLists_.size())
					newClusterIdxsLists_.resize(idx + 1000);

				if (! newClusterIdxsLists_[idx])
					newClusterIdxsLists_[idx] = clusterIdxAllocation_.allocateVector();

				newClusterIdxsLists_[idx]->addElement(clusterIdx);
			}
		}
	}

	assert(numPeaksReadFromDat == newSpectraStats.back().peakWritePos + newSpectraStats.back().numPeaks);

	nextClusterPos_ += newSpectraStats.size();
	nextPeakPos_    += numPeaksReadFromDat;

	assert( nextClusterPos_ <= clusterAlloc_.size() );
	assert( nextPeakPos_    <= peakAlloc_.size() );

	newSpectraStats.clear();

	return numSpectraReadFromDat;
}


// Sets the SQS value needed for an unpaired singleton in order for it to be written in the output.
// This threshold is set according to the spectrum density
void MsClusterDataStorage::setMinSqsForSingleton(float sqsThreshold)
{
	minSqsForSingleton_ = sqsThreshold;
	if (sqsThreshold <= 0.0)
		return;

	assert(datFileManager_.getMaxMOverZ()>0.0 && datFileManager_.getMinMOverZ()>0.0);
	const double spectrumDensity = (datFileManager_.getTotalSpectraToCluster() * 1000.0) / 
		(datFileManager_.getMaxMOverZ() - datFileManager_.getMinMOverZ());

	assert(spectrumDensity>= 0.0);
	double logDiff = log10(spectrumDensity) - log10(50000.0);
	if (logDiff<=0.0)
		return;
	
	minSqsForSingleton_ = (1.25 + 1.15*logDiff) * sqsThreshold;
	
	// don't let it get too high...
	if (minSqsForSingleton_ > 0.45)
		minSqsForSingleton_ = 0.45;
}


void MsClusterDataStorage::writeClustersToOutput(size_t n)
{
	if (verboseLevel_>0)
	{
		cout << "Writing clusters for " << n << " singletons";
		cout << "...";
	}

	size_t numSkipped = 0;
	size_t numLowSqs    = 0;
	size_t numClustersWritten =0;

	if (! params_->gotMergeArchives)
	{
		numClustersWritten=writeRegularClustersToOutput(n,  numSkipped, numLowSqs);
	}
	else
		numClustersWritten=writeMergedArchiveClustersToOutput(n,  numSkipped, numLowSqs);

	totalSpectraWritten_ += n;

	if (verboseLevel_>0)
	{
		cout << " written into " << numClustersWritten << " clusters ";
		if (numLowSqs)
			cout << " [ removed w/low SQS " << numLowSqs <<  " ]";
		cout << "   Cumulative " << totalClustersWritten_ << " : " << totalSpectraWritten_;
		cout.flush();
	}

	// close all open output files that have an index that is -5 below the lowest mass outputted
	// there is no chance that they will have anything written to them
	clusterOutputter_.closeAllWithMassBelow(clusterAlloc_[0].getClusterMOverZ()-5.0);
	if (verboseLevel_>0)
	{
		cout << " ..done." << endl << endl;
		cout.flush();
	}
}




// writes the first n clusters in clusterAlloc_ to the appropriate files
size_t MsClusterDataStorage::writeRegularClustersToOutput(size_t n, size_t& numSkipped, size_t& numLowSqs)
{
	ostringstream ossHeader;
	//ostringstream ossMembers;
	SingleSpectrumHeader clusterHeader;
	SingleSpectrumHeader* header;
	ostringstream clustInfo;

	numSkipped  = 0;
	numLowSqs   = 0;
	size_t numWritten=0;

	// title csv datasetIdx, spectraFileIndexInList, scanNumber, peakArea, peakId, rt, mz, simToConsensus, pValuecharge, peptideStr
	for (size_t i=0; i<n; i++)
	{
		Cluster& cluster = clusterAlloc_[i];
        if (! cluster.getIndInPlay())
		{
			numSkipped++;
			continue;
		}

		if (cluster.getClusterSize()<=1 && minSqsForSingleton_>0.0)
		{
			const float sqs = cluster.getHeader()->getSqs();
			if (sqs>0.0 && sqs<minSqsForSingleton_)
			{
				numLowSqs++;
				cluster.setIndInPlay(false);
				continue;
			}
		}

		assert( cluster.getNumSingletonsIncludingSelf() == cluster.getClusterSize());
		assert( cluster.getHeader() );

		const char* existingTitle = NULL;
		string maximalPeptide = std::string(); // the peptide with maximal counts (if exists)
//		const size_t numSingletonLines = makeClusterMembersStringAndSelectTitleAndPeptide(
//			cluster.getClusterIdx(), ossMembers, existingTitle, maximalPeptide);
		const vector<string> singletonLines = makeClusterMembersStringAndSelectTitleAndPeptide_CSV(
		        cluster.getClusterIdx(), existingTitle, maximalPeptide);

		// set charge according to maximal peptide
		if (maximalPeptide.length()>0)
		{
			Peptide pep;
			pep.parseFromString(config_, maximalPeptide);
			int newCharge = config_->determine_charge(pep, cluster.clusterMOverZ_);
			if (newCharge>0)
				cluster.clusterCharge_ = newCharge;
		}


		clustInfo.str(std::string());
		// make title (look for existing title)
		ossHeader.str(std::string());

		// NP3 add header info
		//ossHeader << "ID/Run" << "\tFile\t" << "N/Scan" << "\t" << "RTmean"  << "\t" << "Mz" << "\t"
		//          << "Chg/Sim" << "\tRTmin/P" << "\tRTmax" << "\tClustInt" << endl;
		if (params_->gotUseInputTitles && existingTitle)
		{
            clustInfo << existingTitle << "," << runningOutputIdx_;

			//ossHeader << existingTitle;
		}
		else
            clustInfo << nameStr_ << "," << runningOutputIdx_;
			//ossHeader << nameStr_ << "." << runningOutputIdx_++;

		// make consensus
		if (cluster.getClusterSize()>1 &&
			cluster.getSingletonVector() &&
			cluster.getSingletonVector()->getSize()>0)
		{
			makeConsensus(&cluster, &clusterHeader);
			header = &clusterHeader;
			assert( cluster.fullSanityCheck(cluster.getClusterSize()) );
		}
		else
		{
			header =const_cast<SingleSpectrumHeader*>(cluster.header_);
		}

		// NP3 GOT OUTPUT MGF skipping low number of peaks < minNumPeaks
		// TODO save removed scans due to low minimum number of peaks
		if (cluster.getNumPeaks() < params_->minNumPeaks)
		{
			numSkipped++;
			continue;
		}

		header->setOriginalNumPeaks(cluster.getNumPeaks());
		//header->setTitle(ossHeader.str());
		// NP3 new title
        header->setTitle(nameStr_ + "." + to_string(runningOutputIdx_));

		runningOutputIdx_++;

		// make rest of cluster header line after the consensus
		//ossHeader << "\t" << cluster.getNumSingletonsIncludingSelf() << "\t"
		//		  << cluster.clusterMOverZ_ << "\t" << cluster.clusterCharge_;
		// NP3 GOT output more cluster info to .clust
        clustInfo << "," << cluster.getNumSingletonsIncludingSelf() << "," << std::to_string(cluster.getClusterRT()) << "," <<
        		to_string(cluster.clusterMOverZ_) << "," << cluster.clusterCharge_ << "," << to_string(cluster.getClusterRTMin()) << "," <<
        		to_string(cluster.getClusterRTMax()) << "," << to_string(cluster.getClusterTotalPrecursorIntensity());

//		ossHeader << "\t" << cluster.getNumSingletonsIncludingSelf() << "\t"
//        				  << cluster.getClusterRT()  << "\t"
//        				  << cluster.clusterMOverZ_ << "\t" << cluster.clusterCharge_
//        				  << "\t" << cluster.getClusterRTMin() << "\t" << cluster.getClusterRTMax() << "\t"
//        				  << cluster.getClusterTotalPrecursorIntensity();
		if (cluster.getHeader()->getPeptideStr().length()>1)
		{
			assert( cluster.clusterCharge_ > 0 );
            clustInfo << "," << cluster.getHeader()->getPeptideStr();
			//ossHeader << "\t" << cluster.getHeader()->getPeptideStr();
		}
		else if (maximalPeptide.length()>0 && cluster.clusterCharge_>0)
		{
            clustInfo << "," << maximalPeptide;
			//ossHeader << "\t" << maximalPeptide << endl;
		} else {
            clustInfo << ",";
		}
		//ossHeader << endl;
		// NP3 add header info of clust members
		//ossHeader << "Process" << "\t" << "FileIdx" << "\t" << "Scan"  << "\t" << "RTmean" << "\t" << "Mz"
		//		  << "\tSim" << "\tP-value" << endl;
		//ossHeader << ossMembers.str();

		//NP3 for each member concate clust info into the same row
		for (int j = 0; j < singletonLines.size(); j++) {
            ossHeader << clustInfo.str() << "," << singletonLines[j] << endl;
		}

		assert(header->getOriginalNumPeaks() > 0);
		cluster.setHeader(header);
		clusterOutputter_.outputCluster(cluster, ossHeader.str(), true, (header->getPeptideStr().length() == 0) );
		totalClustersWritten_++;
		numWritten++;
		cluster.setIndInPlay(false);
	}

	// take care of memory that is specially allocated
	for (size_t i=0; i<n; i++)
	{
		if (clusterAlloc_[i].singletonIdxVector_)
			clusterAlloc_[i].singletonIdxVector_->relinquish();
	
		if (clusterAlloc_[i].singletonDistancePeaks_)
		{
			delete clusterAlloc_[i].singletonDistancePeaks_;
			clusterAlloc_[i].singletonDistancePeaks_ =0;
		}
	}

	return numWritten;
}


// obtains the cluster record from the relevant clust files
bool MsClusterDataStorage::getClusterEntry(const Cluster& cluster, ClusterClustEntry& clustEntry)
{
	const string& searchTitle = cluster.getHeader()->getTitle();
	// obtain the clust records
	map<string, ClusterClustEntry>::const_iterator it = clustEntries_.find(searchTitle);
	if (it == clustEntries_.end())
	{
		const int datasetIdx = cluster.getHeader()->getDatasetIndex();
		const int fileIdx    = cluster.getHeader()->getSpectraFileIndexInList();
		const string& clusterTitle = cluster.getHeader()->getTitle();

		const string& datPath = datFileManager_.getDatFile(fileIdx).getPath();


	#if defined(WIN32) || defined(WIN64)
		const size_t lastSlash = datPath.find_last_of('\\');
		const string clustPath = datPath.substr(0,lastSlash-3) + "clust\\" + datPath.substr(lastSlash+1,datPath.length()- lastSlash - 4) + "clust";
	#else
		const size_t lastSlash = datPath.find_last_of('/');
		const string clustPath = datPath.substr(0,lastSlash-3) + "clust/" + datPath.substr(lastSlash+1,datPath.length()- lastSlash - 4) + "clust";
	#endif

		// read all entries from clust file into clustEntries_

		clustFile_.read(clustPath.c_str(),datasetIdx, fileIdx);
		char* p;
		string title;
	//	int c=0;
	//	cout << endl;
		while (p=clustFile_.getNextEntryPointer(&title))
		{
			ClusterClustEntry e;
			clustFile_.parseEntry(p, e);
			clustEntries_[title]=e;
		//	if (c++<20)
		//		cout << title << "\t" << e.clusterIdx << endl;
		}

		it = clustEntries_.find(searchTitle);
		if (it == clustEntries_.end())
		{
			cout << "Looking for " << searchTitle << " in " << clustPath.c_str() << endl;
			cout << "dataset idx: " << datasetIdx << "\tfileIdx: " << fileIdx << endl;
			error("Could not find entry in clust file for ",searchTitle.c_str());
		}

	}

	
	clustEntry = it->second;
	clustEntries_.erase(searchTitle);
	return true;
}


// writes the first n clusters in clusterAlloc_ to the appropriate files
// this function operates differently since in needs to create the cluster membership records
// by reading the existing ".clust" files (since clusters are actually being clustered)
size_t MsClusterDataStorage::writeMergedArchiveClustersToOutput(size_t n, size_t& numSkipped, size_t& numLowSqs)
{
	ostringstream ossHeader;
	ostringstream ossMembers;
	SingleSpectrumHeader clusterHeader;
	SingleSpectrumHeader* header;

	numSkipped  = 0;
	numLowSqs   = 0;
	size_t numWritten=0;

	for (size_t i=0; i<n; i++)
	{
		Cluster& cluster = clusterAlloc_[i];
		if (! cluster.getIndInPlay())
		{
			numSkipped++;
			continue;
		}
		assert( cluster.getNumSingletonsIncludingSelf() <= cluster.getClusterSize());
		assert( cluster.getHeader() );

		string maximalPeptide = std::string(); // the peptide with maximal counts (if exists)
		string existingTitle  = std::string();

		const size_t numSingletonLines = makeClustMemberStringForMergedClusters(cluster, ossMembers, existingTitle, maximalPeptide);
		assert(existingTitle.size()>0);

		// set charge according to maximal peptide
		if (maximalPeptide.length()>0)
		{
			Peptide pep;
			pep.parseFromString(config_, maximalPeptide);
			int newCharge = config_->determine_charge(pep, cluster.clusterMOverZ_);
			if (newCharge>0)
				cluster.clusterCharge_ = newCharge;
		}


		// make title (look for existing title)
		ossHeader.str(std::string());
		ossHeader << existingTitle;
			
		// make consensus
		if (cluster.getClusterSize()>1 && 
			cluster.getSingletonVector() && 
			cluster.getSingletonVector()->getSize()>0)
		{
			makeConsensus(&cluster, &clusterHeader);
			header = &clusterHeader;
			assert( cluster.fullSanityCheck(cluster.getClusterSize()) );
		}
		else
		{
			header =const_cast<SingleSpectrumHeader*>(cluster.header_);
		}

		header->setOriginalNumPeaks(cluster.getNumPeaks());
		header->setTitle(ossHeader.str());

		// make rest of cluster header line after the consensus
		ossHeader << "\t" << numSingletonLines << "\t" << cluster.clusterMOverZ_ << "\t" << cluster.clusterCharge_;
		if (cluster.getHeader()->getPeptideStr().length()>1)
		{
			assert( cluster.clusterCharge_ > 0 );
			ossHeader << "\t" << cluster.getHeader()->getPeptideStr();
		}
		ossHeader << endl;

		ossHeader << ossMembers.str();

		assert(header->getOriginalNumPeaks() > 0);
		cluster.setHeader(header);
		clusterOutputter_.outputCluster(cluster, ossHeader.str(), true, (header->getPeptideStr().length() == 0) );
		totalClustersWritten_++;
		numWritten++;
		cluster.setIndInPlay(false);
	}

	// take care of memory that is specially allocated
	for (size_t i=0; i<n; i++)
	{
		if (clusterAlloc_[i].singletonIdxVector_)
			clusterAlloc_[i].singletonIdxVector_->relinquish();
	
		if (clusterAlloc_[i].singletonDistancePeaks_)
		{
			delete clusterAlloc_[i].singletonDistancePeaks_;
			clusterAlloc_[i].singletonDistancePeaks_ =0;
		}
	}
	return numWritten;
}





// shifts all relevant memory left (pushing away the clusters on the lower m/z range)
// the sift is done to all storage vectors (peaks, headers, clusters, lists, etc.)
// this (expensive) action is to be taken before a new m/z slice of spectra is added for 
// clustering. Every cluster that is removed is outputed (to DAT and/or MGF)
void  MsClusterDataStorage::shiftStorage(size_t posToMoveToZero, mass_t windowSize)
{
	if (posToMoveToZero==0)
		return;

	// first write clusters that are removed to the output
	writeClustersToOutput(posToMoveToZero);
	
	const size_t numToShift = nextClusterPos_ - posToMoveToZero; // number of items that get shifted
	if (numToShift > 0)
	{
		// this is the number of peaks that get shifted
		const size_t peakStart = (clusterAlloc_[posToMoveToZero].getPeaks() - &peakAlloc_[0]);
		shuntInVector(peakAlloc_, peakStart, nextPeakPos_ - peakStart);
		nextPeakPos_ -= peakStart;

		// move the headers and clusters starting at position n.. n+numToShift-1
		// to positions 0..numToShift-1
		shuntInVector(headerAlloc_, posToMoveToZero, numToShift); 
		shuntInVector(clusterAlloc_, posToMoveToZero, numToShift);
	
		// the singleton idx vector is still kept by the original cluster however, after
		// shunting there are now two pointers to it. This is not good.
		for (size_t i=numToShift; i<nextClusterPos_; i++)
		{
			headerAlloc_[i].setMOverZ(-1.0);
			clusterAlloc_[i].singletonIdxVector_ = NULL; 
			clusterAlloc_[i].clusterSize_ = 0;
			clusterAlloc_[i].indInPlay_   = false;
			clusterAlloc_[i].singletonDistancePeaks_ = NULL;
		}
		nextClusterPos_ = numToShift;
		UpdateClusterIdxsLists(clusterAlloc_[0].getClusterIdx());

		// update peaks_ pointers and header pointers in all shifted clsuter
		for (size_t i=0; i<numToShift; i++)
		{
			clusterAlloc_[i].peaks_ -= peakStart;
			clusterAlloc_[i].setHeader(&headerAlloc_[i]);
		}
	}
	else // simply reset pointers
	{
		for (size_t i=0; i<nextClusterPos_; i++)
		{
			clusterAlloc_[i].singletonIdxVector_ = NULL; 
			clusterAlloc_[i].singletonDistancePeaks_ = NULL;
			clusterAlloc_[i].clusterSize_ = 0;
			clusterAlloc_[i].indInPlay_   = false;
		}
		nextClusterPos_ = 0;
		nextPeakPos_ = 0;
	}
}

// removes all cluster idxs that have a lower
 // index than the one given (these were all shifted out)
void MsClusterDataStorage::UpdateClusterIdxsLists(clusterIdx_t minClusterIdx)
{
	for (size_t i=0; i<newClusterIdxsLists_.size(); i++)
	{
		if (! newClusterIdxsLists_[i] || newClusterIdxsLists_[i]->getSize() == 0)
			continue;
		const size_t size = newClusterIdxsLists_[i]->getSize();
		const clusterIdx_t* idxs = newClusterIdxsLists_[i]->getElementsPtr();

		const size_t pos = newClusterIdxsLists_[i]->getFirstPositionGreaterOrEqual(minClusterIdx);
		if (pos==0)
			continue;

		if (pos>=size)
		{
			newClusterIdxsLists_[i]->relinquish();
			newClusterIdxsLists_[i] = NULL;
			continue;
		}
		newClusterIdxsLists_[i]->shunt(pos);
	}
}

// This struct is used to sort the lines in the cluster membership output
struct ClusterMemberEntry {

	ClusterMemberEntry() : singleton(0), datasetIdx(-1), fileIndexInList(-1), scanNumber(-1), 
							MOverZ(0), charge(0), title(0), peptide(0) {}

	bool operator< (const ClusterMemberEntry& rhs) const
	{
		if (datasetIdx >  rhs.datasetIdx)
			return false;
		if (datasetIdx == rhs.datasetIdx && fileIndexInList >  rhs.fileIndexInList)
			return false;
		if (datasetIdx == rhs.datasetIdx && fileIndexInList == rhs.fileIndexInList && scanNumber >= rhs.scanNumber)
			return false;

		return true;
	}

	const Cluster* singleton;
	int datasetIdx;
	int fileIndexInList;
	int scanNumber;
	mass_t MOverZ;
	short charge;
	const char* title;
	const char* peptide;
};

// NP3 outputting clust file as a CSV table
// make cluster member string. Makes the text entry that should appear in the clust file
// assumes that any singleton with a title at this stage is one from previous generations
// so it parses the title to get the generation information and reports the title
// returns number of lines in the oss
vector<string> MsClusterDataStorage::makeClusterMembersStringAndSelectTitleAndPeptide_CSV(clusterIdx_t idx, const char*& existingTitle, string& maximalPeptide)
{
    vector<string> members_rows = vector<string>();
	ostringstream member_str;

    const clusterIdx_t firstClusterIdx = clusterAlloc_[0].getClusterIdx();
    assert(idx >= firstClusterIdx);

    const Cluster& cluster = clusterAlloc_[idx - firstClusterIdx];
    assert(idx == cluster.getClusterIdx());

    const clusterIdx_t* singltonIdxs = (cluster.singletonIdxVector_ ?
                                        cluster.singletonIdxVector_->getElementsPtr() : 0);
    const size_t numSingletons       = (cluster.singletonIdxVector_ ?
                                        cluster.singletonIdxVector_->getSize() : 0);



    // special (easy) case if there are no additional singletons (no sorting needed)
    if (numSingletons == 0)
    {
		member_str.str(std::string());
        const SingleSpectrumHeader* header = cluster.getHeader();
        const int headerDatasetIndex = header->getDatasetIndex();

		member_str << (headerDatasetIndex >= 0 ? headerDatasetIndex : datasetIdx_)  << "," << header->getSpectraFileIndexInList()
		   		   << "," << std::to_string(header->getScanNumber()) << "," << std::to_string(header->getPeakArea()) << "," << header->getPeakId();
//        oss << (headerDatasetIndex >= 0 ? headerDatasetIndex : datasetIdx_ )<< "\t"
//            << header->getSpectraFileIndexInList() << "\t"
//            << header->getScanNumber() << "\t" << header->getPeakArea() << "\t";

        if (params_->gotCorrectPM)
        {
            const mass_t deltaMz = fabs(header->getOriginalPmWith19() - header->getMOverZ());
            assert( deltaMz<15.0);
            SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
            nonConstHeader->setMOverZ(header->getOriginalPmWith19());
        }

        member_str << "," <<  std::to_string(header->getRetentionTime()) << "," <<  std::to_string(header->getMOverZ())
				   << ",1.0,0," << cluster.clusterCharge_;
        //oss << header->getRetentionTime() << "\t" << header->getMOverZ() << "\t1.0\t0\t" << cluster.clusterCharge_ ;
        if (header->getPeptideStr().length()>1)
            //oss << "\t" << header->getPeptideStr();
            member_str << "," <<  header->getPeptideStr();

        members_rows.push_back(member_str.str());
        //oss << endl << endl;
        return members_rows;
    }

    // there are multiple singletons, they should be sorted according to generation/file/scan
    vector<ClusterMemberEntry> members(numSingletons+1);
    map<string,int> peptideCounts;

    // first collect pointers to the singletons (including the cluster's signelton)
    for (size_t i=0; i<numSingletons; i++)
    {
        const clusterIdx_t singletonIdx = singltonIdxs[i];
        assert(singletonIdx >= firstClusterIdx);
        members[i].singleton = &clusterAlloc_[singletonIdx - firstClusterIdx];
    }
    members[numSingletons].singleton = &cluster;

    for (size_t i=0; i<=numSingletons; i++)
    {
        const SingleSpectrumHeader* singletonHeader = members[i].singleton->getHeader();
        const int headerDatasetIndex =singletonHeader->getDatasetIndex();
        members[i].datasetIdx	   = (headerDatasetIndex >=0 ? headerDatasetIndex : datasetIdx_);
        members[i].fileIndexInList = singletonHeader->getSpectraFileIndexInList();
        members[i].scanNumber	   = singletonHeader->getScanNumber();
        members[i].MOverZ		   = singletonHeader->getMOverZ();
        members[i].charge		   = singletonHeader->getCharge();
        members[i].title		   = (singletonHeader->getTitle().length()>1 ? singletonHeader->getTitle().c_str() : NULL);

        const string& peptideStr = singletonHeader->getPeptideStr();
        if (peptideStr.length()>=1)
        {
            members[i].peptide		   = peptideStr.c_str();
            peptideCounts[peptideStr]++;
        }

        // hack return the orginal m/z if there was correction (so data doesn't look funny)
        if (params_->gotCorrectPM)
        {
            const mass_t deltaMz = fabs(singletonHeader->getOriginalPmWith19() - singletonHeader->getMOverZ());
            assert( deltaMz<15.0);
            members[i].MOverZ = singletonHeader->getOriginalPmWith19();
        }
    }

    sort(members.begin(), members.end());

    // create the string
    for (size_t i=0; i<members.size(); i++)
    {
		member_str.str(std::string());
        // NP3 GOT more cluster info in the .clust file: oss << members[i].datasetIdx << "\t" << members[i].fileIndexInList << "\t" << members[i].scanNumber << "\t" << members[i].MOverZ;
        member_str << members[i].datasetIdx << "," << members[i].fileIndexInList << "," << members[i].scanNumber << ","
        		   << std::to_string(members[i].singleton->getHeader()->getPeakArea()) << "," << members[i].singleton->getHeader()->getPeakId() << ","
        		   << members[i].singleton->getClusterRT() << "," << members[i].MOverZ;
        //        oss << members[i].datasetIdx << "\t" << members[i].fileIndexInList << "\t" << members[i].scanNumber << "\t"
//            << members[i].singleton->getHeader()->getPeakArea() << "\t" << members[i].singleton->getClusterRT() << "\t"  << members[i].MOverZ;
//        // output similarity to consensus and the p-value
        // NP3 GOT sim computation
        //const float similarity = computeSimilarity(cluster.getDistancePeaks(), members[i].singleton->getSingletonPeaks(), Cluster::getPeakIndexToleranceAsInt());
        const float similarity = computeSimilarity(cluster.getDistancePeaks(),
                                                   members[i].singleton->getSingletonPeaks(), Cluster::getPeakIndexToleranceAsInt());
        const float pvalue	   = simModel_->computePValue(cluster.getDistancePeaks()->numPeaks, members[i].singleton->getNumSimilarityPairs(), similarity);

        member_str << "," << setprecision(2) << fixed << similarity;
        //oss << "\t" << setprecision(2) << fixed << similarity << "\t";
        if (pvalue>0.0)
        {
            member_str << "," << setprecision(3) << scientific << pvalue << fixed;
            //oss << setprecision(3) << scientific << pvalue << fixed;
        }
        else
            member_str << ",0";
            //oss << "0";

        member_str << "," <<  members[i].charge;
        //oss << "\t" << members[i].charge;
        if (members[i].peptide)
            member_str << "," << members[i].peptide;
            //oss << "\t" << members[i].peptide;
        //oss << endl;
        members_rows.push_back(member_str.str());
    }
    //oss << endl;

    // if there is a title from previous generations, the earliest one will be here
    existingTitle = members[0].title;

    // choose the peptide that has the maximal counts as the cluster's peptide
    if (peptideCounts.size()> 0)
    {
        int maxCount=0;
        map<string,int>::const_iterator it;
        for (it = peptideCounts.begin(); it != peptideCounts.end(); it++)
        {
            if (it->second>maxCount)
            {
                maxCount = it->second;
                maximalPeptide = it->first;
            }
        }
    }
    return (members_rows);
}


// make cluster member string. Makes the text entry that should appear in the clust file
// assumes that any singleton with a title at this stage is one from previous generations
// so it parses the title to get the generation information and reports the title
// returns number of lines in the oss
size_t MsClusterDataStorage::makeClusterMembersStringAndSelectTitleAndPeptide(clusterIdx_t idx, 
					ostringstream& oss, const char*& existingTitle, string& maximalPeptide) const
{
	oss.str(std::string());
	const clusterIdx_t firstClusterIdx = clusterAlloc_[0].getClusterIdx();
	assert(idx >= firstClusterIdx);

	const Cluster& cluster = clusterAlloc_[idx - firstClusterIdx];
	assert(idx == cluster.getClusterIdx());

	const clusterIdx_t* singltonIdxs = (cluster.singletonIdxVector_ ? 
		cluster.singletonIdxVector_->getElementsPtr() : 0);
	const size_t numSingletons       = (cluster.singletonIdxVector_ ? 
		cluster.singletonIdxVector_->getSize() : 0);

	

	// special (easy) case if there are no additional singletons (no sorting needed)
	if (numSingletons == 0)
	{
		const SingleSpectrumHeader* header = cluster.getHeader();
		const int headerDatasetIndex = header->getDatasetIndex();
		
		oss << (headerDatasetIndex >= 0 ? headerDatasetIndex : datasetIdx_ )<< "\t" 
			<< header->getSpectraFileIndexInList() << "\t"
			<< std::to_string(header->getScanNumber()) << "\t" << std::to_string(header->getPeakArea()) << "\t";
		
		if (params_->gotCorrectPM)
		{
			const mass_t deltaMz = fabs(header->getOriginalPmWith19() - header->getMOverZ());
			assert( deltaMz<15.0);
			SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
			nonConstHeader->setMOverZ(header->getOriginalPmWith19());
		}

		oss << header->getRetentionTime() << "\t" << header->getMOverZ() << "\t1.0\t0\t" << cluster.clusterCharge_ ;
		if (header->getPeptideStr().length()>1)
			oss << "\t" << header->getPeptideStr();

		oss << endl << endl;
		return 1;
	}

	// there are multiple singletons, they should be sorted according to generation/file/scan
	vector<ClusterMemberEntry> members(numSingletons+1);
	map<string,int> peptideCounts;

	// first collect pointers to the singletons (including the cluster's signelton)
	for (size_t i=0; i<numSingletons; i++)
	{
		const clusterIdx_t singletonIdx = singltonIdxs[i];
		assert(singletonIdx >= firstClusterIdx);
		members[i].singleton = &clusterAlloc_[singletonIdx - firstClusterIdx];
	}
	members[numSingletons].singleton = &cluster;

	for (size_t i=0; i<=numSingletons; i++)
	{
		const SingleSpectrumHeader* singletonHeader = members[i].singleton->getHeader();
		const int headerDatasetIndex =singletonHeader->getDatasetIndex();
		members[i].datasetIdx	   = (headerDatasetIndex >=0 ? headerDatasetIndex : datasetIdx_);
		members[i].fileIndexInList = singletonHeader->getSpectraFileIndexInList();
		members[i].scanNumber	   = singletonHeader->getScanNumber();
		members[i].MOverZ		   = singletonHeader->getMOverZ();
		members[i].charge		   = singletonHeader->getCharge();
		members[i].title		   = (singletonHeader->getTitle().length()>1 ? singletonHeader->getTitle().c_str() : NULL);

		const string& peptideStr = singletonHeader->getPeptideStr();
		if (peptideStr.length()>=1)
		{
			members[i].peptide		   = peptideStr.c_str();
			peptideCounts[peptideStr]++;
		}
	
		// hack return the orginal m/z if there was correction (so data doesn't look funny)
		if (params_->gotCorrectPM)
		{
			const mass_t deltaMz = fabs(singletonHeader->getOriginalPmWith19() - singletonHeader->getMOverZ());
			assert( deltaMz<15.0);
			members[i].MOverZ = singletonHeader->getOriginalPmWith19();
		}
	}
	
	sort(members.begin(), members.end());

	// create the string
	for (size_t i=0; i<members.size(); i++)
	{
		// NP3 GOT more cluster info in the .clust file: oss << members[i].datasetIdx << "\t" << members[i].fileIndexInList << "\t" << members[i].scanNumber << "\t" << members[i].MOverZ;
		oss << members[i].datasetIdx << "\t" << members[i].fileIndexInList << "\t" << std::to_string(members[i].scanNumber) << "\t"
				<< std::to_string(members[i].singleton->getHeader()->getPeakArea()) << "\t" <<
				std::to_string(members[i].singleton->getClusterRT()) << "\t"  << std::to_string(members[i].MOverZ);
		// output similarity to consensus and the p-value
		// NP3 GOT sim computation
		//const float similarity = computeSimilarity(cluster.getDistancePeaks(), members[i].singleton->getSingletonPeaks(), Cluster::getPeakIndexToleranceAsInt());
		const float similarity = computeSimilarity(cluster.getDistancePeaks(),
		                                           members[i].singleton->getSingletonPeaks(), Cluster::getPeakIndexToleranceAsInt());
		const float pvalue	   = simModel_->computePValue(cluster.getDistancePeaks()->numPeaks, members[i].singleton->getNumSimilarityPairs(), similarity);
		oss << "\t" << setprecision(2) << fixed << similarity << "\t";
		if (pvalue>0.0)
		{
			oss << setprecision(3) << scientific << pvalue << fixed;
		}
		else
			oss << "0";

		oss << "\t" << members[i].charge;
		if (members[i].peptide)
			oss << "\t" << members[i].peptide;
		oss << endl;
	}
	oss << endl;

	// if there is a title from previous generations, the earliest one will be here
	existingTitle = members[0].title;

	// choose the peptide that has the maximal counts as the cluster's peptide
	if (peptideCounts.size()> 0)
	{
		int maxCount=0;
		map<string,int>::const_iterator it;
		for (it = peptideCounts.begin(); it != peptideCounts.end(); it++)
		{
			if (it->second>maxCount)
			{
				maxCount = it->second;
				maximalPeptide = it->first;
			}
		}
	}
	return (members.size());
}



size_t MsClusterDataStorage::makeClustMemberStringForMergedClusters(const Cluster& cluster, 
											ostringstream& oss, string& existingTitle, string& maximalPeptide, bool verbose)
{
	const clusterIdx_t idx = cluster.getClusterIdx();
	oss.str(std::string());
	const clusterIdx_t firstClusterIdx = clusterAlloc_[0].getClusterIdx();
	assert(idx >= firstClusterIdx);

	const clusterIdx_t* singltonIdxs = (cluster.singletonIdxVector_ ? cluster.singletonIdxVector_->getElementsPtr() : 0);
	const size_t numSingletons       = (cluster.singletonIdxVector_ ? cluster.singletonIdxVector_->getSize() : 0);

	vector<clusterIdx_t> idxsToAdd;
	idxsToAdd.push_back(idx);
	for (size_t i=0; i<numSingletons; i++)
		idxsToAdd.push_back(singltonIdxs[i]);

	int lowestD   = cluster.getHeader()->getDatasetIndex();
	existingTitle = cluster.getHeader()->getTitle();
	map<string,int> peptideCounts;
	vector<SingletonClustEntry> allEntries;

	if (verbose)
	{
		cout << endl << cluster.getHeader()->getTitle() << endl;
		cout << "Idxs " << idxsToAdd.size() << " : ";
		for (size_t i=0; i<idxsToAdd.size(); i++)
			cout << "\t" << idxsToAdd[i];
		cout << endl;
	}
	for (size_t i=0; i<idxsToAdd.size(); i++)
	{
		const Cluster& singleton = *getCluster(idxsToAdd[i]);
		ClusterClustEntry entry;
		getClusterEntry(singleton, entry);
		const int clusterSource = singleton.getHeader()->getArchiveSource();
		assert(clusterSource ==1 || clusterSource ==2);
		const int numToAddToDatasetIdx = (clusterSource == 2 ? archive1_.getNumDatasets() : 0);
		const size_t numEntries = entry.entries.size();
		for (size_t j=0; j<numEntries; j++)
		{
			// modify the entries if they came from clust2 add to the dataset indexes
			entry.entries[j].datasetIdx += numToAddToDatasetIdx;
			allEntries.push_back(entry.entries[j]);
		}

		// select earliest title
		if (singleton.getHeader()->getDatasetIndex()<lowestD)
		{
			lowestD = singleton.getHeader()->getDatasetIndex();
			existingTitle = singleton.getHeader()->getTitle();
		}


		// add peptide counts
		if (singleton.getHeader()->getPeptideStr().length()>0)
			peptideCounts[singleton.getHeader()->getPeptideStr()]+=numEntries;
	}
	sort(allEntries.begin(),allEntries.end());

	for (size_t i=0; i<allEntries.size(); i++)
	{
		SingletonClustEntry& e=allEntries[i];
		oss << e.datasetIdx << "\t" << e.fileIdx << "\t" << e.scanNumber << "\t" << e.mz << "\t" << e.similarityToConsensus << "\t" << e.pValue;
		oss << "\t" << e.charge;
		if (e.peptide.length()>0)
			oss << "\t" << e.peptide;
		oss << endl;
	}
	oss << endl;

	if (verbose)
	{
		cout << "Entry: " << oss.str();
	}

	// choose the peptide that has the maximal counts as the cluster's peptide
	if (peptideCounts.size()> 0)
	{
		int maxCount=0;
		map<string,int>::const_iterator it;
		for (it = peptideCounts.begin(); it != peptideCounts.end(); it++)
		{
			if (it->second>maxCount)
			{
				maxCount = it->second;
				maximalPeptide = it->first;
			}
		}
	}
	return (allEntries.size());
}



size_t calcClusterSizeBin(unsigned int clusterSize)
{
	size_t i;
	for (i=0; i<NUM_SIZES_FOR_REDOING_DISTANCE_PEAKS; i++)
		if (clusterSize<SIZES_FOR_REDOING_DISTANCE_PEAKS[i])
			break;
	return i;
}



// adds addedIdx to the clsuter mainIdx.
// performs all the needed actions including updating the info in the added clusters singletons.
void   MsClusterDataStorage::joinClusters(Cluster* mainCluster, Cluster* addedCluster, int context)
{
	assert(mainCluster && addedCluster && mainCluster->getClusterIdx() < addedCluster->getClusterIdx());
	assert(mainCluster->getIndInPlay() && addedCluster->getIndInPlay());

	const clusterIdx_t mainClusterIdx = mainCluster->clusterIdx_;
	const unsigned int sizeBefore = mainCluster->clusterSize_;
	const size_t addedPosition = getCluserPosition(addedCluster->clusterIdx_);

	assert(addedPosition < nextClusterPos_);

	/*if (abs(mainCluster->getClusterMOverZ() - addedCluster->getClusterMOverZ()) >= MsClusterDataStorage::params_->mzWindow)
	{
		cout << "join mass diff " << mainCluster->getClusterMOverZ() - addedCluster->getClusterMOverZ() << " sim" << endl;
	}*/
	// add other to this cluster's singletons list_
	if (! mainCluster->singletonIdxVector_)
		mainCluster->singletonIdxVector_ = clusterIdxAllocation_.allocateVector();

	if (! mainCluster->singletonDistancePeaks_)
		mainCluster->backupDistancePeaks();

	mainCluster->singletonIdxVector_->addElement(addedCluster->clusterIdx_);
	mainCluster->clusterSize_				 += addedCluster->clusterSize_;
	mainCluster->totalNonPrecursorIntensity_ += addedCluster->totalNonPrecursorIntensity_;

    // NP3 GOT retention time RT mean computation on joining
    mainCluster->clusterTotalPrecursorIntensity_ += addedCluster->clusterTotalPrecursorIntensity_;
    mainCluster->clusterTotalRT_ += addedCluster->clusterTotalRT_;
    mainCluster->clusterRT_ 	 = mainCluster->clusterTotalRT_ / mainCluster->clusterTotalPrecursorIntensity_;
    // NP3 GOT rt width mean computation on joining
	if ((mainCluster->clusterTotalRTMin_ == 0.0 && mainCluster->clusterTotalRTMax_ == 1000000.0) ||
		(addedCluster->clusterTotalRTMin_ == 0.0 && addedCluster->clusterTotalRTMax_ == 1000000.0)) {
		// keep baseline rt range, for blanks
		mainCluster->clusterTotalRTMin_ = 0.0;
		mainCluster->clusterTotalRTMax_ = 1000000.0;
		mainCluster->clusterRTMin_ = 0.0;
		mainCluster->clusterRTMax_ = 1000000.0;
	} else {
		// compute the average
        mainCluster->clusterTotalRTMin_ += addedCluster->clusterTotalRTMin_;
        mainCluster->clusterTotalRTMax_ += addedCluster->clusterTotalRTMax_;
        mainCluster->clusterRTMin_ = mainCluster->clusterTotalRTMin_ / mainCluster->clusterTotalPrecursorIntensity_;
        mainCluster->clusterRTMax_ = mainCluster->clusterTotalRTMax_ / mainCluster->clusterTotalPrecursorIntensity_;
    }

    // NP3 GOT check if rtmean is inside the peak, if not recompute it as the mean of min and max rt
    if (mainCluster->clusterRT_ < mainCluster->clusterRTMin_ || mainCluster->clusterRT_  > mainCluster->clusterRTMax_)
    {
        mainCluster->clusterRT_ = (mainCluster->clusterRTMin_ + mainCluster->clusterRTMax_)/2;
    }

	// update the added cluster
	//addedCluster->indInPlay_		     = false;
	addedCluster->setIndInPlay(false);
	addedCluster->assignedClusterIdx_    = mainClusterIdx;

	// update the assigned cluster index of the other singletons
	if (addedCluster->getClusterSize()>1 && addedCluster->singletonIdxVector_)
	{
		const size_t vecSize    = addedCluster->singletonIdxVector_->getSize();
		const clusterIdx_t* vec = addedCluster->singletonIdxVector_->getElementsPtr();
		for (size_t i=0; i<vecSize; i++)
		{

			mainCluster->singletonIdxVector_->addElement(vec[i]);
			Cluster* singleton = getCluster(vec[i]);
			assert( singleton && ! singleton->getIndInPlay());
			singleton->assignedClusterIdx_ = mainClusterIdx;
		}
		addedCluster->singletonIdxVector_->relinquish();
		addedCluster->singletonIdxVector_ = NULL;
	}

	// do all other updates if there is a needed
	const size_t sizeBinBefore = calcClusterSizeBin(sizeBefore);
	const size_t sizeBinAfter  = calcClusterSizeBin(mainCluster->clusterSize_);

	if (sizeBinAfter > sizeBinBefore &&
		sizeBefore				   < LARGE_CLUSTER_SIZE &&
		addedCluster->clusterSize_ < LARGE_CLUSTER_SIZE)
	{
		// redo distance peaks
		selectDistancePeaksFromMultipleSigneltons(mainCluster);
	}
	else
	{
		if (addedCluster->clusterSize_ > mainCluster->clusterSize_)
		{
			// use larger clsuters distance peaks
			mainCluster->distancePeaks_ = addedCluster->distancePeaks_;
		}
	}

	// check distance peaks
	assert( mainCluster->checkDistancePeaksOk() );
}


bool compIntensity( const DistancePeak& lhs, const DistancePeak& rhs)
{
	return (lhs.intensity > rhs.intensity);
}

/****************************************************************************
This function merges the singletons and selects the distance peaks, it uses
a slightly different method than the one used for selecting the ones in a single
spectrum since it weights intensities according to their frequency in all spectra.
This function deals with at most 255 singletons.

Note: this function does not treat all spectr equally. Beacause the cluster's
distance peaks get re-written each time, they do not reflect the peaks of
the original singleton (but rahter the current spectrum's consensus). This should
not be a problem beacuse there are many more distance peaks from unmodified singletons,
so on average this will have little effect.
*****************************************************************************/
void MsClusterDataStorage::selectDistancePeaksFromMultipleSigneltons(Cluster* cluster)
{
	const size_t peakAreaSize = 2*LARGE_CLUSTER_SIZE*MAX_NUM_PEAKS_FOR_DISTANCE;
	static vector<DistancePeak>  allPeaks;
	static vector<DistancePeak>  tmpPeakArea;
	static vector<unsigned char> peakCounts;
	static vector<float>		 ratios;(256);
	
	if (cluster->clusterSize_ <= 1 || cluster->clusterSize_>255)
		return;

	if (allPeaks.size() < peakAreaSize)
	{
		allPeaks.resize(peakAreaSize);
		tmpPeakArea.resize(peakAreaSize);
		peakCounts.resize(peakAreaSize);
		ratios.resize(256);
	}

	memcpy(&allPeaks[0], cluster->getDistancePeaks()->peaks, sizeof(DistancePeak)*cluster->getDistancePeaks()->numPeaks);
	int totalPeaks = cluster->getDistancePeaks()->numPeaks;

	const AllocatedVector<clusterIdx_t>* const singletonVector = cluster->getSingletonVector();
	const clusterIdx_t* clusterIdxs = singletonVector->getElementsPtr();
	for (size_t i=0; i<singletonVector->getSize(); i++)
	{
		Cluster* singleton = getCluster(clusterIdxs[i]);
		assert(singleton);
		const int numPeaks = singleton->getDistancePeaks()->numPeaks;
		memcpy(&allPeaks[totalPeaks], singleton->getDistancePeaks()->peaks, numPeaks*sizeof(DistancePeak));
		totalPeaks += numPeaks;
	}

	sort(allPeaks.begin(), allPeaks.begin()+totalPeaks);

	// set all counts to 1
	assert( cluster->clusterSize_ < 256 );
	memset(&peakCounts[0], 1, totalPeaks);

	// merge peaks
	tmpPeakArea[0]=allPeaks[0];
	const int maxIntProximity = convertMassToInt(Cluster::getIsoTolerance());
	int prev=0;
	for (int i=1; i<totalPeaks; i++)
	{
		if 	(allPeaks[i].massAsInt - tmpPeakArea[prev].massAsInt < maxIntProximity )
		{
			// join peaks with proportion to their intensities
			const intensity_t intensitySum =(tmpPeakArea[prev].intensity + allPeaks[i].intensity);
			const float		  ratio = tmpPeakArea[prev].intensity/intensitySum;
			const mass_t newMass = ratio * convertIntToMass(tmpPeakArea[prev].massAsInt) + 
				                   (1.0-ratio) * convertIntToMass(allPeaks[i].massAsInt);
			
			tmpPeakArea[prev].intensity = intensitySum;
			tmpPeakArea[prev].massAsInt = convertMassToInt(newMass);
			if (peakCounts[prev]+peakCounts[i]<255)
				peakCounts[prev]+=peakCounts[i];
			else
				peakCounts[prev]=255;
		}
		else
		{
			tmpPeakArea[++prev]= allPeaks[i];
			peakCounts[prev]   = peakCounts[i];
		}
	}
	totalPeaks = prev+1;

	if (totalPeaks < cluster->distancePeaks_.numPeaks)
		cluster->distancePeaks_.numPeaks = totalPeaks;

	// modify the intensity according to the peakWeightTable_
	// that is discount the weight of peaks that have only a few copies
	for (int i=0; i<totalPeaks; i++)
		tmpPeakArea[i].intensity = tmpPeakArea[i].intensity *
				peakWeightTable_.getWeight(static_cast<int>(peakCounts[i]),static_cast<int>(cluster->clusterSize_));

	// select a number of peaks according to their intensity
	sort(tmpPeakArea.begin(), tmpPeakArea.begin()+totalPeaks, compIntensity);

	memcpy(cluster->distancePeaks_.peaks, &tmpPeakArea[0], cluster->distancePeaks_.numPeaks * sizeof(DistancePeak));
	
	// sort according to mass
	sort(cluster->distancePeaks_.peaks, cluster->distancePeaks_.peaks + cluster->distancePeaks_.numPeaks);
	
	// set new adjusted intensity
	cluster->setAdjustedIntensities();
}


/**************************************************************************
Creates the consensus spectrum for the cluster.
Stores the peaks in the local static allocation.
(this destroys the original peak list pointer, so it should be backed up
if this information is still needed).
***************************************************************************/
void MsClusterDataStorage::makeConsensus(Cluster* cluster, 
										 SingleSpectrumHeader* clusterHeader)
{
	if (cluster->clusterSize_ <= 1 || ! cluster->singletonIdxVector_ )
		return;

	const size_t numSingletons = cluster->singletonIdxVector_->getSize();
	const clusterIdx_t* singletonIdxs = cluster->singletonIdxVector_->getElementsPtr();

	if (numSingletons == 0)
		return;

	size_t orgSize;
	if (params_->gotMergeArchives)
	{
		orgSize = cluster->clusterSize_;
		cluster->clusterSize_ = numSingletons + 1;
	}

	if (cluster->clusterSize_ <= MAX_CLUSTER_SIZE_FOR_UPDATE)
	{
		makeConsensusForSmallCluster(cluster, clusterHeader);
	}
	else
	{
		makeConsensusForLargeCluster(cluster, clusterHeader);
	}

	if (params_->gotMergeArchives)
		cluster->clusterSize_ = orgSize;
}


void MsClusterDataStorage::makeConsensusForSmallCluster(Cluster* cluster, 
														SingleSpectrumHeader* clusterHeader,
														bool	indSetConsensusParameters)
{
	const size_t numSingletons = cluster->singletonIdxVector_->getSize();
	const clusterIdx_t* singletonIdxs = cluster->singletonIdxVector_->getElementsPtr();
	static vector<Peak> tmpPeakArea;
	if (tmpPeakArea.size() < 100000)
		tmpPeakArea.resize(100000);

	float maxPrecursorIntensity = cluster->getHeader()->getPrecursorIntensity();
	size_t i;
	for (i=0; i<numSingletons; i++)
	{
		Cluster* singleton = getCluster(singletonIdxs[i]);
		assert( singleton );

		const float precursorIntensity = singleton->getHeader()->getPrecursorIntensity();
		if (precursorIntensity<=0.0)
			break;
		if (precursorIntensity>maxPrecursorIntensity)
			maxPrecursorIntensity=precursorIntensity;
	}
	if (i<numSingletons)
		maxPrecursorIntensity=0.0;

	size_t p=0;
	if (maxPrecursorIntensity>0.0)
	{
		float ratio = (cluster->getHeader()->getPrecursorIntensity() /maxPrecursorIntensity);
		if (ratio<0.2)
			ratio = 0.2;

		for (size_t j=0; j<cluster->numPeaks_; j++)
		{
			Peak& peak = tmpPeakArea[j];
			peak=cluster->peaks_[j];
			peak.intensity *= ratio;
		}							   
	}
	else
		memcpy(&tmpPeakArea[0], cluster->peaks_, cluster->numPeaks_*sizeof(Peak));

	p+= cluster->numPeaks_;

	// these vectors are filled and used for the function adjustMaxPossiblePeakCounts
	// which sets for each peak the maximal number of spectra in which it had the 
	// potential to appear
	vector<MassCount> minMassCounts(numSingletons+1), maxMassCounts(numSingletons+1);
	clusterIdx_t totalMassCounts=0;

	for (size_t i=0; i<numSingletons; i++)
	{
		Cluster* singleton = getCluster(singletonIdxs[i]);
		const size_t numPeaks = singleton->numPeaks_;
		assert (numPeaks>0 && numPeaks<10000);

		if (p+numPeaks >= tmpPeakArea.size())
			tmpPeakArea.resize(tmpPeakArea.size()*2);

		// copy peaks, but multiply the intensity by the ratio of the total
		// intensity of this spectrum comapred to the first main cluster's singleton
		// the idea is that a spectrum with higher total intensity has betterhigher s/n
		// TODO add sqs into ratio computation
		if (maxPrecursorIntensity>0.0)
		{
		    // NP3 fixed the intensity ratio to be between the current singleton and the main cluster consensus, instead of the max of the members
            float ratio = (singleton->getHeader()->getPrecursorIntensity() /maxPrecursorIntensity);
		    //float ratio = (singleton->getHeader()->getPrecursorIntensity()/cluster->getHeader()->getPrecursorIntensity());
			if (ratio<0.2)
				ratio = 0.2;
			for (size_t j=0; j<numPeaks; j++)
			{
				Peak& peak = tmpPeakArea[p+j];
				peak=singleton->peaks_[j];
				peak.intensity *= ratio;
			}
		}
		else
			memcpy(&tmpPeakArea[p], singleton->peaks_, sizeof(Peak)*numPeaks);

		p+=numPeaks;
	
		const clusterIdx_t size = singleton->getHeader()->getClusterSize(); // use this since it is the original size
		minMassCounts[i].mass   = singleton->getPeak(0).mass;
		minMassCounts[i].numSpectra = size;
		maxMassCounts[i].mass  = singleton->getPeak(singleton->getNumPeaks()-1).mass;
		maxMassCounts[i].numSpectra = size;
		totalMassCounts+=size;
	}

	// add info of original cluster
	assert( cluster->getHeader() && cluster->getPeaks() );
	const clusterIdx_t size = cluster->getHeader()->getClusterSize(); // use this since it is the original size
	minMassCounts[numSingletons].mass  = cluster->getPeak(0).mass;
	minMassCounts[numSingletons].numSpectra = size;
	maxMassCounts[numSingletons].mass  = cluster->getPeak(cluster->getNumPeaks()-1).mass;
	maxMassCounts[numSingletons].numSpectra = size;
	totalMassCounts+=size;

	if (! params_->gotMergeArchives)
	{
		if (totalMassCounts != cluster->clusterSize_)
		{
			for (size_t i=0; i<numSingletons; i++)
			{
				Cluster* singleton = getCluster(singletonIdxs[i]);
				cout << singleton->getHeader()->getClusterSize() << "\t";
				singleton->getHeader()->printStats();
			}
		}
		assert(totalMassCounts == cluster->clusterSize_);
	}
	
	// sort peaks and assign them to cluster
	sort(tmpPeakArea.begin(), tmpPeakArea.begin()+p);

	cluster->peaks_    = &tmpPeakArea[0];
	cluster->numPeaks_ = p;

// NP3 debug final cluster consensus
//	if (cluster->getHeader()->getMOverZ() > 293.08 && cluster->getHeader()->getMOverZ() < 293.1) {
//        cout << "@@@@@@@@@@@ AAAA precursor " << cluster->getHeader()->getMOverZ() << endl;
//        cluster->detailedPrint();
//        cluster->printPeaks();
//	}

	// process the peaks
	// counts and maxPosisble fields in peaks should be added (so paramter true is given)
	cluster->joinAdjacentPeaks(config_->getTolerance() ,true);

//    if (cluster->getHeader()->getMOverZ() > 293.08 && cluster->getHeader()->getMOverZ() < 293.1) {
//        cout << "@@@@@@@@@@@ AAAA precursor " << cluster->getHeader()->getMOverZ() << endl;
//        cluster->detailedPrint();
//        cluster->printPeaks();
//    }

	// convert the min/max counts to counts of the number of spectra that have peaks 
	// at a given mass
	sort(minMassCounts.begin(), minMassCounts.end());
	for (size_t i=1; i<minMassCounts.size(); i++)
		minMassCounts[i].numSpectra += minMassCounts[i-1].numSpectra;

	sort(maxMassCounts.begin(), maxMassCounts.end());
	for (size_t i=1; i<maxMassCounts.size(); i++)
		maxMassCounts[i].numSpectra += maxMassCounts[i-1].numSpectra;
	for (size_t i=0; i<maxMassCounts.size(); i++)
		maxMassCounts[i].numSpectra = totalMassCounts - maxMassCounts[i].numSpectra;

	cluster->adjustMaxPossiblePeakCounts(minMassCounts, maxMassCounts);

	// normalize intensities according to peakWeightTable_
	intensity_t totalPeakIntnesity = 0.0;
	for (size_t i=0; i<cluster->numPeaks_; i++)
	{
		Peak& peak = cluster->peaks_[i];
		const float ratio =  peakWeightTable_.getWeight(static_cast<int>(peak.count),
				static_cast<int>(peak.maxPossible));
		peak.intensity *= ratio;
		totalPeakIntnesity += peak.intensity;
	}
	

	// normalize peak intensity according to cluster size
	if (cluster->getClusterSize()>=1 && totalPeakIntnesity>0.0)
	{
		const intensity_t norm = (1000.0 * cluster->getClusterSize())/totalPeakIntnesity;
		for (size_t i=0; i<cluster->numPeaks_; i++)
			cluster->peaks_[i].intensity *= norm;
		
		// don't let this peak get filtered out, might cause problems with first peak mass in header
		if (cluster->peaks_[0].intensity < 0.001)
			cluster->peaks_[0].intensity = 0.001;
	}

	// remove the many weak ones (including a global removal for largish clusters)
	cluster->filterWeakPeaks(config_, 0.0, 0, (cluster->getClusterSize()>6) );

	if (clusterHeader)
	{
		if (clusterHeader->getTitle().length() == 0)
		{
			clusterHeader->setTitle(cluster->getHeader()->getTitle());
		}

		assert(cluster->numPeaks_>0);
		cluster->setHeaderFirstPeakMass(cluster->getPeak(0).mass);
		clusterHeader->setFirstPeakMass(cluster->getPeak(0).mass);
		clusterHeader->setOriginalNumPeaks(cluster->numPeaks_);
		if (indSetConsensusParameters)
			setConsensusParameters(cluster, clusterHeader);
	}

	assert( cluster->sanityCheck() );
}



/* With a large cluster we might  exceed the maximum of 255 peak counts for the fields for peak
   counts use unsigned char.*/
void MsClusterDataStorage::makeConsensusForLargeCluster(Cluster* cluster, SingleSpectrumHeader* clusterHeader)
{
	assert(cluster->getSingletonVector());

	// find maximum sized spectrum if more than half the size - use it as a consensus
	// check if one of the singletons is large enough

	const size_t sizeOfLarge = MAX_CLUSTER_SIZE_FOR_UPDATE/2 + 1;
	const size_t numSingletons = cluster->getSingletonVector()->getSize();
	const clusterIdx_t* singletonIdxs = cluster->getSingletonVector()->getElementsPtr();


	// There is no large cluster. Select a subset of singletons and create a consensus with them
	Cluster tmpCluster;
	SingleSpectrumHeader tmpHeader;
	tmpCluster.copyWithoutSingletonIdxVector(*cluster);
	tmpCluster.singletonIdxVector_ = clusterIdxAllocation_.allocateVector();
	tmpHeader = *(cluster->getHeader());
	tmpCluster.setHeader(&tmpHeader);

	size_t totalSize = tmpCluster.getHeader()->getClusterSize();
	assert(totalSize>0 && totalSize < sizeOfLarge);

	// start adding clusters, until we add one that leads to a size that is too large 
	for (size_t i=0; i<numSingletons; i++)
	{
		const Cluster* singleton = getCluster(singletonIdxs[i]);
		const size_t singletonSize = singleton->getHeader()->getClusterSize();
		// TODO change to continue
		if (totalSize + singletonSize > MAX_CLUSTER_SIZE_FOR_UPDATE)
			break;

		assert(! singleton->getIndInPlay());
		tmpCluster.singletonIdxVector_->addElement(singleton->getClusterIdx());
		totalSize += singletonSize;
	}

	assert( tmpCluster.singletonIdxVector_->getSize() < cluster->singletonIdxVector_->getSize());

	// the setting of the consensus parameters should be done with all the spectra
	tmpCluster.peaks_    = cluster->peaks_;
	tmpCluster.numPeaks_ = cluster->numPeaks_;
	tmpCluster.clusterSize_ = totalSize;
	makeConsensusForSmallCluster(&tmpCluster, const_cast<SingleSpectrumHeader*>(tmpCluster.getHeader()), false);

	cluster->peaks_    = tmpCluster.peaks_;
	cluster->numPeaks_ = tmpCluster.numPeaks_;
	cluster->setHeaderFirstPeakMass(cluster->getPeak(0).mass);

	setConsensusParameters(cluster, clusterHeader);

	assert(cluster->getClusterSize() >= MAX_CLUSTER_SIZE_FOR_UPDATE);
	assert( cluster->sanityCheck() );
}


// This function sets various consensus parameters into the supplied header
// Should be run after the consensus spectrum is created
void MsClusterDataStorage::setConsensusParameters(Cluster* cluster, SingleSpectrumHeader* clusterHeader)
{
	// set other parameters into the cluser header
	if (! clusterHeader)
		clusterHeader = const_cast<SingleSpectrumHeader*>(cluster->getHeader());

	clusterHeader->setClusterSize(cluster->getClusterSize());
	
	const AllocatedVector<clusterIdx_t>* singletonVector = cluster->getSingletonVector();
	const clusterIdx_t* singletonIdxs  = singletonVector->getElementsPtr();
	const size_t		numSingletons  = singletonVector->getSize();

	assert(numSingletons>0);

	// NP3 GOT set RT width
	float  totalRtMin=0.0;
	float  totalRtMax=0.0;

	float  totalRt=0.0;
	double rtWeight=0.0;
	double totalMz=0.0;
	double totalPrecursorIntensity=0.0;
	map<string,int> pepStrCounts;
	map<short, int> chargeCounts;

	// collect singleton pointers (including the cluster's
	vector<const Cluster*> singletons(numSingletons+1,0);
	for (size_t i=0; i<numSingletons; i++)
		singletons[i]=getCluster(singletonIdxs[i]);
	singletons[numSingletons]=cluster;

	// NP3 got set rtmin and max of blanks as baseline
	// count if a blank mz was present
	int blank_mz = 0;
	// collect data from all singletons
	for (size_t i=0; i<=numSingletons; i++)
	{
		const Cluster* singleton = singletons[i];
		const intensity_t precursorIntensity = singleton->getHeader()->getPrecursorIntensity();
		totalMz+=singleton->getHeader()->getMOverZ() * precursorIntensity;
		const float rt = singleton->getHeader()->getRetentionTime();

		if (rt>=0.0)
		{
			totalRt+=rt*precursorIntensity;
			rtWeight+=precursorIntensity;
		}

		// NP3 GOT rt width
		const float rtMin = singleton->getHeader()->getRetentionTimeMin();
		const float rtMax = singleton->getHeader()->getRetentionTimeMax();
		if (rtMin>0.0)
		{
			totalRtMin+=rtMin*precursorIntensity;
		}
		if (rtMax<1000000.0)
		{
			totalRtMax+=rtMax*precursorIntensity;
		}
		if (rtMin == 0.0 && rtMax == 1000000.0) {
			// NP3 set rtmin and max of blanks as baseline
			blank_mz++;
		}

		totalPrecursorIntensity += precursorIntensity;

		if (singleton->getHeader()->getPeptideStr().length()>0)
			pepStrCounts[singleton->getHeader()->getPeptideStr()]++;

		chargeCounts[singleton->getHeader()->getCharge()]++;
	}

	if (totalPrecursorIntensity<=0.0)
	{
		cout << "Problem with cluster " << cluster->getClusterIdx() << " has " << singletons.size() << " singltons." << endl;
		for (size_t i=0; i<=numSingletons; i++)
		{
			const Cluster* singleton = singletons[i];
			const intensity_t precursorIntensity = singleton->getHeader()->getPrecursorIntensity();
			cout << singleton->getHeader()->getDatasetIndex() << "\t" << singleton->getHeader()->getSpectraFileIndexInList()
				<< "\t" << singleton->getHeader()->getScanNumber() << "\tm/z:" << singleton->getHeader()->getMOverZ()
				<< " \tinten: " << precursorIntensity << endl;
		}
	}

	assert(totalPrecursorIntensity>0.0);

	if (rtWeight>0) {
		clusterHeader->setRetentionTime(totalRt / rtWeight);
		// NP3 GOT add rt width
		if (blank_mz > 0) {
			// has a blank mz among its members, set rt as baseline
			clusterHeader->setRetentionTimeMin(0.0);
			clusterHeader->setRetentionTimeMax(1000000.0);
		} else {
			// compute the average rt
            clusterHeader->setRetentionTimeMin(totalRtMin / rtWeight);
            clusterHeader->setRetentionTimeMax(totalRtMax / rtWeight);
        }
	}

	clusterHeader->setPrecursorIntensity(totalPrecursorIntensity);
	clusterHeader->setMOverZ(totalMz / totalPrecursorIntensity);
	cluster->clusterMOverZ_ = clusterHeader->getMOverZ();

	// GoT need change?? mx diff
	// check for errors computing the consensus m/z
	if (fabs(cluster->clusterMOverZ_-singletons[0]->getHeader()->getMOverZ())>10.0)
	{
		double totalMz=0.0;
		double totalPrecursorIntensity=0;
		cout << "PROBLEM WITH CLUSTER " << cluster->getClusterIdx() << " has " << singletons.size() << " singltons." << endl;
		for (size_t i=0; i<=numSingletons; i++)
		{
			const Cluster* singleton = singletons[i];
			const intensity_t precursorIntensity = singleton->getHeader()->getPrecursorIntensity();
			totalMz+=singleton->getHeader()->getMOverZ() * precursorIntensity;
			totalPrecursorIntensity += precursorIntensity;

			cout << i << "\t" << singleton->getHeader()->getMOverZ() << "\t";
			for (size_t j=0; j<NUM_PEAKS_FOR_HEURISTIC; j++)
				cout << "\t" << singleton->getTopPeakIdxs()[j];
			cout << "\t" << precursorIntensity << endl;

			if (singleton->getHeader()->getPeptideStr().length()>0)
				pepStrCounts[singleton->getHeader()->getPeptideStr()]++;

			chargeCounts[singleton->getHeader()->getCharge()]++;
		}
		cout << endl << "Total mz: " << totalMz << "\t Total intensity: " << totalPrecursorIntensity << endl;
		cout << "clusterMz: " << totalMz / totalPrecursorIntensity << endl;
	}
	// NP3 GOT change accepted mz diff in window <= 10 to <= 1000
	assert(fabs(cluster->clusterMOverZ_-singletons[0]->getHeader()->getMOverZ())<=1000.0);

	// select the peptide str
	if (pepStrCounts.size()>0)
	{
		if (pepStrCounts.size()==1)
		{
			clusterHeader->setPeptideStr(pepStrCounts.begin()->first);
		}
		else // go with max
		{
			map<string,int>::const_iterator it;
			const string* maxStr =0;
			int maxCount	   =0;
			for (it = pepStrCounts.begin(); it != pepStrCounts.end(); it++)
				if (it->second > maxCount)
				{
					maxCount = it->second;
					maxStr   = &it->first;
				}
			if (maxStr)
				clusterHeader->setPeptideStr(*maxStr);
		}

		// set charge according to string
		if (clusterHeader->getPeptideStr().length()>0)
		{
			Peptide pep;
			pep.parseFromString(config_, clusterHeader->getPeptideStr());

			int charge = static_cast<int>((pep.get_mass_with_19() / clusterHeader->getMOverZ()) + 0.1);
		//	assert( charge>0 && fabs(clusterHeader->getMOverZ() * charge - (charge+1) - pep.get_mass_with_19())<10.0 );
			clusterHeader->setCharge(charge);
		}
	}

	// select the most abundant charge
	short maxCharge=0;
	int   maxCount=0;
	map<short, int>::const_iterator it;
	for (it = chargeCounts.begin(); it != chargeCounts.end(); it++)
		if (it->second > maxCount)
		{
			maxCount = it->second;
			maxCharge = it->first;
		}

	clusterHeader->setCharge(maxCharge);
	cluster->clusterCharge_ = maxCharge;
}


bool MsClusterDataStorage::testSingletonsWithGoodMzs(const Cluster* cluster, mass_t windowSize) const
{
	mass_t minMz = cluster->getHeader()->getMOverZ();
	mass_t maxMz = minMz;

	if (cluster->getSingletonVector())
	{
		const size_t numSingletons = cluster->getSingletonVector()->getSize();
		const clusterIdx_t* singletonIdxs = cluster->getSingletonVector()->getElementsPtr();
		for (size_t i=0; i<numSingletons; i++)
		{
			const Cluster* singleton = getCluster(singletonIdxs[i]);
			
			assert(singleton && ! singleton->getIndInPlay());
			assert(singleton->getAssignedClusterIdx() == cluster->getClusterIdx());
			assert(singleton->getClusterMOverZ() == singleton->getHeader()->getMOverZ());

			const mass_t mz = singleton->getHeader()->getMOverZ();
			assert(mz>0.0);
			if (mz<minMz)
				minMz=mz;
			if (mz>maxMz)
				maxMz=mz;
		}
	}

	// assume that there is no way to have a correct joining of spectra with such a large m/z difference
	if (maxMz-15.0>minMz)
	{
		cout << endl << "DEBUG THIS..." << endl;
		cout << cluster->getClusterIdx() << "\t" << cluster->getClusterCharge() << "\t" << cluster->getClusterMOverZ() << endl;
		const size_t numSingletons = cluster->getSingletonVector()->getSize();
		const clusterIdx_t* singletonIdxs = cluster->getSingletonVector()->getElementsPtr();
		for (size_t i=0; i<numSingletons; i++)
		{
			const Cluster* singleton = getCluster(singletonIdxs[i]);
			
			assert(singleton && ! singleton->getIndInPlay());
			assert(singleton->getAssignedClusterIdx() == cluster->getClusterIdx());
			assert(singleton->getClusterMOverZ() == singleton->getHeader()->getMOverZ());

			const mass_t mz = singleton->getHeader()->getMOverZ();
			cout << singletonIdxs[i] << "\t" << mz << endl;
		}
	}
	return (maxMz-minMz<15.0);
}








