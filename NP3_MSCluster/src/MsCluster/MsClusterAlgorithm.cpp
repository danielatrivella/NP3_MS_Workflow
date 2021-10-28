#include "MsClusterAlgorithm.h"
#include "MsParameterStruct.h"
#include "MsClusterAuxfuns.h"
#include "MsArchive.h"
#include "../PepNovo/SpectraList.h"
#include "../PepNovo/PeakList.h"
#include "../PepNovo/AnnotationFile.h"
#include "../PepNovo/PepNovo_auxfun.h"


void MsClusterAlgorithm::performClustering(
						   MsParameterStruct* params,
						   const AllScoreModels* model,
						   bool  indGreedyClusterJoin)
{
	model_		= model;
	outDir_		= params->outDir;
	name_		= params->outputName;
	datasetIdx_ = params->datasetIdx;
	batch_		= params->batchIdx;
	window_		= params->mzWindow;
	verboseLevel_ = params->verboseLevel;
	indGreedyClusterJoin_	= indGreedyClusterJoin;

	const double clusteringStartTime = time(NULL);

	simModel_.readSimilarityModel(model->get_config());

	setSimilarityThresholds(params->minSimilarity, params->numRounds);

	data_.allocateMemory(params->memoryGb);

	data_.initialize(params, model_->get_config(), &simModel_);

	size_t iter=0;
	clusterIdx_t totalNumberSpectraRead = 0;

	// main loop
	while (data_.runningClusterIdx_ < data_.getTotalSpectraToCluster())
	{
		assert( data_.nextPeakPos_ <= data_.peakAlloc_.size() );
		assert( data_.nextClusterPos_ <= data_.clusterAlloc_.size() );

		clusterIdx_t numSpectraFreed = 0;
		const clusterIdx_t nextIdxToBeClustered    = data_.runningClusterIdx_;
		const clusterIdx_t numSpectraRead		   = data_.addNewSpectra(window_, &numSpectraFreed);
		const size_t       posNextIdxToBeClustered = data_.getCluserPosition(nextIdxToBeClustered);

		assert( data_.nextPeakPos_ <= data_.peakAlloc_.size() );
		assert( data_.nextClusterPos_ <= data_.clusterAlloc_.size() );

		if (verboseLevel_>0)
		{
			cout << ++iter << "\t#Freed " << numSpectraFreed << "\t#Read " << numSpectraRead
				<< "\t#In memory " << data_.nextClusterPos_ << endl; 
		}

		totalNumberSpectraRead += numSpectraRead;
		if (numSpectraRead == 0)
			break;

		if (verboseLevel_>0)
		{
			cout << "Clustering: [ " << nextIdxToBeClustered << " at " << posNextIdxToBeClustered
				 << " m/z = " << setprecision(2) << fixed << data_.clusterAlloc_[posNextIdxToBeClustered].getClusterMOverZ() << " ]";
			cout << "\t   to\t" << " [ "<< data_.clusterAlloc_[data_.nextClusterPos_-1].getAssignedClusterIdx()
				 << " at " << data_.nextClusterPos_-1 << " m/z = " << data_.clusterAlloc_[data_.nextClusterPos_-1].getClusterMOverZ() << " ]" << endl;
		}

		// fill in-play vector for all clusters that are involved
		// basically checks for all clusters if (clusterIdx == assignedClusterIdx)
		// storing this info as a bit vector means all indexes can be kept in cache simultaneously
		fillInPlayIndicators();

		// clustering is done in rounds (with decreasing similarity thresholds)
		// First round involves many similarity computations that are performed in batches of spectra
		// (sliding window across all the spectra in clusterAlloc_)

		const size_t nextClusterPos = data_.getNextClusterPos();
		for (size_t blockStartPos=posNextIdxToBeClustered; blockStartPos<nextClusterPos;  )
		{
			const mass_t minClusterMz = data_.clusterAlloc_[blockStartPos].getClusterMOverZ() - window_;
			const size_t firstPos = data_.findPosOfClusterWithMzAbove(minClusterMz);
		
			// the spectra in range [pos,lastPos) are going to be compared to all spectra in range
			// [firstPos,lastPos) that appear before them (subject to being in lists, not assigned
			// to different clusters, etc.
			// the first step is to make a large sorted list of all pairs of idxs that need to be compared)
			const size_t plannedLastPos = (blockStartPos + SIMILARITY_BATCH_SIZE > nextClusterPos ? 
										   nextClusterPos : blockStartPos + SIMILARITY_BATCH_SIZE);

			// the actualLastPos might be smaller than the plannedLastPos if the number of similarities
			// that need to be computed is too large (so we don't end up comparing a whole block of
			// SIMILARITY_BATCH_SIZE spectra.
			const size_t actualLastPos = data_.makeListOfSimilarities(firstPos, plannedLastPos, blockStartPos, plannedLastPos, window_);

			data_.computeMinSimilarityForJoiningFromLists( params->maxMixtureProb, firstPos, actualLastPos);
			// compare all new spectra in allocation against all clusters in their window that have
			// a smaller m/z mass
			performFirstRoundSpectraComparisons();

			// if the greedy join was done, we don't need to join since the spectra have already
			// been joined during the similarity computation stage
			if (! indGreedyClusterJoin_)
				joinAllClustersThatPassThreshold(0, blockStartPos, actualLastPos);

			blockStartPos = actualLastPos;
		}

		// perform additional rounds of clustering
		for (int round=1; round<numRounds_; round++)
			performSpectraComparisons(round, posNextIdxToBeClustered, data_.nextClusterPos_);

	}

	// flush remaining clusters
	data_.writeClustersToOutput(data_.nextClusterPos_);

	// report on computations and clusters written
	if (verboseLevel_>=0)
	{
		cout << endl << "Done clustering files." << endl;
		cout << "Time (without dat creation): " << time(NULL)-clusteringStartTime << endl;
		cout << endl << "Final statistics:" << endl;
		cout << "-----------------" << endl;

		cout << endl << "Total similarity computations: " << scientific << setprecision(3)
			<< totalSpectraComps_ << " average per spectrum " << fixed << setprecision(2) 
			<< static_cast<double>(totalSpectraComps_)/static_cast<double>(data_.getTotalSpectraToCluster()) << endl;

		clusterIdx_t totalJoins = 0;
		for (size_t i=0; i<numRounds_; i++)
			totalJoins += numJoinsPerRound_[i];

		cout << endl << "Breakdown of spectra joins:" << endl;
		cout << "---------------------------" << endl;
		cout << "Round\tThresh\t#joins\t% of total" << endl;
		for (size_t i=0; i<numRounds_; i++)
			cout << i+1 << "\t" << similarityThresholds_[i]  << "\t" << numJoinsPerRound_[i]
				<< "\t(" << setprecision(3) << fixed 
				<< static_cast<double>(numJoinsPerRound_[i])/static_cast<double>(totalJoins)
				<< ")" << endl;

		cout << endl;
		data_.clusterOutputter_.reportClusterStats();
		cout.flush();
	}

	data_.clusterOutputter_.closeAll(); // closing collects the file names

	if (! params->gotMergeArchives && params->gotCreateArchive)
	{
		data_.archive1_.writeArchive(data_.clusterOutputter_, params, data_.config_);
	}
}


void MsClusterAlgorithm::mergeTwoArchives(MsParameterStruct* params, const AllScoreModels* model)
{
	MsArchive& archive1 = data_.archive1_;
	MsArchive& archive2 = data_.archive2_;

	if (! archive1.readArchive(params->pathToArchive1))
		error("Could not read archive: ",params->pathToArchive1.c_str());

	if (! archive2.readArchive(params->pathToArchive2))
		error("Could not read archive: ",params->pathToArchive2.c_str());

	if (archive1.getArchiveName() == archive2.getArchiveName())
		error("When merging archives names should be different!, both archives names are: ", archive1.getArchiveName().c_str());

	archive1.generateArchiveDatPaths(params->datPaths1);
	archive2.generateArchiveDatPaths(params->datPaths2);

	archive1.updateGenerationIfNeeded(params->outputName);

	performClustering(params, model);

	archive1.addArchive(params, archive2);

	archive1.writeArchive(data_.clusterOutputter_, params, model->get_config());
}


// fill in-play vector for all clusters that are involved
// basically checks for all clusters if (clusterIdx == assignedClusterIdx)
// storing this info as a bit vector means all indexes can be kept in cache simultaneously
void MsClusterAlgorithm::fillInPlayIndicators()
{
	const size_t lastPos = data_.nextClusterPos_;
	for (size_t i=0; i<lastPos; i++)
	{
		const Cluster& cluster = data_.clusterAlloc_[i];
		const bool inPlay = (cluster.getAssignedClusterIdx() < MAX_CLUSTER_IDX &&
							 cluster.getIndInPlay() &&
							 cluster.getClusterIdx() == cluster.getAssignedClusterIdx());
	}
}



struct ValueIdx {
	ValueIdx() : value(0.0), idx(0) {}
	ValueIdx(float v, clusterIdx_t i) : value(v), idx(i) {}
	const bool operator< (const ValueIdx& rhs) const
	{
		return (value<rhs.value);
	}
	float        value;
	clusterIdx_t idx;
};

/**************************************************************************
Examines all clusters that that are candidates for joining with other clusters
(that have a smaller idx). This is done by examining each clusters list of
highest similarity tragets.
***************************************************************************/
void MsClusterAlgorithm::joinAllClustersThatPassThreshold(int round, size_t firstPos, size_t lastPos)
{
    const float similarityThreshold = similarityThresholds_[round];
    const bool  lastRound			= (round >= similarityThresholds_.size()-1);
    const MsParameterStruct* params = data_.params_;

    for (size_t i=firstPos; i<lastPos; i++)
    {
        Cluster& cluster = data_.clusterAlloc_[i];
        const float* const bestSimilarities = cluster.getBestSimilarityValues();
        const clusterIdx_t* const bestIdxs  = cluster.getBestSimilarityClusterIdxs();
        size_t maxSimilarityPos = 0;
        for (size_t j=1; j<NUM_TOP_SIMILARITIES_TO_SAVE; j++)
            if (bestSimilarities[j]>bestSimilarities[maxSimilarityPos])
                maxSimilarityPos=j;

        if (bestSimilarities[maxSimilarityPos]<similarityThreshold)
            continue;

        Cluster* topCluster = data_.getCluster(bestIdxs[maxSimilarityPos]);
        if (topCluster &&
            (topCluster->getClusterSize() < LARGE_CLUSTER_SIZE ||
             bestSimilarities[maxSimilarityPos] >= similarityForLargeClusters_) )
        {
            data_.joinClusters(topCluster, &cluster, 1);
            continue;
        }

        // the cluster we want to add might be similar enough to some other suboptimal (but not large)
        // clusters. Try and see if there is such a cluster here.

        vector<ValueIdx> aboveThreshold;
        aboveThreshold.clear();
        for (size_t j=0; j<NUM_TOP_SIMILARITIES_TO_SAVE; j++)
            if (bestSimilarities[j]>similarityThreshold)
                aboveThreshold.push_back(ValueIdx(bestSimilarities[j],bestIdxs[j]));

        if (aboveThreshold.size()<=1)
            continue;

        sort(aboveThreshold.begin(),aboveThreshold.end());
        for (size_t j=1; j<aboveThreshold.size(); j++)
        {
            Cluster* otherCluster = data_.getCluster(aboveThreshold[j].idx);
            if (! otherCluster || otherCluster->getClusterSize()>= LARGE_CLUSTER_SIZE)
                continue;
            data_.joinClusters(otherCluster, &cluster, 2);
            break;
        }
    }
}

/**************************************************************************
Performs the comparisons of a block of spectra against all other in its
window.
***************************************************************************/
void MsClusterAlgorithm::performFirstRoundSpectraComparisons()
{	
	const float similarityToJoin = (indGreedyClusterJoin_ ? similarityThresholds_[0] : 1.0);
	const size_t numSimilaities = data_.similarityPairPositions_.size();
	const MsParameterStruct* params = data_.params_;

	//cout << "@@@@@@@@@@@@RT tol = " << Cluster::getRtTolerance() << endl;

	for (size_t i=0; i<numSimilaities; i++)
	{
		const ClusterPositionPair& pair = data_.similarityPairPositions_[i];
		Cluster& lowerPosCluster		= data_.clusterAlloc_[pair.lowerPos];
		Cluster& higherPosCluster		= data_.clusterAlloc_[pair.higherPos];

		// NP3 GOT accepted mz diff in window from < 7.0 to 1000.0
		assert(higherPosCluster.getClusterMOverZ() - lowerPosCluster.getClusterMOverZ() < 1000.0);

		if (! lowerPosCluster.getIndInPlay() || ! higherPosCluster.getIndInPlay())
			   continue;

        // don't compare clusters that are too far away from each other
        //if (abs(lowerPosCluster.getClusterRT() - higherPosCluster.getClusterRT()) > Cluster::getRtTolerance())
        //    continue;
        // NP3 GOT sim computation
		if (!((lowerPosCluster.getClusterRT() >= higherPosCluster.getClusterRTMin() - Cluster::getRtTolerance() &&
			lowerPosCluster.getClusterRT() <= higherPosCluster.getClusterRTMax() + Cluster::getRtTolerance()) ||
			(higherPosCluster.getClusterRT() >= lowerPosCluster.getClusterRTMin() - Cluster::getRtTolerance() &&
			higherPosCluster.getClusterRT() <= lowerPosCluster.getClusterRTMax() + Cluster::getRtTolerance())))
			continue;

        const float similarity = computeSimilarity(lowerPosCluster.getDistancePeaks(),
        										   higherPosCluster.getDistancePeaks(),
        										   Cluster::getPeakIndexToleranceAsInt());

		totalSpectraComps_++;

		if (similarity < minSimilarityForAnything_)
			continue;

		const float minSimilarityToJoinClusters = (higherPosCluster.getMinSimilarityToJoin() > lowerPosCluster.getMinSimilarityToJoin() ?
										       higherPosCluster.getMinSimilarityToJoin() : lowerPosCluster.getMinSimilarityToJoin());

		assert(minSimilarityToJoinClusters<=1.0);

		if (similarity >= similarityToJoin && similarity >= minSimilarityToJoinClusters)
		{
			data_.joinClusters(&lowerPosCluster, &higherPosCluster, 3);
			numJoinsPerRound_[0]++;
		}
		else
		{
			//cout << "Not joined mz " << lowerPosCluster.getClusterMOverZ() << " RT " << lowerPosCluster.getClusterRT() << " and mz " << higherPosCluster.getClusterMOverZ() << " RT " << higherPosCluster.getClusterRT() << " sim " << similarity << " simToJOin " << similarityToJoin << " minSimToJoin " << minSimilarityToJoinClusters << endl;
			lowerPosCluster.updateSimilarityValues(similarity, higherPosCluster.getClusterIdx());
			higherPosCluster.updateSimilarityValues(similarity, lowerPosCluster.getClusterIdx());
		}
	}
}



/**************************************************************************
Performs the comparisons of a block of spectra against all other in its
window.
***************************************************************************/
void MsClusterAlgorithm::performSpectraComparisons(int round, size_t firstPos, size_t lastPos)
{
	assert(round>0);
	const float similarityToJoin = similarityThresholds_[round];
	const float maxMzDiff		 = 2.0 * window_;
	const MsParameterStruct* params = data_.params_;

	for (size_t pos=firstPos; pos<lastPos; pos++)
	{
		Cluster* firstCluster = &data_.clusterAlloc_[pos];
		if (! firstCluster->getIndInPlay())
			continue;

		const clusterIdx_t* bestSimilarityClusterIdxs = firstCluster->getBestSimilarityClusterIdxs();

		for (size_t i=0; i<NUM_TOP_SIMILARITIES_TO_SAVE; i++)
		{
			const clusterIdx_t otherIdx = bestSimilarityClusterIdxs[i];

			if (otherIdx == MAX_CLUSTER_IDX)
				continue;
			
			Cluster* otherCluster = data_.getCluster(otherIdx);
			assert(otherCluster);
			
			const clusterIdx_t otherAssignedClusterIdx = otherCluster->getAssignedClusterIdx();


			// if assigned to this cluster no need to compare
			if (otherAssignedClusterIdx == firstCluster->getClusterIdx())
			{
				firstCluster->updateSimilarityValueInPosition(i, 0.0, MAX_CLUSTER_IDX);
				continue;
			}

			// if the other cluster got assigned, we will compare with the cluster it got assigned to
			if (otherAssignedClusterIdx != otherIdx)
			{
				size_t j;
				for (j=0; j<NUM_TOP_SIMILARITIES_TO_SAVE; j++)
					if (bestSimilarityClusterIdxs[j] == otherAssignedClusterIdx)
						break;

				if (j<NUM_TOP_SIMILARITIES_TO_SAVE) // we already compared to this one
					continue;

				assert(! otherCluster->getIndInPlay());
				otherCluster = data_.getCluster(otherAssignedClusterIdx);
				assert(otherCluster);
			}

            //if (abs(firstCluster->getClusterRT() - otherCluster->getClusterRT()) > Cluster::getRtTolerance())
            //    continue;
            // NP3 GOT rt width comparision
			if (!((firstCluster->getClusterRT() >= otherCluster->getClusterRTMin() - Cluster::getRtTolerance() &&
				   firstCluster->getClusterRT() <= otherCluster->getClusterRTMax() + Cluster::getRtTolerance()) ||
				  (otherCluster->getClusterRT() >= firstCluster->getClusterRTMin() - Cluster::getRtTolerance() &&
				   otherCluster->getClusterRT() <= firstCluster->getClusterRTMax() + Cluster::getRtTolerance())))
				continue;

			// don't compare clusters that are too far away from each other
			if (firstCluster->getClusterMOverZ() - maxMzDiff > otherCluster->getClusterMOverZ() ||
				otherCluster->getClusterMOverZ() - maxMzDiff > firstCluster->getClusterMOverZ())
				continue;

			assert(firstCluster->getIndInPlay() && otherCluster->getIndInPlay());

			// NP3 GOT sim computation
			const float similarity = computeSimilarity(firstCluster->getDistancePeaks(),
													   otherCluster->getDistancePeaks(),
													   Cluster::getPeakIndexToleranceAsInt());

			totalSpectraComps_++;
			
			firstCluster->updateSimilarityValueInPosition(i, similarity, otherCluster->getClusterIdx());
			otherCluster->updateSimilarityValuesCheckIfClusterThere(similarity, firstCluster->getClusterIdx());
		}

        // look for max similarity value
		clusterIdx_t clusterIdxToAssignTo=MAX_CLUSTER_IDX;
		const float maxSimilarityToOther = data_.getMaxSimilarityForCluster(firstCluster, clusterIdxToAssignTo);
		if (maxSimilarityToOther >= similarityToJoin)
		{
			Cluster* bestOtherCluster = data_.getCluster(clusterIdxToAssignTo);

			//if (abs(firstCluster->getClusterRT() - bestOtherCluster->getClusterRT()) > Cluster::getRtTolerance())
			//	continue;
			// NP3 GOT rt width comparision
			if (!((firstCluster->getClusterRT() >= bestOtherCluster->getClusterRTMin() - Cluster::getRtTolerance() &&
				   firstCluster->getClusterRT() <= bestOtherCluster->getClusterRTMax() + Cluster::getRtTolerance()) ||
				  (bestOtherCluster->getClusterRT() >= firstCluster->getClusterRTMin() - Cluster::getRtTolerance() &&
				   bestOtherCluster->getClusterRT() <= firstCluster->getClusterRTMax() + Cluster::getRtTolerance())))
				continue;

			const float minSimilarityToJoinClusters = (firstCluster->getMinSimilarityToJoin() > bestOtherCluster->getMinSimilarityToJoin() ?
													   firstCluster->getMinSimilarityToJoin() : bestOtherCluster->getMinSimilarityToJoin());

			assert(minSimilarityToJoinClusters<=1.0);
			if (maxSimilarityToOther > minSimilarityToJoinClusters)
			{
				if (bestOtherCluster->getClusterIdx() < firstCluster->getClusterIdx())
				{
					data_.joinClusters(bestOtherCluster, firstCluster, 4);
				}
				else
					data_.joinClusters(firstCluster, bestOtherCluster, 5);
			} else {
			   // cout << "NOT JOINING - minSimiToJoin Clusters " << minSimilarityToJoinClusters << " max sim To Other " << maxSimilarityToOther << " mz " << firstCluster->getClusterMOverZ() << endl;
			}

			numJoinsPerRound_[round]++;
		} else {
            //cout << "NOT JOINING - similarityToJoin Clusters " << similarityToJoin << " > max sim To Other " << maxSimilarityToOther << " mz " << firstCluster->getClusterMOverZ() << endl;
		}
	}
    //joinAllClustersThatPassThreshold(round, firstPos, lastPos);
}

void MsClusterAlgorithm::setSimilarityThresholds(float similarity, int numRounds)
{
	similarity_ = similarity;
	numRounds_  = numRounds;

	minSimilarityForAnything_ = 0.15;

	totalSpectraComps_ = 0;
	numJoinsPerRound_.clear();
	numJoinsPerRound_.resize(numRounds, 0);

	if (similarity >= MAX_SIMILARITY*0.9 || numRounds == 1)
	{
		similarityThresholds_.resize(1,similarity);
		similarityForLargeClusters_ = similarity;
		numRounds_ = 1;
		return;
	}

	similarityThresholds_.resize(numRounds_);
	similarityThresholds_[0] = MAX_SIMILARITY;
	similarityThresholds_[numRounds_-1] = similarity;
	if (numRounds_ == 2)
	{
		similarityForLargeClusters_ = (MAX_SIMILARITY + similarity)*0.5;
		return;
	}
	
	float delta = (MAX_SIMILARITY-similarity) / (static_cast<float>(numRounds_) - 1.0);
	for (size_t i=1; i<numRounds-1; i++)
		similarityThresholds_[i] = similarityThresholds_[i-1]-delta;

	similarityForLargeClusters_ = similarityThresholds_[numRounds/2];
}








