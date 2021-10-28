#include "MsClusterBenchmark.h"
#include "MsClusterAlgorithm.h"
#include "MsParameterStruct.h"
#include "MsClusterAuxfuns.h"
#include "../PepNovo/PmcsqsAdditionalFunctions.h"

void performBenchmark(AllScoreModels* model, MsParameterStruct* params)
{
	if (params->gotBenchmarkSqs)
	{
		benchmarkSqs(model, params->list.c_str());
		return;
	}

	if (params->gotBenchmarkSimilarity)
	{
		MsClusterAlgorithm alg;
		alg.benchmarkSimilarity(model, params->nonClusterList, params->clusterList, params->mzWindow, params->simType);
		return;
	}

	if (params->gotBenchmarkConsensus)
	{
		MsClusterAlgorithm alg;
		alg.benchmarkConsensus(model, params->clusterList, params->minSimilarity, params->numRounds);
		return;
	}

	if (params->gotBenchmarkPairs)
	{
		MsClusterAlgorithm alg;
		alg.benchmarkPairs(model, params);
		return;
	}

	if (params->gotBenchmarkPm)
	{
		PMCSQS_Scorer* nonConstPMCSQS = const_cast<PMCSQS_Scorer*>(model->get_pmcsqs_ptr());
		nonConstPMCSQS->benchmarkPm(model->get_config(), params->list.c_str());
		return;
	}

	if (params->gotBenchmarkSimilarityHistogram)
	{
		MsClusterAlgorithm alg;
		alg.benchmarkSimilarityHistogram(model, params);
		return;
	}
}



/**************************************************************************
creates a recall-precision plot for the similarity computations.
Two distributions are compared:
- The nonclusterList of files is used to generate a distribution of the maximum
similarity a spectrum can have to all others (of a different peptide).
- The clusterList is used to to generate a distribution of the inner cluster
similarity (comparing spectra of the same peptide/charge).
***************************************************************************/
void MsClusterAlgorithm::benchmarkSimilarity(const AllScoreModels* model,
											 const string& nonClusterList,
											 const string& clusterList,
											 mass_t window,
											 int simType)
{
	model_	 = model;
	window_	 = window;

	Cluster::setTolerances(model_->get_config()->getTolerance());


	vector<double> sameClusterSimilarities;
	vector< vector<double> > differentClusterSimilarities;

	cout << "Using sim-type " << simType << endl;

	cout << endl << "1. Computing same cluster similarities:" << endl;
	cout <<         "***************************************" << endl << endl;
	computeSameClusterSimilarities(clusterList, sameClusterSimilarities, simType);
	cout << "Histogram of similarities:" << endl << fixed << setprecision(3);
	double cdf=0.0;
	for (size_t i=0; i<sameClusterSimilarities.size(); i++)
	{
		cdf+=sameClusterSimilarities[i];
		cout << i*0.01 << "\t" << sameClusterSimilarities[i] << "\t" << cdf << "\t" << log(1.0-cdf) << endl;
	}
	cout << endl << endl;
	cout << endl << "2. Computing between cluster similarities:" << endl;
	cout <<         "******************************************" << endl << endl;
	computeDifferentClusterSimilarities(nonClusterList, differentClusterSimilarities, simType);
	cout << "Histogram of similarities:" << endl << fixed << setprecision(3);
	cdf=0.0;
	for (size_t i=0; i<differentClusterSimilarities[0].size(); i++)
	{
		cdf+=differentClusterSimilarities[0][i];
		cout << i*0.01 << "\t" << differentClusterSimilarities[0][i] << "\t" << cdf << "\t" << log(1.0-cdf) << endl;
	}
	cout << endl << endl;

	cout << endl << "3. ROC curves:" << endl;
	cout <<         "*************" << endl;
	cout << "\t";
	const int maxExp = static_cast<int>(differentClusterSimilarities.size());
	int n=1;
	for (int i=0; i<maxExp; i++)
	{
		cout << n << "\t\t";
		n*=10;
	}
	cout << endl;
	cout << "Sim";
	for (int i=0; i<maxExp; i++)
		cout << "\tPrec\tRecall";
	cout << endl;

	vector<double> cumDiffs(maxExp,0.0);
	vector<double> cumSames(maxExp,0.0);
	vector< vector<double> > precisions(maxExp);
	vector< vector<double> > recalls(maxExp);

	assert(sameClusterSimilarities.size() == differentClusterSimilarities[0].size());

	vector<double> factors(20,1.0);
	factors[0] = 0.0;
	factors[1] = 10.0/11.0;
	factors[2] = 100.0/101.0;
	factors[3] = 1000.0/1001.0;
	for (size_t i=0; i<sameClusterSimilarities.size(); i++)
	{
		cout << fixed << setprecision(2) << i*0.01 << "\t";
		for (int d=0; d<maxExp; d++)
		{
			cumSames[d] += sameClusterSimilarities[i];
			cumDiffs[d] += differentClusterSimilarities[d][i];

			double sameAbove = 1.0 - cumSames[d];
			double diffAbove = 1.0 - cumDiffs[d];
			double p,r;
			if (sameAbove<=0.0)
			{
				if (i<10)
				{
					p=0.0;
					r=1.0;
				}
				else
				{
					p=1.0;
					r=0.0;
				}
			}
			else
			{
				p = (sameAbove-diffAbove*factors[d]) /(sameAbove+diffAbove);
				if (p<0.0)
					p=0.0;
				r = sameAbove;
			}
			if (p<0.0)
				p=0.0;
			if (p>1.0)
				p=1.0;
			if (r<0.0)
				r=0.0;
			if (r>1.0)
				r=1.0;
			cout << fixed << setprecision(4) << p << "\t" << r << "\t";
			precisions[d].push_back(p);
			recalls[d].push_back(r);
		}
		cout << endl;
	}
	cout << "AUC=\t";
	for (int d=0; d<maxExp; d++)
		cout << computeRocAuc(precisions[d], recalls[d]) << "\t\t";
	cout << endl;
	cout << endl << "----------------------" << endl;
}


void MsClusterAlgorithm::computeSameClusterSimilarities(const string& clusterList,
														vector<double>& distribution,
														int simType) const
{
	cout << "Scanning input files...";
	const Config* config = model_->get_config();
	SpectraAggregator sa;
	sa.initializeFromTextFile(clusterList.c_str(), config);

	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();
	cout << "... done." << endl;

	vector<Cluster> clusters(sl.getNumHeaders());
	map<Annotation, vector<size_t> > annotationMap;

	cout << "Found " << sl.getNumHeaders() << " spectra." << endl;
	cout << "Reading spectra... ";
	size_t tenth = sl.getNumHeaders()/10;
	// create singleton clusters for all spectra
	for (size_t i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
		PeakList pl;
		if (pl.readPeaksToLocalAllocation(sa, header) == 0)
			continue;
		pl.initializePeakList(config, true);

		clusters[i].createNewCluster(i, header, pl.getPeaks(), pl.getNumPeaks());
		Annotation ann;
		ann.peptideStr = header->getPeptideStr();
		ann.charge	   = header->getCharge();

		assert(ann.peptideStr.length()>0 && ann.charge>0);

		annotationMap[ann].push_back(i);

		if (i>0 && i % tenth == 0)
		{
			cout << i/tenth << " ";
			cout.flush();
		}
	}
	cout << "  ...done." << endl;

	// compare all spectra from the same peptide/charge
	distribution.clear();
	distribution.resize(101,0);

	cout << "computing similarities... ";
	map<Annotation, vector<size_t> >::const_iterator it;
	size_t s=0;
	for (it=annotationMap.begin(); it != annotationMap.end(); it++)
		s++;
	tenth = s/10;
	size_t counter=0;
	size_t numPeps=0;
	

	const double startTime = time(NULL);
	for (it=annotationMap.begin(); it != annotationMap.end(); it++)
	{
		counter++;
		if (counter % tenth == 0)
		{
			cout << counter/tenth << " ";
			cout.flush();
		}
		const vector<size_t>& idxs = it->second;
		if (idxs.size() < 2)
			continue;

		numPeps++;
		for (size_t i=1; i<idxs.size(); i++)
			for (size_t j=0; j<i; j++)
			{
				float similarity = -1.0;

				if (simType == 0) 
				{
					/*similarity = computeSimilarity(
								clusters[idxs[i]].getDistancePeaks(),
												clusters[idxs[j]].getDistancePeaks(),
												Cluster::getPeakIndexToleranceAsInt());*/
					// NP3 GOT sim computation
					similarity = computeSimilarity(
                    								clusters[idxs[i]].getDistancePeaks(),
                    								clusters[idxs[j]].getDistancePeaks(),
                    												Cluster::getPeakIndexToleranceAsInt());
				}
				else if (simType == 1)
				{
					similarity=sharedPeakProportion(
								clusters[idxs[i]].getDistancePeaks(),
												clusters[idxs[j]].getDistancePeaks(),
												Cluster::getPeakIndexToleranceAsInt());
				}
				else if (simType== 2)
				{
					similarity=computeSimilarityMixed(
								clusters[idxs[i]].getDistancePeaks(),
												clusters[idxs[j]].getDistancePeaks(),
												Cluster::getPeakIndexToleranceAsInt());
				}
		
				if (similarity >=0.0)
				{
					const size_t idx = static_cast<size_t>(100.0 * similarity + 0.5);
					distribution[idx]++;
				}
			}
	}
	cout << " done." << endl;

	double total=0;
	for (size_t i=0; i<distribution.size(); i++)
		total+=distribution[i];

	cout << "computed inner cluster distribution using " << numPeps 
		<< " unique peptides/charges (" << fixed << setprecision(0) << total << " pairwise similarities)." << endl;

	cout << endl << "*** 1\tTime for computing distances = " <<time(NULL) - startTime << endl << endl;


	// normalize
	for (size_t i=0; i<distribution.size(); i++)
		distribution[i]/=total;
}




void MsClusterAlgorithm::computeDifferentClusterSimilarities(
										 const string& nonClusterList,
										 vector< vector<double> >& distributions,
										 int simType) const
{
	const int maxExp = 7;
	cout << "Scanning input files...";
	const Config* config = model_->get_config();
	SpectraAggregator sa;
	sa.initializeFromTextFile(nonClusterList.c_str(), config);

	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();
	cout << "... done." << endl;

	vector<Cluster> clusters(sl.getNumHeaders());
	map<Annotation, size_t> annotationMap;

	cout << "Found " << sl.getNumHeaders() << " spectra." << endl;
	cout << "Reading spectra... ";
	size_t tenth = sl.getNumHeaders()/10;

	// create singleton clusters for all spectra
	size_t numPeps=0;
	for (size_t i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
		Annotation ann;
		ann.peptideStr = header->getPeptideStr();
		ann.charge	   = header->getCharge();
		if (ann.charge == 0 || ann.peptideStr.length() == 0)
		{
			cout << header->getTitle() << "\t" << ann.peptideStr << "\t" << ann.charge << endl;
			continue;
		}

		assert(ann.peptideStr.length()>0 && ann.charge>0);

		if (annotationMap.find(ann) != annotationMap.end())
			continue;

		PeakList pl;
		if (pl.readPeaksToLocalAllocation(sa, header) == 0)
			continue;
		pl.initializePeakList(config, true);
		clusters[numPeps].createNewCluster(numPeps, header, pl.getPeaks(), pl.getNumPeaks());

		annotationMap[ann]=numPeps++;

		if (i>0 && i % tenth == 0)
		{
			cout << i/tenth << " ";
			cout.flush();
		}
	}
	cout << " ...done." << endl;


	cout << "Read " << numPeps << " unique peptides/charges (from a total of " 
		<< sl.getNumHeaders() << ")" << endl;

	cout << "Computing similarities (using a window of " << window_ << " m/z)..." << endl;
	clusters.resize(numPeps);
	sort(clusters.begin(),clusters.end());

	distributions.clear();
	distributions.resize(maxExp, vector<double>(101,0.0));

	tenth = numPeps/10;
	double numPairs=0;
	vector< vector<float> > maxSims(7);
	double numSameOrg = 0;
	const double startTime = time(NULL);
	for (size_t i=0; i<numPeps; i++)
	{
		string orgI = clusters[i].getHeader()->getTitle().substr(0,6);
		float maxSimilarity=-1.0;
		size_t j;
		for (j=i+1; j<numPeps; j++)
		{
			mass_t massDiff = (clusters[j].getClusterMOverZ()-clusters[i].getClusterMOverZ());
			if (massDiff<0.1)
				continue;
			if (massDiff>3.2)
				break;
			if (j-i>10000)
				break;

			string orgJ = clusters[j].getHeader()->getTitle().substr(0,6);
			if (orgI == orgJ)
			{
				numSameOrg++;
				continue;
			}

			float similarity=-1.0;
			
			if (simType == 0) 
			{
			    // NP3 GOT sim computation
				similarity = computeSimilarityMixed(
							clusters[i].getDistancePeaks(),
							clusters[j].getDistancePeaks(),
							Cluster::getPeakIndexToleranceAsInt());
			}
			else if (simType == 1)
			{
				similarity=sharedPeakProportion(
							clusters[i].getDistancePeaks(),
							clusters[j].getDistancePeaks(),
							Cluster::getPeakIndexToleranceAsInt());
			}
			else if (simType == 2)
			{
			    // NP3 GOT sim computation
				similarity=computeSimilarity(
							clusters[i].getDistancePeaks(),
							clusters[j].getDistancePeaks(),
							Cluster::getPeakIndexToleranceAsInt());
			}

			if (similarity>= 0)
			{
				maxSims[0].push_back(similarity);
				const size_t idx = static_cast<size_t>(100.0 * similarity + 0.49999);
				distributions[0][idx]++;
			}

		}

		
		numPairs+= (j-i-1);

		if (i>0 && i % tenth == 0)
		{
			cout << i/tenth << " " ;
			cout.flush();
		}
	}
	cout << " ... done." << endl;
	cout << "Skipped " << numSameOrg << " pairs from same org." << endl;
	cout << endl << "*** 2 \tTime for computing distances = " <<time(NULL) - startTime << endl << endl;

	double total=0;
	for (size_t i=0; i<distributions[0].size(); i++)
		total+=distributions[0][i];

	for (size_t i=0; i<distributions[0].size(); i++)
		distributions[0][i] /= total;

	cout << "computed between cluster distribution using " << numPeps 
		<< " unique peptides/charges  (" << fixed << setprecision(0) << total << " pairwise similarities)." << endl;

	cout << "H = 0, total = " << total << endl;
	// compute distributions
	for (int h=1; h<maxExp; h++)
	{
		maxSims[h].clear();
		const size_t n = maxSims[h-1].size();
		size_t idx = 0;
		while (idx<n)
		{
			const size_t maxI = (idx + 10 < n ? 10 : n - idx);
			float maxSim = -1.0;
			for (size_t i=0; i<maxI; i++)
				if (maxSims[h-1][idx+i]>maxSim)
					maxSim = maxSims[h-1][idx+i];
			maxSims[h].push_back(maxSim);
			idx += 10;
		}
		int d=h;
		size_t p=1;
		while (maxSims[d].size()<1000 && d>=0)
		{
			d--;
			p*=10;
		}
		if (d<0)
			error("Bad histogram!");

		vector<float>& histVec = maxSims[d];
		const size_t maxI=histVec.size() - p;
		double total=0;
		for (size_t i=0; i<maxI; i++)
		{
			float maxSim = -1.0;
			for (size_t j=0; j<p; j++)
				if (histVec[i + j] > maxSim)
					maxSim = histVec[i+j];
			if (maxSim>0)
			{
				const size_t idx = static_cast<size_t>(100.0 * maxSim + 0.5);
				distributions[h][idx]++;
				total++;
			}
		}

		for (size_t i=0; i<distributions[h].size(); i++)
			distributions[h][i]/=total;
		
		cout << "H = " <<  h << ", total = " << total << endl;
	}
}



// performs clutering on a small subset of the clusters
// doesn't use any of the data structures expect for clusterAlloc_
// doesn't output anything
void MsClusterAlgorithm::simulateClusteringOnSet(const vector<size_t>& positions)
{
	for (int round=0; round<numRounds_; round++)
	{
		for (size_t i=0; i<positions.size()-1; i++)
		{
			const size_t pos1 = positions[i];
			if (! data_.clusterAlloc_[pos1].getIndInPlay())
				continue;

			for (size_t j=i+1; j<positions.size(); j++)
			{
				const size_t  pos2 = positions[j];
				if (! data_.clusterAlloc_[pos2].getIndInPlay())
					continue;

                // NP3 GOT sim computation
				const float similarity = computeSimilarity(
								data_.clusterAlloc_[pos1].getDistancePeaks(),
								data_.clusterAlloc_[pos2].getDistancePeaks(),
								Cluster::getPeakIndexToleranceAsInt());

				if (similarity == -1.0)
				{
					cout << endl << "With " << pos1 << " - " << pos2 << endl;
					cout << "Size: " << data_.clusterAlloc_[pos1].getClusterSize() << " + " <<
						data_.clusterAlloc_[pos2].getClusterSize() << endl;
				}

				if (similarity > similarityThresholds_[round])
				{

				//	cout << round << " sim: " << similarity << "  " << pos1 << " gets " << pos2 << endl;
					data_.joinClusters(&data_.clusterAlloc_[pos1], &data_.clusterAlloc_[pos2], 6);
				}
			}
		}
	}
}



void MsClusterAlgorithm::benchmarkConsensus(const AllScoreModels* model,
											const string& clusterList,
											float similarity,
											int   numRounds)
{
	model_	 = model;
	window_	 = 2.5;
	data_.runningClusterIdx_ = 0;
	Cluster::setTolerances(model_->get_config()->getTolerance());

	const size_t maxIdx = computeMzIndex(5000.0, Cluster::getPeakIndexTolerance(), 0.0)+2;
	data_.newClusterIdxsLists_.resize(maxIdx,0);
	for (size_t i=0; i<data_.newClusterIdxsLists_.size(); i++)
		data_.newClusterIdxsLists_[i] = data_.clusterIdxAllocation_.allocateVector();
	data_.peakWeightTable_.initWeights(64, 0.15);
	data_.config_ = model->get_config();


	cout << "Scanning input files...";
	const Config* config = model_->get_config();
	SpectraAggregator sa;
	sa.initializeFromTextFile(clusterList.c_str(), config);

	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();
	cout << "... done." << endl;

	setSimilarityThresholds(similarity, numRounds);
	data_.allocateMemory(0.7);

	if (data_.clusterAlloc_.size() < sl.getNumHeaders())
	{
		cout << "Need room for " << sl.getNumHeaders() << " spectra." << endl;
		cout << "Allocated room for " << data_.clusterAlloc_.size() << " spectra." << endl;
		error("Not enough memory allocated!");
	}

	vector<Cluster>& clusters = data_.clusterAlloc_;

	vector<float>	explainedIntensity(sl.getNumHeaders(),0.0);
	vector<int>		numBy(sl.getNumHeaders(),0);
	map<Annotation, vector<size_t> > annotationMap;

	cout << "Found " << sl.getNumHeaders() << " spectra." << endl;
	cout << "Reading spectra... ";
	const size_t tenth = sl.getNumHeaders()/10;
	const int bIdx = config->get_frag_idx_from_label("b");
	const int yIdx = config->get_frag_idx_from_label("y");

	// create singleton clusters for all spectra
	size_t numSpectraRead=0;
	for (size_t i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
		PeakList pl;
		if (pl.readPeaksToBuffer(sa, header, &data_.peakAlloc_[data_.nextPeakPos_])==0)
			continue;

		pl.initializePeakList(config, true);
		data_.nextPeakPos_ += pl.getNumPeaks();
	
		const Peak* spectrumPeaks = pl.getPeaks();
		clusters[numSpectraRead].createNewCluster(numSpectraRead, header, spectrumPeaks, pl.getNumPeaks());
		clusters[numSpectraRead].setConfig(config);

		Annotation ann;
		ann.peptideStr = header->getPeptideStr();
		ann.charge	   = header->getCharge();

		assert(ann.peptideStr.length()>0 && ann.charge>0);

		annotationMap[ann].push_back(numSpectraRead);
		if (i>0 && i % tenth == 0)
		{
			cout << i/tenth << " ";
			cout.flush();
		}

		AnnotatedSpectrum as(pl);	
		as.annotate_spectrum(as.get_true_mass_with_19());
		explainedIntensity[numSpectraRead] = as.get_explianed_intensity();
		numBy[numSpectraRead] = as.get_num_observed_frags(bIdx) + as.get_num_observed_frags(yIdx);
		numSpectraRead++;

	}
	cout << " ...done." << endl;
	data_.nextClusterPos_ = numSpectraRead;


	
	// collect stats
	const size_t tsize = 201;
	double         numSpectraInclusters=0;
	vector<double> numClustersPerSize(tsize,0);
	vector<double> numSetsPerSize(tsize,0);

	vector<double> sumConsensusExplained(tsize,0);
	vector<double> sumConsensusBy(tsize,0);
	vector<double> sumBestExplained(tsize,0);
	vector<double> sumBestBy(tsize,0);
	vector<double> sumSqsExplained(tsize,0);
	vector<double> sumSqsBy(tsize,0);
	vector<double> numConsensus(tsize,0);

	vector<double> sumAvgSingletonExplained(tsize,0);
	vector<double> sumAvgSingletonBy(tsize,0);
	vector<double> numSingleton(tsize,0);



	map<Annotation, vector<size_t> >::iterator it;
	for (it=annotationMap.begin(); it != annotationMap.end(); it++)
	{
		const vector<size_t>& idxs = it->second;
		const size_t setSize = idxs.size();
		size_t setBin = setSize;
		if (setBin>30 && setBin<100)
		{
			setBin /=10;
			setBin *=10;
		}
		else if (setBin>100)
		{
			setBin=200;
		}

		//cout << "SET of " << setSize << endl;
		simulateClusteringOnSet(idxs);
		
		int numSpectra=0;
		int numClusters=0;
		for (size_t i=0; i<setSize; i++)
		{
			Cluster& cluster = data_.clusterAlloc_[idxs[i]];
			if (cluster.getIndInPlay())
				numClusters++;

			if (cluster.getIndInPlay() && cluster.getClusterSize()>1 && 
				cluster.getSingletonVector() && cluster.getSingletonVector()->getSize()>0)
			{
				size_t clusterSize =  cluster.getClusterSize();
				numSpectra+=clusterSize;

				size_t clusterBin = clusterSize;
				if (clusterBin>30 && clusterBin<100)
				{
					clusterBin /=10;
					clusterBin *=10;
				}
				else if (clusterBin>100)
				{
					clusterBin=200;
				}

				numSingleton[clusterBin]++;
				sumAvgSingletonExplained[clusterBin]+=explainedIntensity[idxs[i]];
				sumAvgSingletonBy[clusterBin]+=numBy[idxs[i]];

				float bestExplainedIntensity = explainedIntensity[idxs[i]];
				float bestBy				 = numBy[idxs[i]];
			//	float bestSqs				 = cluster.getHeader()->getSqs();


				data_.makeConsensus(&cluster);

				AnnotatedSpectrum as(cluster);
				as.annotate_spectrum(as.get_true_mass_with_19());

				numConsensus[clusterBin]++;
				sumConsensusExplained[clusterBin]+=as.get_explianed_intensity();
				sumConsensusBy[clusterBin]+=as.get_num_observed_frags(bIdx) + 
											 as.get_num_observed_frags(yIdx);

				const clusterIdx_t* singletons = cluster.getSingletonVector()->getElementsPtr();
				const size_t numSingletons     = cluster.getSingletonVector()->getSize();

				for (size_t j=0; j<numSingletons; j++)
				{
					numSingleton[clusterBin]++;
					sumAvgSingletonExplained[clusterBin]+=explainedIntensity[singletons[j]];
					sumAvgSingletonBy[clusterBin]+=numBy[singletons[j]];

					if (explainedIntensity[singletons[j]]>bestExplainedIntensity )
						bestExplainedIntensity = explainedIntensity[singletons[j]];
					if (numBy[singletons[j]]>bestBy)
						bestBy = numBy[singletons[j]];
				}

				sumBestExplained[clusterBin]+=bestExplainedIntensity;
				sumBestBy[clusterBin]+=bestBy;
			}
		}

		numSpectraInclusters+=numSpectra;
		numClustersPerSize[setBin]+=numClusters;
		numSetsPerSize[setBin]++;
	}

	// report
	cout << "Spectra in clusters: " << static_cast<int>(numSpectraInclusters) << " from " << numSpectraRead << " ( "
		<< fixed << setprecision(4) << numSpectraInclusters/static_cast<double>(numSpectraRead) << ")" <<endl;

	cout << endl << "Average number of clusters per peptide: " << endl;
	cout << "Size\t#sets\tAvg clusters" << endl;
	for (size_t i=0; i<tsize; i++)
		if (numSetsPerSize[i]>0)
			cout << i << "\t" << static_cast<int>(numSetsPerSize[i]) << "\t" 
				 << numClustersPerSize[i]/static_cast<double>(numSetsPerSize[i]) << endl;
	
	cout << endl << "Comparison of expalined intensity and #b/y: " << endl;
	cout << "Size\t#clusts\tclustEI\tclustBY\tsingEI\tsingBY\tdelEI\tdelBY\tdBstEI\tdBstBY" << endl;
	for (size_t i=0; i<tsize; i++)
	{
		if (numConsensus[i]>0)
		{
			double clusEI = sumConsensusExplained[i] / numConsensus[i];
			double clusBY = sumConsensusBy[i] / numConsensus[i];
			double singEI = sumAvgSingletonExplained[i] / numSingleton[i];
			double singBY = sumAvgSingletonBy[i] / numSingleton[i];
			double bestEI = sumBestExplained[i] / numConsensus[i];
			double bestBY = sumBestBy[i]	   / numConsensus[i];
			cout << i << "\t" << static_cast<int>(numConsensus[i])  << "\t" << clusEI << "\t" << clusBY
				<< "\t" << singEI << "\t" << singBY << "\t" << (clusEI-singEI) << "\t" << (clusBY-singBY) 
				<< "\t" << clusEI-bestEI << "\t" << clusBY-bestBY << endl;
		}
	}
}






void MsClusterAlgorithm::benchmarkPairs(const AllScoreModels* model,
										const MsParameterStruct* params)
{
	model_	 = model;
	Cluster::setTolerances(model_->get_config()->getTolerance());

	if (params->spectraListToLoad.length() == 0 )
		error("Must supply list of files to load!");

	// the benchmark tests the identifications at different peptide density levels
	const Config* config = model->get_config();

	SpectraAggregator loadSa;
	vector<Cluster> loadClusters;
	map<Annotation, vector<size_t> > loadAnns;
	readAnnotatedSpectraIntoClusters(params->spectraListToLoad, config, loadSa, loadClusters, loadAnns);
	
	vector<int> squareCounts(NUM_PEAKS_FOR_HEURISTIC,0);
	vector<int> cornerCounts(NUM_PEAKS_FOR_HEURISTIC,0);
	vector< vector<int> > cellCounts(NUM_PEAKS_FOR_HEURISTIC, vector<int>(NUM_PEAKS_FOR_HEURISTIC,0));
	int numPairs=0;
	

	map<Annotation, vector<size_t> >::const_iterator it;
	for (it=loadAnns.begin(); it != loadAnns.end(); it++)
	{
		if (it->second.size()<2)
			continue;

		const vector<size_t>& idxs = it->second;
		for (size_t i=1; i<idxs.size(); i++)
		{
			for (size_t j=0; j<i; j++)
			{
			    // NP3 GOT sim computation
				float similarity = computeSimilarity(loadClusters[idxs[i]].getDistancePeaks(),
													 loadClusters[idxs[j]].getDistancePeaks(),
													 Cluster::getPeakIndexToleranceAsInt());
				if (similarity < 0.3)
					continue;

				const unsigned int* top1 = loadClusters[idxs[i]].getTopPeakIdxs();
				const unsigned int* top2 = loadClusters[idxs[j]].getTopPeakIdxs();

				// test different options for overlap (allow idxs to be +- 1 away)
				vector< vector<bool> > matches(NUM_PEAKS_FOR_HEURISTIC, vector<bool>(NUM_PEAKS_FOR_HEURISTIC));
				for (size_t k=0; k<NUM_PEAKS_FOR_HEURISTIC; k++)
					for (size_t l=0; l<NUM_PEAKS_FOR_HEURISTIC; l++)
						matches[k][l] = (top1[k]==top2[l]-2 || top1[k]==top2[l]-1 || 
										 top1[k]==top2[l] || top1[k]==top2[l]+1 || top1[k]==top2[l]+2);

				numPairs++;
				for (size_t k=0; k<NUM_PEAKS_FOR_HEURISTIC; k++)
					for (size_t l=0; l<NUM_PEAKS_FOR_HEURISTIC; l++)
						if (matches[k][l])
							cellCounts[k][l]++;

				// square
				size_t k;
				for (k=0; k<NUM_PEAKS_FOR_HEURISTIC; k++)
				{
					bool good=false;
					for (size_t l=0; l<=k; l++)
						if (matches[k][l] || matches[l][k])
							good=true;
					if (good)
						break;
				}
				for ( ; k<NUM_PEAKS_FOR_HEURISTIC; k++)
					squareCounts[k]++;

				// corner
				for (k=0; k<NUM_PEAKS_FOR_HEURISTIC; k++)
				{
					bool good=false;
					for (size_t l=k; l<NUM_PEAKS_FOR_HEURISTIC; l++)
						if (matches[k][l] || matches[l][k])
							good=true;
					if (good)
						break;
				}
				for ( ; k<NUM_PEAKS_FOR_HEURISTIC; k++)
					cornerCounts[k]++;
			}
		}
	}

	

	// report results
	cout << "Statistics from " << numPairs << " pairs: " << endl;
	cout << fixed << setprecision(4);
	for (size_t i=0; i<NUM_PEAKS_FOR_HEURISTIC; i++)
		cout << "\t" << i+1;
	cout << endl;
	for (size_t i=0; i<NUM_PEAKS_FOR_HEURISTIC; i++)
	{
		cout << i+1;
		for (size_t j=0; j<NUM_PEAKS_FOR_HEURISTIC; j++)
			//cout << "\t" << cellCounts[i][j];
			cout << "\t" << static_cast<float>(cellCounts[i][j])/static_cast<float>(numPairs);
		cout << endl;
	}
	cout << endl << "Squares: " << endl;
	for (size_t i=0; i<NUM_PEAKS_FOR_HEURISTIC; i++)
		cout << i+1 << " : " << static_cast<float>(squareCounts[i])/static_cast<float>(numPairs) << endl;

	cout << endl << "Corners: " << endl;
	for (size_t i=0; i<NUM_PEAKS_FOR_HEURISTIC; i++)
		cout << i+1 << " : " << static_cast<float>(cornerCounts[i])/static_cast<float>(numPairs) << endl;


}


struct PepEntry {
	bool operator< (const PepEntry& rhs) const
	{
		return (mz<rhs.mz);
	}
	string pep;
	int charge;
	mass_t mz;
	size_t idx;
};	


void  MsClusterAlgorithm::benchmarkSimilarityHistogram(const AllScoreModels* model,
													   const MsParameterStruct* params)
{
	model_	 = model;
	Cluster::setTolerances(model_->get_config()->getTolerance());

	if (params->spectraListToLoad.length() == 0 )
		error("Must supply list of files to load!");

	// the benchmark tests the identifications at different peptide density levels
	const Config* config = model->get_config();

	SpectraAggregator sa;
	vector<Cluster> clusters;
	map<Annotation, vector<size_t> > anns;
	readAnnotatedSpectraIntoClusters(params->spectraListToLoad, config, sa, clusters, anns);

	vector<PepEntry> sortedEntries;

	for (size_t i=0; i<clusters.size(); i++)
	{
		const Cluster& clust = clusters[i];
		PepEntry pe;
		pe.pep = clust.getHeader()->getPeptideStr();
		pe.charge = clust.getHeader()->getCharge();
		pe.idx = i;
		pe.mz  = clust.getHeader()->getMOverZ();
		sortedEntries.push_back(pe);
	}

	sort(sortedEntries.begin(),sortedEntries.end());

	size_t minIdx=0;
	size_t maxIdx=0;

	vector<float> sameSims;
	vector<float> diffSims;

	diffSims.clear();
	sameSims.clear();

	cout << setprecision(4);
	for (size_t i=2; i<sortedEntries.size()-3; i++)
	{
		while (minIdx>1 && sortedEntries[i].mz - sortedEntries[minIdx].mz > 5.7)
			minIdx++;
		while (maxIdx<clusters.size()-1 && sortedEntries[maxIdx].mz - sortedEntries[i].mz < 5.7)
			maxIdx++;

		// look for a same peptide
		vector<size_t> pepIdxs;
		for (size_t j=minIdx; j<=maxIdx; j++)
		{
			if (i==j)
				continue;
			if (sortedEntries[j].pep == sortedEntries[i].pep)
				pepIdxs.push_back(j);
		}
	
		if (pepIdxs.size() == 0)
			continue;

		if (i== minIdx)
			continue;

		const size_t compIdx = pepIdxs[static_cast<size_t>(myRandom() * pepIdxs.size())];
		// NP3 GOT sim computation
		const float sim = computeSimilarity(clusters[sortedEntries[i].idx].getDistancePeaks(),
											clusters[sortedEntries[compIdx].idx].getDistancePeaks(), Cluster::getPeakIndexToleranceAsInt());

		assert(sim>=0 && sim<=1.0);
		sameSims.push_back(sim);
		
		
	/*	size_t otherIdx;
		if (minIdx>0 && sortedEntries[minIdx].pep != sortedEntries[i].pep)
		{
			otherIdx = minIdx;
		}
		else if (maxIdx<clusters.size() && sortedEntries[maxIdx].pep != sortedEntries[i].pep)
		{
			otherIdx = maxIdx;
		}
		else
			continue;*/
        // NP3 GOT sim computation
		const float sim2 = computeSimilarity(clusters[sortedEntries[i].idx].getDistancePeaks(),
											 clusters[sortedEntries[minIdx].idx].getDistancePeaks(), Cluster::getPeakIndexToleranceAsInt());

		assert(sim2>=0.0 && sim2<=1.0);
		diffSims.push_back(sim2);
	//	cout << i << "\t" << diffSims.size() << "\t" << sim << "\t" << sim2 << endl;

		
	}


	cout << "computed: " << sameSims.size() << " and " << diffSims.size() << endl;

	vector<double> sameBins(26,0.0), diffBins(26,0.0);
	for (size_t i=0; i<sameSims.size(); i++)
	{
		size_t b_idx=static_cast<size_t>(25.0*sameSims[i]);
		if (b_idx>25)
		{
			cout << i << " Same: " << sameSims[i] << " ==> " << b_idx << endl;
			exit(0);
		}
		sameBins[b_idx]++;
	}
	for (size_t i=0; i<diffSims.size(); i++)
	{
		size_t b_idx=static_cast<size_t>(25.0*diffSims[i]);
		if (b_idx>25)
		{
			cout << i << " Diff: " << diffSims[i] << " ==> " << b_idx << endl;
			exit(0);
		}
		diffBins[b_idx]++;
	}

	vector<double> cdfSame(25,0.0), cdfDiff(25,0.0);
	for (size_t i=0; i<101; i++)
	{
		sameBins[i]/= static_cast<double>(sameSims.size());
		cdfSame[i]+=sameBins[i];
		if (i>0)
			cdfSame[i]+=cdfSame[i-1];

		diffBins[i]/= static_cast<double>(diffSims.size());
		cdfDiff[i]+=diffBins[i];
		if (i>0)
			cdfDiff[i]+=cdfDiff[i-1];

		cout << setprecision(5);
		cout << i*0.04 << "\t" << sameBins[i] << "\t" << diffBins[i] << "\t" << cdfSame[i] << "\t" << cdfDiff[i] << endl;
	}


}
