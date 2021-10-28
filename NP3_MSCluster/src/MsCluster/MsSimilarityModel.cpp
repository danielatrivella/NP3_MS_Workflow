#include "MsSimilarityModel.h"
#include "MsClusterAuxfuns.h"
#include "../PepNovo/SpectraList.h"
#include "../PepNovo/AnnotatedSpectrum.h"

struct mzIdx {
	bool operator< (const mzIdx& rhs) const
	{
		return (mz < rhs.mz);
	}
	mass_t mz;
	size_t idx;
};

void MsSimilarityModel::trainSimilarityModel(const Config* config, const char* annotatedMgfList)
{
	Cluster::setTolerances(config->getTolerance());

	cout << "Training similarity models with tolerance: " << config->getTolerance() << endl;

    SpectraAggregator sa;
	sa.initializeFromTextFile(annotatedMgfList, config);
	sa.setDatasetIdxAccordingToFileName();
	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();

	vector<Cluster> clusters(sl.getNumHeaders());
	const size_t maxIdx = computeMzIndex(4500.0, Cluster::getPeakIndexTolerance(), 0.0);
	vector< vector<size_t> > idxLists(maxIdx);

	Peak* peaks = new Peak[10000];
	size_t numBadClusters=0;
	size_t numInPlay=0;

	cout << "Examining: " << sl.getNumHeaders() << endl;

	for (size_t i=0; i<sl.getNumHeaders(); i++)
	{
		clusters[i].setIndInPlay(false);
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
		if (header->getPeptideStr().length()<2)
			continue;

		if (clusters[i].readPeaksToBuffer(sa, header, peaks) < 7)
			continue;
		clusters[i].initializePeakList(config, true);
		if (! clusters[i].createNewCluster(i, header, peaks, clusters[i].getNumPeaks()))
		{
			numBadClusters++;
			continue;
		}

		clusters[i].setIndInPlay(true);
		numInPlay++;
	}
	sort(clusters.begin(), clusters.end());
	cout << "Read and sorted " << clusters.size() << " spectra, " << numBadClusters << " are bad (" << numInPlay << " are in play)." << endl;

	for (size_t i=0; i<clusters.size(); i++)
	{
		clusters[i].setClusterIdx(i);
		for (size_t j=0; j<NUM_PEAKS_FOR_HEURISTIC; j++)
		{
			unsigned int idx = clusters[i].getTopPeakIdxs()[j];
			if (idx<0 || idx>=maxIdx)
			{
				cout << "bad idx: " << idx << " (" << j << ")" << endl;
				exit(0);

			}
			assert(idx>=0 && idx<maxIdx);
			if (idx>0)
				idxLists[idx].push_back(i);
		}
	}

	vector< vector<double> > simCounts(MAX_NUM_PEAKS_FOR_DISTANCE+2, vector<double>(101, 0));
	size_t minClusterIdx=0;
	size_t maxClusterIdx=0;
	for (size_t i=0; i<clusters.size(); i++)
	{
		const Cluster& clust = clusters[i];
		const mass_t mz      = clusters[i].getClusterMOverZ();
		const int orgIdx = clust.getHeader()->getDatasetIndex();
		const int numDistancePeaks = clust.getDistancePeaks()->numPeaks;
		const string& pepStr = clust.getHeader()->getPeptideStr();

		assert(mz>0.0);
		const mass_t minMz = mz - 17.0;
		const mass_t maxMz = mz + 17.0;

		cout << ">>\t" << i << "\t" << mz;

		while (minClusterIdx<clusters.size() && clusters[minClusterIdx].getClusterMOverZ()<minMz)
			minClusterIdx++;
		while (maxClusterIdx<clusters.size() && clusters[maxClusterIdx].getClusterMOverZ()<maxMz)
			maxClusterIdx++;
		cout << "\t" << minClusterIdx << " - " << maxClusterIdx << "\t";

		vector<size_t> clusterIdxs;
		clusterIdxs.push_back(0);
		for (size_t j=0; j<NUM_PEAKS_FOR_HEURISTIC; j++)
		{
			const unsigned int idx = clusters[i].getTopPeakIdxs()[j];
			if (idx<2)
				continue;

			for (size_t k=idx-2; k<=idx+2; k++)
			{
				const vector<size_t>& list = idxLists[k];
				if (list.size() == 0)
					continue;

				vector<size_t>::const_iterator it = lower_bound(list.begin(),list.end(),minClusterIdx);
				while (it != list.end() && *it<maxClusterIdx)
				{
					const size_t clustIdx = *it++;
					const Cluster& otherClust = clusters[clustIdx];
					if ((otherClust.getHeader()->getDatasetIndex() != orgIdx &&
						fabs(mz-otherClust.getClusterMOverZ())>7.0)&&
						abs(otherClust.getDistancePeaks()->numPeaks - numDistancePeaks) <= ALLOWED_DIFF_IN_NUM_DISTANCE_PEAKS &&
						otherClust.getHeader()->getPeptideStr() != pepStr)
					{
						clusterIdxs.push_back(clustIdx);
					}
				}
			}
		}
		sort(clusterIdxs.begin(),clusterIdxs.end());

		size_t count=0;
		float maxSim = 0.0;
		for (size_t j=1; j<clusterIdxs.size(); j++)
		{
			if (clusterIdxs[j] == clusterIdxs[j-1])
				continue;
			const size_t otherIdx = clusterIdxs[j];
			assert(otherIdx < clusters.size());
			//  NP3 sim computation
            //const float sim = computeSimilarity(clust.getDistancePeaks(), clusters[otherIdx].getDistancePeaks(), Cluster::getPeakIndexToleranceAsInt());
            const float sim = computeSimilarity(clust.getDistancePeaks(), clusters[otherIdx].getDistancePeaks(), Cluster::getPeakIndexToleranceAsInt());
			const size_t bin = static_cast<size_t>(sim * 100);
			assert(bin<=101);
			const int otherNumPeaks = clusters[otherIdx].getDistancePeaks()->numPeaks;
			simCounts[numDistancePeaks][bin]++;
			simCounts[otherNumPeaks][bin]++;

			// smoothing
			if (bin<90)
			{
				simCounts[numDistancePeaks][bin+1]+=0.5;
				simCounts[otherNumPeaks][bin+1]+=0.5;
				simCounts[numDistancePeaks-1][bin+1]+=0.25;
				simCounts[otherNumPeaks-1][bin+1]+=0.25;
				simCounts[numDistancePeaks+1][bin+1]+=0.25;
				simCounts[otherNumPeaks+1][bin+1]+=0.25;
			}
			if (bin>1)
			{
				simCounts[numDistancePeaks][bin-1]+=0.5;
				simCounts[otherNumPeaks][bin-1]+=0.5;
				simCounts[numDistancePeaks-1][bin-1]+=0.25;
				simCounts[otherNumPeaks-1][bin-1]+=0.25;
				simCounts[numDistancePeaks+1][bin-1]+=0.25;
				simCounts[otherNumPeaks+1][bin-1]+=0.25;
			}
			simCounts[numDistancePeaks-1][bin]+=0.5;
			simCounts[otherNumPeaks-1][bin]+=0.5;
			simCounts[numDistancePeaks+1][bin]+=0.5;
			simCounts[otherNumPeaks+1][bin]+=0.5;

			if (sim>maxSim)
				maxSim = sim;
			count++;
		}
		cout <<  "\t" << count << "\t" << maxSim << endl;
	}

	
	simCdfs_.clear();
	simCdfs_.resize(MAX_NUM_PEAKS_FOR_DISTANCE+1, vector<double>(101,1.0 - 1E-12));
	vector<double> sums(MAX_NUM_PEAKS_FOR_DISTANCE+1,0);
	ofstream ofs("dist_new.txt");
	size_t minIdx=0;
	for (size_t i=0; i<=MAX_NUM_PEAKS_FOR_DISTANCE; i++)
	{
		for (size_t j=0; j<101; j++)
			sums[i] += simCounts[i][j];

		if (sums[i]==0)
			continue;
		
		if (minIdx==0)
			minIdx=i;

		// NP3 GOT precision from 3 to 4
		ofs << i << "\t" << sums[i] << scientific << setprecision(4);
		double s=0;
		for (size_t j=0; j<simCounts[i].size(); j++)
		{
			s += simCounts[i][j];
			simCdfs_[i][j] = 1.0 - (s/sums[i]);
			ofs << "\t" << simCdfs_[i][j];
		}
		ofs << endl;
	}
	for (size_t i=0; i<=minIdx; i++)
		simCdfs_[i]=simCdfs_[minIdx+1];

	writeSimilarityModel(config);

	demoSimitalrityValues();
}

void MsSimilarityModel::writeSimilarityModel(const Config* config) const
{
	string path = config->get_resource_dir() + "/" + config->get_model_name() + "_simcdf.txt";
	ofstream ofs(path.c_str());
	if (! ofs.good())
		error("Could not open simialrity file for writing: ",path.c_str());

	ofs << scientific << setprecision(4);
	for (size_t i=0; i<simCdfs_.size() && i<= MAX_NUM_PEAKS_FOR_DISTANCE; i++)
	{
		ofs << i << "\t" << simCdfs_[i].size();
		for (size_t j=0; j<simCdfs_[i].size(); j++)
			ofs << "\t" << simCdfs_[i][j];
		ofs << endl;
	}
	ofs.close();
}

bool MsSimilarityModel::readSimilarityModel(const Config* config)
{
	string path = config->get_resource_dir() + "/" + config->get_model_name() + "_simcdf.txt";
	ifstream ifs(path.c_str());
	if (! ifs.good())
		error("Could not open simialrity file for reading: ",path.c_str());

	simCdfs_.clear();
	simCdfs_.resize(MAX_NUM_PEAKS_FOR_DISTANCE+1, vector<double>(101,1E-12));

	char buffer[8192];
	while (ifs.good())
	{
		ifs.getline(buffer,8192);
		istringstream iss(buffer);
		size_t n=0, size=0;
		iss >> n >> size;
		if (n>0 && n<=MAX_NUM_PEAKS_FOR_DISTANCE)
		{
			simCdfs_[n].resize(size);
			for (size_t i=0; i<size; i++)
				iss >> simCdfs_[n][i];
		}
	}
	ifs.close();

	//demoSimitalrityValues();

	return true;
}


float MsSimilarityModel::computeMinSimilarityAllowed(size_t numPeaks, size_t numPairs, double maxContaminationProb) const
{
	const double minCdfVal = 1.0 - exp(log(1.0 - maxContaminationProb)/static_cast<double>(numPairs));
	const vector<double>& cdfVals = simCdfs_[numPeaks];
	size_t index = 30;
	if (cdfVals[index]>minCdfVal)
	{
		while (index<100 && cdfVals[index]>minCdfVal)
			index++;
	}
	else
	{
		while (index>0 && cdfVals[index]<minCdfVal)
			index--;
		index++;
	}
	return (static_cast<float>(index)*0.01);
}

float MsSimilarityModel::computePValue(size_t numPeaks, size_t numPairs, float similarity) const
{
	assert(similarity>=0.0 && similarity<=1.0);
	//assert(similarity<=1.0);
	const double cdf = simCdfs_[numPeaks][static_cast<size_t>(similarity*100.0)];
	return 1.0 - exp(numPairs * log(1.0 - cdf));
}


void MsSimilarityModel::demoSimitalrityValues() const
{
	const double probs[]={0.25,0.1,0.05,0.01};
	const size_t n_probs=sizeof(probs)/sizeof(double);
	const size_t sizes[]={1,10,100,1000,10000,100000};
	const size_t n_sizes=sizeof(sizes)/sizeof(size_t);

	for (size_t i=0; i<n_probs; i++)
	{
		cout << endl << "Mixure prob: " << setprecision(3) << probs[i] << endl;
		
		cout << "n_peaks" << setprecision(2);
		for (size_t j=0; j<n_sizes; j++)
			cout << "\t" << sizes[j];
		cout << endl;
		for (size_t n =10; n<=MAX_NUM_PEAKS_FOR_DISTANCE; n+=3)
		{
			cout << n;
			for (size_t j=0; j<n_sizes; j++)
				cout << "\t" << computeMinSimilarityAllowed(n, sizes[j], probs[i]);
			cout << endl;
		}
		cout << endl;
	}
//	exit(0);
}


void MsSimilarityModel::makeSingleDistributionTables(const Config* config, const char* annotatedMgfList)
{
	SpectraAggregator sa;
	sa.initializeFromTextFile(annotatedMgfList, config);
	sa.setDatasetIdxAccordingToFileName();
	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();

	vector<Cluster> clusters(sl.getNumHeaders());
	

	Peak* peaks = new Peak[10000];
	size_t numBadClusters=0;
	size_t numInPlay=0;
	for (size_t i=0; i<sl.getNumHeaders(); i++)
	{
		clusters[i].setIndInPlay(false);
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
		if (header->getPeptideStr().length()<2)
			continue;

		if (clusters[i].readPeaksToBuffer(sa, header, peaks) < 7)
			continue;
		clusters[i].initializePeakList(config, true);
		if (! clusters[i].createNewCluster(i, header, peaks, clusters[i].getNumPeaks()))
		{
			numBadClusters++;
			continue;
		}

		clusters[i].setIndInPlay(true);
		numInPlay++;
	}
	sort(clusters.begin(), clusters.end());
	cout << "Read and sorted " << clusters.size() << " spectra, " << numBadClusters << " are bad (" << numInPlay << " are in play)." << endl;

	
	vector<double> sameCounts(51, 0);
	vector<double> diffCounts(51, 0);
	size_t minClusterIdx=0;
	size_t maxClusterIdx=0;
	for (size_t i=0; i<clusters.size(); i++)
	{
		const Cluster& clust = clusters[i];
		const mass_t mz      = clusters[i].getClusterMOverZ();
		const string& pepStr = clust.getHeader()->getPeptideStr();

		assert(mz>0.0);
		const mass_t minMz = mz - 10.0;
		const mass_t maxMz = mz + 10.0;

		cout << ">>\t" << i << "\t" << mz;

		while (minClusterIdx<clusters.size() && clusters[minClusterIdx].getClusterMOverZ()<minMz)
			minClusterIdx++;
		while (maxClusterIdx<clusters.size() && clusters[maxClusterIdx].getClusterMOverZ()<maxMz)
			maxClusterIdx++;
		cout << "\t" << minClusterIdx << " - " << maxClusterIdx << "\t";

		size_t nSame=0, nDiff=0;
		float maxSame = 0.0, maxDiff=0;

		for (size_t idx=minClusterIdx; idx<maxClusterIdx; idx++)
		{
			if (idx == i)
				continue;

			const string& otherStr = clusters[idx].getHeader()->getPeptideStr();
			if (otherStr == pepStr)
			{
			    // GOT NP3 sim computation
				//const float sim = computeSimilarity(clust.getDistancePeaks(), clusters[idx].getDistancePeaks(), Cluster::getPeakIndexToleranceAsInt() );
				const float sim = computeSimilarity(clust.getDistancePeaks(),
				                                    clusters[idx].getDistancePeaks(), Cluster::getPeakIndexToleranceAsInt());
				const size_t bin = static_cast<size_t>(sim * 50);
				sameCounts[bin]++;
				nSame++;
				if (sim>maxSame)
					maxSame = sim;
				continue;
			}

			if (fabs(clusters[idx].getHeader()->getMOverZ()-mz)<7.0)
				continue;
			// NP3 GOT sim computation
			//const float sim = computeSimilarity(clust.getDistancePeaks(), clusters[idx].getDistancePeaks(), Cluster::getPeakIndexToleranceAsInt());
			const float sim = computeSimilarity(clust.getDistancePeaks(),
			                                    clusters[idx].getDistancePeaks(), Cluster::getPeakIndexToleranceAsInt());

			const size_t bin = static_cast<size_t>(sim * 50);
			diffCounts[bin]++;

			nDiff++;
			if (sim>maxDiff)
				maxDiff=sim;
		}
		cout <<  "\t" << nSame << "\t" << maxSame << "\t" << nDiff << "\t" << maxDiff << endl;
	}
	cout << scientific << setprecision(4) << endl;

	double totalSame=0, totalDiff =0;
	for (size_t i=0; i<51; i++)
	{
		totalSame+=sameCounts[i];
		totalDiff+=diffCounts[i];
	}

	cout << "Total same: " << totalSame << endl;
	cout << "Total diff: " << totalDiff << endl;

	for (size_t i=0; i<51; i++)
		cout << i*0.02 << "\t" << sameCounts[i]/totalSame << "\t" << diffCounts[i]/totalDiff << endl;
	
}




