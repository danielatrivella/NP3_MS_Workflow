/*!
@file MsClusterAlgorithm.h
\brief contains the class MsClusterAlgorithm
*/

#ifndef __MSCLUSTERALGORITHM_H__
#define __MSCLUSTERALGORITHM_H__

#include "MsClusterDataStorage.h"
#include "MsSimilarityModel.h"
#include "../PepNovo/AllScoreModels.h"

struct MsParameterStruct; // fwd dclr

/*! @class MsClusterAlgorithm
*/
class MsClusterAlgorithm {
public:

	void performClustering(MsParameterStruct* params,
						   const AllScoreModels* model,
						   bool	  indGreedyClusterJoin = true);


	void mergeTwoArchives(MsParameterStruct* params, const AllScoreModels* model);


	void benchmarkSimilarity(const AllScoreModels* model,
							 const string& nonClusterList,
							 const string& clusterList,
							 mass_t window = 2.0,
							 int simType = 0);


	/*! Checks how well the heuristics work for selecting paris of spectra for distance computation 

	The function expects params to include a list of mgf files \c spectraListToLoad (these are identified peptides - 
	should be at least two spectra per peptide).

	@param model The score models class.
	@param params The structure with the run's parameters.*/
	void benchmarkPairs(const AllScoreModels* model,
						const MsParameterStruct* params);

	void benchmarkSimilarityHistogram(const AllScoreModels* model,
					    			  const MsParameterStruct* params);

	void benchmarkConsensus(const AllScoreModels* model,
							const string& clusterList,
							float similarity,
							int   numRounds);



private:

	const AllScoreModels* model_;

	MsSimilarityModel	  simModel_;

	string outDir_;
	string name_;
	int	   datasetIdx_;
	int	   batch_;
	float  availableGb_;
	float  similarity_;		// the main similarity parameter given by user
	mass_t window_;
	int	   numRounds_;
	int	   verboseLevel_;
	bool   indGreedyClusterJoin_;

	vector<float> similarityThresholds_;    // the threshold to use each round

	float	similarityForLargeClusters_;// for large clusters require a higher similarity so they don'e get corrupted 
									   // a large cluster is one with at least LARGE_CLUSTER_SIZE spectra

	float   minSimilarityForAnything_; // if a pair of spectra do not have at least this much similarity,
									   // don't remember their result



	vector<clusterIdx_t> numJoinsPerRound_;// statistics on the number of joins performed during the clustering

	longInt8_t			 totalSpectraComps_;


	MsClusterDataStorage data_;


	// functions for clustering
	void performFirstRoundSpectraComparisons();
	void performSpectraComparisons(int round,  size_t firstPos, size_t lastPos);
	void fillInPlayIndicators();
	void joinAllClustersThatPassThreshold(int round, size_t firstPos, size_t lastPos);


	// functions for similarity benchmarks
	void computeSameClusterSimilarities(const string& clusterList,
										vector<double>& distribution,
										int simType=0) const;
	void computeDifferentClusterSimilarities(const string& nonClusterList,
											 vector< vector<double> >& distributions,
											 int simType=0) const;
	void setSimilarityThresholds(float similarity, int numRounds);

	// performs clutering on a small subset of the clusters
	// doesn't use any of the data structures expect for clusterAlloc_
	// doesn't output anything
	void simulateClusteringOnSet(const vector<size_t>& positions);

};

#endif





