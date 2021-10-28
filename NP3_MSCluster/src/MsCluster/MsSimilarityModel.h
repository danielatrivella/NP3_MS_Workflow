#ifndef __MSSIMILARITYMODEL_H__
#define __MSSIMILARITYMODEL_H__

/*! @file MsSimilarityModel.h
	\brief Holds the class \c MsSimilarityModel for computing the desired similarity threshold for a spectrum and a list of spectra
*/


#include "Cluster.h"

class Config;

class MsSimilarityModel {
public:
	bool readSimilarityModel(const Config* config);
	void writeSimilarityModel(const Config* config) const;
	void trainSimilarityModel(const Config* config, const char* annotatedMgfList);

	void makeSingleDistributionTables(const Config* config, const char* annotatedMgfList);

	float computeMinSimilarityAllowed(size_t numPeaks, size_t numPairs, double maxContaminationProb) const;
	float computePValue(size_t numPeaks, size_t numPairs, float simialrity) const;

private:

	vector< vector<double> > simCdfs_; /// #peaks/0..100 (rounded similarity)

	void demoSimitalrityValues() const;

};



#endif

