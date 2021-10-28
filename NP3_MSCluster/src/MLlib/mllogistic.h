#ifndef __MLLOGISTIC_H__
#define __MLLOGISTIC_H__

#include "mlscoremodel.h"

class MlTrainingContainer; // fwd dclr
class MlFeatureSet; // fwd dclr

class MlLogisticModel : public MlScoreModel {
public:
	void  calcScores(const MlSample* sample, vector<float>& scores) const;
	double calcScoreForClass(const MlSample* sample, int label=0) const;
	double calcErrorRate(const MlDataSet* mld, bool verbose = false) const;
	double calcErrorRate(const MlDataSet* mld, const vector<size_t>& idxs, bool verbose =false) const
	{ return calcErrorRateWithLogLikelihood(mld, idxs, verbose, 0); }

	void  trainModel(MlTrainingContainer* mltd);

	bool readModel(const char* path);
	bool readModel(ifstream& ifs);
	bool writeModel(ostream& os = cout) const;

	void outputModelAnalysis(const MlDataSet* mld, 
								 const vector<size_t>& idxs, 
								 const MlFeatureSet* featureSet, 
								 ostream& os = cout) const;

private:
	vector<value_t> weights_;

	void learnCG(MlTrainingContainer* params);

	double	calcErrorRateWithLogLikelihood(const MlDataSet* mld, 
			    						const vector<size_t>& idxs,
										bool  verbose,
										double *logLikelihood) const;
};



#endif


