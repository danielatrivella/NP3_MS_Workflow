#ifndef __MLMAXIMUMENTROPY_H__
#define __MLMAXIMUMENTROPY_H__

#include "mlscoremodel.h"

class MlTrainingContainer;


class MlMaximumEntropyModel : public MlScoreModel {
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
	vector< vector<double> > weights_;

	void learnLMBFGS(MlTrainingContainer* params);

	void learnCG(MlTrainingContainer* params);

	double calcErrorRateWithLogLikelihood(const MlDataSet* mld, 
			    						const vector<size_t>& idxs,
										bool  verbose,
										double *logLikelihood) const;
};

// utiltiy functions that can by optimzation algorithms

void computeGradient(const MlDataSet* ds,
					 const vector<size_t>& idxs,
					 const vector< vector<double> >& w,
					 const vector< vector<double> >& wtx, //
					 double lambda,
					 vector< vector<double> >& g);

double computePosterior(const MlDataSet* ds,
					  const vector<size_t>& idxs,
					  const vector< vector<double> >& wtx,
					  const vector< vector<double> >& qtx,
					  double wtw,
					  double qtw,
					  double qtq,
					  double lambda,
					  double eta);

// for when the weights are exacatly correct
double computePosterior(const MlDataSet* ds,
					  const vector<size_t>& idxs,
					  const vector< vector<double> >& wtx,
					  double wtw,
					  double lambda);

double lineSearch(const MlDataSet* ds,
				   const vector<size_t>& idxs,
				   const vector< vector<double> >& w,    // the point x
				   const vector< vector<double> >& wtx,
				   const vector< vector<double> >& qtx,
				   const vector< vector<double> >& g,    // the gradient
				   const vector< vector<double> >& q,    // q=p, direction along which we want to minimize f(x)
				   value_t lambda);




#endif
