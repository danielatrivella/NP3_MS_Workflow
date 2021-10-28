#ifndef __MLSCOREMODEL_H__
#define __MLSCOREMODEL_H__

/*! @file mlscoremodel.h
	\brief Holds the abstract class for all score models (i.e., models that assign a score to a feature
	vector such as Logistic Regression, RankBoost, etc.)
*/

#include "mldata.h"
#include "mltrainingcontainer.h"

class FeatureSet; // fwd dclr


size_t getModelTypeFormLabel(const char* s);
void   printModelTypes(ostream& os = cout);


/*! @class MlScoreModel
	\brief This is the base class for all scoring models: logistic regression,
	rankboost, etc. It provides the common interface and funcitonality
	that is required in order to work with the rest of the library.
*/
class MlScoreModel {
	friend class MlModel;
public:

	virtual void  calcScores(const MlSample* sample, vector<float>& scores) const =0;
	virtual double calcScoreForClass(const MlSample* sample, int label=0) const =0;
	virtual double calcErrorRate(const MlDataSet* mld, bool verbose) const =0;
	virtual double calcErrorRate(const MlDataSet* mld, const vector<size_t>& idxs, bool verbose) const =0;
	virtual void  trainModel(MlTrainingContainer* mltd) =0;

	virtual bool readModel(ifstream& ifs) =0;
	virtual bool writeModel(ostream& os = cout) const =0;

	virtual void outputModelAnalysis(const MlDataSet* mld, 
							 const vector<size_t>& idxs,
							 const MlFeatureSet* featureSet,
							 ostream& os = cout) const =0;

	bool writeModel(const char* path) const;



	bool getIndWasInitialized() const { return indWasInitialized_; }

	static size_t getIndModelIsBinary(size_t typeIndex);
	
	static MlScoreModel* createDerivedModelObject(size_t typeIndex);

	static MlScoreModel* readAndAllocateModel(const char* path); // reads the model file and returns a pointer of
																// the appropriate derived class type
	static MlScoreModel* readAndAllocateModel(ifstream& ifs);

protected:

	MlScoreModel() : indWasInitialized_(false), trainingType_(MAX_SIZE_T) {}

	bool indWasInitialized_;
	size_t trainingType_;

};


class MlTrainingTerminationDecider {
public:
	MlTrainingTerminationDecider(double tolerance = 1e-6) : logLikelihoodTolerance_(tolerance), 
								bestTrainLL_(MIN_FLOAT), bestTestLL_(MIN_FLOAT), bestTestRound_(0) {}

	bool shouldTerminateTraining(size_t round, double trainLL, double testLL = MIN_FLOAT);

private:
	double logLikelihoodTolerance_;
	double bestTrainLL_;
	double bestTestLL_;
	size_t bestTestRound_;
};



#endif


