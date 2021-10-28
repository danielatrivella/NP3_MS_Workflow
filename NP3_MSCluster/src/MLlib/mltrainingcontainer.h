#ifndef __MLTRAININGCONTAINER_H__
#define __MLTRAININGCONTAINER_H__

#include "mldata.h"
#include "mloperatorsearch.h"

// A structure that holds all possible parameters for a model being trained
class MlTrainingContainer {
public:
	MlTrainingContainer() : 
		trainingAlgorithm_(SMT_LOGISTIC_CG), indClearWeights_(true), indDataWasInitialized_(false),
		indHadInternalError_(false), indNormalTermination_(false), indReadPreviousModel_(true), 
		inputDir_(""), outputDir_(""), inputName_("Model"), outputName_("Model"),
		maxNumIterations_(0), verboseLevel_(10), numGenerationRounds_(MAX_UINT), maxNumSplits_(10), 
		testRatio_(0.0), perplexityDelta_(1E-8), lambda_(0.0), lmBfgsMemorySize_(5),
		bhattacharyyaDistance_(0.1), minRelativeInformationGain_(0.075), labelToConvertTo0_(0),
		trainingDataSet_(0), testDataSet_(0), auxilaryData_(0) 
		{trainingSetPaths_.clear(); testSetPaths_.clear(); }

	void  printParameters(ostream& os=cout) const;

	void performTrainingTestSplit(float ratio);

	// reads and processes the data (including split to training and test sets)
	bool  initialize(bool verbose=false);

	size_t getTrainingAlgorithm() const			  { return trainingAlgorithm_; }
	void   setTrainingAlgorithm(size_t t)		  { trainingAlgorithm_ = t; }
	size_t getFeatureGenerationType()         { return featureGenerationType_; }
	void   setFeatureGenerationType(size_t t) { featureGenerationType_=t; }
	bool   getIndClearWeights() const		  { return indClearWeights_; }
	void   setIndClearWeights(bool b)		  { indClearWeights_ = b; }
	bool   getIndDataWasInitialized()	      { return indDataWasInitialized_; }
	void   setIndDataWasInitialized(bool b)   { indDataWasInitialized_ = b; }
	bool   getIndHadInternalError() const	  { return indHadInternalError_; }
	void   setIndHadInternalError(bool b)     { indHadInternalError_ = b; }
	bool   getIndNormalTermination() const	  { return indNormalTermination_; }
	void   setIndNormalTermination(bool b)    { indNormalTermination_ = b; }
	bool   getIndReadPreviousModel() const	  { return indReadPreviousModel_; }
	void   setIndReadPreviousModel(bool b)    { indReadPreviousModel_ = b; }
	const vector<string>& getTrainingSetPaths() const { return trainingSetPaths_; }
	void  addTrainingSetPath(const char* p)   { trainingSetPaths_.push_back(p); }
	const vector<string>& getTestSetPaths()     const { return testSetPaths_; }
	void  addTestSetPath(const char* p)      { testSetPaths_.push_back(p); }
	const string& getInputDir()        const { return inputDir_; }
	void  setInputDir(const char* p)		  { inputDir_ = p; }
	const string& getOutputDir()  const { return outputDir_; }
	void  setOutputDir(const char* p)       { outputDir_ = p; }
	const string& getInputName()        const { return inputName_; }
	void  setInputName(const char* p)		  { inputName_ = p; }
	const string& getOutputName()  const { return outputName_; }
	void  setOutputName(const char* p)       { outputName_ = p; }
	size_t getMaxNumIterations()        const { return maxNumIterations_; }
	void   setMaxNumIterations(size_t n)      { maxNumIterations_ = n; }
	size_t getVerboseLevel()            const { return verboseLevel_; }
	void   setVerboseLevel(size_t n)          { verboseLevel_ = n; }
	size_t getNumGenerationRounds()		const { return numGenerationRounds_; }
	void   setNumGenerationRounds(size_t g)   { numGenerationRounds_ = g; }
	size_t getMaxNumSplits()			const { return maxNumSplits_; }
	void   setMaxNumSplits(size_t g)		  { maxNumSplits_ = g; }
	double getLambda()					const { return lambda_; }
	void   setLambda(double l)				  { lambda_ = l; }
	size_t getLmBfgsMemorySize()	    const { return lmBfgsMemorySize_; }
	void   setLmBfgsMemorySize(size_t i)	  { lmBfgsMemorySize_ = i; }
	float  getTestRatio()               const { return testRatio_; }
	void   setTestRatio(float r)			  { testRatio_ = r; }
	double getPerplexityDelta()         const { return perplexityDelta_; }
	void   setPerplexityDelta(double d)       { perplexityDelta_ = d; }
	double getBhattacharyyaDistance()	const { return bhattacharyyaDistance_; }
	void   setBhattachryyaDistance(double d)  { bhattacharyyaDistance_ = d; }
	double getMinRelativeInformationGain()	const { return minRelativeInformationGain_; }
	void   setMinRelativeInformationGain(double d)  { minRelativeInformationGain_ = d; }
	int	   getlabelToConvertTo0() const { return labelToConvertTo0_; }
	void   setLabelToConvertTo0(int label) { labelToConvertTo0_ = label; }

	const  MlDataSet* getTrainingDataSet() const  { return trainingDataSet_; }
	MlDataSet* getTrainingDataSet()				  { return trainingDataSet_; }
	void   setTrainingDataSet(MlDataSet* mld)     { trainingDataSet_ = mld; }
	const  MlDataSet* getTestDataSet()  const     { return testDataSet_; }
	MlDataSet* getTestDataSet()				      { return testDataSet_; }
	void   setTestDataSet(MlDataSet* mld)         { testDataSet_ = mld; }
	MlOperatorSearchData* getAuxilarySearchData()		    { return auxilaryData_; }
	void   setAuxilarySearchData(MlOperatorSearchData* ad) { auxilaryData_ = ad; }

	void  setDesiredClassWeights(const vector<weight_t>& weights) { desiredClassWeights_ = weights; }
	const vector<weight_t>&   getDesiredClassWeights() const { return desiredClassWeights_; }
	void  clearDesiredClassWeights() { desiredClassWeights_.clear(); }

	const vector<size_t>& getTrainingIdxs() const { return trainingIdxs_; }
	const vector<size_t>& getTestIdxs()     const { return testIdxs_; }

	void removeTrainingAndTestIdxs(const vector<size_t>& idxsToRemove);

	string getInputPath() const;
	string getOutputPath() const;

private:
	size_t trainingAlgorithm_; // 
	size_t featureGenerationType_; // e.g., FGT_REGRESSION, FGT_RANKING

	bool  indClearWeights_;
	bool  indDataWasInitialized_;
	bool  indHadInternalError_;
	bool  indNormalTermination_;
	bool  indReadPreviousModel_;

	vector<string> trainingSetPaths_;
	vector<string> testSetPaths_;
	string inputDir_;
	string outputDir_;
	string inputName_;
	string outputName_;
	size_t maxNumIterations_;
	size_t verboseLevel_;
	size_t numGenerationRounds_;
	size_t maxNumSplits_;
	
	float  testRatio_;
	double perplexityDelta_;
	double lambda_; // regularization
	size_t lmBfgsMemorySize_; // the number of gradients to remember in LMBFGS training
	double bhattacharyyaDistance_;
	double minRelativeInformationGain_; // The relative decrase in entropy that needs to be observed to make a split worthwhile
	int	   labelToConvertTo0_; // if this is set to a non 0 value, goes to training and test sets and converts
							   // the labels to 0

	MlDataSet*   trainingDataSet_;
	MlDataSet*	 testDataSet_;

	vector<double> desiredClassWeights_; // should re-weight the instances to acheive these class weights

	vector<size_t> trainingIdxs_;
	vector<size_t> testIdxs_;

	MlOperatorSearchData* auxilaryData_;

	bool  initalizeDataSets(bool verbose = true);
};




#endif




