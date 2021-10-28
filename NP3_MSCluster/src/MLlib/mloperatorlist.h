#ifndef __MLOPERATORLIT_H__
#define __MLOPERATORLIT_H__

#include "mlfeature.h"
#include "mldata.h"
#include "mloperatorsearch.h"


typedef enum OPERATOR_TYPE { OT_DROP, OT_INDICATOR, OT_NORMALIZATION, OT_FUNCTION,  
							 OT_SPLIT, OT_CONDITIONAL, OT_NUM_OPERATORS } OPERATOR_TYPE;


typedef enum UnaryOperatorType  { // unary operators, need only one feature to work with
							 UOT_SELF, // output f(x)=x
							 UOT_BOOL, // output f(x)=1 if the feature was found in the sample
							 UOT_LOG,  // output f(x)=log(x) for x>0, otherwise minLogValue
							 UOT_EXP,  // output f(x)=e^x
							 UOT_SQR,  // output f(x)=x^2
							 UOT_SQRT, // output f(x)=x^1/2
							 UOT_NEG,  // output f(x)=-x
							 UOT_ABS,  // output f(x)=|x|
							 UOT_SCALE, // normalize x so according to mean and variance
							 UOT_SPLIT,// split x into several features according to values
									   // if add bool flag is set, each bin also gets a boolean indicator
							 UOT_NUM_OPERATORS,
} UnaryOperatorType;

const value_t minLogValue = -20.0;
const char* const unaryOperatorLabels[]={"SELF", "BOOL",  "LOG",  "EXP", "SQR","SQRT","NEG", "ABS", "SCALE",  "SPLIT"};
const size_t numUnaryOperatorLabels = sizeof(unaryOperatorLabels)/sizeof(char*);

typedef enum ConditionalOperatorType {
							 COT_AND,     // performs AND of a list of features
							 COT_OR,	  // performs OR of a list of features
							 COT_NOT_AND, // returns 1 unless all features are present
							 COT_NOT_OR,  // returns 1 only if all featureas are not present
							 COT_NUM_OPERATORS,
} ConditionalOperatorType;

const char* const conditionalOperatorLabels[]={"AND",  "OR",  "NOT_AND", "NOT_OR"};
const size_t numConditionalOperatorLabels = sizeof(conditionalOperatorLabels)/sizeof(char*);


// these are the values that get outputted in the result feature if the conditional
// operator evaluates to "true". Usually there will be a feature that ouputs 1 in addition
// to the type of value(s) listed below.
typedef enum ConditionalValueType {
									CVT_BOOL,	// output 1.0
									CVT_SUM,	// output the sum of values being operated on
									CVT_PROD,   // output the product of the values being operated on
									CVT_MAX,    // output the max val
									CVT_MIN,	// output tht minimal val
									CVT_AVG,	//
									CVT_COUNT,	// count how many features had a value (good only for OR)
									CVT_MEDIAN,
									CVT_FEATURE, // output the value of a specific feature
									CVT_NUM_VALUES,
} ConditionalOperatorValueType;

const char* const conditionalValueLabels[]={"BOOL","SUM","PROD","MAX","MIN","AVG","COUNT","MEDIAN","FEATURE"};
const size_t numConditionalValueLabels = sizeof(conditionalValueLabels)/sizeof(char*);


const size_t unaryOperatorsForRegression[]  = {UOT_SELF, UOT_LOG, UOT_EXP, UOT_SQR, UOT_SQRT, UOT_ABS};
const size_t numUnaryOperatorsForRegression =  sizeof(unaryOperatorsForRegression)/sizeof(size_t);
const size_t unaryOperatorsForRanking[]     = {UOT_SELF, UOT_BOOL};
const size_t numUnaryOperatorsForRanking    =  sizeof(unaryOperatorsForRanking)/sizeof(size_t);

// If a feature doesn't display at least this number of values, non of the real
// vlaued operatores (e.g., SQR,LOG, etc) will be applied to it
const size_t minNumValuesToBeConsideredReal = 7; 





/*****************************************************************************
An operator performs some type of manipulation of the feature values. There
are several type of operators that can be used:
Syntax     : Description
-------------------------
Drop X             - Drop the feature X
Identity X Y       - if the sample has the feature X, assign 1.0 to the feature y
Function X Y func  - if there is a feature X, perform func on it and write result to feature Y
			 (func comes from a list like sqr(x), sqrt(x), log(x), abs(x), etc., see below)
Normalization X Y mu sigma - If the feature x has a value scale it z=(x-mu)/sigma and write result in Y
Split X <T1..Tn> <Y1..Yn+1> <B1..Bn+1>: If X has a value, assign to a bin i, write result to a Yi bin
									    if Bi is a good address, write a boolean to there
Conditional <X1 X2 Xn> Y cond result-type : perform a boolean condition on X1..X2 example (AND OR)
								  write result in Y (result can be a dunction on the variables
								  like MIN, MAX, SUM, COUNT, etc. see defintions below

Operators are applied to the sample in the order they were read.
The operators are not implemented as derivied classes to save the run-time binding

******************************************************************************/



struct IndicatorOperator {
	IndicatorOperator() : sourceIdx(MAX_UINT), targetIdx(MAX_UINT) {}
	size_t sourceIdx;
	size_t targetIdx;
};

struct FunctionOperator {
	FunctionOperator() : sourceIdx(MAX_UINT), targetIdx(MAX_UINT), type(MAX_UINT) {}
	size_t sourceIdx;
	size_t targetIdx;
	size_t type; // UnaryOperatorType
};

struct NormalizationOperator {
	NormalizationOperator() : sourceIdx(MAX_UINT), targetIdx(MAX_UINT), mu(MAX_FLOAT), sigma(MIN_FLOAT) {}
	size_t sourceIdx;
	size_t targetIdx;
	value_t mu;
	value_t sigma;
};

// a split operator creates n>1 variables (bins) to which the values get sent
struct SplitOperator {
	SplitOperator() : sourceIdx(MAX_UINT) {}
	size_t sourceIdx;
	vector<value_t> thresholds;      // dim n
	vector<size_t> indexesForBinValues;  // dim n+1
	vector<size_t> indexesForBinIndicators;   // dim n+1
};

struct ConditionalOperator {
	ConditionalOperator() : targetIdx(MAX_UINT), indexForBool(MAX_UINT), conditionType(MAX_UINT), resultType(MAX_UINT) {}
	size_t targetIdx;
	size_t indexForBool; 
	size_t conditionType; // ConditionalOperatorType
	size_t resultType;    // ConditionalValueType
	vector<size_t> sourceIdxs;
};

class MlTrainingContainer; // fwd dclr
class MlOperatorSearchData; // fwd dclr

class MlOperatorList {
	friend class MlModel;
public:

	MlOperatorList(): indWasInitialized_(false), operatorFeatureSpaceSize_(0) {}

	bool	getIndWasInitialized() const { return indWasInitialized_; }
	void	setIndWasInitialized(bool b) { indWasInitialized_ = b; }

	bool readOperatorList(const char* path);
	bool readOperatorList(ifstream& ifs);

	bool writeOperatorList(const char* path) const;
	bool writeOperatorList(ostream& os) const; 

	void applyOperatorList(MlSample& sample) const;

	void applyOperatorList(MlDataSet* dataSet, const vector<size_t>& sampleIdxs) const;

	size_t getNumOperators() const { return executionOrder_.size(); }

	// chooses operators that can be useful in separating classes
	// this is done before an initial model is created
	// does not include splits
	// this function is useful for regressions, not for boosting
	void selectInitialUnaryOperators(MlTrainingContainer& mtc, 
									 const MlFeatureSet& featureSet,
									 bool verbose=false);

	// chooses split operators that can be usesful in separating the classes
	// this is done before an intial model is created
	// verbose = 0, minimal chatter
	// verbose = 1, one line summary for each feature
	// verbose > 1, more details
	void selectInitialSplitOperators(MlTrainingContainer& mtc, 
									 const MlFeatureSet& featureSet,
									 bool  addBooleanIndicators,
									 size_t verboseLevel=0);

	void createNewOperators(const MlDataSet& ds,
							const vector<size_t>& designatedIndexes,
							size_t maxNumOperators);

	static size_t determineFeatureGenerationType(size_t trainingType);


private:

	bool				indWasInitialized_;

	size_t				operatorFeatureSpaceSize_; // the maximum index that can be accessed/written to 
	
	mutable vector<value_t>		expandedVector_;



	void createSplitOperator(size_t originalFeatureIndex,
							 const vector<value_t>& thresholds,
							 MlOperatorSearchData* auxilaryData,
							 bool  addIndicators);

	vector<IdxPair>		executionOrder_; // holds pairs type/index in array of operators that need to be applied

	vector<size_t>                drops_;
	vector<IndicatorOperator>     indicators_;
	vector<FunctionOperator>      functions_;
	vector<NormalizationOperator> normalizations_;
	vector<SplitOperator>		  splits_;
	vector<ConditionalOperator>	  conditionals_;

};





#endif

