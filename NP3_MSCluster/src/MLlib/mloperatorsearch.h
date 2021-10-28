#ifndef __MLOPERATORSEARCHDATA_H__
#define __MLOPERATORSEARCHDATA_H__

/************************************************************************

This class contains data structures and functions that are used in the
process of operator selection (during the initial training).

*************************************************************************/

#include "mloperatorlist.h"



/********************************************************************************
Data structure that contains all temporary info from samples that is needed
in the operator generation phase.
*********************************************************************************/
class MlOperatorSearchData {
	friend class MlOperator;
	friend class MlOperatorList;
	friend class MlModel;
public:
	MlOperatorSearchData() : numClasses_(0), numFeatures_(0), totalSampleWeight_(0.0) {}

	size_t getNumClasses() const { return numClasses_; }
	size_t getNumBasicFeatures() const { return numFeatures_; }
	weight_t getTotalSampleWeight() const { return totalSampleWeight_; }
	
	// reads the information from the dataset into the various data structures in
	// the search data. If read multiple times, it only adds on data from new features
	void readInfoFromDataSet(const MlDataSet* mld, const vector<size_t>& idxs);

	void determineNumberOfValues();

	// puts current values in featureLists_
	void copyFeatureValuesFromDataSet(const MlDataSet* mld, 
									  const vector<size_t>& idxs,
									  size_t startFeatureIdx =0);

	void evaluateUnaryOperatorsOnFeature(size_t featureIdx, 
										 vector<double>& distances,
										 size_t featureSearchType = FGT_REGRESSION) const;

	double findOptimalSplit(size_t featureIdx,
						  vector<value_t>& thresholds,
						  double gainRegularization = -0.5,
						  double minGainForSplit	= 0.1,
						  size_t maxNumSplits		= 5,
						  size_t maxTableLength		= 100,
						  bool   verbose			= true) const;
	
	void printDetailedSummary(MlDataSet* ds,
							  const MlFeatureSet& fs,
							  ostream& os = cout);
							 
private:

	size_t numClasses_;
	size_t numFeatures_;

	weight_t	totalSampleWeight_;

	vector<weight_t> classWeights_; // the weight of each class

	vector<weight_t> sampleWeights_;    // the weight of each sample

	vector< vector<weight_t> > weightsOfFeatures_; // for each feature and class holds the weight of samples
												  // that have that feature

	vector<weight_t> totalFeatureWeights_; // the weights of the samples that have each feature

	vector< vector< vector<IdxVal> > > featureLists_;  // Feaures X classes X list
													  // feature and class holds a sorted 
											         // list of the indexes and vlaues of the samples 
													// in which it occurs

	vector<size_t> numFeatureValues_;	// how many values does the feature display (MAX_UINT means this feature can have an infinite number of values)

	vector<size_t> featureParents_; // list of feature indexes that are necessary for this feature, e.g., the
							 // the index of the variable being split (if the feature holds one of the bins), 
							 // the indexes of variables in a conditional statement (if the feature is the results), etc.


	ListStorage featureIntersectingLists_;

	
	void   countNumValues(size_t featureIdx, const vector<value_t>& thresholds, vector<size_t>& counts) const;

	void createCumulativeWeightTableForFeature(size_t featureIdx, 
											   vector< vector< weight_t > >& cumulativeWeightsTable,
											   vector< value_t >& tableColumnValues,
											   size_t maxTableLength = 100,
											   bool   verbose=false) const;

	void createEntropyTableForFeature(
								const vector< vector<weight_t> >& cumulativeWeightsTable,
								vector< vector<weight_t> >&       entropyTable,
								bool							  verbose = false) const;

};




#endif


