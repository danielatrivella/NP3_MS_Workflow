#ifndef __MLDATA_H__
#define __MLDATA_H__

#include "mlincludes.h"
#include "mlauxfun.h"

/*! \file mldata.h
	\breif Holds the general classes for the data (\c MlSample and \c MlDataSet).
*/

class MlFeatureSet; //fwd dclr

/*! @struct MlSample
	\brief Holds the feature vector for a sinlge instance.

	There are two sets of feature-value vectors in each sample:
	1. \c pairs - these are the regular pairs of features and values as they appear in the instance space.
	2. \c derivedPairs - these are feature-value pairs of special features that are created by applying
	operators during the model's training.

	The reason there is a separation into these two sets is to make it easy to extend models with add-on
	features without needing to retrain the models or change indexes of features, etc. This complicates
	things with the implementation of the score models since they now need to train on two sets of features.
*/
struct MlSample {
	MlSample() :  nextOperatorToApply(0) , label(-1), weight(1.0), groupIndex(-1), rankInGroup(-1), 
					tag1(-1), tag2(-1), tag3(-1), floatTag(0.0) {}

	value_t getValue(size_t featureIdx) const;

	void addPair(size_t index, value_t value) { pairs.push_back(IdxVal(index,value)); }

	void addDerivedPair(size_t index, value_t value) { derivedPairs.push_back(IdxVal(index,value)); }

	
	void clearPairs() { pairs.clear(); }

	void print(ostream& os = cout) const;

	bool checkConsistency() const
	{
		size_t i;
		for (i=1; i<pairs.size(); i++)
			if (pairs[i-1].index >= pairs[i].index)
				return false;

		for (i=1; i<derivedPairs.size(); i++)
			if (derivedPairs[i-1].index >= derivedPairs[i].index)
				return false;
		
		return true;
	}

	void clear()
	{
		label  = -1;
		weight = 1.0;
		groupIndex  = -1; 
		rankInGroup = -1;  
		tag1 = -1;
		tag2 = -1; 
		tag3 = -1;  
		floatTag = 0.0;
		pairs.clear();
		derivedPairs.clear();
	}

	/*! \brief writes all the sample information to a stream (except for the group scores)
	*/
	void writeToStream(ostream& os, bool indWriteAdditionalInfo = false) const;

	/*! \brief reads all the sample information from a stream
	*/
	bool readFromStream(istream& is);


	size_t	 nextOperatorToApply; /*! holds the index of the next operator that should be applied to the pairs	
								      assumes that all operators below this index were already applied. */
	int		 label;
	weight_t weight;

	int groupIndex;       /// for testing/training purposes
	int rankInGroup;     /// for testing/training purposes
	int tag1,tag2,tag3; /// extra information that can be stored for debugging
	float floatTag;    /// a value not used as feature

	vector<IdxVal> pairs;		  /// The feature-value pairs that involve features from the instnace space
	vector<IdxVal> derivedPairs; /*! The feature-value pairs that invovle features that are derived by
									 applying operators. */
	vector<double> groupScores; ///
};


/*! @class MlDataSet
	\brief A container that holds the samples.
*/
struct MlDataSet {
public:
	MlDataSet() : totalWeight_(0), numClasses_(0), numPairs_(0), numDerivedPairs_(0), 
				maxFeatureIndex_(0), maxDerivedFeatureIndex_(0),
				indStatsValid_(false), indOutStreamOpen_(false), outStream_(0), buffer_(0), bufferSize_(0), bufferPos_(0) 
	{ classWeights_.clear(); samples_.clear(); }
	~MlDataSet();

	size_t				  getNumSamples()  const { return samples_.size(); }
	
	void  addSamples(const vector<MlSample>& otherSamples);
	void  addSamples(int label, const vector<MlSample>& otherSamples); // adds the samples and sets the labels to the new value

	void  addSample(const MlSample& sample) { samples_.push_back(sample); indStatsValid_=false; }
	const vector<MlSample>& getSamples()        const { return samples_; }
	const MlSample&		    getSample(size_t i) const { return samples_[i]; }
	vector<MlSample>&		getSamples()		{ return samples_; }
	MlSample&				getSample(size_t i) { return samples_[i]; }
	
	weight_t			  getTotalWeight() const {  return totalWeight_; }

	size_t				  getNumClasess()  const { return numClasses_; }

	size_t				  getNumPairs()	   const { return numPairs_; }
	size_t				  getNumDerviedPairs() const { return numDerivedPairs_; }

	size_t				  getNumBasicFeatures() const { return maxFeatureIndex_+1; }

	void				  setMaxFeatureIndedx(size_t i) { maxFeatureIndex_ = i; }
	void				  setMaxDerivedFeatureIndex(size_t i) { maxDerivedFeatureIndex_ = i; }

	const vector<weight_t>& getClassWeights() const { return classWeights_; }

	void	writeWholeDataFile(const char* dataFile);

	void	readDataFile(const char* dataFile, 
						 size_t maxNumSamplesToRead=MAX_SIZE_T, 
						 bool randomSelection=true,
						 bool verbose=false);

	void	initializeOutputFile(const char* path);

	void	writeSampleToBuffer(const MlSample& sample);

	void	closeOutputFile();

	void	reweight(const vector<weight_t>& ratios);

	void	randomlyReduce(double weightRatio);

	void	printDatasetStatistics(ostream& os = cout);
	void	reserve(size_t numSamples) { samples_.reserve(numSamples); }
	void    convertClassLabelsForBinary(int additionalLabelToTreatAsZero=0); // converts all labels >0 to 1
																			 // if given another label, it will convert it to 0 also
	void	tallyClassStatistics();

	void	printFeatureStatistics(size_t featureIdx) const;

	void	outputFeatureReports(const vector<size_t>& idxs, 
								const MlFeatureSet* featureSet, 
								ostream& os = cout) const;

	void	clear();

	
protected:

	vector<weight_t> classWeights_;
	vector<MlSample> samples_;

	weight_t totalWeight_;
	size_t	 numClasses_;

	size_t numPairs_;
	size_t numDerivedPairs_;
	size_t maxFeatureIndex_;
	size_t maxDerivedFeatureIndex_;

	bool   indStatsValid_;

	// used for on the fly writing of the dataset
	bool	indOutStreamOpen_;
	FILE*	outStream_;
	char*	buffer_;
	size_t	bufferSize_;
	size_t	bufferPos_;

	void	flushOutputBuffer();
	
};


#endif

