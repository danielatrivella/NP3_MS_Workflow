#ifndef __MLRANKBOOST_H__
#define __MLRANKBOOST_H__

/*! \file  MlRankBoost.h
	\brief Holds all the data structures and classes involved in the implementation of the RankBoost algorithm.
*/

#include "mlscoremodel.h"
#include "mlfeature.h"



/*! @OrderedPair
	\breif holds an ordered pairs of sample indexes (which is used for the feedback function Phi).

	The typical interpretation is that the sample idx1 should be ranked ahead of idx0.
*/
struct OrderedPair {
	OrderedPair() : idx0(MAX_SIZE_T), idx1(MAX_SIZE_T), tag(MIN_INT), weight(1.0) {}
	OrderedPair(size_t i0, size_t i1) : idx0(i0), idx1(i1), tag(MIN_INT), weight(1.0) {}
	OrderedPair(size_t i0, size_t i1, weight_t w) : idx0(i0), idx1(i1), tag(MIN_INT), weight(w) {}
	OrderedPair(size_t i0, size_t i1, int t, weight_t w) : idx0(i0), idx1(i1), tag(t), weight(w) {}

	bool operator< (const OrderedPair& rhs) const
	{
		return ( idx0 < rhs.idx0 ||
			    (idx0 == rhs.idx0 && idx1 < rhs.idx1));
	}
	
	size_t   idx0, idx1;
	int		 tag;
	weight_t weight;
};




/*! @class MlRankBoostDataset
	\brief Holds all the data structures involved in the training of a RankBoost model. This class 
	expands on MlDataSet to offer the common interface and functionality.
*/
class MlRankBoostDataset : public MlDataSet {

public:

	MlRankBoostDataset() : totalPhiWeight_(0.0) {}

	void addPairToPhiVector(size_t idx0, size_t idx1, weight_t w=1.0) { phi_.push_back(OrderedPair(idx0,idx1,w)); totalPhiWeight_+=w; }

	void addPairToPhiVector(const OrderedPair& orderedPair) { phi_.push_back(orderedPair); totalPhiWeight_ += orderedPair.weight; }

	/*! \brief reads and initializes the the data files for training a RankBoost model.

	The data being read is assumed to consist only of basic feature values. The data needs
	to be initialized (and operators need to be applied to each sample to get the full set
	of feature values.

	@param samplePath - the path to the sample files (the format is the same as MlDataSet). 
	The function assumes there is also a file with the same path and a suffix .phi which is 
	the feedback function file for this data.
	@param maxSampleIdx - the maximal number of samples to read (read 0..maxSampleIdx), the default
	is not to restrict the number of samples.

	@return number of samples read successfully
	*/
	size_t readRankBoostDataSet(const char* samplePath, size_t maxSampleIdx = MAX_SIZE_T);

	
	
	/*! \brief writes an MlRankBoostDataset to file for future processing.

	Only the basic features and the feedback function are written.

	@param samplePath - the path where the file goes (file name should end with something like .txt or .dat, etc.)
	The function will also write a file with the suffix phi to that area (this is the feedback function).
	*/
	void   writeRankBoostDataSet(const char* samplePath);



	const vector<OrderedPair>&		 getPhi() const { return phi_; }

	size_t getNumberOfGroups() const { return numberOfGroups_; }
	void   setNumberOfGroups(size_t n) { numberOfGroups_ = n; }

private:


	double totalPhiWeight_;

	size_t numberOfGroups_; /// during training samples can be placed in groups (in which case the phi usually holds samples from the same group)

	vector<OrderedPair> phi_; /// all pairs of samples with weight > 0 (i.e., idx1 should be ranked higher than idx0)
};


/*! \brief Struct used to define a bin in a RankBoost feature function.

values x<=limit are considered in the bin (only the upper limit is used).
*/
struct ValueBin {
	ValueBin() : upperLimit(MAX_VALUE_T), score(0.0) {}
	ValueBin(value_t v, float s): upperLimit(v), score(s) {}
	ValueBin(value_t v) : upperLimit(v), score(0.0) {}

	bool operator< (const ValueBin& rhs) const
	{
		return (upperLimit < rhs.upperLimit);
	}

	value_t upperLimit; /// the highest value that is to be considred as belonging in the bin
	float   score;
};


/*! @class FeatureFunction
	\breif This struct holds a single  feature function used in RankBoost.
	
	Each feature function gives a score to one of the dimensions in the feature vector (one of the 
	feautres of the model's feature set). Each feature function can be:
	- A binary feature f:X->{0,1}. In this case there is no \c abstainScore_, and it is asumed that the
	  score given when x=0 is 0. If x!=0, the score that is given is \c binaryScore_.
	- A real feature f:X->R. In this case there is also a vector of bins designating the score that
	  should be given in each range of values of x. If the feature does not appear in the list, then
	  the sample will get the abstain score.
*/
class FeatureFunction {
public:
	FeatureFunction() : indHasNonZeroValues_(false), indBinaryFunction_(true), indDerivedFeature_(false),
						featureIdx_(MAX_SIZE_T), abstainScore_(0.0), binaryScore_(0.0) {}

	void addBin(value_t maxValue, float score);

	float getAbstainScore() const { return abstainScore_; }

	float getScore(value_t val) const
	{
		if (indBinaryFunction_)
		{
			return (val ? binaryScore_ : 0.0);
		}
		
		assert(bins_.size()>0);
		const vector<ValueBin>::const_iterator it = lower_bound(bins_.begin(), bins_.end(), ValueBin(val));
		assert(it != bins_.end());
		return it->score;
	}

	bool readFeatureFunction(istream& is);

	bool writeFeatureFunction(ostream& os) const;

private:
	
	bool  indHasNonZeroValues_; /// Indicator if this feature can change the score somehow 
	bool  indBinaryFunction_;   /// Indicator if this function is boolean (returns \c binaryScore_ if x != 0, 0 otherwise).
	bool  indDerivedFeature_;   /// The feature for which this function is created is a dervied feature (in the derived feature set)
	size_t featureIdx_;         /// The index of the feautre for which this function was created (regular feture or derived)
	float abstainScore_;        /// score given if a sample does not have this feature (used if indIsBool_ == false)
	float binaryScore_;		    /// score given if indIsBool_ == true and f(x)=1 (if f(x)==0, the score is assumed to be 0).
	vector<ValueBin> bins_;     /// The threshold function for this feature (used for real non-boolean features)
};



/*!
	@class MlRankBoostModel
	\brief holds the feature functions of the RankBoost model.

	This class is an MlScoreModel however, since there is no meaning to classes in the ranking scenario,
	the functions are a bit different. \c calcScores and \c calcScoreForClass are not to be used. Instead
	there is the function \c computeRankingScore.
*/
class MlRankBoostModel : public MlScoreModel {
public:

	// should not be used
	void  calcScores(const MlSample* sample, vector<float>& scores) const { error("should not use calcScores in RankBoost!"); } 	
	
	// should not be used
	double calcScoreForClass(const MlSample* sample, int label=0) const { error("should not use calcScoreForClass in RankBoost!"); return 0.0; } 

	double calcErrorRate(const MlDataSet* mld, bool verbose = false) const;

	double calcErrorRate(const MlDataSet* mld, const vector<size_t>& idxs, bool verbose =false) const;
	
	void  trainModel(MlTrainingContainer* mltd);

	bool readModel(const char* path);

	bool readModel(ifstream& ifs);

	bool writeModel(ostream& os = cout) const;

	void outputModelAnalysis(const MlDataSet* mld, 
							 const vector<size_t>& idxs, 
							 const MlFeatureSet* featureSet, 
							 ostream& os = cout) const;

	float computeRankingScore(const MlSample& sample) const;

private:

	vector<FeatureFunction> featureFunctions_; /// The feature functions created by RankBoost
};



#endif


