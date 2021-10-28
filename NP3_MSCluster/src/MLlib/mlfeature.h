#ifndef __MLFEATURE_H__
#define __MLFEATURE_H__

/*! \file mlfeature.h
	\brief Holds the \c MlFeature class and the \c FeatureSet classes.
*/

#ifndef __INCLUDES_H__
#include "mlincludes.h"
#endif


/*! \brief
	Holds the data describing a single feature (i.e., a dimension in the vector space describing the data. 
*/
class MlFeature {
	friend class MlFeatureSet;
public:
	MlFeature() :  name_("Uninitialized"), indHasNonZeroValues_(false), indIsBool_(false) {}

	MlFeature(const string& n) : name_(n) {}

	void addSuffixToFeatureName(const char* suffix);

	string getName() const    { return name_; }
	void   setName(string& s) { name_ = s; }
	bool   getIndIsBool() const { return indIsBool_; }
	void   setIndIsBool(bool b) { indIsBool_=b; }
	
	const vector<size_t>& getGroups() const { return groups_; }
	void  addGroups(const vector<size_t>& groupIdxs);

private:

	string name_;
	
	bool indHasNonZeroValues_;
	bool indIsBool_;

	vector<size_t> groups_;		/// each feature can be tagged as belonging to multiple groups
};




/*! \brief Manages the features that are used in the vector space describing the data.

	There are two types of features that are involved:
	1. basic features - These features are the dimensions of the sample vectors as they appear 
	in the originial input data that is read from the data files.
	2. derived features - These features are created by applying operators to existing features
	(basic features or derived features with a lower index).

	The reason features are separated into two types is to allow models to be easily extendable, i.e.,
	add features to an existing model without having to retrain the whole models or alter the feature
	indexes. This requires the model functions to process both types of features.
*/
class MlFeatureSet {
	friend class MlModel;
public:

	MlFeatureSet() : indWasInitialized_(false), path_(std::string()) {}

	void readFeatureFile(const char* path, bool verbose=false);
	void writeFeatureFile(const char* path, bool verbose=false) const;

	void readFeatureFile(ifstream& ifs);
	void writeFeatureFile(ostream& os) const;

	size_t	getNumBasicFeatures() const		 { return basicFeatures_.size(); }
	size_t  getNumDerivedFeatrues() const	 { return derivedFeatures_.size(); }
	bool	getIndWasInitialized() const { return indWasInitialized_; }
	void	setIndWasInitialized(bool b) { indWasInitialized_=b; }
	const vector<MlFeature>& getFeatures() const { return basicFeatures_; }
	vector<MlFeature>& getFeatures()			 { return basicFeatures_; }
	const MlFeature& getFeature(size_t i) const { return basicFeatures_[i]; }
	void initializeWithDefaultFeatures(size_t numFeatures);

private:

	bool indWasInitialized_;

	string path_;

	vector<MlFeature> basicFeatures_;	/// Features that make up the data's original instance space
	vector<MlFeature> derivedFeatures_;	/// Features created by applying operators 

	map<string, size_t>		 groupMapping_;	/// Maps each group name to a unique index
	vector<string>			 groupNames_;	/// The names of groups features can belong to
	vector< vector<size_t> > groupMembers_; /// The list of features that belong to each group

};

#endif



