#ifndef __MLMODEL_H__
#define __MLMODEL_H__

#include "mlscoremodel.h"
#include "mlfeature.h"
#include "mloperatorlist.h"

/*******************************************************************
This is the container class for the ML library.
An ML Model has three main componenets:

1) list of features. This is the list of N' feature names that are
   the expansion of the original N features. The list of features is
   stored with the "_fet.txt" suffix.

2) List of operators (that works on the original space of N features).
   The operators exapnd each feature vector to N'>=N dimensions.
   The list of operators is stored with the suffix "_opr.txt".

3) A scoring model, such as logistic regression, ranking, etc.
   The model is designed to score instances in the expanded (N') 
   feature space. Models are stored with the "_scr.txt" suffix

Alternativley, all 3 model files can be concatenated and stored as 
one file with the suffix "_mod.txt".
********************************************************************/


class MlModel {
public:

	MlModel() : inputPath_(""), scoreModel_(0), featureSet_(0), indLocalFeatureSet_(false) {}
	~MlModel();

	const string& getInputPath() const { return inputPath_; }
	void  setInputPath(string& s)      { inputPath_ = s; }

	// performs all steps involved in model creation, including the gneration of operators and the model training and writing
	void  createModel(MlTrainingContainer& mtc);

	const MlOperatorList& getOperatorList() const { return operatorList_; }
	const MlFeatureSet*	  getFeatureSet() const   { return featureSet_; }
	const MlScoreModel*   getScoreModel() const   { return scoreModel_; }

	bool  readModel(const char* path, const char* resourcePathDir = NULL);
	bool  readModel(ifstream& fis, const char* resourcePathDir = NULL);

	bool  writeModel(const char* path) const;
	bool  writeModel(ostream& fos) const;

	void outputScoreModelAnalysisToLog(MlTrainingContainer* params, const char* logFilePath) const;

private:
	
	string inputPath_; // the path for the model files

	MlOperatorList operatorList_;
	MlFeatureSet*  featureSet_;		/// always a pointer to a hashed version
	MlScoreModel*  scoreModel_;

	bool		   indLocalFeatureSet_; /// indicates if the MlFeatureSet was allocated for this object using new

	static map<string, MlFeatureSet*> hashFeatureSets_; /*! each feature set that is read is hashed so if the
														   same file needs to be read again, it doesn't need to be opened. */

	/*! \brief checks if the feature set is in the hash, otherwise reads from the file.
	*/
	void readFeatureSet(const string& path);

	// adds operators and adjusts feature names
	// if updateDataForNewFeatures==true then the sample vectors are reflected to 
	void addOperators(const MlOperatorList& operatorsToAdd, 
					  MlTrainingContainer& mtc,
					  bool updateDataForNewFeatures);
};





#endif
