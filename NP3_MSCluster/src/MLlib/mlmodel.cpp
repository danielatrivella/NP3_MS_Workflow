#include "mlmodel.h"

/*! This is a hash for the feature sets. Since many models will use the same feature files,
	it is more ecnomical and consistent to read and maintain only one version of the file.*/
map<string, MlFeatureSet*> MlModel::hashFeatureSets_;

MlModel::~MlModel()
{
	if (indLocalFeatureSet_)
		delete featureSet_;
}

void MlModel::readFeatureSet(const string& path)
{
	map<string, MlFeatureSet*>::const_iterator it = hashFeatureSets_.find(path);
	if (it != hashFeatureSets_.end())
	{
		featureSet_ = it->second;
	}
	else
	{
		MlFeatureSet* newFeatureSet = new MlFeatureSet;
		newFeatureSet->readFeatureFile(path.c_str());
		hashFeatureSets_[path] = newFeatureSet;
		featureSet_ = newFeatureSet;
		indLocalFeatureSet_ = true;
	}
}

// first tries the the "_mod.txt" file, otherwise tries the
// the three files: "_opr.txt", "_fet.txt", and "_scr.txt"
bool  MlModel::readModel(const char* path, const char* resourcePathDir)
{
	inputPath_ = path;

	const string modelFile = inputPath_ + "_mod.txt";
	ifstream ifs(modelFile.c_str());

	if (ifs.is_open())
	{
		bool retVal = readModel(ifs, resourcePathDir);
		ifs.close();
		return retVal;
	}
	else // try seperate files
	{
		const string featureFile = inputPath_ + "_fet.txt";
		if (checkIfFileExists(featureFile.c_str()))
		{
			readFeatureSet(featureFile);
		}
			
		const string operatorFile = inputPath_ + "_opr.txt";
		operatorList_.readOperatorList(operatorFile.c_str());

		const string scoreFile = inputPath_ + "_scr.txt";
		scoreModel_ = MlScoreModel::readAndAllocateModel(scoreFile.c_str());
	}

	return (featureSet_ && featureSet_->getIndWasInitialized() &&
		    operatorList_.getIndWasInitialized() &&
			scoreModel_->getIndWasInitialized() );
}

// assumes that the relevant files are in a successive order in the stream
// before each type of file should be preceded by a line with either "OPERATORS",
// "FEATURES", or "SCOREMODEL".
bool MlModel::readModel(ifstream& ifs, const char* resourcePathDir)
{
	char buffer[256];
	while (ifs.good() && ifs.getline(buffer,256) && ifs.gcount()>0 && buffer[0] =='#')
		;
	
	if (! strncmp(buffer,"FEATURES",8))
	{
		istringstream iss(buffer+8);
		string path;
		iss >> path;

		if (path.length()>2)
		{
			if (resourcePathDir)
				path = std::string(resourcePathDir) + "/" + path;
			readFeatureSet(path);
		}
		else
		{
			featureSet_ = new MlFeatureSet;
			indLocalFeatureSet_ = true;
			featureSet_->readFeatureFile(ifs);
		}
		while (ifs.good() && ifs.getline(buffer,256) && ifs.gcount()==0 && buffer[0] =='#')
			;
	}

	if (! strcmp(buffer,"OPERATORS"))
	{
		operatorList_.readOperatorList(ifs);
		while (ifs.good() && ifs.getline(buffer,256) && ifs.gcount()>0 && buffer[0] =='#')
			;
	}

	if (! strcmp(buffer,"SCOREMODEL"))
		scoreModel_ = MlScoreModel::readAndAllocateModel(ifs);

	return (featureSet_ &&
			featureSet_->getIndWasInitialized() &&
		    operatorList_.getIndWasInitialized() &&
			scoreModel_->getIndWasInitialized() );
}


bool  MlModel::writeModel(const char* path) const
{
	const string modelFile = static_cast<string>(path) + "_mod.txt";
	ofstream ofs(modelFile.c_str());
	
	bool retVal=writeModel(ofs);
	ofs.close();
	return retVal;
}


bool  MlModel::writeModel(ostream& os) const
{
	bool retVal=true;

	if (featureSet_->getIndWasInitialized())
	{
		if (featureSet_->path_.length()>0)
		{
			os << "FEATURES " << featureSet_->path_ << endl;
		}
		else
		{
			os << "FEATURES" << endl;
			featureSet_->writeFeatureFile(os);
		}
	}
	else
		retVal=false;

	if (operatorList_.getIndWasInitialized())
	{
		os << "OPERATORS" << endl;
		if (! operatorList_.writeOperatorList(os))
			retVal = false;
	}
	else
		retVal=false;

	if (scoreModel_->getIndWasInitialized())
	{
		os << "SCOREMODEL" << endl;
		if (! scoreModel_->writeModel(os))
			retVal = false;
	}
	else
		retVal=false;

	return retVal;
}

// performs all steps involved in model creation, including the gneration of operators and the model training and writing
void  MlModel::createModel(MlTrainingContainer& mtc)
{

	mtc.initialize();

	inputPath_ = mtc.getInputPath();

	if (mtc.getInputName() != "Model" && mtc.getIndReadPreviousModel())
	{
		readModel(inputPath_.c_str());
	}
	else
	{
		if (! featureSet_)
		{
				featureSet_ = new MlFeatureSet;
				indLocalFeatureSet_ = true;
		}
		
		if (! featureSet_->getIndWasInitialized())
			featureSet_->initializeWithDefaultFeatures(mtc.getTrainingDataSet()->getNumBasicFeatures());
	}

	bool addIndicators = (mtc.getFeatureGenerationType() == FGT_REGRESSION);

/*	mtc.getAuxilarySearchData()->printDetailedSummary(mtc.getTrainingDataSet(),featureSet_);
	exit(0);

	MlOperatorList initialSplitSet;
	initialSplitSet.selectInitialSplitOperators(mtc, featureSet_, addIndicators, 1);
	addOperators(initialSplitSet, mtc, true);

	//operatorList_.selectInitialUnaryOperators(mtc, featureSet_, true);
	//operatorList_.selectInitialSplitOperators();
	//operatorList_.createNewOperators();

*/

	if (! scoreModel_)
		scoreModel_ = MlScoreModel::createDerivedModelObject(mtc.getTrainingAlgorithm());

	scoreModel_->trainModel(&mtc);

	writeModel(mtc.getOutputPath().c_str());
}





// adds operators and adjusts feature names
void MlModel::addOperators(const MlOperatorList& operatorsToAdd, 
						   MlTrainingContainer& mtc,
						   bool updateDataForNewFeatures)
{
	MlOperatorSearchData* auxilaryData = mtc.getAuxilarySearchData();

	assert( auxilaryData != 0 );

	for (size_t opIdx=0; opIdx<operatorsToAdd.executionOrder_.size(); opIdx++)
	{
		const size_t newOpType = operatorsToAdd.executionOrder_[opIdx].idx1;
		const size_t newOpPos  = operatorsToAdd.executionOrder_[opIdx].idx2;
		
		switch (newOpType)
		{
			case OT_DROP:
				operatorList_.drops_.push_back(operatorsToAdd.drops_[newOpPos]);
				operatorList_.executionOrder_.push_back(IdxPair(newOpType,operatorList_.drops_.size()));
				break;

			case OT_INDICATOR:
			{
				IndicatorOperator newOp = operatorsToAdd.indicators_[newOpPos];
				MlFeature& oldFeature = featureSet_->basicFeatures_[newOp.sourceIdx];
				if (oldFeature.getIndIsBool()) // no need to add, it is already a boolean value
					break;
				
				if (newOp.sourceIdx != newOp.targetIdx)
				{
					newOp.targetIdx=featureSet_->basicFeatures_.size();
					MlFeature newFeature = oldFeature;
					newFeature.setIndIsBool(true);
				//	newFeature.setIndConsiderConditional(false);
				//	newFeature.setIndConsiderUnary(false);
					newFeature.addSuffixToFeatureName("IND");
					featureSet_->basicFeatures_.push_back(newFeature);
				}
				else
				{
					oldFeature.addSuffixToFeatureName("IND");
					oldFeature.setIndIsBool(true);
				}
				operatorList_.executionOrder_.push_back(IdxPair(newOpType,operatorList_.indicators_.size()));
				operatorList_.indicators_.push_back(newOp);
				break;
			}

			case OT_NORMALIZATION:
			{
				NormalizationOperator newOp = operatorsToAdd.normalizations_[newOpPos];
				MlFeature& oldFeature = featureSet_->basicFeatures_[newOp.sourceIdx];

				if (newOp.sourceIdx != newOp.targetIdx)
				{
					newOp.targetIdx=featureSet_->basicFeatures_.size();
					MlFeature newFeature = oldFeature;
					newFeature.addSuffixToFeatureName("NORM");
					featureSet_->basicFeatures_.push_back(newFeature);
				}
				else
					oldFeature.addSuffixToFeatureName("NORM");

				operatorList_.executionOrder_.push_back(IdxPair(newOpType,operatorList_.normalizations_.size()));
				operatorList_.normalizations_.push_back(newOp);
				break;
			}

			case OT_FUNCTION:
			{
				FunctionOperator newOp = operatorsToAdd.functions_[newOpPos];
				MlFeature& oldFeature = featureSet_->basicFeatures_[newOp.sourceIdx];
				if (newOp.sourceIdx != newOp.targetIdx)
				{
					newOp.targetIdx=featureSet_->basicFeatures_.size();
					MlFeature newFeature = oldFeature;
					newFeature.addSuffixToFeatureName(unaryOperatorLabels[newOp.type]);
					featureSet_->basicFeatures_.push_back(newFeature);
				}
				else
					oldFeature.addSuffixToFeatureName(unaryOperatorLabels[newOp.type]);

				operatorList_.executionOrder_.push_back(IdxPair(newOpType,operatorList_.functions_.size()));
				operatorList_.functions_.push_back(newOp);
				break;
			}

			case OT_SPLIT:
			{
				SplitOperator newOp = operatorsToAdd.splits_[newOpPos];
				const MlFeature oldFeature = featureSet_->basicFeatures_[newOp.sourceIdx];

				for (size_t i=0; i<newOp.indexesForBinValues.size(); i++)
					if (newOp.indexesForBinValues[i]<MAX_UINT)
					{
						newOp.indexesForBinValues[i]=featureSet_->basicFeatures_.size();
						MlFeature newFeature = oldFeature;
						ostringstream oss;
						oss << "F"<< newOp.sourceIdx << "_SPLIT_VAL_" << i+1 << "/" << newOp.indexesForBinValues.size();
						string name = oss.str();
						newFeature.setName(name);
						featureSet_->basicFeatures_.push_back(newFeature);
					}
				for (size_t i=0; i<newOp.indexesForBinIndicators.size(); i++)
					if (newOp.indexesForBinIndicators[i]<MAX_UINT)
					{
						newOp.indexesForBinIndicators[i]=featureSet_->basicFeatures_.size();
						MlFeature newFeature = oldFeature;
						ostringstream oss;
						oss << "F" << newOp.sourceIdx << "_SPLIT_IND_" << i+1 << "/" << newOp.indexesForBinIndicators.size();
						string name = oss.str();
						newFeature.setName(name);
				//		newFeature.setIndConsiderConditional(false);
				//		newFeature.setIndConsiderUnary(false);
						newFeature.setIndIsBool(true);
						featureSet_->basicFeatures_.push_back(newFeature);
					}
				operatorList_.executionOrder_.push_back(IdxPair(newOpType,operatorList_.splits_.size()));
				operatorList_.splits_.push_back(newOp);
				break;
			}

			case OT_CONDITIONAL:
			{
				ConditionalOperator newOp = operatorsToAdd.conditionals_[newOpPos];
				MlFeature newFeature;
	
				// the groups of the new feature are the union of the sources' features
				for (size_t i=0; i<newOp.sourceIdxs.size(); i++)
					newFeature.addGroups(featureSet_->basicFeatures_[newOp.sourceIdxs[i]].getGroups());

				newOp.targetIdx = featureSet_->basicFeatures_.size();

				// create name
				ostringstream oss;
				oss << conditionalOperatorLabels[newOp.conditionType];
				oss << "_F:" << makeRangeString(newOp.sourceIdxs);
				string conditionName = oss.str();
				if (newOp.resultType != CVT_FEATURE)
				{
					oss << "_{" << conditionalValueLabels[newOp.resultType] << "}";	
				}
				else
					oss << "_{F" << newOp.sourceIdxs[0] << "}";

				string featureName = oss.str();
				newFeature.setName(featureName);
				if (newOp.resultType==CVT_BOOL)
					newFeature.setIndIsBool(true);

				featureSet_->basicFeatures_.push_back(newFeature);
				
				// add bool indicator feature
				if (newOp.resultType != CVT_BOOL && newOp.indexForBool < MAX_UINT)
				{
					newOp.indexForBool = featureSet_->basicFeatures_.size();
					MlFeature boolFeature = newFeature;
					string boolName = conditionName + "_{IND}";
					boolFeature.setName(boolName);
					boolFeature.setIndIsBool(true);
				//	boolFeature.setIndConsiderConditional(false);
				//	boolFeature.setIndConsiderUnary(false);
					featureSet_->basicFeatures_.push_back(boolFeature);
				}

				operatorList_.executionOrder_.push_back(IdxPair(newOpType,operatorList_.conditionals_.size()));
				operatorList_.conditionals_.push_back(newOp);
				break;
			}

			default:
				error("Unrecognized operator type: ",newOpType);

		};
	}

	// update the feature space size
	operatorList_.operatorFeatureSpaceSize_ = featureSet_->basicFeatures_.size();
	auxilaryData->numFeatures_ = featureSet_->basicFeatures_.size();

	if (updateDataForNewFeatures)
	{
		if (mtc.getTrainingIdxs().size()>0)
		{
			assert(mtc.getTrainingDataSet() != 0);
			operatorList_.applyOperatorList(mtc.getTrainingDataSet(), mtc.getTrainingIdxs());
			mtc.getTrainingDataSet()->setMaxFeatureIndedx(featureSet_->basicFeatures_.size()-1);
		}

		if (mtc.getTestIdxs().size()>0)
		{
			assert(mtc.getTestDataSet() != 0);
			operatorList_.applyOperatorList(mtc.getTestDataSet(), mtc.getTestIdxs());
			mtc.getTestDataSet()->setMaxFeatureIndedx(featureSet_->basicFeatures_.size()-1);
		}

		auxilaryData->readInfoFromDataSet(mtc.getTrainingDataSet(), mtc.getTrainingIdxs() );
	}
}

void MlModel::outputScoreModelAnalysisToLog(MlTrainingContainer* params, const char* logFile) const
{
	ofstream fs(logFile);

	if (! fs.is_open())
		error("Couldn't open log file for writing: ",logFile);
	scoreModel_->outputModelAnalysis(params->getTrainingDataSet(), params->getTrainingIdxs(),
				 featureSet_, fs);
	fs.close();
}

