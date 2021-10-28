#include "mlscoremodel.h"
#include "mllogistic.h"
#include "mlmaximumentropy.h"
#include "MlRankBoost.h"


size_t getModelTypeFormLabel(const char* s)
{
	for (size_t i=0; i<SMT_NUM_SCORE_MODEL_TYPES; i++)
		if (! strncmp(s,SCORE_MODEL_NAMES[i],strlen(SCORE_MODEL_NAMES[i])))
			return i;

	return MAX_UINT; // no model was found
}

void   printModelTypes(ostream& os)
{
	for (size_t i=0; i<SMT_NUM_SCORE_MODEL_TYPES; i++)
		os << i << "\t" << SCORE_MODEL_NAMES[i] << endl;
}




MlScoreModel* MlScoreModel::createDerivedModelObject(size_t typeIndex)
{
	MlScoreModel* newModel=NULL;
	if (typeIndex == SMT_LOGISTIC_CG || typeIndex == SMT_LOGISTIC_LMBFGS)
	{
		newModel = new MlLogisticModel;

	}
	else if (typeIndex == SMT_MAXIMUM_ENTROPY_LMBFGS || 
			 typeIndex == SMT_MAXIMUM_ENTROPY_CG_FR ||
			 typeIndex == SMT_MAXIMUM_ENTROPY_CG_PR)
	{
		newModel = new MlMaximumEntropyModel;
	}
	else if (typeIndex == SMT_RANKBOOST_ADA ||
			 typeIndex == SMT_RANKBOOST_LOGIT ||
			 typeIndex == SMT_RANKBOOST_SMOOTH)
	{
		newModel = new MlRankBoostModel;
	}

	
	if (newModel)
	{
		newModel->trainingType_ = typeIndex;
		return newModel;
	}
	return 0;
}


size_t MlScoreModel::getIndModelIsBinary(size_t typeIndex)
{
	if (typeIndex == SMT_LOGISTIC_CG || typeIndex == SMT_LOGISTIC_LMBFGS)
		return true;

	return false;
}

// reads the model file and returns a pointer of
// the appropriate derived class type
MlScoreModel* MlScoreModel::readAndAllocateModel(const char* path)
{
	ifstream ifs(path);
	if (! ifs.good())
		return 0;

	MlScoreModel* retPointer = readAndAllocateModel(ifs);
	ifs.close();
	return retPointer;
}


// reads the model file and returns a pointer of
// the appropriate derived class type
// does this by peeking at the line in the model that describes the type,
// creating the approapriate model object, and calling its read function
MlScoreModel* MlScoreModel::readAndAllocateModel(ifstream& ifs)
{
	char buffer[256];

	streampos pos=0;
	while (ifs.good())
	{
		pos = ifs.tellg();
		ifs.getline(buffer,256);
		if (ifs.gcount()>0 && buffer[0] != '#')
			break;
	}
	ifs.seekg(pos); // go back

	char modelTypeStr[64];
	if (sscanf(buffer,"%s",modelTypeStr) != 1)
		return 0;

	const size_t typeIndex = getModelTypeFormLabel(modelTypeStr);
	assert(typeIndex<MAX_SIZE_T);

	MlScoreModel* scoreModel = createDerivedModelObject(typeIndex);
	if (! scoreModel)
		return 0;

	scoreModel->readModel(ifs);
	return scoreModel;
}


bool MlScoreModel::writeModel(const char* path) const
{
	ofstream ofs(path);
	if (! ofs.is_open())
		error("Couldn't open model file for writing: ",path);

	if (! writeModel(ofs))
		error("Could not write model file correctly: ",path);

	ofs.close();
	return true;
}




bool MlTrainingTerminationDecider::shouldTerminateTraining(size_t round, double trainLL, double testLL)
{
	if (testLL == MIN_FLOAT) // use log likelihood
	{
		if (trainLL - bestTrainLL_ < logLikelihoodTolerance_)
		{
			cout << "Training log likelihhod diff is : " << scientific << trainLL-bestTrainLL_ << endl;
			cout << "Training should be stopped..." << endl;
			return true;
		}
		bestTrainLL_ = trainLL;
	}
	else // use test error
	{
		size_t numRoundsWithNoImprovment =10;
		if (round>100)
			numRoundsWithNoImprovment = 10 + static_cast<size_t>(0.05*(round-100));
		if (round>200)
			numRoundsWithNoImprovment = 15 + static_cast<size_t>(0.025*(round-200));
		if (numRoundsWithNoImprovment>35)
			numRoundsWithNoImprovment = 35;

		if (round - bestTestRound_ > numRoundsWithNoImprovment && testLL < bestTestLL_)
		{
			cout << "Testing likelihhod hasn't improved in " << numRoundsWithNoImprovment << " rounds..." << endl;
			cout << "Training should be stopped..." << endl;
			return true;
		}

		if (testLL > bestTestLL_)
		{
			bestTestLL_ = testLL;
			bestTestRound_ = round;
		}
	}
	return false;
}




