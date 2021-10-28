#include "mltrainingcontainer.h"
#include "mloperatorlist.h"
#include "mlscoremodel.h"

string MlTrainingContainer::getInputPath() const
{ 
	string inputPath = inputDir_;
	if (inputDir_ != "" && inputDir_[inputDir_.length()-1] != '/')
		inputPath += "/";
	
	inputPath += inputName_;
	return inputPath; 
}

string MlTrainingContainer::getOutputPath() const 
{ 
	string outputPath = outputDir_;
	if (outputDir_ != "" && outputDir_[outputDir_.length()-1] != '/')
		outputPath += "/";
	
	outputPath += outputName_;
	return outputPath; 
}

void MlTrainingContainer::printParameters(ostream& os) const
{
	
	os << "Training set paths: " << endl;
	for (size_t i=0; i<trainingSetPaths_.size(); i++)
		os << i+1 << " " << trainingSetPaths_[i] << endl;
	os << "Test set paths: " << endl;
	for (size_t i=0; i<testSetPaths_.size(); i++)
		os << i+1 << " " << testSetPaths_[i] << endl;
	if (inputDir_ != "")
		os << "Input dir:  " << inputDir_ << endl;
	if (inputDir_ != "" && inputName_ != "")
		os << "Input name: " << inputName_ << endl;
	if (outputDir_ != "")
		os << "Output dir: " << outputDir_ << endl;
	if (outputName_ != "")
		os << "Output dir: " << outputName_ << endl;
	os << "Maximal number of iterations : " << maxNumIterations_ << endl;
	os << "Vebose level					: " << verboseLevel_ << endl;
	os << "Regularization (Lambda)		: " << lambda_ << endl;
	if (testRatio_>0)
		os << "Ratio of samples used for testing: " << testRatio_ << endl;
	os << "Perplexity delta at which training is terminated: " << perplexityDelta_ << endl;
	os << "Number of training samples used: " << trainingIdxs_.size() << endl;
	os << "Number of test  samples used:    " << testIdxs_.size() << endl;
}


void MlTrainingContainer::performTrainingTestSplit(float ratio)
{
	if (! trainingDataSet_)
		error("Must read dataset before splitting!");

	testRatio_   = ratio;
	testDataSet_ = trainingDataSet_;

	trainingIdxs_.clear();
	testIdxs_.clear();

	const size_t n = trainingDataSet_->getNumSamples();

	chooseKFromN(static_cast<size_t>((1.0-testRatio_)*n+1),n,trainingIdxs_);
	if (trainingIdxs_.size()<2)
		error("Insufficient number of samples for training");

	if (trainingIdxs_.size()<n)
	{
		testIdxs_.reserve(n-trainingIdxs_.size());
		size_t s=0;
		while (s<trainingIdxs_[0])
			testIdxs_.push_back(s++);

		s++;
		size_t r=1;
		while (r<trainingIdxs_.size())
		{
			while (s<trainingIdxs_[r])
				testIdxs_.push_back(s++);
			s++;
			r++;
		}

		s=trainingIdxs_[r-1]+1;
		while (s<n)
			testIdxs_.push_back(s++);
	}
	if (testIdxs_.size() + trainingIdxs_.size() != n)
	{
		cout << "Training: " << trainingIdxs_.size() << endl;
		cout << "Test	 : " << testIdxs_.size() << endl;
	}
	assert(trainingIdxs_.size() + testIdxs_.size() == n);
}


bool  MlTrainingContainer::initialize(bool verbose)
{
	if (outputName_ == "" && inputName_ != "")
		outputName_ = inputName_;

	if (inputDir_ != "" && outputDir_ == "")
		outputDir_ = inputDir_;

	// set training type
	featureGenerationType_ =  MlOperatorList::determineFeatureGenerationType(trainingAlgorithm_);

	// read data files if needed
	if (! indDataWasInitialized_)
	{
		initalizeDataSets(verbose);

		// init search data
		if (numGenerationRounds_ < MAX_UINT)
		{
			auxilaryData_ = new MlOperatorSearchData;
			auxilaryData_->readInfoFromDataSet(trainingDataSet_, trainingIdxs_);
			auxilaryData_->determineNumberOfValues();
		}
	}


	return true;
}


// reads and processes the data (including split to training and test sets)
bool  MlTrainingContainer::initalizeDataSets(bool verbose)
{
	if (trainingSetPaths_.size()>0)
	{
		if (! trainingDataSet_)
		{
			trainingDataSet_ = new MlDataSet;
		}
		else
			trainingDataSet_->clear();

		for (size_t i=0; i<trainingSetPaths_.size(); i++)
			trainingDataSet_->readDataFile(trainingSetPaths_[i].c_str());

		if (MlScoreModel::getIndModelIsBinary(trainingAlgorithm_))
		{
			trainingDataSet_->convertClassLabelsForBinary(labelToConvertTo0_);
		}
		else
			trainingDataSet_->tallyClassStatistics();

		if (desiredClassWeights_.size()>0)
		{
			if (desiredClassWeights_.size() != trainingDataSet_->getNumClasess())
				error("Trying to re-weight with different number of classes!");
			trainingDataSet_->reweight(desiredClassWeights_);
			trainingDataSet_->printDatasetStatistics();
		}
		else if (verbose)
			trainingDataSet_->printDatasetStatistics();
	}
	else
		error("No training file was supplied!");

	if (testSetPaths_.size()>0)
	{
		if (! testDataSet_)
		{
			testDataSet_ = new MlDataSet;
		}
		else	
			testDataSet_->clear();

		for (size_t i=0; i<testSetPaths_.size(); i++)
			testDataSet_->readDataFile(testSetPaths_[i].c_str());
		
		if (MlScoreModel::getIndModelIsBinary(trainingAlgorithm_))
		{
			testDataSet_->convertClassLabelsForBinary(labelToConvertTo0_);
		}
		else
			testDataSet_->tallyClassStatistics();

		if (verbose)
			testDataSet_->printDatasetStatistics();

		// mark all training sample indexes as belonging to trianing
		// and all test sample indexes as belonging to test
		trainingIdxs_.resize(trainingDataSet_->getNumSamples());
		for (size_t i=0; i<trainingIdxs_.size(); i++)
			trainingIdxs_[i]=i;

		testIdxs_.resize(testDataSet_->getNumSamples());
		for (size_t i=0; i<testIdxs_.size(); i++)
			testIdxs_[i]=i;
	}
	else if (testRatio_>0)
	{
		// the split takes care of assigning indexes to the test and train
		performTrainingTestSplit(testRatio_);
	}
	else
	{
		// mark all samples sample indexes as belonging to trianing
		// there are no test indexes
		trainingIdxs_.resize(trainingDataSet_->getNumSamples());
		for (size_t i=0; i<trainingIdxs_.size(); i++)
			trainingIdxs_[i]=i;
	}



	// check that datasets are consistent
	if (trainingDataSet_ && testDataSet_ && (trainingDataSet_ != testDataSet_) )
	{
		if (trainingDataSet_->getNumClasess() != testDataSet_->getNumClasess())
			error("Training and test datasets must have same number of classes!");
	}

	indDataWasInitialized_ = true;

//	trainingDataSet_->randomlyReduce(0.1);
//	trainingDataSet_->writeWholeDataFile("ds1.txt");
//	exit(0);

	return true;
}


/******************************************************************************
Scans the training and test Idxs, removes any sample index that is in
the list of idxs to be romved. Assumes the training and test samples belong to
the same dataset.
*******************************************************************************/
void MlTrainingContainer::removeTrainingAndTestIdxs(const vector<size_t>& idxsToRemove)
{
	if (testDataSet_ && testDataSet_ != trainingDataSet_)
		error("removeTrainingAndTestIdxs: testDataSet_ and trainingDataSet_ must point to the same object!");

	const size_t numSamples = trainingDataSet_->getNumSamples();
	vector<bool> inds(numSamples,false);

	for (size_t i=0; i<idxsToRemove.size(); i++)
	{
		if (idxsToRemove[i]>=numSamples)
			error("bad sample idx in idxsToRemove!");
		inds[idxsToRemove[i]]=true;
	}

	const size_t numPasses = (testDataSet_ ? 2 : 1);
	for (size_t pass=0; pass<numPasses; pass++)
	{
		vector<size_t>& idxs = (pass == 0 ? trainingIdxs_ : testIdxs_);
		for (size_t i=0; i<idxs.size(); i++)
			if (inds[idxs[i]])
				idxs[i]=MAX_UINT;
		sort(idxs.begin(),idxs.end());

		while (idxs.size()>0 && idxs[idxs.size()-1] == MAX_UINT)
			idxs.pop_back();
	}
}
