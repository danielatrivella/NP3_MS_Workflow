#include "mllogistic.h"

void  MlLogisticModel::trainModel(MlTrainingContainer* params)
{
	if (params->getTrainingAlgorithm() == SMT_LOGISTIC_CG)
	{
		learnCG(params);
		return;
	}

	error("Unrecognized model training type for logistic regression: ",params->getTrainingAlgorithm());
}

/**********************************************************************
Conguate gradient algorithm.
Follows the pseudo code of Daume 04 (megam)
Label 0 == class -1
Label 1 == class +1
***********************************************************************/
void MlLogisticModel::learnCG(MlTrainingContainer* params)
{
	params->initialize();

	const MlDataSet* trainingDataSet = params->getTrainingDataSet();
	const MlDataSet* testDataSet     = params->getTestDataSet();
	const vector<size_t>& trainingIdxs = params->getTrainingIdxs();
	const vector<size_t>& testIdxs	   = params->getTestIdxs();

	const bool performTest = (testDataSet && testIdxs.size()>0);

	if (trainingDataSet->getNumClasess() != 2)
		error("learnCG accepts only datasets with 2 classes, yor data has ",
			trainingDataSet->getNumClasess());

	const double  lambda  = params->getLambda();
	const double  perplexityDelta  = params->getPerplexityDelta();
	const size_t numFeatures = trainingDataSet->getNumBasicFeatures(); // F
	const size_t numTraining = trainingIdxs.size();    // N
	const size_t numTest	 = testIdxs.size();
	size_t  reportFrequency = params->getVerboseLevel();

	// data structures used for training
	vector<value_t>& w = weights_;
	vector<value_t>  wtx(numTraining, 0.0);
	vector<value_t>  u(numFeatures,      0.0);
	
	vector<double>  g(numFeatures,   0.0); 
	vector<double>  gPrime(numFeatures, 0.0);

	vector<value_t>	 sigmaWtx(numTraining,0);
	vector<value_t>	 sigmaMinusWtx(numTraining,0);

	vector<float> trainingProbs(numTraining);
	vector<float> testProbs(numTest);
	vector<value_t> bestW(numFeatures);

	// initialize weights
	if (params->getInputPath().length() > 1)
	{
		const string modelFile = params->getInputPath() + "_scr.txt";
		if (readModel(modelFile.c_str()))
			params->setIndClearWeights(false);
	}

	if (params->getIndClearWeights())
		weights_.clear();

	weights_.resize(numFeatures,     0.0);

	const vector<MlSample>& trainingSamples = trainingDataSet->getSamples();
	double previousPerplexity = MAX_FLOAT;
	float  bestTestError=1.0;
	size_t bestTestRound=0;
	float  bestTrainingError=1.0;
	size_t bestTrainingRound=0;

	MlTrainingTerminationDecider terminateDecider;
	size_t round;
	for (round=0; round<params->getMaxNumIterations(); round++)
	{
		bool terminateTraining = false;

		// compute sigma(wtx)
		for (size_t i=0; i<numTraining; i++)
			sigmaWtx[i] = 1.0 / (1.0 + exp(-wtx[i]));
		
		// set g' <- -lamda*w
		if (lambda > 0.0)
		{
			const double minusLambda= -lambda;
			for (size_t i=0; i<numFeatures; i++)
				gPrime[i] = w[i]*minusLambda;
		}
		else
			memset(static_cast<void*>(&gPrime[0]),0,numFeatures*sizeof(double));

		// for i=1 to N
		//	 g' <- g' + sigma(-yn xtw)ynxn
		// end for
		for (size_t i=0; i<numTraining; i++)
		{
			const size_t idx=trainingIdxs[i];
			const MlSample& sample = trainingSamples[idx];
			const value_t sigTimesY = sample.weight * 
				(sample.label == 0 ? -sigmaWtx[i] : 1.0 - sigmaWtx[i]);
			const vector<IdxVal>& pairs = sample.pairs;
			for (size_t j=0; j<pairs.size(); j++)
				gPrime[pairs[j].index] += static_cast<double>(sigTimesY * pairs[j].value);
		}
		
		if (round > 0)
		{
			// compute beta
			double numerator=0, denominator=0;
			for (size_t i=0; i<numFeatures; i++)
			{
				const double gPrimeMinusG = gPrime[i]-g[i];
				numerator  += gPrime[i]*gPrimeMinusG;
				denominator+= u[i]*gPrimeMinusG;
			}
			if (denominator==0.0)
			{
				params->setIndHadInternalError(true);
				terminateTraining = true;
				break;
			}

			const	value_t beta = static_cast<value_t>(numerator/denominator);
		
			// update u
			for (size_t i=0; i<numFeatures; i++)
				u[i] = static_cast<value_t>(gPrime[i]) - beta*u[i];
		}
		else // u=g' for the first round
			for (size_t i=0; i<numFeatures; i++)
				u[i]=static_cast<value_t>(gPrime[i]);
	
		// compute z
		double Zdenominator=0.0;
		if (lambda > 0.0)
			Zdenominator = lambda * computeDotProduct(u,u);

		for (size_t i=0; i<numTraining; i++)
		{
			const MlSample& sample = trainingSamples[trainingIdxs[i]];
			const value_t utx=computeDotProduct(u,sample.pairs);
			Zdenominator += static_cast<double>(sample.weight*sigmaWtx[i]*(1.0-sigmaWtx[i])*utx*utx);
		}
		
		if (Zdenominator<=0.0)
		{
			params->setIndHadInternalError(true);
			terminateTraining = true;
			break;
		}

		const value_t z = computeDotProduct(gPrime,u)/static_cast<value_t>(Zdenominator);

		// update w <- w + zu
		value_t deltaW=0;
		for (size_t i=0; i<numFeatures; i++)
		{
			const value_t zu = z*u[i];
			w[i] += zu;
			deltaW += (zu<0.0 ? -zu : zu);
		}

		// update wtx
		for (size_t i=0; i<numTraining; i++)
			wtx[i] = computeDotProduct(trainingSamples[trainingIdxs[i]].pairs, w);

		// g <- g'
		memcpy(&g[0],&gPrime[0],numFeatures*sizeof(double));


		// compute errors
		{
			double trainingLogLikelihood=0.0, testLogLikelihood=0.0;
			const double trainingError = calcErrorRateWithLogLikelihood(trainingDataSet, trainingIdxs,
																	false, &trainingLogLikelihood);
			double testError=1.0;
			if (performTest)
				testError = calcErrorRateWithLogLikelihood(testDataSet, testIdxs, false, &testLogLikelihood);

			if (round>10)
				reportFrequency=2;
			if (round>25)
				reportFrequency=5;
			if (round>50)
				reportFrequency=10;
		
			if (reportFrequency>0 && round % reportFrequency == 0)
			{
				cout << round << "\t" << setprecision(5) << trainingLogLikelihood << "\t" << fixed << setprecision(5) << trainingError;
				if (performTest)
					cout <<"\t\t" << testLogLikelihood << "\t" << fixed << setprecision(5)<< testError;
				cout << endl;
			}
	
			
			if (performTest)
			{
				if (testError<=bestTestError)
				{
					bestTestRound=round;
					bestTestError=testError;
					memcpy(&bestW[0],&w[0],numFeatures*sizeof(value_t)); // copy weights
				}
			}
			
			if (trainingError<=bestTrainingError)
			{
				bestTrainingRound=round;
				bestTrainingError=trainingError;
				if (! performTest)
					memcpy(&bestW[0],&w[0],numFeatures*sizeof(value_t)); // copy weights
			}
			terminateTraining = terminateDecider.shouldTerminateTraining(round, 
									trainingLogLikelihood, testLogLikelihood);
		}

		if (terminateTraining)
			break;
	}

	if (! params->getIndHadInternalError())
	{
		params->setIndNormalTermination(true);
	}
	else
		cout << "Warning: encountered mathemtical error while training!" << endl;

	weights_ = bestW;

	cout << "Terminated after " << round << " rounds." << endl;
	cout << "Best training error  " << fixed << setprecision(8) << bestTrainingError << " (round " << bestTrainingRound << ")" << endl;
	if (performTest)
	cout << "Best test error      "  << bestTestError     << " (round " << bestTestRound << ")" << endl;

	indWasInitialized_ = true;

/*	if (params->getOutputPath().length()>1)
	{
		const string modelFile = params->getOutputPath() + "_src.txt";
		MlScoreModel::writeModel(modelFile.c_str());
	}
	else
		writeModel();*/
}




double	MlLogisticModel::calcScoreForClass(const MlSample* sample, int label) const
{
	assert( label == 0 || label == 1);

	double wtx=0.0;
	for (size_t i=0; i<sample->pairs.size(); i++)
		wtx += static_cast<double>(weights_[sample->pairs[i].index] * sample->pairs[i].value);

	const double e = exp(wtx);
	if (label == 1)
		return (e/(1.0+e));

	return (1.0/(1.0+e));
}

void MlLogisticModel::calcScores(const MlSample* sample, vector<float>& scores) const
{
	double wtx=0.0;
	for (size_t i=0; i<sample->pairs.size(); i++)
		wtx += static_cast<double>(weights_[sample->pairs[i].index] * sample->pairs[i].value);
	scores.resize(2);
	scores[0]=1.0/(1.0+exp(wtx));
	scores[1]=1.0-scores[0];
}





// assumes probaiblities were calculated for class 1
double	MlLogisticModel::calcErrorRateWithLogLikelihood(const   MlDataSet* mld, 
													 const   vector<size_t>& idxs,
													 bool    verbose,
													 double* logLikelihood) const
{
	assert( mld->getNumSamples()>= idxs.size() );

	if (idxs.size() ==0)
		return 0.0;

	double weightCorrect=0.0;
	double weightTotal=0.0;
	if (logLikelihood)
		*logLikelihood = 0.0;

	for (size_t i=0; i<idxs.size(); i++)
	{
		const MlSample& sample = mld->getSample(idxs[i]);
		const double prob0 = calcScoreForClass(&sample, 0);
		const double p = (sample.label == 0 ? prob0 : 1.0 - prob0);
		
		weightTotal += sample.weight;
		if (p>0.5)
			weightCorrect+=sample.weight;
		
		if (logLikelihood)
			*logLikelihood += (sample.weight * log(p));

		if (verbose)
			cout << idxs[i] << "\t"  << prob0 << "\t" << 1.0-prob0 << endl;
	}

	if (logLikelihood)
		*logLikelihood /= weightTotal;
	
	return (1.0 - (weightCorrect/weightTotal));
}

double	MlLogisticModel::calcErrorRate(const MlDataSet* mld, bool verbose) const
{
	const vector<MlSample>& samples = mld->getSamples();
	double weightCorrect=0.0;
	double weightTotal=0.0;
	for (size_t i=0; i<samples.size(); i++)
	{
		const MlSample& sample = samples[i];
		const float prob0 = calcScoreForClass(&sample, 0);
		const float p = (sample.label == 1 ? prob0 : 1.0-prob0);
		
		weightTotal+= sample.weight;
		if (p>0.5)
			weightCorrect+=sample.weight;

		if (verbose)
			cout << i << "\t"  << prob0 << "\t" << 1.0-prob0 << endl;
	}
	return (1.0 - (weightCorrect/weightTotal));
}


bool MlLogisticModel::readModel(const char* path)
{
	ifstream ifs(path);
	if (! ifs.good())
		return false;

	bool retVal = readModel(ifs);
	ifs.close();
	return retVal;
}




bool MlLogisticModel::readModel(ifstream& ifs)
{
	char buffer[256];
	while (ifs.good())
	{
		ifs.getline(buffer,256);
		if (ifs.gcount()>0 && buffer[0] != '#')
			break;
	}
		
	unsigned int numWeights=0;
	if (sscanf(buffer,"LOGISTIC_CG %u",&numWeights) != 1 &&
		sscanf(buffer,"LOGISTIC_BFGS %u",&numWeights) != 1)
	{
		cout << "Bad line in model file:" << endl << buffer << endl;
		return false;
	}

	weights_.resize(numWeights,0);
	while (ifs.good())
	{
		ifs.getline(buffer,256);
		
		if (! strncmp(buffer,"END_LOGISTIC",10))
			break;

		if (ifs.gcount() == 0 || buffer[0] != 'F')
			continue;

		size_t index;
		float  weight;
		istringstream iss(buffer+1);
		iss >> index >> weight;
		if (iss.fail())
		{
			if (strlen(buffer)<3)
				continue;

			cout << "Bad line in model file:" << endl << buffer << endl;
			return false;
		}
		if (index>weights_.size())
			error("Bad feature index in line: ",buffer);

		weights_[index]=weight;
	}
	return true;
}



bool	MlLogisticModel::writeModel(ostream& os) const
{
	os << SCORE_MODEL_NAMES[trainingType_] << "\t" << weights_.size() << endl;
	for (size_t i=0; i<weights_.size(); i++)
	{
		if (weights_[i] == 0.0)
			continue;

		char buffer[20];
		sprintf(buffer,"%e",weights_[i]);
		os << "F" << i << "\t" << buffer << endl;
	}
	os << "END_LOGISTIC" << endl;
	return true;
}


/*********************************************************************************
Ouputs a table that can be used to weigh the contribution of each of the fetures
**********************************************************************************/
void MlLogisticModel::outputModelAnalysis(const MlDataSet* mld, 
										  const vector<size_t>& idxs,
										  const MlFeatureSet* featureSet,
										  ostream& os) const
{
//	assert(weights_.size() == featureSet->getNumBasicFeatures() &&
//		   weights_.size() == mld->getNumBasicFeatures());

	mld->outputFeatureReports(idxs, featureSet, os);

	const size_t numFeatures = mld->getNumBasicFeatures();

	vector< vector<double> > totalContributionOfFeatures(2, vector<double>(numFeatures,0.0));
	vector< vector<double> > totalWeightOfFeatures(2, vector<double>(numFeatures,0.0));
	vector<double> totalWeightOfSamples(2,0);
	vector<int>	   numSamplesRead(2,0);

	for (size_t i=0; i<idxs.size(); i++)
	{
		const MlSample& sam = mld->getSample(idxs[i]);
		
		totalWeightOfSamples[sam.label]+=sam.weight;
		numSamplesRead[sam.label]++;
		for (size_t j=0; j<sam.pairs.size(); j++)
		{
			const IdxVal pair = sam.pairs[j];
			totalWeightOfFeatures[sam.label][pair.index] += sam.weight;
			totalContributionOfFeatures[sam.label][pair.index] += (weights_[j]*pair.value);
		}
	}

	// report
	os << setprecision(4) << fixed;
	os << endl << "Dataset contains:" << endl;
	os <<         "-----------------" << endl;
	os << "Class 0, " << numSamplesRead[0] << " samples with weight " << totalWeightOfSamples[0] << endl;
	os << "Class 1, " << numSamplesRead[1] << " samples with weight " << totalWeightOfSamples[1] << endl;
	os << endl << setprecision(3) << fixed;
	os << "Idx\tAlpha\t\tW0\tCon0\tRelCon0\t\tW1\tCon1\tRelCon1\t   FeatureName" << endl;
	os << "---\t-----\t\t--\t----\t-------\t\t--\t----\t-------\t   -----------" << endl;
	for (size_t i=0; i<numFeatures && i<weights_.size(); i++)
	{
		if (weights_[i] == 0.0)
			continue;

		os << i << "\t" << scientific << setprecision(4) << weights_[i] << "\t";
		os << setprecision(4) << fixed;
		os << totalWeightOfFeatures[0][i]/totalWeightOfSamples[0] << "\t"
		   			  << totalContributionOfFeatures[0][i]/totalWeightOfSamples[0] << "\t"
					  << totalContributionOfFeatures[0][i]/totalWeightOfFeatures[0][i] << "\t\t";
		os << totalWeightOfFeatures[1][i]/totalWeightOfSamples[1] << "\t"
		   << totalContributionOfFeatures[1][i]/totalWeightOfSamples[1] << "\t"
		   << totalContributionOfFeatures[1][i]/totalWeightOfFeatures[1][i] << "\t  ";
		os << featureSet->getFeature(i).getName() << endl;
	}


}
