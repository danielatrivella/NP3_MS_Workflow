#include "mlmaximumentropy.h"



void computeGradient(const MlDataSet* ds,
					 const vector<size_t>& idxs,
					 const vector< vector<double> >& w,
					 const vector< vector<double> >& wtx, //
					 double lambda,
					 vector< vector<double> >& g)
{
	const size_t numClasses = w.size();
	const size_t numFeatures = w[0].size();

	// g<-lamda*w
	if (lambda != 0.0)
	{
		for (size_t c=0; c<numClasses; c++)
			for (size_t i=0; i<numFeatures; i++)
				g[c][i]=-lambda*w[c][i];
	}
	else
		for (size_t c=0; c<numClasses; c++)
			memset(&g[c][0],0,sizeof(double)*numFeatures);

	for (size_t i=0; i<idxs.size(); i++)
	{
		const MlSample& sample = ds->getSample(idxs[i]);
		const size_t label = static_cast<size_t>(sample.label);

		double z= wtx[0][i];
		for (size_t c=1; c<numClasses; c++)
			z=addInLogSpace(z,wtx[c][i]);

		for (size_t c=0; c<numClasses; c++)
		{
			double t = (label == c ? 1.0 : 0.0);
			t -= exp(wtx[c][i]-z);
			for (size_t j=0; j<sample.pairs.size(); j++)
			{
				const size_t& featureIndex = sample.pairs[j].index;
				const value_t& featureValue = sample.pairs[j].value;
				g[c][featureIndex] += t*featureValue*sample.weight;
			}
		}
	}
}


double computePosterior(const MlDataSet* ds,
					  const vector<size_t>& idxs,
					  const vector< vector<double> >& wtx,
					  const vector< vector<double> >& qtx,
					  double wtw,
					  double qtw,
					  double qtq,
					  double lambda,
					  double eta)
{
	const size_t numClasses = wtx.size();
	
	double p=0.0;
	if (lambda != 0.0)
		p = -0.5 *lambda *( wtw + eta*eta*qtq + 2*eta*qtw);
	
	for (size_t i=0; i<idxs.size(); i++)
	{
		const MlSample& sample = ds->getSample(idxs[i]);
		const size_t label = static_cast<size_t>(sample.label);
		double s=MIN_FLOAT;
		for (size_t c=0; c<numClasses; c++)
		{
			const double chi = wtx[c][i] - eta*qtx[c][i];
			s=addInLogSpace(s,chi);
			if (c==label)
				p+=chi*sample.weight;
		}
		p-=s*sample.weight;
	}
	return p;
}

double computePosterior(const MlDataSet* ds,
					  const vector<size_t>& idxs,
					  const vector< vector<double> >& wtx,
					  double wtw,
					  double lambda)
{
	const size_t numClasses = wtx.size();
	
	double p=0.0;
	if (lambda != 0.0)
		p = -0.5 *lambda * wtw ;
	
	for (size_t i=0; i<idxs.size(); i++)
	{
		const MlSample& sample = ds->getSample(idxs[i]);
		const size_t label = static_cast<size_t>(sample.label);
		double s=MIN_FLOAT;
		for (size_t c=0; c<numClasses; c++)
		{
			s=addInLogSpace(s,wtx[c][i]);
			if (c==label)
				p+=wtx[c][i]*sample.weight;
		}
		p-=s*sample.weight;
	}
	return p;
}



double lineSearch(const MlDataSet* ds,
				   const vector<size_t>& idxs,
				   const vector< vector<double> >& w,    // the point x
				   const vector< vector<double> >& wtx,
				   const vector< vector<double> >& qtx,
				   const vector< vector<double> >& g,    // the gradient
				   const vector< vector<double> >& q,    // q=p, direction along which we want to minimize f(x)
				   double lambda)
{
	const double TOLX = numeric_limits<double>::epsilon();
	const double ALF  = 1e-4;
	double wtw=0.0, qtw=0.0;
	double qtg = computeDotProduct(q,g);
	double qtq = computeDotProduct(q,q);
	double scale = 100.0/sqrt(qtq);
	if (scale>1.0)
		scale=1.0;

	qtq *= (scale * scale);
	double slope = scale * qtg;
	if (slope<0.0)
		return 0.0;

	if (lambda != 0.0)
	{
		wtw = computeDotProduct(w,w);
		qtw = computeDotProduct(q,w)*scale;
	}
	const double fOld  = computePosterior(ds, idxs, wtx, qtx, wtw, qtw, qtq, lambda, 0.0);
	
	double test=0.0;
	for (size_t c=0; c<q.size(); c++)
		for (size_t i=0; i<q[c].size(); i++)
		{
			double maxArg = fabs(w[c][i]);
			if (maxArg<1.0)
				maxArg=1.0;
			double temp=(scale*(fabs(q[c][i])/maxArg));
			if (temp>test)
				test=temp;
		}


	double alamMin = TOLX/test;
	double alam = 1.0;
	double alam2 = 0.0;
	double f2 = fOld;
	while (true)
	{
	//	double slopeAlam = slope * alam;
		if (alam<alamMin)
			return (alamMin*(scale<1.0 ? scale : 1.0));

		double f=computePosterior(ds, idxs, wtx, qtx, wtw, qtw, qtq, lambda, alam*scale);
		if (f>= fOld + ALF * alam * scale * slope)
		{
			return (alam*scale);
		}

		double tmpAlam;
		if (fabs(alam-1.0)<1e-20)
		{
			tmpAlam = -slope / (2.0 *(f - fOld - slope));
		}
		else
		{
			double r1 = f - fOld - alam*slope;
			double r2 = f2 - fOld - alam2*slope;
			double a = (r1/(alam*alam) - r2/(alam2*alam2))/(alam-alam2);
			double b = (-(alam2*r1)/(alam*alam)+(alam*r2)/(alam2*alam2))/(alam-alam2);
			if (fabs(a)<1e-22)
			{
				tmpAlam=-slope/(2.0*b);
			}
			else
			{
				double disc=b*b-3*a*slope;
				if (disc<0)
				{
					tmpAlam = 0.5 * alam;
				}
				else if (b<=0.0)
				{
					tmpAlam = (sqrt(disc)-b)/(3.0*a);
				}
				else
					tmpAlam = -slope / (b + sqrt(disc));
			}
			if (alam * 0.5 < tmpAlam)
				tmpAlam = alam *0.5;
		}
		
		alam2 = alam;
		f2 = f;
		alam = 0.1 * alam;
		if (alam<tmpAlam)
			alam = tmpAlam;
		
	}
	return (alam*scale);//
}



/*******************************************************************

********************************************************************/
void MlMaximumEntropyModel::learnLMBFGS(MlTrainingContainer* params)
{
	params->initialize();

	const MlDataSet* trainingDataSet = params->getTrainingDataSet();
	const MlDataSet* testDataSet     = params->getTestDataSet();
	const vector<size_t>& trainingIdxs = params->getTrainingIdxs();
	const vector<size_t>& testIdxs	   = params->getTestIdxs();
	const size_t numSamples = trainingIdxs.size();
	const bool performTest = (testDataSet && testIdxs.size()>0);

	if (trainingDataSet->getNumClasess() < 2)
		error("learnLMBFGS accepts only datasets with 2 or more classes, yor data has ",
			trainingDataSet->getNumClasess());

	const double  lambda  = params->getLambda();
	const size_t memorySize = params->getLmBfgsMemorySize(); 
	const size_t reportFrequency = params->getVerboseLevel();
	const double perplexityDelta  = params->getPerplexityDelta();
	const size_t numClasses  = trainingDataSet->getNumClasess();
	const size_t numFeatures = trainingDataSet->getNumBasicFeatures(); // F
	const size_t numTraining = trainingIdxs.size();    // N

	const size_t numRestarts=3;

	// data structures used for training
	vector< vector<double> >& w = weights_; // class X features
	vector< vector<double> > wOld(numClasses, vector<double>(numFeatures,0.0)); // class X features
	vector< vector<double> > wtx(numClasses, vector<double>(numTraining, 0.0));  // class X samples
	vector< vector<double> > qtx(numClasses, vector<double>(numTraining, 0.0));  // class X samples
	vector< vector<double> > q(numClasses, vector<double>(numFeatures,0.0));       // class X features
	vector< vector<double> > g(numClasses, vector<double>(numFeatures,0.0));       // class X features
	vector< vector<double> > gOld(numClasses, vector<double>(numFeatures,0.0));  // class X features

	vector< vector<float> > trainingProbs(numClasses, vector<float>(numTraining));
	vector< vector<float> > testProbs(numClasses, vector<float>(numTraining));
	vector< vector<double> > bestW(numClasses, vector<double>(numFeatures));

	// initialize weights
	if (params->getInputPath().length() > 1)
	{
		const string modelFile = params->getInputPath() + "_scr.txt";
		if (readModel(modelFile.c_str()))
			params->setIndClearWeights(false);
	}

	if (params->getIndClearWeights())
		weights_.clear();

	weights_.resize(numClasses, vector<double>(numFeatures,0.0));

	double previousPerplexity = MAX_FLOAT;
	float  bestTestError=1.0;
	size_t bestTestRound=0;
	float  bestTrainingError=1.0;
	size_t bestTrainingRound=0;

	bool terminateTraining = false;
	size_t totalRounds=0;
	size_t megaRound=0;
	for ( megaRound=0; megaRound<numRestarts; megaRound++)
	{
		// first round
		computeGradient(trainingDataSet, trainingIdxs, w, wtx, lambda, g);

		const double gtg = computeDotProduct(g,g);
		const double denominator = 1.0 / sqrt(gtg);
		for (size_t c=0; c<numClasses; c++)
			for (size_t i=0; i<numFeatures; i++)
				q[c][i]=g[c][i]*denominator;

		// qtx <- qTx
		for (size_t c=0; c<numClasses; c++)
			for (size_t i=0; i<numSamples; i++)
			{
				const MlSample& sample = trainingDataSet->getSample(trainingIdxs[i]);
				qtx[c][i]=computeDotProduct(q[c],sample.pairs);
			}

		// eta <- lineSearch(...)
		double eta = lineSearch(trainingDataSet, trainingIdxs, w, wtx, qtx, g, q, lambda);
		//cout << "eta = " << eta << endl;

		// update wtx <- wtx + eta*qtx
		for (size_t c=0; c<numClasses; c++)
			for (size_t i=0; i<wtx[c].size(); i++)
				wtx[c][i]+=eta*qtx[c][i];

		// update wOld<- w ; w <- w + eta *q ; gOld<-g
		for (size_t c=0; c<numClasses; c++)
		{
			memcpy(&wOld[c][0],&w[c][0],sizeof(double)*w[c].size());
			memcpy(&gOld[c][0],&g[c][0],sizeof(double)*g[c].size());
			for (size_t i=0; i<numFeatures; i++)
				w[c][i]+= eta*q[c][i];
		}


		// initialize memory
		vector< vector< vector<double> > > memoryU(memorySize, vector< vector<double> >(numClasses));
		vector< vector< vector<double> > > memoryD(memorySize, vector< vector<double> >(numClasses));
		vector< double > memoryAlpha(memorySize);
		size_t nextMemPosition=0;
		size_t numMemPushes=0;
		
		// iterate until convergence
		size_t round=1;
		while (round<10000)
		{
			// compute errors and report round results
			{
				double trainingLogLikelihood=0.0, testLogLikelihood=0.0;
				const double trainingError = calcErrorRateWithLogLikelihood(trainingDataSet, trainingIdxs,
																		false, &trainingLogLikelihood);
				double testError=1.0;
				if (performTest)
					testError = calcErrorRateWithLogLikelihood(testDataSet, testIdxs, false, &testLogLikelihood);

				if (reportFrequency>0 && round % reportFrequency == 0)
				{
					cout << round << "\t" << scientific << setprecision(5) << trainingLogLikelihood << "\t" << fixed << setprecision(5) << trainingError;
					if (performTest)
						cout <<"\t" << scientific << testLogLikelihood << "\t" << fixed << setprecision(5)<< testError;
					cout << endl;
				}
				
				if (performTest)
				{
					if (testError<=bestTestError)
					{
						bestTestRound=round;
						bestTestError=testError;
						for (size_t c=0; c<numClasses; c++)
							memcpy(&bestW[c][0],&w[c][0],numFeatures*sizeof(double)); // copy weights
					}
				}
				
				if (trainingError<=bestTrainingError)
				{
					bestTrainingRound=round;
					bestTrainingError=trainingError;
					if (! performTest)
					{
						for (size_t c=0; c<numClasses; c++)
							memcpy(&bestW[c][0],&w[c][0],numFeatures*sizeof(double)); // copy weights
					}
				}		
			}

			// Train new round

			computeGradient(trainingDataSet, trainingIdxs, w, wtx, lambda, g);

			double alpha=0.0;
			double sigma=0.0;
			double utu=0.0;

			// write u=g'-g and d=w'-w onto memory, use them to compute alpha and sigma
			vector< vector<double> >& u = memoryU[nextMemPosition];
			vector< vector<double> >& d = memoryD[nextMemPosition];
			for (size_t c=0; c<numClasses; c++)
			{
				const size_t numFeatures = g[c].size();
				u[c].resize(numFeatures);
				d[c].resize(numFeatures);
				for (size_t i=0; i<numFeatures; i++)
				{
					const double gDiff = g[c][i]-gOld[c][i];
					const double wDiff = w[c][i]-wOld[c][i];
					u[c][i]=gDiff;
					d[c][i]=wDiff;
					alpha += gDiff*wDiff;
					utu += gDiff*gDiff;
				}
			}
			sigma = alpha / utu;
			memoryAlpha[nextMemPosition]=alpha;

			// update memory position
			nextMemPosition++;
			if (nextMemPosition == memorySize)
				nextMemPosition = 0;
			numMemPushes++;

			// q<-g
			for (size_t c=0; c<numClasses; c++)
				memcpy(&q[c][0],&g[c][0],g[c].size()*sizeof(double));
			
			// determine memory evaluation order 1..M (M is the newest)
			vector<size_t> memOrder;
			if (numMemPushes<=memorySize)
			{
				for (size_t i=0; i<numMemPushes; i++)
					memOrder.push_back(i);
			}
			else
			{
				for (size_t i=0; i<memorySize; i++)
					memOrder.push_back((i+nextMemPosition) % memorySize);
			}

			vector<double> beta(memOrder.size(),0.0);
			for (int i=memOrder.size()-1; i>=0; i--)
			{
				const size_t m = memOrder[static_cast<size_t>(i)];
				const double alpha = memoryAlpha[m];
				
				const vector< vector<double> >& dM = memoryD[m];
				double& betaM = beta[m];

				// compute beta[m] = (memory_d[m] dot g)/alpha[m]
				for (size_t c=0; c<dM.size(); c++)
					for (size_t i=0; i<dM[c].size(); i++)
						betaM += dM[c][i]*g[c][i];
				betaM/=alpha;
				
				// q <- q - beta[m]*memory_u[m]
				const vector< vector<double> >& uM = memoryU[m]; 
				for (size_t c=0; c<q.size(); c++)
					for (size_t i=0; i<q[c].size(); i++)
						q[c][i] -= betaM * uM[c][i];

			}

			// q <- sigma*q
			for (size_t c=0; c<q.size(); c++)
				for (size_t i=0; i<q[c].size(); i++)
					q[c][i]*=sigma;


			for (size_t i=0; i<memOrder.size(); i++)
			{
				const size_t m = memOrder[static_cast<size_t>(i)];
				const vector< vector<double> >& uM = memoryU[m];
				const vector< vector<double> >& dM = memoryD[m]; 
				const double betaM = beta[m];
				const double oneOverAlpha = 1.0 / memoryAlpha[m];
				double umq = computeDotProduct(uM,q);
				for (size_t c=0; c<numClasses; c++)
					for (size_t j=0; j<q[c].size(); j++)
					{
						const double dq = dM[c][j] * (betaM - umq*oneOverAlpha);
						umq += uM[c][j]*dq;
						q[c][j] += dq;
					}
		
			}

			// q<- -q
			for (size_t c=0; c<numClasses; c++)
				for (size_t i=0; i<q[c].size(); i++)
					q[c][i]=-q[c][i];
			
			// qtx = q*X
			for (size_t i=0; i<trainingIdxs.size(); i++)
			{
				const MlSample& sample = trainingDataSet->getSample(trainingIdxs[i]);
				for (size_t c=0; c<numClasses; c++)
					qtx[c][i]=computeDotProduct(q[c],sample.pairs);
			}
			

			bool needToRestart=false;
			eta = lineSearch(trainingDataSet, trainingIdxs, w, wtx, qtx, g, q, lambda);
			if (eta<= 0.0)
			{
				// restart ?
				needToRestart = true;
			}

			// update wOld<- w ; w <- w + eta *q ; gOld<- g
			for (size_t c=0; c<numClasses; c++)
			{
				memcpy(&wOld[c][0],&w[c][0],sizeof(double)*w[c].size());
				memcpy(&gOld[c][0],&g[c][0],sizeof(double)*g[c].size());
				for (size_t i=0; i<numFeatures; i++)
					w[c][i]+= eta*q[c][i];
			}

			for (size_t c=0; c<numClasses; c++)
				for (size_t i=0; i<numSamples; i++)
					wtx[c][i]+=eta*qtx[c][i];

			round++;
			totalRounds++;
			if (terminateTraining || needToRestart)
				break;
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

	cout << "W=" << endl;
	printVector(weights_);
	cout << endl;

	cout << "Terminated after " << totalRounds << " rounds (" << megaRound << " restarts)" << endl;
	cout << "Best training error  " << fixed << setprecision(8) << bestTrainingError << " (round " << bestTrainingRound << ")" << endl;
	if (performTest)
	cout << "Best test error      "  << bestTestError     << " (round " << bestTestRound << ")" << endl;

	indWasInitialized_ = true;

	//this->calcErrorRateWithPerplexity(trainingDataSet, trainingIdxs, true, NULL);
}



void  MlMaximumEntropyModel::trainModel(MlTrainingContainer* params)
{
	size_t alg = params->getTrainingAlgorithm();
	if (alg == SMT_MAXIMUM_ENTROPY_LMBFGS)
	{
		learnLMBFGS(params);
		return;
	}

	if (alg == SMT_MAXIMUM_ENTROPY_CG_FR || alg == SMT_MAXIMUM_ENTROPY_CG_PR)
	{
		learnCG(params);
		return;
	}

	error("Unrecognized model training type for maximum entropy: ",alg);
}


double	MlMaximumEntropyModel::calcScoreForClass(const MlSample* sample, int label) const
{
	assert( label>=0 && label < static_cast<int>(weights_.size()));

	vector<double> wtx(weights_.size(),0.0);
	for (size_t c=0; c<weights_.size(); c++)
		for (size_t i=0; i<sample->pairs.size(); i++)
			wtx[c] += static_cast<double>(weights_[c][sample->pairs[i].index] * sample->pairs[i].value);

	double numerator = exp(wtx[static_cast<size_t>(sample->label)]);
	double Z=0.0;
	for (size_t c=0; c<wtx.size(); c++)
		Z+=exp(wtx[c]);
	
	return (numerator/Z);
}

void MlMaximumEntropyModel::calcScores(const MlSample* sample, vector<float>& scores) const
{
	vector<double> wtx(weights_.size(),0.0);
	for (size_t c=0; c<weights_.size(); c++)
		for (size_t i=0; i<sample->pairs.size(); i++)
			wtx[c] += static_cast<double>(weights_[c][sample->pairs[i].index] * sample->pairs[i].value);

	vector<double> e(weights_.size());
	double Z=0.0;
	for (size_t c=0; c<weights_.size(); c++)
	{
		e[c]=exp(wtx[c]);
		Z+=e[c];
	}
	
	scores.resize(weights_.size());
	for (size_t c=0; c<weights_.size(); c++)
		scores[c]=e[c]/Z;
}





// assumes probaiblities were calculated for class 1
double	MlMaximumEntropyModel::calcErrorRateWithLogLikelihood(const   MlDataSet* mld, 
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

	vector<float> scores;
	for (size_t i=0; i<idxs.size(); i++)
	{
		const MlSample& sample = mld->getSample(idxs[i]);
		calcScores(&sample, scores);

		size_t maxIdx=0;
		for (size_t j=1; j<scores.size(); j++)
			if (scores[j]>scores[maxIdx])
				maxIdx=j;

		weightTotal += sample.weight;
		if (static_cast<size_t>(sample.label) == maxIdx)
			weightCorrect+=sample.weight;
		
		if (logLikelihood)
			*logLikelihood += (sample.weight * log(scores[static_cast<size_t>(sample.label)]));

		if (verbose)
		{
			cout << idxs[i] << "\t" << sample.label;
			for (size_t c=0; c<scores.size(); c++)
				cout << "\t" << scores[c];
			cout << endl;
		}
	}

	if (logLikelihood)
		*logLikelihood /= weightTotal;

	return (1.0 - (weightCorrect/weightTotal));
}

double	MlMaximumEntropyModel::calcErrorRate(const MlDataSet* mld, bool verbose) const
{
	const vector<MlSample>& samples = mld->getSamples();
	double weightCorrect=0.0;
	double weightTotal=0.0;
	vector<float> scores(weights_.size());
	for (size_t i=0; i<samples.size(); i++)
	{
		const MlSample& sample = samples[i];
		calcScores(&sample, scores);

		size_t maxIdx=0;
		for (size_t j=1; j<scores.size(); j++)
			if (scores[j]>scores[maxIdx])
				maxIdx=j;

		weightTotal += sample.weight;
		if (static_cast<size_t>(sample.label) == maxIdx)
			weightCorrect+=sample.weight;

		if (verbose)
		{
			cout << i;
			for (size_t c=0; c<scores.size(); c++)
				cout << "\t" << scores[c];
			cout << endl;
		}
	}
	return (1.0 - (weightCorrect/weightTotal));
}


bool MlMaximumEntropyModel::readModel(const char* path)
{
	ifstream ifs(path);
	if (! ifs.good())
		return false;

	bool retVal = readModel(ifs);
	ifs.close();
	return retVal;
}




bool MlMaximumEntropyModel::readModel(ifstream& ifs)
{
	char buffer[256];
	while (ifs.good())
	{
		ifs.getline(buffer,256);
		if (ifs.gcount()>0 && buffer[0]!='#')
			break;
	}
		
	unsigned int numClasses=0;
	if (sscanf(buffer,"MAXIMUM_ENTROPY %u",&numClasses) != 1)
	{
		cout << "Bad line in model file:" << endl << buffer << endl;
		return false;
	}

	weights_.resize(numClasses);
	for (size_t c=0; c<numClasses; c++)
	{
		ifs.getline(buffer,256);
		unsigned int numWeights=0;
		if (sscanf(buffer,"%u",&numWeights) != 1)
		{
			cout << "Bad line in model file:" << endl << buffer << endl;
			return false;
		}
		weights_[c].resize(numWeights,0.0);
		while (ifs.good())
		{
			ifs.getline(buffer,256);
			if (! strncpy(buffer,"END_",4))
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
			if (index>weights_[c].size())
				error("Bad feature index in line: ",buffer);

			weights_[c][index]=weight;
		}
	}
	return true;
}



bool	MlMaximumEntropyModel::writeModel(ostream& os) const
{
	os << "MAXIMUM_ENTROPY " << weights_.size() << endl;
	for (size_t c=0; c<weights_.size(); c++)
	{
		os << weights_[c].size() << endl;
		for (size_t i=0; i< weights_[c].size(); i++)
		{
			if (weights_[c][i] == 0.0)
				continue;

			char buffer[20];
			sprintf(buffer,"%e",weights_[c][i]);
			os << "F" << i << "\t" << buffer << endl;
		}
		os << "END_" << c << endl;
	}
	return true;
}

void MlMaximumEntropyModel::outputModelAnalysis(const MlDataSet* mld, 
							 const vector<size_t>& idxs,
							 const MlFeatureSet* featureSet,
							 ostream& os) const
{
	cout << "Analysis? Need to write function that does it..." << endl;
}

