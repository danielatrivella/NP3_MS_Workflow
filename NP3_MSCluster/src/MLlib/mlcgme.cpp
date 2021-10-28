#include "mlmaximumentropy.h"
#include <string.h>

struct CGSearchData {
	friend class MlMaximumEntropyModel;
public:
	CGSearchData(MlTrainingContainer* params, vector< vector<double> >& weights) : 
				params_(params), trainingDataSet_(params->getTrainingDataSet()),
				trainingIdxs_(params->getTrainingIdxs()),
				lambda_(0.0), numClasses_(0), numFeatures_(0), numSamples_(0), w_(weights)
				{ init(); }


	void init();


	// uses current w, computes q
	void computeGradientOfW(); 

	// uses current q
	void updateQtx();

	// updates the point w, recomputes wtx
	void udpateW(double x, const vector< vector<double> >& dir);
	
	// this way ax is always updated with fa, etc.
	double func(double x=0.0, const char* msg = 0);

	// barckets the minimum of the function.
	// When function converges, the minimum should be bracketed by ax,bx, and cx
	// such that bx is between ax and cx, and f(bx) < f(ax) and f(bx)< f(cx)
	void bracket(double a, double b);

	void minimize(double tol = 1e-5);

	double linmin();

	const vector< vector<double> >& getGradient() const { return q_; }

	double	ax,bx,cx;		// for bracketing
	double  fa,fb,fc;      // for bracketing

	double  xmin; // result of minimize
	double  fmin;
	
private:

	MlTrainingContainer* params_;
	const MlDataSet*      trainingDataSet_;
	const vector<size_t>& trainingIdxs_;

	double lambda_;		    // regularization parameter
	size_t numClasses_;
	size_t numFeatures_;
	size_t numSamples_;

	vector< vector<double> >& w_;   // the point p in the CG algorithm
	vector< vector<double> > wtx_; // used to store the dot product wTx
	vector< vector<double> > q_;   // the gradient (xi in the algorithm, the direction of the search)
	vector< vector<double> > qtx_;

	double wtw_, qtw_, qtq_; // for internal computations
};


void CGSearchData::init()
{
	ax=bx=cx=fa=fb=fc=xmin=fmin=NON_FLOAT;

	numClasses_  = trainingDataSet_->getNumClasess();
	numFeatures_ = trainingDataSet_->getNumBasicFeatures(); 
	numSamples_  = trainingIdxs_.size();   

	if (w_.size()>0)
	{
		if (w_.size() != numClasses_)
			error("weights do not have same dim as numClasses!");

		for (size_t c=0; c<numClasses_; c++)
			if (w_[c].size()<numFeatures_)
				w_[c].resize(numFeatures_,0.0);

		wtx_.resize(numClasses_);
		for (size_t c=0; c<numClasses_; c++)
		{
			wtx_[c].resize(numSamples_);
			for (size_t i=0; i<numSamples_; i++)
			{
				const MlSample& sample = trainingDataSet_->getSample(trainingIdxs_[i]);
				wtx_[c][i]=computeDotProduct(w_[c],sample.pairs);
			}
		}
	}
	else
	{
		w_.resize( numClasses_, vector<double>(numFeatures_,0.0));
		wtx_.resize( numClasses_, vector<double>(numSamples_,0.0));
	}

	q_.resize( numClasses_, vector<double>(numFeatures_,0.0));
	qtx_.resize( numClasses_, vector<double>(numSamples_,0.0));

	qtq_ = 0.0;
	qtw_ = 0.0;
	wtw_ = 0.0;
}



double CGSearchData::func(double step, const char *msg)
{
	double posterior;

	if (step != 0.0)
	{
		posterior = computePosterior(trainingDataSet_, trainingIdxs_, wtx_, qtx_, wtw_, qtw_, qtq_, lambda_, step); 
	}
	else 
		posterior = computePosterior(trainingDataSet_, trainingIdxs_, wtx_, wtw_,  lambda_);
	

//	if (msg)
//		cout << "\t" << setprecision(8) << step << "\t" << msg << "\t" << posterior << endl;
	
	return -posterior;
}

// uses current w, computes q and qtx
void CGSearchData::computeGradientOfW()
{
	computeGradient(trainingDataSet_, trainingIdxs_, w_, wtx_, lambda_, q_);

/*	for (size_t c=0; c<numClasses_; c++)
		for (size_t i=0; i<numFeatures_; i++)
			q_[c][i] = -q_[c][i];*/
}


void CGSearchData::updateQtx()
{
	for (size_t c=0; c<numClasses_; c++)
		for (size_t i=0; i<numSamples_; i++)
		{
			const MlSample& sample = trainingDataSet_->getSample(trainingIdxs_[i]);
			qtx_[c][i]=computeDotProduct(q_[c],sample.pairs);
		}

	if (lambda_ !=0.0)
	{
		qtq_ = computeDotProduct(q_,q_);
		qtw_ = computeDotProduct(q_,w_);
	}
}

// updates the point w, recomputes wtx
void CGSearchData::udpateW(double x, const vector< vector<double> >& dir)
{
	for (size_t c=0; c<numClasses_; c++)
		for (size_t i=0; i<numFeatures_; i++)
			w_[c][i] -= x*dir[c][i];

	for (size_t c=0; c<numClasses_; c++)
		for (size_t i=0; i<numSamples_; i++)
		{
			const MlSample& sample = trainingDataSet_->getSample(trainingIdxs_[i]);
			wtx_[c][i] = computeDotProduct(w_[c],sample.pairs);
		}

	if (lambda_ != 0.0)
	{
		wtw_ = computeDotProduct(w_,w_);
		qtw_ = computeDotProduct(q_,w_);
	}
}




// barckets the minimum of the function.
// When function converges, the minimum should be bracketed by ax,bx, and cx
// such that bx is between ax and cx, and f(bx) < f(ax) and f(bx)< f(cx)
void CGSearchData::bracket(double a, double b)
{
	const double GOLD = 1.1618034, GLIMIT = 100.0, TINY =1e-20;
	ax = a;
	bx = b;
	if (fa == NON_FLOAT)
		fa = func(ax,"b(a)");

	if (fb == NON_FLOAT)
		fb = func(bx,"b(b)");

	if (fb > fa) // swap
	{
		double t;
		t = ax; ax = bx; bx = t;
		t = fa; fa = fb; fb = t;
	}
	cx =(bx + GOLD * ( bx - ax)); // first guess
	fc = func(cx,"b(c)");

	double fu;
	while (fb > fc)
	{
		double r=(bx - ax)*(fb - fc);
		double q=(bx - cx)*(fb - fa);
		double t = fabs(q-r); // sign ... max
		if (t<TINY)
			t=TINY;
		if (q<r)
			t=-t;
		double u=bx - ((bx - cx)*q - (bx - ax)*r) / (2.0*t);
		double ulim = bx + GLIMIT * (cx - bx);

		if ((bx - u)*(u- cx) > 0.0)
		{
			fu=func(u,"B(u1)");

			if (fu < fc) // got min between b and c
			{
				ax = bx;
				bx = u;
				fa = fb;
				fb = fu;
				return;
			}
			else if ( fu > fb) // got min between a and u
			{
				cx = u;
				fc = fu;
				return;
			}

			u = cx + GOLD * (cx - bx);
			fu = func(u,"B(u2)");
		}
		else if ((cx - u)*(u-ulim) > 0.0)
		{
			fu=func(u,"B");
			if (fu < fc)
			{
				double t= u+GOLD*(u-cx);
				bx = cx; 
				cx = u; 
				u = t;
				fb = fc; 
				fc = fu;
				fu=func(u,"B(u3)");
			}
		}
		else if ((u-ulim)*(ulim-cx) >= 0.0)
		{
			u = ulim;
			fu= func(u,"B(u4)");
		}
		else
		{
			u = cx + GOLD * (cx - bx);
			fu= func(u,"B(u5)");
		}

		ax = bx; bx = cx; cx = u;
		fa = fb; fb = fc; fc = fu;
	}
}

void CGSearchData::minimize(double tol)
{
	const size_t MAXIT = 100;
	const double CGOLD = 0.3819660;
	const double ZEPS = numeric_limits<double>::epsilon()*0.001;
	double a,b,d=0.0,etemp,fu,fv,fw,fx,e=0.0;
	double p,q,r,tol1,tol2,u,v,w,x,xm;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;

	// check if this value needs to be computed, or is there a stored value that mathces it
	fw=fv=fx = (fb != NON_FLOAT ? fb : func(x,"m"));

	for (size_t iter=0; iter<MAXIT; iter++)
	{
		xm=0.5*(a+b);
		tol1 = tol*fabs(x) + ZEPS;
		tol2 = 2.0 * tol1;
		if (fabs(x-xm)<=(tol2 - 0.5*(b-a)))
		{
			fmin = fx;
			xmin = x;
			return;
		}
		if (fabs(e)>tol1)
		{
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q>0.0)
				p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p>= q*(b-x))
			{
				e =  (x>xm ? a-x : b-x);
				d = CGOLD *e;
			}
			else
			{
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
				{
					d=fabs(tol1);
					if (xm>x)
						d=-d;
				}
			}
		}
		else
		{
			e = (x>xm ? a-x : b-x);
			d = CGOLD * e;
		}

		double t = fabs(tol1);
		if (d<0)
			t=-t;
		u = (fabs(d) >= tol1 ? x+d : x+t);
		fu = func(u,"M(u)");

		// special terminate if funciton values are vary similar
		bool terminate =  (2.0*fabs(fu-fx)/(fabs(fu)+fabs(fx)) < 1.0e-7);



		if (fu<=fx)
		{
			if (u>=x)
			{
				a=x;
			}
			else
				b=x;

			v=w; w=x; x=u;
			fv = fw; fw = fx; fx = fu;
		}
		else
		{
			if (u<x)
			{
				a=u;
			}
			else
				b=u;

			if (fu <= fw || w == x)
			{
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			}
			else if (fu <= fv || v == x || v == w)
			{
				v=u;
				fv=fu;
			}
		}

		if (terminate)
		{
			xmin = x;
			fmin = fx;
			return;
		}
	}
	error("too many iterations in minimize!");
}



double CGSearchData::linmin()
{
	// try find a tight bracket that is similar to previous xmin
	ax = bx = cx = fa = fb = fc = NON_FLOAT;
	double a=0.0,b=0.1;
	if (xmin != NON_FLOAT)
	{
		if (xmin>=0.0)
		{
			ax = 0.0;
			fa = fmin;
			b = 3.0 * xmin;
		}
		else
		{
			a = -3.0 * xmin;
			b = 0.0;
			bx = 0.0;
			fb = fmin;
		}
	}
	

	bracket(a, b);
	minimize();


//	cout << "FMIN: " << fmin << endl;

	udpateW(xmin,q_);
	computeGradientOfW();

	return fmin;
}


/************************************************************************
*************************************************************************/
void MlMaximumEntropyModel::learnCG(MlTrainingContainer* params)
{
	params->initialize();

	const MlDataSet* trainingDataSet = params->getTrainingDataSet();
	const MlDataSet* testDataSet     = params->getTestDataSet();
	const vector<size_t>& trainingIdxs = params->getTrainingIdxs();
	const vector<size_t>& testIdxs	   = params->getTestIdxs();
	const bool performTest = (testDataSet && testIdxs.size()>0);

	if (trainingDataSet->getNumClasess() < 2)
		error("learnCG accepts only datasets with 2 or more classes, yor data has ",
			trainingDataSet->getNumClasess());

	const size_t reportFrequency = params->getVerboseLevel();
	const size_t numClasses  = trainingDataSet->getNumClasess();
	const size_t numFeatures = trainingDataSet->getNumBasicFeatures(); // F
	const size_t numTraining = trainingIdxs.size();    // N
	
	const size_t numRestarts=1;
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
	CGSearchData cgsd(params, weights_);

	cgsd.lambda_ =  params->getLambda();

	vector< vector<float> > trainingProbs(numClasses, vector<float>(numTraining));
	vector< vector<float> > testProbs(numClasses, vector<float>(numTraining));
	vector< vector<double> > bestW(numClasses, vector<double>(numFeatures));
	
	double  bestTestError=1.0;
	size_t  bestTestRound=0;
	double  bestTrainingError=1.0;
	double  bestTrainingLikelihood = MIN_FLOAT;
	size_t  bestTrainingRound=0;

	size_t totalRounds=0;
	size_t megaRound=0;

	const size_t maxIterations = (params->getMaxNumIterations() >0 ? params->getMaxNumIterations() : 10000);
	const double EPS = 1.0e-8;
	const double GTOL = 10.e-8;
	const double ftol = 3.0e-8;

	const size_t resetFrequency = ( numFeatures > 20 ? numFeatures/2  : 10);

	double fret,dgg,gg;

	vector< vector<double> >& p = weights_;
	vector< vector<double> >& xi = cgsd.q_;
	vector< vector<double> > g(numClasses, vector<double>(numFeatures,0.0));
	vector< vector<double> > h(numClasses, vector<double>(numFeatures,0.0));
	

	MlTrainingTerminationDecider terminateDecider;
	bool terminate = false;
	for ( megaRound=0; megaRound<numRestarts; megaRound++)
	{
		double fp = cgsd.func(0.0,"I");
		cgsd.computeGradientOfW();
		for (size_t c=0; c<numClasses; c++)
			for (size_t i=0; i<numFeatures; i++)
			{
				g[c][i] = -xi[c][i];
				xi[c][i]=h[c][i]=g[c][i];
			}
		cgsd.updateQtx();

		for (size_t round=0; round<maxIterations; round++)
		{
			totalRounds++;
			fret = cgsd.linmin();
			if (2.0 * fabs(fret-fp) <= ftol*(fabs(fret)+fabs(fp)+EPS))
				break;

			fp = fret;

			double test=0.0;
			double oneOverDen=1.0/(fp>1.0 ? fp : 1.0);
			for (size_t c=0; c<numClasses; c++)
				for (size_t i=0; i<numFeatures; i++)
				{
					double m = fabs(p[c][i]);
					if (m<1.0)
						m=1.0;
					double temp = fabs(xi[c][i])*m*oneOverDen;
					if (temp>test)
						test=temp;
				}

			if (test<GTOL)
			{
				terminate = true;
				break;
			}

			dgg=0.0;
			gg=0.0;

			if (params->getTrainingAlgorithm() == SMT_MAXIMUM_ENTROPY_CG_FR)
			{
				for (size_t c=0; c<numClasses; c++)
					for (size_t i=0; i<numFeatures; i++)
					{
						gg += g[c][i]*g[c][i];
						dgg += xi[c][i] * xi[c][i];
					}
			}
			else
			{
				for (size_t c=0; c<numClasses; c++)
					for (size_t i=0; i<numFeatures; i++)
					{
						gg += g[c][i]*g[c][i];
						dgg += (xi[c][i]+g[c][i]) * xi[c][i];
					}
			}

			if (gg<1e-30)
			{
				terminate = true;
				break;
			}

			// convergence can be quicker if the gradient is reset (to ignore previous CG directions)
			if (round>0 && resetFrequency>0 && round % resetFrequency == 0)
			{
				for (size_t c=0; c<numClasses; c++)
					for (size_t i=0; i<numFeatures; i++)
					{
						g[c][i] = -xi[c][i];
						xi[c][i]=h[c][i]=g[c][i];
					}
			}
			else
			{
				double gam=dgg/gg;
				for (size_t c=0; c<numClasses; c++)
					for (size_t i=0; i<numFeatures; i++)
					{
						g[c][i]=-xi[c][i];
						xi[c][i]=h[c][i]=g[c][i]+gam*h[c][i];
					}
			}

			// update qtx according to the new direction (since xi=q)
			cgsd.updateQtx();

			// test current model status output progress
			{
				double trainingLogLikelihood=MIN_FLOAT, testLogLikelihood=MIN_FLOAT;
				const double trainingError = calcErrorRateWithLogLikelihood(trainingDataSet, trainingIdxs,
																		false, &trainingLogLikelihood);
				float testError=MAX_FLOAT;
				if (performTest)
					testError = calcErrorRateWithLogLikelihood(testDataSet, testIdxs, false, &testLogLikelihood);

				if (reportFrequency>0 && round % reportFrequency == 0)
				{
					cout << round << "\t"  << setprecision(6) << trainingLogLikelihood << "\t" << fixed << setprecision(5) << trainingError;
					if (performTest)
						cout <<"\t\t" << setprecision(6) << testLogLikelihood << "\t" << fixed << setprecision(5)<< testError;
					cout << endl;
				}
				
				if (performTest)
				{
					if (testError<=bestTestError)
					{
						bestTestRound=round;
						bestTestError=testError;
						for (size_t c=0; c<numClasses; c++)
							memcpy(&bestW[c][0],&weights_[c][0],numFeatures*sizeof(double)); // copy weights
					}

					
				}
				
				if (trainingError<=bestTrainingError)
				{
					bestTrainingRound=round;
					bestTrainingError=trainingError;
					if (! performTest)
					{
						for (size_t c=0; c<numClasses; c++)
							memcpy(&bestW[c][0],&weights_[c][0],numFeatures*sizeof(double)); // copy weights
					}
				}

				terminate = terminateDecider.shouldTerminateTraining(totalRounds, 
									trainingLogLikelihood, testLogLikelihood);
			}
			if (terminate)
				break;
		}
		if (terminate)
			break;
	}
	
	if (! params->getIndHadInternalError())
	{
		params->setIndNormalTermination(true);
	}
	else
		cout << "Warning: encountered mathemtical error while training!" << endl;

	weights_ = bestW;
	
/*	cout << "W=" << endl;
	printVector(weights_);
	cout << endl;*/

	cout << "Terminated after " << totalRounds << " rounds (" << megaRound << " restarts)" << endl;
	cout << "Best training error  " << fixed << setprecision(8) << bestTrainingError << " (round " << bestTrainingRound << ")" << endl;
	if (performTest)
	cout << "Best test error      "  << bestTestError     << " (round " << bestTestRound << ")" << endl;

	indWasInitialized_ = true;
}

