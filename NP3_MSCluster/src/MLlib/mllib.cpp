#include "mllib.h"
#include "mllogistic.h"


// --training-set ..\PepNovo\DataFiles\LTQ_TRYP_SQS_train_0_2_1.txt --max-iterations 2000 --verbose-level 10 --test-ratio 0.3 --input-dir NewModels --input-name LTQ_TRYP_SQS_ORG --output-name SQS_021 --generation-rounds 0


void printUsage(char *message)
{
	cout << "MLlib - Machine Learning library. Created by Ari Frank" << endl << endl;

	cout << "For model training:" << endl;
	cout << "-------------------" << endl;
	cout << "--algorithm <LOGISTIC_CG|...> (default = LOGISTIC_CG), available model types:" << endl;
	printModelTypes();
	cout << endl;
	cout << "--input-dir  <path> (where input model files are located)" << endl;
	cout << "--input-name <name> (name of input model files (without suffixes like \"_mod.txt\")" << endl;
	cout << "--output-dir <path> (where ouput model files should go)" << endl;
	cout << "--output-name <name> (name of output model files (without suffixes like \"_mod.txt\")" << endl;
	cout << "--training-set <path> (path of data files to be read, can use multiple instances of --training-set)" << endl;
	cout << "--test-set     <path> (optional, path of data files to be read, can use multiple instances of --test-set)" << endl;
	cout << "--test-ratio   <0-0.5> (optional)" << endl;
	cout << "--max-iterations <N> (optional, default 10000)" << endl;
	cout << "--perplexity-delta <0-0.1> (optional, default 1E-6)" << endl;
	cout << "--lambda <X> (regularization parameter, optional, default X=0)" << endl;
	cout << "--lmbfgs-memory <X> (optional, 0<=X<=100, default X=5)" << endl;
	cout << "--verbose-level <N> (optional - reports progress every N iterations of trianing, default 10, for no reports use N=0)" << endl;
	
	cout << endl;
	cout << "Additional flags for model expansion (feature generation):" << endl;
	cout << "----------------------------------------------------------" << endl;
	cout << "--generation-rounds <0-10> - max number of generation rounds (feault N=0, only initial feature expansion is performed)" << endl;
	cout << "--bhattachryya <D>    (optional) - the minmal distance D for which a unary transformation will be performed" << endl;
	cout << endl;

	cout << endl;
	cout << "For prediction:" << endl;
	cout << "---------------" << endl;
	cout << "--calsify <model file> <data path>" << endl;
	cout << endl;



	if (message)
		cout << "ERROR: " << message << endl;

#ifdef WIN32
	system("pause");
#endif
	exit(0);
}



double classifyDataFile(const char* modelPath, const char* dataFile, bool verbose)
{
	MlModel mlModel;

	mlModel.readModel(modelPath);
	const MlScoreModel* mlScoreModel = mlModel.getScoreModel();

	if (! mlScoreModel)
		error("couldn't read score model, path=",modelPath);

	MlDataSet mld;
	mld.readDataFile(dataFile);

	const double errorRate = mlScoreModel->calcErrorRate(&mld, verbose);
	return (1.0 - errorRate);
}



int main(int argc, char **argv)
{
	char modelPath[256], dataPath[256]; // for classify
	MlTrainingContainer params;
	MlModel			    model;

	params.setTrainingAlgorithm(SMT_LOGISTIC_CG);

	if (argc <= 1)
		printUsage();

	seedRandom(112233U);
	int i=0;
	while (++i<argc)
	{
		if (! strcmp(argv[i],"--classify"))
		{
			if (++i==argc)
				printUsage("must supply model path and data file!");

			if (sscanf(argv[i++],"%s",modelPath) != 1)
				printUsage("must supply model path and data file!");

			if (sscanf(argv[i++],"%s",dataPath) != 1)
				printUsage("must supply model path and data file!");

			double accuracy = classifyDataFile(modelPath, dataPath);
			if (accuracy>=0)
				cout << "Accuracy: " << accuracy << endl;

			#ifdef WIN32
				system("pause");
			#endif
			return 0;
		}


		if (! strcmp(argv[i],"--algorithm"))
		{
			if (++i==argc)
				printUsage("must supply algorithm type!");

			size_t algorithmType = getModelTypeFormLabel(argv[i]);
			if (algorithmType == MAX_UINT)
				error("unknown algorithm name: ",argv[i]);

			params.setTrainingAlgorithm(algorithmType);
			continue;	
		}

		if (! strcmp(argv[i],"--training-set"))
		{
			if (++i==argc)
				printUsage("must supply training set path!");
			params.addTrainingSetPath(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--test-set"))
		{
			if (++i==argc)
				printUsage("must supply test set path!");
			params.addTestSetPath(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--test-ratio"))
		{
			if (++i==argc)
				printUsage("must supply test ratio!");
			float r=0.0;
			if (sscanf(argv[i],"%f",&r) != 1 || r < 0.0 || r > 0.5)
				printUsage("must supply a test set ratio 0-0.5!");
			params.setTestRatio(r);
			continue;
		}

		if (! strcmp(argv[i],"--max-iterations"))
		{
			if (++i==argc)
				printUsage("must supply maximal number of iterations!");
			size_t n=10000;
			if (sscanf(argv[i],"%u",&n) != 1 || n<5)
				printUsage("must supply maximal number of iterations N>=5!");
			params.setMaxNumIterations(n);
			continue;
		}

		if (! strcmp(argv[i],"--perplexity-delta"))
		{
			if (++i==argc)
				printUsage("must supply perplexity difference to terminate training!");

			double pd;
			if (sscanf(argv[i],"%lf",&pd) != 1 || pd <0 || pd > 0.1)
				printUsage("must supply perplexity difference to terminate training (value between 0-0.1)");
			params.setPerplexityDelta(pd);
			continue;
		}

		if (! strcmp(argv[i],"--lambda"))
		{
			if (++i==argc)
				printUsage("must supply lambda value!");

			double lambda=0;
			if (sscanf(argv[i],"%lf",&lambda) != 1 || lambda <0)
				printUsage("must supply lambda>= 0");
			params.setLambda(lambda);
			continue;
		}

		if (! strcmp(argv[i],"--lmbfgs-memory"))
		{
			if (++i==argc)
				printUsage("must supply lmbfgs size!");

			size_t memorySize=MAX_UINT;
			if (sscanf(argv[i],"%u",&memorySize) != 1 || memorySize>100)
				printUsage("must supply lmbfgs size 0-100");
			params.setLmBfgsMemorySize(memorySize);
			continue;
		}


		if (! strcmp(argv[i],"--verbose-level"))
		{
			if (++i==argc)
				printUsage("must supply verbose level!");
			size_t v=0;
			if (sscanf(argv[i],"%u",&v) != 1)
				printUsage("must supply verbose level!");
			params.setVerboseLevel(v);
			continue;
		}

		if (! strcmp(argv[i],"--input-dir"))
		{
			if (++i==argc)
				printUsage("must supply path to output model");
			params.setInputDir(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--input-name"))
		{
			if (++i==argc)
				printUsage("must supply path to initial model");
			params.setInputName(argv[i]);
			continue;
		}


		if (! strcmp(argv[i],"--output-dir"))
		{
			if (++i==argc)
				printUsage("must supply path to output model");
			params.setOutputDir(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--output-name"))
		{
			if (++i==argc)
				printUsage("must supply path to output model");
			params.setOutputName(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--generation-rounds"))
		{
			if (++i==argc)
				printUsage("must supply number of feature generation rounds");
			size_t g=0;
			if (sscanf(argv[i],"%u",&g) != 1 || g>10)
				printUsage("must supply number of generation rounds 0-10!");
			params.setNumGenerationRounds(g);
			continue;
		}

		if (! strcmp(argv[i],"--bhattacharyya"))
		{
			if (++i==argc)
				printUsage("must supply distance>=0");
			double d;
			if (sscanf(argv[i],"%lf",&d) != 1 || d<0.0)
				printUsage("must supply distance>=0!");
			params.setBhattachryyaDistance(d);
		}


		cout << "Unkown flag: " << argv[i] << endl;
		printUsage();
	}
	
	// if we reached this point, we need to create a model
	params.initialize( params.getVerboseLevel()>0 );
	model.createModel(params);
	

#ifdef WIN32
	system("pause");
#endif
	return 0;
}

