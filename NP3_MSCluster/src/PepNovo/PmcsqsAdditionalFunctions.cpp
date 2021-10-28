#include "PmcsqsAdditionalFunctions.h"
#include "../Common/BufferedOutputFile.h"
#include "../Common/auxfun.h"

void filterDataSet(AllScoreModels* model, 
				   const string& list, 
				   const string& outDir, 
				   float sqsThreshold,
				   const string& suffix)
{
	Config* config = model->get_config();
	PMCSQS_Scorer* pmcsqsModel = const_cast<PMCSQS_Scorer*>(model->get_pmcsqs_ptr());
	if (sqsThreshold>0.0 && ! pmcsqsModel->getIndInitializedSqs())
		error("Sqs model not initialized!, need a valid sqs model if using a filtering threshold!");

	if (sqsThreshold>0.0)
	{
		cout << "Filtering spectra with SQS threshold of " << sqsThreshold << endl;
	}
	else
		cout << "Warning: no SQS threhsold supplied, all spectra will be outputted!" << endl;

	vector<string> paths;
	readListOfPaths(list.c_str(), paths);
	
	cout << "Read " << paths.size() << " paths to spectra files." << endl;
	
	if (paths.size() == 0)
		return;

	size_t peakBufferSize = 10000;
	Peak*  peakBuffer = new Peak[peakBufferSize];
	
	size_t totalValidSpectraFound  = 0; // spectra that have the right mslevel, a minimum number of peaks, etc.
	size_t totalGoodSpectraWritten = 0; // spectra that passed the sqs threshold

	char* mgfBuffer = new char[131072]; // should be enough for any MGF
	
	for (size_t i=0; i<paths.size(); i++)
	{
		size_t numValidSpectraFound  = 0;
		size_t numGoodSpectraWritten = 0;
		
		cout << i << "\tExtracting from: " << paths[i] << endl;
		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(paths[i].c_str(), config);

		SpectraList sl(sa);
		sl.selectAllAggregatorHeaders();

		cout << "\tFound " << sl.getNumHeaders() << " spectra...";
		if (sl.getNumHeaders() <= 0)
			continue;

		BufferedOutputFile outputMgf;
		string fileName;
		getFileNameWithoutExtension(paths[i].c_str(), fileName);
		const string outputMgfPath = outDir + "/" + fileName + "_" + suffix + ".mgf";
		outputMgf.open(outputMgfPath, 1<<20);
	
		for (size_t j=0; j<sl.getNumHeaders(); j++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(j);
			if (header->getOriginalNumPeaks()>= peakBufferSize)
			{
				delete [] peakBuffer;
				peakBufferSize = header->getOriginalNumPeaks()*2;
				peakBuffer = new Peak[peakBufferSize];
			}

			PeakList pl;
			pl.setPeaksPtr( peakBuffer );
			pl.readPeaksToBuffer(sa, header, peakBuffer);
			pl.initializePeakList(config, true);

	        // NP3 GOT change few peaks from 5 to 1 @@@ || header->getOriginalNumPeaks()<15
			if (pl.getNumPeaks()<1) // don't bother with spectra with too few peaks
				continue;

			numValidSpectraFound++;
			if (pmcsqsModel && sqsThreshold>0.0)
			{
				size_t maxCharge=0;
				const float sqs = pmcsqsModel->calculateSqsScore(config, pl, &maxCharge);
				if (sqs<sqsThreshold)
					continue;
				header->setSqs(sqs);
			}
			
			pl.setHeaderFileIndexInList(i);
			const size_t n=pl.outputToMgfBuffer(mgfBuffer);
			outputMgf.writeToBuffer(mgfBuffer, n);	
			numGoodSpectraWritten++;
		}

		outputMgf.close();
		cout << " wrote " << numGoodSpectraWritten << "/" << numValidSpectraFound << endl;
		totalGoodSpectraWritten += numGoodSpectraWritten;
		totalValidSpectraFound  += numValidSpectraFound;
	}

	cout << endl << "Done." << endl;
	cout << "Wrote " << totalGoodSpectraWritten << "/" << totalValidSpectraFound
		<< " spectra that passed sqs threshold " << sqsThreshold << endl;
}

				   
void selectSpectra(AllScoreModels* model, 
				   const string& list, 
				   const string& outDir,
				   int selectCharge)
{

	Config* config = model->get_config();
	PMCSQS_Scorer* pmcsqsModel = const_cast<PMCSQS_Scorer*>(model->get_pmcsqs_ptr());
	
	vector<string> paths;
	readListOfPaths(list.c_str(), paths);

	ostringstream coss;
	coss<<"c"<<selectCharge;
	
	cout << "Read " << paths.size() << " paths to spectra files." << endl;
	
	if (paths.size() == 0)
		return;

	size_t peakBufferSize = 10000;
	Peak*  peakBuffer = new Peak[peakBufferSize];
	
	size_t totalValidSpectraFound  = 0; // spectra that have the right mslevel, a minimum number of peaks, etc.
	size_t totalGoodSpectraWritten = 0; // spectra that passed the sqs threshold

	char* mgfBuffer = new char[131072]; // should be enough for any MGF
	
	for (size_t i=0; i<paths.size(); i++)
	{
		size_t numValidSpectraFound  = 0;
		size_t numGoodSpectraWritten = 0;
		
		cout << i << "\tExtracting from: " << paths[i] << endl;
		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(paths[i].c_str(), config);

		SpectraList sl(sa);
		sl.selectAllAggregatorHeaders();

		cout << "\tFound " << sl.getNumHeaders() << " spectra...";
		if (sl.getNumHeaders() <= 0)
			continue;

		BufferedOutputFile outputMgf;
		string fileName;
		getFileNameWithoutExtension(paths[i].c_str(), fileName);
		const string outputMgfPath = outDir + "/" + fileName + "_" + coss.str() + ".mgf";
		outputMgf.open(outputMgfPath, 1<<20);
	
		for (size_t j=0; j<sl.getNumHeaders(); j++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(j);
			if (header->getOriginalNumPeaks()>= peakBufferSize)
			{
				delete [] peakBuffer;
				peakBufferSize = header->getOriginalNumPeaks()*2;
				peakBuffer = new Peak[peakBufferSize];
			}

			PeakList pl;
			pl.setPeaksPtr( peakBuffer );
			pl.readPeaksToBuffer(sa, header, peakBuffer);
			pl.initializePeakList(config, true);

	        // NP3 GOT change few peaks from 5 to 1 @@@
			if (pl.getNumPeaks()<1) // don't bother with spectra with too few peaks
				continue;

			numValidSpectraFound++;
			if (pmcsqsModel)
			{
				size_t maxCharge=0;
				const float sqs = pmcsqsModel->calculateSqsScore(config, pl, &maxCharge);
				if (maxCharge != static_cast<size_t>(selectCharge))
					continue;
				SingleSpectrumHeader* ssh = const_cast<SingleSpectrumHeader*>(header);
				ssh->setCharge(selectCharge);
			}
			
			pl.setHeaderFileIndexInList(i);
			const size_t n=pl.outputToMgfBuffer(mgfBuffer);

            // GOT BASELINE @@@
			/*if (n > 0)
			{*/
            outputMgf.writeToBuffer(mgfBuffer, n);
            numGoodSpectraWritten++;
			//}
		}

		outputMgf.close();
		cout << " wrote " << numGoodSpectraWritten << "/" << numValidSpectraFound << endl;
		totalGoodSpectraWritten += numGoodSpectraWritten;
		totalValidSpectraFound  += numValidSpectraFound;
	}

	cout << endl << "Done." << endl;
	cout << "Wrote " << totalGoodSpectraWritten << "/" << totalValidSpectraFound
		<< " spectra with predicted charge " << selectCharge << endl;
}


struct CrapSpectrumPointer {
	float sqs;
	int fileIdx;
	int scan;
};


void makeCrapDb(AllScoreModels* model, const string& list, const string& name, 
				const string& outDir, float sqsThresh, int numCrapSpectra)
{
	Config* config = model->get_config();
	PMCSQS_Scorer* pmcsqsModel = const_cast<PMCSQS_Scorer*>(model->get_pmcsqs_ptr());
	if (! pmcsqsModel->getIndInitializedSqs())
		error("Sqs model not initialized!, need a valid sqs model if using a filtering threshold!");

	cout << "Selecting spectra with SQS threshold below " << sqsThresh << endl;
	
	vector<string> paths;
	readListOfPaths(list.c_str(), paths);
	
	cout << "Read " << paths.size() << " paths to spectra files." << endl;
	
	if (paths.size() == 0)
		return;

	size_t peakBufferSize = 10000;
	Peak*  peakBuffer = new Peak[peakBufferSize];
	
	size_t totalValidSpectraFound  = 0; // spectra that have the right mslevel, a minimum number of peaks, etc.
	size_t totalGoodSpectraWritten = 0; // spectra that passed the sqs threshold

	char* mgfBuffer = new char[65536]; // should be enough for any MGF
	
	random_shuffle(paths.begin(), paths.end());

	BufferedOutputFile outputMgf;
	string mgfPath = outDir + "/" + name + "_0.mgf";
	int fileCounter=0;
	outputMgf.open(mgfPath, 1<<20);

	
	int totalCrapWrittenToOutput = 0;
	int numCrapWrittenToCurrentOuput = 0;

	for (size_t i=0; i<paths.size(); i++)
	{
		int numForThisFile=0;
		cout << i << "\tExtracting from: " << paths[i] << endl;
		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(paths[i].c_str(), config);

		SpectraList sl(sa);
		sl.selectAllAggregatorHeaders();

		cout << "\tFound " << sl.getNumHeaders() << " spectra...";
		if (sl.getNumHeaders() <= 0)
			continue;

		for (size_t j=0; j<sl.getNumHeaders(); j++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(j);
			if (header->getOriginalNumPeaks()>1e6)
				continue;

			if (header->getOriginalNumPeaks()>= peakBufferSize)
			{
				delete [] peakBuffer;
				peakBufferSize = header->getOriginalNumPeaks()*2;
				peakBuffer = new Peak[peakBufferSize];
			}

			PeakList pl;
			pl.setPeaksPtr( peakBuffer );
			pl.readPeaksToBuffer(sa, header, peakBuffer);
			pl.initializePeakList(config, true);

	        // NP3 GOT change few peaks from 5 to 1 @@@
			if (pl.getNumPeaks()<1) // don't bother with spectra with too few peaks
				continue;

			size_t maxCharge=0;
			const float sqs = pmcsqsModel->calculateSqsScore(config, pl, &maxCharge);
			if (sqs>sqsThresh)
				continue;
			header->setSqs(sqs);

			if (numCrapWrittenToCurrentOuput >= 20000)
			{
				numCrapWrittenToCurrentOuput=0;
				fileCounter++;
				outputMgf.close();
				ostringstream oss;
				oss << outDir << "/" << name << "_" << fileCounter << ".mgf";
				outputMgf.open(oss.str(), 1<<20);
			}
			
			const size_t n=pl.outputToMgfBuffer(mgfBuffer);
			outputMgf.writeToBuffer(mgfBuffer, n);

			numCrapWrittenToCurrentOuput++;
			numForThisFile++;
			if (++totalCrapWrittenToOutput >= numCrapSpectra)
				break;	
		}

		cout << " Wrote " << numForThisFile << " (" << totalCrapWrittenToOutput << ")." << endl;

		if (totalCrapWrittenToOutput >= numCrapSpectra)
			break;
	}

	outputMgf.close();
	cout << endl << "Done." << endl;
	cout << "Wrote " << totalCrapWrittenToOutput
		<< " spectra that have sqs below threshold " << sqsThresh << endl;
	cout << "These were written into " << fileCounter+1 << " files in " << outDir << endl;
}


void benchmarkSqs(AllScoreModels* model, const string& list)
{
	Config* config = model->get_config();
	PMCSQS_Scorer* pmcsqsModel = const_cast<PMCSQS_Scorer*>(model->get_pmcsqs_ptr());
	if (! pmcsqsModel->getIndInitializedSqs())
		error("Sqs model not initialized!, need a valid sqs model if using a filtering threshold!");

	const int maxCharge = pmcsqsModel->getMaximalChargeWithModels();
	const int numSizes  = pmcsqsModel->getNumSizes();
	
	vector<string> paths;
	readListOfPaths(list.c_str(), paths);
	
	cout << "Read " << paths.size() << " paths to spectra files." << endl;
	
	if (paths.size() == 0)
		return;

	vector< vector<double> > sqsScores(numSizes, vector<double>(101,0));
	vector< vector< vector<double> > > chargeAssigns(numSizes);
	for (int i=0; i<numSizes; i++)
		chargeAssigns[i].resize(maxCharge+1, vector<double>(maxCharge+1,0));

	size_t peakBufferSize = 10000;
	Peak*  peakBuffer = new Peak[peakBufferSize];
	int totalHeaders = 0;
	double totalChargeAssigns = 0;
	double totalSqsAssigns = 0;
	for (size_t i=0; i<paths.size(); i++)
	{
		int numForThisFile=0;
		cout << i << "\tExtracting from: " << paths[i];
		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(paths[i].c_str(), config);

		SpectraList sl(sa);
		sl.selectAllAggregatorHeaders();

		cout << "\tFound " << sl.getNumHeaders() << " spectra..." << endl;
		if (sl.getNumHeaders() <= 0)
			continue;
		totalHeaders += sl.getNumHeaders();

		for (size_t j=0; j<sl.getNumHeaders(); j++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(j);
			if (header->getOriginalNumPeaks()>1e6)
				continue;

			if (header->getOriginalNumPeaks()>= peakBufferSize)
			{
				delete [] peakBuffer;
				peakBufferSize = header->getOriginalNumPeaks()*2;
				peakBuffer = new Peak[peakBufferSize];
			}

			PeakList pl;
			pl.setPeaksPtr( peakBuffer );
			pl.readPeaksToBuffer(sa, header, peakBuffer);
			pl.initializePeakList(config, true);

	        // NP3 GOT change few peaks from 5 to 1 @@@
			if (pl.getNumPeaks()<1) // don't bother with spectra with too few peaks
				continue;

			size_t maxCharge=0;
			const float sqs = pmcsqsModel->calculateSqsScore(config, pl, &maxCharge);
			size_t idx = static_cast<size_t>(sqs*100.0 + 0.4999);
			if (idx>=100)
				idx=99;

			const size_t sizeIndex = pmcsqsModel->getSqsSizeIndex(pl.getHeader()->getMOverZ());
			sqsScores[sizeIndex][idx]++;
			totalSqsAssigns++;

			int charge = pl.getHeader()->getCharge();
			if (charge>0)
			{
				chargeAssigns[sizeIndex][charge][maxCharge]++;
				totalChargeAssigns++;
			}
		}
	}
	
	//
	cout << "Computed SQS scores for " << totalSqsAssigns << " spectra from " << totalHeaders 
		<< " a total of header (some spectra had to few peaks)" << endl;

	// create reports
	cout << "Histogram of SQS scores: " << endl;
	cout << "------------------------ " << endl;

	for (int sizeIndex=0; sizeIndex<numSizes; sizeIndex++)
	{
		double totalForSize=0;
		for (size_t i=0; i<sqsScores[sizeIndex].size(); i++)
			totalForSize += sqsScores[sizeIndex][i];

		if (totalForSize == 0)
			continue;

		cout << "Size " << sizeIndex << ":" << endl;
		cout << "SQS\tCount\tCDF" << endl;
		cout << "---\t-----\t---" << endl;
		double cdf=0.0;
		for (size_t i=0; i<=100; i++)
		{
			cdf+=sqsScores[sizeIndex][i];
			double v = (i==100 ? 1.0 : (i+1.0)*0.01);
			if (cdf>0)
				cout << setprecision(3) << v << "\t" << fixed << setprecision(0) << cdf << "\t"
					<< setprecision(5) << cdf/totalForSize << endl;
			size_t add=0;
			if (i>=9)
				add++;
			if (i>=19)
				add+=3;
			if (i>=49)
				add+=5;
			for (size_t j=1; j<=add && i+j<sqsScores[sizeIndex].size(); j++)
				cdf+=sqsScores[sizeIndex][i+j];
			i+=add;
		}
		cout << endl;
	}

	cout << "Histogram for all sizes combined: " << endl;
	cout << "--------------------------------- " << endl;
	double totalSqs=0;
	for (size_t sizeIndex = 0; sizeIndex<numSizes; sizeIndex++)
		for (size_t i=0; i<sqsScores[sizeIndex].size(); i++)
			totalSqs += sqsScores[sizeIndex][i];
		
	cout << "SQS\tCount\tCDF" << endl;
	cout << "---\t-----\t---" << endl;
	double cdf=0.0;
	for (size_t i=0; i<=100; i++)
	{
		for (size_t sizeIndex = 0; sizeIndex<numSizes; sizeIndex++)
			cdf+=sqsScores[sizeIndex][i];
		double v = (i==100 ? 1.0 : (i+1.0)*0.01);
		if (cdf>0)
			cout << setprecision(3) << v << "\t" << fixed << setprecision(0) << cdf << "\t"
				<< setprecision(5) << cdf/totalSqs << endl;
		size_t add=0;
		if (i>=9)
			add++;
		if (i>=19)
			add+=3;
		if (i>=49)
			add+=5;
		for (size_t sizeIndex = 0; sizeIndex<numSizes; sizeIndex++)
			for (size_t j=1; j<=add && i+j<sqsScores[sizeIndex].size(); j++)
				cdf+=sqsScores[sizeIndex][i+j];
		i+=add;
	}
	cout << endl;

	cout << "Found " << totalChargeAssigns << " spectra with charge assignments." << endl;
	if (totalChargeAssigns<=0)
		return;
	
	vector< vector<double> > chargeCounts(numSizes, vector<double>(maxCharge+1,0));
	for (size_t sizeIndex=0; sizeIndex<numSizes; sizeIndex++)
		for (int charge1=0; charge1<=maxCharge; charge1++)
			for (int charge2=0; charge2<=maxCharge; charge2++)
				chargeCounts[sizeIndex][charge1]+=chargeAssigns[sizeIndex][charge1][charge2];
	
	double totalGoodAssigns=0;
	for (size_t sizeIndex=0; sizeIndex<numSizes; sizeIndex++)
	{
		double sum=0.0;
		for (int charge=1; charge<=maxCharge; charge++)
			sum+=chargeCounts[sizeIndex][charge];

		cout << "SIZE " << sizeIndex <<  endl;
		cout << "\tWeight";
		
		for (int charge=1; charge<=maxCharge; charge++)
			cout << "\tCh. " << charge;
		cout << endl;
		double numCorrect=0.0;
		for (int charge1=1; charge1<=maxCharge; charge1++)
		{
			cout << "Ch. " << charge1 << "\t" << setprecision(3) << chargeCounts[sizeIndex][charge1]/sum;
			for (int charge2=1; charge2<=maxCharge; charge2++)
				cout << "\t" << chargeAssigns[sizeIndex][charge1][charge2]/chargeCounts[sizeIndex][charge1];
			cout << endl;
			numCorrect+=chargeAssigns[sizeIndex][charge1][charge1];
		}
		cout << endl << "Accuracy: " << numCorrect/sum << endl << endl;
		totalGoodAssigns += numCorrect;
	}
	cout << "Total charge assignment accuracy: " << totalGoodAssigns/totalChargeAssigns << endl;
}

