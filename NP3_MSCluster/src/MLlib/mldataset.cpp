#include "mldata.h"
#include "mlauxfun.h"
#include "mlfeature.h"
#include "../Common/auxfun.h"
#include <assert.h>

value_t MlSample::getValue(size_t featureIdx) const
{
	size_t i;
	for (i=0; i<pairs.size(); i++)
		if (pairs[i].index>=featureIdx)
			break;
	
	if (i<pairs.size() && pairs[i].index == featureIdx)
		return pairs[i].value;
	return NON_FLOAT;
}

void MlSample::writeToStream(ostream& os, bool indWriteAdditionalInfo) const
{
	os << ">\t" << label << "\t" << scientific << weight;
	if (indWriteAdditionalInfo)
		os << "\t" << groupIndex << "\t" << rankInGroup << "\t" << tag1 << "\t" <<  tag2 
		   << "\t" << tag3 << "\t" << floatTag;
	os << endl << fixed << setprecision(6) << pairs.size();
	for (size_t i=0; i<pairs.size(); i++)
		os << "\t" << pairs[i].index << "\t" << pairs[i].value;

	os << endl << derivedPairs.size();
		for (size_t i=0; i<derivedPairs.size(); i++)
			os << "\t" << derivedPairs[i].index << "\t" << derivedPairs[i].value;
	os << endl;
}



bool MlSample::readFromStream(istream& is)
{
	static char buffer[128];
	
	while (! is.eof() && is.getline(buffer,128) && buffer[0] != '>') ;
	istringstream iss(buffer+1);
	iss >> label >> weight >> groupIndex >> rankInGroup >> tag1 >>  tag2 
		>> tag3  >> floatTag;

	size_t numPairs = MAX_SIZE_T, numDerivedPairs = MAX_SIZE_T;

	is >> numPairs;
	assert (numPairs<MAX_SIZE_T);
	pairs.resize(numPairs);
	for (size_t i=0; i<numPairs; i++)
		is >> pairs[i].index >> pairs[i].value;

	is >> numDerivedPairs;
	assert(numDerivedPairs<MAX_SIZE_T);
	derivedPairs.resize(numDerivedPairs);
	for (size_t i=0; i<numDerivedPairs; i++)
		is >> derivedPairs[i].index >> derivedPairs[i].value;

	return true;
}


MlDataSet::~MlDataSet()
{
	if (outStream_ && bufferPos_>0)
		closeOutputFile();

	if (buffer_)
		delete [] buffer_;
}

void MlDataSet::clear()
{
	if (indOutStreamOpen_)
	{
		flushOutputBuffer();
		closeOutputFile();
	}
	
	totalWeight_=0.0;
	numClasses_=0;
	numPairs_=0;
	numDerivedPairs_=0;
	maxFeatureIndex_=0;
	indStatsValid_=false; 
	indOutStreamOpen_=false; 
	outStream_=0;
	bufferPos_=0; 
	
	classWeights_.clear(); 
	samples_.clear();
}


void MlDataSet::addSamples(const vector<MlSample>& otherSamples)
{
	const size_t originalSize = samples_.size();
	samples_.resize(originalSize + otherSamples.size());

	for (size_t i=0; i<otherSamples.size(); i++)
		samples_[originalSize+i]=otherSamples[i];

	indStatsValid_ = false;
}


// adds the samples and sets the labels to the new value
void  MlDataSet::addSamples(int label, const vector<MlSample>& otherSamples) 
{
	const size_t originalSize = samples_.size();
	samples_.resize(originalSize + otherSamples.size());

	for (size_t i=0; i<otherSamples.size(); i++)
	{
		samples_[originalSize+i]= otherSamples[i];
		samples_[originalSize+i].label = label;
	}

	indStatsValid_ = false;
}


// tries to reweight data so the samples have the same ratio distribution
void MlDataSet::reweight(const vector<weight_t>& ratios)
{
	if (! indStatsValid_)
		tallyClassStatistics();

	if (ratios.size() < numClasses_)
	{
		cout << "Warning: supplied ratios for " << ratios.size() << " classes, dataset has " << numClasses_ << endl;
		cout << "	      reweighting not performed!" << endl;
		return;
	}

	vector<double> adjustedRatios(numClasses_,0.0);
	vector<double> currentRatios(numClasses_,0.0);
	vector<double> multFactors(numClasses_,0.0);

	double sum=0;
	for (size_t i=0; i<ratios.size(); i++)
	{
		if (ratios[i]<0.0)
			error("Negative ratio supplied to reweighting!");
		sum+=ratios[i];
	}

	if (sum<=0.0)
		error("ratio weights less than or equal to zero!");

	for (size_t i=0; i<numClasses_; i++)
	{
		adjustedRatios[i]=ratios[i]/sum;
		currentRatios[i]=classWeights_[i]/totalWeight_;
		if (currentRatios[i]>0.0)
			multFactors[i] = adjustedRatios[i] / currentRatios[i];
	}

	for (size_t i=0; i<samples_.size(); i++)
		samples_[i].weight *= multFactors[static_cast<size_t>(samples_[i].label)];

	tallyClassStatistics();
}


/*************************************************
Reduces the samples so the weight is only the
fraction of the originial
**************************************************/
void MlDataSet::randomlyReduce(double weightFraction)
{
	if (! indStatsValid_)
		tallyClassStatistics();

	vector<size_t> idxs;
	randomPermutation(samples_.size(),idxs);

	double w =0.0;
	double target = totalWeight_ * weightFraction;

	if (target>totalWeight_)
		return;

	vector<MlSample> newSamples;
	size_t i=0;
	while (w<target)
	{
		size_t idx = idxs[i++];
		newSamples.push_back(samples_[idx]);
		w+=samples_[idx].weight;
	}
	
	samples_.clear();
	samples_ = newSamples;

	tallyClassStatistics();
}


void MlDataSet::initializeOutputFile(const char* path)
{
	if (indOutStreamOpen_ && outStream_)
		closeOutputFile();

	outStream_=fopen(path,"w");
	if (! outStream_)
	{
		cout << "Error: couldn't open out stream: " << path << endl;
		exit(1);
	}

	if (buffer_)
	{
		bufferPos_=0;
		return;
	}


	bufferSize_=1048576;
	buffer_ = new char[bufferSize_];
	bufferPos_ =0;

	indOutStreamOpen_ = true;
}

void MlDataSet::writeSampleToBuffer(const MlSample& sample)
{
	const size_t numSamplePairs = sample.pairs.size();
	const size_t maxSampleSize = numSamplePairs*24 + 120; // size in chars

	if (bufferSize_ - bufferPos_ < maxSampleSize)
		flushOutputBuffer();

	char* pos = buffer_ + bufferPos_;

	if (sample.weight != 1.0)
	{
		pos+=sprintf(pos,"%d $$$WEIGHT %.3f",sample.label, sample.weight);
	} else
		pos+=sprintf(pos,"%d",sample.label);

	for (size_t i=0; i<numSamplePairs; i++)
	{
		const IdxVal& pair = sample.pairs[i];

		if (pair.value == 1.0)
		{
			pos+=sprintf(pos," F%d 1.0",static_cast<int>(pair.index));
		}
		else
			pos+=sprintf(pos," F%d %g",static_cast<int>(pair.index), pair.value);
	}
	*pos++='\n';
	bufferPos_ = (pos-buffer_);
}

void MlDataSet::closeOutputFile()
{
	if (indOutStreamOpen_ && outStream_ && bufferPos_>0)
		flushOutputBuffer();

	bufferPos_ = 0;
	if (indOutStreamOpen_ && outStream_)
		fclose(outStream_);
	outStream_ = 0;
	indOutStreamOpen_ = false;
}

void MlDataSet::flushOutputBuffer()
{
	if (! indOutStreamOpen_ || ! outStream_)
		return;

	size_t numBytesWritten = fwrite(buffer_,1,bufferPos_,outStream_);
	if (numBytesWritten != bufferPos_)
	{
		cout << "Error: could not empty data buffer to stream!" << endl;
		exit(1);
	}
	bufferPos_ = 0;
}


// performs a two pass scan:
// 1. first scan read number of lines, tokens, resolve feature names.
// 2. allocate space and read data.
void MlDataSet::readDataFile(const char *dataFile, 
						   size_t maxNumSamplesToRead,
						   bool randomSelection,
						   bool verbose)
{
	const size_t previousNumSamples = samples_.size();
	vector<size_t> featureCounts;

	ifstream stream(dataFile);
	if (! stream.is_open())
		error("couldn't open input file for reading: ",dataFile);

	char* buffer= new char[131072];
	
	// first pass
	while (! stream.eof())
	{
		stream.getline(buffer,131072);
		const size_t len = stream.gcount();
		if (len<=1 || buffer[0] == '#' || buffer[0] == '\n' || buffer[0] == '\r' )
			continue;

		// assume every feature name is prefixed with F (e.g., F0, F112, F1456, etc.)
		// count number of 'F'
		size_t n=0;
		size_t i;
		for (i=0; i<len; i++)
			if (buffer[i] == 'F')
				n++;
			
		featureCounts.push_back(n);
	}
	stream.close();

	// check if we need to be selective in choice of samples
	vector<bool> sampleFlags;
	sampleFlags.clear();
	if (featureCounts.size() > maxNumSamplesToRead)
	{
		sampleFlags.resize(featureCounts.size(),false);
		if (randomSelection)
		{
			if (verbose)
				cout << "Randomly reducing input file being read to " << 
						maxNumSamplesToRead << "/" << featureCounts.size() << endl;
			vector<size_t> idxs;
			chooseKFromN(maxNumSamplesToRead, featureCounts.size(), idxs);
			size_t i;
			for (i=0; i<idxs.size(); i++)
				sampleFlags[idxs[i]]=true;
		}
		else
		{
			if (verbose)
				cout << "Reducing input file being read to first " << 
						maxNumSamplesToRead << " samples (from a total of " << featureCounts.size() << ")" << endl;
			size_t i;
			for (i=0; i<maxNumSamplesToRead; i++)
				sampleFlags[i]=true;
		}
	}

	// second pass
	stream.clear();
	stream.open(dataFile);
	if (! stream.is_open())
		error("couldn't open input file for reading: ",dataFile);

	samples_.resize(previousNumSamples+ featureCounts.size());
	size_t samCounter=0;
	while (! stream.eof())
	{
		stream.getline(buffer,131072);
		const size_t len = stream.gcount();
		
		if (len<=1 || buffer[0] == '#' || buffer[0] == '\n' || buffer[0] == '\r' )
			continue;

		for (size_t i=0; i<len; i++)
			if (buffer[i]=='F')
				buffer[i]=' ';

		if (sampleFlags.size() > 0 && ! sampleFlags[samCounter])
		{
			samCounter++;
			continue;
		}

		const size_t n = featureCounts[samCounter];
		MlSample& sample = samples_[previousNumSamples + samCounter];
		samCounter++;

		if (n == 0)
			continue;

		numPairs_ += n;
		sample.pairs.resize(n);
		istringstream iss(buffer);

		// see if there is a number and a weight
		bool gotWeight = false;
		int label;
		double w;
		if (sscanf(buffer,"%d $$$WEIGHT %lf",&label,&w) == 2)
		{
			gotWeight = true;
		}

		iss >> sample.label;
		if (iss.fail())
		{
			cout << "BAD LINE: " << buffer << endl;
			error("error parsing input file, sample ",samCounter);
		}	
		
		if (gotWeight) // is this a weight string
		{	
			string wStr;
			iss >> wStr >> sample.weight;
			if (iss.fail())
				error("error parsing weight, sample ",samCounter);
		}
		

		for (size_t i=0; i<n; i++)
		{
			iss >> sample.pairs[i].index >> sample.pairs[i].value;
			if (iss.fail())
				error("error parsing input file, sample ",samCounter, " feature ",i);
		}

		if (sample.pairs[n-1].index > maxFeatureIndex_)
			maxFeatureIndex_ = sample.pairs[n-1].index;
	}

	delete [] buffer;

	indStatsValid_ = false;
	if (verbose)
		printDatasetStatistics();
}


void MlDataSet::convertClassLabelsForBinary(int additionalLabelToTreatAsZero)
{
	for (size_t i=0; i<samples_.size(); i++)
	{
		if (samples_[i].label>0 && samples_[i].label != additionalLabelToTreatAsZero)
		{
			samples_[i].label = 1;
		}
		else
			samples_[i].label = 0;

		if (samples_[i].label<0)
			error("negative label for sample ",i);
	}
	tallyClassStatistics();
}


void MlDataSet::tallyClassStatistics()
{
	numPairs_=0;
	numDerivedPairs_=0;
	maxFeatureIndex_=0;
	maxDerivedFeatureIndex_=0;

	vector<double> weights;
	weights.resize(32,0);
	double total=0;
	size_t i;
	for (i=0; i<samples_.size(); i++)
	{
		const MlSample& sample = samples_[i];
		if (! sample.checkConsistency())
		{
			cout << "Sample " << i << " not consistent!" << endl;
			sample.print();
			error();
		}

		if (sample.label<0)
		{
			sample.print();
			error("Sample with negative label being added! : ",sample.label);
		}

		if (sample.label >= static_cast<int>(weights.size()))
			weights.resize(static_cast<size_t>(sample.label*2),0);

		const double w = static_cast<double>(sample.weight);
		total+=w;
		weights[sample.label]+=w;
		numPairs_ += sample.pairs.size();
		numDerivedPairs_ += sample.derivedPairs.size();

		if (sample.pairs.size()>0 && sample.pairs.back().index > maxFeatureIndex_)
			maxFeatureIndex_ = sample.pairs.back().index;

		if (sample.derivedPairs.size()>0 && sample.derivedPairs.back().index > maxDerivedFeatureIndex_)
			maxDerivedFeatureIndex_ = sample.derivedPairs.back().index;
	}

	size_t maxClass=0;
	for (size_t c=0; c<weights.size(); c++)
		if (weights[c]>0)
			maxClass = c;

	numClasses_ = maxClass+1;
	totalWeight_ = static_cast<weight_t>(total);
	classWeights_.resize(numClasses_);
	for (i=0; i<numClasses_; i++)
		classWeights_[i]=static_cast<weight_t>(weights[i]);

	indStatsValid_ = true;
}


void	MlDataSet::writeWholeDataFile(const char* dataFile)
{
	initializeOutputFile(dataFile);
	
	size_t i;
	for (i=0; i<samples_.size(); i++)
		writeSampleToBuffer(samples_[i]);

	closeOutputFile();
}


void	MlDataSet::printDatasetStatistics(ostream& os)
{
	if (! indStatsValid_)
		tallyClassStatistics();

	os << "DataSet contains:" << endl;
	os << samples_.size() << "\tsamples" << endl;
	os << maxFeatureIndex_+1 << "\tfeatures (maximal index is " << maxFeatureIndex_ << ")" << endl;
	os << numPairs_ << " index-value pairs, avg of " << fixed << setprecision(3) << numPairs_ / static_cast<float>(samples_.size())
	   << " per sample." << endl;
	os << numDerivedPairs_ << " derived index-value pairs, avg of " << fixed << setprecision(3) << numDerivedPairs_ / static_cast<float>(samples_.size())
	   << " per sample." << endl;

	os << "Samples belong to " << numClasses_ << " classes:" << endl;
	size_t i;
	for (i=0; i<numClasses_; i++)
		cout << "Class " << i <<"\t weight " << classWeights_[i] << " (" << classWeights_[i]/totalWeight_ << ")" << endl;

	
	
}


void	MlDataSet::printFeatureStatistics(size_t featureIdx) const
{
	vector<value_t> vals;
	for (size_t i=0; i<samples_.size(); i++)
	{
		value_t v = samples_[i].getValue(featureIdx);
		if (v != NON_FLOAT)
			vals.push_back(v);
	}

	sort(vals.begin(),vals.end());

	cout << "N=" << vals.size() << endl;
	if (vals.size()<20)
	{
		for (size_t i=0; i<vals.size(); i++)
			cout << vals[i] << "\t";
		cout << endl;
	}
	else
	{
		cout << "First: ";
		for (size_t i=0; i<10; i++)
			cout << vals[i] << "\t";
		cout << endl;
		cout << "Last: ";
		for (size_t i=vals.size()-10; i<vals.size(); i++)
			cout << vals[i] << "\t";
		cout << endl;

		create_histogram(vals, 10, vals[0], vals[vals.size()-1]);
	}
}

void MlDataSet::outputFeatureReports(const vector<size_t>& idxs, 
									const MlFeatureSet* featureSet, 
									ostream& os) const
{
	const size_t numFeatures = featureSet->getNumBasicFeatures();
	vector< vector< vector<value_t> > > vals(2, vector< vector<value_t> >(numFeatures));
	vector< vector< double > > weights(2, vector<weight_t>(numFeatures,0.0));
	vector< vector< double > > weightsSum(2, vector<weight_t>(numFeatures,0.0));
	vector< double > totalWeights(2,0.0);
	for (size_t i=0; i<idxs.size(); i++)
	{
		const MlSample& sam = samples_[idxs[i]];
		for (size_t j=0; j<sam.pairs.size(); j++)
		{
			totalWeights[sam.label]+=sam.weight;
			if (sam.pairs[j].value != NON_FLOAT)
			{
				const size_t  index = sam.pairs[j].index;
				const value_t value = sam.pairs[j].value;
				vals[sam.label][index].push_back(value);
				weights[sam.label][index] += sam.weight;
				weightsSum[sam.label][index] += sam.weight*value;
			}
		}
	}

	os << "Report on samples in dataset:" << endl;
	os << "-----------------------------" << endl << endl;
	os << "(label weight weighted-sum #values histogram)" << endl;

	for (size_t i=0; i<numFeatures; i++)
	{
		os << i << " ]\t" << featureSet->getFeature(i).getName() << endl;
		for (int label=0; label<numClasses_; label++)
		{
			if (weights[label][i] == 0.0)
			{
				os << label << " : \t0.0" << endl;
			}
			else
			{
				os << label << " : \t" << setprecision(4) << weights[label][i]/totalWeights[label] << "\t"
					<< weightsSum[label][i] / weights[label][i] << "\t";
			
				os << vals[label][i].size();
				if (vals[label][i].size() > 10)
				{
					int ninth = (vals[label][i].size()/9 + 1);
					sort(vals[label][i].begin(),vals[label][i].end());
					for (size_t k=ninth; k<vals[label][i].size(); k+= ninth)
						os << "\t" << vals[label][i][k];
				}
				os << endl;
			}	
		}
		os << endl;
	}
}



void MlSample::print(ostream& os) const
{
	os << "LABEL  " << label << endl;
	os << "WEIGHT " << weight << endl;
	os << "INDEX-VALUE PAIRS: "<< pairs.size() << endl;
	size_t i;
	for (i=0; i<pairs.size(); i++)
		os << i << "\t" << pairs[i].index << "\t" << pairs[i].value << endl;
}



