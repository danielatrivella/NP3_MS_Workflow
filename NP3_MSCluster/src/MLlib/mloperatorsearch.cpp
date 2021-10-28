#include "mloperatorsearch.h"
#include "mlauxfun.h"

/*******************************************************************
Reads the samples and processes the data into appropriate lists
to make the processing more efficient.
Performs a two pass scan of the data.
All indexes stored in featureLists_ and elsewhere are relative to the
supplied vector idxs. That is, sample 0 in featureLists_ is idxs[0],
sample 1 is idxs[1], etc.

If run multiple times, it only adds on data from new features

********************************************************************/
void MlOperatorSearchData::readInfoFromDataSet(const MlDataSet* mld, const vector<size_t>& idxs)
{
	size_t previousNumFeatures=0;
	if (numClasses_ == 0)
	{
		numClasses_ = mld->getNumClasess() ;
		numFeatures_ = mld->getNumBasicFeatures();
		
		classWeights_.clear();
		classWeights_.resize(numClasses_,0.0);
		featureLists_.clear();
		featureLists_.reserve(numFeatures_*5);
		weightsOfFeatures_.clear();
		weightsOfFeatures_.reserve(numFeatures_*5);
		totalFeatureWeights_.clear();
		totalFeatureWeights_.reserve(numFeatures_*5);
	}
	else
		if (mld->getNumClasess()>numClasses_)
			error("new data has more classes than previous!");

	if (mld->getNumBasicFeatures()>numFeatures_)
	{
		previousNumFeatures = numFeatures_;
		numFeatures_ = mld->getNumBasicFeatures();
	}

	featureLists_.resize(numFeatures_, vector< vector<IdxVal> >(numClasses_));
	
	const vector<MlSample>& samples = mld->getSamples();
	const size_t numIdxs = idxs.size();

	if (previousNumFeatures == 0 || sampleWeights_.size() != numIdxs)
	{
		sampleWeights_.clear();
		sampleWeights_.resize(numIdxs,0.0); 
	}

	// first pass - count number of occerences in each class, feature
	vector< vector<size_t> > numOccurences(numFeatures_, vector<size_t>(numClasses_,0));
	double totalWeight =0;
	for (size_t i=0; i<numIdxs; i++)
	{
		const size_t sampleIdx = idxs[i];
		const MlSample& sample = samples[sampleIdx];

		if (sample.label<0 || sample.label>=static_cast<int>(numClasses_))
			error("Bad label in sample ",sampleIdx);

		const size_t label  = static_cast<size_t>(sample.label);
		const double weight = static_cast<double>(sample.weight);
		sampleWeights_[i]   =  sample.weight;
		totalWeight += weight;
		classWeights_[label]+=weight;
		const vector<IdxVal>& featurePairs = sample.pairs;
		for (size_t j=0; j<featurePairs.size(); j++)
			numOccurences[featurePairs[j].index][label]++;	
	}
	totalSampleWeight_ = totalWeight;

	// allocate memory
	for (size_t i=previousNumFeatures; i<numFeatures_; i++)
		for (size_t j=0; j<numClasses_; j++)
			if (numOccurences[i][j]>0)
				featureLists_[i][j].resize(numOccurences[i][j]);

	
	copyFeatureValuesFromDataSet(mld,idxs,previousNumFeatures);

	// fill in the weightsOfFeatures_
	weightsOfFeatures_.resize(numFeatures_, vector<weight_t>(numClasses_,0.0));
	totalFeatureWeights_.resize(numFeatures_, 0.0);
	for (size_t i=previousNumFeatures; i<numFeatures_; i++)
		for (size_t j=0; j<numClasses_; j++)
		{
			for (size_t k=0; k<featureLists_[i][j].size(); k++)
				weightsOfFeatures_[i][j]+= sampleWeights_[featureLists_[i][j][k].index];
			totalFeatureWeights_[i]+=weightsOfFeatures_[i][j];
		}

}

// assumes that featureLists_ has the correct sizes
// simply iterates over the  feature-index/value in each sample
// and copies the values and relative sample index into featureLists_
void MlOperatorSearchData::copyFeatureValuesFromDataSet(const MlDataSet* mld, 
														const vector<size_t>& idxs,
														size_t startFeatureIdx)
{
	const vector<MlSample>& samples = mld->getSamples();
	const size_t numIdxs = idxs.size();
	
	vector< vector<size_t> > listIdxs(numFeatures_, vector<size_t>(numClasses_,0));
	for (size_t i=0; i<numIdxs; i++)
	{
		const size_t sampleIdx = idxs[i];
		const MlSample& sample = samples[sampleIdx];

		const size_t label = static_cast<size_t>(sample.label);
		const vector<IdxVal>& featurePairs = sample.pairs;
		for (size_t j=0; j<featurePairs.size(); j++)
		{
			const size_t& featureIdx = featurePairs[j].index;
			if (featureIdx<startFeatureIdx)
				continue;

			size_t& listPointer = listIdxs[featureIdx][label];

			// the stored index is relative to idxs not the absolute sample index!
			featureLists_[featureIdx][label][listPointer++]=IdxVal(i,featurePairs[j].value); 																					
		}
	}
}


void MlOperatorSearchData::determineNumberOfValues()
{
	if (featureLists_.size()<1)
		error("Must first set featureLists!");

	numFeatureValues_.resize(numFeatures_,0);

	for (size_t featureIdx=0; featureIdx<numFeatures_; featureIdx++)
	{
		if (numFeatureValues_[featureIdx]>0)
			continue;

		vector<value_t> values;
		for (size_t c=0; c<featureLists_[featureIdx].size(); c++)
		{
			const vector<IdxVal>& listOfValues = featureLists_[featureIdx][c];
			for (size_t i=0; i<listOfValues.size(); i++)
			{
				size_t j;
				for (j=0; j<values.size(); j++)
					if (values[j] == listOfValues[i].value)
						break;

				if (j==values.size())
					values.push_back(listOfValues[i].value);

				if (values.size()>minNumValuesToBeConsideredReal)
					break;
			}
			if (values.size()>minNumValuesToBeConsideredReal)
				break;
		}

		numFeatureValues_[featureIdx] = (values.size()>minNumValuesToBeConsideredReal ?
										 MAX_UINT : values.size());
	}
}


/******************************************************************************

Processing features:
- real valued features (with more than max_values) are evaulated for unary trasnformations
  according to the bhataccharya distance. This doesn't include splits

- splits are evaluated according to the weighted information gain (to be defined below)
  after the optimal number of splits is determined, each variable might take boolean or
  transformed value, depending on the bhatacharrya distance it deplays

- Bhattacharyya distance for more than two classes is computed as the weighted pair-wise
  distance.


*******************************************************************************/



// assumes lists are sorted, does not include split or scale
void generateOperatorValues(value_t x, 
							const size_t* unaryList, 
							size_t numUnary, 
						    vector<value_t>& values)
{
	values.resize(numUnary);
	size_t pos=0;
	
	if (unaryList[pos] == UOT_SELF)
	{
		values[pos++]=x;
		if (pos == numUnary)
			return;
	}

	if (unaryList[pos] == UOT_BOOL)
	{
		values[pos++]=1.0;
		if (pos == numUnary)
			return;
	}

	if (unaryList[pos] == UOT_LOG)
	{
		if (x<=0.0)
		{
			values[pos++]=minLogValue;
		}
		else
		{
			const value_t logVal = log(x);
			values[pos++]=(logVal>minLogValue ? logVal : minLogValue);
		}
		if (pos == numUnary)
			return;
	}

	if (unaryList[pos] == UOT_EXP)
	{
		values[pos++]=exp(x);
		if (pos == numUnary)
			return;
	}

	if (unaryList[pos] == UOT_SQR)
	{
		values[pos++]=(x*x);
		if (pos == numUnary)
			return;
	}

	if (unaryList[pos] == UOT_SQRT)
	{
		values[pos++]=sqrt(x);
		if (pos == numUnary)
			return;
	}

	if (unaryList[pos] == UOT_NEG)
	{
		values[pos++]=-x;
		if (pos == numUnary)
			return;
	}

	if (unaryList[pos] == UOT_ABS)
	{
		values[pos++]=fabs(x);
		if (pos == numUnary)
			return;
	}
	error("Should not have reached this point!");
}


double BhattacharyyaDistance(double m1, double sd1, double m2, double sd2)
{
	if (m1 == m2 && sd1 == sd2)
		return 0.0;

	if (sd1 == 0.0 || sd2 == 0.0)
	{
		if (sd1 == 0.0 && sd2 == 0.0)
			return MAX_FLOAT;

		if (sd1 == 0.0)
			sd1 = sd2 / 10.0;

		if (sd2 == 0.0)
			sd2 = sd1 / 10.0;
	}

	double t1=sd1*sd1+sd2*sd2;
	double t2=m1-m2;
	return (0.25*(t2*t2)/t1 + 0.5*log(t1/(2.0*sd1*sd2)));
}



void MlOperatorSearchData::evaluateUnaryOperatorsOnFeature(size_t featureIdx, 
														 vector<double>& distances,
														 size_t featureEvaluationType) const
{
	const size_t* unaryList=0;
	size_t  numUnary=0;
	
	if (featureEvaluationType == FGT_REGRESSION)
	{
		unaryList = unaryOperatorsForRegression;
		numUnary  = numUnaryOperatorsForRegression;
	}
	else if (featureEvaluationType == FGT_RANKING)
	{
		unaryList = unaryOperatorsForRanking;
		numUnary  = numUnaryOperatorsForRanking;
	}
	else
		error("Unrecognized feature evaluation type: ",featureEvaluationType);

	// compute statistics for each class
	// TODO: for very large lists, might consider sampling
	vector< vector< MeanSdStats > > meanSdStats(numClasses_, vector<MeanSdStats>(numUnary));
	vector< vector< double > > means(numClasses_, vector<double>(numUnary));
	vector< vector< double > > sds(numClasses_, vector<double>(numUnary));

	value_t minVal = MAX_FLOAT;
	value_t maxVal = MIN_FLOAT;
	for (size_t c=0; c<numClasses_; c++)
	{
		vector<value_t> generatedValues;
		const vector<IdxVal>& sampleValuePairs = featureLists_[featureIdx][c];
		for (size_t i=0; i<sampleValuePairs.size(); i++)
		{
			const double w = static_cast<double>(sampleWeights_[sampleValuePairs[i].index]); // weight of sample
			generateOperatorValues(sampleValuePairs[i].value, unaryList, numUnary, generatedValues);
			for (size_t j=0; j<numUnary; j++)
				meanSdStats[c][j].addWX(w, generatedValues[j]);

			// check if the log operator should be invalidated because some samples have zero or negative vlaues
			if (sampleValuePairs[i].value < minVal)
				minVal = sampleValuePairs[i].value;

			if (sampleValuePairs[i].value > maxVal)
				maxVal = sampleValuePairs[i].value;
		}		
		for (size_t i=0; i<numUnary; i++)
			meanSdStats[c][i].calcMeanAndSd(means[c][i],sds[c][i]);
	}

	// compute distances
	distances.clear();
	distances.resize(UOT_NUM_OPERATORS, 0);
	if (numClasses_ == 2)
	{
		for (size_t i=0; i<numUnary; i++)
			distances[unaryList[i]]=BhattacharyyaDistance(means[0][i],sds[0][i],means[1][i],sds[1][i]);
	}
	else // compute weighted distances
	{
		weight_t totalPairWeight=0.0;
		vector<weight_t> weights(numClasses_);
		for (size_t c=0; c<numClasses_; c++)
			weights[c]=weightsOfFeatures_[featureIdx][c];

		for (size_t i=0; i<numUnary; i++)
			distances[unaryList[i]]=0.0;

		for (size_t c1=0; c1<numClasses_-1; c1++)
			for (size_t c2=c1+1; c2<numClasses_; c2++)
				for (size_t i=0; i<numUnary; i++)
				{
					const weight_t pairWeight = weights[c1]+weights[c2];
					totalPairWeight += pairWeight;
					distances[unaryList[i]]+= pairWeight * 
							BhattacharyyaDistance(means[c1][i],sds[c1][i],means[c2][i],sds[c2][i]);
				}

		for (size_t i=0; i<numUnary; i++)
			distances[unaryList[i]]/=totalPairWeight;
	}

	// invalidate certain distances based on un-desired properties of the min and max
	if (minVal<=0.0)
		distances[UOT_LOG]=0.0;
		
	if (minVal<0)
		distances[UOT_SQRT]=0.0;

	if (maxVal>3.0 && maxVal-minVal>4.0)
	{
		distances[UOT_EXP]=0.0;
	}
}




/**********************************************************************************
Creates a C X V table, where cell (c,v) holds the weight of all samples from class
c that have a value <= v . If a sample does not have the feature v, it is not considered
for the table weights. If there are too many values, then they are binned together into
|V|=~maxTableLength bins.
***********************************************************************************/
void MlOperatorSearchData::createCumulativeWeightTableForFeature(size_t featureIdx, 
											   vector< vector< weight_t > >& cumulativeWeightsTable,
											   vector< value_t >& tableColumnValues,
											   size_t maxTableLength,
											   bool verbose) const
{
	if (featureLists_.size()<=featureIdx || featureLists_[featureIdx].size()<numClasses_)
		error("trying to split a feature for which there are no lists: ",featureIdx);

	vector< vector<ValWeight> > sortedLists(numClasses_);

	// create intial sorted lists for each class <v1,w1>, <v2,w2>, ...  where v1<v2<...
	size_t totalLen=0;
	for (size_t c=0; c<numClasses_; c++)
	{
		const vector<IdxVal>& valList = featureLists_[featureIdx][c];
		vector<ValWeight>& sortedList = sortedLists[c];
		sortedList.resize(valList.size());

		for (size_t i=0; i<valList.size(); i++)
		{
			sortedList[i].value  = valList[i].value;
			sortedList[i].weight = sampleWeights_[valList[i].index];
		}
		condenseValWeightVector(sortedList,true);
		totalLen+=sortedList.size();
	}

	// merge sortedLists into one list so we can see if there is a need to condense
	// the number of values we consider
	vector<ValWeight> allVals(totalLen);
	size_t pos=0;
	for (size_t c=0; c<numClasses_; c++)
	{
		memcpy(&allVals[pos],&sortedLists[c][0],sortedLists[c].size()*sizeof(ValWeight));
		pos+=sortedLists[c].size();
	}
	condenseValWeightVector(allVals,true);

	// check if need to condense further
	if (allVals.size()>maxTableLength)
	{
		const weight_t delta = totalFeatureWeights_[featureIdx] / (maxTableLength * 1.05);

		size_t pos=0;
		while (pos<allVals.size())
		{
			weight_t w=0.0;
			size_t f=pos;
			while (w<delta && f<allVals.size())
				w+=allVals[f++].weight;

			if (f==allVals.size())
				f--;

			for (size_t i=pos; i<f; i++)
				allVals[i].value = allVals[f].value;

			pos = f;
			if (pos == allVals.size()-1)
				break;
		}
		condenseValWeightVector(allVals,false);
	}


	// copy table column values
	tableColumnValues.resize(allVals.size());
	for (size_t i=0; i<allVals.size(); i++)
		tableColumnValues[i]=allVals[i].value;
	tableColumnValues[tableColumnValues.size()-1] = MAX_FLOAT;

	// fill in table using the sorted lists
	const size_t numValues = allVals.size();
	cumulativeWeightsTable.clear();
	cumulativeWeightsTable.resize(numClasses_, vector<weight_t>(numValues,0));

	for (size_t c=0; c<numClasses_; c++)
	{
		const vector<ValWeight>& list = sortedLists[c];
		size_t pos=0;
		for (size_t v=0; v<numValues; v++)
		{
			if (v>0)
				cumulativeWeightsTable[c][v]=cumulativeWeightsTable[c][v-1];
		
			while (pos<list.size() && list[pos].value <= tableColumnValues[v])
				cumulativeWeightsTable[c][v]+=list[pos++].weight;
		}
	}

	if (verbose)
	{
		cout << endl << "Cumulative Weights Table:" << endl;
		size_t len=14;
		cout << setprecision(4) << endl;
		for (size_t i=0; i<numValues; i+=len)
		{
			const size_t maxCol = (i+len > numValues ? numValues : i+len);
			for (size_t j=i; j<maxCol; j++)
				cout << "\t" << tableColumnValues[j];
			cout << endl;
			for (size_t j=i; j<maxCol; j++)
				cout << "\t" << j;
			cout << endl;
			for (size_t c=0; c<numClasses_; c++)
			{
				cout << c;
				for (size_t j=i; j<maxCol; j++)
					cout << "\t" << (cumulativeWeightsTable[c][j]/weightsOfFeatures_[featureIdx][c]);
				cout << endl;
			}
			cout << endl << endl;
		}
	}
}


/*************************************************************************
Creates a V X V table (only cells i,j: i<=j)
where cell i,j holds the weighted entroy W(i:j)*H(C|i:j).
W(i:j) is the relative weight of the samples in bins i..j and H(C|i:j) is 
the entropy of the class variable in those bins
**************************************************************************/
void MlOperatorSearchData::createEntropyTableForFeature(
								const vector< vector<weight_t> >& cumulativeWeightsTable,
								vector< vector<weight_t> >& entropyTable,
								bool verbose) const
{
	assert( cumulativeWeightsTable.size() == numClasses_ );
	const size_t tableSize = cumulativeWeightsTable[0].size();
	const weight_t oneOverLog2 = 1.0 / log(2.0);

	weight_t w=0.0;
	for (size_t c=0; c<numClasses_; c++)
		w+=cumulativeWeightsTable[c][tableSize-1];
	const weight_t oneOverTotalWeight = 1.0 / w;

	// resize the entropy table if needed
	if (tableSize>entropyTable.size())
	{
		const size_t oldSize = entropyTable.size();
		entropyTable.resize(tableSize);
		for (size_t i= oldSize; i<tableSize; i++)
			entropyTable[i].resize(i+1);
	}

	vector<weight_t> probOfClass(numClasses_);
	for (size_t i=0; i<tableSize; i++)
		for (size_t j=0; j<=i; j++)
		{
			vector<weight_t> weightPerClass(numClasses_);
			weight_t rangeWeight=0.0;
			for (size_t c=0; c<numClasses_; c++)
			{
				weight_t w=cumulativeWeightsTable[c][i];
				if (j>0)
					w-=cumulativeWeightsTable[c][j-1];
				probOfClass[c]=w;
				rangeWeight+=w;
			}

			weight_t entropy = 0.0;
			for (size_t c=0; c<numClasses_; c++)
			{
				probOfClass[c] /= rangeWeight;
				if (probOfClass[c]>0.0)
					entropy -= (probOfClass[c] * log (probOfClass[c]));
			}
			entropy *= oneOverLog2;
			entropyTable[i][j] = (rangeWeight * oneOverTotalWeight) * entropy; // weighted entropy
		}

	if (verbose)
	{
		cout << endl << "Entropy Table: " << endl;
		for (size_t i=0; i<tableSize; i++)
			cout << "\t" << i;
		cout << endl;
		cout << setprecision(4);
		for (size_t i=0; i<tableSize; i++)
		{
			cout << i;
			for (size_t j=0; j<=i; j++)
				cout << "\t" << entropyTable[i][j];
			cout << endl;
		}
	}

}



struct EntropyDpCell {
	EntropyDpCell() : weightedEntropy(MAX_FLOAT), leftMostSplitIdx(MAX_UINT) {}

	weight_t weightedEntropy;
	size_t   leftMostSplitIdx;
};

/*************************************************************************
Considers all ways to split a node in up to *maxNumSplit* bins so that
the information gain is maximized. IG(Split) = n^reg * [H(C) - sum_i=1..n W(i)H(C|i)]
where reg is the regularization of the gain (-0.5), H(C) is the entropy of the
class variable without the split. W(i) is the weight of the samples in the i'th bin
and H(C|i) is the entropy of the class label of the samples in the i'th bin

The function constructs a dynamic programming table for detrmining the weighted entropy in
a bin between values i and j (if there are too many values then they are initially
condensed into a smaller number of bins e.g., 100 and the table is constructed for
that value). Dynamic programming is then used to find the optimal split with k bins
based on the optimal split with k-1 bins.

If the number of values is V, finding the best split into at most k bins is
done in time and space that is O(kV+V^2)
**************************************************************************/
double MlOperatorSearchData::findOptimalSplit(
						  size_t featureIdx,
						  vector<value_t>& splitThresholds, 
						  double gainRegularization,
						  double minGainForSplit,
						  size_t maxNumSplits,
						  size_t maxTableLength,
						  bool   verbose) const
{
	splitThresholds.clear();

	if (maxNumSplits<1)
		return 0.0;

	vector< vector<weight_t> > cumulativeWeightsTable;
	vector< value_t >		   tableColumnValues;
	vector< vector<weight_t> > weightedEntropy;
	
	createCumulativeWeightTableForFeature(featureIdx, cumulativeWeightsTable, 
										  tableColumnValues, maxTableLength);

	createEntropyTableForFeature(cumulativeWeightsTable, weightedEntropy);

	const size_t numValues = tableColumnValues.size();
	if (maxNumSplits>= numValues)
		maxNumSplits = numValues-1;

	// create dp table: cell k,i has the lowest weighted entropy using k splits with left most 
	// split giving i to the left side
	vector< vector<EntropyDpCell> > dpTable(maxNumSplits+1, vector<EntropyDpCell>(numValues));

	// table filling starts at k=0, there is only one option for a split in each cell
	for (size_t i=0; i<numValues; i++)
	{
		dpTable[0][i].leftMostSplitIdx = MAX_UINT;	// there is no internal split
		dpTable[0][i].weightedEntropy = weightedEntropy[numValues-1][i];
	}

	for (size_t k=1; k<=maxNumSplits; k++)
	{
		const size_t lastLeftIdx = numValues-k;
		for (size_t i=0; i<lastLeftIdx; i++)
		{
			// init with the split that gives the bin i to the left split and all the rest
			// is done using k-1 splits for bins i+1..V
			EntropyDpCell& cell = dpTable[k][i];
			cell.leftMostSplitIdx = i;
			cell.weightedEntropy = weightedEntropy[i][i] + dpTable[k-1][i+1].weightedEntropy;
			
			for (size_t j=i+1; j<lastLeftIdx; j++)
			{
				const weight_t w = weightedEntropy[j][i] + dpTable[k-1][j+1].weightedEntropy;
				if (w<cell.weightedEntropy)
				{
					cell.weightedEntropy = w;
					cell.leftMostSplitIdx = j;
				}
			}
		}
	}

	if (0)
	{
		cout << endl << "DP-Table:" << endl << fixed << setprecision(4);
		for (size_t i=0; i<dpTable[0].size(); i++)
			cout << "\t" << i;
		cout << endl << endl;

		for (size_t k=0; k<dpTable.size(); k++)
		{
			cout << k;
			for (size_t i=0; i<dpTable[0].size(); i++)
				if (dpTable[k][i].weightedEntropy<MAX_FLOAT)
				{
					cout << "\t" << dpTable[k][i].weightedEntropy;
				}
				else
					cout << "\t  -";

			cout << endl;

			for (size_t i=0; i<dpTable[0].size(); i++)
				if (dpTable[k][i].leftMostSplitIdx<MAX_UINT)
				{
					cout << "\t" << dpTable[k][i].leftMostSplitIdx;
				}
				else
					cout << "\t  -";
			cout << endl << endl;
		}
	}

	// find the best splits for k=0..maxNumSplits
	if (verbose)
	{
		cout << endl << "Split results:" <<endl;
		cout << "#splits\tH(C)\t\tIG\t\tnorm-IG\t\tthreshes" << endl;
	}

	vector< vector<size_t> > bestSplits(maxNumSplits+1);
	vector<weight_t>		 bestEntropies(maxNumSplits+1, MAX_FLOAT);
	vector<weight_t>		 informationGains(maxNumSplits+1, 0.0);
	vector<weight_t>		 regularizedInformationGains(maxNumSplits+1, 0.0);
	for (size_t i=0; i<=maxNumSplits; i++)
	{
		bestEntropies[i]=dpTable[i][0].weightedEntropy;
		size_t leftIdx=0;
		for (size_t j=i; j>0; j--)
		{
			size_t newLeftIdx = dpTable[j][leftIdx].leftMostSplitIdx;
			if (newLeftIdx<MAX_UINT)
			{
				bestSplits[i].push_back(newLeftIdx);
				leftIdx = newLeftIdx+1;
			}
		}
		informationGains[i]=bestEntropies[0]-bestEntropies[i];
		regularizedInformationGains[i] = informationGains[i] * pow(i+1.0,gainRegularization);
		if (verbose)
		{
			cout << fixed << i << setprecision(6) << "\t" << bestEntropies[i] << "\t" << informationGains[i]
				<< "\t" << regularizedInformationGains[i] << "\t[" << setprecision(4);
			for (size_t j=0; j<bestSplits[i].size(); j++)
				cout << " " << bestSplits[i][j] << ": " << tableColumnValues[bestSplits[i][j]];
			cout << " ]" << endl;
		}
	}


	size_t splitIndex=0;
	for (size_t i=1; i<=maxNumSplits; i++)
		if (regularizedInformationGains[i]>regularizedInformationGains[splitIndex])
			splitIndex = i;

	double relativeGain = 1.0 - bestEntropies[splitIndex]/bestEntropies[0];
	if (verbose)
		cout << "Relative gain: " << relativeGain << endl;

	// select best split accordind to regularization
	if (relativeGain>= minGainForSplit)
	{
		splitThresholds.resize(bestSplits[splitIndex].size());
		for (size_t i=0; i<bestSplits[splitIndex].size(); i++)
			splitThresholds[i]=tableColumnValues[bestSplits[splitIndex][i]];
	}
	return relativeGain;
}


//
void MlOperatorSearchData::printDetailedSummary(MlDataSet* ds,
												const MlFeatureSet& fs,
												ostream& os)
{
	const vector<MlFeature>& features = fs.getFeatures();
	ds->printDatasetStatistics(os);

	if (features.size() != this->featureLists_.size())
		error("feature dimensions do not match between search data and feature set!");



	// print columns
	os << endl;
	os << "Index\tWeight\t#vals\tIG\tBhatt";
	for (size_t c=0; c<numClasses_; c++)
		os << "\tC" << c <<" W\tC" << c <<" Mu\tC" << c <<" sd";
	os << "\tFeature name" << endl << setprecision(4);

	for (size_t featureIdx=0; featureIdx<features.size(); featureIdx++)
	{
		cout << featureIdx << "\t" << totalFeatureWeights_[featureIdx]/totalSampleWeight_ << "\t";
		if (numFeatureValues_[featureIdx] < MAX_UINT)
		{
			cout << numFeatureValues_[featureIdx] << "\t";
		}
		else
			cout << "Inf\t";

		// information gain of the class conditioned on having the variable 
		if (totalFeatureWeights_[featureIdx]>0.0 &&
			totalFeatureWeights_[featureIdx]<totalSampleWeight_)
		{
			double H=0.0;
			double Hcond = 0.0;

			// compute unconditional entropy of class variable (still split into two cases
			// based on presence or missing of the feature
			for (size_t c=0; c<numClasses_; c++)
			{
				weight_t relativeClassWeight   = classWeights_[c]/totalSampleWeight_;
				weight_t relativeFeatureWeight = totalFeatureWeights_[featureIdx]/totalSampleWeight_;
				
				if (relativeClassWeight>0.0)
				{
					if (relativeFeatureWeight>0.0)
					{
						weight_t p=relativeClassWeight*relativeFeatureWeight;
						H-= p * log(p);
					}
					if (relativeFeatureWeight<1.0)
					{
						weight_t p=relativeClassWeight*(1.0-relativeFeatureWeight);
						H-= p * log(p);
					}
				}
			}
			H/= log(2.0);
			for (size_t c=0; c<numClasses_; c++)
			{
				weight_t pFeature = (weightsOfFeatures_[featureIdx][c]/totalSampleWeight_);
				if (pFeature>0.0)
					Hcond -= pFeature * log(pFeature);
			
				weight_t pWithoutFeature = (classWeights_[c] - weightsOfFeatures_[featureIdx][c])/totalSampleWeight_;
				if (pWithoutFeature>0.0)
					Hcond -= pWithoutFeature * log(pWithoutFeature);
			}
			Hcond /= log(2.0);
			
			cout << H-Hcond << "\t";
		}
		else
			cout << "0.0\t";

		// calc mean and sd of classes
		vector<MeanSdStats> meanSdStats(numClasses_);
		vector< double > means(numClasses_, 0.0);
		vector< double > sds(numClasses_, 0.0);
		for (size_t c=0; c<numClasses_; c++)
		{
			const vector<IdxVal>& sampleValuePairs = featureLists_[featureIdx][c];
			for (size_t i=0; i<sampleValuePairs.size(); i++)
			{
				const double w = static_cast<double>(sampleWeights_[sampleValuePairs[i].index]); // weight of sample
				meanSdStats[c].addWX(w, sampleValuePairs[i].value);
			}		
			for (size_t i=0; i<numClasses_; i++)
				meanSdStats[c].calcMeanAndSd(means[c],sds[c]);
		}

		// Bhattacharyya distance
		if (numFeatureValues_[featureIdx] == MAX_UINT)
		{
			double bhatt=0.0;
			if (numClasses_ == 2)
			{
				bhatt=BhattacharyyaDistance(means[0],sds[0],means[1],sds[1]);
			}
			else // compute weighted distances
			{
				weight_t totalPairWeight=0.0;
				for (size_t c1=0; c1<numClasses_-1; c1++)
					for (size_t c2=c1+1; c2<numClasses_; c2++)
					{
						const weight_t pairWeight = weightsOfFeatures_[featureIdx][c1]+weightsOfFeatures_[featureIdx][c2];
						totalPairWeight += pairWeight;
						bhatt+= pairWeight * 
									BhattacharyyaDistance(means[c1],sds[c1],means[c2],sds[c2]);
					}

				bhatt/=totalPairWeight;
			}
			cout << bhatt;
		}
		else
			cout << "  -";

		for (size_t c=0; c<numClasses_; c++)
		{
			cout << "\t" << weightsOfFeatures_[featureIdx][c]/totalFeatureWeights_[c];
			cout << "\t" << means[c];
			if (numFeatureValues_[featureIdx] == MAX_UINT)
			{
				cout << "\t" << sds[c];
			}
			else
				cout << "\t  -";
		}
		cout << "\t" << features[featureIdx].getName() << endl;
	}


}

