#include "mloperatorlist.h"
#include "mltrainingcontainer.h"

size_t  MlOperatorList::determineFeatureGenerationType(size_t algorithmType)
{
	if (algorithmType == SMT_LOGISTIC_CG ||
		algorithmType == SMT_LOGISTIC_LMBFGS ||
		algorithmType == SMT_MAXIMUM_ENTROPY_LMBFGS ||
		algorithmType == SMT_MAXIMUM_ENTROPY_CG_FR ||
		algorithmType == SMT_MAXIMUM_ENTROPY_CG_PR)
		return FGT_REGRESSION;

	if (algorithmType == SMT_RANKBOOST_ADA ||
		algorithmType == SMT_RANKBOOST_LOGIT ||
		algorithmType == SMT_RANKBOOST_SMOOTH)
		return FGT_RANKING;

	error("Unknown model training type:",algorithmType);
	return MAX_UINT;
}



bool MlOperatorList::readOperatorList(const char* path)
{
	ifstream ifs(path);
	if (! ifs.good() || ! readOperatorList(ifs))
		return false;

	ifs.close();
	return true;
}

bool MlOperatorList::readOperatorList(ifstream& ifs)
{
	char buffer[1024];
	while (ifs.good() && ifs.getline(buffer,1024))
		if (ifs.gcount()>0 && buffer[0] != '#')
			break;

	size_t n=0;
	if (sscanf(buffer,"%d",&n) != 1)
		error("expected line with number of operaotrs");

	executionOrder_.clear();
	executionOrder_.reserve(n);
	drops_.clear();
	indicators_.clear();
	normalizations_.clear();
	functions_.clear();
	splits_.clear();
	conditionals_.clear();

	for (size_t i=0; i<n; i++)
	{
		if (! ifs.good())
			return false;

		ifs.getline(buffer,1024);
		if (ifs.gcount()<=0)
			return false;

		istringstream iss(buffer);

		char type=' ';
		iss >> type;
		switch (type) {
			case 'D':
				{
					size_t dropIdx=0;
					iss >> dropIdx;
					if (iss.fail())
						error("Error reading line:",buffer);

					executionOrder_.push_back(IdxPair(OT_DROP,drops_.size()));
					drops_.push_back(dropIdx);
				}
				break;

			case 'I':
				{
					IndicatorOperator iop;
					iss >> iop.sourceIdx >> iop.targetIdx;
					if (iss.fail())
						error("Error reading line:",buffer);
					executionOrder_.push_back(IdxPair(OT_INDICATOR,indicators_.size()));
					indicators_.push_back(iop);
					break;
				}

			case 'F':
				{
					FunctionOperator fop;
					iss >> fop.sourceIdx >> fop.targetIdx;
					string typeStr;
					iss >> typeStr;
					size_t i;
					for (i=0; i<numConditionalValueLabels; i++)
						if (! strcmp(typeStr.c_str(),conditionalValueLabels[i]))
							break;
					if (i == numConditionalValueLabels)
						error("Error reading line:",buffer);
					fop.type = i;
					executionOrder_.push_back(IdxPair(OT_FUNCTION,functions_.size()));
					functions_.push_back(fop);
				}

			case 'N':
				{
					NormalizationOperator nop;
					iss >> nop.sourceIdx >> nop.targetIdx >> nop.mu >> nop.sigma;
					if (iss.fail())
						error("Error reading line:",buffer);
					executionOrder_.push_back(IdxPair(OT_NORMALIZATION,normalizations_.size()));
					normalizations_.push_back(nop);
				}
				break;

			case 'S':
				{
					SplitOperator sop;
					iss >> sop.sourceIdx;

					size_t n;
					iss >> n;
					if (iss.fail())
						error("Error reading line:",buffer);
					sop.thresholds.resize(n);
					sop.indexesForBinValues.resize(n+1);
					sop.indexesForBinIndicators.resize(n+1);
					for (size_t i=0; i<n; i++)
						iss >> sop.thresholds[i];
					for (size_t i=0; i<=n; i++)
					{
						// a fix so the model file can use -1 for MAX_UINT (so it won't hurt my eyes...)
						int index;
						iss >> index;
						sop.indexesForBinValues[i] = (index>=0 ? static_cast<size_t>(index) : MAX_UINT);
					}
					for (size_t i=0; i<=n; i++)
					{
						// a fix so the model file can use -1 for MAX_UINT (so it won't hurt my eyes...)
						int index;
						iss >> index;
						sop.indexesForBinIndicators[i] = (index>=0 ? static_cast<size_t>(index) : MAX_UINT);
					}
					if (iss.fail())
						error("Error reading line:",buffer);
					executionOrder_.push_back(IdxPair(OT_SPLIT,splits_.size()));
					splits_.push_back(sop);
				}
				break;

			case 'C':
				{
					ConditionalOperator cod;
					size_t n;
					iss >> n;
					cod.sourceIdxs.resize(n);
					for (size_t i=0; i<n; i++)
						iss >> cod.sourceIdxs[i];
					iss >> cod.targetIdx;

					// a fix so the model file can use -1 for MAX_UINT (so it won't hurt my eyes...)
					int indexForBool;
					iss >> indexForBool;
					cod.indexForBool = (indexForBool>=0 ? static_cast<size_t>(indexForBool) : MAX_UINT);
					
					string resultStr, conditionStr;
					iss >> conditionStr >> resultStr;
					if (iss.fail())
						error("Error reading line:",buffer);
					size_t idxType, idxResult;
					for (idxType=0; idxType<numConditionalOperatorLabels; idxType++)
						if (! strcmp(conditionStr.c_str(),conditionalOperatorLabels[idxType]))
							break;
					if (idxType == numConditionalOperatorLabels)
						error("Error reading line:",buffer);
					for (idxResult=0; idxResult<numConditionalValueLabels; idxResult++)
						if (! strcmp(resultStr.c_str(),conditionalValueLabels[idxResult]))
							break;
					if (idxResult == numConditionalValueLabels)
						error("Error reading line:",buffer);
					cod.conditionType = idxType;
					cod.resultType    = idxResult;

					executionOrder_.push_back(IdxPair(OT_CONDITIONAL,conditionals_.size()));
					conditionals_.push_back(cod);
				}
				break;
		};
	}
	return true;
}

bool MlOperatorList::writeOperatorList(const char* path) const
{
	ofstream ofs(path);
	if (! ofs.good() || ! writeOperatorList(ofs))
		error("couldn't write operator list: ",path);
	ofs.close();
	return true;
}

/*bool MlOperatorList::writeOperatorList(ostream& os) const
{
	char buffer[256];
	os << operators_.size() << endl;
	for (size_t i=0; i<operators_.size(); i++)
	{
		operators_[i].writeOperator(buffer);
		os << buffer << endl;
	}
	return true;
}*/

bool MlOperatorList::writeOperatorList(ostream& os) const
{
	ostringstream sstream;
	sstream << executionOrder_.size() << endl;
	for (size_t i=0; i<executionOrder_.size(); i++)
	{
		const size_t operatorType = executionOrder_[i].idx1;
		const size_t pos = executionOrder_[i].idx2;
		sstream << fixed;
		switch (operatorType)
		{
			case OT_DROP: 
				sstream << "D\t" << drops_[pos];
				break;

			case OT_INDICATOR:
				sstream << "I\t" << indicators_[pos].sourceIdx << "\t" << indicators_[i].targetIdx;
				break;
				
			case OT_NORMALIZATION:
				sstream << "N\t" << normalizations_[pos].sourceIdx << "\t" << normalizations_[pos].targetIdx
					<< scientific << normalizations_[pos].mu << "\t" << normalizations_[pos].sigma;
				break;

			case OT_FUNCTION:
				sstream << "F\t" << functions_[pos].sourceIdx << "\t" << normalizations_[pos].targetIdx
					<< unaryOperatorLabels[functions_[pos].type];
				break;

			case OT_SPLIT:
				sstream << "S\t" << splits_[pos].sourceIdx << "\t" << splits_[pos].thresholds.size() << scientific;
				for (size_t i=0; i<splits_[pos].thresholds.size(); i++)
					sstream << "\t" << splits_[pos].thresholds[i];

				// a fix so the model file can use -1 for MAX_UINT (so it won't hurt my eyes...)
				for (size_t i=0; i<splits_[pos].indexesForBinValues.size(); i++)
					if (splits_[pos].indexesForBinValues[i]<MAX_UINT)
					{
						sstream << "\t" << splits_[pos].indexesForBinValues[i];
					}
					else
						sstream << "\t-1";

				// a fix so the model file can use -1 for MAX_UINT (so it won't hurt my eyes...)
				for (size_t i=0; i<splits_[pos].indexesForBinIndicators.size(); i++)
					if (splits_[pos].indexesForBinIndicators[i]<MAX_UINT)
					{
						sstream << "\t" << splits_[pos].indexesForBinIndicators[i];
					}
					else
						sstream << "\t-1";
				break;
				
			case OT_CONDITIONAL:
				sstream << "C\t" << conditionals_[pos].sourceIdxs.size();
				for (size_t i=0; i<conditionals_[pos].sourceIdxs.size(); i++)
					sstream << "\t" << conditionals_[pos].sourceIdxs[i];
				sstream << conditionals_[pos].targetIdx << "\t";

				// a fix so the model file can use -1 for MAX_UINT (so it won't hurt my eyes...)
				if (conditionals_[pos].indexForBool<MAX_UINT)
				{
					sstream << "\t" << conditionals_[pos].indexForBool;
				}
				else
					sstream << "\t-1";

				sstream << conditionalOperatorLabels[conditionals_[pos].conditionType] << "\t" 
						<< conditionalValueLabels[conditionals_[pos].resultType];

			default:
				error("Unrecognized operator type: ",operatorType);
		};
		os << sstream.str() << endl;
	}
	return true;
}



/*************************************************************************
This function sequentially applies the operators to the feature values in
the sample. Intermidate results are stored in an expanded array (as opposed
to the sparse feature representation in the sample).

This function takes note of the previous operators that were applied to the
sample, and only applies operators that are new (larger index)

An invariant that is maintained:
- Every entry in the expanded array that is != NON_FLOAT has a corresponding
  feature-value pair in the sample (so if the operator sees such a value in
  the array, it does not add a pair to the sample).
**************************************************************************/
void MlOperatorList::applyOperatorList(MlSample& sam) const
{
	if (sam.pairs.size() == 0)
		return;

	// enlarge expanded vector
	if (expandedVector_.size() < operatorFeatureSpaceSize_)
	{
		expandedVector_.clear();
		expandedVector_.resize(operatorFeatureSpaceSize_, NON_FLOAT); // NON_FLOAT represents NO VALUE
	}

	// copy original feature values
	for (size_t i=0; i<sam.pairs.size(); i++)
		expandedVector_[sam.pairs[i].index]=sam.pairs[i].value;

	// apply operators
	for (size_t i=sam.nextOperatorToApply; i<executionOrder_.size(); i++)
	{
		const size_t operatorType = executionOrder_[i].idx1;
		const size_t position	  = executionOrder_[i].idx2;

		switch (operatorType) {

			case OT_DROP:
				expandedVector_[drops_[position]]=NON_FLOAT;
				break;

			case OT_INDICATOR:
				if (expandedVector_[indicators_[position].sourceIdx] != NON_FLOAT)
				{
					const size_t targetIndex = indicators_[position].targetIdx;
					if (expandedVector_[targetIndex] == NON_FLOAT)
						sam.addPair(targetIndex,1.0);
					expandedVector_[targetIndex] = 1.0;
				}
				break;

			case OT_NORMALIZATION:
				{
					const value_t v = expandedVector_[normalizations_[position].sourceIdx];
					if (v != NON_FLOAT)
					{
						const value_t normVal = (v - normalizations_[position].mu) / normalizations_[position].sigma;
						const size_t  targetIndex = normalizations_[position].targetIdx;
						if (expandedVector_[targetIndex])
							sam.addPair(targetIndex,v);
						expandedVector_[targetIndex] = v;
					}
				}
				break;

			case OT_SPLIT:
				{
					const SplitOperator& sop = splits_[position];
					const value_t v = expandedVector_[sop.sourceIdx];
					if (v != NON_FLOAT)
					{
						size_t i;
						for (i=0; i<sop.thresholds.size(); i++)
							if (v<= sop.thresholds[i])
								break;

						if (sop.indexesForBinValues[i]>=0)
						{
							const size_t index = static_cast<size_t>(sop.indexesForBinValues[i]);
							if (expandedVector_[index] == NON_FLOAT)
								sam.addPair(index,v);
							expandedVector_[index]=v;
						}
						if (sop.indexesForBinIndicators[i]>=0)
						{
							const size_t index = static_cast<size_t>(sop.indexesForBinIndicators[i]);
							if (expandedVector_[index] == NON_FLOAT)
								sam.addPair(index,1.0);
							expandedVector_[index]=1.0;
						}
					}
				}
				break;
			
			case OT_FUNCTION:
				{
					const value_t v = expandedVector_[functions_[position].sourceIdx];
					if (v != NON_FLOAT)
					{
						const size_t functionType = functions_[position].type;
						// apply functions
						value_t result = NON_FLOAT;
						switch (functionType)
						{
						case UOT_SELF: result = v; 
							break;

						case UOT_BOOL: result = 1.0;
							break;
								
						case UOT_LOG: result = (v > 0.0? log(v) : NON_FLOAT);
							break;

						case UOT_EXP: result = exp(v);
							break;
								
						case UOT_SQR: result = v*v;
							break;

						case UOT_SQRT: result = sqrt(v);
							break;

						case UOT_NEG: result = -v;
							break;

						case UOT_ABS: result = abs(v);
							break;
							
						default:
							error("Unrecognized function type: ",functionType);
						};

						const size_t targetIndex = functions_[position].targetIdx;
						if (expandedVector_[targetIndex]==NON_FLOAT)
							sam.addPair(targetIndex,result);
						expandedVector_[targetIndex] = result;
					}
				}
				break;
			
			case OT_CONDITIONAL:
				{
					const ConditionalOperator& cond = conditionals_[position];

					// evaluate condition
					bool  evaluation = false;
					if (cond.conditionType == COT_AND)
					{
						size_t i;
						for (i=0; i<cond.sourceIdxs.size(); i++)
							if (expandedVector_[cond.sourceIdxs[i]] == NON_FLOAT)
								break;
						if (i==cond.sourceIdxs.size())
							evaluation=true;
					}
					else if (cond.conditionType == COT_OR)
					{
						size_t i;
						for (i=0; i<cond.sourceIdxs.size(); i++)
							if (expandedVector_[cond.sourceIdxs[i]] != NON_FLOAT)
								break;
						if (i<cond.sourceIdxs.size())
							evaluation = true;
					}

					if (evaluation) // create result value
					{	
						value_t result=NON_FLOAT;
						if (cond.resultType == CVT_BOOL)
						{
							result =1.0;
						}
						else if (cond.resultType == CVT_FEATURE)
						{
							result = expandedVector_[cond.sourceIdxs[0]]; // index 0 is always the one
						}
						else
						{
							// collect values
							static vector<value_t> values;
							const size_t maxNumValues = cond.sourceIdxs.size();
							if (values.size()< maxNumValues)
								values.resize(maxNumValues);

							size_t numValues=0;
							for (size_t i=0; i<maxNumValues; i++)
							{
								values[numValues]=expandedVector_[cond.sourceIdxs[i]];
								if (values[numValues] != NON_FLOAT)
									numValues++;
							}

							switch (cond.resultType)
							{
								case CVT_SUM:
									result=0.0;
									for (size_t i=0; i<numValues; i++)
										result+=values[i];
									break;

								case CVT_PROD:
									result=1.0;
									for (size_t i=0; i<numValues; i++)
										result*=values[i];
									break;

								case CVT_MAX:
									result = MIN_FLOAT;
									for (size_t i=0; i<numValues; i++)
										if (values[i]>result)
											result=values[i];
									break;
								
								case CVT_MIN:
									result = MAX_FLOAT;
									for (size_t i=0; i<numValues; i++)
										if (values[i]<result)
											result=values[i];
									break;

								case CVT_AVG:
									result = 0.0;
									for (size_t i=0; i<numValues; i++)
										result+=values[i];
									result /= static_cast<value_t>(numValues);
									break;

								case CVT_COUNT:
									result = static_cast<value_t>(numValues);
									break;

								case CVT_MEDIAN:
									{
									vector<value_t> copyOfValues(values.size());
									memcpy(&copyOfValues[0],&values[0],numValues*sizeof(value_t));
									sort(copyOfValues.begin(),copyOfValues.end());
									result = copyOfValues[numValues/2];
									}
									break;

								default:
									error("Unrecognized result type: ",cond.resultType);
							};
						}
						if (expandedVector_[cond.targetIdx] == NON_FLOAT)
							sam.addPair(cond.targetIdx,result);
						expandedVector_[cond.targetIdx]=result;
					}
				}
				break;

			default: error("Unrecognized operator type: ",operatorType);
		};
	}

	// copy updated values and revert cells to NON_FLOAT
	for (size_t i=0; i<sam.pairs.size(); i++)
	{
		const size_t& index = sam.pairs[i].index;
		sam.pairs[i].value = expandedVector_[index];
		expandedVector_[index] = NON_FLOAT;
	}

	sam.nextOperatorToApply = executionOrder_.size();

}


/**************************************************************************
Applies the operators in the list to all the samples
***************************************************************************/
void MlOperatorList::applyOperatorList(MlDataSet* dataSet, const vector<size_t>& sampleIdxs) const
{
	vector<MlSample>& samples = dataSet->getSamples();
	for (size_t i=0; i<sampleIdxs.size(); i++)
		applyOperatorList(samples[sampleIdxs[i]]);
}


/********************************************************************************
// chooses operators that can be useful in separating classes
// this is done before an initial model is created
// does not include splits
// this function is useful for regressions, not for ranking
*********************************************************************************/
void MlOperatorList::selectInitialUnaryOperators(MlTrainingContainer& mtc,
												 const MlFeatureSet& featureSet, 
												 bool verbose)
{
/*	MlOperatorSearchData* auxilaryData = mtc.getAuxilarySearchData();
	const size_t featureGenerationType = mtc.getFeatureGenerationType();
	const double minBhattacharyyaAdvantage = mtc.getBhattacharyyaDistance();
	const vector<size_t>& numFeatureValues = auxilaryData->numFeatureValues_;
	
	const vector<MlFeature>& features = featureSet.getFeatures();

	assert( auxilaryData->getNumBasicFeatures() == features.size());
	if (operatorFeatureSpaceSize_ < features.size() )
		operatorFeatureSpaceSize_ = features.size();

	operators_.clear();
	if (verbose)
	{
		cout << "Feature";
		for (size_t i=0; i<numUnaryOperatorLabels; i++)
			cout << "\t" << unaryOperatorLabels[i];
		cout << endl;
	}

	for (size_t i=0; i<features.size(); i++)
	{
		if (! features[i].getFlagConsiderUnary())
			continue;

		if (numFeatureValues[i]<MAX_UINT)
			continue;

		size_t bestIdx   = UOT_SELF;

		vector<double> distances;
		auxilaryData->evaluateUnaryOperatorsOnFeature(i, distances, featureGenerationType);
		for (size_t j=0; j<distances.size(); j++)
			if (distances[j]>distances[bestIdx])
				bestIdx=j;

		if (verbose)
		{
			cout << i;
			for (size_t j=0; j<distances.size(); j++)
				cout << fixed << setprecision(3) << "\t" << distances[j];
		}

		// find feature's advantage
		double advantage=0.0;
		if (distances[UOT_SELF]<=0.0 && distances[bestIdx]>0.0)
		{
			advantage = MAX_FLOAT;
		}
		else
			advantage = distances[bestIdx] - distances[UOT_SELF];

		if (advantage>=minBhattacharyyaAdvantage)
		{
			MlOperator newOperator;
			newOperator.indUnaryOperator_ = true;
			newOperator.indAddIndicator_ = false;
		    newOperator.type_ = bestIdx;
			newOperator.targetFeatureIndex_ = i;
			operators_.push_back(newOperator);

			// update feature name
			string newName = features[i].getName() + "[" + static_cast<string>(unaryOperatorLabels[bestIdx]) + "]";
		//	features[i].setName(newName);

			if (verbose)
				cout << "\t(" << setprecision(3) << distances[bestIdx]-distances[UOT_SELF] << ")";
		}
		else
		{
			if (verbose)
				cout << "\t   -";
		}
		if (verbose)
			cout << "\t" << features[i].getName() << endl;
	}

	if (operators_.size()>0)
		indWasInitialized_ = true;*/

	// apply the operators to the data, and update the search data
/*	
*/

	
/*	for (size_t i=0; i<features.size(); i++)
	{
		cout << endl << "FEATURE " << i <<"  " << features[i].getName() << endl;
		mtc.getTrainingDataSet()->printFeatureStatistics(i);
	}
	exit(0);*/
	//mtc.getTestDataSet()->printFeatureStatistics(65);
	// update search data
	

//	vector<value_t> threshes;
	//auxilaryData->findOptimalSplit(54, threshes, -0.333 , 0.1, 15, 500);

/*	vector< vector<weight_t> > ctw;
	vector<value_t> tcv;
	vector<vector<float>> entropy;
	auxilaryData->createCumulativeWeightTableForFeature(1,ctw,tcv,10);
	auxilaryData->createEntropyTableForFeature(1,ctw,entropy);*/
}



// chooses split operators that can be usesful in separating the classes
// this is done before an intial model is created
// This function does not create the appropriate features, and should be
// called on a temporary list.
// verbose = 0, minimal chatter
// verbose = 1, one line summary for each feature
// verbose > 1, more details
void MlOperatorList::selectInitialSplitOperators(MlTrainingContainer& mtc,
												 const MlFeatureSet& featureSet,
												 bool  addBooleanIndicators,
												 size_t verboseLevel)
{
	MlOperatorSearchData* auxilaryData = mtc.getAuxilarySearchData();
	const size_t featureGenerationType = mtc.getFeatureGenerationType();
	const double minRelativeInformationGain = mtc.getMinRelativeInformationGain();
	const size_t maxNumSplits				= mtc.getMaxNumSplits();
	const vector<size_t>& numFeatureValues  = auxilaryData->numFeatureValues_;
	const vector<MlFeature>& features       = featureSet.getFeatures();

	assert( auxilaryData->getNumBasicFeatures() == features.size());
	if (operatorFeatureSpaceSize_ < features.size())
		operatorFeatureSpaceSize_ = features.size();

	for (size_t i=0; i<features.size(); i++)
	{
		if (numFeatureValues[i] < 2)
			continue;

		if (verboseLevel>1)
			cout << endl << "FEATURE " << i << " - " << features[i].getName() << endl;
		
		SplitOperator sop;
		double gain = auxilaryData->findOptimalSplit(i, sop.thresholds, -0.333, minRelativeInformationGain, 
									   maxNumSplits, 300, (verboseLevel>1) );

		if (verboseLevel == 1)
		{
			cout << "F" << i << "\t" << features[i].getName() << "\tgain:" << fixed << setprecision(4) << gain;
			cout << " [";
			for (size_t i=0; i<sop.thresholds.size(); i++)
				cout << " " << sop.thresholds[i];
			cout << " ]" << endl;
		}		

		if (sop.thresholds.size()>0)
		{
			// add for now
			sop.sourceIdx = i;
			sop.indexesForBinValues.resize(sop.thresholds.size()+1,MAX_UINT);
			for (size_t i=0; i<=sop.thresholds.size(); i++)
				sop.indexesForBinValues[i]=operatorFeatureSpaceSize_++;
		
			if (addBooleanIndicators)
			{
				sop.indexesForBinIndicators.resize(sop.thresholds.size()+1,MAX_UINT);
				for (size_t i=0; i<=sop.thresholds.size(); i++)
					sop.indexesForBinIndicators[i]=operatorFeatureSpaceSize_++;
			}
			
			executionOrder_.push_back(IdxPair(OT_SPLIT,splits_.size()));
			splits_.push_back(sop);
			
			if (verboseLevel>1)
			{
				cout << endl << "[" << setprecision(4);
				for (size_t j=0; j<sop.thresholds.size(); j++)
					cout << " " << sop.thresholds[j];
				cout << "]" << endl;
			}
		}
	}
}


/***********************************************************************
The creation of the operator involves several actions:
- adding new features
- adding the operator
- updating the SearchData to include the new features
- if the original feature is modified (because of added indicators), then
  the entries for the original feature must be updated too.
************************************************************************/
void MlOperatorList::createSplitOperator(size_t originalFeatureIdx,
							 const vector<value_t>& thresholds,
							 MlOperatorSearchData* auxilaryData,
							 bool  addIndicators)
{
	vector<size_t> binValueCounts;
	//auxilaryData->countNumValues(originalFeatureIdx, thresholds, binValueCounts);
}


void MlOperatorList::createNewOperators(const MlDataSet& ds,
							const vector<size_t>& missIndexes,
							size_t maxNumOperators)
{

}
