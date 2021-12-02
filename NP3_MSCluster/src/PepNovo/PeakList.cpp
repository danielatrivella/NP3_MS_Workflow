#include "PeakList.h"
#include "Isotopes.h"
#include "PepNovo_auxfun.h"


PeakList &PeakList::operator =(const PeakList & rhs)
{
	header_ = rhs.header_;
	config_ = rhs.config_;
	numPeaks_ = rhs.numPeaks_;
	if (localAllocationSize_>0)
	{
		if (numPeaks_ > localAllocationSize_)
		{
			localAllocationSize_ = (numPeaks_ > 10 ? numPeaks_ : 10);
			if (peaks_)
				delete [] peaks_;
			peaks_ = new Peak[numPeaks_];
		}
		memcpy(peaks_, rhs.peaks_, numPeaks_*sizeof(Peak));
	}
	else
		peaks_ = rhs.peaks_;

	totalPeakIntensity_ = rhs.totalPeakIntensity_;
	return *this;
}

void PeakList::copyPeakListLocally(const PeakList& pl) // copies the other peak list, peaks are copied locally
{
	header_ = pl.header_;
	config_ = pl.config_;
	numPeaks_ = pl.numPeaks_;
	
	if (numPeaks_ > localAllocationSize_)
	{
		localAllocationSize_ = (numPeaks_ > 10 ? numPeaks_ : 10);
		if (peaks_)
			delete [] peaks_;
		peaks_ = new Peak[numPeaks_];
	}
	memcpy(peaks_, pl.peaks_, numPeaks_*sizeof(Peak));


	totalPeakIntensity_ = pl.totalPeakIntensity_;
}

// returns number of peaks that were stored
int	PeakList::readPeaksToLocalAllocation(const SpectraAggregator& sa,	
										 const SingleSpectrumHeader* header)
{
	// create a basic PeakList read function 
	header_ = header;
	
	if (localAllocationSize_>0 && peaks_)
		delete [] peaks_;

	localAllocationSize_ = header->getOriginalNumPeaks();
	peaks_ = new Peak[localAllocationSize_];
	if (! peaks_)
	{
		cout << "Error: couldn't allocate memory for spectrum!" << endl;
		exit(1);
	}

	const int numPeaksRead = sa.readPeakList(header, peaks_);
	numPeaks_ = numPeaksRead;

	if (header->getFileType() != IFT_MZXML && 
		( numPeaksRead != header->getOriginalNumPeaks() || 
		  peaks_[0].mass != header->getFirstPeakMass()) )
	{
		cout << "Error reading scan " << header->getScanNumber() << ": " << header_->getTitle()
			<< " in file " << sa.getSpectraFile(header_->getSpectraFileIndexInList()).getFilePath() << endl;
		if (numPeaksRead != header->getOriginalNumPeaks())
			cout << "Num peaks read " << numPeaksRead << ", expecting " <<  header->getScanNumber() << endl;
		if (peaks_[0].mass != header->getFirstPeakMass())
			cout << setprecision(5) << "First peak mass " << peaks_[0].mass 
				 <<", expecting " << header->getFirstPeakMass() << endl;
		cout << "Could possibly be a dos/unix problem with the files, try running dos2unix (or unix2dos)..." << endl;
		cout <<"Skipping spectrum..." << endl;

		const int numPeaksRead = sa.readPeakList(header, peaks_);

		return 0;
	}

	return numPeaks_;
}

// function assumes that the buffer is sufficently large for all peaks being read
// returns number of peaks that were stored
int	PeakList::readPeaksToBuffer(const SpectraAggregator& sa,	
								const SingleSpectrumHeader* header,
								Peak*	peakBuffer)
{
	// create a basic PeakList read function 
	header_ = header;

	if (localAllocationSize_>0 && peaks_) {
        delete[] peaks_;

    }
	localAllocationSize_ = 0;
	
	peaks_ = peakBuffer;
	const int numPeaksRead = sa.readPeakList(header, peaks_);
	numPeaks_ = numPeaksRead;

	return numPeaks_;	
}


// writes the header and peak list in DAT format
// returns the number of bytes written to the buffer
size_t PeakList::writeToDatBuffer(char* buffer, const SingleSpectrumHeader* newHeader) const
{
	const size_t headerSize = newHeader->writeHeaderToDatBuffer(buffer);

	assert( sanityCheck() );

	char* p = buffer + headerSize;

	// write peaks
	for (size_t i=0; i<numPeaks_; i++)
	{
		mass_t* mt = reinterpret_cast<mass_t*>(p);
		*mt = peaks_[i].mass;
		intensity_t* it = reinterpret_cast<intensity_t*>(p+sizeof(mass_t));
		*it = peaks_[i].intensity;
		unsigned char* up = reinterpret_cast<unsigned char*>(p+sizeof(mass_t)+sizeof(intensity_t));
		*up++ = peaks_[i].count;
		*up++ = peaks_[i].maxPossible;
		unsigned short* us = reinterpret_cast<unsigned short*>(up);
		*us++ = peaks_[i].charge;
		
		p = reinterpret_cast<char*>(us);

	/*	if (peaks_[i].maxPossible > newHeader->getClusterSize())
		{
			cout << "Cluster size: " << newHeader->getClusterSize() << endl;
			cout << "Num peaks: " << numPeaks_ << endl;
			cout << "Original : " << newHeader->getOriginalNumPeaks() << endl;
			this->printPeaks();
		}
		assert(peaks_[i].maxPossible <= newHeader->getClusterSize());*/
	}

	unsigned int* ui = reinterpret_cast<unsigned int*>(buffer);
	*ui = static_cast<unsigned int>(p-buffer);

	return (p-buffer);
}


// perfrom binary search to find first index that has mass >= minRange
PeakRange  PeakList::findPeaksInRange(mass_t minRange, mass_t maxRange) const
{
	const int maxPeakIndex = numPeaks_-1;
	int low, high;
    low = 0; high = maxPeakIndex;
    while(low < high) 
	{
		const int mid = (low+high)/2;
		if(minRange < peaks_[mid].mass)
		{
			high = mid-1;
		}
		else 
			low  = mid+1;
    }

	// make sure that low is brought to an element that is less than minRange
	while (low>0 && peaks_[low].mass > minRange)
		--low;

	// move low to the first element larger than minRange
	while (low<numPeaks_ && peaks_[low].mass < minRange)
		++low;

	PeakRange pr;
	if (low == numPeaks_ || peaks_[low].mass > maxRange)
		return pr; // no peaks in range

	pr.low_idx = low;
	pr.high_idx = low;
	while (pr.high_idx<maxPeakIndex)
	{
		if (peaks_[++pr.high_idx].mass > maxRange)
		{
			--pr.high_idx;
			break;
		}
	}
	pr.num_peaks = pr.high_idx - pr.low_idx + 1;

	return pr;
}


int	 PeakList::findPeakWithMaxIntensity(mass_t expectedMass, mass_t tolerance) const
{
	PeakRange pr = findPeaksInRange(expectedMass-tolerance, expectedMass+tolerance);
	if (pr.num_peaks==0)
		return -1;

	if (pr.num_peaks ==1)
		return pr.low_idx;

	const float doubleTolerance = tolerance + tolerance;

	mass_t maxWeightedValue = (doubleTolerance - fabs(peaks_[pr.low_idx].mass - expectedMass)) *
							   peaks_[pr.low_idx].intensity;

	int maxPeakIdx=pr.low_idx;
	
	// find closest peak to exp_pos
	for (int i=1; i<pr.num_peaks; i++)
	{
		const int peakIdx = pr.low_idx + i;
		const mass_t weightedValue = (doubleTolerance - fabs(peaks_[peakIdx].mass - 
								expectedMass)) * peaks_[peakIdx].intensity;

		if (weightedValue>maxWeightedValue)
		{
			maxWeightedValue=weightedValue;
			maxPeakIdx = peakIdx;
		}
		else
			break;
	}
	return maxPeakIdx;
}

void PeakList::createIndexArray(vector<int>& indexArray) const
{
	const mass_t maxPeakMass = peaks_[numPeaks_-1].mass;
	const mass_t maxMassToConsider = 10.0 + (header_->getOriginalPmWith19() > maxPeakMass ?
		header_->getOriginalPmWith19() : maxPeakMass);
	const int arraySize = static_cast<int>(maxMassToConsider + 20.0);
	
	const int maxPeakIndex = numPeaks_-1;
	
	indexArray.clear();
	indexArray.resize(arraySize+1,maxPeakIndex);
	int i=0;
	int m=static_cast<int>(peaks_[0].mass);
	while (i<m)
		indexArray[i++]=0;

	int c=0;
	while (c< maxPeakIndex)
	{
		int curr_m = static_cast<int>(peaks_[c].mass);
		int next_m = curr_m;
		int next_c = c;

		while (next_m == curr_m && next_c<maxPeakIndex)
			next_m=static_cast<int>(peaks_[++next_c].mass);

		if (next_m >= arraySize)
			next_m =  arraySize-1;

		while (i<next_m)
			indexArray[i++]=c;
		
		c=next_c;
	}
}



void PeakList::calculatePeakRanks(vector<int>& peakRanks) const
{	
	vector<IndexIntensityPair> pairs;
	pairs.resize(numPeaks_);

	int i;
	for (i=0; i<numPeaks_; i++)
	{
		pairs[i].index=i;
		pairs[i].intensity=peaks_[i].intensity;
	}

	sort(pairs.begin(), pairs.end());
	peakRanks.resize(pairs.size());

	for (i=0; i<pairs.size(); i++)
		peakRanks[pairs[i].index]=i+1;
}


void PeakList::calculateLogLocalRanks(mass_t windowSize, vector<float>& logLocalRanks) const
{

	const mass_t halfWindowSize = windowSize * 0.5;	
	logLocalRanks.resize(numPeaks_);	
	
	int i;
	for (i=0; i<numPeaks_; i++)
	{
		const PeakRange pr= findPeaksInRange(peaks_[i].mass - halfWindowSize,
											 peaks_[i].mass + halfWindowSize);
		int above=0;
		int j;
		for (j=pr.low_idx; j<=pr.high_idx && j<numPeaks_; j++)
			if (peaks_[j].intensity>peaks_[i].intensity)
				above++;

		logLocalRanks[i] = log(1.0 + static_cast<float>(above));
	}	
}


void PeakList::normalizePeakIntensities()
{
	
	// TODO:
	// remove the intensity of pm+20 / 2 and pm+2/2 from the total intensity
	// if this is an issue, normalization might have to be done twice,
	// before and after parent mass correction

	if (totalPeakIntensity_<=0.0)
	{
		totalPeakIntensity_=0.0;
		int i;
		for (i=0; i<numPeaks_; i++)
			totalPeakIntensity_+=peaks_[i].intensity;
	}

	const float normalizationValue = 1000.0 / totalPeakIntensity_;
	totalNormalizedPeakIntensity_ = 1000.0;
	int i;
	for (i=0; i<numPeaks_; i++)
		peaks_[i].intensity *= normalizationValue;
}


void PeakList::initializePeakList(const Config* config, bool indFilterSpectrum)
{
	config_ = config;
	joinAdjacentPeaks(config->getTolerance());
	
	// TODO : filterPeaks, does it get the corrected PM in normal operation?
	// are we missing out on b/y pairs because of this?
	if (indFilterSpectrum)
	{
		const mass_t pmWith19 = (header_->getPmWith19()>0 ? header_->getPmWith19() : header_->getOriginalPmWith19());
		filterWeakPeaks(config, pmWith19);
	}

	totalPeakIntensity_=0;
	if (header_->getClusterSize()<=1)
	{
		for (size_t i=0; i<numPeaks_; i++)
		{
			totalPeakIntensity_+=peaks_[i].intensity;
			peaks_[i].count = 1;
			peaks_[i].maxPossible = 1;
		}
	}
	else
		for (size_t i=0; i<numPeaks_; i++)
			totalPeakIntensity_+=peaks_[i].intensity;
	

	normalizePeakIntensities();
	
	setHeaderFirstPeakMass(peaks_[0].mass);
}


void PeakList::computeLogIntensities( vector<float>& logIntensities)
{
	logIntensities.clear();
	logIntensities.resize(numPeaks_,0);

	int i;
	for (i=0; i<numPeaks_; i++) {
		// NP3 GOT scale intensities
		if (config_->get_scale_factor() == 0.0) {
			logIntensities[i] = log(1.0 + peaks_[i].intensity);
		} else {
			logIntensities[i] = pow(peaks_[i].intensity, config_->get_scale_factor());
		}
	}
}


// TODO what about charge 2 isotopes, etc. if tolerance is low enough, we can check for them too
void PeakList::calculateIsotopicLevels(mass_t tolerance, 
									   vector<float>& isotopicLevels) const
{
	const mass_t isotopicTolerance = (tolerance<0.2 ? tolerance : 0.2);
	const int lastPeakIndex = numPeaks_-1;

	isotopicLevels.clear();
	isotopicLevels.resize(numPeaks_,0);
	int i;
	for (i=0; i<lastPeakIndex; i++)
	{	
		if (isotopicLevels[i]>0)
			continue;

		// look for +1 peak
		int idx1 = findPeakWithMaxIntensity(peaks_[i].mass + MASS_PROTON, isotopicTolerance);
		if (idx1<0)  
			continue;

		const float oneOverIntensity = 1.0 / peaks_[i].intensity;
		float ratio1 = peaks_[idx1].intensity * oneOverIntensity;

		// ignore strong +1
		if ( ratio1 > 3.0)
			continue;

		// examine ratios
		vector<float> expectedPeakRatios, observedPeakRatios, relativeRatios;
		vector<int> isotopicIndexs;

		observedPeakRatios.resize(6);
		observedPeakRatios[0]=1.0;
		observedPeakRatios[1]= ratio1;

		isotopicIndexs.resize(6);
		isotopicIndexs[0]=i;
		isotopicIndexs[1]=idx1; // peak index at +1 position

		// find additional peaks
		int j;
		for (j=2; j<=5; j++)
		{
			int idx = findPeakWithMaxIntensity(peaks_[i].mass + j*MASS_ISO, isotopicTolerance);
			if (idx<0)
				break;
			observedPeakRatios[j] = peaks_[idx].mass * oneOverIntensity;
			isotopicIndexs[j]=idx;
		}
		const int lastIsotopicPosition = j-1;

		// get expected iso ratios
		calc_expected_iso_ratios(peaks_[i].mass,expectedPeakRatios,j);

		// calc ratios between observed and expected		
		relativeRatios.clear();
		relativeRatios.resize(6,0);
		relativeRatios[0]=1.0;
		for (j=1; j<=lastIsotopicPosition; j++)
			if (expectedPeakRatios[j]>0)
				relativeRatios[j]=observedPeakRatios[j] / expectedPeakRatios[j];

		float levelDiscount=1.0; // as we move farther from the monoisotopic peak, the rates down
		for (j=1; j<= lastIsotopicPosition; j++)
		{
			float isotopicLevel;

			if (relativeRatios[j]>= 0.75 && relativeRatios[j]<=1.333)
			{
				isotopicLevel=2.0;
			}
			else if (relativeRatios[j] >= 0.5 && relativeRatios[j] <=2.0)
			{
				isotopicLevel=1.3333;
			}
			else if (relativeRatios[j] >= 0.3333 && relativeRatios[j] <=3.0)
			{
				isotopicLevel=0.6666;
			}
			else if (relativeRatios[j] >= 0.25 && relativeRatios[j] <= 4.0)
			{
				isotopicLevel=0.3333;
			}
			else
				break;
			
			isotopicLevels[isotopicIndexs[j]] = isotopicLevels[isotopicIndexs[j-1]] + 
												levelDiscount * isotopicLevel;

			levelDiscount *= 0.5;
		}
	}
}

/*
void PeakList::calculateMonoisotpicRanks(const vector<float>& isotopicLevels,
											   vector<int>&	  monoisotopicRanks) const
{
	vector<IndexIntensityPair> monoisotopicPairs, isotopicPairs;
	int i;
	for (i=0; i<numPeaks_; i++)
	{
		IndexIntensityPair iip(i,peaks_[i].intensity);
		if (isotopicLevels[i] <= 0)
		{
			monoisotopicPairs.push_back(iip);
		}
		else
			isotopicPairs.push_back(iip);

	}

	sort(monoisotopicPairs.begin(),monoisotopicPairs.end());
	sort(isotopicPairs.begin(),isotopicPairs.end());

	monoisotopicRanks.resize(numPeaks_);
	
	for (i=0; i<monoisotopicPairs.size(); i++)
		monoisotopicRanks[i]=monoisotopicPairs[i].index;

	for (i=0; i<isotopicPairs.size(); i++)
		monoisotopicRanks[monoisotopicPairs.size()+i]=isotopicPairs[i].index;
}
*/

void PeakList::selectStrongPeakIndexes(int numStrongPeaksPerLocalWindow,
									   const vector<float>& logLocalRanks,
									   const vector<float>& isotopicLevels,
									   vector<int>& strongIndexes) const
{
	const float thresholdLogLevel = log(static_cast<float>(numStrongPeaksPerLocalWindow));

	strongIndexes.clear();
	strongIndexes.reserve(numPeaks_);
	int i;
	for (i=0 ;i<numPeaks_; i++)
		if (logLocalRanks[i]<=thresholdLogLevel && isotopicLevels[i]<=0.0)
			strongIndexes.push_back(i);
}


// There are fixed magic numbers in this function that should not be changed
// (even though actual tolerances or values in config file might differ)
// changing these numbers without retraining all models might give unanticipated results
void PeakList::calculateLogRandomProbabilities(vector<float>& logIntensities,
											   vector<float>& logRandomProbabilities) const
{
	if (numPeaks_<2)
		return;

	const float oneOverSqrt2pi = 1.0 / sqrt(2.0*3.1415927);
	const float logBias = log(1.2); // to fine tune the zero probabilites
	const mass_t peakWindowSize = 0.6; // this is fixed and independent of the tolerance!!!
	const mass_t margin      = 25.0;   // fixed and independent of values in config file
	const mass_t windowSize = 100.0;   // fixed and independent of values in config file	
	const mass_t minPeakMass = peaks_[0].mass;
	const mass_t maxPeakMass = peaks_[numPeaks_-1].mass;
	const mass_t visibleRange = (maxPeakMass - minPeakMass);

	logRandomProbabilities.resize(numPeaks_);	
	int i;
	for (i=0; i<numPeaks_; i++)
	{
		const mass_t peakMass = peaks_[i].mass;
		const mass_t relativePosition = (peakMass-minPeakMass)/visibleRange;
		const mass_t leftWindow  = margin + relativePosition * windowSize;
		const mass_t rightWindow = margin + windowSize - leftWindow;
		const PeakRange pr = findPeaksInRange(peakMass - leftWindow, peakMass + rightWindow);
		const float peakWindowProb = peakWindowSize /(leftWindow + rightWindow);

		// some freak cases have 0 peak counts (only in unix)
		const int numPeaksInRange = (pr.num_peaks>0 ? pr.num_peaks : 1);
		const float zeroProbability = pow((1.0 - peakWindowProb),numPeaksInRange); 

		if (numPeaksInRange<5)
		{
			logRandomProbabilities[i] = log(1.0-zeroProbability) + logBias;
		}
		else // compute special probability based on peak densitiy model
		{
			vector<float> windowLogIntensities;
			int j;
			for (j=pr.low_idx; j<=pr.high_idx; j++)
				windowLogIntensities.push_back(logIntensities[j]);
		
			float mean=0,sd=1;
			calc_mean_sd(windowLogIntensities,&mean,&sd);
			const float e = (logIntensities[i] - mean)/sd;
			if (e<0)
			{
				logRandomProbabilities[i] = log(1 - zeroProbability) + logBias;
			}
			else
			{
				const float normalizedValue = (oneOverSqrt2pi/ sd) * exp(-0.5*e*e);
				const float normalizationConstant = (1.0 - zeroProbability) / (oneOverSqrt2pi / sd); //
				logRandomProbabilities[i] = log(normalizedValue*normalizationConstant);
			}
		} 
	}
}


void PeakList::joinAdjacentPeaks(mass_t tolerance, bool indAddCounts)
{
	// NP3 GOT join adjacent peaks changed < 3 to < 2
	if (numPeaks_ < 2)
		return;

	const mass_t maxProximity = tolerance * 0.5;
	static vector<Peak> tmpPeakArea;
	if (tmpPeakArea.size() < numPeaks_)
	{
		int size = numPeaks_ * 2;
		if (size<5000)
			size=5000;
		tmpPeakArea.resize(size);
	}
	
	
	tmpPeakArea[0]=peaks_[0];
	int startIdx=1; // skip the first peak since it needs to be kept with the same mass
				// for debug purposes (the header's first peak mass needs to be the
				// same as what is read in the files
    // NP3 GOT change < to <=
	while (startIdx<numPeaks_ && peaks_[startIdx].mass <= tmpPeakArea[0].mass + maxProximity)
	{
		tmpPeakArea[0].intensity   += peaks_[startIdx].intensity;

		if (indAddCounts)
		{
			// NP3 fixed an upper bound for the peaks counts and maxPossible to prevent it from going above the char type limit = 255
			if (tmpPeakArea[0].count + peaks_[startIdx].count > 255)
				tmpPeakArea[0].count = 255;
			else
				tmpPeakArea[0].count += peaks_[startIdx].count;

			if (tmpPeakArea[0].maxPossible + peaks_[startIdx].maxPossible > 255)
				tmpPeakArea[0].maxPossible = 255;
			else
				tmpPeakArea[0].maxPossible += peaks_[startIdx].maxPossible;
		}
		startIdx++;
	}

	if (startIdx<numPeaks_)
		tmpPeakArea[1]=peaks_[startIdx];
	
	int prev = 1;
	for (int i=startIdx+1; i<numPeaks_; i++)
	{
        // NP3 GOT change < to <=
		if 	(peaks_[i].mass - tmpPeakArea[prev].mass <= maxProximity &&
			 peaks_[i].charge == tmpPeakArea[prev].charge )
		{
			// join peaks with proportion to their intensities

			const intensity_t intensitySum=(tmpPeakArea[prev].intensity + peaks_[i].intensity);
			const mass_t ratio = tmpPeakArea[prev].intensity/intensitySum;
			const mass_t newMass = ratio *tmpPeakArea[prev].mass + (1.0-ratio)*peaks_[i].mass;
			
			tmpPeakArea[prev].intensity = intensitySum;
			tmpPeakArea[prev].mass		= newMass;
			if (indAddCounts)
			{
			    // NP3 fixed an upper bound for the peaks counts and maxPossible to prevent it from going above the char type limit = 255
				if (tmpPeakArea[prev].count + peaks_[i].count > 255)
					tmpPeakArea[prev].count = 255;
				else
					tmpPeakArea[prev].count += peaks_[i].count;

                if (tmpPeakArea[prev].maxPossible + peaks_[i].maxPossible > 255)
                    tmpPeakArea[prev].maxPossible = 255;
                else
                    tmpPeakArea[prev].maxPossible += peaks_[i].maxPossible;
			}
		}
		else
			tmpPeakArea[++prev]=peaks_[i];
	}
	// NP3 GOT increment prev to not skip the last peak - only doing it when startIdx is not at the end of the spectra
	// this means that prev should be > 1, e.g. the spectra must have more than one peak
    if (startIdx<numPeaks_)
	    prev++;
	if (prev != numPeaks_)
	{
		// need to copy over peaks
		numPeaks_ = prev;
		memcpy(peaks_, &tmpPeakArea[0], numPeaks_ * sizeof(Peak));
	}
}


void PeakList::filterWeakPeaks(const Config* config, mass_t pmWith19, int peakDensity, bool removeGloballyWeak)
{
	if (numPeaks_<10)
		return;

	const mass_t maxMassToConsider = 10.0 + (pmWith19 > peaks_[numPeaks_-1].mass ?
		pmWith19 : peaks_[numPeaks_-1].mass);

	const mass_t windowSize		  = 0.5 * config->get_local_window_size();
	const int	 numPeaksInWindow = (peakDensity>0 ? peakDensity : config->get_max_number_peaks_per_local_window());
	const mass_t tolerance		  = config->getTolerance();

	vector<bool> keep_peaks(numPeaks_, false);	
	vector<Peak> new_peaks;

	int maximalPeakIndex = numPeaks_ -1;

	// keep first peak and last peak
	keep_peaks[0]=true;
	if (peaks_[maximalPeakIndex].mass<maxMassToConsider)
		keep_peaks[maximalPeakIndex]=true;
	
	int min_idx=1;
	int max_idx=1;

	// check the rest of the peaks
	int i;
	for (i=1; i<maximalPeakIndex; i++)
	{
		mass_t peak_mass=peaks_[i].mass;
		mass_t min_mass=peaks_[min_idx].mass;
		mass_t max_mass=peaks_[max_idx].mass;

		if (peaks_[i].mass > maxMassToConsider)
			break;

		// advance min/max pointers
		while (peak_mass-min_mass > windowSize)
			min_mass=peaks_[++min_idx].mass;

		while (max_idx < maximalPeakIndex && max_mass - peak_mass <= windowSize)
			max_mass=peaks_[++max_idx].mass;

		if (max_mass - peak_mass > windowSize)
			max_idx--;

		// this peak might already be marked for keeping (isotpoic peak)
		if (keep_peaks[i])
			continue;

		// if there are less than the maximum number of peaks in the window, keep it.
		if (max_idx-min_idx < numPeaksInWindow)
		{
			keep_peaks[i]=true;
			continue;
		}

		// check if this is one of the top peaks in the window
		int higher_count=0;
		for (int j=min_idx; j<=max_idx; j++)
			if (peaks_[j].intensity > peaks_[i].intensity)
				higher_count++;

		if (higher_count < numPeaksInWindow)
		{
			keep_peaks[i]=true;
		}
	}


	if (pmWith19>0)
	{
		// look for b/y pairs
	//	mass_t pm_with_20 = (originalPmWith19_>0 ? correctedPmWith19_ : originalPmWith19_) + 
	//						 MASS_PROTON ;
		mass_t pm_with_20 = pmWith19 + MASS_PROTON;

		mass_t pm_with_20_upper  = pm_with_20 + tolerance;
		mass_t pm_with_20_lower = pm_with_20 - tolerance;

		int f_idx =0;
		int b_idx = numPeaks_-1;
		while (f_idx<numPeaks_ && b_idx>=0)
		{
			if (! keep_peaks[f_idx])
			{
				f_idx++;
				continue;
			}

			while (b_idx>=0 && peaks_[f_idx].mass + peaks_[b_idx].mass > pm_with_20_upper )
				b_idx--;

			if (b_idx<0)
				break;

			mass_t mass_sum = peaks_[f_idx].mass + peaks_[b_idx].mass;
			if (mass_sum > pm_with_20_lower && mass_sum < pm_with_20_upper)
			{
				keep_peaks[f_idx]=true;
				keep_peaks[b_idx]=true;
			}
			f_idx++;
		}
	}

	// copy peaks
	int j=0;
	for (size_t i=0; i<numPeaks_; i++)
		if (keep_peaks[i] && peaks_[i].intensity>=0.001)
			peaks_[j++]=peaks_[i];
	numPeaks_ = j;

	// filter very low intensity peaks (without a window)
	if (removeGloballyWeak)
	{
		vector<intensity_t> intens(numPeaks_);
		for (size_t i=0; i<numPeaks_; i++)
			intens[i]=peaks_[i].intensity;
		sort(intens.begin(), intens.end());

		j=1;
		const intensity_t minIntensity = intens[(2*numPeaks_)/3] * 0.001;
		for (size_t i=1; i<numPeaks_; i++)
			if (peaks_[i].intensity > minIntensity)
				peaks_[j++]=peaks_[i];
		numPeaks_ = j;
	}
	
	
}

void  PeakList::markAllPossibleIsotopicPeaks(mass_t tolerance, 
											 vector<bool>& isotopicIndicators) const
{
	const mass_t isotopicTolerance = (tolerance< 0.25) ? tolerance : 0.25;
	const mass_t maxIsotopicDiff = MASS_ISO + isotopicTolerance;
	const mass_t minIsotopicDiff = MASS_ISO - isotopicTolerance;

	isotopicIndicators.clear();
	isotopicIndicators.resize(numPeaks_, false);

	const int lastPeakIndex = numPeaks_-1;

	int i;
	for (i=0; i<lastPeakIndex; i++)
	{	
		// look for +1 peak
		if (peaks_[i].intensity <=0)
		{
			isotopicIndicators[i]=true;
			continue;
		}
		
		if (peaks_[i+1].mass - peaks_[i].mass>1.25)
			continue;

		int forwardIndex=i+1;
		while (forwardIndex < numPeaks_)
		{
			const mass_t massDiff = peaks_[forwardIndex].mass - peaks_[i].mass;
			if (massDiff > maxIsotopicDiff)
				break;

			if (massDiff > minIsotopicDiff)
				isotopicIndicators[i+1]=true;
			++forwardIndex;
		}
	}
}


// A fast heuristic function for marking strong peaks in the spectrum
// Uses a heursitic approach which makes it suitbale for applications
// like MSCluster
bool PeakList::heuristicMarkStrogenstPeaks(mass_t windowSize,
										   int numPeaksPerWindow,
										   vector<bool>& strongPeakIndicators) const
{
	strongPeakIndicators.clear();
	if (numPeaks_<=5)
	{
		strongPeakIndicators.resize(numPeaks_, true);
		return false;
	}

	// filter low intensity noise
	// and mark those that are good peaks
	const mass_t halfWindowSize = 0.5 * windowSize;
	const int maxPeakIndex = numPeaks_ -1;

	strongPeakIndicators.resize(numPeaks_,false);
	strongPeakIndicators[0]=true;
	strongPeakIndicators[maxPeakIndex]=true;

	int startWindowIndex =0;
	while (startWindowIndex<maxPeakIndex)
	{
		const mass_t maxWindowMass = peaks_[startWindowIndex].mass + windowSize;

		int endWindowIndex=startWindowIndex;
		while (endWindowIndex<maxPeakIndex && peaks_[endWindowIndex].mass<maxWindowMass)
			++endWindowIndex;


		if (endWindowIndex - startWindowIndex > numPeaksPerWindow)
		{
			const int numPeaksInWindow = endWindowIndex - startWindowIndex+1;
			vector<IndexIntensityPair> pairs;
			pairs.clear();
			int i;
			for (i=0; i<numPeaksInWindow; i++)
			{
				const int peakIdx = i+startWindowIndex;
				pairs.push_back(IndexIntensityPair(peakIdx,peaks_[peakIdx].intensity));
			}

			sort(pairs.begin(),pairs.end());

			for (i=0; i<numPeaksPerWindow; i++)
				strongPeakIndicators[pairs[i].index]=true;	
		}
		else 
		{
			int i;
			for (i=startWindowIndex; i<=endWindowIndex; i++)
				strongPeakIndicators[i]=true;
		}

		// advance half a window
		const mass_t midMass = peaks_[startWindowIndex].mass + halfWindowSize;
		++startWindowIndex;
		while (startWindowIndex<maxPeakIndex && peaks_[startWindowIndex].mass<midMass)
			startWindowIndex++;

	}
	return true;
}


bool PeakList::selectStrongPeakIndexesForPmcsqs(const vector<float>& isotopicLevels, 
												vector<bool>& strongPeakIndicators) const
{
	// Uses fixed values (window size 120, number peaks 3) that are not influenced 
	// by the model file. This is done to avoid unexpexted consequences (like erroneous
	// quality scores) just because somebody fiddled with the parameters.
	if ( ! this->heuristicMarkStrogenstPeaks(120, 3, strongPeakIndicators))
		return false;

	// also mark the top 20 peaks (non isotopic)
	vector<IndexIntensityPair> pairs;
	pairs.resize(numPeaks_);
	int i;
	for (i=0; i<numPeaks_; i++)
		pairs[i]=IndexIntensityPair(i,peaks_[i].intensity);
	
	sort(pairs.begin(),pairs.end());

	const int halfNumPeaks = numPeaks_/2;
	const int maxToMark = (halfNumPeaks<20 ? halfNumPeaks : 20);
	for (i=0; i<halfNumPeaks; i++)
		strongPeakIndicators[pairs[i].index]=true;

	for (i=0; i<numPeaks_; i++)
		strongPeakIndicators[i] = (strongPeakIndicators[i] && isotopicLevels[i]==0);

	return true;
}


bool PeakList::sanityCheck() const
{
	if (numPeaks_ <= 0) {
		cout << "num peaks " << numPeaks_ << endl;
		return false;
	}

	if (peaks_[0].intensity<0.0 || peaks_[0].mass<=0.0) {
		cout << "p0 int " << peaks_[0].intensity << " p0 mass " << peaks_[0].mass << endl;
		return false;
	}

	for (size_t i=1; i<numPeaks_; i++) {
        if (peaks_[i].mass <= peaks_[i - 1].mass || peaks_[i].intensity < 0.0) {
			cout << "problem with peak " << i << endl;
            printPeaks();
            return false;
        }
    }
	return true;
}

bool PeakList::fullSanityCheck(size_t clusterSize) const
{
    // NP3 GOT changed < 3 to 1
	if (numPeaks_ < 1) {
		cout << "num peaks " << numPeaks_ << endl;
		return false;
	}

	if (peaks_[0].intensity<0.0 || peaks_[0].mass<=0.0 || (peaks_[0].count > peaks_[0].maxPossible)) {
		cout << "p0 int " << peaks_[0].intensity << " p0 mass " << peaks_[0].mass << endl;
		return false;
	}


	for (size_t i=1; i<numPeaks_; i++)
		if (peaks_[i].mass <= peaks_[i-1].mass || peaks_[i].intensity<0.0 ||
			(peaks_[i].count > peaks_[i].maxPossible) ||
			 (clusterSize>0 && peaks_[i].count > clusterSize))
			return false;

	
	return true;
}

void PeakList::adjustPeakCounts(size_t clusterSize )
{
	for (size_t i=0; i<numPeaks_; i++)
	{
		if (peaks_[i].maxPossible>clusterSize)
			peaks_[i].maxPossible = clusterSize;
		if (peaks_[i].count > peaks_[i].maxPossible)
			peaks_[i].count = peaks_[i].maxPossible;
	}
}

void PeakList::printPeaks() const
{
	// NP3 GOT precision from 3 to 5
	cout << setprecision(5) << fixed;
	int i;
	for (i=0; i<numPeaks_; i++)
		cout << i << "\t" << peaks_[i].mass << "\t" << peaks_[i].intensity << "\t"
			 << static_cast<int>(peaks_[i].count) << "\t" << static_cast<int>(peaks_[i].maxPossible)
			 << "\t" << peaks_[i].charge << endl;
}

void PeakList::printSomePeaks() const
{
	// NP3 GOT precision from 3 to 5
	cout << header_->getTitle() << endl;
	cout << "num peaks: " << numPeaks_ << endl << fixed << setprecision(5);
	for (size_t i=0; i<4; i++)
	{
		cout << i << "\t";
		peaks_[i].print();
	}

	for (size_t i=numPeaks_-4; i<numPeaks_; i++)
	{
		cout << i << "\t";
		peaks_[i].print();
	}
	cout << endl;
}


void PeakList::clear()
{
	if (localAllocationSize_>0 && peaks_)
		delete [] peaks_;

	header_=NULL;
	numPeaks_=0;
	peaks_=NULL;
	localAllocationSize_=0;
}


// outputs an mgf spectrum BEGIN IONS ... END IONS to a buffer
// returns the number of bytes written (assumes there is enough room)
size_t	PeakList::outputToMgfBuffer(char* buffer, const SingleSpectrumHeader* newHeader) const
{
	static char beginLine[]={"BEGIN IONS\n"};
	char* p = const_cast<char*>(buffer);
	
	// write header
	memcpy(p, beginLine, sizeof(beginLine));
	p+= sizeof(beginLine)-1;
	if (newHeader)
	{
		p+= newHeader->writeHeaderToMgfBuffer(p);
	}
	else
		p+= header_->writeHeaderToMgfBuffer(p);
	// GOT BASELINE INTENSITY
	/*if (newHeader)
    {
        if (newHeader->getPrecursorIntensity() < 1000.0)
            return 0;
        p+= newHeader->writeHeaderToMgfBuffer(p);
    }
    else
    {
     if (header_->getPrecursorIntensity() < 1000.0)
        return 0;
       p+= header_->writeHeaderToMgfBuffer(p);
    }*/

	// write peaks
	ostringstream oss;
	oss << fixed << setprecision(NUM_SIG_DIGITS);

	// to avoid problems that could arise because we changed the peak intensities
	// (e.g., scoring by InsPecT?) we change the intensities back to reflect
	// the original non-normalized peak intensities
	
	if (totalPeakIntensity_ > totalNormalizedPeakIntensity_)
	{
		assert(totalNormalizedPeakIntensity_ > 0.0);
		const float weight = totalPeakIntensity_ / totalNormalizedPeakIntensity_;
		for (size_t i=0; i<numPeaks_; i++)
		{
			oss << peaks_[i].mass << " " << peaks_[i].intensity*weight;
			if (peaks_[i].charge>0)
				oss << " " << peaks_[i].charge;
			oss << endl;
		}
	}
	else
		for (size_t i=0; i<numPeaks_; i++)
		{
			oss << peaks_[i].mass << " " << peaks_[i].intensity;
			if (peaks_[i].charge>0)
				oss << " " << peaks_[i].charge;
			oss << endl;
		}

	oss << "END IONS" << endl << endl;

	const size_t len = oss.str().length();
	memcpy(p, oss.str().c_str(), len);
	p+= len;
	
	return (p-buffer);
}



