#include "SpectraList.h"
#include "PeakList.h"
#include "PepNovo_auxfun.h"

void SpectraList::selectAllAggregatorHeaders()
{
	headerLocations_.clear();

	const int numSpectraFiles = spectraAggregator_.getNumSpectraFiles();

	int totalNumHeaders=0;
	int fileIndex;
	for (fileIndex=0; fileIndex<numSpectraFiles; fileIndex++)
		totalNumHeaders += spectraAggregator_.getSpectraFile(fileIndex).getNumHeaders();

	headerLocations_.reserve(totalNumHeaders);

	for (fileIndex=0; fileIndex<numSpectraFiles; fileIndex++)
	{
		const int numHeaders = spectraAggregator_.getSpectraFile(fileIndex).getNumHeaders();
		for (int i=0; i<numHeaders; i++)
			headerLocations_.push_back(SpectrumLocationPair(fileIndex,i));
	}
	
}



size_t SpectraList::removeExcludedScans(const ScanListManager& slm)
{
	vector<SpectrumLocationPair> newHeaderLocations;
	const size_t numHeadersBefore = headerLocations_.size();
	size_t numExcluded = 0;
	for (size_t i=0; i<numHeadersBefore; i++)
	{
		const SpectrumLocationPair& p = headerLocations_[i];
		const SingleSpectrumHeader* header = spectraAggregator_.getSpectraFile(p.fileIndex).getHeader(p.indexInVector);
		assert(header);
		const ScanEntry se(header->getDatasetIndex(), header->getSpectraFileIndexInList(), header->getScanNumber(), header->getMOverZ());
		if (slm.checkScanList(se))
		{
			numExcluded++;
			continue;
		}
		newHeaderLocations.push_back(headerLocations_[i]);
	}

	headerLocations_ = newHeaderLocations;
	return numExcluded;
}




size_t SpectraList::keepInclusionScans(const ScanListManager& sem)
{
	vector<SpectrumLocationPair> newHeaderLocations;
	const size_t numHeadersBefore = headerLocations_.size();
	size_t numExcluded = 0;
	for (size_t i=0; i<numHeadersBefore; i++)
	{
		const SpectrumLocationPair& p = headerLocations_[i];
		const SingleSpectrumHeader* header = spectraAggregator_.getSpectraFile(p.fileIndex).getHeader(p.indexInVector);
		assert(header);
		const ScanEntry se(header->getDatasetIndex(), header->getSpectraFileIndexInList(), header->getScanNumber(), header->getMOverZ());
		if (! sem.checkScanList(se))
		{
			numExcluded++;
			continue;
		}
		newHeaderLocations.push_back(headerLocations_[i]);
	}

	headerLocations_ = newHeaderLocations;
	return numExcluded;
}


void SpectraList::selectHeaders(mass_t minMz, mass_t maxMz, int minCharge, int maxCharge, 
								float minSqs, float maxSqs)
{
	headerLocations_.clear();

	const int numSpectraFiles = spectraAggregator_.getNumSpectraFiles();

	int totalNumHeaders=0;
	int fileIndex;
	for (fileIndex=0; fileIndex<numSpectraFiles; fileIndex++)
		totalNumHeaders += spectraAggregator_.getSpectraFile(fileIndex).getNumHeaders();

	headerLocations_.reserve(totalNumHeaders);

	for (fileIndex=0; fileIndex<numSpectraFiles; fileIndex++)
	{
		const SpectraFile& file = spectraAggregator_.getSpectraFile(fileIndex);
		const int numHeaders = file.getNumHeaders();
		int i;

		for (i=0; i<numHeaders; i++)
		{
			const SingleSpectrumHeader* header = file.getHeader(i);
			if (header->getMOverZ() >= minMz     && header->getMOverZ() <= maxMz &&
				header->getCharge() >= minCharge && header->getCharge() <= maxCharge &&
				(header->getSqs()<0 || (header->getSqs()	>= minSqs	 && header->getSqs() <= maxSqs)) )
			{
					headerLocations_.push_back(SpectrumLocationPair(fileIndex,i));
			}
		}
	}	
}

void SpectraList::randomlyReduceListToSize(size_t newSize)
{
	if (newSize>= headerLocations_.size())
		return;

	vector<size_t> idxs;
	chooseKFromN(newSize, headerLocations_.size(), idxs);

	vector<SpectrumLocationPair> selectedPairs;
	selectedPairs.resize(newSize);

	for (size_t i=0; i<newSize; i++)
		selectedPairs[i] = headerLocations_[idxs[i]];

	headerLocations_ = selectedPairs;
}


bool SpectraList::checkFilesReadOk() const
{
	cout << "Checking " << getNumHeaders() << " spectra..." << endl;
	size_t numErrors=0;
	for (size_t i=0; i<getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* ssh = getSpectrumHeader(i);
		PeakList pl;
		pl.readPeaksToLocalAllocation(spectraAggregator_, ssh);		
		
		pl.getHeader()->printStats(0, cout, false);
		cout << "\t" << ssh->getOriginalNumPeaks() << " : " << pl.getNumPeaks();
		if (pl.getNumPeaks() == ssh->getOriginalNumPeaks())
		{
			cout << "\t+";
		}
		else
		{
			cout << "\t-";
			numErrors++;
		}

		if (pl.sanityCheck())
		{
			cout << " +";
		}
		else
		{
			cout << " -";
			numErrors++;
		}
		cout << endl;
	}

	cout << "Num errors found: " << numErrors << endl;
	return (numErrors == 0);
}

