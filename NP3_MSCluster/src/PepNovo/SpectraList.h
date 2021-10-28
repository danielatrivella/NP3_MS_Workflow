#ifndef __SPECTRALIST_H__
#define __SPECTRALIST_H__

#include "SpectraAggregator.h"
#include "ScanList.h"

struct SpectrumLocationPair {

	SpectrumLocationPair() : fileIndex(-1), indexInVector(-1) {}
	SpectrumLocationPair(int f, int i) : fileIndex(f), indexInVector(i) {}

	int fileIndex;
	int indexInVector;

	const bool operator< (const SpectrumLocationPair& rhs) const {
		return ( fileIndex < rhs.fileIndex ||
				(fileIndex == rhs.fileIndex && indexInVector<rhs.indexInVector));
	}
};



class SpectraList {
public:
	SpectraList(const SpectraAggregator& spectraAggregator) : spectraAggregator_(spectraAggregator) {
		headerLocations_.clear();
	}

	void selectAllAggregatorHeaders();

	size_t removeExcludedScans(const ScanListManager& slm);
	size_t keepInclusionScans(const ScanListManager& slm);

	void selectHeaders(	mass_t minMz=0,		mass_t maxMz=POS_INF,
						int	minCharge=0,	int maxCharge=POS_INF,
						float minSqs=0,		float maxSqs=POS_INF );

	void randomlyReduceListToSize(size_t newSize);


	int	getNumHeaders() const { return headerLocations_.size(); }

	const SingleSpectrumHeader* getSpectrumHeader(size_t i) const
	{
		const SpectrumLocationPair& slp = headerLocations_[i];
		return spectraAggregator_.getSpectraFile(slp.fileIndex).getHeader(slp.indexInVector);
	}


	// a debug function, reads files and checks peak counts
	bool	checkFilesReadOk() const; 

private:
	const SpectraAggregator& spectraAggregator_;
	vector<SpectrumLocationPair> headerLocations_;
};



#endif


