#ifndef __PEAKLIST_H__
#define __PEAKLIST_H__

#include "BasicDataStructs.h"
#include "SpectraAggregator.h"

////////////////////////////////////////////////////////////////////////////
/// This is the basic class for all spectra in PepNovo/MSCluster. This allows
/// for the use of common function in the different applications, despite the
/// big differences in the types of the inherted spectrum classes.
/// It holds a pointer to the spectrum's header, number of peaks
/// and a pointer to the peaks in the spectrum. There is an optional
/// array of peaks that can be used if the peaks are not allocated from a central
/// peak list.
/// A note about object distruction. There is no new type allocation performed  
/// by this object, so none of the pointers should be deleted upon destruction.
class PeakList {
public:
	PeakList() : header_(0), config_(0), numPeaks_(0), peaks_(0), localAllocationSize_(0),
				totalPeakIntensity_(0.0), totalNormalizedPeakIntensity_(0.0) {};

	~PeakList() { if (localAllocationSize_>0 && peaks_) delete [] peaks_; }

	PeakList& operator =(const PeakList &); // if memory is already allocated in the local allocation
											// the peaks are written to the local peak area

	void copyPeakListLocally(const PeakList& pl); // copies the other peak list, peaks are copied locally

	const SingleSpectrumHeader* getHeader()   const { return header_; }
	void  setHeader(const SingleSpectrumHeader* header) { header_ = header; }
	const Config*			    getConfig()   const { return config_; }
	void  setConfig(const Config* config) {config_ = config; }
	Peak const*			  getPeaks()     const { return peaks_; }
	Peak const&			  getPeak(int i) const { return peaks_[i]; }
	int					  getNumPeaks()  const { return numPeaks_; }
	void				  setNumPeaks(int i) { numPeaks_ = i; }
	void				  setPeaksPtr(Peak* newPtr) { peaks_ = newPtr; }

	mass_t      getPeakMass(int i) const { return peaks_[i].mass; }
	intensity_t getPeakIntensity(int i) const { return peaks_[i].intensity; }
	intensity_t getTotalPeakIntensity() const { return totalPeakIntensity_; }
	intensity_t getTotalNormalizedPeakIntensity() const { return totalNormalizedPeakIntensity_; }
	int			getLocalAllocationSize() const { return localAllocationSize_; }

	// returns number of peaks that were stored
	int	readPeaksToLocalAllocation(const SpectraAggregator& sa,	
								   const SingleSpectrumHeader* header);

	// function assumes that the buffer is sufficently large for all peaks being read
	// returns number of peaks that were stored
	int	readPeaksToBuffer(const SpectraAggregator& sa,	
						  const SingleSpectrumHeader* header,
						  Peak*	peakBuffer);

	// writes the header and peak list in DAT format
	// returns the number of bytes written to the buffer
	size_t writeToDatBuffer(char* buffer, const SingleSpectrumHeader* newHeader) const;

	// performs simple filtering that can be done on the PeakList item itself
	// joinAdjacentPeaks
	// filterWeakPeaks
	// normalizePeakIntensities
	void initializePeakList(const Config* config, bool indFilterSpectrum);

	void computeLogIntensities(vector<float>& logIntensities);

	void createIndexArray(vector<int>& indexArray) const;
	
	void calculatePeakRanks(vector<int>& peakRanks) const;

	void calculateLogLocalRanks(mass_t windowSize, vector<float>& logLocalRanks) const;

	void calculateIsotopicLevels(mass_t tolerance, vector<float>& isotopicLevels) const;

	void selectStrongPeakIndexes(int numStrongPeaksPerLocalWindow,
									   const vector<float>& logLocalRanks,
									   const vector<float>& isotopicLevels,
									   vector<int>& strongIndexes) const;

	void calculateLogRandomProbabilities(vector<float>& logIntensities,
										 vector<float>& logRandomProbabilities) const;

	bool selectStrongPeakIndexesForPmcsqs(const vector<float>& isotopicLevels, 
										  vector<bool>& strongPeakIndicators) const;

	void markAllPossibleIsotopicPeaks(mass_t tolerance, 
									  vector<bool>& isotopicIndicators) const;

	

	PeakRange  findPeaksInRange(mass_t minRange, mass_t maxRange) const;
	int		findPeakWithMaxIntensity(mass_t peakMass, mass_t tolerance) const;
	void	findIsotopicEnvelope(int p_idx, vector<float>& iso_intens, mass_t iso_tolerance, int charge) const;

	void	printPeaks() const;
	void	printSomePeaks() const;
	bool    sanityCheck() const;
	bool	fullSanityCheck(size_t clusterSize = 0) const;

	void adjustPeakCounts(size_t clusterSize = 1);

	// outputs an mgf spectrum BEGIN IONS ... END IONS to a buffer
	// returns the number of bytes written (assumes there is enough room)
	size_t	outputToMgfBuffer(char* buffer, const SingleSpectrumHeader* newHeader = 0) const;

	void	clear();

	// an overside that ignores the constness of header_
	void setHeaderFileIndexInList(int fileIndex)
	{
		SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header_);
		nonConstHeader->setSpectraFileIndexInList(fileIndex);
	}

	// filters the peaks in the spectrum using a sliding window
	// gives the option to use a different peaks density (number peaks per 200 Da window)
	// than is designated in the config.
	void filterWeakPeaks(const Config* config, mass_t pmWith19, 
						 int peakDensity =0, bool removeGloballyWeak = false); 


protected:
	const SingleSpectrumHeader* header_;
	const Config* config_;

	int numPeaks_;	// the number of peaks stored in peaks_
	Peak *peaks_;	// this is the pointer that should be accesed by the rest of the program
	int localAllocationSize_;    // holds the maximal number of peaks that can be stored locally
								 // should be 0 if the peaks are not written to a local allocation.
								 // If the localAllocationSize is greater than 0, it is assumed that
								 // peaks_ is a pointer that was returned by new, so upon destruction
								 // of this object, the memory is deleted.

	intensity_t totalPeakIntensity_; // the sum of the peak intensities (before normalization)
	intensity_t totalNormalizedPeakIntensity_;

	bool heuristicMarkStrogenstPeaks(mass_t windowSize,
								 int numPeaksPerWindow,
								 vector<bool>& strongPeakIndicators) const;

	
	void joinAdjacentPeaks(mass_t tolerance, bool indAddCounts = false);  
	
	 
	
	void normalizePeakIntensities();

	// an overide that ignores the constness of header_
	void setHeaderFirstPeakMass(mass_t m)
	{
		SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header_);
		nonConstHeader->setFirstPeakMass(m);
	}

	

	

};




#endif

