#include "Cluster.h"
#include "MsClusterAuxfuns.h"
#include "../PepNovo/PepNovo_auxfun.h"
#include "../PepNovo/SpectraList.h"
#include "../PepNovo/AnnotatedSpectrum.h"

mass_t Cluster::precursorPPM_ = MAX_FLOAT;
mass_t Cluster::peakIndexTolerance_  = 0.0;
mass_t Cluster::peakTolerance_ = 0.0;
mass_t Cluster::isoTolerance_  = 0.0;
int	   Cluster::peakIndexToleranceAsInt_ = 0;
float  Cluster::rtTolerance_ = 5.0;
float  Cluster::scale_factor_ = 0.0;


void Cluster::setRtTolerance(float rtTol)
{
    rtTolerance_ = rtTol;
}

void Cluster::setScaleFactor(float factor)
{
    scale_factor_ = factor;
}

void Cluster::setTolerances(mass_t fragmentTolerance, float ppms)
{
    precursorPPM_		= ppms;
    peakTolerance_      = fragmentTolerance;
    isoTolerance_       = (Cluster::peakTolerance_ < 0.1 ? Cluster::peakTolerance_ : 0.1 + (Cluster::peakTolerance_-0.1)*0.5);
    peakIndexTolerance_ = Cluster::isoTolerance_ * 0.5;
    peakIndexToleranceAsInt_ = convertMassToInt(peakIndexTolerance_);
}


Cluster::~Cluster()
{
    if (singletonIdxVector_)
        singletonIdxVector_->relinquish();
    if (singletonDistancePeaks_)
        delete singletonDistancePeaks_;
}


void Cluster::copyWithoutSingletonIdxVector(const Cluster& other)
{
    indInPlay_			   = other.indInPlay_;
    clusterCharge_		   = other.clusterCharge_;
    clusterMOverZ_		   = other.clusterMOverZ_;
    totalNonPrecursorIntensity_ = other.totalNonPrecursorIntensity_;
    clusterIdx_			   = other.clusterIdx_;
    assignedClusterIdx_    = other.assignedClusterIdx_;
    clusterSize_		   = other.clusterSize_;
    distancePeaks_		   = other.distancePeaks_;

    // NP3 GOT retention time RT mean computation
    clusterTotalPrecursorIntensity_ = other.clusterTotalPrecursorIntensity_;
    clusterTotalRT_ = other.clusterTotalRT_;
    clusterRT_ 		= other.clusterRT_;
    clusterTotalRTMin_ = other.clusterTotalRTMin_;
    clusterTotalRTMax_ = other.clusterTotalRTMax_;
    clusterRTMin_	= other.clusterRTMin_;
    clusterRTMax_	= other.clusterRTMax_;

    singletonIdxVector_	   = NULL; // no singleton information

    copy(other.topPeakIdxs_, other.topPeakIdxs_ + NUM_PEAKS_FOR_HEURISTIC, topPeakIdxs_);

    lowestSimilarityPos_   = other.lowestSimilarityPos_;

    copy(other.bestSimilarityValues_, other.bestSimilarityValues_ + NUM_TOP_SIMILARITIES_TO_SAVE, bestSimilarityValues_);
    copy(other.bestSimilarityClusterIdxs_, other.bestSimilarityClusterIdxs_ + NUM_TOP_SIMILARITIES_TO_SAVE, bestSimilarityClusterIdxs_);
}



bool Cluster::createNewCluster(clusterIdx_t idx,
                               const SingleSpectrumHeader* header,
                               const Peak* peaks,
                               int   numPeaks)
{
    indInPlay_	   = true;

    clusterIdx_ = idx;
    assignedClusterIdx_ = idx;

    header_ = header;
    peaks_  = const_cast<Peak*>(peaks);

    minSimilarityToJoin_ = 2.0;
    numSimilarityPairs_ = 0;

    // should always be 0 over here
    assert(! singletonDistancePeaks_);

    // NP3 GOT retention time RT weighted mean computation
    clusterTotalPrecursorIntensity_ = header_->getPrecursorIntensity();
    clusterTotalRT_ = header_->getRetentionTime() * clusterTotalPrecursorIntensity_;
    clusterRT_ 		= clusterTotalRT_ / clusterTotalPrecursorIntensity_;
    // NP3 GOT rt width mean computation
    if ( header_->getRetentionTimeMin() == 0.0 && header_->getRetentionTimeMax() == 1000000.0) {
        // keep baseline rt range
        clusterTotalRTMin_ = header_->getRetentionTimeMin();
        clusterTotalRTMax_ = header_->getRetentionTimeMax();
        clusterRTMin_ = clusterTotalRTMin_;
        clusterRTMax_ = clusterTotalRTMax_;
    } else {
        // compute weighted average
        clusterTotalRTMin_ = header_->getRetentionTimeMin() * clusterTotalPrecursorIntensity_;
        clusterTotalRTMax_ = header_->getRetentionTimeMax() * clusterTotalPrecursorIntensity_;
        clusterRTMin_ = clusterTotalRTMin_ / clusterTotalPrecursorIntensity_;
        clusterRTMax_ = clusterTotalRTMax_ / clusterTotalPrecursorIntensity_;
    }

    /* \todo fix first peak mass issues. The problem is that the first peak
        mass should not be re-written after filtering otherwise the peak in the
        header will not be the same as the peak read from the file
    assert( header_->getFirstPeakMass() == peaks_[0].mass ); */

    clusterCharge_ = header_->getCharge();
    clusterMOverZ_ = header_->getMOverZ();
    numPeaks_	   = numPeaks;

    clusterSize_   = (header->getClusterSize() < 1  ? 1 : header->getClusterSize());

    if (singletonIdxVector_)
    {
        singletonIdxVector_->relinquish();
        singletonIdxVector_ = NULL;
    }

    if (totalNonPrecursorIntensity_ <= 0.0 || totalNormalizedPeakIntensity_ <= 0.0)
        normalizePeakIntensities();

    lowestSimilarityPos_=0;
    for (size_t i=0; i<NUM_TOP_SIMILARITIES_TO_SAVE; i++)
    {
        bestSimilarityValues_[i]=MIN_FLOAT;
        bestSimilarityClusterIdxs_[i]=MAX_CLUSTER_IDX;
    }

    if (! selectDistancePeaksAndTopIdxs())
        return false;

    setAdjustedIntensities();

    return true;
}


struct MassTolerancePair {
    MassTolerancePair() : mass(0), tolerance(0) {}
    MassTolerancePair(mass_t m, mass_t t) : mass(m), tolerance(t) {}
    bool operator< (const MassTolerancePair& rhs) const
    {
        return (mass<rhs.mass);
    }
    mass_t mass;
    mass_t tolerance;
};


// this function selects a number of peaks that are usd for computing dot products (tries to
// avoid peaks areound the precursor mass and isotopic peaks
// chooses the idxs of the strongest n=NUM_PEAKS_FOR_HEURISTIC peaks
// (the bin indexes that are converted from the peak masses is what actually gets stored)
bool Cluster::selectDistancePeaksAndTopIdxs()
{
    assert(peakIndexTolerance_>0.0);

    // get list of masses and margins for M+H, M+H-H2O, (M+2H)/2, (M+2H-H2O)/2, etc.
    static vector<MassTolerancePair> invalidMasses;
    invalidMasses.clear();
    const mass_t maxPeakMass = peaks_[numPeaks_-1].mass;
    if (clusterCharge_ > 0)
    {
        const mass_t mPlusH = (clusterMOverZ_ * clusterCharge_ - (clusterCharge_-1)*MASS_PROTON);
        const mass_t mPlusHMinusH2O = mPlusH - MASS_H2O;
        for (short c=clusterCharge_; c>0; c--)
        {
            const mass_t multVal = 1.0 / static_cast<mass_t>(c);
            const mass_t tolerance = (MASS_ISO + peakTolerance_) * multVal ;
            mass_t mz		 =  (mPlusH + (c-1)*MASS_PROTON) * multVal;
            mass_t mzH2O	 =  (mPlusHMinusH2O + (c-1)*MASS_PROTON) * multVal;

            if (mzH2O>maxPeakMass)
                break;

            invalidMasses.push_back( MassTolerancePair(mzH2O, tolerance) );
            invalidMasses.push_back( MassTolerancePair(mz, tolerance) );
        }
    }
    else
    {
        for (short precursorCharge=1; precursorCharge<=3; precursorCharge++)
        {
            const mass_t mPlusH = (clusterMOverZ_ * precursorCharge - (precursorCharge-1)*MASS_PROTON);
            const mass_t mPlusHMinusH2O = mPlusH - MASS_H2O;
            for (short c=clusterCharge_; c>0; c--)
            {
                const mass_t multVal = 1.0 / static_cast<mass_t>(c);
                const mass_t tolerance = (MASS_ISO + peakTolerance_) * multVal ;
                mass_t mz		 =  (mPlusH + (c-1)*MASS_PROTON) * multVal;
                mass_t mzH2O	 =  (mPlusHMinusH2O + (c-1)*MASS_PROTON) * multVal;

                if (mzH2O>maxPeakMass)
                    break;

                invalidMasses.push_back( MassTolerancePair(mzH2O, tolerance) );
                invalidMasses.push_back( MassTolerancePair(mz, tolerance) );
            }
        }
    }
    sort(invalidMasses.begin(), invalidMasses.end());


    const mass_t MAX_ISO_CHARGE1_MASS = MASS_ISO + isoTolerance_;
    const mass_t MIN_ISO_CHARGE1_MASS = MASS_ISO - isoTolerance_;
    const mass_t MAX_ISO_CHARGE2_MASS = MASS_ISO*0.5 + isoTolerance_;
    const mass_t MIN_ISO_CHARGE2_MASS = MASS_ISO*0.5 - isoTolerance_;
    static vector<IndexIntensityPair> pairs;
    pairs.clear();
    totalNonPrecursorIntensity_ = 0.0;
    int idxForInvalid=-1;
    mass_t maxInvalidMass = 0.0;
    mass_t minInvalidMass = 0.0;
    for (size_t i=0; i<numPeaks_; i++)
    {
        while (idxForInvalid < invalidMasses.size() &&
               peaks_[i].mass > maxInvalidMass)
        {
            idxForInvalid++;
            if (idxForInvalid == invalidMasses.size())
            {
                maxInvalidMass = 0.0;
                minInvalidMass = 0.0;
                break;
            }
            else
            {
                maxInvalidMass = invalidMasses[idxForInvalid].mass + invalidMasses[idxForInvalid].tolerance;
                minInvalidMass = invalidMasses[idxForInvalid].mass - invalidMasses[idxForInvalid].tolerance;
            }
        }
        if (peaks_[i].mass < maxInvalidMass &&
            peaks_[i].mass > minInvalidMass)
            continue;

        // don't use very large peaks for distance peaks! -> peaks_[i].mass > 10000.0  removed
        // GOT NP3 nor peaks with mass greater then the cluster mz
        if (peaks_[i].mass > clusterMOverZ_ + MAX_ISO_CHARGE1_MASS + peakTolerance_)
            continue;

        totalNonPrecursorIntensity_+=peaks_[i].intensity;

        // check for iso ratio don't choose peak if it is less than 2/3 of any possible isotopic peak
        // variant (+-1 , +- 0.5)
        bool inValid = false;
        for (size_t j=i+1; j<numPeaks_; j++)
        {
            if (peaks_[j].mass > peaks_[i].mass + MAX_ISO_CHARGE1_MASS)
                break;

            if ( peaks_[j].mass > peaks_[i].mass + MIN_ISO_CHARGE1_MASS ||
                 (peaks_[j].mass < peaks_[i].mass + MAX_ISO_CHARGE2_MASS &&
                  peaks_[j].mass > peaks_[i].mass + MIN_ISO_CHARGE2_MASS) )
            {
                if (peaks_[j].intensity * 0.66 > peaks_[i].intensity)
                {
                    inValid = true;
                    break;
                }
            }
        }
        if (inValid)
            continue;

        if (i>0)
        {
            for (size_t j=i-1; j>0; j--)
            {
                if (peaks_[j].mass < peaks_[i].mass - MAX_ISO_CHARGE1_MASS)
                    break;

                if ( peaks_[j].mass < peaks_[i].mass - MIN_ISO_CHARGE1_MASS ||
                     (peaks_[j].mass > peaks_[i].mass - MAX_ISO_CHARGE2_MASS &&
                      peaks_[j].mass < peaks_[i].mass - MIN_ISO_CHARGE2_MASS) )
                {
                    if (peaks_[j].intensity * 0.66 > peaks_[i].intensity)
                    {
                        inValid = true;
                        break;
                    }
                }
            }

        }
        if (inValid)
            continue;

        pairs.push_back(IndexIntensityPair(i,peaks_[i].intensity));
    }

    sort(pairs.begin(), pairs.end());
    // take first NUM_PEAKS_FOR_HEURISTIC and compute their idxs
    const size_t numTopPeaksToTake = (pairs.size()>= NUM_PEAKS_FOR_HEURISTIC ?
                                      NUM_PEAKS_FOR_HEURISTIC : pairs.size());
    size_t i;
    for (i=0; i<numTopPeaksToTake; i++)
        topPeakIdxs_[i]=computeMzIndex( peaks_[pairs[i].index].mass, peakIndexTolerance_ , 0.0);

    for ( ; i<NUM_PEAKS_FOR_HEURISTIC; i++)
        topPeakIdxs_[i]=0;

    if (totalPeakIntensity_ <= 0.0)
    {
        totalPeakIntensity_ = 0.0;
        for (size_t i=0; i<numPeaks_; i++)
            totalPeakIntensity_ += peaks_[i].intensity;
    }
    const intensity_t neededIntensity = 0.95 * totalPeakIntensity_;
    intensity_t sumIntensity=0.0;
    size_t peakIdx;
    for (peakIdx=0; peakIdx<numPeaks_ && sumIntensity<neededIntensity; peakIdx++)
        sumIntensity += peaks_[peakIdx].intensity;

    size_t maxDistancePeaks = static_cast<size_t>(10.0 + (peaks_[peakIdx-1].mass * 0.0004) * MAX_NUM_PEAKS_FOR_DISTANCE);
    if (maxDistancePeaks > MAX_NUM_PEAKS_FOR_DISTANCE)
        maxDistancePeaks = MAX_NUM_PEAKS_FOR_DISTANCE;

    // get the masses of the maxDistancePeaks_ peak masses (sort and store them in distancePeaks_)
    if (pairs.size() < maxDistancePeaks)
        maxDistancePeaks = pairs.size();

    distancePeaks_.numPeaks = maxDistancePeaks;
    static vector<DistancePeak> tmp;
    tmp.resize(distancePeaks_.numPeaks);
    for (size_t i=0; i<distancePeaks_.numPeaks; i++)
    {
        const Peak& peak		 = peaks_[pairs[i].index];
        tmp[i].massAsInt		 = convertMassToInt(peak.mass);
        tmp[i].intensity		 = peak.intensity;
        tmp[i].adjustedIntensity = 0.0;
    }

    sort(tmp.begin(),tmp.end());
    memcpy(distancePeaks_.peaks, &tmp[0], sizeof(DistancePeak)*distancePeaks_.numPeaks);

    return (checkDistancePeaksOk());
}


/**************************************************************************
This is a debug function that checks that the distance peaks have vlaid values.
***************************************************************************/
bool Cluster::checkDistancePeaksOk() const
{
    const size_t n=distancePeaks_.numPeaks;
    // NP3 GOT change n < 2 to 1
    if (n<1 || n>MAX_NUM_PEAKS_FOR_DISTANCE)
    {
        cout << "Problem with distance peaks: num peaks = " << n << " (" << numPeaks_ << ")" << endl;
        return false;
    }

    const int maxMassAsInt = convertMassToInt(20000.0);
    size_t i;
    for (i=0; i<n; i++)
        if (distancePeaks_.peaks[i].massAsInt < 0 ||
            distancePeaks_.peaks[i].massAsInt > maxMassAsInt ||
            (i>0 && distancePeaks_.peaks[i].massAsInt == distancePeaks_.peaks[i-1].massAsInt))
            break;

    if (i<n)
    {
        cout << "Cluster idx : " << clusterIdx_ << endl;
        cout << "NUM DP: " << n << endl;
        for (size_t i=0; i<n; i++)
        {

            cout << "DP " << i << fixed << setprecision(3) << "\t" << convertIntToMass(distancePeaks_.peaks[i].massAsInt)
                 << "\t" <<  distancePeaks_.peaks[i].adjustedIntensity << endl;
        }

        cout << "Main header : ";
        header_->printStats();
        cout << "Title: " << header_->getTitle() << endl;
        printPeaks();
        cout << endl;
        error("BAD DISTANCE PEAKS!");
    }
    return true;
}

/**************************************************************************
This function is used to set the adjusted peak intensities which are
the values used to compute the vector dot product distance between clusters.
Function assumes that the mass and intensity values are stored in distancePeaks_
(the number of peaks is in numPeaks_). The function changes the adjusted intensity 
values in  distancePeaks_ .
Also computes the squared sum of adjusted intensity
***************************************************************************/
void Cluster::setAdjustedIntensities()
{
    assert(totalNonPrecursorIntensity_>0.0);
    distancePeaks_.sumSqrAdjustedIntensity = 0.0;
    const intensity_t multVal = 1000.0 / totalNonPrecursorIntensity_;
    float t;

    for (size_t i=0; i<distancePeaks_.numPeaks; i++)
    {
        // NP3 got scale sqrt intensity
        if (scale_factor_ == 0.0) {
            t = 1.0 + log(1.0 + multVal*distancePeaks_.peaks[i].intensity);
        } else {
            t = pow(multVal*distancePeaks_.peaks[i].intensity, scale_factor_);
        }

        distancePeaks_.peaks[i].adjustedIntensity = t;
        distancePeaks_.sumSqrAdjustedIntensity += (t*t);
    }
}



/*****************************************************************************
Sets the max possible field for all peaks.
This is done by examining the min/max fields in the supplied clusters
which are assumed to be the peaklists that make up the current cluster
******************************************************************************/
void Cluster::adjustMaxPossiblePeakCounts(const vector<MassCount>& minMassCounts,
                                          const vector<MassCount>& maxMassCounts)
{
    const mass_t lastMinMass  = minMassCounts[minMassCounts.size()-1].mass;
    const mass_t firstMaxMass = maxMassCounts[0].mass;
    size_t minIdx =0;
    size_t maxIdx =0;
    size_t i;
    for (i=0; i<numPeaks_; i++)
    {
        while (minIdx<minMassCounts.size() &&
               peaks_[i].mass > minMassCounts[minIdx].mass)
        {
            minIdx++;
        }
        if (minIdx >= minMassCounts.size()-1)
            break;
        peaks_[i].maxPossible = minMassCounts[minIdx].numSpectra;

        // add this fix because the first peaks can all get joined to the first peak mass
        // (see joinAdjacentPeaks). In this case the count can be higher than the assigned maxPossible
        if (peaks_[i].count > peaks_[i].maxPossible)
            peaks_[i].count = peaks_[i].maxPossible;
    }

    for ( ; i<numPeaks_; i++)
    {
        if (peaks_[i].mass>=firstMaxMass)
            break;
        peaks_[i].maxPossible = clusterSize_;
        // add this fix because the first peaks can all get joined to the first peak mass
        // (see joinAdjacentPeaks). In this case the count can be higher than the assigned maxPossible
        if (peaks_[i].count > clusterSize_)
            peaks_[i].count = clusterSize_;
    }

    for ( ; i<numPeaks_; i++)
    {
        while (peaks_[i].mass<maxMassCounts[maxIdx].mass)
        {
            maxIdx++;
        }
        assert(maxIdx < maxMassCounts.size());
        peaks_[i].maxPossible = maxMassCounts[maxIdx].numSpectra;

        // add this fix because the first peaks can all get joined to the first peak mass
        // (see joinAdjacentPeaks). In this case the count can be higher than the assigned maxPossible
        if (peaks_[i].count > peaks_[i].maxPossible)
            peaks_[i].count = peaks_[i].maxPossible;
    }
}

void Cluster::detailedPrint() const
{
    cout << "IDX: " << clusterIdx_ << " ";
    cout << "Cluster size: " << clusterSize_ << "  Num singltons ";
    if (! singletonIdxVector_)
    {
        cout << "0" << endl;
    }
    else
    {
        cout << singletonIdxVector_->getSize() << ":";
        for (size_t i=0; i<singletonIdxVector_->getSize(); i++)
            cout << " " << singletonIdxVector_->getElementsPtr()[i];
        cout << endl;
    }
    cout << "vec: " << singletonIdxVector_ << endl;
}




void readAnnotatedSpectraIntoClusters(const string& list,
                                      const Config* config,
                                      SpectraAggregator& sa,
                                      vector<Cluster>& clusters,
                                      map<Annotation, vector<size_t> >& annotationMap)
{
    cout << "Reading spectra from: " << list << endl;
    if (getFileExtensionType(list.c_str()) == IFT_TXT)
    {
        sa.initializeFromTextFile(list.c_str(), config);
    }
    else
        sa.initializeFromSpectraFilePath(list.c_str(), config);

    SpectraList sl(sa);
    sl.selectAllAggregatorHeaders();
    cout << "... done." << endl;

    clusters.clear();
    clusters.resize(sl.getNumHeaders());
    annotationMap.clear();

    cout << "Found " << sl.getNumHeaders() << " spectra." << endl;
    cout << "Reading spectra... ";
    size_t tenth = sl.getNumHeaders()/10;
    // create singleton clusters for all spectra
    for (size_t i=0; i<sl.getNumHeaders(); i++)
    {
        const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
        if (header->getPeptideStr().length()<2)
            continue;

        PeakList pl;
        if (pl.readPeaksToLocalAllocation(sa, header) == 0)
            continue;
        pl.initializePeakList(config, true);
        assert(pl.sanityCheck());


        clusters[i].createNewCluster(i, header, pl.getPeaks(), pl.getNumPeaks());
        Annotation ann;
        ann.peptideStr = header->getPeptideStr();
        ann.charge	   = header->getCharge();

        assert(ann.peptideStr.length()>0 && ann.charge>0);

        annotationMap[ann].push_back(i);

        if (i>0 && i % tenth == 0)
        {
            cout << i/tenth << " ";
            cout.flush();
        }
    }
    cout << " ...done." << endl;
}

void Cluster::backupDistancePeaks()
{
    if (! singletonDistancePeaks_)
        singletonDistancePeaks_ = new DistancePeakList;
    *singletonDistancePeaks_ = distancePeaks_;
}










