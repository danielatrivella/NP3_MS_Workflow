#ifndef __SCANLIST_H__
#define __SCANLIST_H__

/*! \file ScanList.h
	holds pointers to the scan co-ordinates that need to be excluded.
*/

#include "PepNovo_includes.h"

/*! \file ScanExclude.h
	\brief reads an exclusion text file that has lines in the format:
	gen	fileIdx	scan#	m/z	charge
*/

struct ScanEntry {
	ScanEntry() : datasetIdx(0), fileIdx(0), scan(0), mz(0.0) {}
	ScanEntry(int g, int f, int s, mass_t m=0.0) : datasetIdx(g), fileIdx(f), scan(s), mz(m) {}

	bool operator< (const ScanEntry& rhs) const
	{
		return (datasetIdx<rhs.datasetIdx ||
			   (datasetIdx == rhs.datasetIdx && fileIdx < rhs.fileIdx) ||
			   (datasetIdx == rhs.datasetIdx && fileIdx == rhs.fileIdx && scan < rhs.scan));
	}


	int datasetIdx;
	int fileIdx;
	int scan;
	mass_t mz;
};

class ScanListManager {
public:
	/*! Reads all exclusion files listed, and keep
		the exclusions that fall within the desired mass range.

		\return number of scans read into exclusion list
	*/
	size_t initialize(const char* list, mass_t minMz=0.0, mass_t maxMz=999999.0);

	bool checkScanList(const ScanEntry& e) const
	{
		map< ScanEntry, int>::const_iterator it = scans_.find(e);
		if (it == scans_.end())
			return false;

		// test for exclusion errors in mass
		assert( it->first.mz <=0.0 || fabs(it->first.mz - e.mz)<7.0 );
		return true;
	}
	
	
	const map< ScanEntry, int>& getScans() const { return scans_; }

private:
	map< ScanEntry, int> scans_;
};


/*! \fn makeScanListFromClust
	\brief makes a scan list from all files given in the list.
	Writes the output to scanFilePath. The scans get sorted so memory size
	needs to be minded.
*/
void makeScanListFromClust(const char* list, const char* scanFilePath);



#endif

