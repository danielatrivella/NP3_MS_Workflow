#include "MsArchiveUpdater.h"
#include "../PepNovo/BasicDataStructs.h"
#include "../Common/auxfun.h"


void MsArchiveUpdater::scanUpdateFile(const string& filePath, const Config* config)
{
	updates_.clear();
	countsInUpdateFile_.clear();
	minIdxLoaded_=0;
	maxIdxLoaded_=0;

	updateFilePath_ = filePath;
	config_ = config;

	char buffer[256];
	ifstream ifs(filePath.c_str());
	if (! ifs.good())
		error("Could not open archive update file for reading: ",filePath.c_str());

	countsInUpdateFile_.resize(10000,0);
	size_t badLines=0;
	size_t goodLines=0;

	while (ifs.getline(buffer,256))
	{
		if (buffer[0]=='#')
			continue;

		if (ifs.gcount()<5)
			continue;

		string title=std::string(), peptide=std::string();
		int charge =0;
		istringstream iss(buffer);
		iss >> title >> charge >> peptide;
		if (peptide.length()>0 && charge>0 && charge<25)
		{
			Peptide pep;
			pep.parseFromString(config, peptide);
			mass_t mz = (pep.get_mass_with_19() + (charge-1)*MASS_PROTON)/static_cast<mass_t>(charge);
			if (mz<10.0 || mz>=10000.0)
			{
				badLines++;
				continue;
			}
			countsInUpdateFile_[static_cast<size_t>(mz)]++;
			goodLines++;
		}
		else
			badLines++;
	}
	ifs.close();
	if (badLines>0)
		cout << "Warning: found " << badLines << " bad lines in the update file!" << endl;
}


void MsArchiveUpdater::loadUpdates(mass_t minMz, mass_t maxMz, size_t maxNumUpdatesToLoad)
{
	const size_t minIdx = (minMz<1.0 ? 0 : static_cast<size_t>(minMz-1.0));
	size_t maxIdx = static_cast<size_t>(maxMz);
	
	assert (minIdx<=maxIdx);
	if (minIdx >= minIdxLoaded_ && maxIdx <= maxIdxLoaded_)
		return;

	updates_.clear();
	size_t numUpdates=0;
	for (size_t i=minIdx; i<=maxIdx; i++)
		numUpdates += countsInUpdateFile_[i];

	while (maxIdx<countsInUpdateFile_.size()-1 && numUpdates+countsInUpdateFile_[maxIdx+1]<maxNumUpdatesToLoad)
	{
		maxIdx++;
		numUpdates += countsInUpdateFile_[maxIdx];
	}
	
	// read files and load all updates with the idxs in the range
	ifstream ifs(updateFilePath_.c_str());
	if (! ifs.good())
		error("Could not open update file for reading: ",updateFilePath_.c_str());

	char buffer[256];
	while (ifs.getline(buffer,256))
	{
		if (buffer[0]=='#' || ifs.gcount()<5)
			continue;

		string title=std::string(), peptide=std::string();
		int charge =0;
		istringstream iss(buffer);
		iss >> title >> charge >> peptide;
		if (peptide.length()>0 && charge>0 && charge<25)
		{
			Peptide pep;
			pep.parseFromString(config_, peptide);
			mass_t mz = (pep.get_mass_with_19() + (charge-1)*MASS_PROTON)/static_cast<mass_t>(charge);
			if (mz<10.0 || mz>=10000.0)
				continue;

			const size_t idx = static_cast<size_t>(mz);
			if (idx>=minIdx && idx <= maxIdx)
				updates_[title]=UpdateEntry(charge,peptide);
		}
	}
	ifs.close();
}



