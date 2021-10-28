#ifndef __MSARCHIVEUPDATER_H__
#define __MSARCHIVEUPDATER_H__

#include "../PepNovo/PepNovo_includes.h"

class Config; // fwd dclr

struct UpdateEntry {
	UpdateEntry() :  charge(0), peptideStr(std::string()) {}
	UpdateEntry(int c, const string& s) : charge(c), peptideStr(s) {}
	int charge;
	string peptideStr;
};


class MsArchiveUpdater {
public:

	const UpdateEntry* getUpdateEntry(const string& title) const
	{
		map<string, UpdateEntry>::const_iterator it = updates_.find(title);
		if (it == updates_.end())
			return static_cast<const UpdateEntry*>(0);
		return &(it->second);
	}

	void scanUpdateFile(const string& filePath, const Config* config);

	void loadUpdates(mass_t minMZ, mass_t maxMZ, size_t maxNumUpdatesToLoad=4000000);


private:

	const Config* config_;
	string updateFilePath_;
	size_t numLinesInUpdateFile_;
	vector<size_t> countsInUpdateFile_;

	map<string, UpdateEntry> updates_;
	size_t minIdxLoaded_, maxIdxLoaded_;

};


#endif



