#ifndef __METALIST_H__
#define __METALIST_H__


/*! \file MetaList.h
	\brief Holds class MetaList
*/

#include "PepNovo_includes.h"

class ScanManagerList; // fwd dclr

struct SinglePath {
	SinglePath() : path(std::string()), datasetIdx(-1), idxInList(-1) {}
	SinglePath(const string& s, int g, int f) : path(s), datasetIdx(g), idxInList(f) {}
	string path;
	int	   datasetIdx;
	int	   idxInList;	
};

/*! \class MetaList
	\biref A list of lists (where each list is considered a differnt generation)
*/
class MetaList
{
public:
	size_t readMetaList(const char* path);

	const vector<SinglePath>& getSinglePaths() const { return singlePaths_; }

	void writeLists(const char* path, const char*) const;

private:
	string metaListPath_;
	vector<SinglePath> singlePaths_;
	vector< vector<string> > dsPaths_;
};

#endif

