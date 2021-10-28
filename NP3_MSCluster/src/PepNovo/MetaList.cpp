#include "MetaList.h"
#include "../Common/auxfun.h"


size_t MetaList::readMetaList(const char *path)
{	
	// read mega list (which is a list of paths to list files)
	// generation indexes equal the position in the list (starting from 0)

	metaListPath_ = std::string(path);

	vector<string> listPaths;
	readListOfPaths(path, listPaths);

	assert(listPaths.size()>0);

	dsPaths_.resize(listPaths.size());
	for (size_t g=0; g<listPaths.size(); g++)
	{
		vector<string> genPaths;
		readListOfPaths(listPaths[g].c_str(), genPaths);

		if (genPaths.size()==0)
		{
			cout << "Error: read no paths from " << listPaths[g] << endl;
		}
		assert(genPaths.size()>0);

		for (size_t f=0; f<genPaths.size(); f++)
		{
			singlePaths_.push_back(SinglePath(genPaths[f],g,f));
			dsPaths_[g].push_back(genPaths[f]);
		}
	}

	cout << "Read meta list: " << metaListPath_ << endl;
	cout << "Found " << singlePaths_.size() << " paths from " << listPaths.size() << " lists." << endl;
	return singlePaths_.size();
}

void MetaList::writeLists(const char* path, const char* suffix) const
{
	for (size_t i=0; i<dsPaths_.size(); i++)
	{
		ostringstream oss;
		oss << path << "_" << i << suffix;
		writeListOfPaths(oss.str().c_str(), dsPaths_[i]);
	}
}

