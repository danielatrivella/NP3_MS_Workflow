#ifndef __ANNOTATIONFILE_H__
#define __ANNOTATIONFILE_H__

#include "PepNovo_includes.h"
#include "Config.h"

struct FileAndScan {

	FileAndScan() : fileIdx(MAX_SIZE_T), scan(MAX_SIZE_T) {}

	bool operator< (const FileAndScan& rhs) const
	{
		return ( fileIdx < rhs.fileIdx ||
			    (fileIdx == rhs.fileIdx && scan < rhs.scan) );
	}
	size_t fileIdx;
	size_t scan;
};

struct Annotation {

	Annotation() : peptideStr(""), charge(0) {}

	bool operator< (const Annotation& rhs) const
	{
		return (peptideStr<rhs.peptideStr ||
			    (peptideStr == rhs.peptideStr && charge < rhs.charge));
	}

	string peptideStr;
	short  charge;
};

class AnnotationFile {
public:
	void readAnnotationFile(const char* annPath,  bool indClear=true);

	void createAnnotatedMgfs(const char* fileList, 
							 const Config* config,
							 string& name,
							 bool onlyAnnotated = true,
							 size_t maxNumPerFile = 100000) const;

private:

	map<FileAndScan, Annotation> annotations_;
};





#endif


