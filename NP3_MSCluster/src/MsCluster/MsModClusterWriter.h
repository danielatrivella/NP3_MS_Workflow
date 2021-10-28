#ifndef __MSMODCLUSTERWRITER_H__
#define __MSMODCLUSTERWRITER_H__

#include "../Common/BufferedOutputFile.h"
#include "MsClusterIncludes.h"

class ModClusterWriter {
public:

	ModClusterWriter() : outPath_(std::string()), fileIdx_(MAX_SIZE_T), spectraCounter_(0), numSpectraPerFile_(0), indInUse_(0) {}

	void init(const string& path, size_t numSpectra = 20000) {
		outPath_ = path;
		fileIdx_ = MAX_SIZE_T;
		spectraCounter_ = 0;
		numSpectraPerFile_ = numSpectra;
		indInUse_		   = true;
		numClusters_	   = 0;
	}

	void writeToBuffer(const char* buffer, size_t len)
	{
		if (! outFile_.isOpen())
		{
			fileIdx_ = (fileIdx_ == MAX_SIZE_T ? 0 : fileIdx_ +1);
			spectraCounter_ = 0;
			ostringstream oss;
			oss << outPath_ << "_" << fileIdx_ << ".mgf";
			outFile_.open(oss.str(), 1000000);
		}
		outFile_.writeToBuffer(buffer, len);
		numClusters_++;
		if (++spectraCounter_ >= numSpectraPerFile_)
			outFile_.close();
	}

	void close()
	{
		if (indInUse_)
			outFile_.close();
		indInUse_ = false;
	}

	size_t getFileIdx() const { return fileIdx_; }
	size_t getNumClusters() const { return numClusters_; }

private:

	string			   outPath_;
	size_t			   fileIdx_;
	size_t			   spectraCounter_;
	size_t			   numSpectraPerFile_;

	bool			   indInUse_;

	BufferedOutputFile outFile_;

	size_t			   numClusters_;
};

#endif

