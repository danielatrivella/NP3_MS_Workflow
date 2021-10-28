#ifndef __BUFFEREDOUTPUTFILE_H__
#define __BUFFEREDOUTPUTFILE_H__

#include "includes.h"

class BufferedOutputFile {
public:
	BufferedOutputFile() : indOpen_(0), fileName_(std::string()), stream_(0), buffer_(0), nextPos_(0),
		bufferSize_(0) {};

	bool open(const string& name, size_t bufferSize, bool binaryMode=false);
	bool close();
	void flush();
	void setBufferSize(size_t t);
	size_t writeToBuffer(const char* src, size_t n);

	bool isOpen() const { return (indOpen_); }
	size_t getNextPos() const { return nextPos_; }
	size_t getBufferSize() const {return bufferSize_; }

	void writeBufferToStream(ostream& os = cout) const;

private:
	bool   indOpen_;
	string fileName_;
	FILE* stream_;
	char* buffer_;
	size_t nextPos_;
	size_t bufferSize_;
};




#endif




