#include "BufferedOutputFile.h"
#include <string.h>


bool BufferedOutputFile::open(const string& name, size_t bufferSize, bool binaryMode)
{
	close();

	fileName_ = name;
	setBufferSize(bufferSize);

	string mode = (binaryMode ? "wb" : "w");
	stream_ = fopen(fileName_.c_str(), mode.c_str());
	if (! stream_)
	{
		cout << "Error: could not open file for writing: " << fileName_ << endl;
		exit(1);
	}

	indOpen_ = true;
	return true;
}

bool BufferedOutputFile::close()
{
	if (nextPos_>0)
		flush();

	if (stream_ && indOpen_)
		fclose(stream_);
	stream_ = NULL;

	if (buffer_)
		delete [] buffer_;
	buffer_ = NULL;
	
	indOpen_ = false;
	return true;
}

void BufferedOutputFile::flush()
{
	if (indOpen_ && nextPos_>0)
	{
		if (fwrite(buffer_, 1, nextPos_, stream_) != nextPos_)
		{
			cout << "Warning: couldn't write all data from buffer to file: " << fileName_ << endl;
		}
	}
		
	nextPos_ = 0;
}

void BufferedOutputFile::setBufferSize(size_t t)
{
	if (buffer_ && bufferSize_ == t)
		return;

	char* newBuffer;
	bufferSize_ = t;
	newBuffer = new char[bufferSize_];

	if (nextPos_>0)
	{
		size_t numToCopy = nextPos_;
		if (numToCopy>bufferSize_)
			numToCopy = bufferSize_;
		memcpy(newBuffer, buffer_, numToCopy);
		nextPos_ = numToCopy;
	}

	if (buffer_)
		delete [] buffer_;

	buffer_ = newBuffer;
}


size_t BufferedOutputFile::writeToBuffer(const char* src, size_t n)
{
	if (! indOpen_)
	{
		cout << "writing to closed BufferedOutputFile!" << endl;
		exit(1);
	}

	if (n>bufferSize_)
		setBufferSize(n+1);

	if (nextPos_ + n >= bufferSize_)
		flush();

	memcpy(buffer_ + nextPos_, src, n);
	nextPos_ += n;
	return n;
}

void BufferedOutputFile::writeBufferToStream(ostream& os) const
{
	for (size_t i=0; i<nextPos_; i++)
		os << buffer_[i];
}

