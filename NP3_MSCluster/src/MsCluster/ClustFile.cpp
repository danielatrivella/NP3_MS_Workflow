#include "ClustFile.h"
#include "../Common/auxfun.h"


ClustFile::~ClustFile() {
	
	if (buffer_)
		delete [] buffer_;
}

bool ClustFile::read(const char* path, int datasetIdx, int fileIdx)
{
	path_       = path;
	datasetIdx_ = datasetIdx;
	fileIdx_    = fileIdx;

	fileSize_ = getFileSize(path);

	if (fileSize_ == 0)
		return false;

	if (bufferSize_ <= fileSize_)
	{
		delete [] buffer_;
	//	cout << endl << "BUFF: " << bufferSize_ << " --> ";
		bufferSize_ = 0;

		buffer_ = new char[fileSize_+1024];
		if (! buffer_)
			error("Could no allocate buffer to read clust file: ",path);
		bufferSize_ = fileSize_+1024;
	//	cout << bufferSize_ << endl;
	}

	FILE* stream=fopen(path,"rb");
	if (! stream)
		error("Could not open clust file: ",path);

	fread(buffer_,1,fileSize_,stream);
	buffer_[fileSize_]='\0';
	buffer_[fileSize_+1]='\0';

	fclose(stream);
	currentPos_ = buffer_;
	bufferEnd_  = buffer_ + fileSize_;
	return true;
}

char* ClustFile::getNextEntryPointer(string* titleString)
{
	
	while (currentPos_< bufferEnd_ && (*currentPos_ == '\n' || *currentPos_ =='#' || *currentPos_=='\t' || *currentPos_ == ' '))
		currentPos_++;

	if (currentPos_>=bufferEnd_)
		return 0;
// for some reason on linux this doesn't parse correctly
#if !(defined(WIN32) || defined(WIN64))
	if (currentPos_>buffer_)
		currentPos_--;
#endif


	char *entry = currentPos_;

	// parse line to know where to advance
	static char titleBuffer[128];
	int	   numLines = 0;
	sscanf(currentPos_,"%s\t%d",titleBuffer,&numLines);
	if (numLines<1)
		error("Could not parse clust file correctly: ",path_.c_str());

	if (titleString)
		*titleString = std::string(titleBuffer);

	for (int i=0; i<=numLines; i++)
	{
		currentPos_ = strchr(currentPos_, '\n');
		if (! currentPos_ && i<numLines)
			error("did not parse file correctly: ",path_.c_str(),titleBuffer);

		while (currentPos_<bufferEnd_ && *currentPos_++ == '\n')
			;
	}
		
	return entry;
}


//HEK2_0_0.9160	2	462.827	0
//0	0	14077	462.839	0.81	3.139e-005	0
//0	1	14740	462.827	0.95	0.000e+000	0


bool ClustFile::parseEntry(char* p, ClusterClustEntry& clust) const
{
	char* endLinePtr = strchr(p,'\n');
	istringstream iss(std::string(p,endLinePtr-p));
	clust.clusterSize = 0;
	iss >> clust.title >> clust.clusterSize >> clust.mz >> clust.charge >> clust.peptide;
	if (clust.clusterSize<1)
		return false;

	clust.entries.resize(clust.clusterSize);
	p=endLinePtr+1;

	for (size_t i=0; i<clust.clusterSize; i++)
	{
		endLinePtr = strchr(p,'\n');
		istringstream iss(std::string(p,endLinePtr-p));
		SingletonClustEntry& ent = clust.entries[i];
		iss >> ent.datasetIdx >> ent.fileIdx >> ent.scanNumber >> ent.mz >> ent.similarityToConsensus >> ent.pValue >> ent.charge;
		if (ent.charge>0)
			iss >> ent.peptide;

		p = endLinePtr +1;
	}

	return true;
}

