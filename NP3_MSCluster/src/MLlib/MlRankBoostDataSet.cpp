#include "MlRankBoost.h"
#include "../Common/auxfun.h"
#include "../Common/BufferedOutputFile.h"


size_t MlRankBoostDataset::readRankBoostDataSet(const char* samplePath, 
												size_t maxSampleIdx)
{

	string fileName, dirPath;
	getFileNameWithoutExtension(samplePath, fileName);
	getDirFromFullPath(samplePath, dirPath);
	
	// read data file, include tag information
	ifstream ifs(samplePath);
	if (! ifs.good())
		error("Couldn't MlRankBoostDataset for reading: ", samplePath);
	
	// count number of samples
	size_t counter=0;
	while (ifs.good())
	{
		if (ifs.get() == '>')
			counter++;
	}
	ifs.close();
	ifs.open(samplePath);
	samples_.resize(counter<maxSampleIdx ? counter : maxSampleIdx);
	counter=0;
	while (! ifs.eof() && counter < maxSampleIdx)
	{
		if (samples_[counter].readFromStream(ifs))
			counter++;
	}
	if (counter<samples_.size())
		samples_.resize(counter);

	// read phi file
	phi_.clear();
	string phiPath = dirPath + "/" + fileName + ".phi";
	ifs.open(phiPath.c_str());
	if (! ifs.good())
	{
		cout << "Warning: could not find feedback file " << phiPath << endl;
	}
	else
	{
		char buffer[64];
		while (! ifs.eof())
		{
			ifs.getline(buffer,64);
			istringstream iss(buffer);
			OrderedPair op;
			iss >> op.idx0 >> op.idx1 >> op.weight;
			if (op.idx0 < maxSampleIdx && op.idx1 < maxSampleIdx)
			{
				assert(op.idx0 < op.idx1);
				if (op.weight < 0.0)
					op.weight = 1.0;

				phi_.push_back(op);
				totalPhiWeight_ += op.weight;
			}
		}
	}

	return (samples_.size());
}



void MlRankBoostDataset::writeRankBoostDataSet(const char* samplePath)
{
	string fileName, dirPath;
	getFileNameWithoutExtension(samplePath, fileName);
	getDirFromFullPath(samplePath, dirPath);

	// write data file, include tag information
	ofstream ofs(samplePath);
	if (! ofs.good())
		error("Couldn't open data file for writing: ",samplePath);

	cout << "Writing: " << samples_.size() << " samples.." << endl;

	for (size_t i=0; i<samples_.size(); i++)
			samples_[i].writeToStream(ofs, true);
	ofs.close();
	
	// write phi file
	string phiPath = dirPath + "/" + fileName + ".phi";
	ofs.open(phiPath.c_str());
	if (! ofs.good())
	{
		error("Could not write feedback file: ",phiPath.c_str());
	}
	else
	{
		for (size_t i=0; i<phi_.size(); i++)
			ofs << phi_[i].idx0 << "\t" << phi_[i].idx1 << "\t" << scientific << setprecision(4) << phi_[i].weight << endl;

		ofs.close();
	}
}


