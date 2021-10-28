#include "ScanList.h"
#include "../Common/auxfun.h"

size_t ScanListManager::initialize(const char *list, 
										mass_t minMz, 
										mass_t maxMz)
{
	vector<string> paths;
	readListOfPaths(list, paths);

	if (paths.size() == 0)
	{
		cout << "Warning: No exclusion files found!" << endl;
		return 0;
	}

	size_t numExclusions = 0;
	for (size_t i=0; i<paths.size(); i++)
	{
		ifstream excludeStream(paths[i].c_str());
		if (! excludeStream.good())
		{
			cout << "Warning: could not open " << paths[i].c_str() << endl;
			continue;
		}
		

		float minExcludeMz=MIN_FLOAT, maxExcludeMz=MAX_FLOAT;
		char buffer[128];
		excludeStream.getline(buffer,128);
		if ( sscanf(buffer,"#RANGE\t%f\t%f",&minExcludeMz, &maxExcludeMz) != 2)
		{
			minExcludeMz=MIN_FLOAT; 
			maxExcludeMz=MAX_FLOAT;
			ScanEntry se;
			if (sscanf(buffer,"%d\t%d\t%d\t%f", &se.datasetIdx, &se.fileIdx, &se.scan, &se.mz)>=3)
			{
				scans_[se] = i;
				numExclusions++;
			}
		}
		else
		{
			if (minMz > maxExcludeMz || maxMz < minExcludeMz)
			{
				excludeStream.close();
				continue;
			}
		}

		while (! excludeStream.eof())
		{
			excludeStream.getline(buffer,128);
			ScanEntry se;
			if (sscanf(buffer,"%d\t%d\t%d\t%f", &se.datasetIdx, &se.fileIdx, &se.scan, &se.mz) < 3)
				continue;

			scans_[se] = i;
			numExclusions++;
		}
		excludeStream.close();
	}

	return numExclusions;
}


struct ScanStruct {
	ScanStruct() : datasetIdx(-1), fileIdx(-1), scan(-1), mz(0.0), charge(0) {}
	bool operator< (const ScanStruct& rhs) const
	{
		return ( datasetIdx <  rhs.datasetIdx ||
				(datasetIdx == rhs.datasetIdx && fileIdx <  rhs.fileIdx) ||
				(datasetIdx == rhs.datasetIdx && fileIdx == rhs.fileIdx && scan < rhs.scan));
	}

	int datasetIdx;
	int fileIdx;
	int scan;
	float mz;
	int charge;
};

void makeScanListFromClust(const char* list, const char* scanFilePath)
{
	vector<string> paths;
	readListOfPaths(list, paths);

	mass_t minMz=MAX_FLOAT;
	mass_t maxMz=0.0;
	vector<ScanStruct> allScans;

	for (size_t i=0; i<paths.size(); i++)
	{
		ifstream clustStream(paths[i].c_str());
		if (! clustStream.good())
			error("Couldn't open file for reading: ", paths[i].c_str());

		int numScans=0;
		char buffer[256];
		while ( ! clustStream.eof())
		{
			clustStream.getline(buffer,256);
			if (clustStream.gcount()<5)
				continue;
			istringstream iss(buffer);
			string title;
			int size=0;
			iss >> title >> size;
			assert(size>0);
			for (size_t j=0; j<size; j++)
			{
				clustStream.getline(buffer,256);
				istringstream iss(buffer);
				ScanStruct s;
				iss >> s.datasetIdx >> s.fileIdx >> s.scan >> s.mz >> s.charge;
				if (s.datasetIdx<0 || s.fileIdx<0 || s.scan<0 || s.mz<=0.0)
				{
					cout << "Bad line: " << buffer << endl;
					cout << "File: " << paths[i] << endl;
					exit(1);
				}
				assert(s.datasetIdx>=0 && s.fileIdx>=0 && s.scan>=0 && s.mz>0.0);
				allScans.push_back(s);
				if (s.mz<minMz)
					minMz = s.mz;
				if (s.mz>maxMz)
					maxMz = s.mz;
				numScans++;
			}
		}
		clustStream.close();
		cout << i << "\t" << numScans << "\t" << paths[i] << endl;
	}

	sort(allScans.begin(),allScans.end());

	// write scan list
	ofstream ofs(scanFilePath);
	if (! ofs.good())
		error("couldn't open scanlist for writing: ",scanFilePath);

	int numSame=0;
	// NP3 GOT precision from 3 to 4
	ofs << fixed << setprecision(4);
	ofs << "#RANGE\t" << minMz << "\t" << maxMz << endl;
	for (size_t i=0; i<allScans.size(); i++)
	{
		if (i>0 && allScans[i].scan == allScans[i-1].scan && allScans[i].fileIdx == allScans[i-1].fileIdx
				&& allScans[i].datasetIdx == allScans[i-1].datasetIdx)
		{
			numSame++;
			if (numSame<10)
				cout << "SAME: " << allScans[i].datasetIdx << "\t" << allScans[i].fileIdx
					 << "\t" << allScans[i].scan << "\t" << allScans[i].mz << endl;
		}
		ofs << allScans[i].datasetIdx << "\t" << allScans[i].fileIdx << "\t"
			<< allScans[i].scan << "\t" << allScans[i].mz << "\t" << allScans[i].charge << endl;
	}

	ofs.close();
	cout << "Wrote " << allScans.size() << " scans to " << scanFilePath << endl;
	cout << "Found " << numSame << " scans that has same g/f/s" << endl;
}






