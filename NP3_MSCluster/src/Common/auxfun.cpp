/************************************************************************//**
@file auxfun.cpp
\brief Holds general auxilary functions that can be useful for the various
projects.
*/
#include "auxfun.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#ifndef WIN32
#include <unistd.h>
#endif

void MeanSdStats::calcMeanAndSd(double& m, double& sd) const
{
	if (sumW<=0)
	{
		m=NON_FLOAT;
		sd=NON_FLOAT;
		return;
	}

	m=sumWX/sumW;
	sd= sqrt((sumWX2 - 2*m*sumWX)/sumW +m*m);
}

double computeRocAuc(vector<double>& precision, vector<double>& recall)
{
	assert(precision.size() == recall.size());
	bool reversed = false;

	// reverse if recall doesn't go in the direction from 0 to 1.0
	if (recall.front()>recall.back())
	{
		reversed = true;
		reverse(precision.begin(), precision.end());
		reverse(recall.begin(), recall.end());
	}

	assert(recall.front() >=0.0 && recall.back()<=1.0);
	assert(precision.front() <= 1.0 && precision.back()>=0.0);


	const size_t size = recall.size();
	double auc = 0.0;

	if (recall[0]>0.0)
	{
		const double avgPrecision = 0.5*(1.0 + precision[0]);
		auc += recall[0] * avgPrecision;
	}

	for (size_t i=1; i<size; i++)
	{
	//	assert(precision[i]<=precision[i-1] && recall[i]>=recall[i-1]);
		const double avgPrecision = 0.5 *(precision[i] + precision[i-1]);
		const double deltaRecall  = recall[i]-recall[i-1];
		auc += deltaRecall * avgPrecision;
	}

	if (recall.back()<1.0)
	{
		const double avgPrecision = 0.5*precision.back();
		const double deltaRecall  = 1.0 -  recall.back();
		auc += fabs(deltaRecall * avgPrecision); 
	}

	// reverse back if needed
	if (reversed)
	{
		reverse(precision.begin(), precision.end());
		reverse(recall.begin(), recall.end());
	}
	return auc;
}


static unsigned int SEED;	

void seedRandom (unsigned int init)   {
	if (init != 0)
	{
		SEED = init;
	}
	else
	{
		time_t ltime;
		unsigned int t=(unsigned int)time( &ltime );

		SEED = t;
	}
}

unsigned int getRandomSeed() { return SEED; }


/* Returns random uniform number */
double myRandom()  
{
  static unsigned int a = 1588635695, m = 4294967291U, q = 2, r = 1117695901;

   SEED = a*(SEED % q) - r*(SEED / q);
   return ((double)SEED / (double)m);
}

void error()
{
	cout << "Error!!!" << endl;
	#ifdef WIN32
		system("pause");
	#endif
	exit(1);

}

void error(const char *msg)
{
	cout << "Error: " << msg << endl;
	#ifdef WIN32
		system("pause");
	#endif
	exit(1);
}

void error(const char *msg1, const char* msg2)
{
	cout << "Error: " << msg1 << msg2 << "." << endl;
	#ifdef WIN32
		system("pause");
	#endif
	exit(1);
}

void error(const char *msg1, size_t num)
{
	cout << "Error: " << msg1 << num << "." << endl;
	#ifdef WIN32
		system("pause");
	#endif
	exit(1);
}

void error(const char *msg1, size_t num1, const char *msg2, size_t num2)
{
	cout << "Error: " << msg1 << num1 << msg2 << num2 << "." << endl;
	#ifdef WIN32
		system("pause");
	#endif
	exit(1);
}

void error(const char *msg1, const char* msg2, const char *msg3)
{
	cout << "Error: " << msg1 << msg2 << msg3 << "." << endl;
	#ifdef WIN32
		system("pause");
	#endif
	exit(1);
}





struct ChoosePair {
	ChoosePair(size_t i,double d) : idx(i),val(d) {};
	bool operator < (const ChoosePair& other) const
	{
		return val<other.val;
	}
	size_t idx;
	double val;
};



// chooses k numbers from 0,...,n-1 (unique)
void chooseKFromN(size_t k, size_t n, vector<size_t>& idxs)
{
	if (k>n)
	{
		cout << "Error: choose " << k << " from " << n << " !" << endl;
		exit(1);
	}

	idxs.clear();
	idxs.resize(k);
	vector<ChoosePair> pairs;

	size_t i;
	for (i=0; i<n; i++)
		pairs.push_back(ChoosePair(i,myRandom()));
	
	sort(pairs.begin(),pairs.end());

	for (i=0; i<k; i++)
		idxs[i]=pairs[i].idx;

	sort(idxs.begin(),idxs.end());
}

void randomPermutation(size_t k, vector<size_t>& perm)
{
	perm.clear();
	perm.resize(k);
	if (k == 0)
		return;
	
	vector<ChoosePair> pairs;
	for (size_t i=0; i<k; i++)
		pairs.push_back(ChoosePair(i,myRandom()));
	
	sort(pairs.begin(),pairs.end());

	for (size_t i=0; i<k; i++)
		perm[i]=pairs[i].idx;
}




bool checkDirPathExists( const char* userPath)
{
#if defined(WIN32) || defined(WIN64)
	string tmpFile = std::string(userPath) + "\\tmpTMPtmp.XXX";
#else
	string tmpFile = std::string(userPath) + "/tmpTMPtmp.XXX";
#endif
	FILE* stream=fopen(tmpFile.c_str(),"w");
	if (! stream)
		return false;

	fclose(stream);
	unlink(tmpFile.c_str());
	return true;
//   struct stat statBuf;
//   return ( stat( userPath, &statBuf ) == 0 && (statBuf.st_mode & S_IFDIR));
   
}


void createDirIfDoesNotExist( const char* dirPath, int verboseLevel)
{
	if (! checkDirPathExists(dirPath))
	{
		char cmd[256];
		sprintf(cmd,"mkdir \"%s\"",dirPath);
		const int retVal = system(cmd);
		if (! checkDirPathExists(dirPath))
		{
			cout << "Error executing command: " << cmd << endl;
			error("Could not create directory: ",dirPath);
		}
		if (verboseLevel)
			cout << "Made directory: " << dirPath << endl;
        return;
	}
}

void removeFilesInList(const char* listPath)
{
	vector<string> paths;
	readListOfPaths(listPath, paths);
	for (size_t i=0; i<paths.size(); i++)
		remove(paths[i].c_str());
}

bool copyFile(const char* source, const char* dest)
{
	ostringstream oss;
#if defined(WIN32) || defined(WIN64)
	oss << "copy \"" << source << "\" \"" << dest << "\"";
#else
	oss << "cp \"" << source << "\" \"" << dest << "\"";
#endif
//	cout << "cmd: " << oss.str() << endl;
	int retVal = system(oss.str().c_str());
	return (retVal == 0);
}

// returns the file size in bytes
size_t getFileSize(const char* sFileName)
{
  ifstream f;
  f.open(sFileName, ios_base::binary | ios_base::in);
  if (!f.good() || f.eof() || !f.is_open()) { return 0; }
  f.seekg(0, ios_base::beg);
  ifstream::pos_type begin_pos = f.tellg();
  f.seekg(0, ios_base::end);
  return static_cast<size_t>(f.tellg() - begin_pos);
}

size_t getFirstLineInFile(const char* filePath, char* buffer, size_t maxCharsToRead)
{
	ifstream f(filePath);
	if (!f.good() || f.eof() || !f.is_open())
		return 0;

	f.getline(buffer,maxCharsToRead);
	size_t numRead = f.gcount();
	f.close();
	return numRead;
}



bool checkIfFileExists(const char* fullPath)
{
	ifstream f(fullPath);
	if (f.good())
	{
		f.close();
		return true;
	}

	f.close();
	return false;
}


void getFileNameWithoutExtension(const char* fullPath, string& fileName)
{
	const string fullPathString = std::string(fullPath);
	const int last_pos = strlen(fullPath)-1;
	int dot_pos = last_pos;
	while (dot_pos>=0 && fullPath[dot_pos] != '.')
		dot_pos--;

	if (dot_pos <=0)
		dot_pos = last_pos;

	int slash_pos=last_pos;
	while (slash_pos>=0 && (fullPath[slash_pos] != '/') && (fullPath[slash_pos] != '\\'))
		slash_pos--;

	fileName = fullPathString.substr(slash_pos+1,dot_pos-slash_pos-1);
}

void getDirFromFullPath(const char* fullPath, string& dirPath)
{
	const string fullPathString = std::string(fullPath);
	const int last_pos = fullPathString.length()-1;

	int slash_pos=last_pos;
	while (slash_pos>=0 && (fullPath[slash_pos] != '/') && (fullPath[slash_pos] != '\\'))
		slash_pos--;

	if (slash_pos <= 0)
	{
		dirPath =  std::string(".");
	}
	else
		dirPath = fullPathString.substr(0, slash_pos);
}

/*! \fn stripPathOfTrailingSymbols
	\brief removes any trailiing '/' or '\' from the path
*/
string stripPathOfTrailingSymbols(const char* path)
{
	const string fullPathString = std::string(path);
	int last_pos = fullPathString.length()-1;
	while (last_pos>=0 && (fullPathString[last_pos] == '/' || fullPathString[last_pos] == '\\'))
		last_pos--;

	if (last_pos <= 0)
	{
		return  std::string(".");
	}
	else
		return fullPathString.substr(0, last_pos+1);
}


size_t unzipSingleFile(const string& zipPath, vector<string>& unzippedNames)
{
	static bool seeded = false;
	if (! seeded)
		seedRandom();

	unzippedNames.clear();
	string fileName;
	getFileNameWithoutExtension(zipPath.c_str(), fileName);
	fileName += ".ZIPRES_";
	ostringstream oss;
	oss << fileName << time(NULL) % 10000;
	string tmpName = oss.str();
	const string unzipCommand = "unzip -o -d . " + zipPath + " > " + tmpName;
	int systemRetVal = system(unzipCommand.c_str());
	if (systemRetVal != 0)
	{
		unlink(tmpName.c_str());
		systemRetVal = system("sleep 100");
		systemRetVal = system(unzipCommand.c_str());
		if (systemRetVal)
		{
			cout << "Warning: Could not perform unzip command correctly: " << unzipCommand.c_str() << endl;
			return 0;
		}
	}

	systemRetVal = system("sleep 4"); // give the zip file a chance to close (when system is busy?)

	FILE* stream = fopen(tmpName.c_str(),"r");
	if (! stream)
		error("Could not read results of unzip command: ",unzipCommand.c_str());

	char line[256];
	while (fgets(line,256,stream))
	{
		istringstream iss(line);
		string firstArg, fileName;
		iss >> firstArg >> fileName;
		if (firstArg == "inflating:")
		{
			unzippedNames.push_back(fileName);
			cout << "Unzipped " << fileName << endl;
		}
	}

	
	unlink(tmpName.c_str());
	return unzippedNames.size();
}









/*************************************************************
   finds all the permutaitons of n elements, repeated elements
   are allowed and do not create redundant permutations.
**************************************************************/
void generate_all_permutations(const vector<int>& org_vector, 
							   vector< vector<int> >& permutations)
{
	
	vector<int> counts, symbols;
	permutations.clear();

	if (org_vector.size() == 0)
		return;

	counts.clear();
	symbols.clear();

	// create vector with symbols and their counts
	symbols.push_back(org_vector[0]);
	counts.push_back(1);

	for (size_t i=1; i<org_vector.size(); i++)
	{
		
		size_t j;
		for (j=0; j<counts.size(); j++)
		{
			if (org_vector[i] == symbols[j])
			{
				counts[j]++;
				break;
			}
		}

		if (j == counts.size())
		{
			symbols.push_back(org_vector[i]);
			counts.push_back(1);
		}
	}

	vector<int> next_sym_idx,perm;
	int n = org_vector.size(); // total number of elements
	int k = counts.size(); // total number of element types
	next_sym_idx.resize(n,0);
	perm.resize(n,-1);
	int d=0;

	while (true)
	{
		while (next_sym_idx[d]<k && counts[next_sym_idx[d]] == 0)
			next_sym_idx[d]++;

		if (next_sym_idx[0]==k)
			break;

		if (next_sym_idx[d] >= k)
		{
			next_sym_idx[d]=0;
			d--;
			counts[next_sym_idx[d]]++;
			next_sym_idx[d]++;
			continue;
		}

		// add symbol
		perm[d]=symbols[next_sym_idx[d]];
		counts[next_sym_idx[d]]--;
		d++;

		if (d == n)
		{
			permutations.push_back(perm);
	//		int k;
	//		for (k=0; k<perm.size(); k++)
	//			cout << perm[k] << " ";
	//		cout << endl;

			d--;
			counts[next_sym_idx[d]]++;
			next_sym_idx[d]++;
		}
	}
}

// returns the minimal x for which the cumulative probability
// P(X<x)>= target_prob, assuming X~bin(n,p)
int get_min_number_from_binomial_prob(int n, double p, double target_prob)
{
	const double one_minus_p = 1.0 - p;
	double pow_1_minus_p = pow(one_minus_p,n);

	double sum_prob = pow_1_minus_p;
	double bin_coef = 1.0;
	double pow_val  = pow_1_minus_p;

//	cout << 0 << " " <<  pow_val << " " << pow_val << endl;
	int b=0;
	while (sum_prob<target_prob)
	{
		b++;
		bin_coef *= (double)(n-b+1);
		bin_coef /= (double)b;

		pow_val *= p;
		pow_val /= one_minus_p;

		double prob = bin_coef * pow_val;
		sum_prob += prob;

	//	cout << b << " " << prob << " " << sum_prob << endl;
	}
	return b;
}

void computeBinomailCDFs(int n, float p, vector<float>& cdf)
{
	cdf.clear();
	if (n<=0 || p<=0.0 || p>=1.0)
		return;

	cdf.resize(n+1,0.0);
	const double pOverOneMinusP = p/(1.0-p);
	double binomialCoeff = 1.0;				// n choose k
	double powerValue    = pow(1.0-p, n);   // p^k * (1-p)^n-k
	cdf[0] = powerValue;

	for (int k=1; k<=n; k++)
	{
		binomialCoeff *= static_cast<double>(n-k+1.0);
		binomialCoeff /= static_cast<double>(k);
		powerValue	  *= pOverOneMinusP;
		cdf[k]= cdf[k-1] + static_cast<float>(binomialCoeff * powerValue);
	}
}





/******************************************************************************
splits string according to a given delimeter char.
*******************************************************************************/
void split_string(const string& str, vector<string>& results, char delim)
{
	const char *str_buff = str.c_str();
	const int   max_pos  = str.length();

	results.clear();
	int last_pos = 0;
	int pos=0;

	for (pos=0; pos<max_pos; pos++)
	{
		if (str_buff[pos]==delim || pos == max_pos-1)
		{
			if (pos>last_pos)
				results.push_back(str.substr(last_pos,pos-last_pos));
			last_pos=pos+1;
			continue;
		}
	}
}


// gets indexes and outputs a string of the form: 1-4,7,9-12
string makeRangeString(const vector<size_t>& orgIdxs)
{
	vector<size_t> idxs = orgIdxs;
	sort(idxs.begin(),idxs.end());

	ostringstream oss;
	if (idxs.size()==0)
		return oss.str();

	vector<bool> skips(idxs.size(),false);
	for (size_t i=1; i<idxs.size()-1; i++)
		if (idxs[i-1]+1 == idxs[i] && idxs[i]+1   == idxs[i+1])
			skips[i]=true;

	bool inSkip=false;
	for (size_t i=0; i<idxs.size(); i++)
	{
		if (! skips[i])
		{
			oss << idxs[i];
			inSkip=false;
			if (i<idxs.size()-1 && ! skips[i+1])
				oss << ",";
		}
		else
		{
			if (inSkip==false)
			{
				inSkip=true;
				oss << "-";
			}
		}
	}
	return oss.str();
}


size_t readListOfPaths(const char* listPath, vector<string>& paths, bool verbose)
{
	FILE* stream=fopen(listPath,"r");
	if (! stream)
	{
		cout << "Error: couldn't open list file: " << listPath << endl;
		exit(1);
	}

	paths.clear();

	
	char buffer[512], pathBuffer[512];

	if (! fgets(buffer,512,stream))
		return 0;

	// check if first line is a number
	int startIdx = -1;
	bool gotStartIdx = false;
	if (buffer[0] != '#' && buffer[0]>='0' && buffer[0]<='9')
	{
		gotStartIdx = (sscanf(buffer,"%d",&startIdx) == 1 && startIdx>=0);
	}
	
	if (! gotStartIdx)
	{
		startIdx = 0;
		sscanf(buffer,"%s",pathBuffer);
		if (strlen(pathBuffer)>=2)
			paths.push_back(static_cast<string>(pathBuffer));
	}
	

	// read the rest
	while ( fgets(buffer,512,stream) )
	{
		if (buffer[0]=='#')
			continue;

		// remove trailing
		sscanf(buffer,"%s",pathBuffer);
		if (strlen(pathBuffer)>=2)
			paths.push_back(static_cast<string>(pathBuffer));
	}
	fclose(stream);

	if (verbose)
		cout << "Read " << paths.size() << " paths, start idx is " << startIdx << endl;
	return (static_cast<size_t>(startIdx));
}

void   writeListOfPaths(const char* listPath, const vector<string>& paths)
{
	FILE* stream=fopen(listPath,"w");
	if (! stream)
	{
		cout << "Error: couldn't open list file: " << listPath << endl;
		exit(1);
	}

	for (size_t i=0; i<paths.size(); i++)
		fprintf(stream,"%s\n",paths[i].c_str());
	fclose(stream);
}


string getTimeString()
{
	static time_t rawtime;
	time( &rawtime );
	string timeStr = std::string(ctime(&rawtime));
	return timeStr.substr(0,timeStr.length()-1);
}

unsigned int getDateAsInt()
{
	static time_t rawtime;
	tm * ptm;

	time( &rawtime );
	ptm = localtime(&rawtime);
	return ((ptm->tm_year + 1900) * 10000 + (ptm->tm_mon +1 ) * 100 + ptm->tm_mday);
}




