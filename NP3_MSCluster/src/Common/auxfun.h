#ifndef __AUXFUN_H__
#define __AUXFUN_H__

#include "includes.h"

bool copyFile(const char* source, const char* dest);
bool checkDirPathExists( const char* userPath );
void createDirIfDoesNotExist( const char* dirPath, int verboseLevel=0);
void removeFilesInList(const char* listPath);
void getDirFromFullPath(const char* fullPath, string& dirPath);
size_t getFileSize(const char* sFileName);

// return first token until "\t"
inline
string getFirstToken(const char* line, size_t maxLength)
{
	for (size_t i=0; i<maxLength; i++)
		if (line[i] == '\t')
			return std::string(line,i);
	return (std::string());
}


/*! \fn readListOfPaths
	Reads a list of paths form a text file.

	If the first line is a number, it is assumed that this is the first file's index in a larger list.
	All file indexes of this list should be relative to that number.

	@param listPath The path to the text file with full paths to the files.
	@param paths The vector of file paths that were read.
	@param  verbose Give status report (default false).

	@return The first files index (0 if no number is given in the first line of the file).
*/
size_t readListOfPaths(const char* listPath, vector<string>& paths, bool verbose=false);

void   writeListOfPaths(const char* listPath, const vector<string>& paths);

/*! \fn unzipFile
	Unzips a single zip file.

	The function unzips the zip file into the same directory as the zip file.

	@param	zipPath The path to the zip file.
	@param  unzippedNames The names of the extracted files.

	@return number of extracted files.*/
size_t unzipSingleFile(const string& zipPath, vector<string>& unzippedNames);



// returns number of bytes read
size_t getFirstLineInFile(const char* filePath, char* buffer, size_t maxCharsToRead);



/*! \fn getFileNameWithoutExtension
	\brief removes the prefix path and the extension after the last "." in the file name
*/
void getFileNameWithoutExtension(const char* fullPath, string& fileName);

/*! \fn checkIfFileExists
	\brief returns true if file can be opened for reading.
*/
bool checkIfFileExists(const char* fullPath);


/*! \fn stripPathOfTrailingSymbols
	\brief removes any trailiing '/' or '\' from the path
*/
string stripPathOfTrailingSymbols(const char* path);

string getTimeString();
unsigned int getDateAsInt();


/*! \fn computeRocAuc
	\brief Computes an ROC curve's area under the curve.

	Cumputes the area under the ROC curve. Assumes that both recall and
	precision vectors have the same dimensionality and cover the range [0-1]. Also assumes
	that both vectors are monotonically increasing (or decreasing).

	@param precision vector of precision values (in range [0-1]).
	@param recall vector of recal values (in range [0-1]).

	@return The area under the curve  (values in range [0-1]).
*/
double computeRocAuc(vector<double>& precision, vector<double>& recall);

/****************************************************************
Efficiently computes the mean and sd of a variable
*****************************************************************/
struct MeanSdStats {
	MeanSdStats() : sumW(0.0), sumWX(0.0), sumWX2(0.0) {};
	void calcMeanAndSd(double& m, double& sd) const;
	void addWX(double w, double x) { sumW+=w; double wx=w*x; sumWX+=wx; sumWX2+=wx*x; }

	double sumW;
	double sumWX;
	double sumWX2;
};



void seedRandom (unsigned int init = 0);

unsigned int getRandomSeed();

/* Returns random uniform number */
double myRandom();

void error();
void error(const char *msg);
void error(const char* msg1, const char* msg2);
void error(const char* msg1, size_t num);
void error(const char* msg1, size_t num1, const char* msg2, size_t num2);
void error(const char* msg1, const char* mssg2, const char* msg3);

void chooseKFromN(size_t k, size_t n, vector<size_t>& idxs);
void randomPermutation(size_t k, vector<size_t>& perm);

string makeRangeString(const vector<size_t>& idxs);



void generate_all_permutations(const vector<int>& org_vector, 
							   vector< vector<int> >& permutations);

// returns the minimal x for which the cumulative probability
// P(X<x)>= target_prob, assuming X~bin(n,p)
int get_min_number_from_binomial_prob(int n, double p, double target_prob);


void computeBinomailCDFs(int n, float p, vector<float>& cdf);



void split_string(const string& str, vector<string>& results, char delim ='\t');

template<class T>
void create_histogram(vector<T>& vals, size_t num_bins, T min_val,
					  T max_val , ostream& os =cout)
{
	T bin_size = (max_val-min_val)/ (T)(num_bins);
	T one_over_bin = 1.0 / bin_size;
	vector<int> counts;
	vector<float> percents;
	counts.resize(num_bins,0);

	for (size_t i=0; i<vals.size(); i++)
	{
		if (vals[i]<min_val)
		{
			counts[0]++;
			continue;
		}

		int bin_idx = num_bins-1;
		if (vals[i]<max_val)
			bin_idx = (int)(one_over_bin*(vals[i]-min_val));
		
		counts[bin_idx]++;
	}

	T v = min_val;
	int tc=0;
	for (size_t i=0; i<num_bins; i++)
	{
		os << setw(4) <<  setprecision(2) << left << v << " - ";
		v+= bin_size;
		os <<  setw(4) << left << v  << "\t" <<  setw(6) << right << counts[i] << "  ";
		os << setw(6) << left << setprecision(4) << (float)counts[i]/(float)vals.size() << endl;
		tc+= counts[i];
	}

	os << "Total:       " << setw(6) << right << tc << "  " << setw(6) << setprecision(4) << left << (float)tc/(float)vals.size() << endl;
}


template<class T>
void create_histogram(const vector<T>& vals, const vector<T>& separator_vals, 
					  vector<int>& counts, ostream& os =cout)
{

	int i;
	
	vector<float> percents;
	counts.resize(separator_vals.size()+1,0);

	for (i=0; i<vals.size(); i++)
	{
		int j;
		for (j=0; j<separator_vals.size(); j++)
			if (vals[i]<separator_vals[j])
				break;
		
		counts[j]++;
	}

/*	cout << "VALS: " << vals.size() << endl;

	
	int tc=0;

	for (i=0; i<counts.size(); i++)
	{
		T v = (i==0) ? 0 : separator_vals[i-1];
		os << setw(4) <<  setprecision(2) << left << v << " - ";
		v = ( i == counts.size() -1 ) ? 0 : separator_vals[i];
		os <<  setw(4) << left << v  << "  " <<  setw(6) << right << counts[i] << "  ";
		os << setw(6) << left << setprecision(4) << (float)counts[i]/(float)vals.size() << endl;
		tc+= counts[i];
	}

	os << "Total:       " << setw(6) << right << tc << "  " << setw(6) << setprecision(4) << left << (float)tc/(float)vals.size() << endl;
	*/
}



template<class T>
void calc_mean_sd(const vector<T>& v, T *mean, T *sd)
{
	T m=0,var=0;
	int i;

	if (v.size() == 0)
	{
		*mean=0;
		*sd=0;
		return;
	}

	if (v.size() == 1)
	{
		*mean=v[0];
		*sd=0;
	}

	for (i=0; i<v.size(); i++)
		m+=v[i];

	m/=v.size();

	for (i=0; i<v.size(); i++)
		var+=(v[i]-m)*(v[i]-m);

	var /= v.size();

	*mean=m;
	*sd = sqrt(var);
}


template<class T>
void calc_mean_sd_from_counts(const vector<T>& vals, const vector<int>& counts,
							  double *mean, double *sd)
{
	double m=0,var=0;
	int i;

	if (vals.size() == 0)
	{
		*mean=0;
		*sd=0;
		return;
	}

	if (vals.size() != counts.size())
	{
		cout << "Error: values and counts should have same dimension!" << endl;
		exit(1);
	}

	int total_counts=0;
	for (i=0; i<vals.size(); i++)
	{
		m+=vals[i]*counts[i];
		total_counts+=counts[i];
	}

	m/=(double)total_counts;

	for (i=0; i<vals.size(); i++)
		var+=(vals[i]-m)*(vals[i]-m)*counts[i];

	var /= (double)total_counts;

	*mean=m;
	*sd = sqrt(var);
}




template<class T>
void mergeSortedVectors(const vector<T>& a, const vector<T>&b, vector<T>& m)
{
	typename vector<T>::const_iterator it_a = a.begin();
	typename vector<T>::const_iterator it_b = b.begin();

	m.clear();
	m.reserve(a.size()+b.size());

	while (it_a != a.end() && it_b != b.end())
	{
		if (*it_a < *it_b)
		{
			m.push_back(*it_a++);
			continue;
		}

		if(*it_a > *it_b)
		{
			m.push_back(*it_b++);
			continue;
		}

		m.push_back(*it_a++);
		it_b++;
	}

	while (it_a != a.end())
		m.push_back(*it_a++);

	while (it_b != b.end())
		m.push_back(*it_b++);
}

template<class T>
void intersectSortedVectors(const vector<T>& a, const vector<T>&b, vector<T>& intersec)
{
	intersec.clear();
	intersec.reserve(a.size()>b.size() ? b.size() : a.size());
	
	typename vector<T>::const_iterator it_a = a.begin();
	typename vector<T>::const_iterator it_b = b.begin();

	while (it_a != a.end() && it_b != b.end())
	{
		if (*it_a < *it_b)
		{	
			it_a++;
			continue;
		}

		if(*it_a > *it_b)
		{
			it_b++;
			continue;
		}

		intersec.push_back(*it_a++);
		it_b++;
	}
}


template<class T>
void set_minus(const vector<T>& a, const vector<T>& b, vector<T>& diff)
{
	int i;
	diff.clear();
	for (i=0; i<a.size(); i++)
	{
		int j;
		for (j=0; j<b.size(); j++)
			if (a[i]==b[j])
				break;
		if (j==b.size())
			diff.push_back(a[i]);
	}
}

template<class T>
void set_overlap(const vector<T>& a, const vector<T>& b, vector<T>& overlap)
{
	int i;
	overlap.clear();
	for (i=0; i<a.size(); i++)
	{
		int j;
		for (j=0; j<b.size(); j++)
			if (a[i]==b[j])
				break;
		if (j<b.size())
			overlap.push_back(a[i]);
	}
}

struct PermutePair {
	bool operator< (const PermutePair& other) const
	{
		return randVal<other.randVal;
	}
		int    orgIdx;
		double randVal;
};

template<class T>
void permute_vector(vector<T>& vec)
{
	

	vector<PermutePair> idxPairs;
	idxPairs.resize(vec.size());
	int i;
	for (i=0; i<vec.size(); i++)
	{
		idxPairs[i].orgIdx=i;
		idxPairs[i].randVal = myRandom();
	}

	sort(idxPairs.begin(),idxPairs.end());
	vector<T> new_vec;
	new_vec.resize(vec.size());
	for (i=0; i<vec.size(); i++)
		new_vec[i]=vec[idxPairs[i].orgIdx];
	
	vec = new_vec;
}


template<typename T>
T addInLogSpace(const T& logX, const T& logY)
{
	const T diff = logX-logY;
	if (diff>32.0)
		return logX;
	if (diff<-32.0)
		return logY;

	if (logX>logY)
		return (logX + log(1.0 + exp(-diff)));

	return (logY + log(1.0 + exp(diff))); 
}

template<typename T>
void printVector(const vector<T>& v)
{
	if (v.size() == 0)
	{
		cout << "empty" << endl;
		return;
	}

	cout << setprecision(4);
	for (size_t i=0; i<v.size(); i+=10)
	{
		if (i % 10 == 0)
			cout << i <<":";
		for (size_t j=i; j<i+10 && j<v.size(); j++)
			cout << "\t" << v[j];
		cout << endl;
	}
}

template<typename T>
void printVector(const vector< vector<T> >& v)
{
	if (v.size() == 0)
	{
		cout << "empty" << endl;
		return;
	}

	size_t maxIdx = 0;
	for (size_t i=1; i<v.size(); i++)
		if (v[i].size()>v[maxIdx].size())
			maxIdx=i;

	cout << setprecision(4);
	for (size_t i=0; i<v[maxIdx].size(); i++)
	{
		cout << i <<":";
		for (size_t c=0; c<v.size(); c++)
		{
			if (i<v[c].size())
			{
				cout << "\t" << v[c][i];
			}
			else
				cout << "\t ";		
		}
		cout << endl;
	}
}



// moves elements to start of vector
template<class T>
void shuntInVector(vector<T>& v, size_t orgIdx,  size_t numElements, bool resizeInd = false)
{
	if (orgIdx + numElements > v.size())
		numElements = v.size() - orgIdx;
	
	copy(v.begin()+orgIdx, v.begin()+(orgIdx+numElements), v.begin());

	if (resizeInd)
		v.resize(numElements);
}


#endif 
