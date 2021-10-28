#ifndef __DBFILE_H__
#define __DBFILE_H__


/*! \file DbFile.h
	\brief holds the class that reads a fasta file and stores it (with tags)
*/

#include "../Common/includes.h"
#include <set>

struct TagWithLocations {
	TagWithLocations() : tag(0), tagLength(0), numLocations(0), locationsOffset(0) {}

	TagWithLocations(const char* t, unsigned int tl) : tag(t), tagLength(tl), numLocations(0), locationsOffset(0) {}

	bool operator< (const TagWithLocations& rhs) const
	{
		if (tagLength < rhs.tagLength)
			return true;
		if (tagLength > rhs.tagLength)
			return false;
		return (strncmp(tag, rhs.tag, tagLength)<0);
	}

	const char* tag;				   /// the tag itself (pointer to location in some char array)
	unsigned int tagLength;	   /// number of amino acids 
	mutable unsigned int numLocations; /// how many places this tag appears in the fasta file
	mutable unsigned int locationsOffset;	   /// pointer to offset in the list of locations
};






class DbFile {
public:

	DbFile() : 	fasta_(std::string()), path_(std::string()), partIdx_(0), minTagLength_(0), maxTagLength_(0) {}

	bool readFastaFile(const string& fasta, size_t partIdx = 0);

	bool writeFastaFile(const string& path) const;

	bool createTags(unsigned int minTagLength, unsigned int maxTagLength);

	bool writeDbFile(const string& path) const;

	bool readDbFile(const string& path);

	size_t getNumProteins() const { return proteinStartPointers_.size(); }

	void setPartIdx(int partIdx) { partIdx_= partIdx; }
	
	const char* getSequence(size_t idx) const {
		if (idx>=proteinStartPointers_.size())
			return 0;
		return (&aaArray_[0] + proteinStartPointers_[idx]);
	}


	string getName(size_t idx)
	{
		if (idx>=proteinNames_.size())
			return (std::string());
		return (std::string(&nameArray_[proteinNames_[idx]]));
	}

	size_t getNumSequences() const { return proteinStartPointers_.size(); }

	void computeProteinIdxAndOffsetFromLocation(unsigned int location, size_t& proteinIdx, size_t& offset) const;

	bool getTagWithLocations(TagWithLocations& twl) const
	{
		set<TagWithLocations>::const_iterator it = allTags_.find(twl);
		if (it == allTags_.end())
		{
			twl.tag = 0;
			twl.numLocations = 0;
			return false;
		}
		twl = *it;
		return true;
	}

	const char* getAaPointer(unsigned int location) const { return &aaArray_[location]; }
	const char* getFirstAAPtr() const { return &aaArray_[0]; }
	const char* getLastAAPtr()  const { return (&aaArray_[0] + aaArray_.size()); }

	unsigned int getLocation(const char* ptr) const {
		assert( ptr>=&aaArray_[0] );
		assert (ptr <= (&aaArray_[0] + aaArray_.size()) );
		return (ptr-&aaArray_[0]); 
	}

	size_t getLocationPointers(const TagWithLocations& twl, vector<const char*>& ptrs) const;

	unsigned int getProteinStartPointer(size_t i) const { return proteinStartPointers_[i]; }
	size_t getAaArraySize() const { return aaArray_.size(); }

	string makeLocationString(const char* location) const;

	bool checkAaPtr(const char* ptr) const { return (ptr >= &aaArray_[0] && ptr < (&aaArray_[0] + aaArray_.size())); }
	void printStats() const;
	void printTags() const;


private:
	string fasta_;
	string path_;
	int partIdx_;
	unsigned int minTagLength_, maxTagLength_;

	vector<char> nameArray_;
	vector<char> aaArray_;
	
	vector<unsigned int> proteinStartPointers_;
	vector<unsigned int> proteinNames_;
	vector<unsigned int> locationArray_;

	set<TagWithLocations> allTags_;
};



void createIntersectionFasta(const vector<string>& fastaFiles,		 /// full paths to input files
							 size_t minLength,						 /// minimal length of required overlap
							 const char* outputFasta,				 /// fasta file with overlap sequences
							 size_t minNumOccurences = MAX_SIZE_T);  /// minimal number of genomes in which sequence needs to appear

#endif



