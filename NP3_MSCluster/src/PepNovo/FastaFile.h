#ifndef __FASTAFILE_H__
#define __FASTAFILE_H__


/*! \file FastaFile.h
	\brief holds the class that reads a fasta file and stores it (with tags)
*/

#include "../Common/includes.h"

struct TagWithLocations {
	TagWithLocations() : tagLength(0), numLocations(0), tag(0), locations(0) {}

	bool operator< (const TagWithLocations& rhs)
	{
		if (tagLength < rhs.tagLength)
			return true;
		if (tagLength > rhs.tagLength)
			return false;
		return (strncmp(tag, rhs.tag, tagLength)<0);
	}

	size_t tagLength;	   /// number of amino acids (no terminating \0)
	size_t numLocations; /// how many places this tag appears in the fasta file
	char* tag;				   /// the tag itself (pointer to location in some char array)
	char* locations;		   /// pointer to location of list of locations in char* array
	
};


class TagStorage {
public:
	TagStorage() : tagLength_(0), numTags_(0), numLocations_(0) {}



private:
	size_t numTags_;
	size_t numLocations_;
	
	vector<TagWithLocations> allTags_;
	vector<char> tagCharArray_;
	vector<char*> locationArray_;
};


class FastaFile {
public:

private:
	

};


#endif

