#ifndef __MLAUXFUN_H__
#define __MLAUXFUN_H__

#include "../Common/includes.h"
#include "../Common/auxfun.h"


// pair of two features
// assumed but not inforced: idx1<=idx2
struct IdxPair {
	IdxPair() : idx1(0), idx2(0) {}
	IdxPair(size_t i1, size_t i2) : idx1(i1), idx2(i2) {}
	bool operator< (const IdxPair& rhs) const
	{
		return (idx1<rhs.idx1 || (idx1 == rhs.idx1 && idx2<rhs.idx2));
	}
	size_t idx1, idx2;
};

struct IdxVal {
	IdxVal() : index(0), value(0.0) {}
	IdxVal(size_t i, value_t v) : index(i), value(v) {}

	bool operator< (const IdxVal& rhs) const
	{
		return (index<rhs.index || (index == rhs.index && value<rhs.value));
	}

	size_t index;
	value_t value;
};

struct ValWeight {
	ValWeight() : value(NON_FLOAT), weight(0.0) {};
	ValWeight(value_t v, weight_t w) : value(v), weight(w) {}

	bool operator< ( const ValWeight& rhs) const
	{
		return (value< rhs.value);
	}

	value_t value;
	weight_t weight;
};





// condenses the vector so each value appears only once
void condenseValWeightVector(vector<ValWeight>& vec, bool needToSort);

const size_t listHeaderStartMarker = 4040404040U; // marks the start of a list (for debug and memory corruption tests)

// Is a pointer to a (usually sorted) list of feature/sample indexes
struct ListHeader {
	bool operator< (const ListHeader& rhs) const
	{
		return (listNumber < rhs.listNumber);
	}

	bool isValid() const { return (*(list-2)== listHeaderStartMarker &&
								   *(list-1)== listNumber); }

	size_t  listNumber;
	size_t  length;
	size_t* list;
};

// Holds a large memory block where lists can be stored.
// Lists are identified by the pair index. Once a list is added
// it might not stay for ever, it might get over written, in which
// case we will nedd to compute it again.

class ListStorage {
public:
	ListStorage() : numWrapAround_(0), numAdditions_(0), storageBegin_(0), 
		storageEnd_(0),  nextEmpty_(0) {}
	ListStorage(size_t k);
	~ListStorage();
								
	bool initialize(size_t k); // size in k=1024 bytes

	// will copy the header info if a valid record is found, otherwise
	// returns NULL. Note that the header returned is not guaranteed to
	// be valid when used if between the time it was given and the time it
	// is accessed, there have been several additional insertions (if the
	// storage area is large, this should not be a problem, use judgement).
	bool getList(ListHeader& header) const;

	// doesn't check if the IdxPair is already there, if it is, it will
	// overwrite the list. Might return false if the list could not be written
	// (maybe too large for storage)
	bool addList(const ListHeader& header);

	void clear() { headerMap_.clear(); numWrapAround_=0; numAdditions_=0; nextEmpty_ = storageBegin_; }

private:
	
	size_t  numWrapAround_;
	size_t	numAdditions_;
	size_t* storageBegin_;
	size_t*	storageEnd_;
	size_t* nextEmpty_;
	map<ListHeader,size_t> headerMap_;
};


double	computeDotProduct(const vector<IdxVal>& a, const vector<IdxVal>& b);

template<class T>
double computeDotProduct(const vector<IdxVal>& a, const vector<T>& b)
{
	assert( a.size()==0 ||  (a[a.size()-1].index < b.size()) );
	double dotProd=0.0;
	for (size_t i=0; i<a.size(); i++)
		dotProd+= static_cast<double>(a[i].value * b[a[i].index]);
	
	return dotProd;
}

template<class T>
double computeDotProduct(const vector<T>& b, const vector<IdxVal>& a)
{
	assert( a.size()==0 ||  (a[a.size()-1].index < b.size()) );
	double dotProd=0.0;
	for (size_t i=0; i<a.size(); i++)
		dotProd+= static_cast<double>(a[i].value * b[a[i].index]);
	
	return dotProd;
}


template<class S, class T>
double computeDotProduct(const vector<S>& a, const vector<T>& b)
{
	assert(a.size() == b.size());
	double dotProd=0.0;
	for (size_t i=0; i<a.size(); i++)
		dotProd+= static_cast<double>(a[i])*static_cast<double>(b[i]);

	return dotProd;
}


template<class S, class T>
double computeDotProduct(const vector< vector<S> >& a, const vector< vector<T> >& b)
{
	assert(a.size() == b.size());
	double dotProd=0.0;
	for (size_t c=0; c<a.size(); c++)
		dotProd+= computeDotProduct(a[c],b[c]);

	return dotProd;
}

template<class T>
void multiplyByScalar(vector<T>& v, T scalar)
{
	for (size_t i=0; i<v.size(); i++)
		v[i]*=scalar;
}

template<class T>
void multiplyByScalar(vector< vector<T> >& v, T scalar)
{
	for (size_t c=0; c<v.size(); c++)
		for (size_t i=0; i<v.size(); i++)
			v[c][i]*=scalar;
}



#endif 
