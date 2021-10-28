#include "mlauxfun.h"
#include <string.h>


// condenses the vector so each value appears only once
// assumes vector is sorted
void condenseValWeightVector(vector<ValWeight>& vec, bool needToSort)
{
	if (needToSort)
		sort(vec.begin(),vec.end());

	// condense list
	size_t fPos=1; // forward position
	size_t bPos=0; // bacwards position

	while (fPos<vec.size())
	{
		while (fPos<vec.size() &&
			   vec[fPos].value == vec[bPos].value)
			vec[bPos].weight += vec[fPos++].weight;
			
		if (fPos == vec.size())
			break;
			
		vec[++bPos]=vec[fPos++];
	}

	while (vec.size()>bPos)
		vec.pop_back();
}




// use doubles to accumulate the results because doubles have epsilon 1E-9
// versus a float's epsilon of 1E-5
double	computeDotProduct(const vector<IdxVal>& a, const vector<IdxVal>& b)
{	
	double dotProd=0.0;
	vector<IdxVal>::const_iterator it_a = a.begin();
	vector<IdxVal>::const_iterator it_b = b.begin();

	while (it_a != a.end() && it_b != b.end())
	{
		if (it_a->value < it_b->value)
		{	
			it_a++;
			continue;
		}

		if(it_a->value > it_b->value)
		{
			it_b++;
			continue;
		}

		dotProd+= static_cast<double>((it_a->value)*(it_b->value));
		++it_a;
		++it_b;
	}
	return dotProd;
}





ListStorage::ListStorage(size_t k)
{
	if (! initialize(k) && ! initialize(k/4))
		error("Couldn't initialize a storage area (or even a quarter of it), size in K = ",k);
}

ListStorage::~ListStorage()
{
	if (storageBegin_ && storageEnd_)
	{
		size_t size = storageEnd_ - storageBegin_;
		memset(storageBegin_,0,size*sizeof(size_t)); // invalidate the lists in case some one is still looking
		delete [] storageBegin_;
	}
}

bool ListStorage::initialize(size_t k)
{
	size_t size = k *1024;

	storageBegin_ = new size_t[size];
	if (! storageBegin_)
		return false;

	storageEnd_ = storageBegin_ + size;
	nextEmpty_	= storageBegin_;

	numWrapAround_ =0;
	numAdditions_ =0;
	headerMap_.clear();

	return true;
}


bool ListStorage::addList(const ListHeader& header)
{
	size_t neededSpace = header.length + 3;

	if (nextEmpty_ + neededSpace >= storageEnd_)
	{
		numWrapAround_++;
		nextEmpty_=0;
	}

	if (nextEmpty_ + neededSpace >= storageEnd_)
		return false;

	// write list
	*nextEmpty_++ = listHeaderStartMarker;
	*nextEmpty_++ = numAdditions_;
	

	memcpy(nextEmpty_, header.list, header.length * sizeof(size_t));

	ListHeader newHeader;
	newHeader.listNumber = numAdditions_;
	newHeader.length  = header.length;
	newHeader.list	  = nextEmpty_;
	headerMap_[newHeader]=numAdditions_++;

	nextEmpty_ += header.length;

	return true;
}

bool ListStorage::getList(ListHeader& header) const
{
	map<ListHeader,size_t>::const_iterator it = headerMap_.find(header);
	if (it == headerMap_.end() || ! it->first.isValid())
		return false;

	header = it->first;

	return true;
}








