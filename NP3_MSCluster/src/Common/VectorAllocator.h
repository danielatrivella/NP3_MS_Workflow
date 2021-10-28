#ifndef __VECTORALLOCATOR_H__
#define __VECTORALLOCATOR_H__

#include "includes.h"

/*********************************************************************************
This class centrally allocates vectors that can change in size.
The scenario is is designed to handle is the case where a program has millions of
vectors that continually change in size (though the total length of all vectors
remains more or less the same). The container allocates and reuses vectors of different
sizes (increase in size with powers of 2). By reusing the allocated memory, we can
avoid many of the alloc and delete operations. Any memory that is allocated can only
be freed when the whole container is destriyed.
**********************************************************************************/

//class MyVectorAllocator;

template<class T>
class MyVectorAllocator;

template<class T>
class AllocatedVector {
	friend class MyVectorAllocator<T>;
public:
	
	const T* const getElementsPtr() const { return elements_; }
	size_t getSize() const { return size_; }
	void   setSize(size_t t) { size_ = t; }
	unsigned int getCapcityExp() const { return capacityExp_; }

	void relinquish()
	{
		allocator_->reclaimVector(this);
	}

	void   addElement(const T& e)
	{
		if (capacityExp_ == 0 || size_ >= (1<<capacityExp_) )
			allocator_->moveToNewArray(this, capacityExp_+1);

		elements_[size_++]=e;
	}

	void  shrinkIfNeeded()
	{
		if (size_ + size_ < (1<<capacityExp_))
		{
			int newSizeExp = 0;
			while ((size_ >> newSizeExp) != 0) 
			{
				newSizeExp++;
			}
			allocator_->moveToNewArray(this, newSizeExp);
		}
	}

	// moves all elements to the left so the element in idx becomes element 0
	// if shrink == true, might copy to new array
	void shunt(size_t idx, bool shrink = false)
	{
		copy(elements_+idx, elements_+size_, elements_);
		size_ = size_ - idx;
		if (shrink)
			shrinkIfNeeded();
	}


	// performs a binary search, returns first position in the element
	// vector that has an element that is greater or equal to the value of e
	size_t getFirstPositionGreaterOrEqual(const T& e) const
	{
		if (size_== 0 || e > elements_[size_-1])
			return MAX_SIZE_T;
		if (e<elements_[0])	
			return 0;

		const size_t dis = distance(elements_, lower_bound(elements_, elements_ + size_, e));
		assert( elements_[dis]>=e);
		return dis;
	}

	void print() const
	{
		cout << size_ << "/" << (1<<capacityExp_);
		for (size_t i=0; i<size_ && i<10; i++)
			cout << "\t" << elements_[i];
		if (size_>10)
			cout << " ...";
		cout << endl;
	}
private:
	T* elements_;
	size_t size_;
	unsigned int capacityExp_;

	MyVectorAllocator<T>* allocator_;

	AllocatedVector() :
		elements_(0), size_(0), capacityExp_(0), allocator_(0) {}
};


template<class T>
class MyVectorAllocator {
public:

	MyVectorAllocator<T>(unsigned int minSizeExp = 1, size_t memBlockSize = 1<<20, size_t numVectorsInBlock = 1<<16) :
	  minSizeExp_(minSizeExp), memBlockSize_(memBlockSize), numVectorsInBlock_(numVectorsInBlock) { availableArrays_.resize(128);}

	~MyVectorAllocator<T>()
	{
		for (size_t i=0; i<arrayMemBlocks_.size(); i++)
			delete [] arrayMemBlocks_[i];

		for (size_t i=0; i<vectorMemBlocks_.size(); i++)
			delete [] vectorMemBlocks_[i];
	}

	AllocatedVector<T>* allocateVector() 
	{ 
		if (availableVectors_.size() == 0)
			allocateMoreVectors(); 
		AllocatedVector<T>* v = availableVectors_.back();
		availableVectors_.pop_back();
		v->elements_ = 0;
		v->size_		= 0;
		v->capacityExp_ = 0;
		return v;
	}

	AllocatedVector<T>* allocateVector(size_t size) 
	{
		AllocatedVector<T>* v=allocateVector();
		size_t sizeExp = 0;
		while ((size >> sizeExp) != 0) 
		{
			sizeExp++;
		}
		v->elements_ = getArray(sizeExp);
		v->size_     = 0;
		v->capacityExp_ = sizeExp;
		return v;
	}

	void moveToNewArray(AllocatedVector<T>* v,  unsigned int newSizeExp)
	{
		if (newSizeExp < minSizeExp_)
			newSizeExp = minSizeExp_;

		T* newArray = getArray(newSizeExp);
		if (v->size_>0)
			copy(v->elements_,  v->elements_ + v->size_, newArray);
		if (v->elements_)
			availableArrays_[v->capacityExp_].push_back(v->elements_); // return this array
		v->elements_ = newArray;
		v->capacityExp_ = newSizeExp;
	}

	void reclaimVector(AllocatedVector<T>* v)
	{
		if (v->elements_)
			availableArrays_[v->capacityExp_].push_back(v->elements_);
		availableVectors_.push_back(v);
	}


	void printReport() const
	{
		cout << "Allocated " << arrayMemBlocks_.size() << " blocks of memory = " << arrayMemBlocks_.size()*memBlockSize_ << " bytes." << endl;
		cout << "Allocated " << vectorMemBlocks_.size() << " vector memory blocks = " << numVectorsInBlock_*vectorMemBlocks_.size() << " vectors." << endl;
		cout << "Available vectors: " << availableVectors_.size() << endl;
		cout << "Available arrays : " << endl;
		for (size_t i=0; i<availableArrays_.size(); i++)
			if (availableArrays_[i].size()>0)
			{
				cout << i << "\t" << availableArrays_[i].size() << "\tblocks" << endl;
			}
		cout << endl;
	}
	
private:

	unsigned int minSizeExp_;
	size_t memBlockSize_;
	size_t numVectorsInBlock_; // size of batch of AllocatedVectors that are created with each memory allocation 

	vector<T*> arrayMemBlocks_;			// bulk memory is allocated with new
	vector< vector< T* > > availableArrays_; // first dim is size exponent, for each size t holds pointers
											// to buffers of 2^t elements of type T
	
	vector< AllocatedVector<T>* > vectorMemBlocks_;  // bulk memory is allocated with new
	vector< AllocatedVector<T>* > availableVectors_; // holds pointers to individual AllocatedVector that can be given

	void allocateMoreVectors()
	{
		AllocatedVector<T>* p = new AllocatedVector<T>[numVectorsInBlock_];
		if (! p)
		{
			cout << "Error: out of memory in MyVectorAllocator, trying to allocate " <<  sizeof(AllocatedVector<T>)*numVectorsInBlock_ <<
				" bytes for " << numVectorsInBlock_ << " vector blocks." << endl;
			exit(1);
		}
		
		vectorMemBlocks_.push_back(p);
		availableVectors_.reserve(availableVectors_.size() + numVectorsInBlock_);

		for (size_t i=0; i<numVectorsInBlock_; i++)
		{
			p[i].allocator_ = this;
			availableVectors_.push_back(&p[i]);
		}
	}

	void allocateMoreArrays(size_t sizeExp)
	{
		const size_t numElementsPerArray = 1<<sizeExp;
		const size_t bytesPerArray = numElementsPerArray*sizeof(T);

		if (2*bytesPerArray>memBlockSize_) // single array allocation
		{
			T* p = new T[numElementsPerArray];
			if (! p)
			{
				cout << "Error: out of memory in MyVectorAllocator, trying to allocate " <<  bytesPerArray <<
				" bytes for " << numElementsPerArray << " elements." << endl;
				exit(1);
			}
			availableArrays_[sizeExp].push_back(p);
			arrayMemBlocks_.push_back(p);

		//	cout << "ALLOC single " << numElementsPerArray << " (" << numElementsPerArray*sizeof(T) << ")" << endl;
			return;
		}
		else
		{
			const size_t numArrays = memBlockSize_/(numElementsPerArray*sizeof(T));
			T* p = new T[numElementsPerArray*numArrays];
			if (! p)
			{
				cout << "Error: out of memory in MyVectorAllocator, trying to allocate " <<  numArrays*bytesPerArray <<
				" bytes for " << numElementsPerArray*numArrays << " elements." << endl;
				exit(1);
			}

			availableArrays_[sizeExp].reserve(availableArrays_[sizeExp].size()+numArrays);
			for (size_t i=0; i<numArrays; i++)
				availableArrays_[sizeExp].push_back(&p[i*numElementsPerArray]);
			arrayMemBlocks_.push_back(p);

		//	cout << "ALLOC " << numArrays << "*" << numElementsPerArray << " (" <<
		//		numElementsPerArray*numArrays*sizeof(T) << ")" << endl;
			return;
		}
	}

	T* getArray(size_t sizeExp)
	{
		if (availableArrays_[sizeExp].size() == 0)
			allocateMoreArrays(sizeExp);

		T* p = availableArrays_[sizeExp].back();
		availableArrays_[sizeExp].pop_back();
		return p;
	}


};


#endif




