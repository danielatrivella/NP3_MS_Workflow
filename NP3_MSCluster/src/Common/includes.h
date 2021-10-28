#ifndef __INCLUDES_H__
#define __INCLUDES_H__


#ifdef WIN32
#pragma warning (disable:4996)
#pragma warning (disable:4018)
#pragma warning (disable:4244)
#pragma warning (disable:4267)
#pragma warning (disable:4305)
#pragma warning (disable:4503)
#endif


#include <limits>
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>
#include <time.h>
#include <assert.h>



using namespace std;
using std::vector;

#if defined(WIN32) && defined(_DEBUG)
	// Improved vector that checks access violations - has overhead use only in DEBUG
	template<class T, class A = allocator<T> >
	class my_vector : public vector<T,A>
	{
	public:
		explicit my_vector(const A& al = A()) : vector<T,A> (al) {}
		explicit my_vector(size_type n, const T& v = T(), const A& al = A()) : vector<T,A>(n, v, al) {}
		
		my_vector(const my_vector& x) : vector<T,A>(x) {}
		my_vector(const_iterator first, const_iterator last, const A& al = A()) : vector<T,A>(first, last, al) {}

		const_reference operator[](size_type pos) const
		{
			return at(pos);
		}
		
		reference operator[](size_type pos)
		{
			return at(pos);
		}
	};
	#define my_vector vector
#endif // defined(WIN32) && defined(_DEBUG)

typedef float		  value_t;
typedef double		  weight_t;

const unsigned long MAX_ULONG = numeric_limits<unsigned long>::max();
const unsigned int  MAX_UINT  = numeric_limits<unsigned int>::max();
const int			MAX_INT   = numeric_limits<int>::max();
const int			MIN_INT   = numeric_limits<int>::min();
const float			MAX_FLOAT = numeric_limits<float>::max();
const float			MIN_FLOAT = -numeric_limits<float>::max();
const double		MAX_DOUBLE = numeric_limits<double>::max();
const double		MIN_DOUBLE = -numeric_limits<double>::max();
const size_t		MAX_SIZE_T = numeric_limits<size_t>::max();
const value_t		MAX_VALUE_T = numeric_limits<value_t>::max();
const value_t		MIN_VALUE_T = numeric_limits<value_t>::min();
const weight_t		MAX_WEIGHT_T = numeric_limits<weight_t>::max();
const weight_t		MIN_WEIGHT_T = numeric_limits<weight_t>::min();

const float NON_FLOAT   (static_cast<float>(1.1231231e36) - numeric_limits<float>::max()); // my float value that is considered 

#define BYTEORDER_LITTLE_ENDIAN

#define NEG_INF -999999999
#define POS_INF  999999999





#endif


