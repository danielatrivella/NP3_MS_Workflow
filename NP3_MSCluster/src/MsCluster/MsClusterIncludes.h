#ifndef __MSCLUSTERICLUDES_H__
#define __MSCLUSTERICLUDES_H__

/*! @file MsClusterIncludes.h
	A central location for definitions and constants used in MsCluster.
*/

#include "../Common/includes.h"
#include "../Common/auxfun.h"
#include "../PepNovo/PepNovo_includes.h"
// GOT change NUM_PEAKS_FOR_HEURISTIC from 4 to 6
// change SIMILARITY_BATCH_SIZE from 800 to -
// change MAX_SIMILARITY from 0.7 to 0.9
// LARGE_CLUSTER_SIZE from 20 to 200
// ALLOWED_DIFF_IN_NUM_DISTANCE_PEAKS from 5 to 8
// MAX_NUM_PEAKS_FOR_DISTANCE from 40 to 25
// NUM_TOP_SIMILARITIES_TO_SAVE from 25 to 100
const size_t NUM_PEAKS_FOR_HEURISTIC       = 6;	 /*! Number of peaks used for deciding if two spectra should be compared.
												     If they don't have at least one peak mass in common in the top
													 \c NUM_PEAKS_FOR_HEURISTIC in common there is no reason to compare them.
													 It is not recmonneded to go below 3, and there is no benefit in going above 6.*/
/*const size_t NUM_TOP_SIMILARITIES_TO_SAVE  = 500;   /// number of distance computations to remember between rounds
const size_t MAX_NUM_PEAKS_FOR_DISTANCE    = 15;   /// Maximal number of peaks used in distance computation
const float MIN_SIMILARITY_FOR_ACCURATE_COMPUTATION = 0.1; /// minimal similarity required for accurate similarity computation (which includes peak intensities)
const size_t ALLOWED_DIFF_IN_NUM_DISTANCE_PEAKS = 12; /// maximum size difference in distance peaks that is tolerated between two spectra
const float	MAX_SIMILARITY				   = 0.85; /// do not require clusters to have a larger similarity than this for joining
const unsigned int LARGE_CLUSTER_SIZE      = 1000;  /// number of spectra for cluster to be considered "big" and rquire higher similarity to add too
const unsigned int MAX_CLUSTER_SIZE_FOR_UPDATE  = 255; /// adding spectra beyond this size does not change the peaks in the consensus spectrum.
const size_t	SIMILARITY_BATCH_SIZE	   = 3000;	/// how many specra should be clustered in each sweep
const size_t	MAX_SIMILARITY_LIST_SIZE   = 20000000; /// maximal number of pairs to store in list of pairs for similarity
const mass_t MAJOR_MZ_INCREMENT_FOR_DAT	   = 25.0; /*! output and first pass dat split files according to this value
												   of Daltions


//const unsigned int SIZES_FOR_REDOING_DISTANCE_PEAKS[]={2,3,4,5,9,16};
const unsigned int SIZES_FOR_REDOING_DISTANCE_PEAKS[]={2,3,4,5,8,15};
const size_t   NUM_SIZES_FOR_REDOING_DISTANCE_PEAKS = sizeof(SIZES_FOR_REDOING_DISTANCE_PEAKS)/sizeof(unsigned int);*/

// GOT INTENSITY BASELINE
const float INTENSITY_BASELINE = 0.0;

const size_t NUM_TOP_SIMILARITIES_TO_SAVE  = 50;   /// number of distance computations to remember between rounds
const size_t MAX_NUM_PEAKS_FOR_DISTANCE    = 30;   /// Maximal number of peaks used in distance computation
const float MIN_SIMILARITY_FOR_ACCURATE_COMPUTATION = 0.1; /// minimal similarity required for accurate similarity computation (which includes peak intensities)
const size_t ALLOWED_DIFF_IN_NUM_DISTANCE_PEAKS = 30; /// maximum size difference in distance peaks that is tolerated between two spectra
const float	MAX_SIMILARITY				   = 0.7; /// do not require clusters to have a larger similarity than this for joining
const unsigned int LARGE_CLUSTER_SIZE      = 5000;  /// number of spectra for cluster to be considered "big" and rquire higher similarity to add too
const unsigned int MAX_CLUSTER_SIZE_FOR_UPDATE   = 255; /// adding spectra beyond this size does not change the peaks in the consensus spectrum.
const size_t	SIMILARITY_BATCH_SIZE	   = 800;	/// how many specra should be clustered in each sweep
const size_t	MAX_SIMILARITY_LIST_SIZE   = 20000000; /// maximal number of pairs to store in list of pairs for similarity
const mass_t MAJOR_MZ_INCREMENT_FOR_DAT	   = 25.0; /*! output and first pass dat split files according to this value
												   of Daltions*/

const unsigned int SIZES_FOR_REDOING_DISTANCE_PEAKS[]={2,3,4,5, 7, 9, 12, 16,25, 30, 50,100, 200, 500}; //TODO porcentagem da area, quanta informação já possui, se tem 70% nao recalcula mais
const size_t   NUM_SIZES_FOR_REDOING_DISTANCE_PEAKS = sizeof(SIZES_FOR_REDOING_DISTANCE_PEAKS)/sizeof(unsigned int);

											

const clusterIdx_t    MAX_CLUSTER_IDX    = numeric_limits<clusterIdx_t>::max();
const longInt8_t	  MAX_INT8			 = numeric_limits<longInt8_t>::max();

#endif



