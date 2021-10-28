#ifndef __MSCLUSTERBENCHMARK_H__
#define __MSCLUSTERBENCHMARK_H__

/*! @file MsClusterBenchmark.h
	\breif a wrapper for all benchmark experiments */

class AllScoreModels;
struct MsParameterStruct;

/*!	A wrapper for all benchmark experiments. All parameters are passed through the struct \c MsParameterStruct* \c params
@param model The scoring model (class containing all relevant models)
@param params The job parameters
*/
void performBenchmark(AllScoreModels* model, MsParameterStruct* params);


#endif

