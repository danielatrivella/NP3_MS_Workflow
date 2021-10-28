#ifndef __MLLIB_H__
#define __MLLIB_H__

#include "mlmodel.h"



void printUsage(char *message=NULL);

double classifyDataFile(const char* modelPath, const char* dataFile, bool verbose=false);

#endif

