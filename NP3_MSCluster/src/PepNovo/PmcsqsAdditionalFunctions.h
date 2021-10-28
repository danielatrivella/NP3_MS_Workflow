#ifndef __PMCSQSADDITONALFUNCTION_H__
#define __PMCSQSADDITONALFUNCTION_H__

/*! \file  PmcsqsAdditionalFunctions.h
	\brief This file contains additional functions that are not class members but
	involve the SQS or PMC models, such as filtering datasets, selecting training
	samples from datasets, etc.
*/

#include "AllScoreModels.h"

/*!
\brief A function that runs quality filtering on a dataset and produces a new set of MGF
files containing the spectra that passed the quality threshold.

The input files are read one at a time. All spectra in each file are read and the SQS
is computed for each one. Spectra with an SQS below the threshold are discarded. An 
MGF file is created with all the remaining spectra (the name is the same as the input
file, but without the extension (which is replaced with "_fil.mgf").

@param model	Pointer to the model container class (AllScoreModels).
@param list		Path to text file with full paths of input spectra files.
@param outDir	Path to directory where output files should be written.
@param sqsThreshold	The quality threshold that is used to filter, takes vlaues
between 0-1.0 (a typical vlaue is 0.1).
*/
void filterDataSet(AllScoreModels* model, const string& list, const string& outDir, 
				   float sqsThreshold, const string& suffix);



void selectSpectra(AllScoreModels* model, 
				   const string& list, 
				   const string& outDir,
				   int selectCharge);

/*!
\fn makeCrapDb(AllScoreModels* model, const string& list, const string& name, 
				const string& outDir, float sqsThresh, int numCrapSpectra)
\brief Creates a set of low quality spectra (useful for SQS model training)

This function is used for creating a better SQS training set. Once you somehow selected
a set of crap spectra from your training data, you can create an SQS model. However this
model might not be that great since there where many peptide containing spectra in your
"crap" set. Use this function to create a new set of "crap" spectra and retrain the SQS model.

@param model	Pointer to score models container (must include a valid SQS model).
@param list		Path to text file with a list of paths to spectra files.
@param name		The prefix of the file names of the created mgf files.
@param outDir	Path to directory where mgf files will be written.
@param sqsThresh The maximal SQS value for a spectrum in order for it to be written as crap.
@param numCrapSpectra The maximal number of spectra that will be written into mgf files.
*/
void makeCrapDb(AllScoreModels* model, const string& list, const string& name, 
				const string& outDir, float sqsThresh, int numCrapSpectra);



/*!
\fn benchmarkSqs
\brief Computes the spectra qulaity score for a list of spectra files and creates a histogram.
If the spectra have their charges specified, it also checks charge assignment accuracy.

@param model  Pointer to the model container class (AllScoreModels); must include a valid SQS model.
@param list		Path to text file with a list of paths to spectra files.
*/
void benchmarkSqs(AllScoreModels* model, const string& list);

#endif

