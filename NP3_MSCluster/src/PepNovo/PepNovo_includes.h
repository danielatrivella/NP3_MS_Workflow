#ifndef __PEPNOVO_INCLUDES_H__
#define __PEPNOVO_INCLUDES_H__

/*! \file PepNovo_includes.h
	\brief Holds various definitions (typdefs and constants) that are used throughout PepNovo's implementation.
*/

#include "../Common/includes.h"

/*! \enum AminoAcids
	This indexes the different types of standard amino acids that are handeled.
*/
typedef enum AminoAcids {N_TERM, C_TERM, Gap,Xle,Ala,Arg,Asn,Asp,Cys,Gln,Glu,Gly,His,
						 Ile,Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val} AminoAcids;

/*! \enum InputFileTypes_IFT
	This indexes the different types of input file types that can be accepted.
*/
typedef enum InputFileTypes_IFT {
	IFT_DTA, IFT_MGF, IFT_MZXML, IFT_DAT, IFT_PKL, IFT_MS2, IFT_TXT, IFT_ZIP, IFT_NUM_FILE_TYPES
} InputFileTypes_IFT;

// Data types for common variables

typedef float  mass_t;
typedef float  intensity_t;
typedef float  score_t;
typedef unsigned int   clusterIdx_t; // used mostly by MsCluster but needed here because used in DatFile
typedef long long int  longInt8_t;  // used mostly by MsCluster but needed here because used in DatFile

typedef map< string , int , less<string> > STRING2INT_MAP;
typedef map< mass_t , mass_t, less<mass_t> > MASS_T_MAP;
typedef map< int , int , less<int> > INT_MAP;
typedef map< mass_t , int   , less<mass_t> > MASS_T2INT_MAP;

const mass_t MASS_PROTON   = 1.00728;
const mass_t MASS_ELECTRON = 0.00055;
const mass_t MASS_ISO      = 1.0033;	 // mass (13C) - mass (12C)
const mass_t MASS_H2O	   = 18.01056; 
const mass_t MASS_NH3	   = 17.02655;
const mass_t MASS_CO	   = 27.99492;
const mass_t MASS_H2ONH3   = 35.03711;
const mass_t MASS_H2OH2O   = 36.02112;
const mass_t MASS_NH3NH3   = 34.0531;
const mass_t MASS_OHHH     = 19.0184;

// NP3 GOT - change NUM_SIG_DIGITS from 3 to 5
const int NUM_SIG_DIGITS = 5;
const int ESI_MASS_SPEC  = 1;


#endif

