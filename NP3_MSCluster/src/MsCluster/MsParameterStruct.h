#ifndef __MSPARAMETERSTRUCT_H__
#define __MSPARAMETERSTRUCT_H__

/*! \file MsParameterStruct.h */

#include "../PepNovo/PepNovo_includes.h"
#include "MsClusterIncludes.h"

/*! Contains all the parameters involved in running a clustering/archive/library jobs.*/
struct MsParameterStruct {

	MsParameterStruct() : commandLine(std::string()), modelName("LTQ_TRYP"), modelDir("Models"), 
		outDir("out"), tmpDir("tmp"), archiveOutDir(std::string()),
		list(std::string()), metaList(std::string()), datList(std::string()), exclusionList(std::string()), 
		inclusionList(std::string()), ptmString(std::string()), outputName(std::string()), 
		outputStub(std::string()), outputNameWithVersion(std::string()), annFile(std::string()), archivePath(std::string()),
		idListPath(std::string()), clusterList(std::string()), nonClusterList(std::string()), spectraListToLoad(std::string()),
		spectraListToTest(std::string()), inputFile(std::string()), intersectionFastaPath(std::string()),
		pathToArchive1(std::string()), pathToArchive2(std::string()),
		minMz(0.0), maxMz(MAX_FLOAT), precursorPPMs(MAX_FLOAT), fragmentTolrance(0.0), rtTolerance(5.0), minNumPeaks(1), scaleFactor(0.0),
		sqsThreshold(0.0), minSimilarity(0.2), maxMixtureProb(0.1), memoryGb(0.5), mzWindow(2.0), 
		datasetIdx(MAX_INT), batchIdx(0), startFileIdx(0), numRounds(5), peakDensity(-1), 
		majorIncrement(MAJOR_MZ_INCREMENT_FOR_DAT), numCrapSpectra(-1), simType(0),
		verboseLevel(1), fileConversionType(0), outputFileSize(20000), minIntersectionLength(7), minNumberForFastaIntersection(MAX_SIZE_T),
		maxResultsForSpectrum(10), maxPvalue(0.05),
		gotAllowArchiveCombine(false), gotAssignCharges(false), gotDatOnly(false), gotKeepDat(false), gotFirstPass(false), gotSecondPass(false), 
		gotBenchmarkPairs(false), gotBenchmarkSimilarityHistogram(false), gotMakeMgf(false), gotBenchmarkSimilarity(false), gotBenchmarkConsensus(false), 
		gotBenchmarkSqs(false), gotBenchmarkLibrary(false), gotFilterOnly(false), gotMakeCrapDb(false),  gotCreateArchive(false), 
		gotCreateArchiveFromMgfs(false), gotCompleteArchiveSearch(false), gotArchiveStatistics(false), gotOutputMgf(false), gotOutputDat(false), 
		gotOutputModified(false), gotLearnLibraryPValues(false), gotCreateLibrary(false), gotCreateLibraryFromDat(false),
		gotLibraryStats(false), gotMakeLibraryWithPeptidesOnly(false), gotNoSqs(false), gotConvert(false), gotCnvertArchive(false),
		gotUseInputTitles(false), gotCorrectPM(false), gotOverwriteLocations(false), gotMakeExclusion(false), gotBenchmarkPm(false), 
		gotTrainSimilarityModel(false), gotMergeArchives(false), gotSearchArchive(false), gotOutputMatchedSpectra(false) {}

	void printParameters(ostream& os = cout) const;
	bool checkIfGotBenchmark() const { return (gotBenchmarkSimilarity || gotBenchmarkConsensus || 
											   gotBenchmarkSqs || gotBenchmarkLibrary || gotBenchmarkPairs ||
											   gotBenchmarkPm || gotBenchmarkSimilarityHistogram); }
	void copyCommandLine(int argc, char** argv);

	string commandLine; /// The command line entered for the run

	string modelName; /// The model used for filtering etc., e.g., LTQ_TRYP

	string modelDir; /// The directory in which the model files sit (if not ./Models)

	string outDir;   /// Where to write the clustering outout files (e.g., *.mgf *.clust)

	string tmpDir;	 /// Where to write the binary dat files with the spectra

	string archiveOutDir; /// Where to write the archive files

	string list;	 /// The list of paths to the spectra input files

	string metaList; /// The list of the list of paths

	string datList;	 /// The list of paths to the dat input files

	string exclusionList; /// The list of paths that hold scans that should be excluded

	string inclusionList; /// list of scans that should be included

	string ptmString; /*! The string designating the PTMs that should be considered (needed if 
						  the spectra are annotated with peptides that have modifications, unless using
						  an archive which specifies all the PTMs that are included.*/

	string outputName;	/// The name given as a prefix to all generated files

	string outputStub;

	string outputNameWithVersion; /// name + "_" + generation + "_" + batch

	string annFile;   	/*! Used for creating annotated mgf files from dat (not likely to be used by oridinary users) */

	string archivePath; /// The path to the main archive file

	string idListPath; /// Path to text file with paths of id files for clusters

	string clusterList, nonClusterList; /// used for benchmakring

	string spectraListToLoad, spectraListToTest; /// used for benchmarking library identifications

	string inputFile; /// path to input file

	string intersectionFastaPath; /// path of file that is the output of the intersection

	string pathToArchive1; /// Path to archive 1 (when merging two archives)

	string pathToArchive2; /// Path to archive 2 (when merging two archives)

	mass_t minMz; /// minimal m/z of spectra to load/process 
	
	mass_t maxMz; /// maximal m/z of spectra to load/process

	float precursorPPMs; /// The ppm tolerance that should be used for high res data

	float fragmentTolrance; /// The tolerance in Da used for fragments

	// NP3 GOT RT tol
	float rtTolerance; // the tolerance in seconds for retention time - clusters with rt abs diff greater than this wont be joined
	int minNumPeaks; // min number of peaks that a cluster must have to be outputted >=1
	float scaleFactor; //the number to scale the fragmented peaks intensities before any dotproduct comparision.
	                    // Valid values: 0 for the natural logarithm (ln) of the intensities; 1 for no scaling;
	                    // and other values of x will raise the intensities to the power of x (e.g. x = 0.5 is the square root scaling). [x] >= 0

	float sqsThreshold; /// The SQS quality threshold for filtering out bad spectra

	float minSimilarity;/// The minimal similarity threshold for joining spectra into a cluster

	double maxMixtureProb; /// The maximal probability allowed for wrongfully joining  spectrum to a cluster

	float memoryGb;	/*! The number of GB of RAM the program can use. This parameter
							is an approximation, the executable might actually use a bit more 
							than allocated */

	mass_t mzWindow;	/*! The m/z window used in the clustering process. The program will
						    only consider joining spectra that whose precursor m/z differs
						    by this value or less. Note that in some conditions such as when 
						    the memory allocated is not large ebough, the actual window being 
						    used might be smaller than this value. */

	int datasetIdx;	/*! In cases where we are adding new data to pervious clustering jobs,
						    this parameters tells the current jobs generation. This parameter will
						    get added to the name of the clusters (in the titles ) and the 
						    outputted file names.*/

	int	batchIdx;		/*! In cases of clustering large sets of input files, the clustering stage
						    can be split into several batches based on the m/z values (as they
						    are reflected in the dat file name e.g., data_0_1000_0.dat). In such cases
						    we can run different m/z values in batches. The batch flag is used to 
						    avoid giving clusters non-unique names in their titles. */

	size_t startFileIdx;	/*! The index of the first input file. This can be non-zero in a case 
							    where a large set of input files was split into batches (and dat
							    files are created separately for each batch). */

	int numRounds;		/*! number of rounds in which the algorithm performs clustering (using a
						    decreasing similarity threshold that reaches the similarity parameter
						    in the last round. */

	int peakDensity;   /// the number of peaks to keep in a 200 Da window, if set to a value>0, should not be lower than 15

	mass_t majorIncrement; /*! the "slice" size for creating the first pass dat files and the output files.
						  Default is defined in \c MAJOR_MZ_INCREMENT_FOR_DAT */

	int numCrapSpectra; /// for creating SQS training sets

	int simType;   /// type of similarity function to use (for benchmarking)

	int verboseLevel;	/// how noisy should the program be (value>0 gives more messages)

	int fileConversionType; /// when converting files, what type should we write them as (IFT_MGF or IFT_DAT)

	size_t outputFileSize; /// number of spectra to write in each ouput file (dat/mgf/clust)

	size_t minIntersectionLength; /// the minimal length of overlap between fasta files (default X=7)

	size_t minNumberForFastaIntersection; /// the minimal number of fasta files that need to include  a segement for it to be kept in the intersection files (default X=#fasta files)

	int maxResultsForSpectrum; // when searhcing archive

	float maxPvalue;	// when searhcing archive

	// indicators on the status of the input parameters
	bool gotAllowArchiveCombine;	/// Indicator that tells if allow the clustering to compare archived spectra
	bool gotAssignCharges;
	bool gotDatOnly;
	bool gotKeepDat;
	bool gotFirstPass;
	bool gotSecondPass;
	bool gotMakeMgf;
	bool gotBenchmarkSimilarity;
	bool gotBenchmarkConsensus;
	bool gotBenchmarkSqs;
	bool gotBenchmarkLibrary;
	bool gotBenchmarkPairs;
	bool gotBenchmarkSimilarityHistogram;
	bool gotFilterOnly;
	bool gotMakeCrapDb;
	bool gotCreateArchive;
	bool gotCreateArchiveFromMgfs;
	bool gotCompleteArchiveSearch; // both ided and unided clusters in library search
	bool gotArchiveStatistics;
	bool gotOutputMgf;
	bool gotOutputDat;
	bool gotOutputModified;
	bool gotLearnLibraryPValues;
	bool gotCreateLibrary;
	bool gotCreateLibraryFromDat;
	bool gotLibraryStats;
	bool gotMakeLibraryWithPeptidesOnly; /// If true extracts only spectra that are annotated with a peptide
	bool gotNoSqs;
	bool gotConvert; /// file conversion to dat or mgf
	bool gotCnvertArchive; // convert archive to mgf
	bool gotUseInputTitles; /// indicator that the spectrum titles should be taken from the input spectra when possible.
	bool gotCorrectPM;	/// indicates if PM correction should be use (in this case the 
	bool gotOverwriteLocations; /// indicates wether the locations written in dat files should be overwritten with the locations according to genration idx/ index in list
	bool gotMakeExclusion;		
	bool gotBenchmarkPm;
	bool gotTrainSimilarityModel;
	bool gotMergeArchives;
	bool gotSearchArchive;
	bool gotOutputMatchedSpectra;

	vector<string> datPaths1, datPaths2; /// for merging archives
};




#endif


