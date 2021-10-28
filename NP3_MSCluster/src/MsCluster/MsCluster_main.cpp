/*
% Copyright 2009, The Regents of the University of California
% All Rights Reserved
%
% Permission to use, copy, modify and distribute any part of this
% program for educational, research and non-profit purposes, without fee,
% and without a written agreement is hereby granted, provided that the
% above copyright notice, this paragraph and the following three paragraphs
% appear in all copies.
%
% Those desiring to incorporate this work into commercial
% products or use for commercial purposes should contact the Technology
% Transfer & Intellectual Property Services, University of California,
% San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910,
% Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu.
%
% IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
% FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
% INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN
% IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY
% OF SUCH DAMAGE.
%
% THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY
% OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
% ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF CALIFORNIA MAKES NO
% REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR
% EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
% THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.

*/

/*!
@file MsCluster_main.cpp 

\brief This is the main file for MsCluster it manages the command line interface and 
calls the high-level functions that perfrom actions like dat creation, 
clustering, benchmarks etc.

*/



#include "MsParameterStruct.h"
#include "MsClusterAuxfuns.h"
#include "MsClusterBenchmark.h"
#include "MsClusterAlgorithm.h"
#include "MsSimilarityModel.h"
#include "DatFileWriter.h"
#include "../PepNovo/AnnotationFile.h"
#include "../PepNovo/PmcsqsAdditionalFunctions.h"
#include "../PepNovo/DbFile.h"
#include "../Common/auxfun.h"


const string build_name = "20101018"; 


/*!
Prints the help message for the user, specifing all command line flags.
*/
void print_help(const char *message) 
{
	if (message)
		printf("***************************************************************************\n\n%s\n",message);

	printf("\nMsCluster v2.00 - Clustering algorithm for billions of spectra.\n");
	printf("Release %s.\nAll rights reserved to the Regents of the University of California.\n\n",build_name.c_str());

	
	printf("Required arguments:\n");
	printf("-------------------\n\n");
	printf("--model <model name> \n\n");
	printf("--list  <path> - path to text file listing input files>\n");
	printf("  OR\n");
	printf("--meta-list <path> - path to list of lists (each list considered a dataset)\n");
	printf("  OR\n");
	printf("--dat-list  <path> - path to text file listing dat files>\n");
	printf("--output-name <XX> - the prefix name given to all clusters and files (should include generation and batch numbers if repeated or split runs, respectively).\n");
	
	printf("\n\nOptional arguments (file I/O): \n");
	printf(    "------------------------------ \n");
	printf("--start-file-index <X> - start numbering the source files from X (useful if splitting a dat creation job).\n");
	printf("--dataset-idx      <X> - default X=0, this number is added to file names and cluster titles (to represent the number of times clustering was run on an incrementally increasing dataset.\n");
	printf("--batch-idx        <X> - default X=0, a clustering job can be split to batches run on different m/z values in parallel.\n");
	printf("--output-file-size <X> - number of clusters written per output file (default X=20000).\n");
	printf("--tmp-dir  <path>      - path to where dat files are written (default ./tmp).\n");
	printf("--out-dir  <path>      - path to where the output files (e.g., *.mgf, *.clust, archive or library) are written.\n");
	printf("--dat-only             - only create dat files, do not cluster (performs first and second passes in tandem).\n");
	printf("--model-dir  < path >  - directory where model files are kept (default ./Models)\n");
	printf("--PTMs   <PTM string>  - separated  by a colons (no spaces) e.g., M+16:S+80:N+1. Should be used if spectra are annotated with modified peptides.\n");	
	printf("--keep-dat	           - do not delete the dat files that were written in the tmp directory.\n");
	printf("--first-pass           - only perform first pass of dat writing (input should be with --list).\n");
	printf("--second-pass          - only perform the second pass of dat writing (input should be with --dat-list).\n");
	printf("--major-increment      - size of m/z slices used for first pass file gneration and ouput files (default is X=25), value must be in whole Da.\n");
	printf("--assign-charges	   - use the charges in the spectra files, otherwise all charges are set to 0.\n");
	printf("--use-input-titles     - use the titles given in the spectra as cluster titles (when possible).\n");
	printf("--output-mgf	       - write the clusters as mgf files (default if --create-archive is not used).\n");
	printf("--convert <type>       - convert the input files to another type (type can be dat or mgf).\n");

	printf("\n\nArchives:\n");
	printf(    "---------\n");
	printf("--create-archive     - Creates the files for an archive (including dat files).\n");
	printf("--convert-archive	 - Convert dat files to mgf (must supply --dat-list and --output-name)\n");
	printf("--merge-archives  <path to archive1> <path to archive2> - merges two archives (and writes them inot a third). The paths point to the main archive files \"_archive.txt\"\n");
	printf("--archive-statistics - Writes statistics about archive, (number of spectra, cluster sizees, etc.)\n");
	printf("--output-modified    - Writes clusters that are new or have been modified.\n");
	
	printf("--create-from-mgfs	 - Create an archive from mgf files, expects to have --list with mgf files, --id-list with ids (only ided mgfs will be included), and --output-name\n");
	
	printf("\n\nSearching archives:\n");
	printf(    "-------------------\n");
	printf("--search-archive  <path to archive>   - searches the files provided in --list  against the archive.\n");
	printf("--output-matched  - Writes spectra that were matched in search.\n");
	printf("--id-list         <path to list file> - paths of id files for the archive clusters, these are external user-generated files.\n");
	printf("--max-pvalue	  <X> - maximal p-value for a reported match (default X=0.05).\n");
	printf("--max-results	  <X> - maximal number results returned for each spectrum (default X=10).\n");
	printf("--complete-search     - if and id-list is provided and this flag is set all clusters are searched in the archive (the deault is to exclude unannotated clusters from the search results)\n");
	
	printf("\n\nAdditional optional arguments (clustering / archive search): \n");
	printf(    "------------------------------------------------- \n");
	printf("--memory-gb   <GG>   - the number of GB of RAM available for this process (Must allocate at least 0.35GB, do not allocate more than 3GB on a 32-bit machine).\n");
	printf("--min-mz       <X> - the minimal precursor m/z to process (default X=0)\n");
	printf("--max-mz       <X> - the maximal precursor m/z to process (default X=10000)\n");
	printf("--sqs          <X> - threshold for quality filtering 0.0<X<1.0 (typical value should be X=0.1)\n");
	printf("--fragment-tolerance <X> - the tolerance in Da for fragment peaks (default 0.34 Da.)\n");
	printf("--precursor-ppm		 <X> - the precursor mass tolerance in ppm for high-resolution data (not used with regular data).\n");
	printf("--correct-pm        - tries to correct the precursor m/z (uses correct value for clustering, but outputs the original m/z)\n");
	printf("--peak-density <X>  - for spectra preprocessing; the number of peaks to keep in a window of 200 Da (default is X=20, but higher values might be needed to avoid loss of identifications with certain database search perograms).\n");
	printf("--assign-charges    - tries to assign clusters with charges according to the charges in the spectra files (default assigns charge 0 to all spectra).\n");
	printf("--window        <X> - default X=2.0 Da. This is the window width for the m/z of the precursor that determines if two spectra will be compared and possibly joined.\n");
	printf("--mixture-prob  <X> - the probability wrongfully adding a spectrum to a cluster (default X=0.05)\n");
	printf("--num-rounds    <X> - determines how many rounds are used for the hierarchical clustering (default X=3).\n");
	printf("--keep-dataset-idx  - for large cluster jobs with pre-processed dat files.\n");
	printf("--verbose-level <X> - for values X>0 returns debug information (default X=0).\n");
	//NP3 parms
	printf("--rt-tolerance <X> - tolerances in seconds for the retention time width of the precursor that determines if two spectra will be compared and possibly joined. X>=0 (default X=0). \n");
	printf("--min-peaks-output <X> - the minimum number of fragment peaks that a spectrum must have to be outputted after the final clustering step. Spectra with less than x peaks will be discarted. x>=1 (default X=1).\n");
	printf("--scale-factor <X> - the number to scale the fragmented peaks intensities before the dotproduct comparisions. Valid values: 0 for the natural logarithm (ln) of the intensities; 1 for no scaling; and other values of x will raise the intensities to the power of x (e.g. x = 0.5 is the square root scaling). [x] >= 0 (default X=0.0).\n");
//	printf("--benchmark-similarity <CLUSTER LIST> <NON-CLUSTER LIST> - takes to text files, the first a list of mgf files for clsuters, the second, a list of mgf files of non-clusters. Creates an ROC curve.\n");
//	printf("--benchmark-pmc - runs benchmarks on the list provided (with --list).\n");
//  printf("--sim-type	<X> - similarity type (0 = hybrid, 1 = only peak count, 2 = only dot-prod)
//	printf("--intersection-fasta <X> - path to fasta that is the intersection of fasta files given with --list.\n");
//	printf("--min-intersection-length <X> - the minimal length of overlap between fasta files (default X=7).\n");
//  printf("--min-number-for-intersection <X>) - the minimal number of fasta files that need to include  a segement for it to be kept in the intersection files (default X=#fasta files).\n");
	printf("\n\nPlease send comments and bug reports to Ari Frank (arf@cs.ucsd.edu).\n\n");
#ifdef WIN32
	system("Pause");
#endif
	exit(1);
}


/*!
\brief The main function entry point.
*/
int main(int argc, char** argv)
{
	MsParameterStruct params;

	seedRandom(112233);		// use a fixed seed for "random" so behavior is the same in different runs

	if (argc <= 0)
		print_help(NULL);

#if defined(WIN32) || defined(WIN64)
	params.modelDir = "Models_Windows";
#endif
	params.copyCommandLine(argc, argv);

	int selectCharge=0;
	bool indKeepDatasetIdx = false;
	bool gotDatasetIdx     = false;
	bool gotCorrectDat	   = false;

	// This portion of the code parses the command line arguments
	// There is probably a more elegant way to do this
	for (int i=1; i<argc; i++)
	{
		if (! strcmp(argv[i],"--model"))
		{
			if (++i == argc)
				print_help("missing model name!");
			params.modelName=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--model-dir"))
		{
			if (++i == argc)
				print_help("missing model dir!");
			params.modelDir=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--file"))
		{
			if (++i == argc)
				print_help("missing path to file!");
			params.inputFile=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--list"))
		{
			if (++i == argc)
				print_help("missing path to list!");
			params.list=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--meta-list"))
		{
			if (++i == argc)
				print_help("missing path to meta-list!");
			params.metaList=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--dat-list"))
		{
			if (++i == argc)
				print_help("missing path to dat list!");
			params.datList=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--keep-dataset-idx"))
		{
			indKeepDatasetIdx = true;
			continue;
		}

		if (! strcmp(argv[i],"--exclusion-list"))
		{
			if (++i == argc)
				print_help("missing path to exclusion list!");
			params.exclusionList=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--make-exclusion"))
		{
			if (++i == argc)
				print_help("missing path to exclusion output file!");
			params.exclusionList = argv[i];
			params.gotMakeExclusion = true;
			continue;
		}

		if (! strcmp(argv[i],"--inclusion-list"))
		{
			if (++i == argc)
				print_help("missing path to inclusion list!");
			params.inclusionList=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--out-dir"))
		{
			if (++i == argc)
				print_help("missing output directory!");
			params.outDir=stripPathOfTrailingSymbols(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--tmp-dir"))
		{
			if (++i == argc)
				print_help("missing dat directory!");
			params.tmpDir=stripPathOfTrailingSymbols(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--convert"))
		{
			if (++i == argc)
				print_help("Missing type for --convert (options are dat or mgf)");

			params.gotConvert = true;
			if (! strcmp(argv[i],"mgf"))
			{ 
				params.fileConversionType = IFT_MGF;
			}
			else if (! strcmp(argv[i],"dat"))
			{
				params.fileConversionType = IFT_DAT;
			}
			else
				error("File conversion type not supported for: ",argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--create-from-mgfs"))
		{
			params.gotCreateArchiveFromMgfs = true;
			continue;
		}

		if (! strcmp(argv[i],"--complete-search"))
		{
			params.gotCompleteArchiveSearch = true;
			continue;
		}

		if (! strcmp(argv[i],"--convert-archive"))
		{
			params.gotCnvertArchive = true;
			continue;
		}

		if (! strcmp(argv[i],"--output-name"))
		{
			if (++i == argc)
				print_help("missing name to assign to clusters and files!");
			params.outputName=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--PTMs"))
		{
			if (++i == argc)
				print_help("missing PTM string!");
			params.ptmString=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--sqs"))
		{
			if (++i == argc)
				print_help("missing sqs threshold!");
			
			params.sqsThreshold = MIN_FLOAT;
			params.sqsThreshold = atof(argv[i]);
			if (params.sqsThreshold<0.0 || params.sqsThreshold>=1.0)
				print_help("--sqs <X>, the treshold should be in the range 0<X<1");
			continue;
		}

		// NP3 GOT RT PARAM
		if (! strcmp(argv[i],"--rt-tolerance"))
        {
            if (++i == argc)
            {
                print_help("missing rt tolerance value!");
            }
			params.rtTolerance = MIN_FLOAT;
			params.rtTolerance = atof(argv[i]);

			if (params.rtTolerance < 0.0)
				print_help("--rt-tolerance <X>, the treshold should be X >= 0");
			continue;
        }
		// NP3 GOT minPeaksOutput PARAM
		if (! strcmp(argv[i],"--min-peaks-output"))
		{
			if (++i == argc)
			{
				print_help("missing min-peaks-output value!");
			}
			params.minNumPeaks = 1;
			params.minNumPeaks = atoi(argv[i]);

			if (params.minNumPeaks <= 0)
				print_help("--min-peaks-output <X>, the treshold should be X > 0");
			continue;
		}
		// NP3 scale factor parameter
        if (! strcmp(argv[i],"--scale-factor"))
        {
            if (++i == argc)
                print_help("missing scale factor number!");

            params.scaleFactor = MIN_FLOAT;
            params.scaleFactor = atof(argv[i]);
            if (params.scaleFactor<0.0)
                print_help("--scale-factor <X>, the number should be X >= 0.0");
            continue;
        }

		if (! strcmp(argv[i],"--fragment-tolerance"))
		{
			if (++i == argc)
				print_help("missing fragment tolerance threshold!");
			
			params.fragmentTolrance = MIN_FLOAT;
			params.fragmentTolrance = atof(argv[i]);
			if (params.fragmentTolrance<=0.0 || params.fragmentTolrance>=1.0)
				print_help("--fragment-tolerance <X>, the treshold should be in the range 0<X<1");
			continue;
		}

		if (! strcmp(argv[i],"--precursor-ppm"))
		{
			if (++i == argc)
				print_help("missing precursor ppm threshold!");
			
			params.precursorPPMs = MIN_FLOAT;
			params.precursorPPMs = atof(argv[i]);
			if (params.precursorPPMs<=0.0)
				print_help("--precursor-ppm <X>, the treshold should be in the range X>0");
			continue;
		}
		

		if (! strcmp(argv[i],"--peak-density"))
		{
			if (++i == argc)
				print_help("missing peak density!");
			params.peakDensity = atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--window"))
		{
			if (++i == argc)
				print_help("missing window size!");
			
			params.mzWindow = MIN_FLOAT;
			params.mzWindow = atof(argv[i]);
			if (params.mzWindow<=0.0)
				print_help("--window <X>, select X>0");
			continue;
		}

		if (! strcmp(argv[i],"--memory-gb"))
		{
			if (++i == argc)
				print_help("missing memory size!");
			
			params.memoryGb = MIN_FLOAT;
			params.memoryGb = atof(argv[i]);
			if (sizeof(size_t) == 4 && params.memoryGb>3.21)
				print_help("--memory-gb <X>, select 0.35<=X<=3.2 for 32bit or 0.35<=X for 64bit");
			continue;
		}

		if (! strcmp(argv[i],"--similarity"))
		{
			if (++i == argc)
				print_help("missing similarity threshold!");
			
			params.minSimilarity = MIN_FLOAT;
			params.minSimilarity = atof(argv[i]);
			if (params.minSimilarity>=1.0)
				print_help("--similarity <X>, select X<1.0");
			continue;
		}

		if (! strcmp(argv[i],"--mixture-prob"))
		{
			if (++i == argc)
				print_help("missing mixture prob threshold!");
			
			params.maxMixtureProb = MIN_FLOAT;
			params.maxMixtureProb = atof(argv[i]);
			if (params.maxMixtureProb >= 0.5)
				print_help("--mixture-prob <X>, select X<0.5");
			continue;
		}

		if (! strcmp(argv[i],"--num-rounds"))
		{
			if (++i == argc)
				print_help("missing snum rounds!");
			
			params.numRounds = 0;
			params.numRounds = atoi(argv[i]);
			if (params.numRounds<1 || params.numRounds>20)
				print_help("--num-rounds <X>, choose 1<=X<=20");
			continue;
		}

		if (! strcmp(argv[i],"--start-file-index"))
		{
			if (++i == argc)
				print_help("missing start file index!");
			
			params.startFileIdx = 0;
			params.startFileIdx = atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--verbose-level"))
		{
			if (++i == argc)
				print_help("missing verbose level!");
			
			params.verboseLevel = atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--assign-charges"))
		{
			params.gotAssignCharges = true;
			continue;
		}

		if (! strcmp(argv[i],"--dataset-idx"))
		{
			if (++i == argc)
				print_help("missing generation!");
			
			params.datasetIdx = 0;
			params.datasetIdx = atoi(argv[i]);
			gotDatasetIdx = true;
			continue;
		}

		if (! strcmp(argv[i],"--generation"))
		{
			if (++i == argc)
				print_help("missing generation!");
			
			params.datasetIdx = 0;
			params.datasetIdx = atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--batch-idx"))
		{
			if (++i == argc)
				print_help("missing batch!");
			
			params.batchIdx = 0;
			params.batchIdx = atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--correct-pm"))
		{
			params.gotCorrectPM = true;
			continue;
		}

		if (! strcmp(argv[i],"--output-mgf"))
		{
			params.gotOutputMgf = true;
			continue;
		}

		if (! strcmp(argv[i],"--output-dat"))
		{
			params.gotOutputDat = true;
			continue;
		}

		if (! strcmp(argv[i],"--output-file-size"))
		{
			if (++i == argc)
				print_help("missing output file size!");
			params.outputFileSize = atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--use-input-titles"))
		{
			params.gotUseInputTitles = true;
			continue;
		}

		if (! strcmp(argv[i],"--overwrite-locations"))
		{
			params.gotOverwriteLocations = true;
			continue;
		}

		if (! strcmp(argv[i],"--dat-only"))
		{
			params.gotDatOnly = true;
			continue;
		}

		if (! strcmp(argv[i],"--keep-dat"))
		{
			params.gotKeepDat = true;
			continue;
		}

		if (! strcmp(argv[i],"--first-pass"))
		{
			params.gotFirstPass = true;
			continue;
		}

		if (! strcmp(argv[i],"--second-pass"))
		{
			params.gotSecondPass = true;
			continue;
		}

		if (! strcmp(argv[i],"--major-increment"))
		{
			if (++i == argc)
				print_help("missing major increment!");
			
			int increment = 0;
			increment = atoi(argv[i]);
			if (increment<1)
				print_help("--major-increment <X>, choose X>=1");

			params.majorIncrement = static_cast<mass_t>(increment);
			continue;
		}

		if (! strcmp(argv[i],"--min-mz"))
		{
			if (++i == argc)
				print_help("missing min mz");
			float mz= 0.0;
			mz=atof(argv[i]);
			params.minMz = mz;
			continue;
		}

		if (! strcmp(argv[i],"--max-mz"))
		{
			if (++i == argc)
				print_help("missing max mz");
			float mz= 0.0;
			mz=atof(argv[i]);
			params.maxMz = mz;
			continue;
		}

		if (!  strcmp(argv[i],"--create-archive"))
		{
			params.gotCreateArchive = true;
			continue;
		}

		if (! strcmp(argv[i],"--merge-archives"))
		{
			params.gotMergeArchives = true;
			params.gotCreateArchive = true;
			if (++i == argc)
				print_help("missing path to archive 1");
			params.pathToArchive1 = argv[i];
			if (++i == argc)
				print_help("missing path to archive 2");
			params.pathToArchive2 = argv[i];
			continue;
		}


		if (! strcmp(argv[i],"--archive-statistics"))
		{
			params.gotArchiveStatistics=true;
			continue;
		}


		if (! strcmp(argv[i],"--output-modified"))
		{
			params.gotOutputModified = true;
			continue;
		}

		if (! strcmp(argv[i],"--search-archive"))
		{
			if (++i == argc)
				print_help("missing archive file!");
			params.archivePath=std::string(argv[i]);
			params.gotSearchArchive = true;
			continue;
		}

		if (! strcmp(argv[i],"--output-matched"))
		{
			params.gotOutputMatchedSpectra = true;
			continue;
		}

		if (! strcmp(argv[i],"--id-list"))
		{
			if (++i == argc)
				print_help("missing id list file");
			params.idListPath = std::string(argv[i]);
			continue;
		}
		
		if (! strcmp(argv[i],"--max-pvalue"))
		{
			if (++i == argc)
				print_help("missing p-value!");
			params.maxPvalue = atof(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--max-results"))
		{
			if (++i == argc)
				print_help("missing max results!");
			params.maxResultsForSpectrum = atoi(argv[i]);
			if (params.maxResultsForSpectrum<1)
				error("--max-results must be > 0 !");
			continue;
		}


		if (! strcmp(argv[i],"--ann-file"))
		{
			if (++i == argc)
				print_help("missing ann file!");
			params.annFile = argv[i];
			params.gotMakeMgf = true;
			continue;
		}

		if (! strcmp(argv[i],"--filter-only"))
		{
			params.gotFilterOnly=true;
			continue;
		}

		if (! strcmp(argv[i],"--make-crap-db"))
		{
			if (++i == argc)
				print_help("missing number of spectra!");
			params.numCrapSpectra=atoi(argv[i]);
			params.gotMakeCrapDb=true;
			continue;
		}

		if (! strcmp(argv[i],"--select-charge"))
		{
			if (++i == argc)
				print_help("missing number of spectra!");
			selectCharge=atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--benchmark-similarity"))
		{
			if (++i == argc)
				print_help("missing cluster list!");
			params.clusterList = argv[i];
			if (++i == argc)
				print_help("missing non-cluster list!");
			params.nonClusterList = argv[i];

			params.gotBenchmarkSimilarity = true;
			continue;
		}

		if (! strcmp(argv[i],"--sim-type"))
		{
			if (++i == argc)
				print_help("missing sim-type!");
			params.simType = atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--benchmark-consensus"))
		{
			if (++i == argc)
				print_help("missing cluster list!");
			params.clusterList = argv[i];

			params.gotBenchmarkConsensus = true;
			continue;
		}

		if (! strcmp(argv[i],"--benchmark-sqs"))
		{
			params.gotBenchmarkSqs = true;
			continue;
		}

		if (! strcmp(argv[i],"--benchmark-library"))
		{
			if (++i == argc)
				print_help("missing list of mgf files to load!");
			params.spectraListToLoad = argv[i];
			if (++i == argc)
				print_help("missing list of mgf files to test!");
			params.spectraListToTest = argv[i];

			params.gotBenchmarkLibrary = true;
			continue;
		}

		if (! strcmp(argv[i],"--benchmark-pairs"))
		{
			if (++i == argc)
				print_help("missing list of mgf files to load!");
			params.spectraListToLoad = argv[i];
			
			params.gotBenchmarkPairs = true;
			continue;
		}

		if (! strcmp(argv[i],"--benchmark-pm"))
		{
			params.gotBenchmarkPm = true;
			continue;
		}

		if (! strcmp(argv[i],"--benchmark-similarity-histogram"))
		{
			if (++i == argc)
				print_help("missing list of mgf files to load!");
			params.spectraListToLoad = argv[i];
			params.gotBenchmarkSimilarityHistogram = true;
			continue;
		}

		if (! strcmp(argv[i],"--intersection-fasta"))
		{
			if (++i == argc)
				print_help("missing path to fasta file!");
			params.intersectionFastaPath = argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--min-intersection-length"))
		{
			if (++i == argc)
				print_help("missing minimal intersection length!");
			params.minIntersectionLength=atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--min-number-for-intersection"))
		{
			if (++i == argc)
				print_help("missing minimal number of fastas in intersection!");
			params.minNumberForFastaIntersection=atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--train-similarity-model"))
		{
			params.gotTrainSimilarityModel = true;
			continue;
		}

		if (! strcmp(argv[i],"--correct-dat"))
		{
			gotCorrectDat = true;
			continue;
		}


	
		if (! strcmp(argv[i],"--help") || ! strcmp(argv[i],"-h"))
			print_help(NULL);
	

		printf("**********************************************************\n");
		printf("\nError: Unkown command line option: %s\n\n",argv[i]);
		printf("**********************************************************\n");
		exit(1); 
	}

	// first check for fasta file intersection, this doesn't invlovle most of the other code
	if (params.intersectionFastaPath.length()>0)
	{
		if (params.list.length() == 0)
			print_help("Must supply list of fasta files (using --list)");
		vector<string> fastas;
		readListOfPaths(params.list.c_str(), fastas);
		createIntersectionFasta(fastas, params.minIntersectionLength, params.intersectionFastaPath.c_str(), params.minNumberForFastaIntersection);
		exit(0);
	}


	// Always need a model because it creates the config (that is always needed)
	if (params.modelName.length() == 0)
		print_help("Must supply model name (--model)!");

	if (! params.gotCreateArchive && ! params.gotOutputDat)
		params.gotOutputMgf = true;

	if (params.gotCreateArchive && (gotDatasetIdx || indKeepDatasetIdx))
		error("When creating or merging archives you cannot use the flags \"--dataset-idx\" or \"--keep-dataset-idx\", dataset idx will be set according to defaults");

	AllScoreModels model;
	Config* config = model.get_config();
	
	if (params.modelDir.length()>0)
		config->set_resource_dir(params.modelDir);

	model.read_model(params.modelName.c_str());
	if (params.ptmString.length()>0)
		config->apply_selected_PTMs(params.ptmString.c_str());

	if (indKeepDatasetIdx)
		config->setKeepOriginalDatasetIdx(true);

	if (params.peakDensity>0)
	{
		if (params.peakDensity<15)
			error("--peak-density must have a value >= 15 !\n");
		config->set_max_number_peaks_per_local_window(params.peakDensity);
	}

	if (params.gotTrainSimilarityModel)
	{
		mass_t tolerance = (params.fragmentTolrance > 0.0 ? params.fragmentTolrance : model.get_config()->getTolerance());
		Cluster::setTolerances(tolerance, params.precursorPPMs);
		MsSimilarityModel msm;
	//	msm.readSimilarityModel(model.get_config());
		msm.trainSimilarityModel(model.get_config(), params.list.c_str());
	//	msm.makeSingleDistributionTables(model.get_config(), params.list.c_str());
		return 0;
	}

	if (gotCorrectDat)
	{
		correctDatError(params.list.c_str(), model.get_config());
		return 0;
	}

	if (params.gotConvert)
	{
		convertSpectra(&params, &model);
		return 0;
	}

	if (params.gotCnvertArchive)
	{
		convertArchive(&params, &model);
		return 0;
	}

	if (params.gotCreateArchiveFromMgfs)
	{
		MsArchive arch;
		arch.createArchiveFromMgf(&params, &model);
		return 0;
	}

	if (params.gotMakeExclusion)
	{
		makeScanListFromClust(params.list.c_str(), params.exclusionList.c_str());
		return 0;
	}

	if (params.gotMakeMgf)
	{
		AnnotationFile ann;
		ann.readAnnotationFile(params.annFile.c_str());
		ann.createAnnotatedMgfs(params.list.c_str(), model.get_config(), params.outputName);
		return 0;
	}

	if (params.checkIfGotBenchmark() )
	{
		performBenchmark(&model, &params);
		cout << "Done with benchmarks." << endl;
		return 0;
	}
 
	if (params.gotFilterOnly)
	{
		createDirIfDoesNotExist(params.outDir.c_str());

		if (params.modelName.length() <=1 || params.list.length()<1)
			error("Must supply: --model, --out-dir, --list, and --sqs-threshold (--output-name optional)");

		string suffix = "fil";
		if (params.outputName.length()>0)
			suffix = params.outputName;
		filterDataSet(&model, params.list, params.outDir, params.sqsThreshold, suffix);
		return 0;
	}

	if (selectCharge>0)
	{
		selectSpectra(&model, params.list, params.outDir, selectCharge);
		return 0;
	}

	if (params.gotMakeCrapDb)
	{
		const string name = (params.outputName.length()>1 ? params.outputName : std::string("crap"));
		const float  sqsThresh = (params.sqsThreshold>0.0 ? params.sqsThreshold : 0.05);
		if (params.numCrapSpectra<=0)
			error("Must supply positive size for crap db!");
		if (params.list.length()==0)
			error("Must supply one of --list!");

		createDirIfDoesNotExist(params.outDir.c_str());
		makeCrapDb(&model, params.list, params.outputName, params.outDir, sqsThresh, params.numCrapSpectra);
		return 0;
	}

	if (params.gotArchiveStatistics)
	{
		MsArchive arch;
		arch.readArchive(params.archivePath);
		arch.printStatistics(config);
		return 0;
	}

	if (params.gotSearchArchive)
	{
		MsArchive arch;
		if (! arch.readArchive(params.archivePath))
			error("Could not read archive: ", params.archivePath.c_str());
		arch.searchArchive(&params, &model);
		return 0;
	}


	// This is a clustering job...

	// NP3 GOT set RT tol to config
	config->set_rt_tolerance(params.rtTolerance);
	config->set_scale_factor(params.scaleFactor);

	if (params.outputName.length() == 0)
		print_help("Must supply job output name (--output-name)!");

	params.outputStub = params.outDir + "/" + params.outputName;

	if (! params.gotMergeArchives)
	{
		if (params.datList.length()==0 && params.list.length()==0 && params.metaList.length()==0)
			print_help("Must supply one of --list, --dat-list, or --meta-list !");

		if (params.gotCreateArchive && params.list.length()==0)
			error("When creating an archive you must supply the original list of spectra files with \"--list\" (you can use it along with \"--dat-list\").\n");
	}

	// when working with archives the increment should always be 1
	if (params.gotCreateArchive || params.archivePath.length()>0 || params.gotOutputDat)
		params.majorIncrement = 1.0; 

	// set datasetIdx to 0 if it wasn't set by a previous archive
	if (params.datasetIdx == MAX_INT)
		params.datasetIdx = 0;

	if (params.verboseLevel>0)
	{
		cout << "MsCluster v2.0 " << endl << "Started job on: " << getTimeString() << endl << endl;
		params.printParameters();
	}

	/*! \todo There should be no second pass. Instead, we should look at the number of peaks/cluster
	          spots are available, and decide if a slice of dat files needs to be split or not (if 
			  they all fit in the memory (with a margin), then there is no reason to split the file...
	*/
	// make DAT first
	bool madeFirstRoundDat=false;
	if ((params.datList.length()==0 && 
		(params.list.length()>0 || params.metaList.length()>0))
		|| params.gotSecondPass)
	{
		// make sure dat dir is present
		createDirIfDoesNotExist(params.tmpDir.c_str(), params.verboseLevel);
		DatFileWriter dfw(&model);
		const double startTime = time(NULL);

		// perform first pass
		if (! params.gotSecondPass)
		{
			params.datList = dfw.convertDataToDatFirstPass(&params);

			if (params.gotFirstPass)
			{
				cout << endl << "Finished with first pass of dat creation." << endl;
				cout << "Elapsed time: " << time(NULL) - startTime << endl;
				cout << "Exiting..." << endl;
				exit(0);
			}
			madeFirstRoundDat=true;
		}

		// perform second pass
		params.datList = dfw.convertDataToDatSecondPass(params.datList, params.tmpDir, params.outputName, 
			params.sqsThreshold, params.verboseLevel);
		if (params.gotSecondPass)
			cout << endl << "Finished with second pass of dat creation." << endl;
			
		if (params.gotDatOnly)
		{
			cout << endl << "Finished creating dat files." << endl;
			cout << "Elapsed time: " << time(NULL) - startTime << endl;
			cout << "Exiting..." << endl;
			exit(0);
		}
		cout << "Time required for dat file creation: " << time(NULL) - startTime << endl;
	}

	// make sure out dir is present
	createDirIfDoesNotExist(params.outDir.c_str(), params.verboseLevel);

	// Perform clsutering
	MsClusterAlgorithm alg;
	if (params.gotMergeArchives)
	{
		alg.mergeTwoArchives(&params, &model);
	}
	else
		alg.performClustering(&params, &model);

	if (madeFirstRoundDat && ! params.gotKeepDat)
		removeFilesInList(params.datList.c_str());

	if (params.verboseLevel>0)
		cout << endl << "Ended job on: " << getTimeString() << endl;
	return 0;
}








