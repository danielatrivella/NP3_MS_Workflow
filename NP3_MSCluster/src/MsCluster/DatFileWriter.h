#ifndef __DATFILEWRITER_H__
#define __DATFILEWRITER_H__

/*! @file DatFileWriter.h
Holds the class DatFileWriter
*/

#include "MsClusterIncludes.h"
#include "../PepNovo/DatFile.h"
#include "../PepNovo/AllScoreModels.h"

struct MsParameterStruct; // fwd dclr

/*! @class DatFileWriter
Manages the process of writing spectra to dat files.

This class's functions are designed for specific tasks such as writing dat files for a lage set of 
spectra files or quick splitting of dat files with large m/z increments to dat files with smaller
increments. If the dat files are to be used simply as a container for spectra (i.e., conceptually 
the same as an mzXML or mgf), then the class DatFile should probably be used.
*/
class DatFileWriter {
public:

	DatFileWriter(AllScoreModels* model) :   model_(model), config_(0), mzIncrement_(20.0), 
		maxMz_(5000.0), maxDatFileSize_(1<<29), maxBinIdx_(0), verboseLevel_(0),
		numSpectraWrittenFirstPass_(0), numSpectraWrittenSecondPass_(0), 
		numSpectraReadFromOriginalFiles_(0), numOriginalPaths_(0) 
		{ init(MAJOR_MZ_INCREMENT_FOR_DAT); }

	//1073741824
	~DatFileWriter();

	/*! Reads spectra files (mgf, mzXML, dat) and writes them into dat files with large (e.g., 25 Da)
	increments. Also performs quality filtering at this stage if requested.


	@param params the parameters struct; should include:
	-datDir Path to the directory where the dat files are to be written.
	-datName The prefix to all file/cluster names that will be written.
	-sqsThreshold Minimal quality score for a spectrum to be kept; values in the range(0,1). If 
	no filtering is to be done, the value 0 should be passed.
	-fileStartIdx The offset that should be added to the file indexes that get assigned to the scans.
	This parameters should only be used when the dat creation is done in parallel (when splitting 
	the data creation of large clustering jobs).
	*/
	string convertDataToDatFirstPass(const MsParameterStruct* params);

	string convertDataToDatSecondPass(const string& datList,
									  const string& datDir, 
									  const string& datName, 
									  float sqsThreshold = 0.0,
									  int verboseLevel=0);
									  
	
	void init(mass_t mzIncrement);

	void closeAllOpenDats();

	string writeDatPaths(bool indSort = true);
	
	void addPeakListToDat(const PeakList& pl, bool ignoreMzIdx = false);

	mass_t getMzIncrement() const   { return mzIncrement_; }
	void   setMzIncrement(mass_t m) { mzIncrement_ = m; }

	mass_t getMaxMz() const		  { return maxMz_; }
	void   setMaxMz(mass_t m)     { maxMz_ = m; }

	size_t getMaxDatFileSize() const   { return maxDatFileSize_; }
	void   setMaxDatFileSize(size_t s) { maxDatFileSize_ = s; }

	const string& getDatDir() const   { return datDir_; }
	void  setDatDir(const string& s)  { datDir_ = s; }
	
	const string& getDatName() const  { return datName_; }
	void  setDatName(const string& s) { datName_ = s; }

	clusterIdx_t getNumSpectraWrittenSecondPass() const { return numSpectraWrittenSecondPass_; }

	void generateDatPaths(vector<string>& datPaths) const
	{
		datPaths.clear();
		for (size_t i=0; i<datFiles_.size(); i++)
			datPaths.push_back(datFiles_[i].getPath());
	}

protected:
	
	AllScoreModels* model_;
	const Config*   config_;

	mass_t mzIncrement_;
	mass_t maxMz_;
	size_t maxDatFileSize_;
	size_t maxBinIdx_;

	int	   verboseLevel_;
	clusterIdx_t numSpectraWrittenFirstPass_;
	clusterIdx_t numSpectraWrittenSecondPass_;
	clusterIdx_t numSpectraReadFromOriginalFiles_;
	int  numOriginalPaths_;
	
	vector<DatFile> datFiles_;
	vector<size_t>  fileCounters_; // running indexes on the files
	vector<size_t>  fileSizes_;	   // counts size of files

	vector<string>  datPaths_;     // the paths of the created DAT files
	string datDir_, datName_;
	
	// writes the spectrum directly from a buffer
	void writeDatSpectrumDirectly(const char* spec, size_t numBytes);

	// reads directly from the dat file and writes it in the DatFileWriter
	// without performing any additional processing except sqs filtering
	// this is used in the second pass
	size_t writeDatFileToOtherDats(const string& datPath, 
								   const string& datDir, 
								   const string& datName,
								   float sqsThreshold,
								   size_t& totalSpectraRead);

};





#endif

