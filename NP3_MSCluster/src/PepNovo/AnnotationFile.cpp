#include "AnnotationFile.h"
#include "SpectraAggregator.h"
#include "SpectraList.h"
#include "PeakList.h"
#include "PepNovo_auxfun.h"
#include "../Common/BufferedOutputFile.h"

void AnnotationFile::readAnnotationFile(const char* path, bool indClear)
{
	FILE* stream = fopen(path,"r");
	if (! stream)
		error("Could not open annotation file for reading: ",path);

	if (indClear)
		annotations_.clear();

	int numAnns=0;
	char line[256];
	char pepBuffer[256];
	while (fgets(line, 256, stream) )
	{
		int fileIdx, scan, charge;
		if (sscanf(line,"%d %d %s %d",&fileIdx,&scan,pepBuffer,&charge) == 4)
		{
			FileAndScan fas;
			Annotation ann;
			fas.fileIdx = fileIdx;
			fas.scan    = scan;
			ann.peptideStr = pepBuffer;
			ann.charge     = charge;

			annotations_[fas]=ann;
			numAnns++;
		}
	}
	fclose(stream);
	cout << "Read " << numAnns << " annotations (" << annotations_.size() << ")" << endl;
}

void AnnotationFile::createAnnotatedMgfs(const char* fileList, 
										 const Config* config,
										 string& name,
										 bool onlyAnnotated,
										 size_t maxNumPerFile) const
{
	SpectraAggregator sa;
	sa.initializeFromTextFile(fileList, config);

	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();

	size_t fileIdx  = 0;
	size_t numScans = 0;
	size_t numMismatches = 0;

	char* mgfBuffer = new char[131072];
	BufferedOutputFile bof;
	
	cout << "found " << sl.getNumHeaders() << " spectrum headers." << endl;
	if (sl.getNumHeaders() == 0)
		error("No spectra!");

	for (size_t i=0; i<sl.getNumHeaders(); i++)
	{
		SingleSpectrumHeader* header = const_cast<SingleSpectrumHeader*>(sl.getSpectrumHeader(i));

		FileAndScan fas;
		fas.fileIdx = header->getSpectraFileIndexInList();
		fas.scan    = (header->getFileType() == IFT_MGF ? header->getIndexInFile() : header->getScanNumber());

		map<FileAndScan, Annotation>::const_iterator it = annotations_.find(fas);
		
		if (it == annotations_.end())
			continue;

		header->setCharge(it->second.charge);
		header->setPeptideStr(it->second.peptideStr);

		Peptide pep;
		pep.parseFromString(config, it->second.peptideStr);
		pep.get_mass_with_19();
		mass_t theoreticalMz = (pep.get_mass_with_19() + (it->second.charge-1)*MASS_PROTON)/it->second.charge;

		if (fabs(theoreticalMz - header->getMOverZ())>10.0)
		{
			cout << "Type: " << (header->getFileType() == IFT_MGF ? "MGF" : "OTHER") << endl;
			cout << "Mass mismatch for ann " << fas.fileIdx << "\t" << fas.scan << endl;
			cout << "Spectrum m/z: " << header->getMOverZ() << " vs. " << it->second.peptideStr << " charge "
				 << it->second.charge << " which gives " << theoreticalMz << endl;
			numMismatches++;
			continue;
		}
		

		PeakList pl;
		pl.readPeaksToLocalAllocation(sa, header);
		pl.initializePeakList(config, true);


		if (! bof.isOpen())
		{
			ostringstream oss;
			oss << name << "_" << fileIdx++ << ".mgf";
			bof.open(oss.str(), 1<<20);
			cout << "Opened " << oss.str() << endl;
		}

		size_t n =pl.outputToMgfBuffer(mgfBuffer);
		bof.writeToBuffer(mgfBuffer, n);
	
		numScans++;
		if (numScans > 0 && numScans % maxNumPerFile == 0)
		{
			
			bof.close();
		}
	}
	bof.close();

	cout << "Wrote " << numScans << " annotated spectra to " << fileIdx << " files." << endl;
	cout << "Found " << numMismatches << " annotation mismatches. " << endl;
}



