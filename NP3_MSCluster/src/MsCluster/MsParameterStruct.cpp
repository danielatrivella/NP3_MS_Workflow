#include "MsParameterStruct.h"

void MsParameterStruct::copyCommandLine(int argc, char** argv)
{
	commandLine = std::string();
	if (argc <= 0)
		return;

	commandLine = std::string(argv[0]);
	for (int i=1; i<argc; i++)
		commandLine += (" " + std::string(argv[i]));
}

void MsParameterStruct::printParameters(ostream& os) const
{
	os << "Parameter values:" << endl;
	os << "-----------------" << endl;
	
	if (outputName.length()>0)
		os << "outputName = " << outputName <<endl;

	if (modelName.length()>0 && modelName != "LTQ_TRYP")
		os << "modelName  = " << modelName << endl;

	if (modelDir.length()>0 && modelDir != "Models")
		os << "modelDir   = " << modelDir << endl;

	if (outDir.length()>0 && outDir != "out")
		os << "outDir     = " << outDir << endl;

	if (tmpDir.length()>0 && tmpDir != "tmp")
		os << "tmpDir     = " << tmpDir << endl;

	if (list.length()>0)
		os << "list	      = " << list << endl;

	if (datList.length()>0)
		os << "datList    = " << datList << endl;

	if (archivePath.length()>0)
		os << "Archive    = " << archivePath << endl;

	if (ptmString.length()>0)
		os << "ptmString  = " << ptmString << endl;

	if (sqsThreshold>0.0)
		os << "sqsThreshold  = " << setprecision(4) << sqsThreshold << endl;

	os <<	  "mixtureProb   = " << setprecision(4) << maxMixtureProb << endl;
	os <<     "minSimilarity = " << setprecision(4) << minSimilarity << endl;
	// NP3 additional parms
	os <<     "rtTolerance = " << setprecision(2) << rtTolerance << endl;
	os <<     "minNumPeaks = " << minNumPeaks << endl;
    os <<     "scaleFactor = " << scaleFactor << endl;



	os << endl << "COMMAND LINE = " << commandLine << endl;
}



