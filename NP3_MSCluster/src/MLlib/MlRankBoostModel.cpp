#include "MlRankBoost.h"

double MlRankBoostModel::calcErrorRate(const MlDataSet* mld, bool verbose) const
{
	return 0.0;
}


double MlRankBoostModel::calcErrorRate(const MlDataSet* mld, const vector<size_t>& idxs, bool verbose) const
{
	return 0.0;
}
	
void  MlRankBoostModel::trainModel(MlTrainingContainer* mltd)
{

}


bool MlRankBoostModel::readModel(const char* path)
{
	return true;
}

bool MlRankBoostModel::readModel(ifstream& ifs)
{
	return true;
}

bool MlRankBoostModel::writeModel(ostream& os) const
{
	return true;
}

void MlRankBoostModel::outputModelAnalysis(const MlDataSet* mld, 
										   const vector<size_t>& idxs, 
										   const MlFeatureSet* featureSet, 
										   ostream& os) const
{

}




