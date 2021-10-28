#include "mlfeature.h"
#include "mlauxfun.h"


// adds a suffix in the form _{XX} to the feature name
void MlFeature::addSuffixToFeatureName(const char* suffix)
{
	if (name_[name_.length()-1] != '}')
	{
		name_ += "_{";
		name_ += static_cast<string>(suffix);
		name_ += "}";
		return;
	}

	// remove the previous suffix
	size_t i=name_.length()-1;
	while (i>0 && name_[i] != '{')
		i--;

	if (i==0)
	{
		name_ += "_{";
		name_ += static_cast<string>(suffix);
		name_ += "}";
	}
	else
	{
		name_ = name_.substr(0,i);
		name_ += static_cast<string>(suffix);
		name_ += "}";
	}
}


void MlFeature::addGroups(const vector<size_t>& newIdxs)
{
	for (size_t i=0; i<newIdxs.size(); i++)
	{
		size_t j;
		for (j=0; j<groups_.size(); j++)
			if (newIdxs[i]==groups_[j])
				break;
		if (j==groups_.size())
			groups_.push_back(newIdxs[i]);
	}
}


void MlFeatureSet::initializeWithDefaultFeatures(size_t numFeatures)
{
	basicFeatures_.clear();
	basicFeatures_.resize(numFeatures);

	for (size_t i=0; i<numFeatures; i++)
	{
		char buff[16];
		sprintf(buff,"F%d",static_cast<int>(i));
		basicFeatures_[i].name_ = static_cast<string>(buff);
	}

	indWasInitialized_ = true;
}


void MlFeatureSet::readFeatureFile(ifstream& ifs)
{
	basicFeatures_.clear();
	char buffer[256];
	while (! ifs.eof())
	{
		ifs.getline(buffer,256);
		if (ifs.gcount()<=3)
			continue;

		if (buffer[0] == '#')
		{
			int n=0;
			if (sscanf(buffer,"#NUM FEATURES %d",&n) == 1)
				basicFeatures_.resize(n);

			continue;
		}

		if (basicFeatures_.size() == 0)
			error("Must supply number of features \"#NUM FEAUTRES xxx\"");

		istringstream iss(buffer);
		size_t index;
		iss >> index;
		if (iss.fail())
			continue;

		if (index>=basicFeatures_.size())
			error("Feature index out of bounds ",index);

		MlFeature& mlf = basicFeatures_[index];

		iss >> mlf.name_;

		if (mlf.name_.length()<1)
			error("bad feature line: ",buffer);

		// read groups
		while (true)
		{
			string g="";
			iss >> g;
			if (g.length()>0)
			{
				// get group index
				map<string,size_t>::const_iterator it=this->groupMapping_.find(g);
				if (it != groupMapping_.end())
				{
					mlf.groups_.push_back(it->second);
					continue;
				}

				groupMapping_[g]=groupNames_.size();
				groupNames_.push_back(g);
				continue;
			}
			break;
		}

		if (index == basicFeatures_.size()-1)
			break; // quit reading file
	}
	
	

	// create lists of group members
	groupMembers_.resize(groupNames_.size());
	for (size_t i=0; i<basicFeatures_.size(); i++)
	{
		size_t j;
		for (j=0; j<basicFeatures_[i].groups_.size(); j++)
			groupMembers_[basicFeatures_[i].groups_[j]].push_back(i);
	}

	indWasInitialized_ = true;
}

void MlFeatureSet::readFeatureFile(const char* path, bool verbose)
{
	basicFeatures_.clear();
	ifstream ifs;
	ifs.open(path,ios::in);

	if (! ifs.good())
		error("couldn't open feature file for reading: ",path);

	path_ = std::string(path);
	readFeatureFile(ifs);
	ifs.close();

	if (verbose)
		cout << "Read " << basicFeatures_.size() << " features from " << path << endl;
}


void MlFeatureSet::writeFeatureFile(ostream& os) const
{
	if (! indWasInitialized_)
		error("trying to write uninitialized feature set!");

	os << "#NUM FEATURES " << basicFeatures_.size() << endl;

	for (size_t i=0; i<basicFeatures_.size(); i++)
	{
		char buffer[256];
		size_t len=sprintf(buffer,"%d\t%s",static_cast<int>(i),basicFeatures_[i].getName().c_str());

		const vector<size_t>& groups = basicFeatures_[i].getGroups();
		char* pos = buffer + len;
		size_t j;
		for (j=0; j<groups.size(); j++)
			pos+=sprintf(pos,"\t%s",this->groupNames_[groups[j]].c_str());

		os << buffer << endl;
	}
}




void MlFeatureSet::writeFeatureFile(const char* path, bool verbose) const
{
	ofstream ofs(path);
	if (! ofs.is_open())
		error("couldn't open feature file for writing: ",path);

	ofs.close();
	if (verbose)
		cout << "Wrote " << basicFeatures_.size() << " features to " << path << endl;
}


