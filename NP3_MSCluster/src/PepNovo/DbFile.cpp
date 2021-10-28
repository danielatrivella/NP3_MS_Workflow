#include "DbFile.h"
#include "../Common/auxfun.h"


const char separator[]={"[:-(]"};


bool DbFile::readFastaFile(const std::string &fasta, size_t partIdx)
{
	// first pass: count number of aas, number of proteins and lengths of names
	// need to make sure the reading can handle long lines
	
	ifstream ifs(fasta.c_str());
	if (! ifs.good())
		error("Couldn't open fasta file for reading: ",fasta.c_str());

	fasta_ = fasta;
	partIdx_ = partIdx;
	
	char* nameBuffer= new char[20000]; // name should not be longer than this
	size_t totalAaLengths=0;
	size_t totalNameLengths=0;
	size_t numSequences=0;
	while (ifs.good())
	{
		ifs.getline(nameBuffer,20000);
		if (ifs.gcount() < 2 || nameBuffer[0] != '>')
			continue;

		totalNameLengths += ifs.gcount();
		numSequences++;

		char c;
		ifs.get(c);
		while (ifs.good() && c != '>')
		{
			if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c<= 'z'))
				totalAaLengths++;
			ifs.get(c);
		}
		ifs.unget();
		totalAaLengths++;
	}

	ifs.close();
	totalAaLengths++;

	// allocate storage area
	
	try
	{
		aaArray_.resize(totalAaLengths);
		nameArray_.resize(totalNameLengths);
		proteinStartPointers_.resize(numSequences);
		proteinNames_.resize(numSequences);
	}
	catch ( ... )
	{
		error("Could not alloacte memory for fasta file, try breaking file into smaller parts!");
	}

	const char* aaArrayBefore = &aaArray_[0];

	// second pass: write data from fasta
	ifs.clear();
	ifs.open(fasta.c_str(),ios::in);
	if (! ifs.good())
		error("could not read fasta file second time: ",fasta.c_str());

	size_t namePos=0;
	size_t aaPos=0;
	size_t idx=0;

	while (ifs.good())
	{
		ifs.getline(nameBuffer,20000);
		if (ifs.gcount() < 2 || nameBuffer[0] != '>')
			continue;

		size_t nameLength = ifs.gcount()-1;
		assert(nameLength>0 && nameLength<100000);

		proteinNames_[idx]=namePos;
		memcpy(&nameArray_[namePos], nameBuffer+1, nameLength);
		namePos+=nameLength;
		nameArray_[namePos-1]='\0';
		while (namePos>1 && (nameArray_[namePos-1] == '\n' || nameArray_[namePos-1] == ' ' || 
							 nameArray_[namePos-1] == '\t' || nameArray_[namePos-1] == '\b'))
		{
			namePos--;
			nameArray_[namePos-1]='\0';
		}

		proteinStartPointers_[idx]=aaPos;
		char c;
		ifs.get(c);
		while (ifs.good() && c != '>')
		{
			if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c<= 'z'))
				aaArray_[aaPos++]=c;
			ifs.get(c);
		}
		ifs.unget();
		aaArray_[aaPos++]='\0';
		assert(aaPos<totalAaLengths);

		idx++;
	}
	ifs.close();

	assert(proteinStartPointers_.size() == proteinNames_.size());

//	printStats();

	delete [] nameBuffer;

	return true;
}


bool DbFile::writeFastaFile(const string& path) const
{
	assert(proteinStartPointers_.size() == proteinNames_.size());
	ofstream ofs(path.c_str());

	for (size_t i=0; i<proteinStartPointers_.size(); i++)
	{
		ofs << ">" << proteinNames_[i] << endl;
		const char* ptr= &aaArray_[0] + proteinStartPointers_[i];
		size_t counter=0;
		while (*ptr != '\0')
		{
			ofs << *ptr++;
			if (++counter == 80)
			{
				counter=0;
				ofs<<endl;
			}
		}
		ofs << endl;
		if (counter>0)
			ofs << endl;
	}
	ofs.close();
	return true;
}

bool DbFile::createTags(unsigned int minTagLength, unsigned int maxTagLength)
{
	allTags_.clear();

	// first pass, find out what tags are there and get counts

	for (unsigned int tagLength=minTagLength; tagLength<=maxTagLength; tagLength++)
	{
		for (size_t i=0; i<proteinStartPointers_.size(); i++)
		{
			// check that protein is not too short
			const char* ptr= &aaArray_[0] + proteinStartPointers_[i];
			size_t j;
			for (j=0; j<tagLength; j++)
			{
				if (*ptr++ == '\0')
					break;
			}
			if (j<tagLength)
				continue;

			// parse all tags from the protein
			ptr = &aaArray_[0] + proteinStartPointers_[i];
			while (1)
			{
				TagWithLocations twl(ptr, tagLength);
				set<TagWithLocations>::iterator it = allTags_.find(twl);
				
				if (it != allTags_.end())
				{
					twl.numLocations = it->numLocations+1;
					(it->numLocations)++;
				}
				else
				{
					twl.numLocations = 1;
					allTags_.insert(twl);
				}

				if (*(ptr+tagLength) == '\0')
					break;
				ptr++;
			}
		}
	}

	// allocate memory for locations
	size_t totalNumLocations = 0;
	for (set<TagWithLocations>::iterator it = allTags_.begin(); it != allTags_.end(); it++)
	{
		(it->locationsOffset) = totalNumLocations; // this is the offset into the array
		totalNumLocations += (it->numLocations+1);
		(it->numLocations) = 0; // this will be used as a counter, so it should be reset
	}
	locationArray_.resize(totalNumLocations, 0);

	// second pass, write the locations of the tags
	size_t numLocationsWritten = 0;
	for (size_t tagLength=minTagLength; tagLength<=maxTagLength; tagLength++)
	{
		for (size_t i=0; i<proteinStartPointers_.size(); i++)
		{
			// check that protein is not too short
			const char* ptr = &aaArray_[0] + proteinStartPointers_[i];
			size_t j;
			for (j=0; j<tagLength; j++)
			{
				if (*ptr++ == '\0')
					break;
			}
			if (j<tagLength)
				continue;

			// parse all tags from the protein
			const char* aaArrayStart = &aaArray_[0];
			ptr = &aaArray_[0] + proteinStartPointers_[i];
			while (1)
			{
				TagWithLocations twl(ptr, tagLength);
				set<TagWithLocations>::iterator it = allTags_.find(twl);

				assert(it != allTags_.end());
				locationArray_[it->locationsOffset + it->numLocations] = (ptr - aaArrayStart);
				it->numLocations++;
				numLocationsWritten++;
				if (*(ptr+tagLength) == '\0')
					break;
				ptr++;
			}
		}
	}

	assert(numLocationsWritten + allTags_.size() == locationArray_.size());
//	printTags();

	return true;
}


bool DbFile::writeDbFile(const std::string &path) const
{
	if (aaArray_.size() >= MAX_UINT || nameArray_.size() >= MAX_UINT ||
		locationArray_.size() >= MAX_UINT)
	{
		error("DB size too large to write! make sure that number of AAs and locations can be expressed using 32 bits!");
	}

	ofstream ofs(path.c_str(),ofstream::binary);
	if (! ofs.good())
		error("Could not open db file for writing: ", path.c_str());

	const unsigned int fastaLength = fasta_.length();
	ofs.write(reinterpret_cast<const char*>(&fastaLength), sizeof(unsigned int));
	ofs.write(fasta_.c_str(), fastaLength);
	ofs.write(reinterpret_cast<const char*>(&partIdx_), sizeof(unsigned int));
	ofs.write(reinterpret_cast<const char*>(&minTagLength_), sizeof(unsigned int));
	ofs.write(reinterpret_cast<const char*>(&maxTagLength_), sizeof(unsigned int));

	
	ofs.write(separator, sizeof(separator));

	const unsigned int nameArrayLength = nameArray_.size();
	ofs.write(reinterpret_cast<const char*>(&nameArrayLength), sizeof(unsigned int));
	ofs.write(&nameArray_[0], nameArrayLength);

	ofs.write(separator, sizeof(separator));

	const unsigned int aaArrayLength = aaArray_.size();
	ofs.write(reinterpret_cast<const char*>(&aaArrayLength), sizeof(unsigned int));
	ofs.write(&aaArray_[0], aaArrayLength);

	ofs.write(separator, sizeof(separator));

	const unsigned int proteinStartPointersLength = proteinStartPointers_.size();
	ofs.write(reinterpret_cast<const char*>(&proteinStartPointersLength), sizeof(unsigned int));
	ofs.write(reinterpret_cast<const char*>(proteinStartPointers_[0]), sizeof(unsigned int) * proteinStartPointersLength);

	ofs.write(separator, sizeof(separator));

	const unsigned int proteinNamesLength = proteinNames_.size();
	ofs.write(reinterpret_cast<const char*>(&proteinNamesLength), sizeof(unsigned int));
	ofs.write(reinterpret_cast<const char*>(proteinNames_[0]), sizeof(unsigned int) * proteinNamesLength);

	ofs.write(separator, sizeof(separator));

	const unsigned int locationArrayLength = locationArray_.size();
	ofs.write(reinterpret_cast<const char*>(&locationArrayLength), sizeof(unsigned int));
	ofs.write(reinterpret_cast<const char*>(locationArray_[0]), sizeof(unsigned int) * locationArrayLength);

	ofs.write(separator, sizeof(separator));

	// locations need to be written one at time because the the char* needs to be converted to an unsigned int
	const unsigned int allTagsSize = allTags_.size();
	ofs.write(reinterpret_cast<const char*>(&allTagsSize), sizeof(unsigned int));

	const char* aaArrayStart = &aaArray_[0];
	unsigned int buffer[4];
	const size_t bufferSize = 4*sizeof(unsigned int);
	
	for (set<TagWithLocations>::const_iterator it = allTags_.begin(); it != allTags_.end(); it++)
	{
		buffer[0] = static_cast<unsigned int>(it->tag - aaArrayStart);
		buffer[1] = it->tagLength;
		buffer[2] = it->numLocations;
		buffer[3] = it->locationsOffset;
		ofs.write(reinterpret_cast<const char*>(buffer), bufferSize);
	}

	ofs.write(separator, sizeof(separator));

	return true;
}

bool DbFile::readDbFile(const std::string &path)
{
	ifstream ifs(path.c_str(), istream::binary);
	if (! ifs.good())
		error("Could not open db file for reading: ",path.c_str());

	try
	{
		char separatorBuffer[sizeof(separator)];

		unsigned int fastaLength;
		ifs.read(reinterpret_cast<char*>(&fastaLength), sizeof(unsigned int));
		assert(fastaLength<5000);
		
		vector<char> buffer(fastaLength+1);
		ifs.read(&buffer[0], fastaLength);
		buffer[fastaLength]='\0';
		fasta_ = std::string(&buffer[0]);

		ifs.read(reinterpret_cast<char*>(&partIdx_), sizeof(unsigned int));
		ifs.read(reinterpret_cast<char*>(&minTagLength_), sizeof(unsigned int));
		ifs.read(reinterpret_cast<char*>(&maxTagLength_), sizeof(unsigned int));

		ifs.read(separatorBuffer, sizeof(separator));
		assert(! strncmp(separator, separatorBuffer, sizeof(separator)) );

		unsigned int nameArrayLength = 0;
		ifs.read(reinterpret_cast<char*>(&nameArrayLength), sizeof(unsigned int));
		nameArray_.resize(nameArrayLength);
		ifs.read(&nameArray_[0], nameArrayLength);

		ifs.read(separatorBuffer, sizeof(separator));
		assert(! strncmp(separator, separatorBuffer, sizeof(separator)) );

		unsigned int aaArrayLength = 0;
		ifs.read(reinterpret_cast<char*>(&aaArrayLength), sizeof(unsigned int));
		aaArray_.resize(aaArrayLength);
		ifs.read(&aaArray_[0], aaArrayLength);

		ifs.read(separatorBuffer, sizeof(separator));
		assert(! strncmp(separator, separatorBuffer, sizeof(separator)) );

		unsigned int proteinStartPointersLength = 0;
		ifs.read(reinterpret_cast<char*>(&proteinStartPointersLength), sizeof(unsigned int));
		proteinStartPointers_.resize(proteinStartPointersLength);
		ifs.read(reinterpret_cast<char*>(proteinStartPointers_[0]), sizeof(unsigned int) * proteinStartPointersLength);

		ifs.read(separatorBuffer, sizeof(separator));
		assert(! strncmp(separator, separatorBuffer, sizeof(separator)) );

		unsigned int proteinNamesLength = 0;
		ifs.read(reinterpret_cast<char*>(&proteinNamesLength), sizeof(unsigned int));
		proteinNames_.resize(proteinNamesLength);
		ifs.read(reinterpret_cast<char*>(proteinNames_[0]), sizeof(unsigned int) * proteinNamesLength);

		ifs.read(separatorBuffer, sizeof(separator));
		assert(! strncmp(separator, separatorBuffer, sizeof(separator)) );

		unsigned int locationArrayLength = 0;
		ifs.read(reinterpret_cast<char*>(&locationArrayLength), sizeof(unsigned int));
		locationArray_.resize(locationArrayLength);
		ifs.read(reinterpret_cast<char*>(locationArray_[0]), sizeof(unsigned int) * locationArrayLength);

		ifs.read(separatorBuffer, sizeof(separator));
		assert(! strncmp(separator, separatorBuffer, sizeof(separator)) );

		// locations need to be written one at time because the the char* needs to be converted to an unsigned int
		unsigned int allTagsSize = allTags_.size();
		ifs.read(reinterpret_cast<char*>(&allTagsSize), sizeof(unsigned int));
		allTags_.clear();

		const char* aaArrayStart = &aaArray_[0];
		unsigned int tmpBuffer[4];
		const size_t bufferSize = 4*sizeof(unsigned int);

		pair<set<TagWithLocations>::const_iterator,bool> ret(allTags_.begin(),true);
		for (size_t i=0; i<allTags_.size(); i++)
		{
			ifs.read(reinterpret_cast<char*>(tmpBuffer), bufferSize);
			TagWithLocations twl(tmpBuffer[0] + aaArrayStart, tmpBuffer[1]);
			twl.numLocations = tmpBuffer[2];
			twl.locationsOffset = tmpBuffer[3];
			ret = allTags_.insert(twl); // would like to put ret.first there, but it won't compile
			assert( ret.second );
		}

		ifs.read(separatorBuffer, sizeof(separator));
		assert(! strncmp(separator, separatorBuffer, sizeof(separator)) );

	}
	catch (...)
	{
		return false;
	}


	return true;
}



void DbFile::computeProteinIdxAndOffsetFromLocation(unsigned int location, size_t& proteinIdx, size_t& offset) const
{
	proteinIdx=0;
	offset=0;
	if (location==0)
		return;

	vector<unsigned int>::const_iterator it = lower_bound(proteinStartPointers_.begin(), proteinStartPointers_.end(), location);
	if (it == proteinStartPointers_.end())
	{
		assert( location >= proteinStartPointers_.back() && location < aaArray_.size());
		proteinIdx=proteinStartPointers_.size()-1;
		offset = location - proteinStartPointers_.back();
		return;
	}

	if (*it == location)
	{
		proteinIdx = it - proteinStartPointers_.begin();
		return;
	}

	proteinIdx = it - proteinStartPointers_.begin() - 1;
	assert(proteinIdx<proteinStartPointers_.size());
	assert(location > proteinStartPointers_[proteinIdx]);
	offset = location - proteinStartPointers_[proteinIdx];
}


void DbFile::printStats() const
{
	cout << "fasta: " << this->fasta_ << endl;
	cout << "Num proteins: " << this->proteinNames_.size() << endl;
	for (size_t i=0; i<proteinNames_.size(); i++)
		cout << i << "\t" << proteinStartPointers_[i] << "\t" << proteinNames_[i] << endl;
	cout << endl;
		
}

void DbFile::printTags() const
{
	cout << endl << "Tags: " << allTags_.size() << endl;
	for (set<TagWithLocations>::const_iterator it=allTags_.begin(); it != allTags_.end(); it++)
	{
		for (unsigned int i=0; i<it->tagLength; i++)
			cout << it->tag[i];
		cout << "\t" << it->numLocations << ":";
		for (unsigned int i=0; i<it->numLocations; i++)
			cout << " " << locationArray_[it->locationsOffset+i];
		cout << endl;
	}
}


size_t DbFile::getLocationPointers(const TagWithLocations& twl, vector<const char*>& ptrs) const
{
	ptrs.resize(twl.numLocations, 0);
	if (twl.numLocations == 0)
		return 0;
	for (size_t i=0; i<twl.numLocations; i++)
	{
		assert( locationArray_[twl.locationsOffset + i] < aaArray_.size() );
		ptrs[i] = &aaArray_[0] + locationArray_[twl.locationsOffset + i];
		assert(i==0 || ptrs[i]>ptrs[i-1]);
	}
	return (ptrs.size());
}

string DbFile::makeLocationString(const char* ptr) const
{
	const char *aa0 = &aaArray_[0];
	const char* aaLast = aa0 + aaArray_.size();
	const unsigned int location = getLocation(ptr);
	size_t proteinIdx=0, offset=0;
	computeProteinIdxAndOffsetFromLocation(location, proteinIdx, offset);
	ostringstream oss;
	oss << partIdx_ << ":" << proteinIdx << ":" << offset;
	return oss.str();
}



struct NameAndSequence {
	NameAndSequence() : name(std::string()), sequence(std::string()) {}
	NameAndSequence(const string& n, const string& s) : name(n), sequence(s) {}

	string name;
	string sequence;
};

void createIntersectionFasta(const vector<string>& fastaFiles, size_t minLength, 
							 const char* outputFasta, size_t minNumOccurences)
{
	if (minNumOccurences > fastaFiles.size())
		minNumOccurences = fastaFiles.size();

	if (fastaFiles.size() <= 1)
	{
		cout << "less than two fasta files were supplied, no intersection created!" << endl;
		return;
	}

	// first fasta is the ref, others are also parsed into tags
	const size_t numFastas = fastaFiles.size();
	vector<DbFile> dbFiles(numFastas);

	for (size_t i=0; i<numFastas; i++)
		dbFiles[i].readFastaFile(fastaFiles[i].c_str(), i);

	for (size_t i=1; i<numFastas; i++)
		dbFiles[i].createTags(minLength, minLength);

	const size_t numTemplateProteins = dbFiles[0].getNumProteins();
	vector<TagWithLocations>   tagWithLocations(numFastas);
	vector<NameAndSequence>    intersectionSequences; // holds all sequences that should appear in the intersection file
	map<const string, vector<size_t> > tagMapping; // maps each tag to the sequences that contain it somwhere in them

	for (size_t protIdx=0; protIdx<numTemplateProteins; protIdx++)
	{
		const char* seq = dbFiles[0].getSequence(protIdx);
		const size_t seqLength = (protIdx < numTemplateProteins -1 ?
			dbFiles[0].getProteinStartPointer(protIdx+1) - dbFiles[0].getProteinStartPointer(protIdx) :
			dbFiles[0].getAaArraySize() -  dbFiles[0].getProteinStartPointer(protIdx) );

		if (seqLength<minLength)
			continue;

		for (size_t startIdx=0; startIdx<=seqLength-minLength; startIdx++)
		{
			// check if tag exists in all other dbs
			const char* const templateSeq = (seq +startIdx);
			if (! dbFiles[0].checkAaPtr(templateSeq))
			{
				cout << "numTemplateProteins: " << numTemplateProteins << endl;
				cout << "Problem:           : " << protIdx << "\t" << startIdx << endl;
				cout << "TemplateSeq        : " << templateSeq << endl;
				cout << dbFiles[0].getProteinStartPointer(protIdx+1) << " ::: " << dbFiles[0].getProteinStartPointer(protIdx) << endl;
				cout << endl;
				exit(0);
			}
			assert(dbFiles[0].checkAaPtr(templateSeq));
			size_t numWithTag=1;
			for (size_t j=1; j<numFastas; j++)
			{
				tagWithLocations[j].tag = templateSeq;
				tagWithLocations[j].tagLength = minLength;
				if (dbFiles[j].getTagWithLocations(tagWithLocations[j]))
					numWithTag++;
			}
			if (numWithTag<minNumOccurences)
				continue;

			// extend sequences according to the seq of the firs fasta (it is the template)
			size_t maxCommonLength = minLength;

			vector< vector<const char*> > seqPointers(numFastas);
			seqPointers[0].resize(1,templateSeq); // look at only one position in the template fasta
			for (size_t i=1; i<numFastas; i++)
				dbFiles[i].getLocationPointers(tagWithLocations[i], seqPointers[i]);

			vector<const char*> bestExtensions(numFastas,0);
			size_t bestWithTag=1;
			bestExtensions[0]=templateSeq;
			for (size_t i=1; i<numFastas; i++)
				if (seqPointers[i].size()>0)
				{
					bestExtensions[i] = seqPointers[i][0];
					bestWithTag++;
				}
			assert( bestWithTag >= minNumOccurences);

			while (templateSeq[maxCommonLength] != '\0')
			{
				size_t numFastasInPlay=0;
				for (size_t i=0; i<numFastas; i++)
					if (seqPointers[i].size()>0)
					{
						bestExtensions[i] = seqPointers[i][0];
						numFastasInPlay++;
					}
				assert(numFastasInPlay>=minNumOccurences);

				// remove pointers that don't exted according to the template
				const char nextTemplateChar = templateSeq[maxCommonLength];
				for (size_t i=1; i<numFastas; i++)
				{
					if (seqPointers[i].size()>0)
					{
						for (size_t j=0; j<seqPointers[i].size(); j++)
						{
							assert(seqPointers[i][j]);
							if (*(seqPointers[i][j]+maxCommonLength) != nextTemplateChar)
								seqPointers[i][j]=0;
						}
						seqPointers[i].erase(remove(seqPointers[i].begin(),seqPointers[i].end(), static_cast<const char* >(0)), seqPointers[i].end());
					}
				}

				// check that we still have enough fastas that extend correctly
				numFastasInPlay=1;
				for (size_t i=1; i<numFastas; i++)
					if (seqPointers[i].size()>0)
						numFastasInPlay++;
				
				if (numFastasInPlay>=minNumOccurences)
				{
					maxCommonLength++;
					for (size_t i=0; i<numFastas; i++)
					{
						bestExtensions[i]=0;
						if (seqPointers[i].size()>0)
						{
							bestExtensions[i]=seqPointers[i][0];
							assert( dbFiles[i].checkAaPtr(bestExtensions[i]) );
						}
					}
				}
				else
					break;
			}

			// verify that the best extensions are ok
			bool allSeqsOk = true;
			for (size_t i=0; i<numFastas; i++)
				if (strncmp(templateSeq, bestExtensions[i], maxCommonLength) || 
					! dbFiles[i].checkAaPtr(bestExtensions[i]))
					allSeqsOk = false;

			if (allSeqsOk == false)
			{
				cout << "Prot  idx   = " << protIdx << endl;
				cout << "Start idx   = " << startIdx << endl;
				printf( "TemplateSe  = %p\n",templateSeq);

				for (size_t k=0; k<numFastas; k++)
				{
					cout << "SEQ: ";
					for (size_t j=0; j<maxCommonLength; j++)
						cout << *(bestExtensions[k]+j);
					cout << endl;
					printf("%p\t%p\t%p\n",bestExtensions[k],dbFiles[k].getFirstAAPtr(),dbFiles[k].getLastAAPtr());
				}
				exit(0);
			}
			assert( allSeqsOk );

			// make seqeuence name
			string locationsStr= dbFiles[0].makeLocationString(templateSeq);
			for (size_t i=1; i<numFastas; i++)
				if (bestExtensions[i])
					locationsStr += std::string(" ") + dbFiles[i].makeLocationString(bestExtensions[i]);

			
			// make sure this is not a substring of any other seqeunce 
			const string newSequence = std::string(templateSeq, maxCommonLength);
			const string firstTagSequence = std::string(templateSeq, minLength);
		
			map<const string, vector<size_t> >::iterator it = tagMapping.find(firstTagSequence);
			bool sequeneAlreadyThere = false;
			if (it != tagMapping.end()) 
			{
				for (size_t i=0; i<it->second.size(); i++)
				{
					const size_t seqIdx =  it->second[i];
					if (seqIdx >= intersectionSequences.size())
						continue;

					const string& existingSeq = intersectionSequences[seqIdx].sequence;
					if (existingSeq.length() >= newSequence.length() && existingSeq.find(newSequence) != string::npos)
					{
						sequeneAlreadyThere = true;
						break;
					}
				}
			}
			if (sequeneAlreadyThere)
				continue;

			// remove any seqeunce contained in this sequence (check all tags that are substrings of the new sequence)
			for (size_t startIdx = 0; startIdx<=newSequence.length()-minLength; startIdx++)
			{
				const string tag = std::string(templateSeq + startIdx, minLength);
				map<const string, vector<size_t> >::iterator it = tagMapping.find(tag);
				if (it == tagMapping.end())
					continue;
				
				vector<size_t>& idxs = it->second;
				for (size_t i=0; i<idxs.size(); i++)
					if (idxs[i]<MAX_SIZE_T &&
						intersectionSequences[idxs[i]].sequence.length()<newSequence.length() &&
						newSequence.find(intersectionSequences[idxs[i]].sequence) != string::npos)
					{
						intersectionSequences[idxs[i]].name=std::string();  // remove
						intersectionSequences[idxs[i]].sequence=std::string();
						idxs[i]=MAX_SIZE_T;
					}
				idxs.erase(remove(idxs.begin(),idxs.end(), MAX_SIZE_T), idxs.end());
			}

			// add it
			const size_t idxInVec = intersectionSequences.size();
			intersectionSequences.push_back(NameAndSequence(locationsStr, newSequence));

			// add pointers to all inner tags
			for (size_t startIdx = 0; startIdx<=newSequence.length()-minLength; startIdx++)
			{
				const string tag = std::string(templateSeq + startIdx, minLength);
				tagMapping[tag].push_back(idxInVec);
			}
		}
	}

	// write fasta file
	ofstream ofs(outputFasta);
	if (! ofs.good())
		error("Could not open intersection fasta for writing: ",outputFasta);

	for (size_t i=0; i<intersectionSequences.size(); i++)
	{
		if (intersectionSequences[i].sequence.length()<minLength)
			continue;
		ofs << ">" << intersectionSequences[i].name << endl;
		ofs << intersectionSequences[i].sequence << endl << endl;
	}
	ofs.close();
}



