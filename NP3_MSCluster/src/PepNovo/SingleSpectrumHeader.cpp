#include "SingleSpectrumHeader.h"
#include "BasicDataStructs.h"
#include "PepNovo_auxfun.h"





bool SingleSpectrumHeader::scanSpectrumHeader(FILE* stream, const Config* config)
{
	if (fileType_ == IFT_DTA)
		return (scanDtaSpectrumHeader(stream, config));

	if (fileType_ == IFT_MGF)
		return (scanMgfSpectrumHeader(stream, config));

	if (fileType_ == IFT_MZXML)
		return (scanMzxmlSpectrumHeader(stream, config));

	if (fileType_ == IFT_MS2)
		return (scanMs2SpectrumHeader(stream, config));

	error("File type not supported for scanning from stream: ",fileType_);

	return false;
}

bool	SingleSpectrumHeader::scanSpectrumHeaderFromBuffer(const char* buffer, const Config *config)
{
	// initialize fields that might not get assigned when reading
	scanNumber_ = MIN_INT;
	spectraFileIndexInList_ = MIN_INT;
	
	if (fileType_ == IFT_DAT)
		return (scanDatSpectrumHeaderFromBuffer(buffer, config));

	error("File type not supported for scanning from buffer: ",fileType_);
	return false;
}




bool SingleSpectrumHeader::scanMzxmlSpectrumHeader(FILE *stream, const Config *config)
{
	error("Trying to scan mzXML header, should not reach this funciton!");
	return false;
}




// assumes the stream is pointing to the begining of the header
bool	SingleSpectrumHeader::scanDatSpectrumHeaderFromBuffer(const char* buffer, const Config *config)
{
	assert(fileType_ == IFT_DAT);

	const unsigned int* ui = reinterpret_cast<const unsigned int*>(buffer);
	unsigned int spectrumSize = *ui++;
	originalNumPeaks_ = *ui++;
	
	const mass_t* mt = reinterpret_cast<const mass_t*>(ui);
	mOverZ_ = *mt++;
	originalPmWith19_ = *mt++;
	pmWith19_		  = *mt++;
	firstPeakMass_    = *mt++;

	const short* sh = reinterpret_cast<const short*>(mt);
	charge_   = *sh++;
	fileType_ = *sh++;
	msLevel_  = *sh++; 

	
	const int* ip = reinterpret_cast<const int*>(sh);
	scanNumber_ = *ip++;
	clusterSize_ = *ip++;
	spectraFileIndexInList_ = *ip++;
	datasetIndex_   		= *ip++;
	
	const float* fp = reinterpret_cast<const float*>(ip);
	retentionTime_ = *fp++;
	precursorIntensity_ = *fp++;
	sqs_ = *fp++;
	// NP3 GOT read rt min and max from dat
	retentionTimeMin_ = *fp++;
	retentionTimeMax_ = *fp++;
	peakArea_ = *fp++;
	
	sh = reinterpret_cast<const short*>(fp);
	// NP3 peak id
	const short peakIdLength = *sh++;
	char* p = reinterpret_cast<char*>(const_cast<short*>(sh));
	if (peakIdLength>0)
	{
		char t = p[peakIdLength];
		p[peakIdLength]='\0';
		peakId_ = p;
		p[peakIdLength] = t;
		p += peakIdLength;
	}
	else
		peakId_ = std::string("fake_");

	sh = reinterpret_cast<const short*>(p);
	const short titleLength = *sh++;
	p = reinterpret_cast<char*>(const_cast<short*>(sh));
	if (titleLength>0)
	{
		char t = p[titleLength];
		p[titleLength]='\0';
		title_ = p;
		p[titleLength] = t;
		p += titleLength;
	}
	else
		title_ = std::string();
	
	sh = reinterpret_cast<const short*>(p);
	const short peptideLength = *sh++;

	p = reinterpret_cast<char*>(const_cast<short*>(sh));
	if (peptideLength>0)
	{
		char t = p[peptideLength];
		p[peptideLength]='\0';
		const char lastChar = p[peptideLength-1];
		if (lastChar == '\n' || lastChar == '\r' || lastChar == '\t' || lastChar == ' ')
			p[peptideLength-1]='\0';
		
		peptideStr_ = p;
		p[peptideLength] = t;
		p += peptideLength;
	}
	else
		peptideStr_ = std::string();

	assert(scanNumber_ == MIN_INT || scanNumber_ >=0);
	assert(originalNumPeaks_>0);
	return true;	
}


bool SingleSpectrumHeader::scanDtaSpectrumHeader(FILE* dtaStream, const Config* config)
{
	char buffer[128];

	assert(fileType_ == IFT_DTA);

	double pm19=0.0;
	int c=0, scan=-1;
	while (fgets(buffer, 128, dtaStream))
	{
		if (buffer[0] == '#' || buffer[0] == ' ' || buffer[0] == '\n')
			continue;

		if (sscanf(buffer,"%lf %d\tscan=%d",&pm19,&c,&scan) >= 2)
			break;
	}

	if (scan>=0)
		scanNumber_ = scan;

	if (feof(dtaStream))
		return false;


	assert(pm19>0.0 && c>0);

	positionInFile_ = ftell(dtaStream);
	precursorIntensity_ = 0.0;
	intensity_t firstPeakIntensity=0.0;

	pmWith19_ = static_cast<mass_t>(pm19);
	charge_ = static_cast<short>(c);
	mOverZ_ = (pmWith19_ + (charge_ -1) * MASS_PROTON) / static_cast<mass_t>(charge_);
	if (! fgets(buffer, 128, dtaStream))
		return false;

	istringstream is(buffer);
	Peak p;
	p.mass = -1.0;
	p.intensity = -1.0;
	is >> p.mass >> p.intensity;

	if (p.mass >0 && p.intensity>=0.0)
	{
		firstPeakMass_ = p.mass;
		firstPeakIntensity = p.intensity;
		precursorIntensity_ = p.intensity;

		originalNumPeaks_ = 1;
		while (fgets(buffer, 128, dtaStream))
		{
			if (buffer[0] < '0' || buffer[0] > '9')
				break;

			istringstream is(buffer);
			Peak p;
			is >> p.mass >> p.intensity;
			precursorIntensity_ += p.intensity;
			originalNumPeaks_++;
		}
		
		if (precursorIntensity_>0.0)
			return true;
	}
	
	return false;
}


bool SingleSpectrumHeader::scanMs2SpectrumHeader(FILE* ms2Stream, const Config* config)
{
	char buffer[128];

	assert(fileType_ == IFT_MS2);

	double pm19=0.0;
	int c=0, scan=-1;
	while (fgets(buffer, 128, ms2Stream))
	{
		if (buffer[0] == '#' || buffer[0] == ' ' || buffer[0] == '\n')
			continue;

		if (buffer[0] == ':')
		{
			size_t len = strlen(buffer);
			if (len>0 && len < 128)
				title_ = std::string(buffer+1, len-2);

			break;
		}
	}

	if (! fgets(buffer, 128, ms2Stream))
		return false;

	if (sscanf(buffer,"%lf %d",&pm19,&c) != 2)
		return false;

	if (feof(ms2Stream))
		return false;


	assert(pm19>0.0 && c>0);

	positionInFile_ = ftell(ms2Stream);
	precursorIntensity_ = 0.0;
	intensity_t firstPeakIntensity=0.0;

	pmWith19_ = static_cast<mass_t>(pm19);
	charge_ = static_cast<short>(c);
	mOverZ_ = (pmWith19_ + (charge_ -1) * MASS_PROTON) / static_cast<mass_t>(charge_);
	if (! fgets(buffer, 128, ms2Stream))
		return false;

	istringstream is(buffer);
	Peak p;
	p.mass = -1.0;
	p.intensity = -1.0;
	is >> p.mass >> p.intensity;

	if (p.mass >0 && p.intensity>=0.0)
	{
		firstPeakMass_ = p.mass;
		firstPeakIntensity = p.intensity;
		precursorIntensity_ = p.intensity;

		originalNumPeaks_ = 1;
		while (fgets(buffer, 128, ms2Stream))
		{
			if (buffer[0] < '0' || buffer[0] > '9')
				break;

			istringstream is(buffer);
			Peak p;
			is >> p.mass >> p.intensity;
			precursorIntensity_ += p.intensity;
			originalNumPeaks_++;
		}
		
		if (precursorIntensity_>0.0)
			return true;
	}
	
	return false;
}




bool SingleSpectrumHeader::scanMgfSpectrumHeader(FILE* mgfStream, const Config* config)
{
	char buffer[256];

	assert(fileType_ == IFT_MGF);

	while (fgets(buffer, 256, mgfStream))
	{
		if (strncmp(buffer,"BEGIN IONS",10))
			continue;
		break;
	}

	positionInFile_ = ftell(mgfStream);
	precursorIntensity_ = 0.0;
	intensity_t firstPeakIntensity=0.0;

	if (positionInFile_<0)
		error("Bad skip position in mgf file! This can often be corrected by running unix2dos (or vice versa if appropriate)");

	charge_ = 0;
	mOverZ_ = -1.0;

	// read header info and first peak
	while (true)
	{
		if( ! fgets(buffer, 256, mgfStream))
			return false;

		if (! strncmp(buffer,"END IONS",7))
		{
			originalNumPeaks_ = 0;
			return true;
		}

		if (! strncmp(buffer,"TITLE=",6) )
		{
			int len = strlen(buffer)-1;
			buffer[len]='\0';
			if (buffer[len-1]=='\r' || buffer[len-1]=='\n' )
				buffer[len-1]='\0';
			
			string titleStr = buffer + 6;
			setTitle(titleStr);

			// see if title includes scan number information.
			// this works only if the title ends with: .xxxx.yyyy.d.dta
			// e.g., MyMSMSData.2000.2000.2.dta
			// if the title has this format then scanNumber is set to xxxx
			if (scanNumber_<0)
			{
				len = title_.length();
				if (len>7 &&
					title_[len-1] == 'a' && title_[len-2] == 't' && title_[len-3]=='d' && 
					title_[len-6] == '.' && title_[len-4]== '.')
				{
					int pos = len-7;
					int numDots = 0;

					while (pos>0 && numDots<2)
					{
						--pos;
						if (title_[pos] == '.')
							++numDots;
					}

					if (numDots == 2)
					{
						string scanString = title_.substr(pos+1,len-7-pos);
						int i;
						for (i=0; i<scanString.length(); i++)
							if (scanString[i] == '.')
							{
								scanString[i]=' ';
								break;
							}
						
						istringstream iss(scanString);
						int scan1=-1, scan2=-1;
						iss >> scan1 >> scan2;
						if (scan1 <= scan2 && scan1>0)
							scanNumber_ = scan1;
					}
				}
				continue;
			}
		}
		else
		if (! strncmp(buffer,"SEQ=",4) )
		{
			peptideStr_ = buffer+4;
			continue;		
		}
		else
		if (! strncmp(buffer,"PEPSEQ=",7) )
		{
			peptideStr_ = buffer+7;
			continue;		
		}
		else
		if (! strncmp(buffer,"SCAN=",5) )
		{
			if (sscanf(buffer+5,"%d",&scanNumber_) != 1)
			{
				cout << "Error: couldn't read scan number from mgf file!" << endl;
				exit(1);
			}
			continue;
		}
		else
		if (! strncmp(buffer,"SCANS=",6) ) // this is the offical MGF field, only the first number is kept
		{
			if (sscanf(buffer+6,"%d",&scanNumber_) != 1)
			{
				cout << "Error: couldn't read scan number!" << endl;
				exit(1);
			}
			continue;
		}
		else
		if (! strncmp(buffer,"RT=",3) )
		{
			if (sscanf(buffer+3,"%f",&retentionTime_) != 1)
			{
				cout << "Error: couldn't read retention_time!" << endl;
				exit(1);
			}
			continue;
		}
		else
		if (! strncmp(buffer,"RTINSECONDS=",12) ) // this is the official MGF field name
		{
			if (sscanf(buffer+12,"%f",&retentionTime_) != 1)
			{
				cout << "Error: couldn't read retention_time!" << endl;
				exit(1);
			}
			continue;
		}
		// GOT NP3 peak profile and area
		else
		if (! strncmp(buffer,"RTMIN=",6) ) // this is the NP3 MGF field name
		{
			if (sscanf(buffer+6,"%f",&retentionTimeMin_) != 1)
			{
				cout << "Error: couldn't read retention_time_min!" << endl;
				exit(1);
			}
			continue;
		}
		else
		if (! strncmp(buffer,"RTMAX=",6) ) // this is the Np3 MGF field name
		{
			if (sscanf(buffer+6,"%f",&retentionTimeMax_) != 1)
			{
				cout << "Error: couldn't read retention_time_max!" << endl;
				exit(1);
			}
			continue;
		}
        else
        if (! strncmp(buffer,"PEAK_AREA=",10) ) // this is the Np3 MGF field name
        {
            if (sscanf(buffer+10,"%f",&peakArea_) != 1)
            {
                cout << "Error: couldn't read peak_area!" << endl;
                exit(1);
            }
            continue;
        }
        if (! strncmp(buffer,"NUM_PEAKS=",10) ) // this is the Np3 MGF field name
        {
            if (sscanf(buffer+10,"%d",&originalNumPeaks_) != 1)
            {
                cout << "Error: couldn't read num_peaks!" << endl;
                exit(1);
            }
            //cout << "num_peaks " << originalNumPeaks_ << endl;
            continue;
        }
		else
		if (! strncmp(buffer,"PEAK_ID=",8) ) // this is the Np3 MGF field name
		{
			int len = strlen(buffer)-1;
			buffer[len]='\0';
			if (buffer[len-1]=='\r' || buffer[len-1]=='\n' )
				buffer[len-1]='\0';

			peakId_ = buffer + 8;

//			if (sscanf(buffer+8,"%s",&peakId_) != 1)
//			{
//				cout << "Error: couldn't read peak_id!" << endl;
//				exit(1);
//			}
			continue;

		}
		else
		if (! strncmp(buffer,"CLUSTER_SIZE=",13) )
		{
			if (sscanf(buffer+13,"%d",&clusterSize_) != 1)
			{
				cout << "Error: couldn't read cluster size!" << endl;
				exit(1);
			}
			// NP3 GOT change clusterSize init when reading mgf
			clusterSize_ = 1;
			continue;
		}
		else
		if (! strncmp(buffer,"PRECURSOR_INTENSITY=",20) )
		{
			if (sscanf(buffer+20,"%f",&precursorIntensity_) != 1)
			{
				cout << "Error: couldn't read cluster size!" << endl;
				exit(1);
			}
			continue;
		}
		else	
		if ( ! strncmp(buffer,"CHARGE=",7))
		{
			int c;
			if (sscanf(buffer,"CHARGE=%d",&c) != 1)
			{
				cout <<  "Error: couldn't read charge!" << endl;
				return false;
			}
			charge_ = static_cast<short>(c);
			//cout << "@@@@@@@@@@@@@ CHARGE " << charge_ << endl;
		}
		else
		if (! strncmp(buffer,"PEPMASS=",8))
		{
			istringstream is(buffer+8);
			is >> mOverZ_;
			
			if (mOverZ_ < 0)
			{
				cout << "Error: reading pepmass:" << mOverZ_ << endl;
				return false;
			}		
		}
		else // is this a peak?
		{
			istringstream is(buffer);
			Peak p;
			p.mass = -1.0;
			p.intensity = -1.0;
			is >> p.mass >> p.intensity;

			if (p.mass >0.0 && p.intensity>0.0)
			{
				firstPeakMass_ = p.mass;
				firstPeakIntensity = p.intensity;
				break;
			}
		}
	}

	if (charge_<0 || mOverZ_<0)
		return false;

	originalPmWith19_ = mOverZ_ * charge_ + MASS_PROTON * (1 - charge_);

	if (originalPmWith19_ < 0)
		originalPmWith19_ = -1.0;

	pmWith19_ = originalPmWith19_;

	originalNumPeaks_=1;

	if (precursorIntensity_ > 0.0)
	{
		while ( fgets(buffer, 256, mgfStream) )
		{
			if (! strncmp(buffer,"END IONS",8) )
				break;
			++originalNumPeaks_;
		}
	}
	else // find precursor intensity by summing up peak mass
	{
		precursorIntensity_ = firstPeakIntensity;
		while ( fgets(buffer, 256, mgfStream) )
		{
			if (! strncmp(buffer,"END IONS",8) )
				break;
			istringstream is(buffer);
			Peak p;
			is >> p.mass >> p.intensity;
			precursorIntensity_ += p.intensity;
			++originalNumPeaks_;
		}
	}

	if (peptideStr_.length()>0)
	{
		// trim white space from peptidStr_
		size_t pos = peptideStr_.length()-1;
		while (pos>=0 && (peptideStr_[pos] == ' ' || peptideStr_[pos] == '\t' ||
						  peptideStr_[pos] == '\r' || peptideStr_[pos] == '\n') )
		{
						 pos--;
		}
		if (pos<peptideStr_.length()-1)
			peptideStr_.erase(pos+1);
	}

	// make sure charge is correct for annotated peptide
	if (peptideStr_.length()>0 && charge_>=0)
	{
		Peptide peptide;
		if (peptide.parseFromString(config, peptideStr_))
		{
			mass_t diff = originalPmWith19_ - peptide.get_mass() - MASS_OHHH;
			if (fabs(diff)>8.0)
			{
				// try and correct charge!
				int c;
				for (c=1; c<=7; c++)
				{
					mass_t newPmWith19 =  mOverZ_ * c + MASS_PROTON * (1 - c);
					diff = fabs(newPmWith19 - peptide.get_mass() - MASS_OHHH);
					if (diff<12.0)
					{
						charge_ = c;
						originalPmWith19_ = newPmWith19;
						return true;
					}
				}
			

				cout << "Warning: MGF sequence mass doesn't add up: " << title_ << " (diff: "
					 << diff << ")" <<  endl;
				cout << "m/z: " << mOverZ_ << " originalPmWith19: " << originalPmWith19_ << endl;
				cout << "Pepitde: " << peptideStr_ << " (" << peptide.get_mass() << ")" << endl;
				cout << "Mass Cys = " << config->get_session_tables().get_aa2mass(Cys) << endl;

				return false;
			}
		}
	}

	// NP3 GOT check if MIN and MAX rt were provided, if not set them with default rt tolerance from config
	if (retentionTimeMin_ == -1 || retentionTimeMax_ == -1)
	{
	    //cout << "@@@@@@@@@@@@@@ NO PEAK PROFILE - Creating it from the rt tol @@@@@@" << endl;
		retentionTimeMin_ = retentionTime_ - config->get_rt_tolerance();
		retentionTimeMax_ = retentionTime_ + config->get_rt_tolerance();
		// peak_width*intensity/2
		peakArea_ = precursorIntensity_;
		peakId_ = "fake_raw_" + scanNumber_;
	}

	return true;
}


void SingleSpectrumHeader::printStats(const Config *config, ostream& os, bool print_endl) const
{
	if (fileType_ == IFT_MGF)
	{
		os << ">> " << spectraFileIndexInList_<< " " << (scanNumber_ > MIN_INT ? scanNumber_ :  indexInFile_) << 
			" " << title_;
	}
	else if (fileType_ == IFT_MZXML)
	{
		os << ">> " << spectraFileIndexInList_<< " " << scanNumber_;
	}
	else if (fileType_ == IFT_DAT)
	{
		os << ">> " << spectraFileIndexInList_ << " " << scanNumber_;
	}	
	else
		os << ">> " << spectraFileIndexInList_ << " " << title_;

	if (peptideStr_.length()>0)
	{
		os << " " << peptideStr_;
//		os << " " << peptide_->get_mass() + MASS_OHHH;
	}

	if (print_endl)
		os << endl;
	
}

size_t  SingleSpectrumHeader::writeHeaderToDatBuffer(char* buffer) const
{
	assert(originalNumPeaks_>0);
	// first write information that is needed for the quick scan
	// write 0 for header size (to be filled in later)
	unsigned int* ui = reinterpret_cast<unsigned int*>(buffer);
	*ui++ = 0;
	*ui++ = originalNumPeaks_;

	mass_t* mt = reinterpret_cast<mass_t*>(ui);
	*mt++ = getMOverZ();

	*mt++ = originalPmWith19_;
	*mt++ = pmWith19_;
	*mt++ =	firstPeakMass_;

	short* sh = reinterpret_cast<short*>(mt);
	*sh++ = charge_;
	*sh++ = static_cast<short>(IFT_DAT);
	*sh++ =	msLevel_;
	
	int* ip = reinterpret_cast<int*>(sh);
	*ip++ = scanNumber_;
	*ip++ = clusterSize_;
	*ip++ = spectraFileIndexInList_;
	*ip++ = datasetIndex_;
	
	float *fp = reinterpret_cast<float*>(ip);
	*fp++ = retentionTime_;
	*fp++ = precursorIntensity_;
	*fp++ = sqs_;
	// NP3 GOT write to dat rt min and max
	*fp++ = retentionTimeMin_;
	*fp++ = retentionTimeMax_;
	*fp++ = peakArea_;

	// NP3 GOT peak ID
	sh = reinterpret_cast<short*>(fp);
	*sh++ = peakId_.length();
	char* p = reinterpret_cast<char*>(sh);
	if (peakId_.length()>0)
	{
		strncpy(p,peakId_.c_str(),peakId_.length());
		p+= peakId_.length();
	}

	sh = reinterpret_cast<short*>(p);
	*sh++ = title_.length();
	p = reinterpret_cast<char*>(sh);
	if (title_.length()>0)
	{
		strncpy(p,title_.c_str(),title_.length());
		p+= title_.length();
	}

	sh = reinterpret_cast<short*>(p);
	*sh++ = peptideStr_.length();
	p = reinterpret_cast<char*>(sh);
	if (peptideStr_.length()>0)
	{
		strncpy(p,peptideStr_.c_str(),peptideStr_.length());
		p+= peptideStr_.length();
	}

	return (p-buffer);
}


// writes header without BEGIN IONS
size_t  SingleSpectrumHeader::writeHeaderToMgfBuffer(char* buffer) const
{
	ostringstream oss;
	oss << fixed << setprecision(NUM_SIG_DIGITS);

	if (title_.length()>0)
	{
		oss << "TITLE=" << title_ << endl;
	}
	else
		oss << "TITLE=spectrum_" << spectraFileIndexInList_ << "_" << scanNumber_ << endl;

	if (peptideStr_.length()>0)
		oss << "SEQ=" << peptideStr_ << endl;

    //  do not output SCANS= because it is inconsistent. Only singletons have a meaningful scan number.
    //	if (scanNumber_ >= 0)
    //		oss << "SCANS=" << scanNumber_ << endl;
    // NP3 GOT changed the SCAN initialization from 0 to 1
	static int scan = 1;
	oss << "SCANS=" << scan++ << endl;

	if (retentionTime_>0.0)
		oss << "RTINSECONDS=" << retentionTime_ << endl;

	// GOT NP3 write to mgf rt min and max
    if (retentionTimeMin_>=0.0)
        oss << "RTMIN=" << retentionTimeMin_ << endl;
    if (retentionTimeMax_>=0.0)
        oss << "RTMAX=" << retentionTimeMax_ << endl;
	
	if (clusterSize_>0)
		oss << "CLUSTER_SIZE=" << clusterSize_ << endl;

	// GOT possible num of peaks
	if (getOriginalNumPeaks() > 0)
	    oss << "NUM_PEAKS=" << originalNumPeaks_ << endl;

	 oss << "SQS=" << sqs_ << endl;

	if (precursorIntensity_>0.0)
		oss << "PRECURSOR_INTENSITY=" << scientific << precursorIntensity_ << endl;
	
	oss << "CHARGE=" << charge_ << "+" << endl;
	
	oss << "PEPMASS=" << fixed << setprecision(NUM_SIG_DIGITS) << mOverZ_ << endl;
		
	const size_t len = oss.str().length();
	memcpy(buffer,oss.str().c_str(), len);

	return len;
}

void SingleSpectrumHeader::printStats(ostream& os, bool print_endl) const
{
	if (fileType_ == IFT_MGF)
	{
		os << ">> " << spectraFileIndexInList_<< " " << (scanNumber_ > MIN_INT ? scanNumber_ :  indexInFile_) << 
			" " << title_;
	}
	else if (fileType_ == IFT_MZXML)
	{
		os << ">> " << spectraFileIndexInList_<< " " << scanNumber_;
	}
	else if (fileType_ == IFT_DAT)
	{
		os << ">> " << spectraFileIndexInList_ << " " << scanNumber_;
	}	
	else if (fileType_ == IFT_DTA)
	{
		os << ">> " << spectraFileIndexInList_ << " " << indexInFile_;
	}
	else
		os << ">> " << spectraFileIndexInList_ << " " << title_;

	if ( peptideStr_.length()>0)
	{
		os << " " << peptideStr_;
	//	os << " " << peptide_->get_mass() + MASS_OHHH;
	}
	os << "\t" << setprecision(4) << fixed << mOverZ_;

	if (print_endl)
		os << endl;	
}





