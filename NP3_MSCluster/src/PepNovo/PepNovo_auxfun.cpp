#include "PepNovo_auxfun.h"


// Parse an int - skip characters until you see digits or -, then read until you 
// see something else.
int parseIntFromXml(char* AttributeString)
{
    char Buffer[256];
    int CharCount;
    //
    if (!AttributeString || !*AttributeString)
    {
        return 0;
    }
    CharCount = 0;
    while ((*AttributeString < '0' || *AttributeString > '9') && *AttributeString != '-')
    {
        if (!*AttributeString || CharCount > 256)
        {
            return 0; // too much non-digit garbage!
        }
        AttributeString++;
    }
    CharCount = 0;
    while (*AttributeString >= '0' && *AttributeString <= '9')
    {
        Buffer[CharCount++] = *AttributeString;
        if (CharCount > 10)
        {
            break;
        }
        AttributeString++;
    }
    Buffer[CharCount] = '\0';
    return atoi(Buffer);
}

mass_t parseMassFromXml(char* AttributeString)
{
    char Buffer[256];
    int CharCount;
    //
    if (!AttributeString || !*AttributeString)
    {
        return 0;
    }
    CharCount = 0;
    while ((*AttributeString < '0' || *AttributeString > '9') && *AttributeString != '-')
    {
        if (!*AttributeString || CharCount > 256)
        {
            return 0; // too much non-digit garbage!
        }
        AttributeString++;
    }
    CharCount = 0;
    while ((*AttributeString >= '0' && *AttributeString <= '9') || *AttributeString == '.')
    {
        Buffer[CharCount++] = *AttributeString;
        if (CharCount > 10)
        {
            break;
        }
        AttributeString++;
    }
    Buffer[CharCount] = '\0';
    return (mass_t)atof(Buffer);
}


mass_t ppm_val(mass_t offset, mass_t total_mass)
{
	return (offset / total_mass) * 1000000.0;
}


void add_to_mass_vector(vector<mass_t>& vec, mass_t val, mass_t tolerance)
{
	int i;
	for (i=0; i<vec.size(); i++)
		if (fabs(vec[i]-val)<tolerance)
			break;
	if (i<vec.size())
		return;
	vec.push_back(val);
	
}



/// parses the type from the file name from the extention, -1 if not recognized
int getFileExtensionType(const char* fileName)
{
	int lastPos = strlen(fileName)-1;

	while (lastPos>0 && 
		   (fileName[lastPos] == '\n' || fileName[lastPos] == '\r' || 
			fileName[lastPos] == '\t' || fileName[lastPos] == '\f') )
		--lastPos;

	if (lastPos>2 &&
		fileName[lastPos-2]=='d' && 
		fileName[lastPos-1]=='t' && 
		fileName[lastPos  ]=='a')
		return IFT_DTA;

	if (lastPos>2 &&
		fileName[lastPos-2]=='m' && 
		fileName[lastPos-1]=='g' && 
		fileName[lastPos  ]=='f')
		return IFT_MGF;

	if (lastPos>4 &&
		fileName[lastPos-4] == 'm' &&
	    fileName[lastPos-3] == 'z' &&
		fileName[lastPos-2]=='X' &&
		fileName[lastPos-1]=='M' && 
		fileName[lastPos  ]=='L')
		return IFT_MZXML;

	if (lastPos>2 &&
		fileName[lastPos-2]=='d' && 
		fileName[lastPos-1]=='a' && 
		fileName[lastPos  ]=='t')
		return IFT_DAT;

	if (lastPos>2 &&
		fileName[lastPos-2]=='m' && 
		fileName[lastPos-1]=='s' && 
		fileName[lastPos  ]=='2')
		return IFT_MS2;

	if (lastPos>2 &&
		fileName[lastPos-2]=='p' && 
		fileName[lastPos-1]=='k' && 
		fileName[lastPos  ]=='l')
		return IFT_PKL;

	if (lastPos>7 &&
		fileName[lastPos-7]=='_' &&
		fileName[lastPos-6]=='d' &&
		fileName[lastPos-5]=='t' &&
		fileName[lastPos-4]=='a' &&
		fileName[lastPos-3]=='.' &&
		fileName[lastPos-2]=='t' && 
		fileName[lastPos-1]=='x' && 
		fileName[lastPos  ]=='t')
		return IFT_DTA;

	if (lastPos>2 &&
		fileName[lastPos-2]=='t' && 
		fileName[lastPos-1]=='x' && 
		fileName[lastPos  ]=='t')
		return IFT_TXT;
	
	if (lastPos>2 &&
		fileName[lastPos-2]=='z' && 
		fileName[lastPos-1]=='i' && 
		fileName[lastPos  ]=='p')
		return IFT_ZIP;

	return -1;
}


