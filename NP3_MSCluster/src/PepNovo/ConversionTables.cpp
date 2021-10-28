#include "ConversionTables.h"


/*************************************************************
// initializes the tables with values for the 
// ordinary amino acids (unmodified).
**************************************************************/
void ConversionTables::init_for_standard_aas()
{
	int i;

	int max_char = 0xFF;

	
	char2aa.clear();
	char2mass.clear();
	aa2mass.clear();
	aa2label.clear();
	aa_positions.clear();
	
	char2aa.resize(max_char); 
	char2mass.resize(max_char);

    for (i=0; i<max_char; i++)
    {
        char2aa[i]=-999999;
        char2mass[i]=-999999;
    }

	aa2mass.resize(Val+1);
	aa2label.resize(Val+1);
	aa2char.resize(Val+1);

	aa_positions.resize(Val+1,0);

	org_aa.resize(Val+1);
	for (i=0; i<=Val; i++)
		org_aa[i]=i;

//	char2aa['^']=N_TERM; aa2char[N_TERM]='^'; char2mass['^']=0.0; aa2mass[N_TERM]=0.0;
//	char2aa['$']=C_TERM; aa2char[C_TERM]='$'; char2mass['$']=0.0; aa2mass[C_TERM]=0.0;
	char2aa['_']=Gap;  aa2char[Gap]='_';  char2mass['_']= 9999.999;   aa2mass[Gap]=9999.999;
	char2aa['X']=Xle;  aa2char[Xle]='X';  char2mass['X']= 113.08406;  aa2mass[Xle]= 113.08406;  
    char2aa['A']=Ala;  aa2char[Ala]='A';  char2mass['A']=  71.03711;  aa2mass[Ala]=  71.03711;
    char2aa['R']=Arg;  aa2char[Arg]='R';  char2mass['R']= 156.10111;  aa2mass[Arg]= 156.10111;
    char2aa['N']=Asn;  aa2char[Asn]='N';  char2mass['N']= 114.04293;  aa2mass[Asn]= 114.04293; 
    char2aa['D']=Asp;  aa2char[Asp]='D';  char2mass['D']= 115.02694;  aa2mass[Asp]= 115.02694;
    char2aa['C']=Cys;  aa2char[Cys]='C';  char2mass['C']= 103.00919;  aa2mass[Cys]= 103.00919;
    char2aa['Q']=Gln;  aa2char[Gln]='Q';  char2mass['Q']= 128.05858;  aa2mass[Gln]= 128.05858;
    char2aa['E']=Glu;  aa2char[Glu]='E';  char2mass['E']= 129.04259;  aa2mass[Glu]= 129.04259;
    char2aa['G']=Gly;  aa2char[Gly]='G';  char2mass['G']=  57.02146;  aa2mass[Gly]=  57.02146;
    char2aa['H']=His;  aa2char[His]='H';  char2mass['H']= 137.05891;  aa2mass[His]= 137.05891;
    char2aa['I']=Ile;  aa2char[Ile]='I';  char2mass['I']= 113.08406;  aa2mass[Ile]= 113.08406;
    char2aa['L']=Leu;  aa2char[Leu]='L';  char2mass['L']= 113.08406;  aa2mass[Leu]= 113.08406;
    char2aa['K']=Lys;  aa2char[Lys]='K';  char2mass['K']= 128.09496;  aa2mass[Lys]= 128.09496;
    char2aa['M']=Met;  aa2char[Met]='M';  char2mass['M']= 131.04049;  aa2mass[Met]= 131.04049;
    char2aa['F']=Phe;  aa2char[Phe]='F';  char2mass['F']= 147.06841;  aa2mass[Phe]= 147.06841;
    char2aa['P']=Pro;  aa2char[Pro]='P';  char2mass['P']=  97.05276;  aa2mass[Pro]=  97.05276;
    char2aa['S']=Ser;  aa2char[Ser]='S';  char2mass['S']=  87.03203;  aa2mass[Ser]=  87.03203;
    char2aa['T']=Thr;  aa2char[Thr]='T';  char2mass['T']= 101.04768;  aa2mass[Thr]= 101.04768;
    char2aa['W']=Trp;  aa2char[Trp]='W';  char2mass['W']= 186.07931;  aa2mass[Trp]= 186.07931;
    char2aa['Y']=Tyr;  aa2char[Tyr]='Y';  char2mass['Y']= 163.06333;  aa2mass[Tyr]= 163.06333;
    char2aa['V']=Val;  aa2char[Val]='V';  char2mass['V']=  99.06841;  aa2mass[Val]=  99.06841;
    
    // chose to make B=D so
    char2aa['B']=Asp;                     char2mass['B']= 115.088;

	for (i=Gap; i<=Val; i++)
		aa2label[i]=aa2char[i];
}


bool ConversionTables::add_optional_PTM_aa(int aa, const string& label, mass_t delta, int position)
{
	if (aa2mass[aa]+delta<12.00)
	{
		cout << "Warning: illegal mass for PTM " << label << " : " << aa2mass[aa]+delta << " (mass of optional PTM should be above 12 Da)" << endl;
		return false;
	}
	org_aa.push_back(aa);
	aa2mass.push_back(aa2mass[aa]+delta);
	aa2label.push_back(label);
	aa_positions.push_back(position);

	return true;
}

bool ConversionTables::make_fixed_mod(int aa, mass_t delta)
{
	if (aa2mass[aa]+delta<12.00)
		return false;

	aa2mass[aa]+=delta;
	if (aa>0 && aa<=Val)
		char2mass[aa2char[aa]] += delta;

	return true;
}


bool ConversionTables::add_optional_PTM_terminal_aa(mass_t delta, 
													int position, 
													const string& ptm_label)
{
	if (position != -1 && position != 1)
	{
		cout << "Error: terminal PTMs must have position +1 or -1!" << endl;
		exit(1);
	}

	if (delta<-50.0)
	{
		cout << "Warning: illegal mass offset for terminal PTM: " << ptm_label << " (offset should not be less than -50 Da)" << endl;
		return false;
	}

	int aa;
	for (aa=Ala; aa<=Val; aa++)
	{
		org_aa.push_back(aa);

		aa2mass.push_back(aa2mass[aa]+delta);
		aa_positions.push_back(position);

		// make label
		string aa_label = aa2char[aa] + ptm_label;

		aa_label[0] = aa2char[aa];
		aa2label.push_back(aa_label);
	}

	return true;
}



