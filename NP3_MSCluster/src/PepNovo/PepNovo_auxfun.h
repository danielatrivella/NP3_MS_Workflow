#ifndef __PEPNOVO_AUXFUN_H__
#define __PEPNOVO_AUXFUN_H__

#include "PepNovo_includes.h"
#include "../Common/auxfun.h"

int parseIntFromXml(char* AttributeString);

mass_t parseMassFromXml(char* AttributeString);

mass_t ppm_val(mass_t offset, mass_t total_mass);

void add_to_mass_vector(vector<mass_t>& vec, mass_t val, mass_t tolerance);

int getFileExtensionType(const char* fileName);

#endif
