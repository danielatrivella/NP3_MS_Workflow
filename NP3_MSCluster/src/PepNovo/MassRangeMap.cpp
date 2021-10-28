#include "Config.h"


void MassRangeMap::clear(mass_t max_mass)
{ 
	ranges.clear(); 
	max_map_mass = max_mass; 

	ranges.insert(MASS_T_MAP::value_type(2*NEG_INF,NEG_INF));
	ranges.insert(MASS_T_MAP::value_type(max_map_mass,max_map_mass+1));
	ranges.insert(MASS_T_MAP::value_type(max_map_mass*10,20*max_map_mass));
}

/****************************************************
Since the initial range being added is always smaller or 
equal to the ranges in the map, there can be 5 cases for 
the added range:
*****************************************************/
void MassRangeMap::add_range(mass_t min_mass, mass_t max_mass)
{
	MASS_T_MAP::iterator it_left_min,it_right_min;

	if (min_mass >= max_mass || max_mass> max_map_mass)
		return;


	// get iterator values

	it_right_min=ranges.lower_bound(min_mass);
	if (it_right_min->first == min_mass)
		return;

	

	if (it_right_min != ranges.begin())
	{
		it_left_min = it_right_min;
		it_left_min--;
	}
	else
		it_left_min = ranges.end();


	// just add it - case 4
	if (it_right_min->first > max_mass && it_left_min->second <min_mass)
	{
		ranges.insert(MASS_T_MAP::value_type(min_mass,max_mass));
		return;
	}

	// case 1 - connects to gapped ranges
	if (it_left_min != ranges.end() && it_left_min->second > min_mass && it_right_min->first<max_mass)
	{
		it_left_min->second = it_right_min->second;
	
		ranges.erase(it_right_min);

		return;
	}

	// case 2 - extends some range to the right
	if (it_left_min->first < min_mass && 
		it_left_min->second >= min_mass && it_left_min->second <max_mass)
	{
		it_left_min->second = max_mass;
		return;
	}

	// case 3 - extends some range to the left
	if (it_right_min->first > min_mass && it_right_min->first <= max_mass)
	{
		mass_t new_max = it_right_min->second;
		ranges.erase(it_right_min);
		ranges.insert(MASS_T_MAP::value_type(min_mass, new_max));
		return;
	}


}


// adds a new set of ranges *shift_size* away from the current existing ranges
void MassRangeMap::add_shifted_ranges(mass_t shift_size)
{
	MASS_T_MAP::const_iterator it;

	MassRangeMap new_map = *this;
	it=ranges.begin();
	it++;
	int c=0;
	for ( ; it != ranges.end(); it++)
	{
		mass_t min_mass = it->first  + shift_size;
		mass_t max_mass = it->second + shift_size;
	
		new_map.add_range(min_mass, max_mass);
	}
	
	*this = new_map;
}



void MassRangeMap::print(ostream& os) const
{
	MASS_T_MAP::const_iterator it;

	for (it=ranges.begin(); it != ranges.end(); it++)
		os << it->first << " " << it->second << endl;
}


void MassRangeMap::read_ranges(istream& is)
{
	char buff[32];
	is.getline(buff,32);
	if (strcmp(buff,"#AA_COMBO_RANGES"))
	{
		cout << "Error: expecting \"#AA_COMBO_RANGES\" : " << buff << endl;
		exit(1);
	}

	max_map_mass=-1;
	ranges.clear();

	is.getline(buff,32);
	istringstream iss(buff);
	iss >> max_map_mass;

	if (max_map_mass<0)
	{
		cout << "Error reading max map mass!" << endl;
		exit(1);
	}

	while( ! is.eof())
	{
		is.getline(buff,32);
		if (! strcmp(buff,"#END_AA_COMBO_RANGES"))
			break;

		istringstream iss(buff);
		mass_t min_mass, max_mass;
		iss >> min_mass >> max_mass;
		ranges.insert(MASS_T_MAP::value_type(min_mass,max_mass));
	}

	was_initialized =1;
}

void MassRangeMap::write_ranges(ostream& os) const
{
	if (! was_initialized)
	{
		cout << "Error: aa combo ranges were not initialized!" << endl;
		exit(1);
	}

	MASS_T_MAP::const_iterator it;
	os << "#AA_COMBO_RANGES"<< endl;
	os << max_map_mass << endl;
	for (it=ranges.begin(); it!= ranges.end(); it++)
		os << fixed << setprecision(4) << it->first << " " 
		<< fixed << setprecision(4) << it->second << endl;
	os << "#END_AA_COMBO_RANGES" << endl;
}

/************************************************************
// initializes the allowed_prefix_masses map
// and the allowed suffix masses map
*************************************************************/
void Config::init_allowed_node_masses(mass_t max_mass)
{
	const int num_repeats = (int)(max_mass / 57.0);
	const mass_t two_tolerance = tolerance * 2;
	const vector<mass_t>& aa2mass = session_tables.get_aa2mass();
	allowed_node_masses.clear(max_mass);
	allowed_node_masses.add_range(-tolerance,tolerance);

	int r;
	for (r =0; r<num_repeats; r++)
	{
	
		int a;
		for (a=0; a<session_aas.size(); a++)
		{
			allowed_node_masses.add_shifted_ranges(aa2mass[session_aas[a]]);	
		}

	}

	allowed_node_masses.set_was_initialized(1);

//	allowed_node_masses.print();
}


