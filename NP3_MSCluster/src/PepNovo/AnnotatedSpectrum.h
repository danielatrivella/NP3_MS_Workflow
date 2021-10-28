#ifndef __ANNOTATEDSPECTRUM_H__
#define __ANNOTATEDSPECTRUM_H__

#include "Spectrum.h"
#include "BasicDataStructs.h"

struct PeakAnnotation {
	string label;
	int breakage_idx;
	int frag_type_idx;
};

/***************************************************************
Finds all the peaks in the given spectrum that support a breakage
at a certain location.
****************************************************************/
void annotate_breakage(Spectrum *spec,
					   mass_t pm_with_19, 
					   int peptide_size_idx,
					   Breakage& breakage);




class AnnotatedSpectrum : public Spectrum {
public:

	AnnotatedSpectrum() {}
	AnnotatedSpectrum(PeakList& pl) { copyFromPeakList(pl); }

	vector<Breakage>& get_non_const_breakages() { return breakages; }
	const vector<Breakage>& get_breakages() const { return breakages; }
	const vector< vector< PeakAnnotation> >& get_peak_annotations() const { return peak_annotations; }

	void extract_annotated_intens_and_masses(vector< vector<intensity_t> >& intens,
											 vector< vector<mass_t> >&      masses) const;


	// annoatates spectrum using the spectrum's config
	void annotate_spectrum(mass_t pm_with_19,
						   int spec_charge =0,  // if 0 use spectrum charge
						   bool reset_annotations = false);

	int get_num_annotated_peaks() const;

	
	bool has_stretch_of_b_or_y(int min_stretch_length = 8, int max_skip =1);

	void print_annotated(ostream& os = cout) const;

	// how many of the expected fragment peaks were observed
	void get_number_observed_frags(const vector<int>& frag_types, int& num_obs, int &num_exp) const;

	int get_num_observed_frags(int frag_idx) const;


	float get_explianed_intensity() const;


	// chooses the charge that gives a good_pm_with_19
	void set_charge_from_seq_and_m_over_z();

	
private:
	vector< Breakage > breakages;
	vector< vector< PeakAnnotation > > peak_annotations;


};

void print_dataset_spectra_by_stats(Config *config, char *mgf_file);

#endif
