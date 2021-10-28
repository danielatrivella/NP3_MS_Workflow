#ifndef __FRAGMENTMODEL_H__
#define __FRAGMENTMODEL_H__

#include "Spectrum.h"
#include "ME_REG.h"



struct FragStats {
	FragStats() : frag_idx(NEG_INF), is_viz(false), has_intensity(false), peak_idx(NEG_INF),
					   mass(NEG_INF), log_intensity(NEG_INF), log_local_rank(NEG_INF), log_global_rank(NEG_INF) {};

	void fill_from_breakage(const Breakage *breakage, Spectrum *spec, int f)
	{
		frag_idx = f;
		if (breakage->is_frag_type_visible(f))
		{
			is_viz=true;
			const int pos = breakage->get_position_of_frag_idx(f);
			if (pos>=0)
			{
				has_intensity=true;
				peak_idx = breakage->fragments[pos].peak_idx;

				const Peak& peak = spec->getPeak(peak_idx);
				mass		    = peak.mass;
				iso_level	    = spec->get_peak_iso_level(peak_idx);
				log_intensity   = spec->get_peak_log_intensity(peak_idx);
				log_local_rank  = spec->getPeakLogLocalRank(peak_idx);
				log_global_rank = log(1.0+ static_cast<float>(spec->get_peak_rank(peak_idx)));
			}
		}
	}

	int frag_idx;
	bool is_viz;
	bool has_intensity;
	int    peak_idx;
	mass_t mass;
	float  iso_level;
	float  log_intensity;
	float  log_local_rank;
	float  log_global_rank;
};

/***************************************************************************
Virtual class.
Holds common elements to derived classes: strong fragment and regular fragment.
****************************************************************************/
class FragmentModel {
public:
	FragmentModel() : config_(NULL), tolerance(NEG_INF), frag_tolerance(NEG_INF), exact_peak_tolerance(NEG_INF), 
		one_over_tolerance(NEG_INF),  ind_has_models(0), model_frag_idx(NEG_INF), model_frag_charge(NEG_INF), 
		inten_log_scaling_factor(0), no_inten_log_scaling_factor(0) {}

	virtual void fill_combo_vectors(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							const vector<BreakageInfo>& infos,
							vector< ME_Regression_Sample > & samples) const =0;

	virtual void fill_single_frag_vector(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							BreakageInfo& info,
							vector< fval >& f_vals) const =0;

	void set_config_and_tolerance(const Config *config)
	{
		config_ = config;
		tolerance = config_->getTolerance();
		frag_tolerance = (tolerance<0.1 ? tolerance : tolerance * 0.6);
		exact_peak_tolerance = (tolerance<0.1 ? 0.5 * tolerance : 0.3 *tolerance);
		one_over_tolerance = 1.0 / tolerance;
	}


	virtual bool read_model(const Config *config, istream& is, bool silent_ind) =0;

	virtual bool write_model (ostream& os) const =0;

	int get_model_frag_idx() const { return model_frag_idx; }

protected:
	const Config *config_;
	mass_t tolerance;
	mass_t frag_tolerance;
	mass_t exact_peak_tolerance;
	mass_t one_over_tolerance;

	int ind_has_models;

	int model_frag_idx;
	int model_frag_charge;

	score_t inten_log_scaling_factor;
	score_t no_inten_log_scaling_factor;

	ME_Regression_Model inten_model;
	ME_Regression_Model no_inten_model;
};


#endif

