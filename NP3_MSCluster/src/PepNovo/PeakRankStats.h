#ifndef __PEAKRANKSTATS_H__
#define __PEAKRANKSTATS_H__

#include "PeakRankModel.h"


//void make_proton_mobility_summary(FileManager& fm);


void center_cleavage_reports(string frag_label, int charge);

void aa_composition_stats(const vector<TrainingPeptide>& all_tps, Config *config);

void fragment_detection_stats(const vector<TrainingPeptide>& all_tps, Config *config);

void n_terminal_cleavage_reports(string frag_label, int cut_idx);

void c_terminal_cleavage_reports(string frag_label, int cut_idx_offset);

void proline_cleavage_reports(string frag_label, int charge);


#endif

