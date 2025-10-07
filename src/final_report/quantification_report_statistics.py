#!/usr/bin/python

import pandas as pd
import numpy as np
from pathlib import Path
import sys


def compute_quantification_report_statistics(clean_table_file, output_path, mz_tolerance=0.025):
	clean_table_file = Path(clean_table_file)
	output_path = Path(output_path)
	
	if not clean_table_file.exists() or not clean_table_file.is_file():
		sys.exit("The provided path to the clean data file does not exists. Quantification report statistics aborted.")
	if not output_path.exists() or not output_path.is_dir():
		sys.exit("The provided quantification report output path does not exists. Quantification report statistics aborted.")
	
	clean_data = pd.read_csv(clean_table_file)
	n = clean_data.shape[0]
	
	print("  - Computing the quantification statistics\n")
	# create dictionary to store the quantification statistics of the job
	quantification_statistics = {'Statistics': [],
	                       'Value': [],
	                       'Description': []}
	
	quantification_statistics['Statistics'].append("Total number of m/zs")
	quantification_statistics['Value'].append(str(n))
	quantification_statistics['Description'].append("Total number of m/zs in the final table, no filter.")
	
	# if any blank sample, compute their statistics
	if "BLANKS_TOTAL" in clean_data.columns:
		number_not_blank_mzs = sum((clean_data.BLANKS_TOTAL == 0))
		quantification_statistics['Statistics'].append("Total number of not blank m/zs")
		quantification_statistics['Value'].append(f"{number_not_blank_mzs} ({number_not_blank_mzs / n * 100:.1f}%)")
		quantification_statistics['Description'].append(
			"Total number of not blank m/zs in the final table, BLANKS_TOTAL == 0.")
		number_blank_mzs = sum((clean_data.BLANKS_TOTAL > 0))
		quantification_statistics['Statistics'].append("Total number of blank m/zs")
		quantification_statistics['Value'].append(f"{number_blank_mzs} ({number_blank_mzs / n * 100:.1f}%)")
		quantification_statistics['Description'].append("Total number of blank m/zs in the final table, BLANKS_TOTAL > 0.")
		
	# if any bed sample, compute their statistics
	if "BED_TOTAL" in clean_data.columns:
		number_not_bed_mzs = sum((clean_data.BED_TOTAL == 0))
		quantification_statistics['Statistics'].append("Total number of not bed m/zs")
		quantification_statistics['Value'].append(f"{number_not_bed_mzs} ({number_not_bed_mzs / n * 100:.1f}%)")
		quantification_statistics['Description'].append(
			"Total number of not bed (culture media) m/zs in the final table, BED_TOTAL == 0.")
		number_bed_mzs = sum((clean_data.BED_TOTAL > 0))
		quantification_statistics['Statistics'].append("Total number of bed m/zs")
		quantification_statistics['Value'].append(f"{number_bed_mzs} ({number_bed_mzs / n * 100:.1f}%)")
		quantification_statistics['Description'].append(
			"Total number of bed (culture media) m/zs in the final table, BED_TOTAL > 0.")
	
	# Compute number of isotope ions
	number_isotope_ions = clean_data.isotope_ion.sum()
	quantification_statistics['Statistics'].append("Number of isotope m/zs")
	quantification_statistics['Value'].append(
		f"{number_isotope_ions} ({number_isotope_ions / n * 100:.1f}%)")
	quantification_statistics['Description'].append(
		"Total number of isotope m/zs that were assigned as [M+1]+ (isotope_ion == 1). And its percentage over the total number of m/zs.")
	
	# UNPD identification statistics for spectra identification rate
	if 'tremolo_UNPD_category_best' in clean_data.columns:
		number_identified_unpd = (clean_data.tremolo_UNPD_category_best != "out").sum()
		quantification_statistics['Statistics'].append("Number of m/zs identified in UNPD and spectral identification rate")
		quantification_statistics['Value'].append(
			f"{number_identified_unpd} ({number_identified_unpd / n * 100:.1f}%)")
		quantification_statistics['Description'].append(
			"Total number of m/zs that were identified against the UNPD using tremolo and passed the curation (tremolo_UNPD_category_best != 'out'). And its percentage over the total number of m/zs - the spectral identification rate.")
		# unique identifications
		number_unique_unpd_curated = clean_data.tremolo_SMILES_best[(clean_data.tremolo_UNPD_category_best != "out")].unique().size
		quantification_statistics['Statistics'].append("Number of unique UNPD identifications")
		quantification_statistics['Value'].append(
			f"{number_unique_unpd_curated}")
		quantification_statistics['Description'].append(
			"Total number of unique UNPD identifications that passed the curation (unique tremolo_SMILES_best with tremolo_UNPD_category_best != 'out') - unique identified molecules")
	
	# number of pre processed raw MS2 ions
	number_raw_ms2_ions = clean_data.numSpectra.sum()
	quantification_statistics['Statistics'].append("Total number of pre processed raw MS2 ions")
	quantification_statistics['Value'].append(
		f"{number_raw_ms2_ions}")
	quantification_statistics['Description'].append(
		"Total number of pre processed raw MS2 ions (sum of numSpectra).")
	
	# number of unique m/zs within tolerance
	number_unique_mzs = np.round(clean_data.mzConsensus / (mz_tolerance)).unique().size
	quantification_statistics['Statistics'].append("Number of unique m/zs withing tolerance")
	quantification_statistics['Value'].append(
		f"{number_unique_mzs}")
	quantification_statistics['Description'].append(
		"Number of unique m/zs withing tolerance (mzConsensus +- mz_tolerance)")
	
	# number of putative isomers
	describe_isomers_mzs = (clean_data.mzConsensus / (mz_tolerance)).value_counts().describe()
	quantification_statistics['Statistics'].append("Average number of putative isomers by m/z")
	quantification_statistics['Value'].append(
		f"{describe_isomers_mzs['mean']:.2} +- {describe_isomers_mzs['std']:.2}")
	quantification_statistics['Description'].append(
		"Mean number of m/zs within tolerance and its standard deviation - putative average number of isomers (same m/z, different retention times)")
	
	# number of putative fragmented clusters
	number_fragmented_clusters = (clean_data.fragmented_clusters < 0).sum()
	number_mzs_w_fragmented_clusters = (clean_data.fragmented_clusters > 0).sum()
	quantification_statistics['Statistics'].append("Number of m/zs that had at least one fragmented cluster")
	quantification_statistics['Value'].append(
		f"{number_mzs_w_fragmented_clusters} ({number_mzs_w_fragmented_clusters / n * 100:.1f}%)")
	quantification_statistics['Description'].append(
		"Number of m/zs that had at least one putative fragmented cluster within the RT and m/z tolerance and percentage overall m/zs - which are clusters from the same ion that kept separated - number of m/zs with mass dissipation.")
	quantification_statistics['Statistics'].append("Number of putative fragmented cluster")
	quantification_statistics['Value'].append(
		f"{number_fragmented_clusters} ({number_fragmented_clusters / n * 100:.1f}%)")
	quantification_statistics['Description'].append(
		"Number of m/zs that were assigned as putative fragmented clusters from the same ion and percentage overall m/zs - number of putative mass dissipation - this could also be counting minor numerical coincidences. This measures the quality of the clustering (very low < 5% mass dissipation is better). For blank m/zs, only its precursor m/z accounts for the match.")
	
	quantification_statistics_table = pd.DataFrame(quantification_statistics)
	quantification_statistics_table.to_csv(output_path / "quantification_statistics.csv", index=False)

