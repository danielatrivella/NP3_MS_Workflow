#!/usr/bin/python

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from itertools import chain


def compute_chemical_report_statistics(clean_table_file, output_path):
	clean_table_file = Path(clean_table_file)
	output_path = Path(output_path)
	
	if not clean_table_file.exists() or not clean_table_file.is_file():
		sys.exit("The provided path to the clean data file does not exists. Chemical report statistics aborted.")
	if not output_path.exists() or not output_path.is_dir():
		sys.exit("The provided chemical report output path does not exists. Chemical report statistics aborted.")
		
	print("  - Computing the chemical statistics\n")
	clean_data = pd.read_csv(clean_table_file)
	n = clean_data.shape[0]
	# if any blank sample, remove blank mzs from this analysis
	if "BLANKS_TOTAL" in clean_data.columns:
		clean_data = clean_data.loc[clean_data.BLANKS_TOTAL == 0, :]
	
	# create dictionary to store the chemical statistics of the job
	chemical_statistics = {'Statistics': [],
	                       'Value': [],
	                       'Description': []}
	
	chemical_statistics['Statistics'].append("Total number of m/zs")
	chemical_statistics['Value'].append(str(n))
	chemical_statistics['Description'].append("Total number of m/zs in the final table, no filter.")
	
	#  number of not blank m/zs
	number_not_blank_mzs = clean_data.shape[0]
	chemical_statistics['Statistics'].append("Number of not blank m/zs")
	chemical_statistics['Value'].append(f"{number_not_blank_mzs} ({number_not_blank_mzs/n*100:.1f}%)")
	chemical_statistics['Description'].append("Total number of m/zs in the final table that do not appear in any blank sample (BLANKS_TOTAL == 0). And its percentage over the total number of m/zs.")
	
	# compute number of putative molecules - protonated m/zs [M+H]+
	number_protonated_mzs = clean_data.protonated_representative.sum()
	chemical_statistics['Statistics'].append("Number of not blank [M+H]+ m/zs")
	chemical_statistics['Value'].append(f"{number_protonated_mzs} ({number_protonated_mzs / number_not_blank_mzs*100:.1f}%)")
	chemical_statistics['Description'].append(
		"Total number of not blank m/zs that were selected as [M+H]+ - protonated (protonated_representative == 1) - putative molecules. And its percentage over the total number of not blank m/zs.")
	
	# filter only [M+H]+
	clean_data = clean_data.loc[clean_data.protonated_representative == 1,:]
	
	# UNPD identification statistics for putative molecules
	# # number of not blank [M+H]+ m/zs identified / UNPD size; potential for chemical novelty of the set (percentage of [M+H]+ not identified
	if 'tremolo_UNPD_category_best' in clean_data.columns:
		total_unpd_fixo = 170602  # total unique UNPD
		number_protonated_identified_unpd = (clean_data.tremolo_UNPD_category_best != "out").sum()
		chemical_statistics['Statistics'].append("Number of [M+H]+ m/zs identified in UNPD")
		chemical_statistics['Value'].append(
			f"{number_protonated_identified_unpd} ({number_protonated_identified_unpd / number_protonated_mzs*100:.1f}%)")
		chemical_statistics['Description'].append(
			"Total number of not blank [M+H]+ m/zs that were identified against the UNPD using tremolo and passed the curation (tremolo_UNPD_category_best != 'out'). And its percentage over the total number of not blank [M+H]+ m/zs.")
		chemical_statistics['Statistics'].append("Number of putative novel [M+H]+ m/zs")
		chemical_statistics['Value'].append(
			f"{(clean_data.shape[0]-number_protonated_identified_unpd)} ({(clean_data.shape[0]-number_protonated_identified_unpd) / number_protonated_mzs * 100:.1f}%)")
		chemical_statistics['Description'].append(
			"Total number of not blank [M+H]+ m/zs that were NOT identified against the UNPD using tremolo (tremolo_UNPD_category_best == 'out'), thus, may represent putative novel compounds not present in the database. And its percentage over the total number of not blank [M+H]+ m/zs.")
		number_protonated_unique_identification_unpd = clean_data.loc[clean_data.tremolo_UNPD_category_best != "out","tremolo_SMILES_best"].unique().size
		chemical_statistics['Statistics'].append("Number of [M+H]+ m/zs unique identifications and UNPD coverage")
		chemical_statistics['Value'].append(
			f"{number_protonated_unique_identification_unpd} ({number_protonated_unique_identification_unpd / total_unpd_fixo*100:.1f}%)")
		chemical_statistics['Description'].append(
			"The number of unique and curated identifications in UNPD for the not blank [M+H]+ m/zs and its percentage over the total number of entries in UNPD (unique tremolo_SMILES_best for tremolo_UNPD_category_best != 'out') - the UNPD coverage by [M+H]+.")
		# chemistry diversity for NP based on superclass annotation, from the annotated and curated data
		total_superclass_unpd = 91 # from the superclass_groupings dictionary in tremolo_UNPD_curate_identification
		# use curated superclass, remove NA and None
		number_unique_superclass = np.unique([x for x in chain.from_iterable(clean_data.curated_superclass.str.split(":", expand=True).to_numpy())
		                                      if x==x and x is not None]).size
		chemical_statistics['Statistics'].append("Chemical diversity of [M+H]+ in UNPD Superclasses")
		chemical_statistics['Value'].append(
			f"{number_unique_superclass} ({number_unique_superclass / total_superclass_unpd*100:.1f}%)")
		chemical_statistics['Description'].append(
			"The unique number of superclasses that got identified by the not blank [M+H]+ m/zs in UNPD. And its percentage over the total number of unique superclasses considered (88 for UNPD using NPClassifier). ")

	# samples redundancy (number of [M+H]+ that appear in more than one sample, mediana+-sd of number of samples that they appear
	clean_data["number_samples"] = (clean_data.loc[:, clean_data.columns.str.endswith("_area")] > 0).sum(1)
	number_protonated_in_redundant_samples = (clean_data.number_samples > 1).sum()
	number_samples_describe = clean_data["number_samples"].describe()
	
	chemical_statistics['Statistics'].append("Number of [M+H]+ m/zs that are redundant")
	chemical_statistics['Value'].append(
		f"{number_protonated_in_redundant_samples} ({number_protonated_in_redundant_samples / number_protonated_mzs*100:.1f}%)")
	chemical_statistics['Description'].append(
		"Number of not blank [M+H]+ m/zs that appear in more than one sample (<sample_code>_area > 0 for at least two samples). And its percentage over the total number of not blank [M+H]+ m/zs.")
	
	chemical_statistics['Statistics'].append("Median number of samples that the [M+H]+ m/zs appear")
	chemical_statistics['Value'].append(
		f"{number_samples_describe['50%']:.1f} +-({number_samples_describe['std']:.1f})")
	chemical_statistics['Description'].append(
		"Median number of samples that the not blank [M+H]+ m/zs appear +- its standard deviation.")
	
	chemical_statistics['Statistics'].append("Mean number of samples that the [M+H]+ m/zs appear")
	chemical_statistics['Value'].append(
		f"{number_samples_describe['mean']:.1f} +-({number_samples_describe['std']:.1f})")
	chemical_statistics['Description'].append(
		"Mean number of samples that the not blank [M+H]+ m/zs appear +- its standard deviation.")
	
	chemical_statistics_table = pd.DataFrame(chemical_statistics)
	chemical_statistics_table.to_csv(output_path/"chemical_statistics_UNPD.csv", index=False)
	
