#!/usr/bin/python

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from itertools import chain
import matplotlib.pyplot as plt
from tremolo_UNPD_curate_identification import superclass_groupings_names # list of superclass groupings names

# clean_table_file must contain the clean count table with peak area quantification
def plot_superclass_samples_distribution(metadata_file, clean_table_file, output_path):
	clean_table_file = Path(clean_table_file)
	output_path = Path(output_path)
	metadata_file = Path(metadata_file)
	
	if not clean_table_file.exists() or not clean_table_file.is_file():
		sys.exit("The provided path to the clean data file does not exists. Superclasses distribution plot aborted.")
	if not output_path.exists() or not output_path.is_dir():
		sys.exit("The provided chemical report output path does not exists. Superclasses distribution plot aborted.")
	if not metadata_file.exists() or not metadata_file.is_file():
		sys.exit("The provided path to the metadata table file does not exists. Superclasses distribution plot aborted.")
	
	
	# read the data
	clean_data = pd.read_csv(clean_table_file)
	metadata = pd.read_csv(metadata_file)
	# if there is a curated identification, proceed for plotting
	if "tremolo_curated_superclass_grouping" in clean_data.columns:
		print("  - Creating the superclass grouping distribution by not blank sample \n")
		# get the columns names containing the count of spectra by peak area without blanks
		samples_area_name = metadata.SAMPLE_CODE[metadata.SAMPLE_TYPE.str.lower() != "blank"].values
		samples_area_col = samples_area_name + "_area"
		
		# group the quantification columns by the superclass grouping and sum the respective rows
		samples_area_by_superclass_grouping = clean_data.groupby("tremolo_curated_superclass_grouping")[samples_area_col].sum()
		# normalize the quantification by superclass
		samples_area_by_superclass_grouping = samples_area_by_superclass_grouping.div(samples_area_by_superclass_grouping.sum(axis=0), axis=1)
		# rename columns with the original samples codes
		samples_area_by_superclass_grouping.columns = samples_area_name
		
		# define the superclasses colors for the plot
		grouping_colors = ['#e8ff00', '#ff8b00', '#ff008b', '#00cc00', '#e800ff', '#5dff00', '#6fffff', '#5d00ff', '#00b9ff', '#002eff', '#cccccc']
		superclass_groupings_colors = dict(zip(superclass_groupings_names, grouping_colors))
		
		fig, ax = plt.subplots(figsize=(20, 10))  # plot size - bigger plot
		
		bottom = pd.Series([0] * len(samples_area_name), index=samples_area_name)
		
		for superclass_group in samples_area_by_superclass_grouping.index:
			distribution_superclass = samples_area_by_superclass_grouping.loc[superclass_group, samples_area_name]
			superclass_color = superclass_groupings_colors.get(superclass_group, 'gray')
			ax.bar(samples_area_name, distribution_superclass, bottom=bottom.values, label=superclass_group, color=superclass_color)
			bottom += distribution_superclass
		
		# Axes e style
		ax.set_ylabel('Normalize distribution by superclass grouping', fontsize=18)
		ax.set_title('Samples Composition by Superclass Grouping (normalized by sample)', fontsize=20)
		ax.set_xticks(range(len(samples_area_name)))
		ax.set_xticklabels(samples_area_name, rotation=45, ha='right', fontsize=16)
		plt.yticks(fontsize=16)
		
		# add Legend
		ax.legend(title='Superclass Grouping', loc='lower center', bbox_to_anchor=(0.5, -0.35),
		          ncol=4, fontsize='14', title_fontsize='14')
		
		plt.tight_layout()
		
		barplot_filepath = output_path / "samples_composition_superclass_grouping_distribution.png"
		plt.savefig(barplot_filepath, dpi=300, bbox_inches='tight')
		#plt.show()


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
	
	if number_not_blank_mzs > 0:
		# compute number of putative molecules - protonated m/zs [M+H]+
		number_protonated_mzs = clean_data.protonated_representative.sum()
		chemical_statistics['Statistics'].append("Number of not blank [M+H]+ m/zs")
		chemical_statistics['Value'].append(f"{number_protonated_mzs} ({number_protonated_mzs / number_not_blank_mzs*100:.1f}%)")
		chemical_statistics['Description'].append(
			"Total number of not blank m/zs that were selected as [M+H]+ - protonated (protonated_representative == 1) - putative molecules. And its percentage over the total number of not blank m/zs.")
		
		if number_protonated_mzs > 0:
			# filter only [M+H]+
			clean_data = clean_data.loc[clean_data.protonated_representative == 1,:]
			
			# UNPD identification statistics for putative molecules
			# # number of not blank [M+H]+ m/zs identified / UNPD size; potential for chemical novelty of the set (percentage of [M+H]+ not identified
			if 'tremolo_UNPD_category_best' in clean_data.columns:
				total_unpd_fixo = 170602  # total unique UNPD
				number_protonated_identified_unpd = (clean_data.tremolo_UNPD_category_best != "out").sum()
				chemical_statistics['Statistics'].append("Number of [M+H]+ m/zs identified in UNPD and spectral identification rate")
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
					"The number of unique and curated identifications in UNPD for the not blank [M+H]+ m/zs and its percentage over the total number of entries in UNPD (unique tremolo_SMILES_best for tremolo_UNPD_category_best != 'out' over 170602 compounds for UNPD total) - the UNPD coverage by [M+H]+.")
				# chemistry diversity for NP based on superclass annotation, from the annotated and curated data
				total_superclass_unpd = 91 # from the superclass_groupings dictionary in tremolo_UNPD_curate_identification
				# use curated superclass, remove NA and None
				number_unique_superclass = np.unique([x for x in chain.from_iterable(clean_data.tremolo_curated_superclass.str.split(":", expand=True).to_numpy())
				                                      if x==x and x is not None and x is not ""]).size
				chemical_statistics['Statistics'].append("Chemical diversity of [M+H]+ in UNPD Superclasses")
				chemical_statistics['Value'].append(
					f"{number_unique_superclass} ({number_unique_superclass / total_superclass_unpd*100:.1f}%)")
				chemical_statistics['Description'].append(
					"The unique number of superclasses that got identified by the not blank [M+H]+ m/zs in UNPD. And its percentage over the total number of unique superclasses considered ("+str(total_superclass_unpd)+" for UNPD using NPClassifier). ")
				total_superclass_npclassifier_grouping = 11
				number_unique_superclass_grouping = np.unique(clean_data.tremolo_curated_superclass_grouping.values).size
				chemical_statistics['Statistics'].append("Chemical diversity of [M+H]+ in UNPD Superclasses grouping")
				chemical_statistics['Value'].append(
					f"{number_unique_superclass_grouping} ({number_unique_superclass_grouping / total_superclass_npclassifier_grouping * 100:.1f}%)")
				chemical_statistics['Description'].append(
					"The unique number of superclasses grouping that got identified by the not blank [M+H]+ m/zs in UNPD. And its percentage over the total number of unique superclasses grouping considered (" + str(
						total_superclass_npclassifier_grouping) + " groups proposed by NP3 from the NPClassifier superclasses). ")
			
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
	
def compute_chemical_identification_report_GNPS_result(clean_table_file, output_path):
	clean_table_file = Path(clean_table_file)
	output_path = Path(output_path)
	
	if not clean_table_file.exists() or not clean_table_file.is_file():
		sys.exit("The provided path to the clean data file does not exists. Chemical and identification report for GNPS result statistics aborted.")
	output_path.mkdir(exist_ok=True)
	if not output_path.exists() or not output_path.is_dir():
		sys.exit("The provided chemical report output path does not exists. Chemical and identification report for GNPS result statistics aborted.")
	
	print("* Computing the chemical and identification statistics for GNPS and UNPDxGNPS best origin *\n")
	clean_data = pd.read_csv(clean_table_file)
	# if any blank sample, remove blank mzs from this analysis
	if "BLANKS_TOTAL" in clean_data.columns:
		clean_data = clean_data.loc[clean_data.BLANKS_TOTAL == 0, :]
	n = clean_data.shape[0]
	if 'gnps_category' not in clean_data.columns or (clean_data.gnps_category == "out").all():
		sys.exit("  - No valid GNPS curated data is available. All gnps_category is out or missing. Aborting chemical identification report for GNPS and best origin.")
	
	# create dictionary to store the chemical and identification statistics of the job for GNPS result
	gnps_statistics = {'Statistics': [],
                       'Value': [],
                       'Description': []}
	
	gnps_statistics['Statistics'].append("Total number of not blank m/zs")
	gnps_statistics['Value'].append(str(n))
	gnps_statistics['Description'].append("Total number of not blank m/zs in the final table (BLANKS_TOTAL == 0).")
	# also create statistics for the UNPDxGNPS - if best origin exist
	# create dictionary to store the chemical and identification statistics of the job for UNPDxGNPS result
	if "best_origin_SMILES" in clean_data.columns:
		best_origin_exists = True
		unpd_gnps_statistics = {'Statistics': [],
		                        'Value': [],
		                        'Description': []}
		
		unpd_gnps_statistics['Statistics'].append("Total number of not blank m/zs")
		unpd_gnps_statistics['Value'].append(str(n))
		unpd_gnps_statistics['Description'].append("Total number of not blank m/zs in the final table (BLANKS_TOTAL == 0).")
	else:
		best_origin_exists = False
	if n > 0:
		# identification statistics for spectra identification rate
		number_identified_gnps = (clean_data.gnps_category != "out").sum()
		gnps_statistics['Statistics'].append("Number of m/zs identified in GNPS and spectral identification rate")
		gnps_statistics['Value'].append(
			f"{number_identified_gnps} ({number_identified_gnps / n * 100:.1f}%)")
		gnps_statistics['Description'].append(
			"Total number of not blank m/zs that were identified against the GNPS libraries (gnps_category != 'out'). And its percentage over the total number of not blank m/zs - the spectral identification rate.")
		if best_origin_exists:
			number_identified_gnps = (clean_data.curated_identification_best_origin == "GNPS").sum()
			number_identified_unpd = (clean_data.curated_identification_best_origin == "UNPD").sum()
			unpd_gnps_statistics['Statistics'].append("Number of m/zs identified in UNPD or GNPS, and spectral identification rate")
			unpd_gnps_statistics['Value'].append(
				f"{number_identified_gnps+number_identified_unpd} ({(number_identified_gnps+number_identified_unpd) / n * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of not blank m/zs that were identified against the UNPD and GNPS libraries - best origin (curated_identification_best_origin != ''). And its percentage over the total number of not blank m/zs - the spectral identification rate.")
			unpd_gnps_statistics['Statistics'].append(
				"Number of m/zs identified in UNPD as best origin, and spectral identification rate")
			unpd_gnps_statistics['Value'].append(
				f"{number_identified_unpd} ({(number_identified_unpd) / n * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of not blank m/zs that were identified against the UNPD and selected as best origin (curated_identification_best_origin == 'UNPD'). And its percentage over the total number of not blank m/zs - the spectral identification rate.")
			unpd_gnps_statistics['Statistics'].append(
				"Number of m/zs identified in GNPS as best origin, and spectral identification rate")
			unpd_gnps_statistics['Value'].append(
				f"{number_identified_gnps} ({(number_identified_gnps) / n * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of not blank m/zs that were identified against the GNPS and selected as best origin (curated_identification_best_origin == 'GNPS'). And its percentage over the total number of not blank m/zs - the spectral identification rate.")
		# unique identifications
		number_unique_gnps_curated = clean_data.gnps_Smiles[
			(clean_data.gnps_category != "out")].unique().size
		gnps_statistics['Statistics'].append("Number of unique GNPS identifications")
		gnps_statistics['Value'].append(
			f"{number_unique_gnps_curated}")
		gnps_statistics['Description'].append(
			"Total number of unique GNPS identifications that passed the curation (unique gnps_Smiles with gnps_category != 'out') - unique identified molecules")
		if best_origin_exists:
			number_unique_unpd_gnps_curated = clean_data.best_origin_SMILES[
				(clean_data.best_origin_SMILES != "")].unique().size
			unpd_gnps_statistics['Statistics'].append("Number of unique UNPD or GNPS identifications")
			unpd_gnps_statistics['Value'].append(
				f"{number_unique_unpd_gnps_curated}")
			unpd_gnps_statistics['Description'].append(
				"Total number of unique UNPD and GNPS identifications that passed the best origin curation (unique best_origin_SMILES with best_origin_SMILES != '') - unique identified molecules")
			number_unique_unpd_curated = clean_data.best_origin_SMILES[
				(clean_data.curated_identification_best_origin == "UNPD")].unique().size
			number_unique_gnps_curated = clean_data.best_origin_SMILES[
				(clean_data.curated_identification_best_origin == "GNPS")].unique().size
			unpd_gnps_statistics['Statistics'].append("Number of unique UNPD best origin identifications")
			unpd_gnps_statistics['Value'].append(
				f"{number_unique_unpd_curated}")
			unpd_gnps_statistics['Description'].append(
				"Total number of unique UNPD identifications selected as best origin in the curation (unique best_origin_SMILES with curated_identification_best_origin == 'UNPD') - unique identified molecules")
			unpd_gnps_statistics['Statistics'].append("Number of unique GNPS best origin identifications")
			unpd_gnps_statistics['Value'].append(
				f"{number_unique_gnps_curated}")
			unpd_gnps_statistics['Description'].append(
				"Total number of unique GNPS identifications selected as best origin in the curation (unique best_origin_SMILES with curated_identification_best_origin == 'GNPS') - unique identified molecules")
		# identification statistics for [M+H]+
		# now compute number of putative molecules - protonated m/zs [M+H]+
		number_protonated_mzs = clean_data.protonated_representative.sum()
		gnps_statistics['Statistics'].append("Number of [M+H]+ m/zs")
		gnps_statistics['Value'].append(
			f"{number_protonated_mzs} ({number_protonated_mzs / n * 100:.1f}%)")
		gnps_statistics['Description'].append(
			"Total number of not blank m/zs that were selected as [M+H]+ - protonated (protonated_representative == 1) - putative molecules. And its percentage over the total number of not blank m/zs.")
		# if any protonated, proceed
		if number_protonated_mzs > 0:
			# filter only [M+H]+
			clean_data = clean_data.loc[clean_data.protonated_representative == 1, :]
			
			# GNPS identification statistics for putative molecules
			# # number of not blank [M+H]+ m/zs identified / GNPS size; potential for chemical novelty of the set (percentage of [M+H]+ not identified
			total_gnps_fixo = 600000 # total unique GNPS
			number_protonated_identified_gnps = (clean_data.gnps_category != "out").sum()
			gnps_statistics['Statistics'].append("Number of [M+H]+ m/zs identified in GNPS")
			gnps_statistics['Value'].append(
				f"{number_protonated_identified_gnps} ({number_protonated_identified_gnps / number_protonated_mzs * 100:.1f}%)")
			gnps_statistics['Description'].append(
				"Total number of not blank [M+H]+ m/zs that were identified against the GNPS and passed the curation (gnps_category != 'out'). And its percentage over the total number of not blank [M+H]+ m/zs.")
			if best_origin_exists:
				number_protonated_identified_gnps_best = (clean_data.curated_identification_best_origin == "GNPS").sum()
				number_protonated_identified_unpd_best = (clean_data.curated_identification_best_origin == "UNPD").sum()
				unpd_gnps_statistics['Statistics'].append("Number of [M+H]+ m/zs identified in UNPD or GNPS")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_identified_unpd_best+number_protonated_identified_gnps_best} ({(number_protonated_identified_unpd_best+number_protonated_identified_gnps_best) / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of not blank [M+H]+ m/zs that were identified against the UNPD or GNPS and passed the best origin curation (curated_identification_best_origin != ''). And its percentage over the total number of not blank [M+H]+ m/zs.")
				unpd_gnps_statistics['Statistics'].append("Number of [M+H]+ m/zs identified in UNPD as best origin")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_identified_unpd_best} ({number_protonated_identified_unpd_best / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of not blank [M+H]+ m/zs that were identified against the UNPD and selected as best origin in the curation (curated_identification_best_origin == 'UNPD'). And its percentage over the total number of not blank [M+H]+ m/zs.")
				unpd_gnps_statistics['Statistics'].append("Number of [M+H]+ m/zs identified in GNPS as best origin")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_identified_gnps_best} ({number_protonated_identified_gnps_best / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of not blank [M+H]+ m/zs that were identified against the GNPS and selected as best origin in the curation (curated_identification_best_origin == 'GNPS'). And its percentage over the total number of not blank [M+H]+ m/zs.")
			# not identified m/zs - novelty
			gnps_statistics['Statistics'].append("Number of putative novel [M+H]+ m/zs in GNPS")
			gnps_statistics['Value'].append(
				f"{(clean_data.shape[0] - number_protonated_identified_gnps)} ({(clean_data.shape[0] - number_protonated_identified_gnps) / number_protonated_mzs * 100:.1f}%)")
			gnps_statistics['Description'].append(
				"Total number of not blank [M+H]+ m/zs that were NOT identified against the GNPS (gnps_category == 'out'), thus, may represent putative novel compounds not present in the database. And its percentage over the total number of not blank [M+H]+ m/zs.")
			if best_origin_exists:
				unpd_gnps_statistics['Statistics'].append("Number of putative novel [M+H]+ m/zs in GNPS and UNPD")
				unpd_gnps_statistics['Value'].append(
					f"{(clean_data.shape[0] - number_protonated_identified_unpd_best - number_protonated_identified_gnps_best)} ({(clean_data.shape[0] - number_protonated_identified_gnps) / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of not blank [M+H]+ m/zs that were NOT identified against the UNPD or the GNPS libraries (curated_identification_best_origin == ''), thus, may represent putative novel compounds not present in the databases. And its percentage over the total number of not blank [M+H]+ m/zs.")
			# dataset coverage
			number_protonated_unique_identification_gnps = clean_data.loc[
				clean_data.gnps_category != "out", "gnps_Smiles"].unique().size
			gnps_statistics['Statistics'].append(
				"Number of [M+H]+ m/zs unique identifications and GNPS coverage")
			gnps_statistics['Value'].append(
				f"{number_protonated_unique_identification_gnps} ({number_protonated_unique_identification_gnps / total_gnps_fixo * 100:.1f}%)")
			gnps_statistics['Description'].append(
				"The number of unique and curated identifications in GNPS for the not blank [M+H]+ m/zs and its percentage over the total number of entries in GNPS (unique gnps_Smiles for gnps_category != 'out' over 600000 compounds for GNPS total) - the GNPS coverage by [M+H]+.")
			if best_origin_exists:
				number_protonated_unique_identification_gnps_best = clean_data.best_origin_SMILES[
					(clean_data.curated_identification_best_origin == "GNPS")].unique().size
				number_protonated_unique_identification_unpd_best = clean_data.best_origin_SMILES[
					(clean_data.curated_identification_best_origin == "UNPD")].unique().size
				unpd_gnps_statistics['Statistics'].append(
					"Number of [M+H]+ m/zs unique best identifications and GNPS coverage")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_unique_identification_gnps_best} ({number_protonated_unique_identification_gnps_best / total_gnps_fixo * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"The number of unique and curated identifications in GNPS selected as the best origin for the not blank [M+H]+ m/zs and its percentage over the total number of entries in GNPS (unique best_origin_SMILES for curated_identification_best_origin == 'GNPS' over 600000 compounds for GNPS total) - the best GNPS coverage by [M+H]+.")
				total_unpd_fixo = 170602  # total unique UNPD
				unpd_gnps_statistics['Statistics'].append(
					"Number of [M+H]+ m/zs unique best identifications and UNPD coverage")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_unique_identification_unpd_best} ({number_protonated_unique_identification_unpd_best / total_unpd_fixo * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"The number of unique and curated identifications in UNPD selected as the best origin for the not blank [M+H]+ m/zs and its percentage over the total number of entries in UNPD (unique best_origin_SMILES for curated_identification_best_origin == 'UNPD' over 170602 compounds for UNPD total) - the best UNPD coverage by [M+H]+.")
		# chemistry diversity for NP based on superclass annotation, from the annotated and curated data
		if 'gnps_curated_superclass' in clean_data.columns:
			total_superclass_npclassifier = 91  # from the superclass_groupings dictionary in tremolo_UNPD_curate_identification
			# use curated superclass, remove NA and None
			number_unique_superclass = np.unique([x for x in clean_data.gnps_curated_superclass.values[clean_data.gnps_category != "out"]
			                                      if x == x and x is not None and x is not ""]).size
			gnps_statistics['Statistics'].append("Chemical diversity of [M+H]+ in GNPS Superclasses")
			gnps_statistics['Value'].append(
				f"{number_unique_superclass} ({number_unique_superclass / total_superclass_npclassifier * 100:.1f}%)")
			gnps_statistics['Description'].append(
				"The unique number of superclasses that got identified by the not blank [M+H]+ m/zs in GNPS. And its percentage over the total number of unique superclasses considered ("+str(total_superclass_npclassifier)+" for GNPS using NPClassifier). ")
			if best_origin_exists and "best_origin_curated_superclass" in clean_data.columns:
				number_unique_superclass_best = np.unique([x for x in chain.from_iterable(
					clean_data.best_origin_curated_superclass[clean_data.curated_identification_best_origin != ""].str.split(":", expand=True).to_numpy())
				                                           if x==x and x is not None and x is not ""]).size
				unpd_gnps_statistics['Statistics'].append("Chemical diversity of [M+H]+ in GNPS and UNPD Superclasses best")
				unpd_gnps_statistics['Value'].append(
					f"{number_unique_superclass_best} ({number_unique_superclass_best / total_superclass_npclassifier * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"The unique number of superclasses that got identified by the not blank [M+H]+ m/zs in GNPS or UNPD for the best origin. And its percentage over the total number of unique superclasses considered (" + str(
						total_superclass_npclassifier) + " for NPClassifier). ")
		if 'gnps_curated_superclass_grouping' in clean_data.columns:
			total_superclass_npclassifier_grouping = 10
			number_unique_superclass_grouping = np.unique(clean_data.gnps_curated_superclass_grouping.values[(clean_data.gnps_category != "out") & (clean_data.gnps_curated_superclass_grouping != "Not_Annotated")]).size
			gnps_statistics['Statistics'].append("Chemical diversity of [M+H]+ in GNPS Superclasses grouping")
			gnps_statistics['Value'].append(
				f"{number_unique_superclass_grouping} ({number_unique_superclass_grouping / total_superclass_npclassifier_grouping * 100:.1f}%)")
			gnps_statistics['Description'].append(
				"The unique number of superclasses grouping that got identified by the not blank [M+H]+ m/zs in GNPS. And its percentage over the total number of unique superclasses grouping considered (" + str(
					total_superclass_npclassifier_grouping) + " groups proposed by NP3 from the NPClassifier superclasses without the not annotated ones). ")
			if best_origin_exists and "best_origin_curated_superclass_grouping" in clean_data.columns:
				number_unique_superclass_grouping_best = np.unique(
					clean_data.best_origin_curated_superclass_grouping.values[(clean_data.curated_identification_best_origin != "") & (clean_data.best_origin_curated_superclass_grouping != "Not_Annotated")]).size
				unpd_gnps_statistics['Statistics'].append("Chemical diversity of [M+H]+ in GNPS and UNPD Superclasses best grouping")
				unpd_gnps_statistics['Value'].append(
					f"{number_unique_superclass_grouping_best} ({number_unique_superclass_grouping_best / total_superclass_npclassifier_grouping * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"The unique number of superclasses grouping that got identified by the not blank [M+H]+ m/zs in UNPD or GNPS libraries. And its percentage over the total number of unique superclasses grouping considered (" + str(
						total_superclass_npclassifier_grouping) + " groups proposed by NP3 from the NPClassifier superclasses without the not annotated ones). ")
	# save results
	# save the gnps statistics
	gnps_statistics_table = pd.DataFrame(gnps_statistics)
	gnps_statistics_table.to_csv(output_path / "chemical_identification_statistics_GNPS.csv", index=False)
	# save the unpd and gnps best statistics
	if best_origin_exists:
		unpd_gnps_statistics_table = pd.DataFrame(unpd_gnps_statistics)
		unpd_gnps_statistics_table.to_csv(output_path / "chemical_identification_statistics_UNPDxGNPS_best_origin.csv", index=False)

