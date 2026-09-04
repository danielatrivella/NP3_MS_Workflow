#!/usr/bin/python

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from itertools import chain

# Get the path of the src folder (two levels up from this file)
project_root = str(Path(__file__).resolve().parents[1])
# Add the src folder to sys.path if it isn't already there
if project_root not in sys.path:
	sys.path.insert(0, project_root)
# import the post_process_analysis func with relative path
from post_process_analysis.drug_discovery_post_analysis import plot_stacked_bar_pandas_df, superclass_colors

total_unpd_unique_SMILES = 183962  # total unique SMILES in UNPD
total_gnps_unique_SMILES = 90253 # total unique SMILES in GNPS from April 2026
gnps_update_date = "April 2026"

# clean_table_file must contain the clean count table with peak area quantification
def plot_superclass_samples_distribution(metadata_file, clean_table_file, output_path,
                                         superclass_grouping_name="tremolo_curated_superclass_grouping"):
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
	
	# fix the metadata column to upper and the sample types to lower
	metadata.columns = metadata.columns.str.upper()
	metadata["SAMPLE_TYPE"] = metadata.SAMPLE_TYPE.str.lower()
	# if there is a curated identification and any not blank sample, proceed for plotting
	if superclass_grouping_name in clean_data.columns and (~metadata.SAMPLE_TYPE.isin(["blank","bed"])).any():
		print("  - Creating the samples composition by superclass grouping distribution of not blank or bed m/z \n")
		# remove blank and beds if any:
		if "BLANKS_TOTAL" in clean_data.columns:
			clean_data = clean_data.loc[clean_data.BLANKS_TOTAL == 0, :]
		if "BEDS_TOTAL" in clean_data.columns:
			clean_data = clean_data.loc[clean_data.BEDS_TOTAL == 0, :]
		
		# get the columns names containing the count of spectra by peak area without blanks or bed
		samples_area_name = metadata.SAMPLE_CODE[~metadata.SAMPLE_TYPE.isin(["blank","bed"])].values
		if clean_table_file.name.find("peak_area") > 0:
			samples_area_col = samples_area_name + "_area"
		else:
			samples_area_col = samples_area_name + "_spectra"
			
		# create for all not protonated mzs then only for protonated
		for protonated in [0, 1]:
			if protonated == 1:
				if not ("protonated_representative" in clean_data.columns):
					print(
						"    - No protonated representative column present, skipping the protonated distribution plot.")
					break
				output_name_plot = superclass_grouping_name + "_protonated"
				# filter only the protonated nodes
				clean_data = clean_data.loc[clean_data.protonated_representative == 1, :]
				if clean_data.shape[0] == 0:
					print(
						"    - No protonated m/z present, skipping the protonated distribution plot.")
					break
				mzs_selected = "[M+H]+ "
			else:
				output_name_plot = superclass_grouping_name
				mzs_selected = ""
		
			# this is equal to the DD post analysis, reuse the plotting function and make it by presence of the mz (peak area > 0)
			# check the m/z that are present in each sample, which have a peak area > 0
			present_mzs_by_sample = clean_data[samples_area_col].fillna(0) > 0
			# group the quantification columns by the superclass grouping and sum the respective rows
			samples_area_by_superclass_grouping = present_mzs_by_sample.groupby(clean_data[superclass_grouping_name]).sum().T
			# and order the superclass grouping columns with the colors order, fill any missing class with 0
			superclass_grouping_columns = list(superclass_colors.keys())
			missing_groups = samples_area_by_superclass_grouping.columns.symmetric_difference(
				superclass_grouping_columns).values
			if len(missing_groups) > 0:
				samples_area_by_superclass_grouping[missing_groups] = 0
			samples_area_by_superclass_grouping = samples_area_by_superclass_grouping.loc[:, superclass_grouping_columns]
			# normalize the quantification by sample
			samples_area_by_superclass_grouping = samples_area_by_superclass_grouping.div(samples_area_by_superclass_grouping.sum(axis=1), axis=0)
			# rename rows with the original samples codes
			samples_area_by_superclass_grouping.index = samples_area_name
			
			# call the stacked bar plot and pass the superclass fixed colors
			plot_stacked_bar_pandas_df(
				output_path / ("samples_composition_"+output_name_plot+"_distribution.png"),
				samples_area_by_superclass_grouping,
				title='Samples Composition by Superclass Grouping Distribution (normalized by sample)',
				xlabel='Samples', ylabel='Percentage of detected ' + mzs_selected + 'm/z',
				colors=list(superclass_colors.values()),
				figsize=(20, 10), label_size=16,
				title_size=20, legend_title="Superclass Grouping",
				legend_bbox_to_anchor=(0.5, -0.35),
				legend_fontsize=15,
				legend_ncol=4)


def compute_chemical_report_statistics(clean_table_file, output_path):
	total_superclass_unpd = 94  # from the superclass_groupings dictionary in tremolo_UNPD_curate_identification
	total_superclass_npclassifier_grouping = 10
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
	number_valid_mzs = clean_data.shape[0]
	chemical_statistics['Statistics'].append("Number of not blank m/zs")
	chemical_statistics['Value'].append(f"{number_valid_mzs} ({number_valid_mzs/n*100:.1f}%)")
	chemical_statistics['Description'].append("Total number of m/zs in the final table that do not appear in any blank sample (BLANKS_TOTAL == 0). And its percentage over the total number of m/zs.")
	
	if number_valid_mzs > 0:
		# if any bed sample, compute their statistics
		if "BEDS_TOTAL" in clean_data.columns:
			number_not_bed_mzs = sum((clean_data.BEDS_TOTAL == 0))
			chemical_statistics['Statistics'].append("Number of not blank and not bed m/zs ")
			chemical_statistics['Value'].append(
				f"{number_not_bed_mzs} ({number_not_bed_mzs / number_valid_mzs * 100:.1f}%)")
			chemical_statistics['Description'].append(
				"Total number of not bed (not culture media) m/zs over the not blank m/zs, BEDS_TOTAL == 0.")
			#number_bed_mzs = sum((clean_data.BEDS_TOTAL > 0))
			#chemical_statistics['Statistics'].append("Total number of bed m/zs not blank")
			#chemical_statistics['Value'].append(
			#	f"{number_bed_mzs} ({number_bed_mzs / number_not_blank_mzs * 100:.1f}%)")
			#chemical_statistics['Description'].append(
			#	"Total number of bed (culture media) m/zs in the not blank m/zs, BEDS_TOTAL > 0.")
			# skip one empty row
			chemical_statistics['Statistics'].append("")
			chemical_statistics['Value'].append("")
			chemical_statistics['Description'].append("")
			# compute number of putative molecules - protonated m/zs [M+H]+ without blanks only - SSMN [M+H]+ #nodes
			if "protonated_representative" in clean_data.columns:
				number_protonated_mzs = clean_data.protonated_representative.sum()
			else:
				number_protonated_mzs = 0
			chemical_statistics['Statistics'].append("Number of not blank [M+H]+ m/zs")
			chemical_statistics['Value'].append(
				f"{number_protonated_mzs} ({number_protonated_mzs / number_valid_mzs * 100:.1f}%)")
			chemical_statistics['Description'].append(
				"Total number of not blank m/zs that were selected as [M+H]+ - protonated (protonated_representative == 1) - putative molecules. And its percentage over the total number of blanks m/zs, beds are included.")
			
			# filter only the not bed mzs
			mz_types = "not blank and not bed"
			clean_data = clean_data.loc[clean_data.BEDS_TOTAL == 0, :]
			number_valid_mzs = clean_data.shape[0]
		else:
			# skip one empty row
			chemical_statistics['Statistics'].append("")
			chemical_statistics['Value'].append("")
			chemical_statistics['Description'].append("")
			mz_types = "not blank"
			
		# compute number of putative molecules - protonated m/zs [M+H]+
		if "protonated_representative" in clean_data.columns:
			number_protonated_mzs = clean_data.protonated_representative.sum()
		else:
			number_protonated_mzs = 0
		chemical_statistics['Statistics'].append("Number of "+mz_types+" [M+H]+ m/zs")
		chemical_statistics['Value'].append(f"{number_protonated_mzs} ({number_protonated_mzs / number_valid_mzs*100:.1f}%)")
		chemical_statistics['Description'].append(
			"Total number of "+mz_types+" m/zs that were selected as [M+H]+ - protonated (protonated_representative == 1) - putative molecules. And its percentage over the total number of "+mz_types+" m/zs.")
		
		if number_protonated_mzs > 0:
			# filter only [M+H]+
			clean_data = clean_data.loc[clean_data.protonated_representative == 1,:]
			
			# UNPD identification statistics for putative molecules
			# # number of "+mz_types+" [M+H]+ m/zs identified / UNPD size; potential for chemical novelty of the set (percentage of [M+H]+ not identified
			if 'tremolo_UNPD_category_best' in clean_data.columns:
				# skip one empty row
				chemical_statistics['Statistics'].append("")
				chemical_statistics['Value'].append("")
				chemical_statistics['Description'].append("")
				#
				total_unpd_fixo = total_unpd_unique_SMILES  # total unique Smiles in UNPD
				number_protonated_identified_unpd = (clean_data.tremolo_UNPD_category_best != "out").sum()
				chemical_statistics['Statistics'].append("Number of [M+H]+ m/zs identified in UNPD curated and spectral identification rate")
				chemical_statistics['Value'].append(
					f"{number_protonated_identified_unpd} ({number_protonated_identified_unpd / number_protonated_mzs*100:.1f}%)")
				chemical_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were identified against the UNPD using tremolo and passed the curation (tremolo_UNPD_category_best != 'out'). And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
				chemical_statistics['Statistics'].append("Number of putative novel [M+H]+ m/zs curated")
				chemical_statistics['Value'].append(
					f"{(number_protonated_mzs-number_protonated_identified_unpd)} ({(number_protonated_mzs-number_protonated_identified_unpd) / number_protonated_mzs * 100:.1f}%)")
				chemical_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were NOT identified against the UNPD using tremolo after the curation (tremolo_UNPD_category_best == 'out'), thus, may represent putative novel compounds not present in the database. And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
				number_protonated_unique_identification_unpd = clean_data.loc[clean_data.tremolo_UNPD_category_best != "out","tremolo_SMILES_best"].unique().size
				chemical_statistics['Statistics'].append("Number of unique curated library annotations of the [M+H]+ m/zs in UNPD and their coverage")
				chemical_statistics['Value'].append(
					f"{number_protonated_unique_identification_unpd} ({number_protonated_unique_identification_unpd / total_unpd_fixo*100:.1f}%)")
				chemical_statistics['Description'].append(
					"The number of unique and curated library annotations in UNPD for the "+mz_types+" [M+H]+ m/zs and its percentage over the total number of unique SMILES in UNPD (unique tremolo_SMILES_best for tremolo_UNPD_category_best != 'out' over "+str(total_unpd_fixo)+" compounds for UNPD total unique SMILES) - the UNPD coverage by [M+H]+.")
				# chemistry diversity for NP based on superclass annotation, from the annotated and curated data
				# use curated superclass, remove NA and None
				number_unique_superclass = np.unique([x for x in chain.from_iterable(clean_data.tremolo_curated_superclass.str.split(":", expand=True).to_numpy())
				                                      if x==x and x is not None and x is not ""]).size
				chemical_statistics['Statistics'].append("Chemical diversity of [M+H]+ m/zs in UNPD Superclasses curated")
				chemical_statistics['Value'].append(
					f"{number_unique_superclass} ({number_unique_superclass / total_superclass_unpd*100:.1f}%)")
				chemical_statistics['Description'].append(
					"The unique number of curated superclasses that got identified by the "+mz_types+" [M+H]+ m/zs in UNPD (column tremolo_curated_superclass). And its percentage over the total number of unique superclasses considered ("+str(total_superclass_unpd)+" for UNPD using NPClassifier). ")
				number_unique_superclass_grouping = np.unique(clean_data.tremolo_curated_superclass_grouping[(clean_data.tremolo_curated_superclass_grouping != "Not_Annotated")].values).size
				chemical_statistics['Statistics'].append("Chemical diversity of [M+H]+ m/zs in UNPD Superclasses curated grouping")
				chemical_statistics['Value'].append(
					f"{number_unique_superclass_grouping} ({number_unique_superclass_grouping / total_superclass_npclassifier_grouping * 100:.1f}%)")
				chemical_statistics['Description'].append(
					"The unique number of curated superclasses grouping that got identified by the "+mz_types+" [M+H]+ m/zs in UNPD (column tremolo_curated_superclass_grouping). And its percentage over the total number of unique superclasses grouping considered (" + str(
						total_superclass_npclassifier_grouping) + " groups proposed by NP3 from the NPClassifier superclasses). ")
			# chemical identification statistics before curation using the best tremolo result
			if 'tremolo_best_position' in clean_data.columns:
				# skip one empty row
				chemical_statistics['Statistics'].append("")
				chemical_statistics['Value'].append("")
				chemical_statistics['Description'].append("")
				# add best library annotations stats
				number_protonated_identified_unpd = (~clean_data.tremolo_best_position.isna()).sum()
				chemical_statistics['Statistics'].append(
					"Number of [M+H]+ m/zs identified in UNPD best and spectral identification rate")
				chemical_statistics['Value'].append(
					f"{number_protonated_identified_unpd} ({number_protonated_identified_unpd / number_protonated_mzs * 100:.1f}%)")
				chemical_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were identified against the UNPD using tremolo and selected as best result (tremolo_best_position != NA). And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
				chemical_statistics['Statistics'].append("Number of putative novel [M+H]+ m/zs best")
				chemical_statistics['Value'].append(
					f"{(number_protonated_mzs - number_protonated_identified_unpd)} ({(number_protonated_mzs - number_protonated_identified_unpd) / number_protonated_mzs * 100:.1f}%)")
				chemical_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were NOT identified against the UNPD using tremolo best result (tremolo_best_position == NA), thus, may represent putative novel compounds not present in the database. And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
				number_protonated_unique_identification_unpd = clean_data.loc[
					(~clean_data.tremolo_best_position.isna()), "tremolo_SMILES_best"].unique().size
				chemical_statistics['Statistics'].append(
					"Number of unique best library annotations of the [M+H]+ m/zs in UNPD and their coverage")
				chemical_statistics['Value'].append(
					f"{number_protonated_unique_identification_unpd} ({number_protonated_unique_identification_unpd / total_unpd_fixo * 100:.1f}%)")
				chemical_statistics['Description'].append(
					"The number of unique best library annotations in UNPD for the "+mz_types+" [M+H]+ m/zs and its percentage over the total number of unique SMILES in UNPD (unique tremolo_SMILES_best for tremolo_best_position != NA over " + str(
						total_unpd_fixo) + " compounds for UNPD total unique SMILES) - the UNPD coverage by [M+H]+.")
				# use superclass clean, remove NA and None
				number_unique_superclass = np.unique([x for x in chain.from_iterable(
					clean_data.tremolo_NPClassifier_superclass_clean_best.str.split(":", expand=True).to_numpy())
				                                      if x == x and x is not None and x is not ""]).size
				chemical_statistics['Statistics'].append(
					"Chemical diversity of [M+H]+ m/zs in UNPD Superclasses best")
				chemical_statistics['Value'].append(
					f"{number_unique_superclass} ({number_unique_superclass / total_superclass_unpd * 100:.1f}%)")
				chemical_statistics['Description'].append(
					"The unique number of best superclasses that got identified by the "+mz_types+" [M+H]+ m/zs in UNPD (column tremolo_NPClassifier_superclass_clean_best not NA or empty). And its percentage over the total number of unique superclasses considered (" + str(
						total_superclass_unpd) + " for UNPD using NPClassifier). ")
				number_unique_superclass_grouping = np.unique(clean_data.tremolo_NPClassifier_superclass_grouping_best[(
							clean_data.tremolo_NPClassifier_superclass_grouping_best != "Not_Annotated")].values).size
				chemical_statistics['Statistics'].append(
					"Chemical diversity of [M+H]+ m/zs in UNPD Superclasses best grouping")
				chemical_statistics['Value'].append(
					f"{number_unique_superclass_grouping} ({number_unique_superclass_grouping / total_superclass_npclassifier_grouping * 100:.1f}%)")
				chemical_statistics['Description'].append(
					"The unique number of best superclasses grouping that got identified by the "+mz_types+" [M+H]+ m/zs in UNPD (column tremolo_NPClassifier_superclass_grouping_best not NA or empty). And its percentage over the total number of unique superclasses grouping considered (" + str(
						total_superclass_npclassifier_grouping) + " groups proposed by NP3 from the NPClassifier superclasses). ")
			
			# skip one empty row
			chemical_statistics['Statistics'].append("")
			chemical_statistics['Value'].append("")
			chemical_statistics['Description'].append("")
			# samples redundancy (number of [M+H]+ that appear in more than one sample, mediana+-sd of number of samples that they appear
			if clean_table_file.name.find("peak_area") > 0:
				clean_data["number_samples"] = (clean_data.loc[:, clean_data.columns.str.endswith("_area")] > 0).sum(1)
			else:
				clean_data["number_samples"] = (clean_data.loc[:, clean_data.columns.str.endswith("_spectra")] > 0).sum(1)
			number_protonated_in_redundant_samples = (clean_data.number_samples > 1).sum()
			number_samples_describe = clean_data["number_samples"].describe()
			
			chemical_statistics['Statistics'].append("Number of [M+H]+ m/zs that are redundant among the samples")
			chemical_statistics['Value'].append(
				f"{number_protonated_in_redundant_samples} ({number_protonated_in_redundant_samples / number_protonated_mzs*100:.1f}%)")
			chemical_statistics['Description'].append(
				"Number of "+mz_types+" [M+H]+ m/zs that appear in more than one sample (<sample_code>_area > 0 for at least two samples). And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
			
			chemical_statistics['Statistics'].append("Median number of samples that the [M+H]+ m/zs appear")
			chemical_statistics['Value'].append(
				f"{number_samples_describe['50%']:.1f} +-({number_samples_describe['std']:.1f})")
			chemical_statistics['Description'].append(
				"Median number of samples that the "+mz_types+" [M+H]+ m/zs appear +- its standard deviation.")
			
			chemical_statistics['Statistics'].append("Mean number of samples that the [M+H]+ m/zs appear")
			chemical_statistics['Value'].append(
				f"{number_samples_describe['mean']:.1f} +-({number_samples_describe['std']:.1f})")
			chemical_statistics['Description'].append(
				"Mean number of samples that the "+mz_types+" [M+H]+ m/zs appear +- its standard deviation.")
			# samples redundancy by number of samples in which they appear (number of [M+H]+ that appear in different count of samples)
			for count_samples_redundancy in np.sort(clean_data.number_samples.unique()):
				number_protonated_in_redundant_samples_count = (clean_data.number_samples == count_samples_redundancy).sum()
				chemical_statistics['Statistics'].append("Number of [M+H]+ m/zs that appear in "+str(count_samples_redundancy)+" samples")
				chemical_statistics['Value'].append(
					f"{number_protonated_in_redundant_samples_count} ({number_protonated_in_redundant_samples_count / number_protonated_mzs * 100:.1f}%)")
				chemical_statistics['Description'].append(
					"Number of "+mz_types+" [M+H]+ m/zs that appear in exactly "+str(count_samples_redundancy)+" samples and thus are redundant among them (<sample_code>_area > 0 for exactly "+str(count_samples_redundancy)+" samples). And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
	
	chemical_statistics_table = pd.DataFrame(chemical_statistics)
	chemical_statistics_table.to_csv(output_path/"chemical_statistics_UNPD.csv", index=False)
	
def compute_chemical_identification_report_GNPS_result(clean_table_file, output_path):
	total_superclass_npclassifier = 94  # from the superclass_groupings dictionary in tremolo_UNPD_curate_identification
	total_superclass_npclassifier_grouping = 10
	total_gnps_fixo = total_gnps_unique_SMILES  # total unique SMILES in GNPS
	gnps_date = gnps_update_date # updated date from which the total_gnps_unique_SMILES was computed
	total_unpd_fixo = total_unpd_unique_SMILES  # total unique Smiles in UNPD
	clean_table_file = Path(clean_table_file)
	output_path = Path(output_path)
	
	if not clean_table_file.exists() or not clean_table_file.is_file():
		sys.exit("The provided path to the clean data file does not exists. Chemical and identification report for GNPS result statistics aborted.")
	output_path.mkdir(exist_ok=True)
	if not output_path.exists() or not output_path.is_dir():
		sys.exit("The provided chemical report output path does not exists. Chemical and identification report for GNPS result statistics aborted.")
	
	print("* Computing the chemical and identification statistics for GNPS and UNPDxGNPS final curation *\n")
	clean_data = pd.read_csv(clean_table_file)
	# if any blank sample, remove blank mzs from this analysis
	if "BLANKS_TOTAL" in clean_data.columns:
		clean_data = clean_data.loc[clean_data.BLANKS_TOTAL == 0, :]
	n = clean_data.shape[0]
	gnps_curated = True
	if 'gnps_category' not in clean_data.columns or (clean_data.gnps_category == "out").all():
		gnps_curated = False
		print("  - No valid GNPS curated data is available. All gnps_category is 'out' or missing. The chemical identification report for GNPS will use the best result not curated.")
	if 'gnps_Smiles' in clean_data.columns and clean_data.gnps_Smiles.isna().all():
		print("  - No valid GNPS identification data is available. All gnps_Smiles are NA. Aborting the chemical identification report for GNPS and final curated library annotation.")
		return 0
	
	# create dictionary to store the chemical and identification statistics of the job for GNPS result
	gnps_statistics = {'Statistics': [],
	                   'Value': [],
	                   'Description': []}
	
	gnps_statistics['Statistics'].append("Total number of not blank m/zs")
	gnps_statistics['Value'].append(str(n))
	gnps_statistics['Description'].append("Total number of not blank m/zs in the final table (BLANKS_TOTAL == 0).")
	# also create statistics for the UNPDxGNPS - if final curated library annotation exist
	# create dictionary to store the chemical and identification statistics of the job for UNPDxGNPS result
	if "curated_lib_annotation_SMILES" in clean_data.columns:
		curated_lib_annotation_exists = True
		unpd_gnps_statistics = {'Statistics': [],
		                        'Value': [],
		                        'Description': []}
		
		unpd_gnps_statistics['Statistics'].append("Total number of not blank m/zs")
		unpd_gnps_statistics['Value'].append(str(n))
		unpd_gnps_statistics['Description'].append("Total number of not blank m/zs in the final table (BLANKS_TOTAL == 0).")
	else:
		print("  - The curated library annotation data is not present (curated_lib_annotation_SMILES column is missing), its report will be skipped.")
		curated_lib_annotation_exists = False
	# also remove BEDs if present, add in the description a formated text with this info
	# if any bed sample, compute their statistics
	if "BEDS_TOTAL" in clean_data.columns:
		number_not_bed_mzs = sum((clean_data.BEDS_TOTAL == 0))
		gnps_statistics['Statistics'].append("Total number of not blank and not bed m/zs ")
		gnps_statistics['Value'].append(
			f"{number_not_bed_mzs} ({number_not_bed_mzs / n * 100:.1f}%)")
		gnps_statistics['Description'].append(
			"Total number of not bed (not culture media) m/zs over the not blank m/zs, BEDS_TOTAL == 0.")
		if curated_lib_annotation_exists:
			unpd_gnps_statistics['Statistics'].append("Total number of not blank and not bed m/zs ")
			unpd_gnps_statistics['Value'].append(
				f"{number_not_bed_mzs} ({number_not_bed_mzs / n * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of not bed (not culture media) m/zs over the not blank m/zs, BEDS_TOTAL == 0.")
		# filter only the not bed mzs
		mz_types = "not blank and not bed"
		clean_data = clean_data.loc[clean_data.BEDS_TOTAL == 0, :]
		n = clean_data.shape[0]
	else:
		mz_types = "not blank"
	if n > 0:
		if gnps_curated:
			# skip one empty row
			gnps_statistics['Statistics'].append("")
			gnps_statistics['Value'].append("")
			gnps_statistics['Description'].append("")
			# identification statistics for spectra identification rate
			number_identified_gnps = (clean_data.gnps_category != "out").sum()
			gnps_statistics['Statistics'].append("Number of m/zs identified in GNPS curated and spectral identification rate")
			gnps_statistics['Value'].append(
				f"{number_identified_gnps} ({number_identified_gnps / n * 100:.1f}%)")
			gnps_statistics['Description'].append(
				"Total number of "+mz_types+" m/zs that were identified against the GNPS libraries (gnps_category != 'out'). And its percentage over the total number of "+mz_types+" m/zs - the spectral identification rate.")
			# unique curated library annotations
			number_unique_gnps_curated = clean_data.gnps_Smiles[
				(clean_data.gnps_category != "out")].unique().size
			gnps_statistics['Statistics'].append("Number of unique GNPS library annotations curated and their coverage")
			gnps_statistics['Value'].append(
				f"{number_unique_gnps_curated} ({number_unique_gnps_curated / total_gnps_fixo * 100:.1f}%)")
			gnps_statistics['Description'].append(
				"Total number of unique GNPS library annotations from "+mz_types+" m/zs that passed the curation (unique gnps_Smiles with gnps_category != 'out') - unique identified molecules and its percentage over the total number of unique SMILES in GNPS from "+gnps_date+" (unique gnps_Smiles for gnps_category != 'out' over " + str(
					total_gnps_fixo) + " compounds for GNPS total unique SMILES) - the GNPS coverage.")
		if curated_lib_annotation_exists:
			# skip one empty row
			unpd_gnps_statistics['Statistics'].append("")
			unpd_gnps_statistics['Value'].append("")
			unpd_gnps_statistics['Description'].append("")
			# identification statistics for spectra identification rate
			number_identified_gnps = (clean_data.curated_lib_annotation_origin == "GNPS").sum()
			number_identified_unpd = (clean_data.curated_lib_annotation_origin == "UNPD").sum()
			unpd_gnps_statistics['Statistics'].append("Number of m/zs identified in UNPD or GNPS curated and spectral identification rate")
			unpd_gnps_statistics['Value'].append(
				f"{number_identified_gnps+number_identified_unpd} ({(number_identified_gnps+number_identified_unpd) / n * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of "+mz_types+" m/zs that were identified against the UNPD and GNPS libraries - final curated library annotation (curated_lib_annotation_origin != ''). And its percentage over the total number of "+mz_types+" m/zs - the spectral identification rate.")
			unpd_gnps_statistics['Statistics'].append(
				"Number of m/zs identified in UNPD as final curated library annotation curated and spectral identification rate")
			unpd_gnps_statistics['Value'].append(
				f"{number_identified_unpd} ({(number_identified_unpd) / n * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of "+mz_types+" m/zs that were identified against the UNPD and selected as origin in the final curated library annotation (curated_lib_annotation_origin == 'UNPD'). And its percentage over the total number of "+mz_types+" m/zs - the spectral identification rate.")
			unpd_gnps_statistics['Statistics'].append(
				"Number of m/zs identified in GNPS in the final curated library annotation and spectral identification rate")
			unpd_gnps_statistics['Value'].append(
				f"{number_identified_gnps} ({(number_identified_gnps) / n * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of "+mz_types+" m/zs that were identified against the GNPS and selected as origin in the final curated library annotation (curated_lib_annotation_origin == 'GNPS'). And its percentage over the total number of "+mz_types+" m/zs - the spectral identification rate.")
			# unique curated library annotations for unpd and gnps together
			number_unique_unpd_gnps_curated = clean_data.curated_lib_annotation_SMILES[
				(clean_data.curated_lib_annotation_SMILES != "") & (~clean_data.curated_lib_annotation_SMILES.isna())].unique().size
			unpd_gnps_statistics['Statistics'].append("Number of unique UNPD or GNPS final curated library annotations and their UNPD and GNPS coverage")
			unpd_gnps_statistics['Value'].append(
				f"{number_unique_unpd_gnps_curated} ({number_unique_unpd_gnps_curated/(total_gnps_fixo+total_unpd_fixo)*100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of unique UNPD and GNPS annotations that passed the final library annotation curation (unique curated_lib_annotation_SMILES with curated_lib_annotation_SMILES != '') - unique identified molecules and their percentage over the total number of unique SMILES in GNPS from "+gnps_date+" plus UNPD ("+str(total_gnps_fixo)+"+"+str(total_unpd_fixo)+" compounds for GNPS and UNPD total unique SMILES) - the final origin GNPS and UNPD coverage by [M+H]+.")
			# unique curated library annotations for unpd and gnps separated
			number_unique_unpd_curated = clean_data.curated_lib_annotation_SMILES[
				(clean_data.curated_lib_annotation_origin == "UNPD")].unique().size
			number_unique_gnps_curated = clean_data.curated_lib_annotation_SMILES[
				(clean_data.curated_lib_annotation_origin == "GNPS")].unique().size
			unpd_gnps_statistics['Statistics'].append("Number of unique UNPD origin from the final curated library annotations and its UNPD coverage")
			unpd_gnps_statistics['Value'].append(
				f"{number_unique_unpd_curated} ({number_unique_unpd_curated/total_unpd_fixo*100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of unique UNPD library annotations selected as final origin in the library annotation curation (unique curated_lib_annotation_SMILES with curated_lib_annotation_origin == 'UNPD') - unique identified molecules and its percentage over the total number of unique SMILES in UNPD (" + str(
						total_unpd_fixo) + " compounds for UNPD total unique SMILES) - the UNPD coverage.")
			unpd_gnps_statistics['Statistics'].append("Number of unique GNPS origin from the final curated library annotations and its GNPS coverage")
			unpd_gnps_statistics['Value'].append(
				f"{number_unique_gnps_curated} ({number_unique_gnps_curated / total_gnps_fixo * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of unique GNPS library annotations selected as final origin in the library annotation curation (unique curated_lib_annotation_SMILES with curated_lib_annotation_origin == 'GNPS') - unique identified molecules and its percentage over the total number of unique SMILES in GNPS from "+gnps_date+" (" + str(
				total_gnps_fixo) + " compounds for GNPS total unique SMILES) - the GNPS coverage.")
		# skip one empty row
		gnps_statistics['Statistics'].append("")
		gnps_statistics['Value'].append("")
		gnps_statistics['Description'].append("")
		# all gnps
		# add identification stats for all gnps results not curated
		# identification statistics for spectra identification rate
		number_identified_gnps_all = (~clean_data.gnps_Smiles.isna()).sum()
		gnps_statistics['Statistics'].append(
			"Number of m/zs identified in GNPS all and spectral identification rate")
		gnps_statistics['Value'].append(
			f"{number_identified_gnps_all} ({number_identified_gnps_all / n * 100:.1f}%)")
		gnps_statistics['Description'].append(
			"Total number of all "+mz_types+" m/zs that were identified against the GNPS libraries (gnps_Smiles != NA). And its percentage over the total number of "+mz_types+" m/zs - the spectral identification rate.")
		# unique all library annotations
		number_unique_gnps_all = clean_data.gnps_Smiles[(~clean_data.gnps_Smiles.isna())].unique().size
		gnps_statistics['Statistics'].append("Number of unique GNPS library annotations from all and its GNPS coverage")
		gnps_statistics['Value'].append(
			f"{number_unique_gnps_all} ({number_unique_gnps_all / total_gnps_fixo * 100:.1f}%)")
		gnps_statistics['Description'].append(
			"Total number of unique GNPS library annotations from all "+mz_types+" m/zs (unique gnps_Smiles that are not NA) - unique identified molecules and its percentage over the total number of unique SMILES in GNPS from "+gnps_date+" (" + str(
				total_gnps_fixo) + " compounds for GNPS total unique SMILES) - the GNPS coverage.")
		if curated_lib_annotation_exists:
			# skip one empty row
			unpd_gnps_statistics['Statistics'].append("")
			unpd_gnps_statistics['Value'].append("")
			unpd_gnps_statistics['Description'].append("")
			# identification statistics for all gnps library annotations
			unpd_gnps_statistics['Statistics'].append(
				"Number of m/zs identified in GNPS all and spectral identification rate")
			unpd_gnps_statistics['Value'].append(
				f"{number_identified_gnps_all} ({number_identified_gnps_all / n * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of all "+mz_types+" m/zs that were identified against the GNPS libraries (gnps_Smiles != NA). And its percentage over the total number of "+mz_types+" m/zs - the spectral identification rate.")
			# unique library annotations gnps all
			unpd_gnps_statistics['Statistics'].append("Number of unique GNPS library annotations from all and its GNPS coverage")
			unpd_gnps_statistics['Value'].append(
				f"{number_unique_gnps_all} ({number_unique_gnps_all / total_gnps_fixo * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of unique GNPS library annotations from "+mz_types+" m/zs among all results (unique gnps_Smiles that are not NA) - unique identified molecules and its percentage over the total number of unique SMILES in GNPS from "+gnps_date+" (" + str(
					total_gnps_fixo) + " compounds for GNPS total unique SMILES) - the GNPS coverage.")
			# identification statistics for all unpd library annotations
			number_identified_unpd_all = (~clean_data.tremolo_SMILES_best.isna()).sum()
			unpd_gnps_statistics['Statistics'].append(
				"Number of m/zs identified in UNPD all and spectral identification rate")
			unpd_gnps_statistics['Value'].append(
				f"{number_identified_unpd_all} ({number_identified_unpd_all / n * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of all "+mz_types+" m/zs that were identified against the UNPD using tremolo (tremolo_SMILES_best != NA). And its percentage over the total number of "+mz_types+" m/zs - the spectral identification rate.")
			# unique library annotations unpd all
			number_unique_identification_unpd = clean_data.loc[
				(~clean_data.tremolo_SMILES_best.isna()), "tremolo_SMILES_best"].unique().size
			unpd_gnps_statistics['Statistics'].append("Number of unique UNPD library annotations from all and its UNPD coverage")
			unpd_gnps_statistics['Value'].append(
				f"{number_unique_identification_unpd} ({number_unique_identification_unpd / total_unpd_fixo * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of unique UNPD library annotations from "+mz_types+" m/zs among all results (unique tremolo_SMILES_best that are not NA) - unique identified molecules and its percentage over the total number of unique SMILES in UNPD (" + str(
					total_unpd_fixo) + " compounds for UNPD total unique SMILES) - the UNPD coverage.")
			
		# identification statistics for [M+H]+
		# now compute number of putative molecules - protonated m/zs [M+H]+
		#skip one empty row
		gnps_statistics['Statistics'].append("")
		gnps_statistics['Value'].append("")
		gnps_statistics['Description'].append("")
		
		number_protonated_mzs = clean_data.protonated_representative.sum()
		gnps_statistics['Statistics'].append("Number of [M+H]+ m/zs")
		gnps_statistics['Value'].append(
			f"{number_protonated_mzs} ({number_protonated_mzs / n * 100:.1f}%)")
		gnps_statistics['Description'].append(
			"Total number of "+mz_types+" m/zs that were selected as [M+H]+ (protonated_representative == 1) - putative molecules. And its percentage over the total number of "+mz_types+" m/zs.")
		if curated_lib_annotation_exists:
			# skip one empty row
			unpd_gnps_statistics['Statistics'].append("")
			unpd_gnps_statistics['Value'].append("")
			unpd_gnps_statistics['Description'].append("")
			
			unpd_gnps_statistics['Statistics'].append("Number of [M+H]+ m/zs")
			unpd_gnps_statistics['Value'].append(
				f"{number_protonated_mzs} ({number_protonated_mzs / n * 100:.1f}%)")
			unpd_gnps_statistics['Description'].append(
				"Total number of "+mz_types+" m/zs that were selected as [M+H]+ (protonated_representative == 1) - putative molecules. And its percentage over the total number of "+mz_types+" m/zs.")
		# if any protonated, proceed
		if number_protonated_mzs > 0:
			# filter only [M+H]+
			clean_data = clean_data.loc[clean_data.protonated_representative == 1, :]
			if gnps_curated:
				# skip one empty row
				gnps_statistics['Statistics'].append("")
				gnps_statistics['Value'].append("")
				gnps_statistics['Description'].append("")
				# GNPS identification statistics for putative molecules
				# # number of not blank or bed [M+H]+ m/zs identified / GNPS size; potential for chemical novelty of the set (percentage of [M+H]+ not identified
				number_protonated_identified_gnps = (clean_data.gnps_category != "out").sum()
				gnps_statistics['Statistics'].append("Number of [M+H]+ m/zs identified in GNPS curated and spectral identification rate")
				gnps_statistics['Value'].append(
					f"{number_protonated_identified_gnps} ({number_protonated_identified_gnps / number_protonated_mzs * 100:.1f}%)")
				gnps_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were identified against the GNPS and passed the curation (gnps_category != 'out'). And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
				# not identified m/zs - novelty
				gnps_statistics['Statistics'].append("Number of putative novel [M+H]+ m/zs in GNPS curated")
				gnps_statistics['Value'].append(
					f"{(number_protonated_mzs - number_protonated_identified_gnps)} ({(number_protonated_mzs - number_protonated_identified_gnps) / number_protonated_mzs * 100:.1f}%)")
				gnps_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were NOT identified against the GNPS or did not pass the curation (gnps_category == 'out'), thus, may represent putative novel compounds not present in the database. And its percentage over the total number of "+mz_types+" [M+H]+ m/zs - the spectral identification rate.")
				# dataset coverage
				number_protonated_unique_identification_gnps = clean_data.loc[
					((clean_data.gnps_category != "out") & (~clean_data.gnps_Smiles.isna())), "gnps_Smiles"].unique().size
				gnps_statistics['Statistics'].append(
					"Number of unique curated library annotations in GNPS of the [M+H]+ m/zs and their GNPS coverage")
				gnps_statistics['Value'].append(
					f"{number_protonated_unique_identification_gnps} ({number_protonated_unique_identification_gnps / total_gnps_fixo * 100:.1f}%)")
				gnps_statistics['Description'].append(
					"The number of unique and curated library annotations in GNPS for the "+mz_types+" [M+H]+ m/zs and its percentage over the total number of unique SMILES in GNPS from "+gnps_date+" (unique gnps_Smiles not NA for gnps_category != 'out' over " + str(
						total_gnps_fixo) + " compounds for GNPS total unique SMILES) - the GNPS coverage by [M+H]+.")
				if 'gnps_curated_superclass' in clean_data.columns:
					# use curated superclass, remove NA and None
					number_unique_superclass = np.unique(
						[x for x in clean_data.gnps_curated_superclass.values[clean_data.gnps_category != "out"]
						 if x == x and x is not None and x is not ""]).size
					gnps_statistics['Statistics'].append(
						"Chemical diversity of [M+H]+ in GNPS Superclasses curated")
					gnps_statistics['Value'].append(
						f"{number_unique_superclass} ({number_unique_superclass / total_superclass_npclassifier * 100:.1f}%)")
					gnps_statistics['Description'].append(
						"The unique number of superclasses that got identified by the " + mz_types + " [M+H]+ m/zs in GNPS. And its percentage over the total number of unique superclasses considered (" + str(
							total_superclass_npclassifier) + " for GNPS using NPClassifier). ")
				if 'gnps_curated_superclass_grouping' in clean_data.columns:
					number_unique_superclass_grouping = np.unique(
						clean_data.gnps_curated_superclass_grouping.values[
							(clean_data.gnps_category != "out") & (
									clean_data.gnps_curated_superclass_grouping != "Not_Annotated")]).size
					gnps_statistics['Statistics'].append(
						"Chemical diversity of [M+H]+ in GNPS Superclasses grouping curated")
					gnps_statistics['Value'].append(
						f"{number_unique_superclass_grouping} ({number_unique_superclass_grouping / total_superclass_npclassifier_grouping * 100:.1f}%)")
					gnps_statistics['Description'].append(
						"The unique number of superclasses grouping that got identified by the " + mz_types + " [M+H]+ m/zs in GNPS. And its percentage over the total number of unique superclasses grouping considered (" + str(
							total_superclass_npclassifier_grouping) + " groups proposed by NP3 from the NPClassifier superclasses without the not annotated ones). ")
			# all gnps M+H
			# add M+H identification stats for all gnps results not curated
			# skip one empty row
			gnps_statistics['Statistics'].append("")
			gnps_statistics['Value'].append("")
			gnps_statistics['Description'].append("")
			# identification statistics for spectra identification rate
			number_protonated_identified_gnps_all = (~clean_data.gnps_SpectrumID.isna()).sum()
			gnps_statistics['Statistics'].append(
				"Number of [M+H]+ m/zs identified in GNPS all and spectral identification rate")
			gnps_statistics['Value'].append(
				f"{number_protonated_identified_gnps_all} ({number_protonated_identified_gnps_all / number_protonated_mzs * 100:.1f}%)")
			gnps_statistics['Description'].append(
				"Total number of all "+mz_types+" [M+H]+ m/zs that were identified against the GNPS libraries (gnps_SpectrumID is not NA). And its percentage over the total number of "+mz_types+" [M+H]+ m/zs - the spectral identification rate.")
			# not identified m/zs - novelty
			gnps_statistics['Statistics'].append("Number of putative novel [M+H]+ m/zs in GNPS all")
			gnps_statistics['Value'].append(
				f"{(number_protonated_mzs - number_protonated_identified_gnps_all)} ({(number_protonated_mzs - number_protonated_identified_gnps_all) / number_protonated_mzs * 100:.1f}%)")
			gnps_statistics['Description'].append(
				"Total number of "+mz_types+" [M+H]+ m/zs that were NOT identified against the GNPS (gnps_SpectrumID is NA), thus, may represent putative novel compounds not present in the database. And its percentage over the total number of "+mz_types+" [M+H]+ m/zs - the spectral identification rate.")
			# unique all library annotations
			number_unique_gnps_all = clean_data.gnps_Smiles[(~clean_data.gnps_Smiles.isna())].unique().size
			gnps_statistics['Statistics'].append("Number of unique GNPS library annotations from all of the [M+H]+ m/zs and GNPS coverage")
			gnps_statistics['Value'].append(
				f"{number_unique_gnps_all} ({number_unique_gnps_all / total_gnps_fixo * 100:.1f}%)")
			gnps_statistics['Description'].append(
				"Total number of unique GNPS library annotations from all "+mz_types+" [M+H]+ m/zs (unique gnps_Smiles that are not NA) - unique identified molecules and its percentage over the total number of unique SMILES in GNPS from "+gnps_date+" (unique gnps_Smiles that are not NA over " + str(
					total_gnps_fixo) + " compounds for GNPS total unique SMILES) - the GNPS coverage.")
			# chemistry diversity for NP based on superclass annotation, from the all annotated and clean data
			if gnps_curated:
				# chemical diversity GNPS all
				# use clean superclass, remove NA and None
				number_unique_superclass = np.unique(
					[x for x in clean_data.gnps_npclassifier_superclass_clean.values[
						~clean_data.gnps_npclassifier_superclass_clean.isna()]
					    if x == x and x is not None and x is not ""]).size
				gnps_statistics['Statistics'].append("Chemical diversity of [M+H]+ in GNPS Superclasses all clean")
				gnps_statistics['Value'].append(
					f"{number_unique_superclass} ({number_unique_superclass / total_superclass_npclassifier * 100:.1f}%)")
				gnps_statistics['Description'].append(
					"The unique number of superclasses that got identified by the "+mz_types+" [M+H]+ m/zs in GNPS. And its percentage over the total number of unique superclasses considered (" + str(
						total_superclass_npclassifier) + " for GNPS using NPClassifier). ")
				# clean super class grouping
				number_unique_superclass_grouping = np.unique(
					clean_data.gnps_npclassifier_superclass_grouping.values[
						(clean_data.gnps_npclassifier_superclass_grouping != "Not_Annotated")]).size
				gnps_statistics['Statistics'].append(
					"Chemical diversity of [M+H]+ in GNPS Superclasses grouping all clean")
				gnps_statistics['Value'].append(
					f"{number_unique_superclass_grouping} ({number_unique_superclass_grouping / total_superclass_npclassifier_grouping * 100:.1f}%)")
				gnps_statistics['Description'].append(
					"The unique number of superclasses grouping that got identified by the "+mz_types+" [M+H]+ m/zs in GNPS. And its percentage over the total number of unique superclasses grouping considered (" + str(
						total_superclass_npclassifier_grouping) + " groups proposed by NP3 from the NPClassifier superclasses without the not annotated ones). ")
			if curated_lib_annotation_exists:
				# skip one empty row
				unpd_gnps_statistics['Statistics'].append("")
				unpd_gnps_statistics['Value'].append("")
				unpd_gnps_statistics['Description'].append("")
				# compute M+H identification stats for the final curated library annotation
				number_protonated_identified_gnps_best = (clean_data.curated_lib_annotation_origin == "GNPS").sum()
				number_protonated_identified_unpd_best = (clean_data.curated_lib_annotation_origin == "UNPD").sum()
				# curated library annotations
				unpd_gnps_statistics['Statistics'].append("Number of [M+H]+ m/zs identified in UNPD or GNPS curated")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_identified_unpd_best+number_protonated_identified_gnps_best} ({(number_protonated_identified_unpd_best+number_protonated_identified_gnps_best) / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were identified against the UNPD or GNPS and passed the final library annotation curation (curated_lib_annotation_origin != ''). And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
				# unpd origin in the final curated library annotation
				unpd_gnps_statistics['Statistics'].append("Number of [M+H]+ m/zs identified in UNPD as final origin in the library annotation curation")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_identified_unpd_best} ({number_protonated_identified_unpd_best / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were identified against the UNPD and selected as final origin in the library annotation curation (curated_lib_annotation_origin == 'UNPD'). And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
				# gnps final origin in the library annotation curation
				unpd_gnps_statistics['Statistics'].append("Number of [M+H]+ m/zs identified in GNPS as final origin in the library annotation curation")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_identified_gnps_best} ({number_protonated_identified_gnps_best / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were identified against the GNPS and selected as final origin in the library annotation curation (curated_lib_annotation_origin == 'GNPS'). And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
				# unpd and gnps novel
				unpd_gnps_statistics['Statistics'].append("Number of putative novel [M+H]+ m/zs in GNPS and UNPD curated")
				unpd_gnps_statistics['Value'].append(
					f"{(number_protonated_mzs - number_protonated_identified_unpd_best - number_protonated_identified_gnps_best)} ({(number_protonated_mzs - number_protonated_identified_unpd_best - number_protonated_identified_gnps_best) / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were NOT identified against the UNPD or the GNPS libraries (curated_lib_annotation_origin == ''), thus, may represent putative novel compounds not present in the databases. And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
				# dataset coverage for final origin in the library annotation curation
				number_protonated_unique_identification_gnps_best = clean_data.curated_lib_annotation_SMILES[
					(clean_data.curated_lib_annotation_origin == "GNPS")].unique().size
				number_protonated_unique_identification_unpd_best = clean_data.curated_lib_annotation_SMILES[
					(clean_data.curated_lib_annotation_origin == "UNPD")].unique().size
				# GNPS coverage
				unpd_gnps_statistics['Statistics'].append(
					"Number of unique final curated library annotations of the [M+H]+ m/zs in GNPS and their coverage")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_unique_identification_gnps_best} ({number_protonated_unique_identification_gnps_best / total_gnps_fixo * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"The number of unique and final curated library annotations in GNPS selected as the final origin for the "+mz_types+" [M+H]+ m/zs and its percentage over the total number of unique SMILES in GNPS from "+gnps_date+" (unique curated_lib_annotation_SMILES for curated_lib_annotation_origin == 'GNPS' over "+str(total_gnps_fixo)+" compounds for GNPS total unique SMILES) - the final curated GNPS coverage by [M+H]+.")
				# UNPD coverage
				unpd_gnps_statistics['Statistics'].append(
					"Number of unique final curated library annotations of the [M+H]+ m/zs in UNPD and their coverage")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_unique_identification_unpd_best} ({number_protonated_unique_identification_unpd_best / total_unpd_fixo * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"The number of unique and final curated library annotations in UNPD selected as the final origin for the "+mz_types+" [M+H]+ m/zs and its percentage over the total number of unique SMILES in UNPD (unique curated_lib_annotation_SMILES for curated_lib_annotation_origin == 'UNPD' over "+str(total_unpd_fixo)+" compounds for UNPD total unique SMILES) - the final curated UNPD coverage by [M+H]+.")
				#
				# chemical diversity of the final curated library annotation
				if "curated_lib_annotation_superclass" in clean_data.columns:
					number_unique_superclass_best = np.unique([x for x in chain.from_iterable(
						clean_data.curated_lib_annotation_superclass[
							(clean_data.curated_lib_annotation_origin != "") &
							(~clean_data.curated_lib_annotation_origin.isna())].str.split(":", expand=True).to_numpy())
					                                           if x == x and x is not None and x is not ""]).size
					unpd_gnps_statistics['Statistics'].append(
						"Chemical diversity of [M+H]+ in GNPS and UNPD Superclasses curated best")
					unpd_gnps_statistics['Value'].append(
						f"{number_unique_superclass_best} ({number_unique_superclass_best / total_superclass_npclassifier * 100:.1f}%)")
					unpd_gnps_statistics['Description'].append(
						"The unique number of superclasses that got identified by the "+mz_types+" [M+H]+ m/zs in GNPS or UNPD after the final curated library annotation. And its percentage over the total number of unique superclasses considered (" + str(
							total_superclass_npclassifier) + " for NPClassifier). ")
				if "curated_lib_annotation_superclass_grouping" in clean_data.columns:
					number_unique_superclass_grouping_best = np.unique(
						clean_data.curated_lib_annotation_superclass_grouping.values[(clean_data.curated_lib_annotation_origin != "") &
						                                                             (~clean_data.curated_lib_annotation_origin.isna()) &
						                                                             (clean_data.curated_lib_annotation_superclass_grouping != "Not_Annotated")]).size
					unpd_gnps_statistics['Statistics'].append("Chemical diversity of [M+H]+ in GNPS and UNPD Superclasses curated best grouping")
					unpd_gnps_statistics['Value'].append(
						f"{number_unique_superclass_grouping_best} ({number_unique_superclass_grouping_best / total_superclass_npclassifier_grouping * 100:.1f}%)")
					unpd_gnps_statistics['Description'].append(
						"The unique number of superclasses grouping that got identified by the "+mz_types+" [M+H]+ m/zs in UNPD or GNPS libraries. And its percentage over the total number of unique superclasses grouping considered (" + str(
							total_superclass_npclassifier_grouping) + " groups proposed by NP3 from the NPClassifier superclasses without the not annotated ones). ")
				# add chemical diversity for all GNPS and UNPD results before curation - no filter
				# add M+H identification stats for all final curated library annotation results not curated - identification no filter
				# skip one empty row
				unpd_gnps_statistics['Statistics'].append("")
				unpd_gnps_statistics['Value'].append("")
				unpd_gnps_statistics['Description'].append("")
				# identification statistics for spectra identification rate GNPS
				unpd_gnps_statistics['Statistics'].append(
					"Number of [M+H]+ m/zs identified in GNPS all and spectral identification rate")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_identified_gnps_all} ({number_protonated_identified_gnps_all / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of all "+mz_types+" [M+H]+ m/zs that were identified against the GNPS libraries (gnps_SpectrumID is not NA). And its percentage over the total number of "+mz_types+" [M+H]+ m/zs - the spectral identification rate.")
				# not identified m/zs - novelty GNPS
				unpd_gnps_statistics['Statistics'].append("Number of putative novel [M+H]+ m/zs in GNPS all")
				unpd_gnps_statistics['Value'].append(
					f"{(number_protonated_mzs - number_protonated_identified_gnps_all)} ({(number_protonated_mzs - number_protonated_identified_gnps_all) / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were NOT identified against the GNPS (gnps_SpectrumID is NA), thus, may represent putative novel compounds not present in the database. And its percentage over the total number of "+mz_types+" [M+H]+ m/zs - the spectral identification rate.")
				# unique all library annotations GNPS
				unpd_gnps_statistics['Statistics'].append(
					"Number of unique GNPS library annotations in all of the [M+H]+ m/zs and GNPS coverage")
				unpd_gnps_statistics['Value'].append(
					f"{number_unique_gnps_all} ({number_unique_gnps_all / total_gnps_fixo * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of unique GNPS library annotations from all "+mz_types+" [M+H]+ m/zs (unique gnps_Smiles that are not NA) - unique identified molecules and its percentage over the total number of unique SMILES in GNPS from "+gnps_date+" (unique gnps_Smiles that are not NA over " + str(
						total_gnps_fixo) + " compounds for GNPS total unique SMILES) - the GNPS coverage.")
				# chemical diversity GNPS all
				if gnps_curated:
					# use clean superclass, remove NA and None
					unpd_gnps_statistics['Statistics'].append("Chemical diversity of [M+H]+ in GNPS Superclasses curated")
					unpd_gnps_statistics['Value'].append(
						f"{number_unique_superclass} ({number_unique_superclass / total_superclass_npclassifier * 100:.1f}%)")
					unpd_gnps_statistics['Description'].append(
						"The unique number of superclasses that got identified by the "+mz_types+" [M+H]+ m/zs in GNPS. And its percentage over the total number of unique superclasses considered (" + str(
							total_superclass_npclassifier) + " for GNPS using NPClassifier). ")
					# clean super class grouping
					unpd_gnps_statistics['Statistics'].append(
						"Chemical diversity of [M+H]+ in GNPS Superclasses curated grouping")
					unpd_gnps_statistics['Value'].append(
						f"{number_unique_superclass_grouping} ({number_unique_superclass_grouping / total_superclass_npclassifier_grouping * 100:.1f}%)")
					unpd_gnps_statistics['Description'].append(
						"The unique number of superclasses grouping that got identified by the "+mz_types+" [M+H]+ m/zs in GNPS. And its percentage over the total number of unique superclasses grouping considered (" + str(
							total_superclass_npclassifier_grouping) + " groups proposed by NP3 from the NPClassifier superclasses without the not annotated ones). ")
				# same stats for UNPD
				# skip one empty row
				unpd_gnps_statistics['Statistics'].append("")
				unpd_gnps_statistics['Value'].append("")
				unpd_gnps_statistics['Description'].append("")
				# add UNPD best library annotations stats
				# identification statistics for spectra identification rate UNPD
				number_protonated_identified_unpd = (~clean_data.tremolo_best_position.isna()).sum()
				unpd_gnps_statistics['Statistics'].append(
					"Number of [M+H]+ m/zs identified in UNPD best and spectral identification rate")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_identified_unpd} ({number_protonated_identified_unpd / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were identified against the UNPD using tremolo and selected as best result (tremolo_best_position != NA). And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
				# not identified m/zs - novelty UNPD
				unpd_gnps_statistics['Statistics'].append("Number of putative novel [M+H]+ m/zs in UNPD best")
				unpd_gnps_statistics['Value'].append(
					f"{(number_protonated_mzs - number_protonated_identified_unpd)} ({(number_protonated_mzs - number_protonated_identified_unpd) / number_protonated_mzs * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"Total number of "+mz_types+" [M+H]+ m/zs that were NOT identified against the UNPD using tremolo best result (tremolo_best_position == NA), thus, may represent putative novel compounds not present in the database. And its percentage over the total number of "+mz_types+" [M+H]+ m/zs.")
				number_protonated_unique_identification_unpd = clean_data.loc[
					(~clean_data.tremolo_SMILES_best.isna()), "tremolo_SMILES_best"].unique().size
				# unique all library annotations UNPD
				unpd_gnps_statistics['Statistics'].append(
					"Number of unique UNPD best library annotations for the [M+H]+ m/zs and their UNPD coverage")
				unpd_gnps_statistics['Value'].append(
					f"{number_protonated_unique_identification_unpd} ({number_protonated_unique_identification_unpd / total_unpd_fixo * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"The number of unique best library annotations in UNPD for the "+mz_types+" [M+H]+ m/zs (unique tremolo_SMILES_best for tremolo_best_position != NA ) and its percentage over the total number of unique SMILES in UNPD (" + str(
						total_unpd_fixo) + " compounds for UNPD total unique SMILES) - the UNPD coverage by [M+H]+.")
				# chemical diversity UNPD
				# use superclass clean, remove NA and None
				number_unique_superclass = np.unique([x for x in chain.from_iterable(
					clean_data.tremolo_NPClassifier_superclass_clean_best.str.split(":", expand=True).to_numpy())
				                                      if x == x and x is not None and x is not ""]).size
				unpd_gnps_statistics['Statistics'].append(
					"Chemical diversity of [M+H]+ m/zs in UNPD Superclasses best")
				unpd_gnps_statistics['Value'].append(
					f"{number_unique_superclass} ({number_unique_superclass / total_superclass_npclassifier * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"The unique number of best superclasses that got identified by the "+mz_types+" [M+H]+ m/zs in UNPD (column tremolo_NPClassifier_superclass_clean_best not NA or empty). And its percentage over the total number of unique superclasses considered (" + str(
						total_superclass_npclassifier) + " for UNPD using NPClassifier). ")
				number_unique_superclass_grouping = np.unique(clean_data.tremolo_NPClassifier_superclass_grouping_best[(
						clean_data.tremolo_NPClassifier_superclass_grouping_best != "Not_Annotated")].values).size
				unpd_gnps_statistics['Statistics'].append(
					"Chemical diversity of [M+H]+ m/zs in UNPD Superclasses best grouping")
				unpd_gnps_statistics['Value'].append(
					f"{number_unique_superclass_grouping} ({number_unique_superclass_grouping / total_superclass_npclassifier_grouping * 100:.1f}%)")
				unpd_gnps_statistics['Description'].append(
					"The unique number of best superclasses grouping that got identified by the "+mz_types+" [M+H]+ m/zs in UNPD (column tremolo_NPClassifier_superclass_grouping_best not NA or empty). And its percentage over the total number of unique superclasses grouping considered (" + str(
						total_superclass_npclassifier_grouping) + " groups proposed by NP3 from the NPClassifier superclasses). ")
	#
	# save results
	# save the gnps statistics
	gnps_statistics_table = pd.DataFrame(gnps_statistics)
	gnps_statistics_table.to_csv(output_path / "chemical_identification_statistics_GNPS.csv", index=False)
	# save the unpd and gnps best statistics
	if curated_lib_annotation_exists:
		unpd_gnps_statistics_table = pd.DataFrame(unpd_gnps_statistics)
		unpd_gnps_statistics_table.to_csv(output_path / "chemical_identification_statistics_UNPDxGNPS_final_curation.csv", index=False)

