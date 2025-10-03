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
	if "curated_superclass_grouping" in clean_data.columns:
		print("  - Creating the superclass grouping distribution by not blank sample \n")
		# get the columns names containing the count of spectra by peak area without blanks
		samples_area_name = metadata.SAMPLE_CODE[metadata.SAMPLE_TYPE.str.lower() != "blank"].values
		samples_area_col = samples_area_name + "_area"
		
		# group the quantification columns by the superclass grouping and sum the respective rows
		samples_area_by_superclass_grouping = clean_data.groupby("curated_superclass_grouping")[samples_area_col].sum()
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
	
