#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26 Sep 2025
@author: Cris

Script to select the best result from the UNPD identification using tremolo
Selects the best tremolo result with the greater mqscore and mzerror < 20, if any, or >= 20 within conditions of number of matched peaks and mqscore range
   Uses a defined scoring and classification system to select the best results
  creates the coluns tremolo_<>_best with the best selected identification result
  creates the coluns tremolo_UNPD_score_best and tremolo_UNPD_category_best with the final scores of the best identification and its classification
  creates the column curated_identification_best_origin to store the current best origin of valid identifications - checked as 'UNPD' for all valid classifications
  cretes the curated_superclass columns and the curated superclasses groupings columns with the counts of occurrence in each group fro the best origin, here UNPD only
Stores the results in the clean table - overwritting it

Args:

	clean table csv
"""

import pandas as pd
import numpy as np
from pathlib import Path

# return the position of the first float with value < value_lt from a list of string with floats
# - which will be the one with the greater mqscore, following the tremolo result ordering -,
# else return the position of the first smaller float rounded in decimal places that is greater than the cutoff value
def position_first_float_list_lt_or_smaller(str_list, value_lt):
	x_values = np.float64(str_list.split(";"))
	x = (x_values < value_lt)
	if any(x):
		return int(np.argmax(x))
	else:
		return np.argmin(np.round(x_values/10)*10)

def get_position_str(x,pos,col):
	if type(x) is str:
		if col == 'tremolo_SMILES':
			x = x.split(",")
		else:
			x = x.split(";")
		return x[int(pos)]
	else:
		return x

superclass_groupings_names = ['Aminoacids_and_Peptides_and_OtherNComp', 'Alkaloids_and_Lactams',
	                               'Terpenes_and_Carotenoids', 'Fatty_Acids_and_Lipids', 'Polyketides',
	                               'Benzenoids', 'Flavonoids_and_Phenolic_derivatives', 'Organohalogen_and_Organometallic',
	                               'Lignans_and_Other_Ocompounds', 'Organic_Acids_and_OthersGenerals', 'Not_Annotated']

# add curate superclass groupings
# groupings of the curated superclass in 10 major groups plus the not annotated ones
# count number of occurrences in each group - a single identification may have 2 superclasses separated by :
# return a dataframe with the columns equal to the 10 major groups with their occurrences counts 0,1,2,...
def group_curated_superclass_toCols(curated_superclass_col):
	# check if there is at least one not NaN value, otherwise skip the curation
	if not curated_superclass_col.isna().all():
		# create dataframe with curated superclass split in columns
		curated_superclass_col = curated_superclass_col.str.split(":", expand=True)
		n = curated_superclass_col.shape[1]
		# for each superclass split, check their groupings and compute their counts by major group
		# create the 10 major superclasses plus not annotated cols and initialize them with 0
		
		curated_superclass_col.loc[:, superclass_groupings_names] = 0
		# simplified classification of the superclass in 10 defined groups
		# Groupingss definition:
		superclass_groupings = {
			'Aminoacids_and_Peptides_and_OtherNComp':
				["Amino acid glycosides", "Mycosporine derivatives", "Oligopeptides",
				"Small peptides", "Aminosugars and aminoglycosides", "Diazotetronic acids and derivatives",
				"Mitomycin derivatives",
				"Organic nitrogen compounds", "Nucleosides, nucleotides, and analogues", "Nucleosides"],
		    'Alkaloids_and_Lactams':
			    ["Alkaloids and derivatives", "Anthranilic acid alkaloids", "Lysine alkaloids", "Guanidine alkaloids",
				"Histidine alkaloids", "Ornithine alkaloids",
				"Peptide alkaloids", "Proline alkaloids", "Pseudoalkaloids", "Tetramate alkaloids",
				"Miscellaneous alkaloids", "Nicotinic acid alkaloids",
				"Serine alkaloids", "Tryptophan alkaloids", "Tyrosine alkaloids", "β-lactams", "γ-lactam-β-lactones"],
			'Terpenes_and_Carotenoids':
				["Diterpenoids", "Meroterpenoids", "Monoterpenoids",
				"Sesterterpenoids", "Triterpenoids", "Steroids", "Sesquiterpenoids",
				"Apocarotenoids", "Carotenoids (C40)", "Carotenoids (C45)", "Carotenoids (C50)"],
			'Fatty_Acids_and_Lipids':
				["Fatty Acids and Conjugates", "Fatty acyl glycosides", "Fatty acyls",
				"Fatty amides", "Fatty esters", "Lipids and lipid-like molecules", "Docosanoids",
				"Eicosanoids", "Glycerophospholipids", "Glycerolipids", "Octadecanoids",
				"Sphingolipids", "Hydrocarbons", "Hydrocarbon derivatives"],
			'Polyketides':
				["Aromatic polyketides", "Cyclic polyketides", "Macrolides", "Linear polyketides",
				"Miscellaneous polyketides", "Polycyclic aromatic polyketides", "Polyethers"],
			'Benzenoids':
				["Benzenoids", "Diarylheptanoids", "Naphthalenes", "Phenanthrenoids",
				"Stilbenoids", "Terphenyls"],
			'Flavonoids_and_Phenolic_derivatives':
				["Alkylresorcinols", "Flavonoids", "Isoflavonoids", "Phenolic acids (C6-C1)",
				"Phenylethanoids (C6-C2)", "Phenylethanoids (C6-C3)", "Phenylpropanoids (C6-C3)",
				"Phloroglucinols", "Phenylpropanoids and polyketides"],
			'Organohalogen_and_Organometallic':
				["Organohalogen compounds", "Organometallic compounds", "Acetylides"],
			'Lignans_and_Other_Ocompounds':
				["Chromanes", "Coumarins", "Organic oxygen compounds",
				"Polyols", "Polyprenols", "Saccharides", "Styrylpyrones", "Xanthones",
				"Lignans", "Lignans, neolignans and related compounds"],
			'Organic_Acids_and_OthersGenerals':
				["Homogeneous non-metal compounds", "Organic 1,3-dipolar compounds", "Organic acids and derivatives",
				"Organoheterocyclic compounds"],
			'Not_Annotated': ["not_annotated"]}
		
		# iterate over the names of the superclasses groupings and count the curated superclasses (n first cols)
		#  that appear in each one of the groupings
		for superclass_group in superclass_groupings_names:
			for i in range(0, n):
				check_superclass_grouping = curated_superclass_col.loc[:, i].isin(superclass_groupings[superclass_group])
				if check_superclass_grouping.any():
					curated_superclass_col.loc[check_superclass_grouping, superclass_group] += 1
		# check the not annotated entries, the ones without any counts in the superclasses groupings
		curated_superclass_col.loc[(curated_superclass_col.loc[:, superclass_groupings_names].sum(1) == 0),
		                           'Not_Annotated'] = 1
		# filter only the superclasses grouping columns
		curated_superclass_col = curated_superclass_col.iloc[:, n:]
		# choose the most representative superclass grouping of each m/z depending on its counts
		curated_superclass_col['curated_superclass_grouping'] = curated_superclass_col.idxmax(axis=1)
		# renamed the superclass columns using the prefix curated_superclass_GR_
		curated_superclass_col.columns = list("curated_superclass_GR_" + curated_superclass_col.columns[:-1]) + [curated_superclass_col.columns[-1]]
	else:
		# not valid superclass, set all as not annotated
		curated_superclass_col = pd.DataFrame(curated_superclass_col)
		curated_superclass_col.loc[:, "curated_superclass_GR_Not_Annotated"] = 1
		curated_superclass_col.loc[:, "curated_superclass_grouping"] = "Not_Annotated"
		curated_superclass_col = curated_superclass_col.iloc[:, 1:]
	# return the superclasses groupings counts
	return curated_superclass_col

def curate_tremolo_unpd_identification(clean_table_file):
	clean_table_file = Path(clean_table_file)
	if not clean_table_file.exists() or not clean_table_file.is_file():
		sys.exit(
			"The provided path to the clean table file for UNPD curation does not exists." +
			" Tremolo-UNPD identification curation aborted.")
	
	clean_table = pd.read_csv(clean_table_file, converters={'msclusterID':str}, low_memory=False)
	
	# check if the tremolo result is present, if not exit without further processing
	if not 'tremolo_UNPD_IDs' in clean_table.columns:
		sys.exit("- The Tremolo-UNPD result is not present in the provided count table, curation skipped.")
	
	print("  * Creating the curated identification result for Tremolo-UNPD *")
	# Create the tremolo identification best result column

	# for each mz, select first tremolo result with tremolo_mzErrorPPM < 20 or the first minimum mzErrorPPM
	# and return its position in the results
	# this will be the better result with greater mqscore and mzerror < 20 or minimum possible
	# filter all tremolo columns results to get the position of the best result and save to the tremolo_<>_best columns
	# create the scoring and category columns and store the valid best results in the curated_identification column
	# then overwrite corresponding clean table
	clean_table["tremolo_best_position"] = [position_first_float_list_lt_or_smaller(x, 20) if type(x) is str else None
													   for x in clean_table.tremolo_mzErrorPPM]

	# apply function to create modified columns with best identification in the selected position
	columns_to_modify = ['tremolo_MQScore', 'tremolo_NPClassifier_superclass', 'tremolo_NPClassifier_class',
						 'tremolo_NPClassifier_pathway', 'tremolo_molecularFormula', 'tremolo_chemicalNames',
						 'tremolo_InChIKey', 'tremolo_SMILES', 'tremolo_UNPD_IDs', 'tremolo_numSharedPeaks',
						 'tremolo_mzErrorPPM']

	for col in columns_to_modify:
		clean_table[f'{col}_best'] = [get_position_str(clean_table.loc[i,col],pos,col)
												if pos==pos else pos
												for i,pos in enumerate(clean_table.tremolo_best_position)]
		# convert to decimal columns
		if col in ['tremolo_mzErrorPPM', 'tremolo_MQScore', 'tremolo_numSharedPeaks']:
			clean_table[f'{col}_best'] = pd.to_numeric(clean_table[f'{col}_best'], errors='coerce')

	# add unpd category and scoring
	# Define conditions and categories
	conditions = [
		{
			"condition": (clean_table['tremolo_mzErrorPPM_best'] <= 20) & (clean_table['tremolo_numSharedPeaks_best'] >= 6) & (
						clean_table['tremolo_MQScore_best'] >= 0.5),
			"category": "Top1|mqs>=0.5"
		},
		{
			"condition": (clean_table['tremolo_mzErrorPPM_best'] <= 20) & (clean_table['tremolo_numSharedPeaks_best'] >= 6) & (
						clean_table['tremolo_MQScore_best'] >= 0.4) & (clean_table['tremolo_MQScore_best'] < 0.5),
			"category": "Top2|0.4<=mqs<0.5"
		},
		{
			"condition": (clean_table['tremolo_mzErrorPPM_best'] <= 20) & (clean_table['tremolo_numSharedPeaks_best'] < 6) & (
						clean_table['tremolo_MQScore_best'] >= 0.5),
			"category": "Top3|sp<6 mqs>=0.5"
		},
		{
			"condition": (clean_table['tremolo_mzErrorPPM_best'] <= 20) & (clean_table['tremolo_numSharedPeaks_best'] < 6) & (
						clean_table['tremolo_MQScore_best'] >= 0.4) & (clean_table['tremolo_MQScore_best'] < 0.5),
			"category": "Top4|sp<6 0.4<=mqs<0.5"
		},
		{
			"condition": (clean_table['tremolo_mzErrorPPM_best'] <= 20) & (
						clean_table['tremolo_numSharedPeaks_best'] >= 6) & (clean_table['tremolo_MQScore_best'] >= 0.2),
			"category": "Top5|mqs>=0.2"
		},
		{
			"condition": (clean_table['tremolo_mzErrorPPM_best'] > 20) & (clean_table['tremolo_numSharedPeaks_best'] >= 6) & (
						clean_table['tremolo_MQScore_best'] >= 0.5),
			"category": "Analog1|mqs>=0.5"
		},
		{
			"condition": (clean_table['tremolo_mzErrorPPM_best'] > 20) & (clean_table['tremolo_numSharedPeaks_best'] >= 6) & (
						clean_table['tremolo_MQScore_best'] >= 0.4) & (clean_table['tremolo_MQScore_best'] < 0.5),
			"category": "Analog2|0.4<=mqs<0.5"
		},
		{
			"condition": (clean_table['tremolo_mzErrorPPM_best'] > 20) & (clean_table['tremolo_numSharedPeaks_best'] < 6) & (
						clean_table['tremolo_MQScore_best'] >= 0.5),
			"category": "Analog3|sp<6 mqs>=0.5"
		},
		{
			"condition": (clean_table['tremolo_mzErrorPPM_best'] > 20) & (clean_table['tremolo_numSharedPeaks_best'] < 6) & (
						clean_table['tremolo_MQScore_best'] >= 0.4) & (clean_table['tremolo_MQScore_best'] < 0.5),
			"category": "Analog4|sp<6 0.4<=mqs<0.5"
		},
		{
			"condition": (clean_table['tremolo_mzErrorPPM_best'] > 20) & (clean_table['tremolo_numSharedPeaks_best'] >= 6) & (
						clean_table['tremolo_MQScore_best'] >= 0.2),
			"category": "Analog5|mqs>=0.2"
		}
	]
	# Apply the conditions, if non is filtered set it as out
	clean_table['tremolo_UNPD_category_best'] = 'out'  # default
	for cond in conditions:
		clean_table.loc[cond["condition"], 'tremolo_UNPD_category_best'] = cond["category"]

	# Mapping from the categories and the scores
	score_mapping = {
		"Top1|mqs>=0.5": 100,
		"Top2|0.4<=mqs<0.5": 80,
		"Top3|sp<6 mqs>=0.5": 60,
		"Top4|sp<6 0.4<=mqs<0.5": 50,
		"Top5|mqs>=0.2": 45,
		"Analog1|mqs>=0.5": 40,
		"Analog2|0.4<=mqs<0.5": 30,
		"Analog3|sp<6 mqs>=0.5": 20,
		"Analog4|sp<6 0.4<=mqs<0.5": 10,
		"Analog5|mqs>=0.2": 5,
		"out": 0  # default pontuations for not mapped category
	}

	# Creates the new column "tremolo_UNPD_score_best" with the mapped score based on the category of the best identification
	clean_table['tremolo_UNPD_score_best'] = clean_table['tremolo_UNPD_category_best'].map(score_mapping)

	# defines the current best source/origin of the identificaiton - only UNPD
	clean_table['curated_identification_best_origin'] = ''
	clean_table.loc[clean_table.tremolo_UNPD_score_best > 0,'curated_identification_best_origin'] = 'UNPD'
	
	# make the curated superclass for UNPD identification
	# get the first superclass from NPClassifier, the one before the pipe "|", if not NA else leave as NA
	clean_table['tremolo_curated_superclass'] = clean_table.tremolo_NPClassifier_superclass_best.apply(lambda x: x.split("|")[0] if x == x else np.NaN)
	# now call function to create the superclass grouping - counting the number of superclass match by group in the case
	# of doubled superclasses separated by ':'
	curated_superclass_groupings = group_curated_superclass_toCols(clean_table['tremolo_curated_superclass'])
	curated_superclass_groupings.columns = "tremolo_"+curated_superclass_groupings.columns
	if curated_superclass_groupings.columns.isin(clean_table.columns.values).all():
		clean_table.drop(curated_superclass_groupings.columns.values, axis=1, inplace=True)
	clean_table = pd.concat([clean_table, curated_superclass_groupings], axis=1)
	
	print("  - Saving the clean table with the UNPD curated identification result and superclasses groupings: ", clean_table_file)
	# save the data with the new columns - overwrite original table
	clean_table.to_csv(clean_table_file, index=False, float_format="%.4f")
	# save the curated columns to the other count table if it exists
	if clean_table_file.name.find("peak_area") > 0:
		clean_table_other_file = clean_table_file.parent / clean_table_file.name.replace("peak_area", "spectra")
	elif clean_table_file.name.find("spectra") > 0:
		clean_table_other_file = clean_table_file.parent / clean_table_file.name.replace("spectra", "peak_area")
	else:
		# no other clean table to save the curation result
		sys.exit(0)
	if clean_table_other_file.exists() and clean_table_other_file.is_file():
		print("  - Saving the clean table with the UNPD curated identification result and superclasses groupings: ",
		      clean_table_other_file)
		new_curated_columns = np.concatenate([["tremolo_best_position"], clean_table.columns.values[
			(clean_table.columns.str.startswith("tremolo_") & clean_table.columns.str.endswith("_best")) |
			(clean_table.columns.str.startswith("tremolo_curated_superclass"))],
		                                      ["curated_identification_best_origin"]])
		clean_table_other = pd.read_csv(clean_table_other_file, converters={'msclusterID':str}, low_memory=False)
		clean_table_other.drop(new_curated_columns, axis=1, inplace=True, errors="ignore") # remove existing new columns
		clean_table_other = clean_table_other.merge(
			clean_table.loc[:, np.concatenate([["msclusterID"], new_curated_columns])], how="left",
			on="msclusterID")
		clean_table_other.to_csv(clean_table_other_file, index=False, float_format="%.4f")


if __name__ == "__main__":
	import sys, os
	if len(sys.argv) > 1:
		# print(sys.argv)
		clean_table_file = sys.argv[1]
	else:
		print("Error: One argument must be supplied to curate the UNPD identification result using tremolo:\n",
			  " 1 - clean_table_file: Path to the clean table containing identification results from UNPD using tremolo.\n")
		sys.exit(1)
	# call the curate function
	curate_tremolo_unpd_identification(clean_table_file)
	