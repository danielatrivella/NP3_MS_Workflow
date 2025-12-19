#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08 Octu 2025
@author: Cris

Script to select the best result from the GNPS identification Library Search or Molecular Networking workflows
Compare the gnps smiles with the best tremolo smiles when present, then categorize the gnps result using mqscore,
mzerror and number of shared peaks. Uses a defined scoring and classification system to curate the results and select
best results.
  creates the coluns gnps_score and gnps_category with the final scores of the best identification and its classification
  creates the column curated_identification_best_origin to store the current best origin of valid identifications - checked as 'UNPD' for all valid classifications
  creates the gnps_curated_superclass columns and the curated superclasses groupings columns with the counts of occurrence in each group for GNPS and the best origin (UNPDxGNPS)
Stores the results in the clean table as new coloumns - overwritting it

Args:

	clean table csv
	output_path to the final result folder of the np3 job, inside the outs dir, named with the output_name,
		it will be used to store somes analysis in the final_reports folder
	metadata path to the metadata table in csv
"""

import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from tremolo_UNPD_curate_identification import group_curated_superclass_toCols
from pathlib import Path
from pca_calculation_ref_plot import pca_calculation_smiles_rcdk_ref_plot
from chemical_report_statistics import compute_chemical_identification_report_GNPS_result, plot_superclass_samples_distribution
import os

# funcao para calcular o coeficiente de Tanimoto
def calculate_tanimoto(smiles1, smiles2):
	# Verifica se os SMILES sao validos
	if pd.isna(smiles1) or pd.isna(smiles2) or smiles1.strip() == '' or smiles2.strip() == '':
		return None
	# Tenta converter os SMILES em moleculas
	try:
		mol1 = Chem.MolFromSmiles(smiles1)
		mol2 = Chem.MolFromSmiles(smiles2)
		# Verifica se os SMILES sao validos
		if mol1 is None or mol2 is None:
			return None
		# Gera os fingerprints
		mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
		fp1 = mfpgen.GetFingerprint(mol1)
		fp2 = mfpgen.GetFingerprint(mol2)

		# Calcula o coeficiente de Tanimoto
		tanimoto = DataStructs.TanimotoSimilarity(fp1, fp2)
		return min(round(tanimoto, 4), 1.0)  # Arredonda o valor para 4 casas decimais
	except Exception as e:
		print(f"Error calculating Tanimoto for SMILES: {smiles1} and {smiles2}. Error: {e}")
		return None
	

def curate_gnps_identification(clean_table_file, output_path, metadata_file):
	clean_table_file = Path(clean_table_file)
	if not clean_table_file.exists() or not clean_table_file.is_file():
		sys.exit(
			"The provided path to the clean table file for GNPS curation does not exists." +
			" GNPS identification curation aborted.")
	
	clean_table = pd.read_csv(clean_table_file, converters={'msclusterID':str}, low_memory=False)
	
	if not 'gnps_Smiles' in clean_table.columns:
		sys.exit("- The GNPS result is not present in the provided count table, curation skipped.")
		
	print("* Creating the curated identification result for GNPS *")
	if 'tremolo_SMILES_best' in clean_table.columns:
		print("* Computing Tanimoto between UNPDxGNPS smiles *")
		# Calculate tanimoto score between the best unpd and the best gnps for each mz
		clean_table['tanimoto_unpd_gnps'] = clean_table.apply(
			lambda row: calculate_tanimoto(row['tremolo_SMILES_best'],
										   row['gnps_Smiles']), axis=1)
		
	# Create column gnps_GoldCategory = GOLD if mzerrorppm <= 20 and mqscore >= 0.9 else out
	clean_table['gnps_GoldCategory'] = 'out'
	clean_table.loc[
		(clean_table.gnps_MZErrorPPM <= 20) & (clean_table.gnps_MQScore >= 0.9), 'gnps_GoldCategory'] = 'GOLD'
	
	# Add GNPS category and scoring
	# Define conditions
	conditions = [
		{
			"condition": (clean_table['gnps_MZErrorPPM'] <= 20) & (clean_table['gnps_SharedPeaks'] >= 6) & (
					clean_table['gnps_MQScore'] >= 0.9),
			"category": "Top1|mqs>=0.9"
		},
		{
			"condition": (clean_table['gnps_MZErrorPPM'] <= 20) & (clean_table['gnps_SharedPeaks'] >= 6) & (
					clean_table['gnps_MQScore'] >= 0.8) & (clean_table['gnps_MQScore'] < 0.9),
			"category": "Top2|0.8<=mqs<0.9"
		},
		{
			"condition": (clean_table['gnps_MZErrorPPM'] <= 20) & (clean_table['gnps_SharedPeaks'] < 6) &
						 (clean_table['gnps_MQScore'] >= 0.9),
			"category": "Top3|sp<6 mqs>=0.9"
		},
		{
			"condition": (clean_table['gnps_MZErrorPPM'] <= 20) & (clean_table['gnps_SharedPeaks'] < 6) & (
					clean_table['gnps_MQScore'] >= 0.8) & (clean_table['gnps_MQScore'] < 0.9),
			"category": "Top4|sp<6 0.8<=mqs<0.9"
		},
		{
			"condition": (clean_table['gnps_MZErrorPPM'] <= 20) & (clean_table['gnps_SharedPeaks'] >= 6) & (
					clean_table['gnps_MQScore'] >= 0.7) & (clean_table['gnps_MQScore'] < 0.8),
			"category": "Top5|mqs>=0.7"
		},
		{
			"condition": (clean_table['gnps_MZErrorPPM'] > 20) & (clean_table['gnps_SharedPeaks'] >= 6) &
						 (clean_table['gnps_MQScore'] >= 0.9),
			"category": "Analog1|mqs>=0.9"
		},
		{
			"condition": (clean_table['gnps_MZErrorPPM'] > 20) & (clean_table['gnps_SharedPeaks'] >= 6) & (
					clean_table['gnps_MQScore'] >= 0.8) & (clean_table['gnps_MQScore'] < 0.9),
			"category": "Analog2|0.8<=mqs<0.9"
		},
		{
			"condition": (clean_table['gnps_MZErrorPPM'] > 20) & (clean_table['gnps_SharedPeaks'] < 6) &
						 (clean_table['gnps_MQScore'] >= 0.9),
			"category": "Analog3|sp<6 mqs>=0.9"
		},
		{
			"condition": (clean_table['gnps_MZErrorPPM'] > 20) & (clean_table['gnps_SharedPeaks'] < 6) & (
					clean_table['gnps_MQScore'] >= 0.8) & (clean_table['gnps_MQScore'] < 0.9),
			"category": "Analog4|sp<6 0.8<=mqs<0.9"
		},
		{
			"condition": (clean_table['gnps_MZErrorPPM'] > 20) & (clean_table['gnps_SharedPeaks'] >= 6) & (
					clean_table['gnps_MQScore'] >= 0.7) & (clean_table['gnps_MQScore'] < 0.8),
			"category": "Analog5|mqs>=0.7"
		}
	]
	# Apply the conditions
	clean_table['gnps_category'] = 'out'  # default
	for cond in conditions:
		clean_table.loc[cond["condition"], 'gnps_category'] = cond["category"]
	
	# Map the categories to the scores
	score_mapping = {
		"Top1|mqs>=0.9": 100,
		"Top2|0.8<=mqs<0.9": 80,
		"Top3|sp<6 mqs>=0.9": 60,
		"Top4|sp<6 0.8<=mqs<0.9": 50,
		"Top5|mqs>=0.7": 45,
		"Analog1|mqs>=0.9": 40,
		"Analog2|0.8<=mqs<0.9": 30,
		"Analog3|sp<6 mqs>=0.9": 20,
		"Analog4|sp<6 0.8<=mqs<0.9": 10,
		"Analog5|mqs>=0.7": 5,
		"out": 0  # not mapped category
	}
	# Create the column gnps_score to store the mapping result
	clean_table['gnps_score'] = clean_table['gnps_category'].map(score_mapping)
	
	# define the best source of identification, if tremolo UNPD present compare their score
	if 'tremolo_UNPD_score_best' in clean_table.columns:
		clean_table.loc[(clean_table.gnps_score > 0) &
						(clean_table.gnps_score >= clean_table.tremolo_UNPD_score_best),
						'curated_identification_best_origin'] = 'GNPS'
	else:
		# no tremolo unpd result, set best source as empty and select the GNPS results
		clean_table['curated_identification_best_origin'] = ''
		clean_table.loc[(clean_table.gnps_score > 0), 'curated_identification_best_origin'] = 'GNPS'
		
	# make the clean and the curated superclass for GNPS identification
	# get the first superclass from NPClassifier, the one before the pipe "|", if not NA else leave as NA
	clean_table['gnps_npclassifier_superclass_clean'] = clean_table.gnps_npclassifier_superclass.apply(
		lambda x: x.split("|")[0] if x == x else np.NaN)
	clean_table['gnps_curated_superclass'] = clean_table['gnps_npclassifier_superclass_clean']
	clean_table.loc[clean_table.gnps_category == "out", "gnps_curated_superclass"] = np.NaN
	# now call function to create the superclass grouping - counting the number of superclass match by group in the case
	# of doubled superclasses separated by ':'
	curated_superclass_groupings = group_curated_superclass_toCols(clean_table['gnps_npclassifier_superclass_clean'])
	curated_superclass_groupings.columns = "gnps_" + curated_superclass_groupings.columns
	# set the clean grouping and then set the curated result to retain only the gnps_category != "out" results
	clean_table['gnps_npclassifier_superclass_grouping'] = curated_superclass_groupings['gnps_curated_superclass_grouping']
	# set the superclasses grouping count columns to 0 where gnps category is out and the grouping to NA
	curated_superclass_groupings.loc[(clean_table.gnps_category == "out"), :] = 0
	curated_superclass_groupings.loc[(clean_table.gnps_category == "out"),
									 'gnps_curated_superclass_grouping'] = "Not_Annotated"
	curated_superclass_groupings.loc[(clean_table.gnps_category == "out"),
									 'gnps_curated_superclass_GR_Not_Annotated'] = 0
	if curated_superclass_groupings.columns.isin(clean_table.columns.values).all():
		clean_table.drop(curated_superclass_groupings.columns.values, axis=1, inplace=True)
	clean_table = pd.concat([clean_table, curated_superclass_groupings], axis=1)
	
	# and now curate the superclass of the best origin, if tremolo result is present
	# extract the best origin smiles and superclass and then group them
	if 'tremolo_curated_superclass' in clean_table.columns:
		clean_table['best_origin_curated_superclass'] = ''
		clean_table['best_origin_SMILES'] = ''
		# get superclass for best origin gnps
		best_origin_gnps = (clean_table.curated_identification_best_origin == 'GNPS')
		clean_table.loc[best_origin_gnps, 'best_origin_curated_superclass'] = clean_table.loc[
			best_origin_gnps, 'gnps_curated_superclass']
		clean_table.loc[best_origin_gnps, 'best_origin_SMILES'] = clean_table.loc[
			best_origin_gnps, 'gnps_Smiles']
		# get super class best origin unpd
		best_origin_unpd = (clean_table.curated_identification_best_origin == 'UNPD')
		clean_table.loc[best_origin_unpd, 'best_origin_curated_superclass'] = clean_table.loc[
			best_origin_unpd, 'tremolo_curated_superclass']
		clean_table.loc[best_origin_unpd, 'best_origin_SMILES'] = clean_table.loc[
			best_origin_unpd, 'tremolo_SMILES_best']
		# create grouping for best_origin_curated_superclass
		curated_superclass_groupings = group_curated_superclass_toCols(clean_table['best_origin_curated_superclass'])
		curated_superclass_groupings.columns = "best_origin_" + curated_superclass_groupings.columns
		# drop previous result if exists
		clean_table.drop(curated_superclass_groupings.columns.values, axis=1, inplace=True, errors="ignore")  # remove existing new columns
		clean_table = pd.concat([clean_table, curated_superclass_groupings], axis=1)
	
	print("  - Saving the clean table with the GNPS curated identification result, superclasses groupings and best origin: ",
		  clean_table_file)
	# save the data with the new columns - overwrite original table
	clean_table.to_csv(clean_table_file, index=False, float_format="%.4f")
	# also save result to the other count table in terms of peak area or spectra
	# save the curated columns to the other count table if it exists
	if clean_table_file.name.find("peak_area") > 0:
		clean_table_other_file = clean_table_file.parent / clean_table_file.name.replace("peak_area", "spectra")
	elif clean_table_file.name.find("spectra") > 0:
		clean_table_other_file = clean_table_file.parent / clean_table_file.name.replace("spectra", "peak_area")
	else:
		# no other clean table to save the curation result, set as empty
		clean_table_other_file = Path('')
	if clean_table_other_file.exists() and clean_table_other_file.is_file():
		print("  - Saving the clean table with the GNPS curated identification result, superclasses groupings and best origin: ",
			clean_table_other_file)
		gnps_columns = ['gnps_GoldCategory', 'gnps_category', 'gnps_score', 'gnps_npclassifier_superclass_clean',
						'gnps_npclassifier_superclass_grouping', 'curated_identification_best_origin']
		if 'tanimoto_unpd_gnps' in clean_table.columns:
			gnps_columns = ['tanimoto_unpd_gnps'] + gnps_columns
		new_curated_columns = np.concatenate([gnps_columns,
											  clean_table.columns.values[
												  clean_table.columns.str.startswith("best_origin_") |
												  clean_table.columns.str.startswith("gnps_curated_superclass")]])
		clean_table_other = pd.read_csv(clean_table_other_file, converters={'msclusterID': str}, low_memory=False)
		clean_table_other.drop(new_curated_columns, axis=1, inplace=True, errors="ignore") # remove existing new columns
		clean_table_other = clean_table_other.merge(
			clean_table.loc[:, np.concatenate([["msclusterID"], new_curated_columns])], how="left",
			on="msclusterID")
		clean_table_other.to_csv(clean_table_other_file, index=False, float_format="%.4f")
		
	# call superclass grouping plot for GNPS and UNDPxGNPS, when metadata file is provided
	output_path = Path(output_path)
	chemical_report_path = (output_path / "final_reports" / "chemical_report")
	chemical_report_path.mkdir(exist_ok=True, parents=True)
	if metadata_file != "":
		plot_superclass_samples_distribution(metadata_file, clean_table_file, chemical_report_path,
											 superclass_grouping_name="gnps_curated_superclass_grouping")
		if 'best_origin_curated_superclass_grouping' in clean_table.columns:
			plot_superclass_samples_distribution(metadata_file, clean_table_file, chemical_report_path,
												 superclass_grouping_name="best_origin_curated_superclass_grouping")
	# call create report table
	compute_chemical_identification_report_GNPS_result(clean_table_file, chemical_report_path)
	# call PCA for GNPS and UNPDxGNPS
	data_reference_path = Path(os.path.dirname(__file__) , "Chemical_space_data",
							   "descriptors_reference_unpd_drugbank_allo_rev_natural_pubmedID_clean_top24.csv")
	pca_calculation_smiles_rcdk_ref_plot(data_reference_path, clean_table_file, output_path , output_path.name,
										 data_type="GNPS")
	
	
if __name__ == "__main__":
	import sys
	metadata_file = ""
	if len(sys.argv) > 2:
		# print(sys.argv)
		clean_table_file = sys.argv[1]
		output_path = sys.argv[2]
		if len(sys.argv) > 3:
			metadata_file = sys.argv[3]
	else:
		print("Error: Two argument must be supplied to curate the GNPS identification result:\n",
			  " 1 - clean_table_file: Path to the clean table containing identification results from GNPS joined.;\n",
			  " 2 - output_path: the final clustering result, inside the outs dir, named with the output_name.\n",
			  " 3 - metadata_file: path to the metadata table CSV file of the NP3 job. This is necessary to plot the",
			  "distribution of the superclasses grouping by sample, it may be missing if this plot is not desired ("
			  "leave as empty string).\n")
		sys.exit(1)
	# call the curate function
	curate_gnps_identification(clean_table_file, output_path, metadata_file)
	