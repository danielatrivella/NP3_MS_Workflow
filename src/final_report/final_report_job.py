from chemical_report_statistics import compute_chemical_report_statistics, plot_superclass_samples_distribution
from pca_calculation_ref_plot import pca_calculation_smiles_rcdk_ref_plot,pca_calculation_mz_ref_plot
from quantification_report_statistics import compute_quantification_report_statistics
from molecular_networking_report_statistics import compute_mn_report_statistics
from pathlib import Path
import sys

# creates the final report by calling the chemical, quantification and pca scripts
# the output_path should be the final result folder inside the outs folder, named with the output_name
# inside it will be created the final_report folder
def run_final_report(metadata_file, clean_table_path, output_path, output_name, mz_tolerance):
	output_path = Path(output_path)
	if not output_path.exists() or not output_path.is_dir():
		sys.exit("The provided output path does not exists. Final report aborted.")
	clean_table_path = Path(clean_table_path)
	if not clean_table_path.exists() or not clean_table_path.is_file():
		sys.exit("The provided path to the clean data file does not exists. Final report aborted.")
	metadata_file = Path(metadata_file)
	if not metadata_file.exists() or not metadata_file.is_file():
		sys.exit("The provided path to the metadata table file does not exists. Final report aborted.")
	
	print("\n*** Creating the NP3 final report ***\n\n")
	# create the final report folder
	final_report_path = (output_path / "final_reports")
	final_report_path.mkdir(exist_ok=True)
	if not final_report_path.exists() or not final_report_path.is_dir():
		sys.exit("The final report path could not be created. Final report aborted.")
	
	print("** Creating the NP3 final chemical report **\n")
	# create the chemical report subfolder inside the final_report
	chemical_report_path = final_report_path / "chemical_report"
	chemical_report_path.mkdir(exist_ok=True)
	if not chemical_report_path.exists() or not chemical_report_path.is_dir():
		print("The final chemical report path could not be created. Chemical report statistics aborted.")
	else:
		# call function to create the chemical report statistics
		compute_chemical_report_statistics(clean_table_path, chemical_report_path)
		# call function to create the superclasses grouping distribution by sample
		plot_superclass_samples_distribution(metadata_file, clean_table_path, chemical_report_path)
		# create the chemical space identification subfolder inside the chemical_report folder
		chemical_space_identification_path = chemical_report_path / "chemical_space_identifications"
		chemical_space_identification_path.mkdir(exist_ok=True)
		if not chemical_space_identification_path.exists() or not chemical_space_identification_path.is_dir():
			print("The final chemical space identification path could not be created. PCA for identifications aborted.")
		else:
			# call function to create the pca of chemical space identification using the UNPD result
			# get the directory of the current script and point to the np3 pca reference table
			np3_chemical_space_reference_path = Path(__file__).resolve().parent / "Chemical_space_data" / "descriptors_reference_unpd_drugbank_allo_rev_natural_pubmedID_clean_top24.csv"
			pca_calculation_smiles_rcdk_ref_plot(np3_chemical_space_reference_path, clean_table_path,
												 chemical_space_identification_path, output_name, data_type="UNPD")
		
		# create the chemical space identification subfolder inside the chemical_report folder
		chemical_space_mzs_path = chemical_report_path / "chemical_space_mzs"
		chemical_space_mzs_path.mkdir(exist_ok=True)
		if not chemical_space_mzs_path.exists() or not chemical_space_mzs_path.is_dir():
			print("The final chemical space mzs path could not be created. PCA for m/z quantifications aborted.")
		else:
			# call creation of the chemical space for the mzs quantification
			pca_calculation_mz_ref_plot(clean_table_path, chemical_space_mzs_path, output_name)
	
	print("\n** Creating the NP3 final quantification report **\n")
	# create the quantification report subfolder inside the final_report
	quantification_report_path = final_report_path / "quantification_report"
	quantification_report_path.mkdir(exist_ok=True)
	if not quantification_report_path.exists() or not quantification_report_path.is_dir():
		print("The final quantification report path could not be created. Quantificaiton report aborted.")
	else:
		# call function to create the quantification reports
		compute_quantification_report_statistics(clean_table_path, quantification_report_path, mz_tolerance)
	
	print("** Creating the NP3 final molecular networking report **\n")
	# create the molecular networking report subfolder inside the final_report
	mn_report_path = final_report_path / "molecular_networking_report"
	mn_report_path.mkdir(exist_ok=True)
	if not mn_report_path.exists() or not mn_report_path.is_dir():
		print("The final molecular networking report path could not be created. Molecular networking report aborted.")
	else:
		# call function to create the molecular networking reports for each
		# network present in the output_mn_path
		output_mn_path = (output_path / "molecular_networking")
		compute_mn_report_statistics(output_mn_path, mn_report_path)


if __name__ == "__main__":
	if len(sys.argv) > 5:
		metadata_path = sys.argv[1]
		clean_table_path = sys.argv[2]
		output_path = sys.argv[3]
		output_name = sys.argv[4]
		mz_tolerance = float(sys.argv[5])
	else:
		print("Error: Five arguments must be supplied to created the final report of the NP3 result:\n"
			"  1 - metadata_path: Path to the job samples metadata table with the original samples codes and types (.csv);\n"
			"  2 - clean_table_path: Path to the clean table with the final list of consensus spectra and UNPD identification if any (.csv);\n"
			"  3 - output_path: Path to the final result folder, inside the outs folder, named with the output_name;\n"
			"  4 - output_name: The name of the job;\n"
		    "  5 - mz_tolerance: m/z tolerance in daltons for the precursor m/z.\n")
		sys.exit(1)
	# call final report creation
	run_final_report(metadata_path, clean_table_path, output_path, output_name, mz_tolerance)
	