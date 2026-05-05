#!/usr/bin/python

# @Crisfbazz
# adapted from https://github.com/Wang-Bioinformatics-Lab/NextflowModules/blob/a31f16297ea5d7f996edb70897179b2e190b9778/bin/library_search/library_search_wrapper.py
# changed temp folder location to the system root /tmp for executing internally in CNPEM Marvin system - not changing main root folder files

import sys
import getopt
import os
import json
import argparse
import uuid
import pandas as pd
from collections import defaultdict

def search_files(spectrum_file, library_file, temp_folder, tempresults_folder, path_to_convert, path_to_main_exec,
    min_cosine=0.7,
    min_matched_peaks=6, 
    top_k_results=1, 
    pm_tolerance=2.0, 
    ion_tolerance=0.5, 
    analog_search=0, 
    max_shift_mass=0.5,
    filter_precursor=1,
    filter_window=1):
    
    parameter_filename = os.path.join(temp_folder, str(uuid.uuid4()) + ".params")
    output_parameter_file = open(parameter_filename, "w")

    #Search Criteria
    output_parameter_file.write("MIN_MATCHED_PEAKS_SEARCH={}\n".format(min_matched_peaks))
    output_parameter_file.write("TOP_K_RESULTS={}\n".format(top_k_results))
    output_parameter_file.write("search_peak_tolerance={}\n".format(ion_tolerance))
    output_parameter_file.write("search_parentmass_tolerance={}\n".format(pm_tolerance))
    output_parameter_file.write("ANALOG_SEARCH={}\n".format(analog_search))
    output_parameter_file.write("MAX_SHIFT_MASS={}\n".format(1999))

    #Filtering Criteria
    output_parameter_file.write("FILTER_PRECURSOR_WINDOW={}\n".format(filter_precursor))
    output_parameter_file.write("MIN_PEAK_INT={}\n".format(0))
    output_parameter_file.write("WINDOW_FILTER={}\n".format(filter_window))
    output_parameter_file.write("FILTER_LIBRARY={}\n".format(1))

    #Scoring Criteria
    output_parameter_file.write("SCORE_THRESHOLD={}\n".format(min_cosine))

    #Output
    output_parameter_file.write("RESULTS_DIR={}\n".format(tempresults_folder))

    output_parameter_file.write("NODEIDX={}\n".format(0))
    output_parameter_file.write("NODECOUNT={}\n".format(1))

    output_parameter_file.write("EXISTING_LIBRARY_MGF={}\n".format(library_file))

    # Formatting and converting query files
    all_query_spectra_list = []
    
    fileName, fileExtension = os.path.splitext(os.path.basename(spectrum_file))
    output_filename = ""

    if spectrum_file.find("mzXML") != -1 or spectrum_file.find("mzxml") != -1 or spectrum_file.find("mzML") != -1:
        output_filename = os.path.join(temp_folder, fileName + ".pklbin")
        cmd = "%s %s %s" % (path_to_convert, spectrum_file, output_filename)
        print(cmd)
        os.system(cmd)
    else:
        output_filename = os.path.join(temp_folder, os.path.basename(spectrum_file))
        cmd = "cp %s %s" % (spectrum_file, output_filename)
        os.system(cmd)

    #Input
    faked_output_filename = os.path.join(temp_folder, os.path.basename(spectrum_file))
    all_query_spectra_list.append(faked_output_filename)

    output_parameter_file.write("searchspectra=%s\n" % (" ".join(all_query_spectra_list)))
    output_parameter_file.close()

    cmd = "{} ExecSpectralLibrarySearchMolecular {} -ccms_input_spectradir {} -ccms_results_prefix {} -ll 9".format(path_to_main_exec, \
        parameter_filename, temp_folder, os.path.join(tempresults_folder,"tempresults"))
    print(cmd)
    os.system(cmd)

    #Removing the spectrum
    try:
        os.remove(output_filename)
    except:
        print("Can't remove", output_filename)


def main():
    parser = argparse.ArgumentParser(description='Running library search parallel')
    parser.add_argument('spectrum_file', help='spectrum_file')
    parser.add_argument('library_file', help='library_file')
    parser.add_argument('result_folder', help='output folder for results')
    parser.add_argument('convert_binary', help='conversion binary')
    parser.add_argument('librarysearch_binary', help='librarysearch_binary')

    parser.add_argument('--pm_tolerance', default=0.5, help='pm_tolerance')
    parser.add_argument('--fragment_tolerance', default=0.5, help='fragment_tolerance')
    parser.add_argument('--library_min_cosine', default=0.7, help='library_min_cosine')
    parser.add_argument('--library_min_matched_peaks', default=6, help='library_min_matched_peaks')
    parser.add_argument('--topk', default=1, help='topk')

    parser.add_argument('--analog_search', default=0, help='Turn on analog search, 0 or 1', type=int)

    # This is good for bookkeeping in GNPS2 if you need the full path
    parser.add_argument('--full_relative_query_path', default=None, help='This is the original full relative path of the input file')
    
    parser.add_argument('--filter_precursor', default=1, type=int, help='Filter precursor peaks (1) or not (0)')
    parser.add_argument('--filter_window', default=1, type=int, help='Filter peaks in a window (1) or not (0)')


    args = parser.parse_args()

    temp_folder = "/tmp/gnps_library_search/temp"
    try:
        os.makedirs(temp_folder, exist_ok=True)
    except:
        print("folder error: ", temp_folder)

    tempresults_folder = "/tmp/gnps_library_search/tempresults"
    try:
        os.makedirs(tempresults_folder, exist_ok=True)
    except:
        print("folder error: ", tempresults_folder)

    print(args)

    # performing the search
    search_files(args.spectrum_file, args.library_file,
        temp_folder, tempresults_folder,
        args.convert_binary, args.librarysearch_binary,
        min_cosine=args.library_min_cosine,
        min_matched_peaks=args.library_min_matched_peaks, 
        top_k_results=args.topk, 
        pm_tolerance=args.pm_tolerance,
        ion_tolerance=args.fragment_tolerance,
        analog_search=args.analog_search,
        filter_precursor=args.filter_precursor,
        filter_window=args.filter_window
    )
    
    lib_mgf_basename = os.path.basename(args.library_file)
    
    if not os.path.exists(args.result_folder):
        os.makedirs(args.result_folder, exist_ok=True)
    
    # Reformatting the output
    output_results_file = os.path.join(args.result_folder, os.path.basename(args.spectrum_file) + "_" + lib_mgf_basename + "_gnps.tsv")
    
    if not os.path.exists(os.path.join(tempresults_folder, "tempresults")):
        sys.exit("ERROR: The GNPS library search results file does not exist. Something went wrong. Aborting.")
    results_df = pd.read_csv(os.path.join(tempresults_folder, "tempresults"), sep="\t")

    # Fixing Results, by basename
    results_df["SpectrumFile"] = results_df["SpectrumFile"].apply(lambda x: os.path.basename(x))
    results_df["LibraryName"] = lib_mgf_basename
    # We should rewrite if we have the full path
    if args.full_relative_query_path is not None:
        results_df["SpectrumFile"] = args.full_relative_query_path

        # Create a safe filename from a full path using from args.full_relative_query_path
        safe_filename = args.full_relative_query_path.replace("/", "_").replace(".", "_").replace(" ", "_")

        # lets hash safe_filename to make it shorter
        safe_filename = str(uuid.uuid3(uuid.NAMESPACE_DNS, safe_filename)) + ":" + safe_filename[-20:]

        # fixing the output filename
        output_results_file = os.path.join(args.result_folder, safe_filename + "_" + lib_mgf_basename + "_gnps.tsv")

    results_df.to_csv(output_results_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
