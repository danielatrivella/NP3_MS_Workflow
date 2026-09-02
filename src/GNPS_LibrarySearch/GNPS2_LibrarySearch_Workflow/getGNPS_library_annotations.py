#!/usr/bin/python

# @Cris
# Adapted from https://github.com/Wang-Bioinformatics-Lab/NextflowModules/blob/a31f16297ea5d7f996edb70897179b2e190b9778/bin/library_search/getGNPS_library_annotations.py
# This is intended to be executed only offline and to retrieve all the data from an enriched library summary

import sys
import os
import pandas as pd
import ming_gnps_library
import requests
from collections import defaultdict
import argparse
import urllib.parse
from tqdm import tqdm
from functools import lru_cache 

# GLOBALS
# TODO: Use LRU Caching Mechanism
spectrum_id_cache = {}


@lru_cache()
def _get_gnps_library_spectrum(spectrum_id):
    # TODO: We can better sanity check the spectrum_id itself

    try:
        gnps_library_spectrum = ming_gnps_library.get_library_spectrum(spectrum_id)

        #Making sure not an error
        gnps_library_spectrum["annotations"][0]["Compound_Name"]
    except KeyboardInterrupt:
        raise
    except:
        gnps_library_spectrum = None
        pass

    return gnps_library_spectrum

def _prep_library_dict(library_summary_df):
    library_dict = {}
    for row in library_summary_df.to_dict(orient="records"):
        library_dict[row["spectrum_id"]] = row
    return library_dict


def _enrich_librarysummary_annotations(output_result_dict, library_dict=None):
    spectrum_id = output_result_dict["SpectrumID"]

    if library_dict is not None and spectrum_id in library_dict:
        library_spectrum = library_dict[spectrum_id]
    else:
        library_spectrum = {}

    output_result_dict["Compound_Name"] = str(library_spectrum.get("Compound_Name", library_spectrum.get("NAME", ""))).replace("\t", "").replace("\"", "")
    output_result_dict["Ion_Source"] = str(library_spectrum.get("Ion_Source", "")).replace("\t", "")
    output_result_dict["Instrument"] = str(library_spectrum.get("Instrument", "")).replace("\t", "")
    output_result_dict["LibMZ"] = library_spectrum.get("precursormz", "")
    output_result_dict["Adduct"] = str(library_spectrum.get("Adduct", "")).replace("\t", "")
    output_result_dict["Charge"] = str(library_spectrum.get("Charge", "")).replace("\t", "")
    output_result_dict["Smiles"] = str(library_spectrum.get("Smiles", "")).replace("\t", "")
    output_result_dict["LibScan"] = str(library_spectrum.get("scan", "")).replace("\t", "")
    
    #  check if the library summary contains additional GNPS annotations (previous enriched online) and
    #  aggregate the info, otherwise leave empty
    output_result_dict["INCHI"] = (str(library_spectrum.get("INCHI")).replace("\t", "") if library_spectrum.get("INCHI") is not None else "N/A")
    output_result_dict["INCHI_AUX"] = (str(library_spectrum.get("INCHI_AUX")).replace("\t", "") if library_spectrum.get("INCHI_AUX") is not None else "N/A")
    output_result_dict["Library_Class"] = (str(library_spectrum.get("Library_Class")).replace("\t", "") if library_spectrum.get("Library_Class") is not None else "N/A")
    output_result_dict["tags"] = (str(library_spectrum.get("tags")).replace("\t", "") if library_spectrum.get("tags") is not None else "N/A")
    output_result_dict["IonMode"] = (str(library_spectrum.get("IonMode")).replace("\t", "") if library_spectrum.get("IonMode") is not None else "N/A")
    output_result_dict["PI"] = (str(library_spectrum.get("PI")).replace("\t", "") if library_spectrum.get("PI") is not None else "N/A")
    output_result_dict["Data_Collector"] = (str(library_spectrum.get("Data_Collector")).replace("\t", "") if library_spectrum.get("Data_Collector") is not None else "N/A")
    output_result_dict["Precursor_MZ"] = (str(library_spectrum.get("Precursor_MZ")).replace("\t", "") if library_spectrum.get("Precursor_MZ") is not None else "N/A")
    output_result_dict["ExactMass"] = (str(library_spectrum.get("ExactMass")).replace("\t", "") if library_spectrum.get("ExactMass") is not None else "N/A")
    output_result_dict["CAS_Number"] = (str(library_spectrum.get("CAS_Number")).replace("\t", "") if library_spectrum.get("CAS_Number") is not None else "N/A")
    output_result_dict["Pubmed_ID"] = (str(library_spectrum.get("Pubmed_ID")).replace("\t", "") if library_spectrum.get("Pubmed_ID") is not None else "N/A")
    output_result_dict["Organism"] = (str(library_spectrum.get("Organism")).replace("\t", "") if library_spectrum.get("Organism") is not None else "N/A")
    output_result_dict["Compound_Source"] = (str(library_spectrum.get("Compound_Source")).replace("\t", "") if library_spectrum.get("Compound_Source") is not None else "N/A")
    
    output_result_dict["UpdateWorkflowName"] = (str(library_spectrum.get("UpdateWorkflowName")).replace("\t", "") if library_spectrum.get("UpdateWorkflowName") is not None else "N/A")
    output_result_dict["LibraryQualityString"] = (str(library_spectrum.get("LibraryQualityString")).replace("\t", "") if library_spectrum.get("LibraryQualityString") is not None else "N/A")
    
    # get molecular formula, np classyfire and np classifier if present
    output_result_dict["molecular_formula"] = (str(library_spectrum.get("molecular_formula")).replace("\t", "") if library_spectrum.get("molecular_formula") is not None else "N/A")
    output_result_dict["InChIKey"] = (str(library_spectrum.get("InChIKey")).replace("\t", "") if library_spectrum.get("InChIKey") is not None else "N/A")
    output_result_dict["InChIKey-Planar"] = (str(library_spectrum.get("InChIKey-Planar")).replace("\t", "") if library_spectrum.get("InChIKey-Planar") is not None else "N/A")
    output_result_dict["superclass"] = (str(library_spectrum.get("superclass")).replace("\t", "") if library_spectrum.get("superclass") is not None else "N/A")
    output_result_dict["class"] = (str(library_spectrum.get("class")).replace("\t", "") if library_spectrum.get("class") is not None else "N/A")
    output_result_dict["subclass"] = (str(library_spectrum.get("subclass")).replace("\t", "") if library_spectrum.get("subclass") is not None else "N/A")
    output_result_dict["npclassifier_superclass"] = (str(library_spectrum.get("npclassifier_superclass")).replace("\t", "") if library_spectrum.get("npclassifier_superclass") is not None else "N/A")
    output_result_dict["npclassifier_class"] = (str(library_spectrum.get("npclassifier_class")).replace("\t", "") if library_spectrum.get("npclassifier_class") is not None else "N/A")
    output_result_dict["npclassifier_pathway"] = (str(library_spectrum.get("npclassifier_pathway")).replace("\t", "") if library_spectrum.get("npclassifier_pathway") is not None else "N/A")
    
    # checking all the values, and change to N/A if empty
    for key in output_result_dict:
        if output_result_dict[key] == "" or output_result_dict[key] == "nan":
            output_result_dict[key] = "N/A"

    return output_result_dict


def enrich_output(input_filename, output_filename, 
                    topk=None, 
                    library_summary_df=None, 
                    filtertostructures=False):
    library_dict = {}
    if library_summary_df is not None:
        try:
            library_dict = _prep_library_dict(library_summary_df)
        except:
            pass

    if not os.path.exists(input_filename):
        #open(output_filename, "w").close()
        sys.exit("ERROR: Input file does not exist: "+input_filename)
    
    try:
        input_results_df = pd.read_csv(input_filename, sep="\t", low_memory=False)
    except pd.errors.EmptyDataError:
        print("WARNING: Input file is empty, there was no result in the library search. Skipping the annotation.\n")
        sys.exit()
    except:
        #open(output_filename, "w").close()
        sys.exit("ERROR: Input file is not a valid tsv file")


    # Here we will try to filter to topk
    if topk is not None:
        try:
            input_results_df["MQScore"] = input_results_df["MQScore"].astype(float)
            input_results_df = input_results_df.sort_values("MQScore", ascending=False)
            input_results_df = input_results_df.groupby('FileScanUniqueID').head(topk).reset_index(drop=True) 
        except:
            pass

    # Counting number of hits per filename
    number_hits_per_query = defaultdict(lambda: 0)
    for result_obj in input_results_df.to_dict(orient="records"):
        number_hits_per_query[result_obj["FileScanUniqueID"]] += 1

    output_list = []
    for result_obj in tqdm(input_results_df.to_dict(orient="records")):
        # Reading exsting data
        spectrum_id = result_obj["LibrarySpectrumID"]
        score = result_obj["MQScore"]
        filename = result_obj["SpectrumFile"]
        libfilename = result_obj["LibraryName"]
        scan = result_obj["#Scan#"]
        TIC_Query = result_obj["UnstrictEvelopeScore"]
        RT_Query = result_obj["p-value"]
        SpecCharge = result_obj["Charge"]
        SpecMZ = result_obj["SpecMZ"]
        MZErrorPPM = result_obj["mzErrorPPM"]
        SharedPeaks = result_obj["LibSearchSharedPeaks"]
        MassDiff = result_obj["ParentMassDiff"]

        # Here we will start to write the output dictionary
        output_result_dict = {}

        output_result_dict["SpectrumID"] = (spectrum_id)
        output_result_dict["#Scan#"] = (scan)
        output_result_dict["SpectrumFile"] = (filename)
        output_result_dict["LibraryName"] = (libfilename)
        output_result_dict["MQScore"] = (score)
        output_result_dict["TIC_Query"] = (TIC_Query)
        output_result_dict["RT_Query"] = (RT_Query)
        output_result_dict["MZErrorPPM"] = (MZErrorPPM)
        output_result_dict["SharedPeaks"] = (SharedPeaks)
        output_result_dict["MassDiff"] = (MassDiff)
        
        output_result_dict["SpecMZ"] = (SpecMZ)
        output_result_dict["SpecCharge"] = (SpecCharge)
        output_result_dict["FileScanUniqueID"] = (result_obj["FileScanUniqueID"])
        output_result_dict["NumberHits"] = (number_hits_per_query[result_obj["FileScanUniqueID"]])

        # Here we are going to do the enrichment
        # always with forceoffline, online enrichment is disabled in the np3 flow
        output_result_dict = _enrich_librarysummary_annotations(output_result_dict, library_dict=library_dict)
        
        output_list.append(output_result_dict)

    # Here we can filter based upon the structure criteria
    if filtertostructures is True:
        # Filtering only if the length of Smiles and InchI are small
        output_list = [x for x in output_list if len(x["Smiles"]) > 5 or len(x["INCHI"]) > 5]

    pd.DataFrame(output_list).to_csv(output_filename, sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser(description='Pulling down GNPS identifications from the library summary with annotations. Only offline enrichment.')
    parser.add_argument("input_filename")
    parser.add_argument("output_filename")
    parser.add_argument("--topk", default=None, type=int, help="Top K results per query, default no filter")
    parser.add_argument("--librarysummary", default=None, type=str, help="Library Summary, importnat for non-GNPS libraries")
    parser.add_argument("--filtertostructures", default="0", type=str, help="Filter to structures only if 1")
    #parser.add_argument("--forceoffline", default="No", type=str, help="Force offline so we can avoid using a bunch of API calls, Yes or No")

    args = parser.parse_args()

    input_result_filename = args.input_filename
    output_result_filename = args.output_filename

    # Trying to load the library summary
    try:
        library_summary_df = pd.read_csv(args.librarysummary, sep="\t", low_memory=False)
    except:
        library_summary_df = None

    enrich_output(input_result_filename, output_result_filename, topk=args.topk, 
                    library_summary_df=library_summary_df,
                    filtertostructures=(args.filtertostructures == "1"))

if __name__ == "__main__":
    main()
