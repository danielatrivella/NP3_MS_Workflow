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


def _enrich_gnps_annotation(output_result_dict):
    spectrum_id = output_result_dict["SpectrumID"]

    gnps_library_spectrum = _get_gnps_library_spectrum(spectrum_id)

    if gnps_library_spectrum is None:
        return output_result_dict

    output_result_dict["Compound_Name"] = (gnps_library_spectrum["annotations"][0]["Compound_Name"].replace("\t", ""))
    output_result_dict["Ion_Source"] = (gnps_library_spectrum["annotations"][0]["Ion_Source"].replace("\t", ""))
    output_result_dict["Instrument"] = (gnps_library_spectrum["annotations"][0]["Instrument"].replace("\t", ""))
    output_result_dict["Compound_Source"] = (gnps_library_spectrum["annotations"][0]["Compound_Source"].replace("\t", ""))
    output_result_dict["PI"] = (gnps_library_spectrum["annotations"][0]["PI"].replace("\t", ""))
    output_result_dict["Data_Collector"] = (gnps_library_spectrum["annotations"][0]["Data_Collector"].replace("\t", ""))
    output_result_dict["Adduct"] = (gnps_library_spectrum["annotations"][0]["Adduct"].replace("\t", ""))
    output_result_dict["Precursor_MZ"] = (gnps_library_spectrum["annotations"][0]["Precursor_MZ"].replace("\t", ""))
    output_result_dict["ExactMass"] = (gnps_library_spectrum["annotations"][0]["ExactMass"].replace("\t", ""))
    output_result_dict["Charge"] = (gnps_library_spectrum["annotations"][0]["Charge"].replace("\t", ""))
    output_result_dict["CAS_Number"] = (gnps_library_spectrum["annotations"][0]["CAS_Number"].replace("\t", ""))
    output_result_dict["Pubmed_ID"] = (gnps_library_spectrum["annotations"][0]["Pubmed_ID"].replace("\t", ""))
    output_result_dict["Smiles"] = (gnps_library_spectrum["annotations"][0]["Smiles"].replace("\t", ""))
    output_result_dict["INCHI"] = (gnps_library_spectrum["annotations"][0]["INCHI"].replace("\t", ""))
    output_result_dict["INCHI_AUX"] = (gnps_library_spectrum["annotations"][0]["INCHI_AUX"].replace("\t", ""))
    output_result_dict["Library_Class"] = (gnps_library_spectrum["annotations"][0]["Library_Class"].replace("\t", ""))
    output_result_dict["IonMode"] = (gnps_library_spectrum["annotations"][0]["Ion_Mode"].replace("\t", ""))

    output_result_dict["Organism"] = (gnps_library_spectrum["spectruminfo"]["library_membership"])
    output_result_dict["LibMZ"] = (gnps_library_spectrum["annotations"][0]["Precursor_MZ"])

    if gnps_library_spectrum["annotations"][0]["Library_Class"] == "1":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-GOLD")
        output_result_dict["LibraryQualityString"] = ("Gold")
    elif gnps_library_spectrum["annotations"][0]["Library_Class"] == "2":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-SILVER")
        output_result_dict["LibraryQualityString"] = ("Silver")
    elif gnps_library_spectrum["annotations"][0]["Library_Class"] == "3":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Bronze")
    elif gnps_library_spectrum["annotations"][0]["Library_Class"] == "4":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Insilico")
    elif gnps_library_spectrum["annotations"][0]["Library_Class"] == "5":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Insilico")
    elif gnps_library_spectrum["annotations"][0]["Library_Class"] == "10":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Challenge")
    else:
        print("Invalid Library Class", gnps_library_spectrum["annotations"][0]["Library_Class"])


    tag_list = [ (tag["tag_desc"] + "[" + tag["tag_type"] + "]") for tag in gnps_library_spectrum["spectrum_tags"]]
    tag_string = "||".join(tag_list).replace("\t", "")

    output_result_dict["tags"] = (tag_string)

    return output_result_dict

def _enrich_librarysummary_annotations(output_result_dict, library_dict=None):
    spectrum_id = output_result_dict["SpectrumID"]

    if library_dict is not None and spectrum_id in library_dict:
        library_spectrum = library_dict[spectrum_id]
    else:
        library_spectrum = {}

    output_result_dict["Compound_Name"] = str(library_spectrum.get("compound_name", library_spectrum.get("NAME", ""))).replace("\t", "")
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


# Here we will enrich the smiles
def _enrich_annotations(output_result_dict):
    # Calculating inchi
    if len(output_result_dict["Smiles"]) > 5 and len(output_result_dict["INCHI"]) < 5:
        try:
            inchi_url = "https://structure.gnps2.org/inchi?smiles={}".format(urllib.parse.quote_plus(output_result_dict["Smiles"]), 
                                urllib.parse.quote_plus(output_result_dict["INCHI"]))
            r = requests.get(inchi_url)
            r.raise_for_status()
            output_result_dict["INCHI"] = r.text
        except:
            output_result_dict["INCHI"] = "N/A"
    
    # Calculating smiles
    if len(output_result_dict["Smiles"]) < 5 and len(output_result_dict["INCHI"]) > 5:
        try:
            smiles_url = "https://structure.gnps2.org/smiles?inchi={}".format(urllib.parse.quote_plus(output_result_dict["INCHI"]), 
                                urllib.parse.quote_plus(output_result_dict["Smiles"]))
            r = requests.get(smiles_url)
            r.raise_for_status()
            output_result_dict["Smiles"] = r.text
        except:
            output_result_dict["Smiles"] = "N/A"

    # Calculating molecular formula
    if len(output_result_dict["Smiles"]) > 5:
        try:
            formula_url = "https://structure.gnps2.org/formula?smiles={}".format(output_result_dict["Smiles"])
            r = requests.get(formula_url)
            r.raise_for_status()
            output_result_dict["molecular_formula"] = r.text
        except:
            output_result_dict["molecular_formula"] = "N/A"
    else:
        output_result_dict["molecular_formula"] = "N/A"
        
    # Calculating inchi key
    if len(output_result_dict["Smiles"]) < 5 and len(output_result_dict["INCHI"]) < 5:
        output_result_dict["InChIKey"] = "N/A"
        output_result_dict["InChIKey-Planar"] = "N/A"
    else:
        try:
            inchikey_url = "https://structure.gnps2.org/inchikey?smiles={}&inchi={}".format(urllib.parse.quote_plus(output_result_dict["Smiles"]), 
                                urllib.parse.quote_plus(output_result_dict["INCHI"]))
            r = requests.get(inchikey_url)
            r.raise_for_status()
            output_result_dict["InChIKey"] = r.text
            output_result_dict["InChIKey-Planar"] = r.text.split("-")[0]
        except:
            output_result_dict["InChIKey"] = "N/A"
            output_result_dict["InChIKey-Planar"] = "N/A"

    # Getting Classyfire
    if len(output_result_dict["InChIKey"]) > 5:
        try:
            classyfire_url = "https://classyfire.gnps2.org/entities/{}.json".format(output_result_dict["InChIKey"])
            r = requests.get(classyfire_url, timeout=0.5)
            r.raise_for_status()
            classification_json = r.json()
            output_result_dict["superclass"] = classification_json["superclass"]["name"]
            output_result_dict["class"] = classification_json["class"]["name"]
            output_result_dict["subclass"] = classification_json["subclass"]["name"]
        except:
            output_result_dict["superclass"] = "N/A"
            output_result_dict["class"] = "N/A"
            output_result_dict["subclass"] = "N/A"
    else:
        output_result_dict["superclass"] = "N/A"
        output_result_dict["class"] = "N/A"
        output_result_dict["subclass"] = "N/A"

    # Getting NP Classifier
    if len(output_result_dict["Smiles"]) > 5:
        try:
            npclassifier_url = "https://npclassifier.gnps2.org/classify?smiles={}".format(output_result_dict["Smiles"])
            r = requests.get(npclassifier_url, timeout=10)
            r.raise_for_status()
            classification_json = r.json()

            output_result_dict["npclassifier_superclass"] = "|".join(classification_json["superclass_results"])
            output_result_dict["npclassifier_class"] = "|".join(classification_json["class_results"])
            output_result_dict["npclassifier_pathway"] = "|".join(classification_json["pathway_results"])
        except:
            output_result_dict["npclassifier_superclass"] = "N/A"
            output_result_dict["npclassifier_class"] = "N/A"
            output_result_dict["npclassifier_pathway"] = "N/A"
    else:
        output_result_dict["npclassifier_superclass"] = "N/A"
        output_result_dict["npclassifier_class"] = "N/A"
        output_result_dict["npclassifier_pathway"] = "N/A"

    # Adding a USI
    try:
        output_result_dict["library_usi"] = "mzspec:GNPS:GNPS-LIBRARY:{}".format(output_result_dict["SpectrumID"])
    except:
        output_result_dict["library_usi"] = "No USI"

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
        #else:
        #    if "CCMSLIB" in str(spectrum_id):
        #        output_result_dict = _enrich_gnps_annotation(output_result_dict)
        #    else:
        #        # checking if in library summary
        #        if library_summary_df is not None:
        #            output_result_dict = _enrich_librarysummary_annotations(output_result_dict, library_dict=library_dict)
            # Doing further enrichment
        #    try:
        #        output_result_dict = _enrich_annotations(output_result_dict)
        #    except:
        #        pass
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
