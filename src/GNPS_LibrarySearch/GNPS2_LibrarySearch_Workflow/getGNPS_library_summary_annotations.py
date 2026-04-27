#!/usr/bin/python

# @Cris
# Adapted from https://github.com/Wang-Bioinformatics-Lab/NextflowModules/blob/a31f16297ea5d7f996edb70897179b2e190b9778/bin/library_search/getGNPS_library_annotations.py
# implemented parallel routine to retrieve the identifications info from the GNPS library and adapted to be used to enrich a library summary complete table
# also computing the index of the last/most recent annotation from the retrieve spectrum ID and getting the most updated data (annotations are not always sorted - sorting inplace)
# this script is intended to be executed by the NP3 adm to add the GNPS2 online info to the np3 local library summary

import sys
import os
import pandas as pd
import ming_gnps_library
import requests
from collections import defaultdict
import argparse
import urllib.parse
#from tqdm import tqdm
from functools import lru_cache 

from tqdm.contrib.concurrent import thread_map
#import multiprocessing.dummy as threading

from requests.exceptions import HTTPError

# GLOBALS
# TODO: Use LRU Caching Mechanism
spectrum_id_cache = {}


@lru_cache()
def _get_gnps_library_spectrum(spectrum_id):
    # TODO: We can better sanity check the spectrum_id itself

    try:
        gnps_library_spectrum = ming_gnps_library.get_library_spectrum(spectrum_id)

        #Making sure not an error
        gnps_library_spectrum["annotations"][-1]["Compound_Name"]
    except KeyboardInterrupt:
        raise
    except HTTPError as http_err:
        print(f"Spectrum enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
        gnps_library_spectrum = None
        pass
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
    spectrum_id = output_result_dict["spectrum_id"]

    gnps_library_spectrum = _get_gnps_library_spectrum(spectrum_id)

    if gnps_library_spectrum is None:
        return output_result_dict

    # compute the index of the last/most recent annotation, the annotations are not sorted by date in all entries
    idx_last_annotation = pd.to_datetime([ann["create_time"] for ann in gnps_library_spectrum["annotations"]]).argmax()

    # retrieving the last annotations update, not the oldest and first one in the list
    output_result_dict["Compound_Name"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Compound_Name"].replace("\t", ""))
    output_result_dict["Ion_Source"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Ion_Source"].replace("\t", ""))
    output_result_dict["Instrument"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Instrument"].replace("\t", ""))
    output_result_dict["Compound_Source"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Compound_Source"].replace("\t", ""))
    output_result_dict["PI"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["PI"].replace("\t", ""))
    output_result_dict["Data_Collector"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Data_Collector"].replace("\t", ""))
    output_result_dict["Adduct"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Adduct"].replace("\t", ""))
    output_result_dict["Precursor_MZ"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Precursor_MZ"].replace("\t", ""))
    output_result_dict["ExactMass"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["ExactMass"].replace("\t", ""))
    output_result_dict["Charge"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Charge"].replace("\t", ""))
    output_result_dict["CAS_Number"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["CAS_Number"].replace("\t", ""))
    output_result_dict["Pubmed_ID"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Pubmed_ID"].replace("\t", ""))
    output_result_dict["Smiles"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Smiles"].replace("\t", ""))
    output_result_dict["INCHI"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["INCHI"].replace("\t", ""))
    output_result_dict["INCHI_AUX"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["INCHI_AUX"].replace("\t", ""))
    output_result_dict["Library_Class"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Library_Class"].replace("\t", ""))
    output_result_dict["IonMode"] = (gnps_library_spectrum["annotations"][idx_last_annotation]["Ion_Mode"].replace("\t", ""))

    output_result_dict["Organism"] = (gnps_library_spectrum["spectruminfo"]["library_membership"])
    # set the LibMZ equal to the MGF precursormz
    output_result_dict["LibMZ"] = (output_result_dict["precursormz"])

    if gnps_library_spectrum["annotations"][idx_last_annotation]["Library_Class"] == "1":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-GOLD")
        output_result_dict["LibraryQualityString"] = ("Gold")
    elif gnps_library_spectrum["annotations"][idx_last_annotation]["Library_Class"] == "2":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-SILVER")
        output_result_dict["LibraryQualityString"] = ("Silver")
    elif gnps_library_spectrum["annotations"][idx_last_annotation]["Library_Class"] == "3":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Bronze")
    elif gnps_library_spectrum["annotations"][idx_last_annotation]["Library_Class"] == "4":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Insilico")
    elif gnps_library_spectrum["annotations"][idx_last_annotation]["Library_Class"] == "5":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Insilico")
    elif gnps_library_spectrum["annotations"][idx_last_annotation]["Library_Class"] == "10":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Challenge")
    else:
        print("Invalid Library Class", gnps_library_spectrum["annotations"][idx_last_annotation]["Library_Class"])

    tag_list = [ (tag["tag_desc"] + "[" + tag["tag_type"] + "]") for tag in gnps_library_spectrum["spectrum_tags"]]
    tag_string = "||".join(tag_list).replace("\t", "")

    output_result_dict["tags"] = (tag_string)
    return output_result_dict

def isNA(x):
    return x!=x

# Here we will enrich the smiles, the Smiles column should not be na
# changed columns: "Smiles","INCHI","molecular_formula","InChIKey","InChIKey-Planar","superclass","class","subclass","npclassifier_superclass","npclassifier_class","npclassifier_pathway","library_usi"
def _enrich_annotations(output_result_dict):
    # Calculating inchi
    if isNA(output_result_dict["INCHI"]) or (len(output_result_dict["Smiles"]) > 5 and len(output_result_dict["INCHI"]) < 5):
        try:
            inchi_url = "https://structure.gnps2.org/inchi?smiles={}".format(urllib.parse.quote_plus(output_result_dict["Smiles"]), 
                                urllib.parse.quote_plus(output_result_dict["INCHI"]))
            r = requests.get(inchi_url)
            r.raise_for_status()
            output_result_dict["INCHI"] = r.text
        except HTTPError as http_err:
            print(f"Smiles enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            output_result_dict["INCHI"] = pd.NA
        except:
            output_result_dict["INCHI"] = pd.NA
    
    # Calculating smiles, in case of error keep the smiles and remove the INCHI
    if not isNA(output_result_dict["INCHI"]) and (len(output_result_dict["Smiles"]) < 5 and len(output_result_dict["INCHI"]) > 5):
        try:
            smiles_url = "https://structure.gnps2.org/smiles?inchi={}".format(urllib.parse.quote_plus(output_result_dict["INCHI"]), 
                                urllib.parse.quote_plus(output_result_dict["Smiles"]))
            r = requests.get(smiles_url)
            r.raise_for_status()
            if r.text is not None and r.text != "" and r.text == r.text:
                output_result_dict["Smiles"] = r.text
            else:
                output_result_dict["INCHI"] = pd.NA
        except HTTPError as http_err:
            print(f"Inchi enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            #output_result_dict["Smiles"] = pd.NA
            output_result_dict["INCHI"] = pd.NA
        except:
            #output_result_dict["Smiles"] = pd.NA
            output_result_dict["INCHI"] = pd.NA

    # Calculating molecular formula
    if len(output_result_dict["Smiles"]) > 5:
        try:
            formula_url = "https://structure.gnps2.org/formula?smiles={}".format(output_result_dict["Smiles"])
            r = requests.get(formula_url)
            r.raise_for_status()
            output_result_dict["molecular_formula"] = r.text
        except HTTPError as http_err:
            print(f"Molecular formula enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            output_result_dict["molecular_formula"] = pd.NA
        except:
            output_result_dict["molecular_formula"] = pd.NA
    else:
        output_result_dict["molecular_formula"] = pd.NA
        
    # Calculating inchi key
    if isNA(output_result_dict["INCHI"]) or (len(output_result_dict["Smiles"]) < 5 and len(output_result_dict["INCHI"]) < 5):
        output_result_dict["InChIKey"] = pd.NA
        output_result_dict["InChIKey-Planar"] = pd.NA
    else:
        try:
            inchikey_url = "https://structure.gnps2.org/inchikey?smiles={}&inchi={}".format(urllib.parse.quote_plus(output_result_dict["Smiles"]), 
                                urllib.parse.quote_plus(output_result_dict["INCHI"]))
            r = requests.get(inchikey_url)
            r.raise_for_status()
            output_result_dict["InChIKey"] = r.text
            output_result_dict["InChIKey-Planar"] = r.text.split("-")[0]
        except HTTPError as http_err:
            print(f"InChIKey enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            output_result_dict["InChIKey"] = pd.NA
            output_result_dict["InChIKey-Planar"] = pd.NA
        except:
            output_result_dict["InChIKey"] = pd.NA
            output_result_dict["InChIKey-Planar"] = pd.NA

    # Getting Classyfire "superclass","class","subclass"
    if not isNA(output_result_dict["InChIKey"]) and len(output_result_dict["InChIKey"]) > 5:
        try:
            classyfire_url = "https://classyfire.gnps2.org/entities/{}.json".format(output_result_dict["InChIKey"])
            r = requests.get(classyfire_url, timeout=0.5)
            r.raise_for_status()
            classification_json = r.json()
            output_result_dict["superclass"] = classification_json["superclass"]["name"]
            output_result_dict["class"] = classification_json["class"]["name"]
            output_result_dict["subclass"] = classification_json["subclass"]["name"]
        except HTTPError as http_err:
            print(f"Classyfire enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            output_result_dict["superclass"] = pd.NA
            output_result_dict["class"] = pd.NA
            output_result_dict["subclass"] = pd.NA
        except:
            output_result_dict["superclass"] = pd.NA
            output_result_dict["class"] = pd.NA
            output_result_dict["subclass"] = pd.NA
    else:
        output_result_dict["superclass"] = pd.NA
        output_result_dict["class"] = pd.NA
        output_result_dict["subclass"] = pd.NA

    # Getting NP Classifier "npclassifier_superclass","npclassifier_class","npclassifier_pathway"
    if not isNA(output_result_dict["Smiles"]) and len(output_result_dict["Smiles"]) > 5:
        try:
            npclassifier_url = "https://npclassifier.gnps2.org/classify?smiles={}".format(output_result_dict["Smiles"])
            r = requests.get(npclassifier_url, timeout=10)
            r.raise_for_status()
            classification_json = r.json()

            output_result_dict["npclassifier_superclass"] = "|".join(classification_json["superclass_results"])
            output_result_dict["npclassifier_class"] = "|".join(classification_json["class_results"])
            output_result_dict["npclassifier_pathway"] = "|".join(classification_json["pathway_results"])
        except HTTPError as http_err:
            print(f"NPclassifier enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            output_result_dict["npclassifier_superclass"] = pd.NA
            output_result_dict["npclassifier_class"] = pd.NA
            output_result_dict["npclassifier_pathway"] = pd.NA
        except:
            output_result_dict["npclassifier_superclass"] = pd.NA
            output_result_dict["npclassifier_class"] = pd.NA
            output_result_dict["npclassifier_pathway"] = pd.NA
    else:
        output_result_dict["npclassifier_superclass"] = pd.NA
        output_result_dict["npclassifier_class"] = pd.NA
        output_result_dict["npclassifier_pathway"] = pd.NA

    # Adding a USI
    try:
        output_result_dict["library_usi"] = "mzspec:GNPS:GNPS-LIBRARY:{}".format(output_result_dict["SpectrumID"])
    except:
        output_result_dict["library_usi"] = "No USI"

    return output_result_dict

def enrich_spectrum_annotation_parallel(result_obj):
    # Reading existing data
    spectrum_id = result_obj["spectrum_id"]
    
    # Here we will start to write the output dictionary
    output_result_dict = result_obj
    
    # Here we are going to do the enrichment
    if "CCMSLIB" in str(spectrum_id):
        output_result_dict = _enrich_gnps_annotation(output_result_dict)
    
    return output_result_dict

def enrich_smiles_annotation_parallel(result_dict):
    # Doing further enrichment
    try:
        output_result_dict = _enrich_annotations(result_dict)
    except:
        print("SMILES enrichment failed")
        output_result_dict = result_dict
    return output_result_dict

def enrich_summary_parallell_threads(library_summary_df, output_filename, num_threads=(os.cpu_count() or 1) * 5):
    # clean the lib summary df, remove not used cols
    library_summary_df = library_summary_df.loc[:, ~library_summary_df.columns.isin(['collision_energy', 'instrument',
                                                                                     'ion_source', 'adduct'])]
    #with threading.Pool(processes=num_threads) as pool:
    #    output_list = pool.map(enrich_annotation_parallel, tqdm(library_summary_df.to_dict(orient="records")))
    # enrich parallell with tqdm interface for multiple threads
    n = library_summary_df.shape[0]
    print("Number of spectra in library: {}".format(n))
    output_list = thread_map(enrich_spectrum_annotation_parallel,
                         library_summary_df.to_dict(orient="records"),
                         max_workers=num_threads)
    output_spectrum_annotation = pd.DataFrame(output_list)
    output_spectrum_annotation.to_csv(output_filename.replace(".tsv", "_spectrum_anns.tsv"), sep="\t", index=False)
    m = output_spectrum_annotation.shape[0]
    print("Number of spectra with annotations: {}".format(m))
    # clear some memory space
    del output_list,library_summary_df
    # get list of unique smiles to retrieve their info
    output_unique_smiles_dict = output_spectrum_annotation.loc[(~output_spectrum_annotation.Smiles.duplicated() &
                                                                ~output_spectrum_annotation.Smiles.isna() &
                                                                (output_spectrum_annotation.Smiles != "")),:].to_dict(orient="records")
    output_list = thread_map(enrich_smiles_annotation_parallel,
                             output_unique_smiles_dict,
                             max_workers=min(12,num_threads))
    output_smiles_annotation = pd.DataFrame(output_list)
    del output_unique_smiles_dict,output_list
    output_smiles_annotation.to_csv(output_filename.replace(".tsv", "_smiles_anns.tsv"), sep="\t", index=False)
    # merge the unique smiles annotation to the spectrum annotation
    smiles_ann_cols = ["Smiles", "INCHI", "molecular_formula", "InChIKey", "InChIKey-Planar", "superclass", "class",
                       "subclass", "npclassifier_superclass", "npclassifier_class", "npclassifier_pathway", "library_usi"]
    output_spectrum_annotation = output_spectrum_annotation.loc[:, ~output_spectrum_annotation.columns.isin(smiles_ann_cols[1:])].merge(
        output_smiles_annotation.loc[:, (smiles_ann_cols)], on="Smiles", how='left')
    m = output_spectrum_annotation.shape[0]
    print("Number of spectra with annotations and Smiles info: {}".format(m))
    # Check if final size matches
    if m != n:
        print("WARNING: Number of spectra in the original library does not match with the final annotated table.")
    # output full annotations
    output_spectrum_annotation.to_csv(output_filename, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description='Pulling down GNPS identifications from online API data - needs internet connection. Enrich a library summary table with online info.')
    parser.add_argument("input_librarysummary", type=str, help="Library Summary to be enriched")
    parser.add_argument("output_filename")
    parser.add_argument("--numthreads", default=(os.cpu_count() or 1) * 5, type=int, help="Number of parallel threads to use")
    
    args = parser.parse_args()

    input_library_filename = args.input_librarysummary
    output_result_filename = args.output_filename

    # Trying to load the library summary
    if not os.path.exists(input_library_filename):
        sys.exit("The library summary file does not exists.")
    try:
        library_summary_df = pd.read_csv(input_library_filename, sep="\t", low_memory=False)
    except:
        library_summary_df = None
        sys.exit("The library summary could not be loaded, it is not a valid tsv file.")
    
    enrich_summary_parallell_threads(library_summary_df, output_result_filename, num_threads=args.numthreads)

if __name__ == "__main__":
    main()
