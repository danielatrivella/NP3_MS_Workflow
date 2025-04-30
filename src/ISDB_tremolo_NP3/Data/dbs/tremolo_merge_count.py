from __future__ import print_function

import csv
import sys
import pandas as pd
import numpy as np

def map_and_index(arr1, arr2):
    """
    Maps elements from arr1 to their indices in arr2,
    returning an array of indices with the order of arr2.

    Args:
        arr1 (np.ndarray): The array whose elements are to be mapped.
        arr2 (np.ndarray): The array defining the desired order of indices.

    Returns:
        np.ndarray: An array of indices, where each element at index `i`
                    corresponds to the index of arr2[i] in arr1.
    """
    mapping = {val: idx for idx, val in enumerate(arr1)}
    return np.array([mapping[val] for val in arr2])

csv.field_size_limit(sys.maxsize)

if len(sys.argv) < 3:
    print(" Incorrect number of arguments")
    print(" 1. result file; 2. max results to merge and 3. and 4. count output files")
    sys.exit()

result_file = sys.argv[1]
max_results = int(sys.argv[2])

with open(result_file, 'rt') as f:
    reader = csv.reader(f, delimiter=',')
    header = False
    result_smile = {}
    result_id = {}
    result_cn = {}
    result_inchikey = {}
    result_score = {}
    result_ppmError = {}
    result_sharedPeaks = {}
    result_mf = {}
    result_mw = {}
    result_cas = {}
    scan_count = {}
    # NPClassifier columns result
    result_npc_superclass = {}
    result_npc_class = {}
    result_npc_pathway = {}
    result_npc_isglycoside = {}
    result_cf_subclass = {}
    # NPAtlas columns result
    result_npa_id = {}
    result_npa_compound_names = {}
    result_npa_origin_type = {}
    result_npa_genus = {}
    result_npa_origin_species = {}
    for row in reader:
        if not header:
            scan_id_pos = row.index('msclusterID')
            unpd_id_pos = row.index('CompoundName')
            smiles_id_pos = row.index('SMILES')
            inchikey_id_pos = row.index('InChIKey')
            cas_id_pos = row.index('CAS')
            cn_id_pos = row.index('chemicalNames')
            mf_id_pos = row.index('molecularFormula')
            mw_id_pos = row.index('molecularWeight')
            score_id_pos = row.index('MQScore')
            ppmError_id_pos = row.index('mzErrorPPM')
            sharPeaks_id_pos = row.index('LibSearchSharedPeaks')
            # NPC columns index
            npc_superclass_pos = row.index('NPClassifier_class')
            npc_class_pos = row.index('NPClassifier_superclass')
            npc_pathway_pos = row.index('NPClassifier_pathway')
            npc_isglycoside_pos = row.index('NPClassifier_isglycoside')
            cf_subclass_pos = row.index('ClassyFire_subclass')
            # NPA columns index
            npa_id_pos = row.index('NPAtlas_id')
            npa_compound_names_pos = row.index('NPAtlas_compound_names')
            npa_origin_type_pos = row.index('NPAtlas_origin_type')
            npa_genus_pos = row.index('NPAtlas_genus')
            npa_origin_species_pos = row.index('NPAtlas_origin_species')
            header = True
        else:
            # obtain scan results
            if row[scan_id_pos] in result_smile:
                # Handle multiple molecules by node until max
                if scan_count[row[scan_id_pos]] >= max_results:
                    continue
                result_smile[row[scan_id_pos]] += "," + row[smiles_id_pos]
                result_cn[row[scan_id_pos]] += ";" + row[cn_id_pos].replace(";", ":")
                result_inchikey[row[scan_id_pos]] += ";" + row[inchikey_id_pos]
                result_mf[row[scan_id_pos]] += ";" + row[mf_id_pos]
                result_mw[row[scan_id_pos]] += ";" + row[mw_id_pos]
                result_cas[row[scan_id_pos]] += ";" + row[cas_id_pos]
                result_id[row[scan_id_pos]] += ";" + row[unpd_id_pos]
                result_score[row[scan_id_pos]] += ";" + row[score_id_pos]
                result_ppmError[row[scan_id_pos]] += ";" + row[ppmError_id_pos]
                result_sharedPeaks[row[scan_id_pos]] += ";" + row[sharPeaks_id_pos]
                result_npc_superclass[row[scan_id_pos]] += ";" + row[npc_superclass_pos]
                result_npc_class[row[scan_id_pos]] += ";" + row[npc_class_pos]
                result_npc_pathway[row[scan_id_pos]] += ";" + row[npc_pathway_pos]
                result_npc_isglycoside[row[scan_id_pos]] += ";" + row[npc_isglycoside_pos]
                result_cf_subclass[row[scan_id_pos]] += ";" + row[cf_subclass_pos]
                result_npa_id[row[scan_id_pos]] += ";" + row[npa_id_pos]
                result_npa_compound_names[row[scan_id_pos]] += ";" + row[npa_compound_names_pos]
                result_npa_origin_type[row[scan_id_pos]] += ";" + row[npa_origin_type_pos]
                result_npa_genus[row[scan_id_pos]] += ";" + row[npa_genus_pos]
                result_npa_origin_species[row[scan_id_pos]] += ";" + row[npa_origin_species_pos]
                scan_count[row[scan_id_pos]] += 1
            else:
                result_smile[row[scan_id_pos]] = row[smiles_id_pos]
                result_cn[row[scan_id_pos]] = row[cn_id_pos].replace(";", ":")
                result_inchikey[row[scan_id_pos]] = row[inchikey_id_pos]
                result_mf[row[scan_id_pos]] = row[mf_id_pos]
                result_mw[row[scan_id_pos]] = row[mw_id_pos]
                result_cas[row[scan_id_pos]] = row[cas_id_pos]
                result_id[row[scan_id_pos]] = row[unpd_id_pos]
                result_score[row[scan_id_pos]] = row[score_id_pos]
                result_ppmError[row[scan_id_pos]] = row[ppmError_id_pos]
                result_sharedPeaks[row[scan_id_pos]] = row[sharPeaks_id_pos]
                result_npc_superclass[row[scan_id_pos]] = row[npc_superclass_pos]
                result_npc_class[row[scan_id_pos]] = row[npc_class_pos]
                result_npc_pathway[row[scan_id_pos]] = row[npc_pathway_pos]
                result_npc_isglycoside[row[scan_id_pos]] = row[npc_isglycoside_pos]
                result_cf_subclass[row[scan_id_pos]] = row[cf_subclass_pos]
                result_npa_id[row[scan_id_pos]] = row[npa_id_pos]
                result_npa_compound_names[row[scan_id_pos]] = row[npa_compound_names_pos]
                result_npa_origin_type[row[scan_id_pos]] = row[npa_origin_type_pos]
                result_npa_genus[row[scan_id_pos]] = row[npa_genus_pos]
                result_npa_origin_species[row[scan_id_pos]] = row[npa_origin_species_pos]
                scan_count[row[scan_id_pos]] = 1


print("  Merging the tremolo result to the count tables files")
for count_file in sys.argv[3:]:
    # make the merging with pandas, replace previous identification result
    count_table = pd.read_csv(count_file, low_memory=False)
    # convert the msclusterIDs to str to allow mapping them with the keys of the dictionaries of the tremolo data
    msclusterIDs_mapping = map_and_index(count_table.msclusterID.astype(str).values, list(result_id.keys()))

    if "tremolo_UNPD_IDs" in count_table:
        print("    - Overwriting previous tremolo identification result")
        count_table = count_table.drop(["tremolo_UNPD_IDs","tremolo_SMILES","tremolo_chemicalNames","tremolo_molecularFormula",
                          "tremolo_molecularWeight","tremolo_CAS","tremolo_MQScore","tremolo_mzErrorPPM",
                          "tremolo_numSharedPeaks","tremolo_NPClassifier_superclass","tremolo_NPClassifier_class",
                          "tremolo_ClassyFire_subclass","tremolo_NPClassifier_pathway",
                          "tremolo_NPClassifier_isglycoside","tremolo_NPAtlas_id","tremolo_NPAtlas_compound_names",
                          "tremolo_NPAtlas_origin_type","tremolo_NPAtlas_genus","tremolo_NPAtlas_origin_species"],
                         axis=1, errors='ignore')
    # set the tremolo values to the count_table columns
    count_table.loc[msclusterIDs_mapping,"tremolo_UNPD_IDs"] = list(result_id.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_SMILES"] = list(result_smile.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_chemicalNames"] = list(result_cn.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_InChIKey"] = list(result_inchikey.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_molecularFormula"] = list(result_mf.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_molecularWeight"] = list(result_mw.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_CAS"] = list(result_cas.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_MQScore"] = list(result_score.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_mzErrorPPM"] = list(result_ppmError.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_numSharedPeaks"] = list(result_sharedPeaks.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_NPClassifier_superclass"] = list(result_npc_superclass.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_NPClassifier_class"] = list(result_npc_class.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_ClassyFire_subclass"] = list(result_cf_subclass.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_NPClassifier_pathway"] = list(result_npc_pathway.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_NPClassifier_isglycoside"] = list(result_npc_isglycoside.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_NPAtlas_id"] = list(result_npa_id.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_NPAtlas_compound_names"] = list(result_npa_compound_names.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_NPAtlas_origin_type"] = list(result_npa_origin_type.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_NPAtlas_genus"] = list(result_npa_genus.values())
    count_table.loc[msclusterIDs_mapping,"tremolo_NPAtlas_origin_species"] = list(result_npa_origin_species.values())
    # write output file
    print("    - Outputting ", count_file)
    count_table.to_csv(count_file, index=False)
