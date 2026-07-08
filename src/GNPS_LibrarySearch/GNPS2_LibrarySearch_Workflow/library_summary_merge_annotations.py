import argparse
import pandas as pd
import os, sys
import numpy as np

from rdkit import Chem, RDLogger
from rdkit import RDLogger
# Disable RDKit logging msgs
RDLogger.DisableLog('rdApp.*')

# validate smiles string, MolFromSmiles returns None if the string is invalid
def is_valid_smiles(smiles_string):
    if pd.isna(smiles_string):
        return False
    mol = Chem.MolFromSmiles(smiles_string)
    return mol is not None


def join_library_summary_annotations(library_summary_file, library_summary_annotations_file, output_filename):
    # Trying to load the library summary
    if not os.path.exists(library_summary_file):
        sys.exit("The library summary file does not exists.")
    try:
        library_summary = pd.read_csv(library_summary_file, sep="\t", low_memory=False)
    except:
        sys.exit("The library summary could not be loaded, it is not a valid tsv file.")
    # Trying to load the library summary with annotations, as provided by the NP3 repository
    if not os.path.exists(library_summary_annotations_file):
        sys.exit("The library summary with annotations file does not exists.")
    try:
        library_summary_annotations = pd.read_csv(library_summary_annotations_file, sep="\t", low_memory=False)
    except:
        sys.exit("The library summary with annotations could not be loaded, it is not a valid tsv file.")
    
    # get the annotations columns that are not in the library summary original file
    annotation_cols_id = np.concatenate([["spectrum_id"],
                                        library_summary_annotations.columns[
                                            ~library_summary_annotations.columns.isin(library_summary.columns)].values])
    
    # check if all spectrum_id of the library summary are present in the annotations
    if not library_summary.spectrum_id.isin(library_summary_annotations.spectrum_id).all():
        print("WARNING: There are new Spectrum IDs in the library summary which do not appear in the annotations file.",
              "A total of", str(sum(~library_summary.spectrum_id.isin(library_summary_annotations.spectrum_id))),
              "Spectrum IDs are new in the downloaded library and will be missing the GNPS enriched annotations.\n",
              "Please contact the dev team (adm) to update the library summary annotations.")
    
    # merge the library summary with the annotations columns
    library_summary_enriched = library_summary.merge(library_summary_annotations.loc[:, annotation_cols_id],
                                                     on='spectrum_id', how='left')

    # check consistency of merge, no spectrum may be lost here
    if library_summary_enriched.shape[0] != library_summary.shape[0]:
        sys.exit("The library summary enrichment with the summary annotations failed, some Spectrum IDs went missing after the match.")
    
    # validate SMILES list, remove invalid SMILES from the "Smiles" column and then
    # # for the invalid ones check if there is a valid Smiles string in the "smiles" column
    invalid_smiles = 0
    print("- Validating the Smiles string and removing invalid SMILES...")
    check_Smiles = pd.DataFrame({'Smiles':library_summary_enriched.Smiles.unique(), 'Smiles_valid':np.nan})
    for i in range(check_Smiles.shape[0]):
        check_Smiles_is_valid = is_valid_smiles(check_Smiles.loc[i, "Smiles"])
        if not check_Smiles_is_valid:
            check_Smiles.loc[i, "Smiles_valid"] = np.nan
            invalid_smiles += 1
        else:
            check_Smiles.loc[i, "Smiles_valid"] = check_Smiles.loc[i, "Smiles"]
    # merge the checked smiles with the library summary
    library_summary_enriched = library_summary_enriched.merge(check_Smiles, on='Smiles', how='left')
    # remove invalid Smiles
    library_summary_enriched.loc[:,'Smiles'] = library_summary_enriched.loc[:,'Smiles_valid']
    library_summary_enriched.drop(['Smiles_valid'], axis=1, inplace=True)
    print("  - Removed a total of", invalid_smiles-1, " unique invalid SMILES string.") # remove one from NA
    # for the invalid Smiles, check if there is a valid SMILES in the "smiles" column
    # this is a second try for the invalid Smiles
    check_Smiles = pd.DataFrame({'smiles': library_summary_enriched.loc[(library_summary_enriched.Smiles.isna() & ~library_summary_enriched.smiles.isna()),
                                                                        'smiles'].unique(),
                                 'smiles_valid': np.nan})
    for i in range(check_Smiles.shape[0]):
        check_Smiles_is_valid = is_valid_smiles(check_Smiles.loc[i,"smiles"])
        if check_Smiles_is_valid:
            check_Smiles.loc[i,"smiles_valid"] = check_Smiles.loc[i,"smiles"]
    # if any smiles from the invalid ones was validated, replace them in the Smiles column where it had NAs
    if (~check_Smiles.smiles_valid.isna()).any():
        library_summary_enriched = library_summary_enriched.merge(check_Smiles, on='smiles', how='left')
        check_Smiles_is_valid = (~library_summary_enriched.smiles_valid.isna() & library_summary_enriched.Smiles.isna())
        library_summary_enriched.loc[check_Smiles_is_valid, 'Smiles'] = library_summary_enriched.loc[check_Smiles_is_valid, 'smiles_valid']
        library_summary_enriched.drop(['smiles_valid'], axis=1, inplace=True)
    
    revalidated_smiles = library_summary_enriched.loc[check_Smiles_is_valid, 'Smiles'].unique().size
    print("  - Recovered a total of ",revalidated_smiles, " unique SMILES string.")
    print("  - Total number of spectra present in the GNPS library summary =",library_summary_enriched.shape[0])
    print("  - Total number of spectra with a valid SMILES =", (~library_summary_enriched.Smiles.isna()).sum(),
          "(",round((~library_summary_enriched.Smiles.isna()).sum()/library_summary_enriched.shape[0]*100),"%)")
    print("  - Total number of unique SMILES present in the GNPS library summary =",
          library_summary_enriched.Smiles.unique().size - 1)  # remove one NA
    
    # store the library summary enriched and validated
    library_summary_enriched.to_csv(output_filename, sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser(description='Merge GNPS library summary with enriched summary data provided by NP3 adm.')
    parser.add_argument("librarysummary", type=str, help="Library Summary table in tsv format")
    parser.add_argument("librarysummary_annotations", type=str, help="Library Summary table with annotations in tsv format")
    parser.add_argument("output_filename")
    args = parser.parse_args()

    join_library_summary_annotations(args.librarysummary, args.librarysummary_annotations, args.output_filename)

if __name__ == "__main__":
    main()