import argparse
import pandas as pd
import os, sys
import numpy as np

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