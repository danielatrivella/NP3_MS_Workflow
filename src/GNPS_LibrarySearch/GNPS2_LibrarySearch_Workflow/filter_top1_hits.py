#!/usr/bin/python

# @ Cris
# code adapted from https://github.com/Wang-Bioinformatics-Lab/NextflowModules/blob/a31f16297ea5d7f996edb70897179b2e190b9778/bin/library_search/filter_top1_hits.py
# code adapted to select the top1 hits using not only the MQScore,
# # but also the LibraryQualityString and the number of shared peaks to resolve ties

import pandas as pd
import argparse
import sys


def main():
    # Parsing the arguments
    parser = argparse.ArgumentParser(description='Formatting data results')
    parser.add_argument('input_full_results_file', help='input_full_results_file')
    parser.add_argument('output_top1_file', help='output_top1_file')

    args = parser.parse_args()

    try:
        results_df = pd.read_csv(args.input_full_results_file, sep="\t", low_memory=False)
    except pd.errors.EmptyDataError:
        print("WARNING: Input result file is empty, there was no result in the library search. Skipping the top1 filter.\n")
        sys.exit()
    except:
        #open(output_filename, "w").close()
        sys.exit("ERROR: Input result file is not a valid tsv file")

    # change the LibraryQualityString column to a categorical column to allow correct sorting
    # where Challenge is Unknown Identity, open to community to help annotate - it is still an experimental data, which is better than insilico
    libraryQualityString_order = ['Insilico', 'Challenge', 'Bronze', 'Silver', 'Gold']
    results_df['LibraryQualityString'] = pd.Categorical(results_df['LibraryQualityString'],
                                                        categories=libraryQualityString_order, ordered=True)
    
    # grouping by filename and scan, sorting by MQScore, SharedPeaks and LibraryQualityString
    results_df = results_df.sort_values(by=['MQScore', 'SharedPeaks', 'LibraryQualityString'], ascending=[False, False, False])

    # getting the top1 hits
    top1_df = results_df.groupby(['SpectrumFile', '#Scan#']).head(1)

    # Outputting
    top1_df.to_csv(args.output_top1_file, sep="\t", index=False)

if __name__ == "__main__":
    main()
