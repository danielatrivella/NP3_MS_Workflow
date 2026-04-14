#!/usr/bin/python


import sys
import getopt
import os
import pandas as pd
from collections import defaultdict
import argparse
import uuid
import glob

def _parse_file(input_filename):
    # checking extension
    if input_filename.endswith(".tsv"):
        df = pd.read_csv(input_filename, sep="\t")
    elif input_filename.endswith(".csv"):
        df = pd.read_csv(input_filename, sep=",")

    return df


def main():
    # Parsing the arguments
    parser = argparse.ArgumentParser(description='Merging Library Results Files')
    parser.add_argument('input_folder', help='input_folder')
    parser.add_argument('output_file', help='output_file')
    parser.add_argument('--topk', default=1, help='topk')

    args = parser.parse_args()

    all_results_files = glob.glob(os.path.join(args.input_folder, "*"))

    all_results_list = []
    for i, results_file in enumerate(all_results_files):
        temp_df = _parse_file(results_file)

        if len(temp_df) > 0:
            all_results_list.append(temp_df)
            
    if len(all_results_list) == 0:
        return
    elif len(all_results_list) == 1:
        # If there is nothing to concat, skip concatentation
        all_results_df = all_results_list[0]
    else:
        # merging results
        all_results_df = pd.concat(all_results_list, ignore_index=True)

    # Get Topk results by score per file and scan
    all_results_df["key"] = all_results_df["SpectrumFile"].astype(str) + ":" + all_results_df["#Scan#"].astype(str) 
    all_results_df = all_results_df.sort_values(by=["MQScore"], ascending=False)

    # Grouping by key and getting topk
    all_results_df = all_results_df.groupby("key").head(int(args.topk))

    # Flatten
    all_results_df = all_results_df.drop(columns=["key"])

    # writing results
    all_results_df.to_csv(args.output_file, sep="\t", index=False)

if __name__ == "__main__":
    main()
