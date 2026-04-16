#!/usr/bin/python


import sys
import getopt
import os
import pandas as pd
import argparse
import glob


def main():
    # Parsing the arguments
    parser = argparse.ArgumentParser(description='Formatting data results')
    parser.add_argument('input_full_results_file', help='input_full_results_file')
    parser.add_argument('output_top1_file', help='output_top1_file')

    args = parser.parse_args()

    results_df = pd.read_csv(args.input_full_results_file, sep="\t")

    # grouping by filename and scan, sorting by MQScore
    results_df = results_df.sort_values(by=['MQScore'], ascending=[False])

    # getting the top1 hits
    top1_df = results_df.groupby(['SpectrumFile', '#Scan#']).head(1)

    # Outputting
    top1_df.to_csv(args.output_top1_file, sep="\t", index=False)

if __name__ == "__main__":
    main()
