#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 09 2021

@author: np3
"""
import pandas as pd


def str2bool(v):
    return v.lower() in ('true', '1', 't')

def diagnostic_joined_isomers(count_table_path, overwrite_table=True):
    print("\n* Computing the diagnostic of joined isomers *")
    count_table = pd.read_csv(count_table_path)
    count_table['numSamples'] = 0
    count_table['numSamplesJoinedIsomers'] = 0
    count_table['numRealPeaks'] = 0
    count_table['numFakePeaks'] = 0
    # compute the number of samples where each mz appears, the number of samples with joined isomers (more than one peakIds),
    # the number of real peaks and the number of fake peaks
    for i in range(count_table.shape[0]):
        # split the peaksIds and also split the numeric ID from the sample code
        peakIds_sample = [x.split('_',1) for x in count_table.peakIds[i].split(';')]
        n_peaks = len(peakIds_sample)
        n_samples = len(pd.unique([peakId[1] for peakId in peakIds_sample]))
        # get the samples where there is a real peak
        samples_realPeaks = pd.Series([peakId[1] for peakId in peakIds_sample if peakId[0] != 'fake'], dtype=str)
        n_samples_realPeaks = samples_realPeaks.unique().size
        n_realPeaks = len(samples_realPeaks)
        n_fakePeaks = n_peaks - n_realPeaks
        # get the number of samples which there is only one peak (no isomer was joined),
        # and then compute the number of samples where there was at least one joined isomer
        n_samples_noJoinedIsomer = samples_realPeaks[~samples_realPeaks.duplicated(False)].size
        n_samples_joinedIsomers = n_samples_realPeaks-n_samples_noJoinedIsomer
        count_table.loc[i,
                        ['numSamples', 'numSamplesJoinedIsomers',
                         'numRealPeaks', 'numFakePeaks']] = [n_samples, n_samples_joinedIsomers,
                                                             n_realPeaks, n_fakePeaks]

    if 'BLANKS_TOTAL' in count_table.columns:
        not_blank_mz = (count_table.BLANKS_TOTAL == 0)
        print("\n- Removed msclusterIDs that appears in blank samples")
    else:
        not_blank_mz = (count_table.msclusterID >= 0)

    n_mzs_joinedIsomers = (count_table[not_blank_mz].numSamplesJoinedIsomers > 0).sum()
    n_mzs_isolatedPeaks = ((count_table[not_blank_mz].numSamples == 1) &
                           (count_table[not_blank_mz].numSamplesJoinedIsomers == 0)).sum()
    n_mzs_noJoinedIsomer = count_table[not_blank_mz].shape[0]-n_mzs_joinedIsomers
    n_mzs_fakePeaks = (count_table[not_blank_mz].numFakePeaks > 0).sum()
    n = count_table.shape[0] - (~not_blank_mz).sum()
    print("\n- Percentage of msclusterIDs whit joined isomers: %0.1f%%" %
          (n_mzs_joinedIsomers / n * 100))
    print("  - Joined isomers are counted when there is more than one real peakId from the same sample code,",
          "fake peaks are ignored")
    print("- Percentage of msclusterIDs with NO joined isomer: %0.1f%%" % (
            n_mzs_noJoinedIsomer / n * 100))
    print("  - Percentage of msclusterIDs with NO joined isomer and and",
          "that appears in only one peakId from a single sample code (isolated MS1 peaks):",
          "%0.1f%%" % (n_mzs_isolatedPeaks / n * 100))
    print("  - Percentage of msclusterIDs with NO joined isomer and NO isolated MS1 peak: %0.1f%%"
          % ((n_mzs_noJoinedIsomer - n_mzs_isolatedPeaks) / n * 100))

    print("- Percentage of msclusterIDs with fake peaks: %0.1f%%" % ((n_mzs_fakePeaks) / n * 100))
    print("")
    if overwrite_table:
        count_table.to_csv(count_table_path, index=False)

    return count_table


# parse args
if __name__ == "__main__":
    import argparse
    # Parser for arguments
    parser = argparse.ArgumentParser()

    # Ok
    parser.add_argument('-c', '--count_file_path',
                        help='Path to the .csv clean count table without correlation (header in the first row)',
                        type=str, required=True)
    parser.add_argument('-o', '--overwrite_file',
                        help="Boolean, if True overwrites the count_file_path table adding new columns",
                        type=str2bool, default=True)
    args = parser.parse_args()

    diagnostic_joined_isomers(args.count_file_path, args.overwrite_file)

