import argparse
import os

import numpy as np
import pandas as pd

from _utils import clean_peaks
from cosine import cosine_similarity
from entropy import entropy_similarity
from file_io import batch_process_queries, iterate_gnps_lib_mgf

BUFFER_SIZE = 2000  # number of matched rows to write at once
QRY_BATCH_SIZE = 2500  # number of query spectra to process at once


def main_batch(gnps_lib_mgf, qry_file,
               algorithm='cos', analog_search=False, analog_max_shift=200.,
               pm_tol=0.02, frag_tol=0.05,
               min_score=0.7, min_matched_peak=3,
               rel_int_threshold=0.0, prec_mz_removal_da=1.5,
               peak_transformation='sqrt', max_peak_num=100,
               unmatched_penalty_factor=0.6,
               qry_batch_size=QRY_BATCH_SIZE):
    """
    Main function to search GNPS library

    algorithm: str. cos, entropy, rev_cos, rev_entropy
    peak_transformation: str. 'sqrt' or 'none', only applied on cosine similarity
    """

    if algorithm in ['cos', 'rev_cos']:
        search_eng = cosine_similarity
    elif algorithm in ['entropy', 'rev_entropy']:
        search_eng = entropy_similarity
    else:
        raise ValueError(f"Invalid algorithm: {algorithm}")

    search_kwargs = {
        'tolerance': frag_tol,
        'min_matched_peak': min_matched_peak,
        'sqrt_transform': True if peak_transformation == 'sqrt' else False,
        'penalty': unmatched_penalty_factor if 'rev' in algorithm else 0.0  # penalty for unmatched peaks, only applied on reverse search
    }

    # Some preprocessing
    qry_file_name = os.path.basename(qry_file)
    lib_mgf_basename = os.path.splitext(os.path.basename(gnps_lib_mgf))[0]
    qry_basename = os.path.basename(qry_file).replace("/", "_").replace(".", "_").replace(" ", "_")
    out_path = os.path.join('search_results', f"{qry_basename}_{lib_mgf_basename}_matches.tsv")

    min_matched_peak = max(min_matched_peak, 1)

    # Initialize list for batch writing
    matches_buffer = []

    bad_ref_indices = set()  # ref spectra that has peaks less than min_matched_peak
    for batch_specs, batch_prec_mzs in batch_process_queries(qry_file, min_matched_peak, qry_batch_size):

        # iterate GNPS library
        for gnps_idx, spec in enumerate(iterate_gnps_lib_mgf(gnps_lib_mgf)):

            if gnps_idx in bad_ref_indices:
                continue

            # ref peaks, number of peaks check
            if len(spec['peaks']) < min_matched_peak:
                bad_ref_indices.add(gnps_idx)
                continue

            # precursor mz check
            if analog_search:
                v = np.where(np.abs(batch_prec_mzs - spec['PEPMASS']) <= analog_max_shift)[0]
            else:
                # exact search
                v = np.where(np.abs(batch_prec_mzs - spec['PEPMASS']) <= pm_tol)[0]

            if len(v) == 0:
                continue

            # clean ref peaks
            ref_peaks = clean_peaks(spec['peaks'],
                                    spec['PEPMASS'],
                                    rel_int_threshold=rel_int_threshold,
                                    prec_mz_removal_da=prec_mz_removal_da,
                                    max_peak_num=max_peak_num)
            if len(ref_peaks) < min_matched_peak:
                bad_ref_indices.add(gnps_idx)
                continue

            for i in v:
                qry_spec = batch_specs[i]

                # clean peaks
                if not qry_spec.peaks_cleaned:
                    qry_spec.peaks = clean_peaks(qry_spec.peaks,
                                                 qry_spec.precursor_mz,
                                                 rel_int_threshold=rel_int_threshold,
                                                 prec_mz_removal_da=prec_mz_removal_da,
                                                 max_peak_num=max_peak_num)
                    qry_spec.peaks_cleaned = True
                if len(qry_spec.peaks) < min_matched_peak:
                    continue

                # calculate similarity score
                if analog_search:
                    score, n_matches = search_eng(qry_spec.peaks, ref_peaks, **search_kwargs, shift=qry_spec.precursor_mz - spec['PEPMASS'])
                else:
                    score, n_matches = search_eng(qry_spec.peaks, ref_peaks, **search_kwargs, shift=0.0)

                # filter by minimum score and minimum matched peaks
                if score < min_score or n_matches < min_matched_peak:
                    continue

                # store matched rows
                matches_buffer.append({
                    '#Scan#': qry_spec.scan,
                    'SpectrumFile': qry_file_name,
                    'Annotation': '',
                    'OrigAnnotation': '',
                    'Protein': '',
                    'dbIndex': '',
                    'numMods': '',
                    'matchOrientation': '',
                    'startMass': '',
                    'Charge': qry_spec.charge,
                    'MQScore': round(score, 4),
                    'p-value': qry_spec.rt,  # RT value
                    'isDecoy': '',
                    'StrictEnvelopeScore': '',
                    'UnstrictEvelopeScore': round(qry_spec.tic),  # TIC value
                    'CompoundName': '',
                    'Organism': '',
                    'FileScanUniqueID': f'{qry_file_name}_{qry_spec.scan}',
                    'FDR': '',
                    'LibraryName': '',
                    'mzErrorPPM': round((qry_spec.precursor_mz - spec['PEPMASS']) / spec['PEPMASS'] * 1e6, 4),
                    'LibMetaData': '',
                    'Smiles': '',
                    'Inchi': '',
                    'LibSearchSharedPeaks': n_matches,
                    'Abundance': '',
                    'ParentMassDiff': round(qry_spec.precursor_mz - spec['PEPMASS'], 4),
                    'SpecMZ': qry_spec.precursor_mz,
                    'ExactMass': '',
                    'LibrarySpectrumID': spec['SPECTRUMID'] if spec['SPECTRUMID'] != '' else f'scans_{spec["SCANS"]}',
                })

                # Write buffer to file
                if len(matches_buffer) >= BUFFER_SIZE:
                    write_batch_results(matches_buffer, out_path)
                    matches_buffer = []

    # Write remaining matches
    if matches_buffer:
        write_batch_results(matches_buffer, out_path)

    return


def write_batch_results(match_rows, out_path):
    """Write batch results to file"""

    df = pd.DataFrame(match_rows)
    # If file doesn't exist, write with header
    if not os.path.exists(out_path):
        df.to_csv(out_path, sep='\t', index=False)
    else:
        # Append without header
        df.to_csv(out_path, sep='\t', index=False, mode='a', header=False)


if __name__ == "__main__":
    argparse = argparse.ArgumentParser(description='Search GNPS library')
    argparse.add_argument('--gnps_lib_mgf', type=str, help='GNPS library MGF file path')
    argparse.add_argument('--qry_file', type=str, help='Query file path')
    argparse.add_argument('--algorithm', type=str, default='cos',
                          help='Algorithm: cos, rev_cos, entropy, rev_entropy')
    argparse.add_argument('--analog_search', type=str, default="0", help='Turn on analog search, 0 or 1')
    argparse.add_argument('--analog_max_shift', type=float, default=200, help='Analog search max shift')
    argparse.add_argument('--pm_tol', type=float, default=0.02, help='Precursor m/z tolerance')
    argparse.add_argument('--frag_tol', type=float, default=0.05, help='Fragment m/z tolerance')
    argparse.add_argument('--min_score', type=float, default=0.7, help='Minimum score')
    argparse.add_argument('--min_matched_peak', type=int, default=3, help='Minimum matched peaks')
    # argparse.add_argument('--rel_int_threshold', type=float, default=0.0, help='Relative intensity threshold')
    # argparse.add_argument('--prec_mz_removal_da', type=float, default=1.5, help='Precursor m/z removal')
    argparse.add_argument('--peak_transformation', type=str, default='sqrt',
                          help='Peak transformation, sqrt or none')
    # argparse.add_argument('--max_peak_num', type=int, default=50, help='Maximum number of peaks')
    argparse.add_argument('--unmatched_penalty_factor', type=float, default=0.6, help='Penalty factor for reverse spectral search')

    args = argparse.parse_args()

    main_batch(args.gnps_lib_mgf, args.qry_file,
               algorithm=args.algorithm, analog_search=True if args.analog_search == "1" else False,
               analog_max_shift=args.analog_max_shift,
               pm_tol=args.pm_tol, frag_tol=args.frag_tol,
               min_score=args.min_score, min_matched_peak=args.min_matched_peak,
               # rel_int_threshold=args.rel_int_threshold, prec_mz_removal_da=args.prec_mz_removal_da,
               peak_transformation=args.peak_transformation,
               # max_peak_num=args.max_peak_num,
               unmatched_penalty_factor=args.unmatched_penalty_factor)

    ############################
    # main_batch('/Users/shipei/Documents/projects/reverse_search/play/test/test_ref.mgf',
    #            '/Users/shipei/Documents/projects/reverse_search/play/test/test_qry.mgf',
    #            algorithm='cos',
    #            min_score=0.4, min_matched_peak=1)
