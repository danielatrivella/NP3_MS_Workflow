# # @Crisfbazz
# adapted from https://github.com/Wang-Bioinformatics-Lab/NextflowModules/blob/a31f16297ea5d7f996edb70897179b2e190b9778/bin/library_search/gnps_index.py
# not support mzml spectrum files, remove dependency to pymzml

import numba
import numpy as np
from numba.typed import List as ListNumba
import os
import concurrent.futures
#import pymzml
import math
import time
import argparse

SHIFTED_OFFSET = 6000
MINMATCHES = 6
TOPPRODUCTS = 30000
ADJACENT_BINS = np.array([-1, 0, 1], dtype=np.int64)

SpectrumTuple = numba.types.Tuple([
    numba.float64[:],  # mz
    numba.float64[:],  # intensity
    numba.float32,  # precursor_mz
    numba.int32  # precursor_charge
])


class FenwickTree:
    def __init__(self, size):
        self.size = size
        self.tree = np.zeros(size + 1, dtype=np.float64)

    def add(self, position):
        # position = self.size - position - 1
        while position < self.size:
            self.tree[position] += 1
            position = position | (position + 1)

    def query(self, position):
        sum = 0
        # position = self.size - position - 1
        while (position >= 0):
            sum += self.tree[position]
            position = (position & (position + 1)) - 1
        
        return sum

def filter_window_peaks_optimized(
    mz_array: np.ndarray,
    intensity_array: np.ndarray,
    window_size: float = 50,
    peaks_to_keep_in_window: int = 6
) -> (np.ndarray, np.ndarray):
    """
    For each peak, calculate the number of peaks within its m/z window 
    (current_mz - window_size, current_mz + window_size) 
    that have larger intensity. Keep only the peaks that have less 
    than `peaks_to_keep_in_window` such peaks in their window.
    Assumes the input mz_array is sorted.
    """
    
    if len(mz_array) == 0:
        return mz_array, intensity_array

    n = len(mz_array)
    
    rev_sorted_intensities = np.argsort(-intensity_array)
    fenwick_tree = FenwickTree(n)
    res = []
    
    for i in range(n):
        current_mz = mz_array[rev_sorted_intensities[i]]
        current_intensity = intensity_array[rev_sorted_intensities[i]]
        lower_mz_bound = current_mz - window_size
        start_idx = np.searchsorted(mz_array, lower_mz_bound, side='left')
        left = fenwick_tree.query(rev_sorted_intensities[i]) - fenwick_tree.query(start_idx)
        end_idx = np.searchsorted(mz_array, current_mz + window_size, side='right')
        right = fenwick_tree.query(end_idx) - fenwick_tree.query(rev_sorted_intensities[i])
        
        if left + right <= peaks_to_keep_in_window:
            res.append((current_mz, current_intensity))
        
        # Update Fenwick tree with the current peak
        fenwick_tree.add(rev_sorted_intensities[i])


    mz_array, intensity_array = zip(*res)
    # sort the results by m/z
    index = np.argsort(mz_array)
    mz_array = np.array(mz_array, dtype=np.float64)[index]
    intensity_array = np.array(intensity_array, dtype=np.float64)[index]
    return np.array(mz_array, dtype=np.float64), np.array(intensity_array, dtype=np.float64)


def smart_filter_window_peaks_optimized(
    mz_array: np.ndarray,
    intensity_array: np.ndarray,
    window_size: float = 50,
    peaks_to_keep_in_window: int = 6
) -> (np.ndarray, np.ndarray):
    """
    For each peak, calculate the number of peaks within its m/z window 
    (current_mz - window_size, current_mz + window_size) 
    that have larger intensity. Keep only the peaks that have less 
    than `peaks_to_keep_in_window` such peaks in their window.
    Assumes the input mz_array is sorted.
    """
    
    if len(mz_array) == 0:
        return mz_array, intensity_array
    
    #---------------- this is to improve speed, it will cause change in result compared to legacy GNPS ----------------
    mz_array, intensity_array = filter_window_peaks(
        mz_array, intensity_array, window_size=window_size//4, peaks_to_keep_in_window=peaks_to_keep_in_window
    )

    n = len(mz_array)
    # Create a boolean mask to mark peaks to keep
    keep_mask = np.zeros(n, dtype=bool)

    for i in range(n):
        current_mz = mz_array[i]
        current_intensity = intensity_array[i]

        #---------------correct way to do it:----------------
        # count_larger_in_window = np.sum(
        #     (mz_array > (current_mz - window_size)) &
        #     (mz_array < (current_mz + window_size)) &
        #     (intensity_array > current_intensity)
        # )
        
        #----------------gnps legacy way:----------------
        mask = ((mz_array > (current_mz - window_size)) &
            (mz_array < (current_mz + window_size)) &
            (intensity_array > current_intensity))
        #------------------------------------------------
        
        # count unique values in the intensity array:
        count_larger_in_window = len(np.unique(intensity_array[mask]))
        

        if count_larger_in_window < peaks_to_keep_in_window:
            keep_mask[i] = True # Mark this peak to be kept

    return mz_array[keep_mask], intensity_array[keep_mask]

def filter_window_peaks(
    mz_array: np.ndarray,
    intensity_array: np.ndarray,
    window_size: float = 50,
    peaks_to_keep_in_window: int = 6
) -> (np.ndarray, np.ndarray):
    """
    Sliding‐window peak filter.  For each peak i, consider all peaks j
    with |mz[j] - mz[i]| <= window_size.  Keep peak i if its intensity
    is among the top `peaks_to_keep_in_window` in that local window.

    Returns filtered (mz, intensity), sorted by mz ascending.
    """
    if len(mz_array) == 0:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)

    # 1) find the maximum mass to size our buckets
    max_mass = np.max(mz_array)

    # 2) create buckets
    num_buckets = int(int(max_mass) // window_size) + 1
    buckets = [list() for _ in range(num_buckets)]

    # 3) assign each peak to its bucket
    for mass, intensity in zip(mz_array, intensity_array):
        idx = int(int(mass) // window_size)
        buckets[idx].append((mass, intensity))

    # 4) sort each bucket by intensity descending
    for bucket in buckets:
        bucket.sort(key=lambda peak: peak[1], reverse=True)

    # 5) collect up to max_peaks per bucket
    filtered = list()
    for bucket in buckets:
        if len(bucket) > peaks_to_keep_in_window:
            bucket = bucket[:peaks_to_keep_in_window]
        filtered.extend(bucket)

    # 6) sort the final list by mass ascending
    filtered.sort(key=lambda peak: peak[0])
    # 7) convert to numpy arrays
    mz_filtered = np.array([peak[0] for peak in filtered], dtype=np.float64)
    int_filtered = np.array([peak[1] for peak in filtered], dtype=np.float64)
    return mz_filtered, int_filtered


def filter_around_precursor(
    mz_array: np.ndarray,
    intensity_array: np.ndarray,
    precursor_mz,
    window_size = 17.0
    ) -> (np.ndarray, np.ndarray):
    """
    Filter peaks around the precursor m/z value.
    Peaks within 50 Da of the precursor m/z are kept.
    Returns filtered (mz, intensity), sorted by mz ascending.
    """
    if len(mz_array) == 0:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)
    
    # Create a mask for peaks within the window size of the precursor m/z
    mask = np.abs(mz_array - precursor_mz) >= window_size
    
    mz_filtered = mz_array[mask]
    int_filtered = intensity_array[mask]
    if len(mz_filtered) == 0:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)
    else:
        return mz_filtered.astype(np.float64), int_filtered.astype(np.float64)
        

def filter_peaks_optimized(mz_array, intensity_array, precursor_mz, precursor_charge,
                           filter_precursor=True, filter_window=True, filter_precursor_radius=17.0,
                           filter_window_size=50.0, peaks_to_keep_in_window=6):
    """Only apply sqrt transform and L2 normalization without filtering"""
    if len(mz_array) == 0:
        return (
            np.array([], dtype=np.float64),
            np.array([], dtype=np.float64),
            precursor_mz,
            precursor_charge
        )

    # Sort by m/z first
    sorted_idx = np.argsort(mz_array)
    mz_sorted = mz_array[sorted_idx].astype(np.float64)
    int_sorted = intensity_array[sorted_idx].astype(np.float64)
    
    if filter_precursor:
        # Filter around precursor m/z
        mz_sorted, int_sorted = filter_around_precursor(
            mz_sorted, int_sorted, precursor_mz, window_size=filter_precursor_radius
        )
    
    if filter_window:
        mz_sorted, int_sorted = smart_filter_window_peaks_optimized(
            mz_sorted, int_sorted,
            window_size=filter_window_size,
            peaks_to_keep_in_window=peaks_to_keep_in_window
        )
    

    int_sorted = np.sqrt(int_sorted)
    # Apply sqrt transform and L2 normalization
    norm = np.linalg.norm(int_sorted)
    if norm > 0:
        int_sorted /= norm

    return (mz_sorted, int_sorted, np.float32(precursor_mz), np.int32(precursor_charge)) # precursor mass is converted to float 32 to match legacy GNPS result.

def parse_mgf_file(path, filter_precursor=True, filter_window=True,
                     filter_precursor_radius=17.0, filter_window_size=50.0,
                     peaks_to_keep_in_window=6):
    results=[]
    with open(path,'r') as f:
        in_ions=False
        current={"mz":[],"int":[],"pepmass":None,"charge":None}
        for line in f:
            line=line.strip()
            if not line:
                continue
            if line=="BEGIN IONS":
                in_ions=True
                current={"mz":[],"int":[],"pepmass":None,"charge":None}
            elif line=="END IONS":
                if current["pepmass"] and current["charge"]:
                    pepmass_val=float(current["pepmass"])
                    cstr = current["charge"].rstrip("+")
                    if cstr.isdigit():
                        charge_val=int(cstr)
                    else:
                        charge_val=1
                    # Always create mz and intensity arrays, even if empty
                    mz_arr = np.array(current["mz"], dtype=np.float64)
                    in_arr = np.array(current["int"], dtype=np.float64)
                    filtered = filter_peaks_optimized(mz_arr, in_arr, pepmass_val, charge_val,
                                                       filter_precursor=filter_precursor,
                                                       filter_window=filter_window,
                                                       filter_precursor_radius=filter_precursor_radius,
                                                       filter_window_size=filter_window_size,
                                                       peaks_to_keep_in_window=peaks_to_keep_in_window)
                    # Append to results regardless of mz/int being empty
                    results.append(filtered)
                in_ions=False
            else:
                if in_ions:
                    if "=" in line:
                        key,val=line.split("=",1)
                        key=key.lower()
                        if key=="pepmass":
                            val_parts=val.split()
                            current["pepmass"]=val_parts[0]
                        elif key=="charge":
                            current["charge"]=val
                    else:
                        parts=line.split()
                        if len(parts)==2:
                            try:
                                mz_v=float(parts[0])
                                in_v=float(parts[1])
                                current["mz"].append(mz_v)
                                current["int"].append(in_v)
                            except:
                                pass
    return results

'''
def parse_mzml_file(path, filter_precursor=True, filter_window=True,
                     filter_precursor_radius=17.0, filter_window_size=50.0,
                     peaks_to_keep_in_window=6):
    """Parse mzML file and filter peaks based on precursor m/z and intensity.
    Args:
        path (str): Path to the mzML file.
        filter_precursor (bool): Whether to filter based on precursor m/z.
        filter_window (bool): Whether to apply window filtering.
        filter_precursor_radius (float): Radius for filtering around precursor m/z.
        filter_window_size (float): Size of the window for filtering.
        peaks_to_keep_in_window (int): Number of peaks to keep in the window.
    """
    
    run = pymzml.run.Reader(path)
    results=[]
    for spec in run:
        if spec.ms_level==2:
            if not spec.selected_precursors:
                continue
            mz_arr=np.array(spec.mz,dtype=np.float64)
            in_arr=np.array(spec.i,dtype=np.float64)
            precursor_mz=float(spec.selected_precursors[0]['mz'])
            precursor_charge=spec.selected_precursors[0].get('charge',1)
            if not isinstance(precursor_charge,int):
                try:
                    precursor_charge=int(precursor_charge)
                except:
                    precursor_charge=1
            filtered=filter_peaks_optimized(mz_arr,in_arr,precursor_mz,precursor_charge,
                                              filter_precursor=filter_precursor,
                                                filter_window=filter_window,
                                                filter_precursor_radius=filter_precursor_radius,
                                                filter_window_size=filter_window_size,
                                                peaks_to_keep_in_window=peaks_to_keep_in_window)
            
            if len(filtered[0])>0:
                results.append(filtered)
    return results
'''

def parse_one_file(file_path, filter_precursor=True, filter_window=True,
                        filter_precursor_radius=17.0, filter_window_size=50.0,
                        peaks_to_keep_in_window=6):
    # get file extension for absolute path
    # if path.startswith("/"):
    #     ext = os.path.splitext(path)[1]
    # else:
    ext = os.path.splitext(os.path.basename(file_path))[1]
    # ext to lower
    ext = ext.lower()
    print(f"Parsing file {file_path} with extension {ext}")
    if ext==".mgf":
        return parse_mgf_file(file_path, filter_precursor=filter_precursor,
                                filter_window=filter_window,
                                filter_precursor_radius=filter_precursor_radius,
                                filter_window_size=filter_window_size,
                                peaks_to_keep_in_window=peaks_to_keep_in_window)
    #elif ext==".mzml":
    #    return parse_mzml_file(file_path, filter_precursor=filter_precursor,
    #                            filter_window=filter_window,
    #                            filter_precursor_radius=filter_precursor_radius,
    #                            filter_window_size=filter_window_size,
    #                            peaks_to_keep_in_window=peaks_to_keep_in_window)
    else:
        raise ValueError(f"Unsupported file extension: {ext}. Supported extensions are .mgf and .mzml.")

def parse_files_in_parallel(file_paths, threads=1, filter_precursor=True, filter_window=True,
                            filter_precursor_radius=17.0, filter_window_size=50.0,
                            peaks_to_keep_in_window=6):
    if threads<=1 or len(file_paths)<=1:
        all_spectra=[]
        for fp in file_paths:
            parsed=parse_one_file(fp, filter_precursor=filter_precursor,
                                  filter_window=filter_window,
                                  filter_precursor_radius=filter_precursor_radius,
                                  filter_window_size=filter_window_size,
                                  peaks_to_keep_in_window=peaks_to_keep_in_window)
            all_spectra.extend(parsed)
        return all_spectra
    else:
        all_spectra=[]
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
            futs=[ex.submit(parse_one_file,fp, filter_precursor=filter_precursor,
                             filter_window=filter_window,
                             filter_precursor_radius=filter_precursor_radius,
                             filter_window_size=filter_window_size,
                             peaks_to_keep_in_window=peaks_to_keep_in_window) for fp in file_paths]
            for fut in concurrent.futures.as_completed(futs):
                all_spectra.extend(fut.result())
        return all_spectra


@numba.njit
def create_index(spectra, is_shifted, tolerance, shifted_offset):
    """Create index structure for fast lookups"""
    entries = ListNumba()
    for spec_idx in range(len(spectra)):
        spec = spectra[spec_idx]
        mz_arr = spec[0]
        precursor_mz = spec[2]

        for peak_idx in range(len(mz_arr)):
            mz = mz_arr[peak_idx]
            intensity = spec[1][peak_idx]
            
            if is_shifted:
                bin_val = np.int64(round(
                    (precursor_mz - mz + shifted_offset) / tolerance
                ))
            else:
                bin_val = np.int64(round(mz / tolerance))

            entries.append((bin_val, spec_idx, peak_idx, mz, intensity))
            
    entries.sort()
    return entries


# @numba.njit
# def find_bin_range(entries, target_bin):
#     """Binary search for bin ranges"""
#     left = 0
#     right = len(entries)
#     while left < right - 1:
#         mid = (left + right) // 2
#         if entries[mid][0] <= target_bin:
#             left = mid
#         else:
#             right = mid
#     start = left

#     right = len(entries)
#     while left < right - 1:
#         mid = (left + right) // 2
#         if entries[mid][0] <= target_bin:
#             left = mid
#         else:
#             right = mid
#     return start, right

@numba.njit
def find_bin_range(entries, target_bin):
    """Binary search for bin ranges"""
    left = 0
    right = len(entries)
    while left < right:
        mid = (left + right) // 2
        if entries[mid][0] < target_bin:
            left = mid + 1
        else:
            right = mid
    start = left

    right = len(entries)
    while left < right:
        mid = (left + right) // 2
        if entries[mid][0] <= target_bin:
            left = mid + 1
        else:
            right = mid
    return start, left

@numba.njit
def decode_ij(k, N):
    a = 1.0
    b = -(2 * N - 1.0)
    c = 2.0 * k
    disc = b * b - 4 * a * c
    sqrt_disc = math.sqrt(disc)
    i = int(math.floor((2 * N - 1 - sqrt_disc) / 2))
    total_before_i = i * (2 * N - i - 1) // 2
    j = i + (k - total_before_i) + 1
    return i, j


@numba.njit
def compute_all_pairs(spectra, shared_entries, shifted_entries, tolerance, threshold, scoring_func, query_start, query_end):
    results = ListNumba()
    n_spectra = len(spectra)

    for query_idx in range(query_start, query_end+1):
        query_spec = spectra[query_idx]
        upper_bounds = np.zeros(n_spectra, dtype=np.float64)
        match_counts = np.zeros(n_spectra, dtype=np.int32)

        # Process both shared and shifted peaks
        for peak_idx in range(len(query_spec[0])):
            mz = query_spec[0][peak_idx]
            intensity = query_spec[1][peak_idx]
            precursor_mz = query_spec[2]

            # Shared peaks processing
            shared_bin = np.int64(round(mz / tolerance))
            shifted_bin = np.int64(round((precursor_mz - mz + SHIFTED_OFFSET) / tolerance))

            # Check both shared and shifted entries
            for entries, bin_val in [(shared_entries, shared_bin),
                                     (shifted_entries, shifted_bin)]:
                for delta in ADJACENT_BINS:
                    target_bin = bin_val + delta
                    start, end = find_bin_range(entries, target_bin)
                    pos = start
                    while pos < end and entries[pos][1] <= query_idx:
                        pos += 1
                    # Find matches in this bin
                    while pos < end and entries[pos][0] == target_bin:
                        spec_idx = entries[pos][1]
                        upper_bounds[spec_idx] += intensity * entries[pos][4]
                        match_counts[spec_idx] += 1
                        pos += 1

        # Collect candidates using threshold parameter
        candidates = ListNumba()
        for spec_idx in range(query_idx + 1, n_spectra):
            if (upper_bounds[spec_idx] >= threshold and
                    match_counts[spec_idx] >= MINMATCHES):
                candidates.append((spec_idx, upper_bounds[spec_idx]))

        # Sort candidates by upper bound score
        candidates.sort(key=lambda x: -x[1])

        # Process top candidates for exact matching
        exact_matches = ListNumba()
        for spec_idx, _ in candidates[:TOPPRODUCTS * 2]:
            target_spec = spectra[spec_idx]
            score, shared, shifted, num_matches = scoring_func(spectra[query_idx], target_spec,tolerance)
            if score >= threshold:
                exact_matches.append((spec_idx, score, shared, shifted))

        # Sort and store top results
        exact_matches.sort(key=lambda x: -x[1])
        results.append((query_idx, exact_matches[:TOPPRODUCTS]))

    return results

@numba.njit(fastmath=True)
def calculate_exact_score_GNPS(query_spec, target_spec, TOLERANCE):
    """Numba-optimized cosine scoring with shift handling"""
    q_mz = query_spec[0]
    q_int = query_spec[1]
    q_pre = query_spec[2]
    q_charge = query_spec[3]

    t_mz = target_spec[0]
    t_int = target_spec[1]
    t_pre = target_spec[2]

    # Calculate precursor mass difference (assuming charge=1)
    precursor_mass_diff = (q_pre - t_pre)*q_charge
    allow_shift = True
    fragment_tol = TOLERANCE

    # Pre-allocate arrays for matches (adjust size as needed)
    max_matches = len(q_mz) * 2  # Estimate maximum possible matches
    scores_arr = np.zeros(max_matches, dtype=np.float64)
    idx_q = np.zeros(max_matches, dtype=np.int32)
    idx_t = np.zeros(max_matches, dtype=np.int32)
    match_count = 0

    # print("log", f"Processing query spectrum with {len(q_mz)} peaks against target with {len(t_mz)} peaks")
    # For each peak in query spectrum
    for q_idx in range(len(q_mz)):
        q_mz_val = q_mz[q_idx]
        q_int_val = q_int[q_idx]

        # For each possible shift (charge=1)
        num_shifts = 1
        if allow_shift and abs(precursor_mass_diff) >= fragment_tol:
            num_shifts += 1

        for shift_idx in range(num_shifts):
            if shift_idx == 0:
                # No shift
                adjusted_mz = q_mz_val
            else:
                # Apply shift
                adjusted_mz = q_mz_val - precursor_mass_diff

            # Find matching peaks in target using binary search
            start = np.searchsorted(t_mz, adjusted_mz - fragment_tol)
            end = np.searchsorted(t_mz, adjusted_mz + fragment_tol)

            for t_idx in range(start, end):
                if match_count >= max_matches:
                    break  # Prevent overflow
                if abs(t_mz[t_idx] - adjusted_mz)<fragment_tol:
                    # Store match information
                    scores_arr[match_count] = q_int_val * t_int[t_idx]
                    idx_q[match_count] = q_idx
                    idx_t[match_count] = t_idx
                    match_count += 1

    # Sort matches by score descending using argsort
    if match_count == 0:
        return 0.0, 0.0, 0.0, 0.0

    # Get valid matches
    valid_scores = scores_arr[:match_count]
    valid_q_idx = idx_q[:match_count]
    valid_t_idx = idx_t[:match_count]

    # Argsort descending
    sort_order = np.argsort(-valid_scores)

    # Track used peaks
    q_used = np.zeros(len(q_mz), dtype=np.bool_)
    t_used = np.zeros(len(t_mz), dtype=np.bool_)
    total = 0.0
    shared = 0.0
    shifted = 0.0
    num_matches = 0
    # Accumulate top matches
    for i in sort_order:
        q_idx = valid_q_idx[i]
        t_idx = valid_t_idx[i]

        if not q_used[q_idx] and not t_used[t_idx]:
            score = valid_scores[i]
            total += score

            # Determine match type
            mz_diff = abs(q_mz[q_idx] - t_mz[t_idx])
            if mz_diff <= fragment_tol:
                shared += score
            else:
                shifted += score
            num_matches += 1
            q_used[q_idx] = True
            t_used[t_idx] = True

    return total, shared, shifted, num_matches


@numba.njit
def calculate_exact_score_GNPS_multi_charge(query_spec, target_spec, TOLERANCE):
    """Numba-optimized cosine scoring with shift handling based on original GNPS code."""
    q_mz = query_spec[0]
    q_int = query_spec[1]
    q_pre = query_spec[2]
    q_charge = max(query_spec[3], 1)  # Ensure charge is at least 1

    t_mz = target_spec[0]
    t_int = target_spec[1]
    t_pre = target_spec[2]
    t_charge = max(target_spec[3], 1)  # Not used in calculation

    allow_shift = True
    fragment_tol = TOLERANCE

    # Compute precursor mass difference based on query's charge
    precursor_mass_diff = (q_pre - t_pre) * q_charge
    num_shifts = 1
    if allow_shift and abs(precursor_mass_diff) >= fragment_tol:
        num_shifts += q_charge

    # Pre-allocate arrays for matches
    max_matches = len(q_mz) * num_shifts * 10  # Estimate, adjust as needed
    scores_arr = np.zeros(max_matches, dtype=np.float64)
    idx_q = np.zeros(max_matches, dtype=np.int32)
    idx_t = np.zeros(max_matches, dtype=np.int32)
    is_shifted = np.zeros(max_matches, dtype=np.bool_)
    match_count = 0

    for q_idx in range(len(q_mz)):
        q_mz_val = q_mz[q_idx]
        q_int_val = q_int[q_idx]

        for cpi in range(num_shifts):
            # Calculate mass difference for this shift
            if cpi == 0:
                mass_diff_val = 0.0
            else:
                mass_diff_val = precursor_mass_diff / cpi

            adjusted_query_mz = q_mz_val - mass_diff_val
            start = np.searchsorted(t_mz, adjusted_query_mz - fragment_tol)
            end = np.searchsorted(t_mz, adjusted_query_mz + fragment_tol)

            # Collect all target peaks in this range
            for t_idx in range(start, end):
                if match_count >= max_matches:
                    break  # Prevent overflow

                scores_arr[match_count] = q_int_val * t_int[t_idx]
                idx_q[match_count] = q_idx
                idx_t[match_count] = t_idx
                is_shifted[match_count] = cpi != 0
                match_count += 1

    # Handle cases with no matches
    if match_count == 0:
        return 0.0, 0.0, 0.0, 0.0

    # Extract valid entries
    valid_scores = scores_arr[:match_count]
    valid_q_idx = idx_q[:match_count]
    valid_t_idx = idx_t[:match_count]
    valid_shifted = is_shifted[:match_count]

    # Sort matches by descending score
    sort_order = np.argsort(-valid_scores)

    # Track used peaks
    q_used = np.zeros(len(q_mz), dtype=np.bool_)
    t_used = np.zeros(len(t_mz), dtype=np.bool_)
    total = 0.0
    shared = 0.0
    shifted = 0.0
    num_matches = 0
    # Accumulate scores greedily
    for i in sort_order:
        q_i = valid_q_idx[i]
        t_i = valid_t_idx[i]
        if not q_used[q_i] and not t_used[t_i]:
            score = valid_scores[i]
            total += score
            if valid_shifted[i]:
                shifted += score
            else:
                shared += score
            q_used[q_i] = True
            t_used[t_i] = True
            num_matches += 1

    return total, shared, shifted, num_matches

@numba.njit
def calculate_exact_score(query_spec, target_spec, TOLERANCE):
    """Calculate exact cosine similarity with tolerance checks"""
    # Access tuple elements by index
    q_mz = query_spec[0]  # Already sorted by mz
    q_int = query_spec[1]
    q_pre = query_spec[2]

    t_mz = target_spec[0]  # Already sorted by mz
    t_int = target_spec[1]
    t_pre = target_spec[2]

    q_used = np.zeros(len(q_mz), dtype=np.bool_)
    t_used = np.zeros(len(t_mz), dtype=np.bool_)
    total = 0.0
    shared = 0.0
    shifted = 0.0

    # 1. Shared peak matching
    for q_idx in range(len(q_mz)):
        if q_used[q_idx]:
            continue

        mz_q = q_mz[q_idx]
        # Find matches using binary search on sorted mz
        start = np.searchsorted(t_mz, mz_q - TOLERANCE)
        end = np.searchsorted(t_mz, mz_q + TOLERANCE)

        best_score = 0.0
        best_t_idx = -1
        for t_idx in range(start, end):
            if not t_used[t_idx]:
                current_score = q_int[q_idx] * t_int[t_idx]
                if current_score > best_score:
                    best_score = current_score
                    best_t_idx = t_idx

        if best_t_idx != -1:
            total += best_score
            shared += best_score
            q_used[q_idx] = True
            t_used[best_t_idx] = True

    # 2. Shifted peak matching (critical fix)
    # Create and sort shifted mz arrays while tracking original indices
    # ---------------------------------------------------------------
    # For query
    q_shifted = q_pre - q_mz + SHIFTED_OFFSET
    q_shifted_sorted_idx = np.argsort(q_shifted)
    q_shifted_sorted = q_shifted[q_shifted_sorted_idx]

    # For target
    t_shifted = t_pre - t_mz + SHIFTED_OFFSET
    t_shifted_sorted_idx = np.argsort(t_shifted)
    t_shifted_sorted = t_shifted[t_shifted_sorted_idx]
    # ---------------------------------------------------------------

    # Match using sorted shifted mz arrays
    for q_pos in range(len(q_shifted_sorted)):
        q_orig_idx = q_shifted_sorted_idx[q_pos]
        if q_used[q_orig_idx]:
            continue

        mz_q = q_shifted_sorted[q_pos]
        start = np.searchsorted(t_shifted_sorted, mz_q - TOLERANCE)
        end = np.searchsorted(t_shifted_sorted, mz_q + TOLERANCE)

        best_score = 0.0
        best_t_pos = -1
        for t_pos in range(start, end):
            t_orig_idx = t_shifted_sorted_idx[t_pos]
            if not t_used[t_orig_idx]:
                current_score = q_int[q_orig_idx] * t_int[t_orig_idx]
                if current_score > best_score:
                    best_score = current_score
                    best_t_pos = t_pos

        if best_t_pos != -1:
            t_orig_idx = t_shifted_sorted_idx[best_t_pos]
            total += best_score
            shifted += best_score
            q_used[q_orig_idx] = True
            t_used[t_orig_idx] = True

    return total, shared, shifted


def get_query_range(n_spectra, chunk_id, total_chunks):
    if n_spectra == 0 or total_chunks == 0:
        return (0, 0)

    n = n_spectra
    k = total_chunks
    i = chunk_id

    def get_split_point(ratio):
        if ratio >= 1.0:
            return n - 1
        return (n - 1) * (1 - math.sqrt(1 - ratio))

    start_ratio = i / k
    end_ratio = (i + 1) / k

    start_cont = get_split_point(start_ratio)
    end_cont = get_split_point(end_ratio)

    # Calculate start and end indices, ensuring they are within bounds and contiguous
    start = math.ceil(start_cont)
    start = max(0, start)

    if chunk_id == k - 1:
        end = n - 1
    else:
        end = math.ceil(end_cont) - 1  # end is inclusive
        end = min(end, n - 1)

    # Ensure start does not exceed end, adjust if necessary
    if start > end:
        return (0, 0)

    return (start, end)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Nextflow-compatible spectral networking"
    )
    parser.add_argument("--chunk_id", type=int, default=0,
                        help="ID of this chunk (0-based index)")
    parser.add_argument("--total_chunks", type=int, default=1,
                        help="Total number of chunks to divide pairwise space into")
    parser.add_argument("-t", "--input_files", nargs="+", required=True,
                        help="Input MGF/mzML files")
    parser.add_argument("--tolerance", type=float, default=0.01,
                        help="Fragment mass tolerance")
    parser.add_argument("--threshold", type=float, default=0.7,
                        help="Similarity score threshold")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of processing threads")
    parser.add_argument("--alignment_strategy", type=str, default="index_single_charge",
                        help="Choosing which alignment strategy to use, options are: index_single_charge, index_multi_charge")
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Parse spectra
    spectra_list = parse_files_in_parallel(args.input_files, args.threads)
    print(f"Parsed {len(spectra_list)} spectra")

    # Convert to numba-compatible format
    numba_spectra = ListNumba()
    for spec in spectra_list:
        numba_spectra.append((
            spec[0].astype(np.float64),
            spec[1].astype(np.float64),
            np.float64(spec[2]),
            np.int32(spec[3])
        ))
    n_spectra = len(numba_spectra)
    query_start, query_end = get_query_range(n_spectra, args.chunk_id, args.total_chunks)
    print(f"Chunk {args.chunk_id}/{args.total_chunks} processing pairs {query_start} to {query_end}")
    # Build indexes
    
    if args.alignment_strategy == "index_multi_charge":
        print("Building indexes...")
        shared_idx = create_index(numba_spectra, False, args.tolerance, SHIFTED_OFFSET)
        shifted_idx = create_index(numba_spectra, True, args.tolerance, SHIFTED_OFFSET)
        
        scoring_func = calculate_exact_score_GNPS_multi_charge
    elif args.alignment_strategy == "index_single_charge":
        shared_idx = create_index(numba_spectra, False, args.tolerance, SHIFTED_OFFSET)
        shifted_idx = create_index(numba_spectra, True, args.tolerance, SHIFTED_OFFSET)

        scoring_func = calculate_exact_score_GNPS

    # Compute matches
    print("Computing pairs...")
    matches = compute_all_pairs(numba_spectra, shared_idx, shifted_idx,
                                args.tolerance, args.threshold, scoring_func,query_start=query_start,
    query_end=query_end)

    output_filename = f"{args.chunk_id}.params_aligns.tsv"
    
    # Write output in Nextflow-compatible format
    print("Writing results...")
    with open(output_filename, 'w') as f:
        f.write("CLUSTERID1\tCLUSTERID2\tDeltaMZ\tMinRatio\tCosine\tAlignScore2\tAlignScore3\n")
        for query_idx, candidates in matches:
            q_spec = numba_spectra[query_idx]
            for cand in candidates:
                t_spec = numba_spectra[cand[0]]

                cluster1 = query_idx + 1
                cluster2 = cand[0] + 1
                delta_mz = q_spec[2] - t_spec[2]
                cosine = cand[1]

                f.write(f"{cluster1}\t{cluster2}\t"
                        f"{delta_mz:.3f}\t0.000\t{cosine:.4f}\t"
                        f"0.000\t0.000\n")

    print("Done!")


if __name__ == "__main__":
    main()