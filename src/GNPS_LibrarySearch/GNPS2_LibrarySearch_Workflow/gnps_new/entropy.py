"""
This file contains code modified from the matchms project and MSEntropy project
(https://github.com/matchms/matchms)
Copyright matchms Team 2020

(https://github.com/YuanyueLi/MSEntropy)
Copyright Yuanyue Li 2023

Modified by Shipei Xing in 2024

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from typing import Tuple

import numba as nb
import numpy as np


@nb.njit
def find_matches(ref_spec_mz: np.ndarray, qry_spec_mz: np.ndarray,
                 tolerance: float, shift: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
    """Find matching peaks between two spectra."""
    matches_idx1 = np.empty(len(ref_spec_mz) * len(qry_spec_mz), dtype=np.int64)
    matches_idx2 = np.empty_like(matches_idx1)
    match_count = 0
    lowest_idx = 0

    for peak1_idx in range(len(ref_spec_mz)):
        mz = ref_spec_mz[peak1_idx]
        low_bound = mz - tolerance
        high_bound = mz + tolerance

        for peak2_idx in range(lowest_idx, len(qry_spec_mz)):
            mz2 = qry_spec_mz[peak2_idx] - shift
            if mz2 > high_bound:
                break
            if mz2 < low_bound:
                lowest_idx = peak2_idx
            else:
                matches_idx1[match_count] = peak1_idx
                matches_idx2[match_count] = peak2_idx
                match_count += 1

    return matches_idx1[:match_count], matches_idx2[:match_count]


@nb.njit
def collect_peak_pairs(ref_spec: np.ndarray, qry_spec: np.ndarray, min_matched_peak: int,
                       tolerance: float, shift: float = 0.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Find and score matching peak pairs between spectra."""
    if len(ref_spec) == 0 or len(qry_spec) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.float32)

    #     # Exact matching
    matches_idx1, matches_idx2 = find_matches(ref_spec[:, 0], qry_spec[:, 0], tolerance, 0.0)

    # If shift is not 0, perform hybrid search
    if abs(shift) > 1e-6:
        matches_idx1_shift, matches_idx2_shift = find_matches(ref_spec[:, 0], qry_spec[:, 0], tolerance, shift)
        matches_idx1 = np.concatenate((matches_idx1, matches_idx1_shift))
        matches_idx2 = np.concatenate((matches_idx2, matches_idx2_shift))

    if len(matches_idx1) < min_matched_peak:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.float32)

    # Calculate scores for matches
    scores = (ref_spec[matches_idx1, 1] + qry_spec[matches_idx2, 1]).astype(np.float32)

    # Sort by score descending
    sort_idx = np.argsort(-scores)
    return matches_idx1[sort_idx], matches_idx2[sort_idx], scores[sort_idx]


@nb.njit
def score_matches(matches_idx1: np.ndarray, matches_idx2: np.ndarray,
                  ref_spec: np.ndarray, qry_spec: np.ndarray, penalty: float):
    """Calculate entropy similarity score from matching peaks."""

    # Apply weight to intensity
    ref_spec = apply_weight_to_intensity(ref_spec)
    qry_spec = apply_weight_to_intensity(qry_spec)

    # Use boolean arrays for tracking used peaks
    used1 = np.zeros(len(ref_spec), dtype=nb.boolean)
    used2 = np.zeros(len(qry_spec), dtype=nb.boolean)

    # Initialize arrays for matched intensities
    matched_ref_intensities = np.empty(len(matches_idx1), dtype=np.float32)
    matched_qry_intensities = np.empty(len(matches_idx1), dtype=np.float32)

    used_matches = 0

    # Find best non-overlapping matches
    for i in range(len(matches_idx1)):
        idx1 = matches_idx1[i]
        idx2 = matches_idx2[i]
        if not used1[idx1] and not used2[idx2]:
            matched_ref_intensities[used_matches] = ref_spec[idx1, 1]
            matched_qry_intensities[used_matches] = qry_spec[idx2, 1]
            used1[idx1] = True
            used2[idx2] = True
            used_matches += 1

    if used_matches == 0:
        return 0.0, 0

    # Trim arrays to used matches only
    matched_ref_intensities = matched_ref_intensities[:used_matches]
    matched_qry_intensities = matched_qry_intensities[:used_matches]

    # Normalize intensities
    sum_ref = np.sum(ref_spec[:, 1])

    # For qry, add penalty to unmatched peaks
    sum_matched_qry = np.sum(matched_qry_intensities)
    sum_qry = sum_matched_qry + (np.sum(qry_spec[:, 1]) - sum_matched_qry) * (1 - penalty)

    if sum_ref <= 0 or sum_qry <= 0:
        return 0.0, used_matches

    matched_ref_intensities /= sum_ref
    matched_qry_intensities /= sum_qry

    # Calculate entropy similarity
    peak_ab_intensities = matched_ref_intensities + matched_qry_intensities
    entropy_sim = 0.0

    for i in range(used_matches):
        # Calculate joint entropy term
        if peak_ab_intensities[i] > 0:
            entropy_sim += peak_ab_intensities[i] * np.log2(peak_ab_intensities[i])

        # Subtract individual entropy terms
        if matched_ref_intensities[i] > 0:
            entropy_sim -= matched_ref_intensities[i] * np.log2(matched_ref_intensities[i])

        if matched_qry_intensities[i] > 0:
            entropy_sim -= matched_qry_intensities[i] * np.log2(matched_qry_intensities[i])

    entropy_sim /= 2
    return min(float(entropy_sim), 1.0), used_matches


@nb.njit
def apply_weight_to_intensity(peaks: np.ndarray) -> np.ndarray:
    """Apply a weight to the intensity of a spectrum based on spectral entropy."""
    if peaks.size == 0:
        return np.zeros((0, 2), dtype=np.float32)

    # normalize intensity
    peak_sum = np.sum(peaks[:, 1])
    if peak_sum <= 0:
        return peaks

    # Make a copy of peaks to avoid modifying the input
    weighted_peaks = peaks.copy()
    weighted_peaks[:, 1] /= peak_sum

    # Calculate the spectral entropy
    valid_intensities = weighted_peaks[:, 1] > 0
    entropy = -np.sum(weighted_peaks[valid_intensities, 1] *
                      np.log(weighted_peaks[valid_intensities, 1]))

    # Apply the weight if entropy is below threshold
    if entropy < 3:
        weight = 0.25 + 0.25 * entropy
        weighted_peaks[:, 1] = np.power(weighted_peaks[:, 1], weight)
        intensity_sum = np.sum(weighted_peaks[:, 1])
        if intensity_sum > 0:
            weighted_peaks[:, 1] /= intensity_sum

    return weighted_peaks


def entropy_similarity(qry_spec: np.ndarray, ref_spec: np.ndarray,
                       tolerance: float = 0.1,
                       min_matched_peak: int = 1,
                       sqrt_transform: bool = False,  # Unused
                       penalty: float = 0.,
                       shift: float = 0.0,
                       require_sorted:bool = False) -> Tuple[float, int]:
    """
    Calculate similarity between two spectra.

    Parameters
    ----------
    qry_spec: np.ndarray
        Query spectrum.
    ref_spec: np.ndarray
        Reference spectrum.
    tolerance: float
        Tolerance for m/z matching.
    min_matched_peak: int
        Minimum number of matched peaks.
    penalty: float
        Penalty for unmatched peaks. If set to 0, traditional cosine score; if set to 1, traditional reverse cosine score.
    shift: float
        Shift for m/z values. If not 0, hybrid search is performed. shift = prec_mz(qry) - prec_mz(ref)
    require_sorted: bool
        If True, requires spectra sorted by m/z. If False, sorts the spectra by m/z if needed.
    
    Returns
    -------
    score: float
        Similarity score between 0 and 1.
    n_matches: int
        Number of matched peaks.
    """
    tolerance = np.float32(tolerance)
    penalty = np.float32(penalty)
    shift = np.float32(shift)

    if qry_spec.size == 0 or ref_spec.size == 0:
        return 0.0, 0
    
    if require_sorted:
        if not np.all(qry_spec[:, 0][:-1] <= qry_spec[:, 0][1:]):   # True for array of shape [1,2]
            raise ValueError("Query spectrum is not sorted by m/z.")
        if not np.all(ref_spec[:, 0][:-1] <= ref_spec[:, 0][1:]):
            raise ValueError("Reference spectrum is not sorted by m/z.")
    else:
        if not np.all(qry_spec[:, 0][:-1] <= qry_spec[:, 0][1:]):
            qry_spec = qry_spec[np.argsort(qry_spec[:, 0])]
        if not np.all(ref_spec[:, 0][:-1] <= ref_spec[:, 0][1:]):
            ref_spec = ref_spec[np.argsort(ref_spec[:, 0])]

    # normalize the intensity
    ref_spec[:, 1] /= np.sum(ref_spec[:, 1])
    qry_spec[:, 1] /= np.sum(qry_spec[:, 1])

    matches_idx1, matches_idx2, scores = collect_peak_pairs(
        ref_spec, qry_spec, min_matched_peak, tolerance, shift
    )

    if len(matches_idx1) == 0:
        return 0.0, 0

    return score_matches(
        matches_idx1, matches_idx2, ref_spec, qry_spec, penalty
    )


if __name__ == "__main__":

    # Example usage
    peaks1 = np.array([[50, 8.0], [70, 100.0], [80, 50.0], [100, 50.0]], dtype=np.float32)

    peaks2 = np.array([[55, 38.0], [80, 66.0], [90, 999.0]], dtype=np.float32)

    # Example with standard entropy
    score, n_matches = entropy_similarity(peaks1, peaks2, tolerance=0.05, penalty=0)
    print(f"Standard Score: {score:.3f}, Matches: {n_matches}")

    # Example with enhanced reverse entropy
    score, n_matches = entropy_similarity(peaks1, peaks2, tolerance=0.05, penalty=0.6)
    print(f"Reverse Score: {score:.3f}, Matches: {n_matches}")

    # Example with traditional reverse entropy
    score, n_matches = entropy_similarity(peaks1, peaks2, tolerance=0.05, penalty=1)
    print(f"Reverse Score: {score:.3f}, Matches: {n_matches}")

    # Example with enhanced reverse entropy, analog search
    score, n_matches = entropy_similarity(peaks1, peaks2, tolerance=0.05, penalty=0.6, shift=10)
    print(f"Reverse Score: {score:.3f}, Matches: {n_matches}")
