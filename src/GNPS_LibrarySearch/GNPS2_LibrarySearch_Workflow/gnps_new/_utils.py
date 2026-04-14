from dataclasses import dataclass

import numpy as np
from numba import njit


@dataclass
class Spectrum:
    """
    A class to represent a mass spectrum with its associated metadata.
    """
    scan: object
    precursor_mz: float
    rt: float
    charge: object
    tic: float
    peaks: np.ndarray
    peaks_cleaned: bool = False

    # cleaned_peaks: np.ndarray = None

    def __post_init__(self):
        """Validate peaks data"""
        assert self.peaks.ndim == 2, "Peaks must be 2D array"
        assert self.peaks.shape[1] == 2, "Peaks must have shape (n, 2) for (mz, intensity)"

        # Ensure float32 type
        self.peaks = np.asarray(self.peaks, dtype=np.float32)


@njit
def remove_electronic_noise_peaks(peaks: np.ndarray, min_group_size: int = 5) -> np.ndarray:
    """
    Remove electronic noise peaks from MS2 spectrum.

    Parameters
    ----------
    peaks : np.ndarray
        2D array with shape (n_peaks, 2) where first column is m/z and second is intensity
    min_group_size : int, default=5
        Minimum number of peaks with similar intensities to be considered noise

    Returns
    -------
    np.ndarray
        Cleaned spectrum with noise peaks removed
    """
    if len(peaks) < min_group_size:
        return peaks

    sorted_peaks = peaks[np.argsort(peaks[:, 1])]
    max_intensity = sorted_peaks[-1, 1]

    max_start_idx = len(sorted_peaks) - min_group_size

    cutoff_intensity = 0

    for i in range(max_start_idx):
        peak_group = sorted_peaks[i:i + min_group_size, 1]
        ref_intensity = peak_group[0]

        if ref_intensity > max_intensity * 0.60:
            break

        # Check if all intensities are within 0.5% of the max intensity
        if np.all(np.abs(peak_group - ref_intensity) / max_intensity <= 0.005):
            cutoff_intensity = peak_group[-1]

    return peaks[peaks[:, 1] > cutoff_intensity]


@njit
def clean_peaks(peaks: np.ndarray,
                prec_mz: float,
                rel_int_threshold: float = 0.0,
                prec_mz_removal_da: float = 1.5,
                max_peak_num: int = 100):
    """
    Clean MS/MS peaks
    """

    # m/z and intensity must be positive
    peaks = peaks[np.bitwise_and(peaks[:, 0] > 0, peaks[:, 1] > 0)]
    if peaks.size == 0:
        return np.zeros((0, 2), dtype=np.float32)

    # Remove peaks with mz > prec_mz - prec_mz_removal_da
    max_allowed_mz = prec_mz - prec_mz_removal_da
    peaks = peaks[peaks[:, 0] < max_allowed_mz]
    if peaks.size == 0:
        return np.zeros((0, 2), dtype=np.float32)

    # Remove electronic noise peaks
    peaks = remove_electronic_noise_peaks(peaks)
    if peaks.size == 0:
        return np.zeros((0, 2), dtype=np.float32)

    # Remove low intensity peaks
    peaks = peaks[peaks[:, 1] > rel_int_threshold * np.max(peaks[:, 1])]

    if peaks.size == 0:
        return np.zeros((0, 2), dtype=np.float32)

    # Maximum number of peaks
    if max_peak_num > 0 and len(peaks) > max_peak_num:
        # Sort the spectrum by intensity.
        peaks = peaks[np.argsort(peaks[:, 1])[-max_peak_num:]]

    # Sort peaks by m/z
    peaks = peaks[np.argsort(peaks[:, 0])]

    return np.asarray(peaks, dtype=np.float32)
