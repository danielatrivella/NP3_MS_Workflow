import os
import gensim
import time
import numpy as np
import pandas as pd
from matchms.importing import load_from_mgf
from matchms.filtering import normalize_intensities
from matchms.filtering import add_precursor_mz
from matchms.filtering import select_by_relative_intensity
from matchms.similarity.spectrum_similarity_functions import find_matches
from spec2vec import Spec2Vec
from spec2vec import SpectrumDocument

from matchms.Fragments import Fragments
from matchms.typing import SpectrumType


# some code from https://github.com/iomega/spec2vec_gnps_data_analysis/blob/master/notebooks/
# Adapted notebooks 3 and 6
# Compute spec2vec similarities on mass spectra dataset

# code from matchms, added mz_tolerance to keep prec mz within tolerance - prevent empty spectrum
def remove_peaks_around_precursor_mz_keep_prec(spectrum_in: SpectrumType, trim_mz: float = 17, mz_tolerance: float = 0.025) -> SpectrumType:
    """Remove peaks that are within trim_mz (in Da) of
       the precursor mz, excluding the precursor peak +- mz_tolerance.

    Parameters
    ----------
    spectrum_in:
        Input spectrum.
    mz_tolerance:
        Tolerance of mz values that are not allowed to lie
        within the precursor mz. Default is 17 Da.
    """
    if spectrum_in is None:
        return None

    spectrum = spectrum_in.clone()

    precursor_mz = spectrum.get("precursor_mz", None)
    assert precursor_mz is not None, "Precursor mz absent."
    assert isinstance(precursor_mz, (float, int)), ("Expected 'precursor_mz' to be a scalar number.",
                                                    "Consider applying 'add_precursor_mz' filter first.")
    assert trim_mz >= 0, "trim_mz must be a positive scalar."
    assert mz_tolerance >= 0, "mz_tolerance must be a positive scalar."

    mzs, intensities = spectrum.peaks.mz, spectrum.peaks.intensities
    peaks_to_remove = ((np.abs(precursor_mz-mzs) <= trim_mz) & ~(np.abs(precursor_mz-mzs) <= mz_tolerance))
    assert (~peaks_to_remove).sum() > 0, "the filtered spectrum must have at least one fragment peak."
    new_mzs, new_intensities = mzs[~peaks_to_remove], intensities[~peaks_to_remove]
    spectrum.peaks = Fragments(mz=new_mzs, intensities=new_intensities)

    return spectrum

# Pre-processing of spectrum before pairwise comparison - similar to how it is done in MSCluster and NP3
# normalize spectrum, remove peaks around precursor if trim_mz is T, remove very low intensity peaks
def pre_process_spectrum(s, trim_mz=False, mz_tolerance=0.25):
    s = normalize_intensities(s)
    s = add_precursor_mz(s)
    #
    if trim_mz is True:
        s = remove_peaks_around_precursor_mz_keep_prec(s, trim_mz=20, mz_tolerance=mz_tolerance)
    # remove very low intensity peaks
    if len(s.peaks) >= 10:
        s = select_by_relative_intensity(s, intensity_from=np.quantile(s.peaks.intensities, 0.75)*0.001)
    # return spectrum
    return s

# use matchms to compute the number o matched peaks symmetrically in a list of spectra
def compute_peak_matches_symmetric(spectra_list, bin_size):
    matches_matrix = []
    for i in range(len(spectra_list)):
        # set nan values to the lower triangular matrix
        matches_i = [np.nan] * i
        for j in range(i, len(spectra_list)):
            # in the diagonal return the number of peaks
            # in the upper triangular matrix return the number of matched peaks
            if i == j:
                matches_i.append(len(spectra_list[i].peaks))
            else:
                matches_i.append(len(find_matches(spectra_list[i].peaks.mz, spectra_list[j].peaks.mz,
                                                  tolerance=bin_size)))
        matches_matrix.append(matches_i)
    # return the matrix with the number of matched peaks as a np ndarray
    return np.array(matches_matrix)


# Compute the pairwise similarity table using spec2vec for similarity comparison and
# compute number of matched peaks with matchms
# uses model spec2vec_UniqueInchikeys_ratio05_filtered_iter_50 in the spectra comparison
# scale_factor is fixed to be power of 0.5 - the same used in the model
def compute_pairwise_similarity_spec2vec(data_name, path_mgf, output_path, bin_size, trim_mz):
    # get script directory to set spec2vec model path
    # models from https://zenodo.org/records/3978054
    src_dir = os.path.dirname(os.path.abspath(__file__))
    model_file = os.path.join(src_dir, "spec2vec_models", "spec2vec_UniqueInchikeys_ratio05_filtered_iter_50.model")

    # load mgf with consensus spectra
    spectra = list(load_from_mgf(path_mgf))

    print("number of spectra:", len(spectra))
    # pre process spectra
    # apply pre processing steps to the data
    spectra_preprocessed = [pre_process_spectrum(s, trim_mz, bin_size) for s in spectra]

    # convert spectra to documents in spec2vec format
    tstart = time.time()
    spectra_documents = [SpectrumDocument(s, n_decimals=2) for s in spectra_preprocessed]
    tend = time.time()
    print(f"Time to create {len(spectra_documents)} documents: {tend - tstart} s.")

    # load spec2vec model
    model = gensim.models.Word2Vec.load(model_file)

    # compute pairwise comparison using spec2vec
    # Actual score calculation
    # Using Spec2Vec with intensity_weighting_power=0.5.
    # Calculate matrix of all-vs-all similarity scores.
    spec2vec_similarity = Spec2Vec(model, intensity_weighting_power=0.5, allowed_missing_percentage=50)
    tstart = time.time()
    similarity_matrix = np.round(spec2vec_similarity.matrix(spectra_documents, spectra_documents, is_symmetric=True),3)
    # remove lower triangular matrix values and reset diagonal with 1.0
    similarity_matrix[np.tril_indices_from(similarity_matrix)] = np.nan
    np.fill_diagonal(similarity_matrix, 1.0)
    # set negative values to 0 - expected similarity values must be >= 0
    similarity_matrix[similarity_matrix < 0.0] = 0.0
    tend = time.time()

    if tend - tstart < 120:
        print(f"Calculated {similarity_matrix.shape[0]}x{similarity_matrix.shape[1]} scores in {tend - tstart} s.")
    else:
        print(f"Calculated {similarity_matrix.shape[0]}x{similarity_matrix.shape[1]} scores in {(tend - tstart) / 60} min.")

    # also compute matched peaks
    # the number of peaks is directly influenced by the filtering method, this result will be different from the R code
    tstart = time.time()
    matched_peaks_matrix = compute_peak_matches_symmetric(spectra_preprocessed, bin_size)
    tend = time.time()
    if tend - tstart < 120:
        print(f"Calculated {matched_peaks_matrix.shape[0]}x{matched_peaks_matrix.shape[1]} matches in {tend - tstart} s.")
    else:
        print(f"Calculated {matched_peaks_matrix.shape[0]}x{matched_peaks_matrix.shape[1]} matches in {(tend - tstart) / 60} min.")

    # save similarity table in the right format
    scans_number = [s.get("scans") for s in spectra_preprocessed]
    header_index = "parameters sim pairwise - spec2vec - scale_factor:0.5;bin_size:"+str(bin_size)+";trim_mz:"+str(trim_mz)
    # Store similarity matrix
    similarity_matrix = pd.DataFrame(similarity_matrix, index=scans_number, columns=scans_number)
    similarity_matrix.to_csv(os.path.join(output_path, "similarity_table_"+data_name+".csv"),
              sep=',', na_rep="", index_label=header_index)
    # store matches matrix
    matched_peaks_matrix = pd.DataFrame(matched_peaks_matrix, index=scans_number, columns=scans_number)
    matched_peaks_matrix.to_csv(os.path.join(output_path, "similarity_table_matches_"+data_name+".csv"),
              sep=',', na_rep="", index_label=header_index)


if __name__ == "__main__":
    import sys, os
    if len(sys.argv) >= 6:
        # print(sys.argv)
        job_name = sys.argv[1]
        path_mgf = sys.argv[2]
        output_path = sys.argv[3]
        bin_size = float(sys.argv[4])
        trim_mz = bool(sys.argv[5])
    else:
        print("Error: Five arguments must be supplied to created the pairwise similarity table using spec2vec:\n",
            " 1 - job_name: The name of the job being executed;\n",
            " 2 - path_mgf: Path to the MGF file with the spectra to be compared pairwise;\n",
            " 3 - output_path: Path to the output folder to store the resulting similarity table;\n"
            " 4 - bin_size: The bin size to consider two fragmented peaks m/z's the same, used in the peaks matches procedure;\n",
            " 5 - trim_mz: A boolean indicating if the spectra should be trimmed by the precursor mass, "
            "which removes all peaks around the precursor mz +-20 Da.\n")
        sys.exit(1)
    # call compute spec2vec
    compute_pairwise_similarity_spec2vec(data_name=job_name, path_mgf=path_mgf, output_path=output_path,
                                         bin_size=bin_size, trim_mz=trim_mz)
