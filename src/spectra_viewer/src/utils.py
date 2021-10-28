# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 14:59:03 2021

@author: np3
"""

import pandas as pd
import numpy as np

from pyteomics import mgf
from matchms import Spectrum, Spikes
from matchms.filtering import default_filters

import streamlit as st
from abc import ABC, abstractmethod
from random import randint
from random import random

from src.dotprod import *


def rand_list(initial, final, len_list=100, order=True, integer=False):
    """
    Creates a list of random values, with integer or floats, ordered or messy.

    Parameters
    ----------
    initial : int
        Initial (minimum) number for the random values.
    final : int
        Final (maximum) number for the random values.
    len_list : int, optional
        Length of the list. The default is 100.
    order : bool, optional
        If True, list ir ordinated. The default is True
    integer : bool, optional
        If True, the list is of integers. The default is False.

    Returns
    -------
    The list of random values.

    """
    out = []

    # Validations
    # Correcting floats for ints
    if type(initial) == float:
        print("The please use integers next time for initial parameter,\nthe value will be truncated")
        initial = int(initial)
    if type(final) == float:
        print("The please use integers next time for final parameter,\nthe value will be truncated")
        final = int(final)
    if type(len_list) == float:
        print("The please use integers next time for len_list parameter,\nthe value will be truncated")
        len_list = int(abs(len_list))

    # Correcting order
    if final < initial:
        print("Please use the final value as the major value\nI'll correct for you this time")
        aux = initial
        initial = final
        final = aux
    elif final == initial:
        print("Please use different values for initial and final parameters,\nadding 1 at final")
        final += 1

    # Generating random list

    # Range on size of list
    for i in range(len_list):
        val = randint(initial, final)

        # Get the decimal value
        if not integer:
            increment = random()
            increment = round(increment, 4)

            # add for positive, sub for negative
            if val >= 0:
                val += increment
            else:
                val -= increment
        out.append(val)

    if order:
        out.sort()

    return out


def create_peaks(num_peak: object = 80) -> object:
    """
    Creates a fragments peak list with random values

    Parameters
    ----------
    num_peak : int, optional
        number of peaks. The default is 80.

    Returns
    -------
    A pandas data frame with the peak list.

    """
    if type(num_peak) == float:
        num_peak = int(num_peak)

    if type(num_peak) != int or num_peak < 1:
        print("Invalid input, generating standard 80 rows list")
        num_peak = 80

    mz = rand_list(initial=90,
                   final=600,
                   len_list=num_peak,
                   order=True,
                   integer=False)

    ints = rand_list(initial=15,
                     final=100,
                     len_list=num_peak,
                     order=False,
                     integer=True)

    if 100 not in ints:
        ale = randint(0, num_peak - 1)
        ints[ale] = 100

    return pd.DataFrame({"mz": mz, "int": ints})


def frags_to_text(dataframe):
    """
    From pandas dataframe fragment list,
    creates a text to input as example on forms

    Parameters
    ----------
    dataframe : DataFrame
        Fragments dataframe with mz and int.

    Returns
    -------
    the data frame in text.

    """
    frame_str = ""
    for i in dataframe.index:
        cols = list(dataframe.columns)
        for col in cols:
            frame_str += str(dataframe[col][i])
            if col != cols[-1]:
                frame_str += "\t"
            else:
                frame_str += "\n"

    return frame_str.strip()


def text_to_frags(string):
    """
    From a text fragment list,
    creates a pandas dataframe to update the plot

    Parameters
    ----------
    string : a peak list string separated by white space and \n

    Returns
    -------
    A mz, int pandas dataframe
    """

    # Clean head and tails
    string = string.strip()
    # Initialize lists
    mz = []
    intensity = []

    # run the string
    for peak in string.split("\n"):
        if len(peak) > 0:
            mass_int = peak.split()
            mz.append(float(mass_int[0]))
            intensity.append(abs(float(mass_int[1])))

    return pd.DataFrame({"mz": mz, "int": intensity})


def float_or_none(string):
    """
    Try to coerce text to float

    Parameters
    ----------
    string : str
        A float typed as string.

    Returns
    -------
    float or none.

    """
    if string == '':
        return None

    if ',' in string:
        string = string.replace(',', '.')

    try:
        out = float(string)
        return out
    except Exception:
        return None


def text_to_frags_spec(string,
                       id_spec='spectrumA',
                       title='spectrumA',
                       pepmass=None):
    """
    From a text fragment list,
    creates a matchms spectra to update the plot

    Parameters
    ----------
    string : str
        a peak list string separated by white space and \n.
    id_spec : string, optional
        DESCRIPTION. The default is 'spectrumA'.
    title : string, optional
        DESCRIPTION. Title for the spectrum
    pepmass : float, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    spectrum : Matchms spectra
        DESCRIPTION.

    """

    # Clean head and tails
    string = string.strip()
    # Initialize lists
    mz = []
    intensities = []

    # run the string
    for peak in string.split("\n"):
        if len(peak) > 0:
            mass_int = peak.split()
            mz.append(float(mass_int[0]))
            intensities.append(abs(float(mass_int[1])))

    # Sort by mz (if not sorted already)
    mz = np.array(mz)
    intensities = np.array(intensities)

    if not np.all(mz[:-1] <= mz[1:]):
        idx_sorted = np.argsort(mz)
        mz = mz[idx_sorted]
        intensities = intensities[idx_sorted]

    if (pepmass is None) or (type(pepmass) not in [float, int]):
        spectrum = Spectrum(mz=np.array(mz),
                            intensities=np.array(intensities),
                            metadata={'id': id_spec,
                                      'title': title, })
    else:
        spectrum = Spectrum(mz=np.array(mz),
                            intensities=np.array(intensities),
                            metadata={'id': id_spec,
                                      'title': title,
                                      'precursor_mz': pepmass})

    return default_filters(spectrum)


def sp_2_np(spec):
    """
    Converts spectra (matchms) to numpy 2d array
    """
    return Spikes(spec.peaks.mz, spec.peaks.intensities).to_numpy


###################################################
#######        DOT PRODUCT FUNCTIONS        #######
###################################################

def try_get_pep(spec):
    try:
        pep = spec.get('pepmass')[0]
    except Exception:
        pep = spec.get('pepmass')

    if pep is None:
        pep = spec.get('precursor_mz')

    if type(pep) == float:
        return pep
    else:
        return None


def dot_product_sp(spA, spB, bin_size=0.05):
    """
    Calc cosine from spectrumA and spectrumB

    Parameters
    ----------
    spA : matchms Spectrum
    spB : matchms Spectrum
    bin_size : float, optional
        Bin size window to merge peaks. The default is 0.025.

    Returns
    -------
    Cosine between A and B.

    """
    return normDotProduct(spA.peaks.mz,
                          spA.peaks.intensities,
                          spB.peaks.mz,
                          spB.peaks.intensities,
                          bin_size)


def dot_product_shift_sp(spA, spB, bin_size=0.05, join=True):
    """
    Calc cosine from spectrumA and spectrumB
    using Modified cosine by shift

    Parameters
    ----------
    spA  : matchms Spectrum
    spB  : matchms Spectrum
    bin_size : float, optional
        Bin size window to merge peaks. The default is 0.025.
    join : bool, optional
        If True join similar mz's at bin_size window
    Returns
    -------
    Cosine between A and B.
    """
    pepA = try_get_pep(spA)
    pepB = try_get_pep(spB)

    shift = (pepA - pepB)
    if shift is None:
        shift = 0

    # Join adjacent and/or isotopic
    pairs_a = joinAdjPeaksScale(spA.peaks.mz,  # mz list
                                spA.peaks.intensities,  # intensities list
                                bin_size,  # bin size
                                pepA,  # use pepmass to trim
                                1,  # Scale method :: if int >= 1 power of root, if 0 uses ln. So 1 does nothing
                                int(join)  # join isotopic peaks. if 1 True
                                )

    pairs_b = joinAdjPeaksScale(spB.peaks.mz,  # mz list
                                spB.peaks.intensities,  # intensities list
                                bin_size,  # bin size
                                pepB,  # use pepmass to trim
                                1,  # Scale method :: if int >= 1 power of root, if 0 uses ln. So 1 does nothing
                                int(join)  # join isotopic peaks. if 1 True
                                )

    return normDotProductShift(pairs_a[0],
                               pairs_a[1],
                               pairs_b[0],
                               pairs_b[1],
                               bin_size,
                               shift)


def dot_product(dfA, dfB, bin_size=0.05):
    """
    Calc cosine from dfA and dfB peak list

    Parameters
    ----------
    dfA : Pandas Dataframe
        Peaks dataframe with mz and int columns.
    dfB : Pandas Dataframe
        Peaks dataframe with mz and int columns.
    bin_size : float, optional
        Bin size window to merge peaks. The default is 0.02.

    Returns
    -------
    Cosine between A and B.

    """
    mz_a = dfA.mz.tolist()
    int_a = dfA.int.tolist()

    mz_b = dfB.mz.tolist()
    int_b = dfB.int.tolist()

    return normDotProduct(mz_a,
                          int_a,
                          mz_b,
                          int_b,
                          bin_size)


def dot_product_shift(dfA, dfB, pep_mA, pep_mB, bin_size=0.025):
    """
    Calc cosine from dfA and dfB peak list

    Parameters
    ----------
    dfA : Pandas Dataframe
        Peaks dataframe with mz and int columns.
    dfB : Pandas Dataframe
        Peaks dataframe with mz and int columns.
    pep_mA: float
        Parent mass of A
    pep_mB: float
        Parent mass of B
    bin_size : float, optional
        Bin size window to merge peaks. The default is 0.025.

    Returns
    -------
    Cosine between A and B with parent mass delta.

    """
    delta = abs(pep_mA - pep_mB)

    mz_a = dfA.mz.tolist()
    int_a = dfA.int.tolist()

    mz_b = dfB.mz.tolist()
    int_b = dfB.int.tolist()

    return normDotProductShift(mz_a,
                               int_a,
                               mz_b,
                               int_b,
                               bin_size,
                               delta)


###################################################
#######          SCALING FUNCTIONS          #######
###################################################
def log_sp(spectrum, dec=2):
    spec = Spectrum(mz=spectrum.peaks.mz,
                    intensities=np.round(np.log(spectrum.peaks.intensities),
                                         dec),
                    metadata=spectrum.metadata)
    return default_filters(spec)


def sqrt_sp(spectrum, dec=2):
    spec = Spectrum(mz=spectrum.peaks.mz,
                    intensities=np.round(np.sqrt(spectrum.peaks.intensities),
                                         dec),
                    metadata=spectrum.metadata)
    return default_filters(spec)


###################################################
#######        MGF RELATED FUNCTIONS        #######
###################################################

def load_from_mgf_idx(spectra_bin):
    """Load spectrum(s) from mgf file.

    Example:

    .. code-block:: python

        from matchms.importing import load_from_mgf

        file_mgf = "pesticides.mgf"
        spectrum = list(load_from_mgf(file_mgf))

    """

    for pyteomics_spectrum in mgf.IndexedMGF(spectra_bin, convert_arrays=1):

        metadata = pyteomics_spectrum.get("params", None)
        mz = pyteomics_spectrum["m/z array"]
        intensities = pyteomics_spectrum["intensity array"]

        # Sort by mz (if not sorted already)
        if not np.all(mz[:-1] <= mz[1:]):
            idx_sorted = np.argsort(mz)
            mz = mz[idx_sorted]
            intensities = intensities[idx_sorted]

        yield Spectrum(mz=mz, intensities=intensities, metadata=metadata)


def get_unique_keys_sp(spectra):
    """
    Parameters
    ----------
    spectra : pyteomics or matchms spectra
        Spectra collection converted from a mfg file

    Returns
    -------
    unique_keys : list
        list of unique keys candidates from file.
        Returns None if the list is empty
    """

    all_headers = [s.get("params") for s in spectra]
    mgf_size = len(all_headers)
    # Params for the first spectrum
    inter_keys = set(all_headers[0].keys())

    # Get parameters that belongs to all spectra
    for header in all_headers:
        inter_keys = inter_keys & set(header.keys())

    # Get values of those parameters and keep it in a dict of lists
    dict_frame = {}
    for k in inter_keys:
        dict_frame[k] = []
    for header in all_headers:
        for k in inter_keys:
            dict_frame[k].append(header.get(k))

    # Get list of keys that have a unique value
    unique_keys = []
    for k in dict_frame:
        # Try Except  if keys is charge an error is raised on set
        try:
            if len(set(dict_frame[k])) == mgf_size:
                unique_keys.append(k)
        except Exception:
            pass

    # Returns non empty list results
    if len(unique_keys) > 0:
        return unique_keys
    else:
        return None


def get_all_titles(spectra):
    """
    For a mgf read file spectra, get all titles in list

    Parameters
    ----------
    spectra: mgf file read on pyteomics

    Returns
    -------
    List of titles at the spectra
    """
    out = []
    for s in spectra:
        t = s["params"]["title"]
        if t is not None:
            out.append(t)

    out = list(set(out))
    out.sort()
    return out


def spec_to_df(spectra, title):
    """
    For a mgf read file spectra, get all titles in list

    Parameters
    ----------
    spectra: mgf file read on pyteomics
    title: The title of the peak list that you want
    Returns
    -------
    Peak list as dataframe
    """
    mz = spectra.get_spectrum(title)["m/z array"]
    intense = spectra.get_spectrum(title)["intensity array"]

    return pd.DataFrame({"mz": mz, "int": intense})


def pyteomics_2_matchms(spectra, title):
    """
    From indexed reader from pyteomics
    converts to matchms type

    Parameters
    ----------
    spectra: mgf file read on pyteomics
    title: The title of the peak list that you want
    Returns
    -------
    Matchms Spectrum
    """

    metadata = spectra.get_spectrum(title).get('params')
    mz = spectra.get_spectrum(title)['m/z array']
    intensities = spectra.get_spectrum(title)['intensity array']

    # Sort by mz (if not sorted already)
    if not np.all(mz[:-1] <= mz[1:]):
        idx_sorted = np.argsort(mz)
        mz = mz[idx_sorted]
        intensities = intensities[idx_sorted]

    return default_filters(Spectrum(mz=mz, intensities=intensities, metadata=metadata))


def get_pepmass(spectra, title):
    """
    For a mgf read file spectra it, returns the
    pepmass parameter

    Parameters
    ----------
    spectra : mgf file read on pyteomics
    title : The title of the peak list that you want

    Returns
    -------
    float
        Pepmass value.

    """

    params = spectra.get_spectrum(title)
    params = params["params"]

    if "pepmass" in params.keys():
        mass = params["pepmass"]
        if type(mass) == float:
            return mass
        if type(mass) == tuple or type(mass) == list:
            return mass[0]
        else:
            return None
    else:
        return None


def choice_bar_at_side(mgf_file, key=None):
    """


    Parameters
    ----------
    mgf_file : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    if mgf_file is not None:
        # st.sidebar.write("Select your spectra")
        spectra = mgf.IndexedMGF(mgf_file)

        # List all spectra to choose
        all_title = get_all_titles(spectra)
        idx_sel = st.sidebar.selectbox("Select your spectra", all_title, key=key)

        return pyteomics_2_matchms(spectra, idx_sel)
    else:
        return text_to_frags_spec("0 0")


def show_meta_details(spectra, check_key='meta_key_A'):
    # st.write("Spectrum parameters")
    st.write(spectra.metadata)
    # Show molecule if has smiles
    if 'smiles' in spectra.metadata:
        smi_str = 'smiles'
    elif 'SMILES' in spectra.metadata:
        smi_str = 'SMILES'
    elif 'Smiles' in spectra.metadata:
        smi_str = 'Smiles'
    else:
        smi_str = None

    if smi_str is not None:
        smiles = spectra.metadata.get(smi_str)
        st.text("Smiles Detected")
        st.text(f"\t{smiles}")

        show_smi = st.checkbox("Show smiles structure", key=check_key)

        if show_smi:
            try:
                if smiles is not None:
                    from rdkit import Chem
                    from rdkit.Chem import Draw

                    # Draw
                    m = Chem.MolFromSmiles(smiles)
                    st.image(Draw.MolToImage(m))
            except Exception:
                st.write("Sorry, something went wrong")
                st.write("please check if it is a valid smiles")


###################################################
#######        PEAK LISTS  FUNCTIONS        #######
###################################################


def remove_not_shifted_common_pks(peak_array, shared_list, col):
    """
    Parameters
    ----------
    peak_array: np 2d array
        peak array from matchms Spikes function
    shared_list: np array
        from matchms collect_peak_pairs function
    col: int
        col

    Returns
    -------
    np 2d peak array peaks remained from filter

    """
    # Old and new Shapes
    old_shape_x = peak_array.shape[0]
    flat_idx = list(2 * shared_list[:, col]) + list((2 * shared_list[:, col]) + 1)
    list_not_shared_idx = [i for i in range(old_shape_x * 2) if i not in flat_idx]

    peaks_output = np.take(peak_array, list_not_shared_idx).reshape(-1, 2)
    return peaks_output


# Initiate Pages
class Page(ABC):
    @abstractmethod
    def write(self):
        pass
