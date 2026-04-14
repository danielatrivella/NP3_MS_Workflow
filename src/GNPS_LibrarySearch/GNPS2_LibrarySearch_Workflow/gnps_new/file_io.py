import os

import numpy as np
from pyteomics import mzml, mzxml, mgf

from _utils import Spectrum


def batch_process_queries(qry_file, min_peak, batch_size=1000):
    """Process query spectra in batches to reduce memory usage"""
    file_format = os.path.splitext(qry_file)[1].lower()
    batch_specs = []
    batch_prec_mzs = []

    if file_format == '.mzml':
        with mzml.MzML(qry_file) as reader:
            for spectrum in reader:
                try:
                    if spectrum['ms level'] == 2:  # MS2 scans only
                        precursor_mz = float(
                            spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0][
                                'selected ion m/z'])

                        mz_array = np.array(spectrum['m/z array'])
                        intensity_array = np.array(spectrum['intensity array'])
                        tic = np.sum(intensity_array)
                        intensity_array = intensity_array / np.max(intensity_array) * 999  # will be converted to float32

                        # Skip peaks with less than min_peak peaks
                        if sum(intensity_array > 0) < min_peak:
                            continue

                        peaks = np.column_stack((mz_array, intensity_array))

                        scan_number = int(spectrum['index']) + 1

                        rt = float(spectrum['scanList']['scan'][0]['scan start time'])

                        if 'positive scan' in spectrum:
                            polarity = 1
                        elif 'negative scan' in spectrum:
                            polarity = -1
                        else:
                            polarity = 0

                        try:
                            charge = int(
                                spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0][
                                    'charge state'])
                            charge = charge * polarity
                        except:
                            charge = polarity

                        spec = Spectrum(
                            scan=scan_number,
                            precursor_mz=precursor_mz,
                            rt=rt,
                            charge=charge,
                            tic=tic,
                            peaks=peaks)

                        batch_specs.append(spec)
                        batch_prec_mzs.append(precursor_mz)

                        if len(batch_specs) >= batch_size:
                            yield batch_specs, np.array(batch_prec_mzs)
                            batch_specs = []
                            batch_prec_mzs = []
                except:
                    continue

    elif file_format == '.mzxml':
        with mzxml.MzXML(qry_file) as reader:
            for spectrum in reader:
                try:
                    if spectrum['msLevel'] == 2:  # MS2 scans only
                        precursor_mz = float(spectrum['precursorMz'][0]['precursorMz'])

                        mz_array = np.array(spectrum['m/z array'])
                        intensity_array = np.array(spectrum['intensity array'])
                        tic = np.sum(intensity_array)
                        intensity_array = intensity_array / np.max(intensity_array) * 999  # will be converted to float32

                        if sum(intensity_array > 0) < min_peak:
                            continue

                        peaks = np.column_stack((mz_array, intensity_array))

                        scan_number = int(spectrum['num'])

                        rt = float(spectrum['retentionTime'])

                        if spectrum['polarity'] == '+':
                            polarity = 1
                        elif spectrum['polarity'] == '-':
                            polarity = -1
                        else:
                            polarity = 0

                        try:
                            charge = int(spectrum['precursorMz'][0]['precursorCharge'])
                            charge = charge * polarity
                        except:
                            charge = polarity

                        spec = Spectrum(
                            scan=scan_number,
                            precursor_mz=precursor_mz,
                            rt=rt,
                            charge=charge,
                            tic=tic,
                            peaks=peaks)

                        batch_specs.append(spec)
                        batch_prec_mzs.append(precursor_mz)

                        if len(batch_specs) >= batch_size:
                            yield batch_specs, np.array(batch_prec_mzs)
                            batch_specs = []
                            batch_prec_mzs = []
                except:
                    continue

    elif file_format == '.mgf':
        scan_idx = 0
        with mgf.MGF(qry_file, convert_arrays=0, read_charges=False) as reader:
            for spectrum in reader:
                try:
                    scan_idx += 1
                    d = spectrum['params']

                    try:
                        if 'precursor_mz' not in d:
                            precursor_mz = float(d['pepmass'][0])
                        else:
                            precursor_mz = float(d['precursor_mz'])
                    except:
                        continue

                    mz_array = np.array(spectrum['m/z array'])
                    intensity_array = np.array(spectrum['intensity array'])
                    tic = np.sum(intensity_array)
                    intensity_array = intensity_array / np.max(intensity_array) * 999  # will be converted to float32

                    if sum(intensity_array > 0) < min_peak:
                        continue

                    peaks = np.column_stack((mz_array, intensity_array))

                    if 'rtinseconds' in d:
                        rt = float(d['rtinseconds'])
                    elif 'rt' in d:
                        rt = float(d['rt'])
                    else:
                        rt = 0

                    if 'charge' in d:
                        charge = d['charge'][0]
                    else:
                        charge = 0

                    if 'scans' in d:
                        scan = d['scans']
                    else:
                        scan = scan_idx

                    spec = Spectrum(
                        scan=scan,
                        precursor_mz=precursor_mz,
                        rt=rt,
                        charge=charge,
                        tic=tic,
                        peaks=peaks)

                    batch_specs.append(spec)
                    batch_prec_mzs.append(precursor_mz)

                    if len(batch_specs) >= batch_size:
                        yield batch_specs, np.array(batch_prec_mzs)
                        batch_specs = []
                        batch_prec_mzs = []
                except:
                    continue
    else:
        raise ValueError(f"Unsupported file format: {file_format}")

    if batch_specs:  # yield remaining spectra
        yield batch_specs, np.array(batch_prec_mzs)


def read_mgf_spectrum(file_obj):
    """Read a single spectrum block from an open MGF file.
    Args:
        file_obj: An opened file object positioned at the start of a spectrum
    Returns:
        dict: Spectrum information containing all metadata and peaks,
        or None if end of file is reached
    """
    # Initialize spectrum with common metadata fields
    spectrum = {
        'PEPMASS': 0.0,
        # 'CHARGE': '',
        # 'MSLEVEL': '',
        # 'SOURCE_INSTRUMENT': '',
        # 'FILENAME': '',
        # 'SEQ': '',
        # 'IONMODE': '',
        # 'ORGANISM': '',
        # 'NAME': '',
        # 'PI': '',
        # 'DATACOLLECTOR': '',
        # 'SMILES': '',
        # 'INCHI': '',
        # 'INCHIAUX': '',
        # 'PUBMED': '',
        # 'SUBMITUSER': '',
        # 'LIBRARYQUALITY': '',
        'SPECTRUMID': '',
        'SCANS': '',
        'peaks': []
    }

    # Skip any empty lines before BEGIN IONS
    for line in file_obj:
        if line.strip() == 'BEGIN IONS':
            break
    else:  # EOF reached
        return None

    # Read spectrum metadata and peaks
    for line in file_obj:
        line = line.strip()

        if not line:  # Skip empty lines
            continue

        if line == 'END IONS':
            if spectrum['peaks']:
                this_peaks = np.asarray(spectrum['peaks'])
                this_peaks[:, 1] = this_peaks[:, 1] / np.max(this_peaks[:, 1]) * 999
                this_peaks = this_peaks[np.bitwise_and(this_peaks[:, 0] > 0, this_peaks[:, 1] > 0)]
                spectrum['peaks'] = np.asarray(this_peaks, dtype=np.float32)

            return spectrum

        # Handle peak data
        if line and not line.startswith(('BEGIN', 'END')):
            if '=' in line:
                # First try to parse as metadata
                key, value = line.split('=', 1)

                # Handle specific numeric fields
                if key in ['PEPMASS', 'PRECURSORMZ']:
                    try:
                        spectrum[key] = float(value.strip())
                    except:
                        continue
                elif key in spectrum:  # Store other known metadata fields as strings
                    spectrum[key] = value
            else:
                # If no '=' found, treat as peak data
                try:
                    mz, intensity = line.split()
                    spectrum['peaks'].append((float(mz), float(intensity)))
                except:
                    continue

    return None


def iterate_gnps_lib_mgf(mgf_path, buffer_size=8192):
    """Iterate through spectra in an MGF file efficiently using buffered reading.
    Args:
        mgf_path: Path to the MGF file
        buffer_size: Read buffer size in bytes
    Yields:
        dict: Spectrum information containing all metadata and peaks
    """
    with open(mgf_path, 'r', buffering=buffer_size) as f:
        while True:
            spectrum = read_mgf_spectrum(f)
            if spectrum is None:
                break
            yield spectrum


if __name__ == "__main__":

    for batch_specs, prec_mzs in batch_process_queries('/Users/shipei/Documents/test_data/mgf/input_file_larger.mgf', 1, 20):
        print(len(batch_specs), prec_mzs)
        break
