#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 15:23:19 2021

@author: np3
"""

import pandas as pd
import numpy as np
import termplotlib as tpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import shapiro, anderson, norm, normaltest

def treat_dfMS1_list_pw(dfMS1_list_with_MS2, dfMeta):
    '''
    Computes the peak widths of the MS1 peaks.
    Exclude blanks samples using the metadata information

    Parameters
    ----------
    dfMS1_list_with_MS2 : pandas dataframe
        MS1 table.
    dfMeta : pandas dataframe
        Metadata table.

    Returns
    -------
    pandas dataframe
        Table with mz, rtMean and peak width as pw.
    '''
    print("\n* Computing the peak width of the MS1 list with MS2 and removing blanks *\n")
    n = dfMS1_list_with_MS2.shape[0]
    # transform the column names and the samples types to lower case (no case sensitive)
    dfMeta.columns = dfMeta.columns.str.lower()
    dfMeta.columns = dfMeta.columns.str.strip()
    dfMeta["sample_type"] = dfMeta.sample_type.str.lower()    
    # Get blank samples names
    blank_sample_codes = dfMeta[dfMeta.sample_type == "blank"]["sample_code"].to_list()
    # get the blank samples columns names
    blank_cols = [blank_sample_code + '_area' for blank_sample_code in blank_sample_codes]
    # remove MS1 peaks from blank samples
    blank_peaks = (dfMS1_list_with_MS2[blank_cols].sum(1) > 0)
    n_blank_peaks_rm = blank_peaks.sum()
    dfMS1_list_with_MS2 = dfMS1_list_with_MS2[~blank_peaks].reset_index(drop=True)
    # remove peaks with a mz close to a blank mz which BFLAG is TRUE
    n_blank_peaks_rm += dfMS1_list_with_MS2.BFLAG.sum()
    dfMS1_list_with_MS2 = dfMS1_list_with_MS2[(~dfMS1_list_with_MS2.BFLAG)].reset_index(drop=True)
    # peak_width calc
    dfMS1_list_with_MS2["pw"] = (dfMS1_list_with_MS2.rtMax - dfMS1_list_with_MS2.rtMin)
    print("- Removed a total of", n_blank_peaks_rm, "MS1 peaks (%.1f%%) of blank m/zs from the list" % (n_blank_peaks_rm/n*100))
    # return the MS1 peaks mz, rtMean and peak width
    return dfMS1_list_with_MS2[["mz", "rtMean", "pw"]]

# tests if the peak width distribution is normal,
# uses the Shapiro-Wilk and the anderson-darling methods,
# if yes in more than half of the tests, use the confidence interval of 95% for the pw range and the mean for the rt_tol
# else use the median for the rt_tol and the median +- 1.5 IQR for the pw range
def pw_min_max_rtTol_suggestions(pw):
    norm_tests_result = np.array([False]*3)
    print("\n===============")
    print("Normality check: Null hypothesis, or H0, the peak width distribution is Normal - alpha = 5%")
    print("===============\n")
    # only chck for normality if enough data
    if len(pw) > 50:
        # interpret shapiro
        alpha = 0.05
        shapiro_tresult = shapiro(pw)
        print("\n  Shapiro-Wilk Test statistic = %0.3f" % shapiro_tresult.statistic)
        if shapiro_tresult.pvalue > alpha:
            norm_tests_result[0] = True
            print('    p-value = %.2f : data looks normal (fail to reject H0)' % shapiro_tresult.pvalue)
        else:
            print('    p-value = %.2f : data does not look normal (reject H0)' % shapiro_tresult.pvalue)
        # interpret normal test
        normaltest_result = normaltest(pw)
        print("  D'Agostino and Pearson Test statistic = %0.3f" % normaltest_result.statistic)
        if normaltest_result.pvalue > alpha:
            norm_tests_result[1] = True
            print('    p-value = %.2f : data looks normal (fail to reject H0)' % normaltest_result.pvalue)
        else:
            print('    p-value = %.2f : data does not look normal (reject H0)' % normaltest_result.pvalue)
        # interpret anderson
        anderson_dresult = anderson(pw)
        print("  Anderson-Darling Test statistic = %0.3f" % anderson_dresult.statistic)
        alpha_5_index = 2
        sl, cv = anderson_dresult.significance_level[alpha_5_index], anderson_dresult.critical_values[alpha_5_index]
        if anderson_dresult.statistic < anderson_dresult.critical_values[alpha_5_index]:
            norm_tests_result[2] = True
            print('    %.3f: %.3f, data looks normal (fail to reject H0)' % (sl, cv))
        else:
            print('    %.3f: %.3f, data does not look normal (reject H0)' % (sl, cv))
    # final assumption
    if norm_tests_result.sum() > 0:
        print("\n===============")
        print("Failed to reject H0 in at least one of the normality tests.",
              "Let's assume that the peak width distribution is Normal.")
        print("===============")
        print("\nComputing the parameters suggestions as:\n",
              "- Retention time tolerance = mean / 2\n",
              "- Peak width range (minimum and maximum) = 95% confidence interval")
        print("   - The peak width suggestions must be between the Q2.5 and the Q97.5 values,",
              "or they will receive these limit values.")
        # assuming that this list of peak widths are a sample and we want to estimate the population mean and
        # # standard deviation, we will use
        mean, sigma = pw.mean(), pw.std(ddof=1)
        pw_min, pw_max = norm.interval(0.95, loc=mean, scale=sigma)
        pw_max = np.round(min(pw_max, np.quantile(pw, 0.975)), 1)
        pw_min = np.round(max(pw_min, np.quantile(pw, 0.025)), 1)
        rt_tol = mean / 2
    else:
        print("\n===============")
        print("Rejected H0 in all of the normality tests.",
              "Let's NOT assume that the peak width distribution is Normal.")
        print("===============")
        rt_tol = np.median(pw) / 2
        q25_pw = np.quantile(pw, 0.25)
        q75_pw = np.quantile(pw, 0.75)
        iqr = q75_pw - q25_pw
        iqr_min_factor = 0
        pw_limit = min(pw)
        # select the iqr_min_factor until the pw_min is above the minimum value
        iqr_factor_range = [1.5, 1.25, 1, 0.75, 0.5, 0.25, 0]
        for iqr_factor in iqr_factor_range:
            if np.round(q25_pw - iqr * iqr_factor, 1) >= pw_limit:
                iqr_min_factor = iqr_factor
                break
        iqr_max_factor = 0
        pw_limit = max(pw)
        # select the iqr_max_factor until the pw_max is below the maximum value
        for iqr_factor in iqr_factor_range:
            if np.round(q75_pw + iqr * iqr_factor, 1) <= pw_limit:
                iqr_max_factor = iqr_factor
                break
        pw_min = np.round(q25_pw - iqr * iqr_min_factor, 1)
        pw_max = np.round(q75_pw + iqr * iqr_max_factor, 1)
        print("\nComputing the parameters suggestions as:\n",
              "- Retention time tolerance = median / 2\n",
              "- Peak width range (minimum and maximum) = Q25 -",str(iqr_min_factor),
              "* IQR, Q75 +",str(iqr_max_factor),"* IQR")
        print("   - The IQR (interquartile range) factors start at 1.5 and are decreased by 0.25 until the peak width suggestions are",
              "between the minimum and maximum peak width limit values.")
    # return the suggestions of pw and rt tol
    return ((pw_min, pw_max), rt_tol)

def pw_hist_plot(dfpw, n_bins, counts, bin_edges, plot_filename="HistPlot_PW.png"):
    sns.set(rc={'figure.figsize':(15,25)},
            font_scale=1.7)
    sns.set_theme(style="whitegrid", palette="muted")
    fig = sns.histplot(dfpw,bins=n_bins, y="pw")
    # add the bins pw interval
    for i in range(n_bins):
        if i+1 < n_bins:
            txt = "{} - [{},{}[".format(counts[i],
                                        round(bin_edges[i],1),
                                        round(bin_edges[i+1],1))
        else:
            txt = "{} - [{},{}]".format(counts[i],
                                        round(bin_edges[i],1),
                                        round(bin_edges[i+1],1))
        plt.text(x=counts[i],
                 y = (1)*(bin_edges[i]+bin_edges[i+1])/2,
                 fontsize="medium",
                 va="center",
                 s=txt)
    # plot vertical lines in the suggestions in regions of interest
    pw = dfpw.pw.values
    mean, sigma = pw.mean(), pw.std(ddof=1)
    pw_min, pw_max = norm.interval(0.95, loc=mean, scale=sigma)
    pw_max = min(pw_max, np.quantile(pw, 0.975))
    pw_min = max(pw_min, np.quantile(pw, 0.025))
    median = np.median(pw)
    plt.axhline(y=pw_min, linestyle='-',label='CI 95% min', color='skyblue')
    plt.axhline(y=pw_max, linestyle='-', label='CI 95% max', color='blue')
    plt.axhline(y=mean, linestyle='-',label='pw mean', color='gray')
    plt.axhline(y=median, linestyle='-', label='pw median', color='black')
    # compute the iqr factor
    q25_pw = np.quantile(pw, 0.25)
    q75_pw = np.quantile(pw, 0.75)
    iqr = q75_pw - q25_pw
    iqr_factor_range = [1.5, 1.25, 1, 0.75, 0.5, 0.25, 0]
    pw_limit = min(pw)
    iqr_min_factor = 0
    for iqr_factor in iqr_factor_range:
        if np.round(q25_pw - iqr * iqr_factor, 1) >= pw_limit:
            iqr_min_factor = iqr_factor
            break
    iqr_max_factor = 0
    pw_limit = max(pw)
    # select the iqr_max_factor until the pw_max is below the maximum value
    for iqr_factor in iqr_factor_range:
        if np.round(q75_pw + iqr * iqr_factor, 1) <= pw_limit:
            iqr_max_factor = iqr_factor
            break
    pw_min = q25_pw - iqr * iqr_min_factor
    pw_max = q75_pw + iqr * iqr_max_factor
    plt.axhline(y=pw_min, linestyle='-', label='pw Q25 - IQR * '+str(iqr_min_factor), color='pink')
    plt.axhline(y=pw_max, linestyle='-', label='pw Q75 + IQR * '+str(iqr_max_factor), color='magenta')
    plt.axhline(y=q25_pw, linestyle='-', label='pw Q25', color='green')
    plt.axhline(y=q75_pw, linestyle='-', label='pw Q75', color='yellow')
    plt.legend()
    # Add title and x limits
    plt.title("Peak Width distribution", fontsize="large")
    plt.xlim(0, max(counts)*1.2)
    fig.figure.savefig(plot_filename, dpi=150)

def pw_hist_plot_terminal(counts, bin_edges):
    print("\n* The peak width histogram for the MS1 list with MS2 and without blanks *\n")
    fig = tpl.figure()
    fig.hist(counts,
              list(bin_edges),
              orientation="horizontal",
              force_ascii=False)
    fig.show()

def plot_pw_hist_suggest_parms(pathMeta, pathMsYes, n_bins, plot_filename):
    dfMeta = pd.read_csv(pathMeta, low_memory=False)
    dfMS1_list_with_MS2 = pd.read_csv(pathMsYes, low_memory=False)

    # Png name verification
    if not plot_filename.endswith(".png"):
        plot_filename += ".png"
    # add the pw, remove blanks and get the pw histogram
    dfMS1_list_with_MS2Treated = treat_dfMS1_list_pw(dfMS1_list_with_MS2, dfMeta)
    del dfMS1_list_with_MS2
    if dfMS1_list_with_MS2Treated.shape[0] > 0:
        counts, bin_edges = np.histogram(dfMS1_list_with_MS2Treated.pw, bins=n_bins)

        pw_hist_plot(dfpw=dfMS1_list_with_MS2Treated,
                  n_bins=n_bins,
                  counts=counts,
                  bin_edges=bin_edges,
                  plot_filename=plot_filename)
        pw_hist_plot_terminal(counts=counts,
                           bin_edges=bin_edges)

        print("\n* The peak width descriptive statistics for the MS1 list with MS2 and without blanks *\n")
        print(dfMS1_list_with_MS2Treated.pw.describe())
        pw_range, rt_tol = pw_min_max_rtTol_suggestions(dfMS1_list_with_MS2Treated.pw.values)
        print("\n===============")
        print("Suggestion of values for the MS1 processing parameters:")
        print("  -e, --peak_width [X,Y] = [%.1f,%.1f] seconds" % pw_range)
        print("  -t, --rt_tolerance [X] = %.1f seconds" % rt_tol)
        print("===============")
        print("We recommend the user to analyse the pre process final diagnostic result before running the job with",
              "the suggested parameters values. If the rates of no correspondence between the MS2 spectra and the MS1",
              "peaks by sample are above the cutoff, the user must evaluate the problem and pay attention to the",
              "respective rates descriptive statistics values (mean, median, etc) to decide if it is really an issue.")
        print("===============")
        print("The user must take the following relations into account when choosing the parameters values:",
              "\n  - Minimum Peak width: ",
              "\n    - A small value can split short MS1 peaks",
              "\n    - A big value can exclude very short peaks and/or join adjacent isomers",
              "\n  - Maximum Peak width: ",
              "\n    - A small value can split large MS1 peaks",
              "\n    - A big value can join adjacent isomers",
              "\n  - Retention time tolerance: ",
              "\n    - A small value can prevent the match between a MS2 spectrum and its respective MS1 peak",
              "\n    - A big value can assign a MS2 spectrum to a wrong MS1 peak, which can hide a bad MS1",
              "integration processing")
        print("===============\n")
    else:
        print("* There is no MS1 peak in the MS1 list with MS2 from a not blank sample - pre process parameters",
              "suggestion aborted *")


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Wrong number of arguments. Please insert the following args:\n",
              "\t1. <metadata_path>\n",
              "\t2. <table_MS1_with_ms2_path>\n",
              "\t3. [output_png_name] default: MS1_list_with_MS2_noBlank_peak_width_hist.png\n",
              "\t4. [hitogram_number_bins] default: 150")
        sys.exit()
    elif len(sys.argv) < 4:
        pathMeta = sys.argv[1]
        pathMsYes = sys.argv[2]
        plot_filename = "MS1_list_with_MS2_noBlank_peak_width_hist"
        n_bins = 150
    elif len(sys.argv) < 5:
        pathMeta = sys.argv[1]
        pathMsYes = sys.argv[2]
        plot_filename = sys.argv[3]
        n_bins = 150
    else:
        pathMeta = sys.argv[1]
        pathMsYes = sys.argv[2]
        plot_filename = sys.argv[3]
        n_bins = int(sys.argv[4])

    plot_pw_hist_suggest_parms(pathMeta, pathMsYes, n_bins, plot_filename)

