import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from os.path import dirname

def plot_basePeakInt_blank_distribution(clustering_counts_path, bflag_cutoff_factor):
    output_path = dirname(clustering_counts_path)

    # read the clustering counts and check if there is a blank sample
    clustering_counts = pd.read_csv(clustering_counts_path)
    # set plot parameters
    sns.set_theme(style="darkgrid")
    sns.set_context("paper", rc={"font.size": 18, "axes.titlesize": 18, "axes.labelsize": 16, "legend.fontsize": 14,
                                 "xtick.labelsize": 12, "ytick.labelsize": 12})
    #print("\n* Plotting the base peak intensity distribution *\n")
    # only plot the distribution if there is a blank sample
    if clustering_counts.columns.str.match('^BLANKS_TOTAL$').any():
        clustering_counts['BLANKS'] = (clustering_counts.BLANKS_TOTAL > 0)
        basePeakInt_blank_IQR = (clustering_counts.loc[clustering_counts.BLANKS, 'basePeakInt'].quantile(0.75) -
                                 clustering_counts.loc[clustering_counts.BLANKS, 'basePeakInt'].quantile(0.25))

        descriptive_output = "\n* Descriptive analysis of the basePeakInt distribution for BLANKS *\n\n" + \
                             clustering_counts.loc[clustering_counts.BLANKS, 'basePeakInt'].describe().to_string() + \
                             "\n\n- Interquartile rank (IQR) = " + str(basePeakInt_blank_IQR)

        # bflag_cutoff_factor equals zero is the same as the median;
        # also the bflag_cutoff_factor can be zero when the cutoff is disabled
        if bflag_cutoff_factor > 0:
            descriptive_output += "\n- Median + " + str(bflag_cutoff_factor) + " * IQR = "+\
                                  str(clustering_counts.loc[clustering_counts.BLANKS,'basePeakInt'].median() +
                                      bflag_cutoff_factor * basePeakInt_blank_IQR)
        else:
            # print the default cutoff value to ilustrate the computation
            descriptive_output += "\n- Median + " + str(1.5) + " * IQR = " + \
                                  str(clustering_counts.loc[clustering_counts.BLANKS, 'basePeakInt'].median() +
                                      1.5 * basePeakInt_blank_IQR)

        # print and save the basePeakInt summary
        print(descriptive_output)
        with open(output_path + '/basePeakInt_distribution_summary.txt', "w") as text_file:
            text_file.write(descriptive_output)

        # plot the distribution of the basePeakInt entries in log scale
        basePeakInt_displot = sns.displot(data=clustering_counts,
                    x="basePeakInt", multiple="stack", hue="BLANKS", col='BFLAG', log_scale=True, height=10,
                    aspect=1.2)

        # set the bflag factors vertical lines and color, skip 0 factor which is equals de median value
        bflag_cutoff_factors_axvlines = [1.5, 3, 5, 10, 15, 20]
        if bflag_cutoff_factor not in bflag_cutoff_factors_axvlines and not bflag_cutoff_factor == 0:
            bflag_cutoff_factors_axvlines.append(bflag_cutoff_factor)
            bflag_cutoff_factors_axvlines.sort()
        bflag_cutoff_factors_axvlines_color = [(bflag_cutoff_factors_axvlines[i], cl)
                              for i, cl in enumerate(sns.color_palette("husl", len(bflag_cutoff_factors_axvlines)))]

        basePeakInt_blank_median = clustering_counts.loc[clustering_counts.BLANKS, 'basePeakInt'].median()
        # add vertical lines to both grids with the describe statistics used and some bflag factors
        for ax in basePeakInt_displot.axes.flatten():
            ax.axvline(x=clustering_counts.loc[clustering_counts.BLANKS, 'basePeakInt'].mean(), linestyle='-',
                        color='gray',
                        label='mean basePeakInt of BLANKS')
            ax.axvline(x=basePeakInt_blank_median, linestyle='-',
                        color='black',
                        label='median basePeakInt of BLANKS')
            ax.axvline(x=clustering_counts.loc[clustering_counts.BLANKS,'basePeakInt'].quantile(0.25), linestyle='-',
                        color='brown',label='q25th basePeakInt of BLANKS')
            ax.axvline(x=clustering_counts.loc[clustering_counts.BLANKS,'basePeakInt'].quantile(0.75), linestyle='-',
                        color='yellow',label='q75th basePeakInt of BLANKS')
            for bflag_factor, color_line in bflag_cutoff_factors_axvlines_color:
                ax.axvline(x=basePeakInt_blank_median + bflag_factor * basePeakInt_blank_IQR,
                           linestyle='-', color=color_line,
                           label='median + ' + str(bflag_factor) + ' IQR basePeakInt of BLANKS')
            # ax.xaxis.set_tick_params(labelrotation=-30)
    else:
        #  write file if IQR of the hole distribution
        basePeakInt_IQR = (clustering_counts.loc[:, 'basePeakInt'].quantile(0.75) -
                                 clustering_counts.loc[:, 'basePeakInt'].quantile(0.25))

        descriptive_output = "\n* Descriptive analysis of the complete basePeakInt distribution *\n\n" + \
                             clustering_counts.loc[:, 'basePeakInt'].describe().to_string() + \
                             "\n\n- Interquartile rank (IQR) = " + str(basePeakInt_IQR)
        # add the default noise cutoff factor
        descriptive_output += "\n- Median + " + str(1.5) + " * IQR = " + \
                                  str(clustering_counts.loc[:, 'basePeakInt'].median() +
                                      1.5 * basePeakInt_IQR)
        # print and save the basePeakInt summary
        print(descriptive_output)
        with open(output_path + '/basePeakInt_distribution_summary.txt', "w") as text_file:
            text_file.write(descriptive_output)

        # plot the distribution of all spectra in log scale
        basePeakInt_displot = sns.displot(data=clustering_counts,
                    x="basePeakInt", multiple="stack", log_scale=True, height=10, aspect=1.2)

        # set the noise cutoff factors vertical lines and color, skip 0 factor which is equals de median value
        noise_cutoff_factors_axvlines = [1.5, 3, 5, 10, 15, 20]
        if bflag_cutoff_factor not in noise_cutoff_factors_axvlines and not bflag_cutoff_factor == 0:
            noise_cutoff_factors_axvlines.append(bflag_cutoff_factor)
            noise_cutoff_factors_axvlines.sort()
        noise_cutoff_factors_axvlines_color = [(noise_cutoff_factors_axvlines[i], cl)
                                               for i, cl in
                                               enumerate(sns.color_palette("husl", len(noise_cutoff_factors_axvlines)))]

        basePeakInt_median = clustering_counts.loc[:, 'basePeakInt'].median()
        for ax in basePeakInt_displot.axes.flatten():
            ax.axvline(x=clustering_counts['basePeakInt'].mean(), linestyle='-',
                        color='gray',
                        label='mean basePeakInt')
            ax.axvline(x=clustering_counts['basePeakInt'].median(), linestyle='-',
                        color='black',
                        label='median basePeakInt')
            # add the other vertical lines
            ax.axvline(x=clustering_counts['basePeakInt'].quantile(0.25), linestyle='-',
                       color='brown', label='q25th basePeakInt')
            ax.axvline(x=clustering_counts['basePeakInt'].quantile(0.75), linestyle='-',
                       color='yellow', label='q75th basePeakInt')
            for noise_factor, color_line in noise_cutoff_factors_axvlines_color:
                ax.axvline(x=basePeakInt_median + noise_factor * basePeakInt_IQR,
                           linestyle='-', color=color_line,
                           label='median + ' + str(noise_factor) + ' IQR basePeakInt')
            # ax.xaxis.set_tick_params(labelrotation=-30)

    basePeakInt_displot.fig.subplots_adjust(top=0.9)
    basePeakInt_displot.fig.suptitle("Base Peak Intensity Distribution of the Clustering Counts Consensus Spectra")
    plt.legend()
    plt.savefig(output_path + '/basePeakInt_distribution.png', bbox_inches='tight', dpi=100)

    print("\nBase peak intensity distribution saved to basePeakInt_distribution.png file!\nDONE!\n")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        # print(sys.argv)
        clustering_counts_path = sys.argv[1]
        bflag_cutoff_factor = float(sys.argv[2])
    else:
        print(
            "Error: Two arguments must be supplied to plot the basePeakInt distribution showing the blanks and bflag_cutoff:\n"
            "1. The path to the clustering counts table;\n"
            "2. The bflag_cutoff factor that will be applied.\n")
        sys.exit(1)

    plot_basePeakInt_blank_distribution(clustering_counts_path, bflag_cutoff_factor)
