
import sys
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import argparse

# superclass fixed colors

superclass_colors = {'Alkaloids_and_Lactams': '#ff8b00',                    # orange
	                'Aminoacids_and_Peptides_and_OtherNComp': '#e8ff00',   # yellow
	                'Benzenoids': '#5dff00',                               # light green
	                'Fatty_Acids_and_Lipids': '#00cc00',                   # dark green
	                'Flavonoids_and_Phenolic_derivatives': '#6fffff',      # cyan
	                'Lignans_and_Other_Ocompounds': '#00b9ff',             # blue
	                'Organic_Acids_and_OthersGenerals': '#002eff',         # dark blue
	                'Polyketides': '#ff008b',                              # pink
	                'Terpenes_and_Carotenoids': '#cc0000',                 # dark red
	                'Not_Annotated': "#cccccc"}                            # grey

# save a donut plot
# where values is the list of numeric values to plot, the categories is the list of their labels
# and colors is the list with the colors of each value
# donutplot_filepath is the filepath to where the image will be stored and saved, with .png extension
# Create the pie chart with a hollow center (width=0.4 creates the donut ring)
def plot_donut_values(donutplot_filepath, values, categories, title, colors = ["blue", "green"],
                      title_size=16, text_size=14):
	total = sum(values)
	plt.pie(values, labels=categories, colors=colors,
	        autopct=lambda pct: f'{int(round(pct * total / 100))}\n({pct:.1f}%)',
	        startangle=90, pctdistance=0.70, wedgeprops=dict(width=0.6, edgecolor='white'),
	        textprops={'fontsize': text_size})
	plt.axis('equal') # Equal aspect ratio ensures that pie is drawn as a circle
	plt.title(title, fontsize=title_size, pad=30)
	plt.tight_layout()
	plt.savefig(donutplot_filepath, dpi=300, bbox_inches='tight')
	plt.close()
#plt.show()

# save a stacked bar from a pandas dataframe
# where df is the dataframe used to plot a stackedbar, where the  rows will be the x values and
# the columns will be the y values stacked - the column names will be the legend
# the title, xlabel and y label of the plot must be given
# stackedbarplot_filepath is the filepath to where the image will be stored and saved, with .png extension
def plot_stacked_bar_pandas_df(stackedbarplot_filepath, df, title, xlabel, ylabel, figsize=(18,8),
                               label_size=13, title_size=20, legend_title = None,
                               legend_bbox_to_anchor=(0.5, -0.3), legend_fontsize=17, legend_ncol=4, colors=None):
	ax = df.plot(kind='bar', stacked=True, figsize=figsize,color=colors)
	ax.legend([col.replace("_", " ") for col in df.columns.values],
	          title = legend_title, title_fontsize=legend_fontsize,
	          bbox_to_anchor=legend_bbox_to_anchor, loc='upper center', ncol=legend_ncol, fontsize=legend_fontsize)
	plt.title(title, fontsize=title_size, pad=20)
	plt.xlabel(xlabel, fontsize=label_size, labelpad=5)
	plt.ylabel(ylabel, fontsize=label_size)
	plt.xticks(ticks=range(df.shape[0]),
	           labels=[str(i).replace('_area', '') for i in df.index],
	           rotation=45, ha='right', fontsize=label_size)
	plt.tight_layout()
	plt.savefig(stackedbarplot_filepath, dpi=300, bbox_inches='tight')
	plt.close()
# plt.show()

# post analysis routine for drug discovery research
# plot the m/z occurrence and superclass grouping by sample and save them to the output_path together with
# the respective tables with one sample by row;
# use all samples to compute the metrics (except the ones selected for removal)
# # and only use the samples present in the metadata for the plots and for selecting the topk novelty samples
# where
# # metadata_path is the path to the metadata of the job to be post analysed, it may contain the complete metadata or just a selection of the original samples, which will be used to filter the final plots and created tables by sample
# # clean_counts_path is the path to the clean counts table of the job to be post analysed, same of the metadata
# # output_path Path to the output directory where the plots and tables of the post analysis will be stored
# # topk is the number of samples to select and filter in the plots with top values of exclusive and not annotated m/z - top novelty
# # rm_blanks, rm_beds and rm_controls receive a True or False to allow removing blanks, culture media or control samples and m/z from the metrics computation
# # use_protonated a True of False defining if only the protonated m/z should be used in the plots (filter the table with protonated_representative==1).
# # plots args for sizing personalization and legend size and position:
# # mzs_barplot (exclusive and redudant m/z distribution with library annotation by sample),
# # donutplots of library annotation proportion and exclusive/redundant m/z proportion in the selected samples
# # mzs_barplot_legend_bbox Sets the anchoring coordinates relative to the plot area. (0, 0) is the bottom-left corner of the plot. (1, 1) is the top-right corner of the plot. Values smaller than 0 or greater than 1 will place the legend completely outside the plot area.
# # mzs_barplot_colors colors for each sample or None to use default coloring
def post_dd_analysis_plots(metadata_path, clean_counts_path, output_path, topk = None, rm_blanks = True,
                           rm_beds = True, rm_controls = True, use_protonated=False, donutplots_title_size=16,
                           donutplots_text_size=14, donutplot_libAnnotations_colors=['#ff8b00', '#6372b4', "#c6c6c6"],
						   donutplot_mzs_distr_colors=["#0072c3", "#42be65"],
                           mzs_barplot_figsize=(18, 8), mzs_barplot_label_size=13, mzs_barplot_title_size=20,
                           mzs_barplot_legend_bbox=(0.5, -0.2), mzs_barplot_legend_fontsize=17,
                           mzs_barplot_legend_ncol=4, mzs_barplot_colors=['#0072c3', '#FF8C00', '#42be65', '#FDDA0D'],
                           superclass_barplot_figsize=(20, 10), superclass_barplot_label_size=16,
                           superclass_barplot_title_size=20, superclass_barplot_legend_bbox=(0.5, -0.25),
                           superclass_barplot_legend_fontsize=15, superclass_barplot_legend_ncol=4):
	print("* Post Processing Drug Discovery analysis for NP3 results *\n")
	# check the paths
	clean_counts_path = Path(clean_counts_path)
	metadata_path = Path(metadata_path)
	output_path = Path(output_path)
	if not clean_counts_path.exists() or not clean_counts_path.is_file():
		sys.exit("The provided path to the clean data file does not exists. Post processing drug discovery analysis aborted.")
	if not metadata_path.exists() or not metadata_path.is_file():
		sys.exit("The provided path to the metadata table file does not exists. Post processing drug discovery analysis aborted.")
	if not output_path.exists() or not output_path.is_dir():
		print("The provided output path directory does not exists. It will be created:", output_path)
		output_path.mkdir(exist_ok=True, parents=True)
		if not output_path.exists() or not output_path.is_dir():
			sys.exit("The provided output path could not be created and does not exists. Post processing drug discovery analysis aborted.")
	
	# extract the output name from the clean table and concatenate it with the metadata name
	output_name = clean_counts_path.name.replace("_peak_area_clean_ann.csv", "")
	output_name = output_name+"_"+metadata_path.name.replace(".csv", "")
	
	clean_counts_df = pd.read_csv(clean_counts_path, low_memory=False)
	n = clean_counts_df.shape[0]
	metadata_df = pd.read_csv(metadata_path)
	# fix the metadata column to upper and the sample types to lower
	metadata_df.columns = metadata_df.columns.str.upper()
	metadata_df["SAMPLE_TYPE"] = metadata_df.SAMPLE_TYPE.str.lower()
	# if no topk was informed, set it to the size of the provided metadata
	if topk is None or topk <= 0 or topk > metadata_df.shape[0]:
		topk = metadata_df.shape[0]
		print("- Top k novelty samples set to the number of samples in the provided metadata table = ", topk,
		      "\n  - Select the k samples containing more exclusive and not annotated m/z for plotting")
	
	# check the colors parms
	if mzs_barplot_colors is not None and mzs_barplot_colors != [None] and len(mzs_barplot_colors) < 4:
		print("WARNING: the provided list of colors in the mzs_barplot_colors parameter has a length smaller than 4 and thus will not be used. The default coloring will be used instead. ")
		mzs_barplot_colors = ['#0072c3', '#FF8C00', '#42be65', '#FDDA0D']
	if donutplot_libAnnotations_colors is not None and  donutplot_libAnnotations_colors != [None] and len(donutplot_libAnnotations_colors) < 3:
		print("WARNING: the provided list of colors in the donutplot_libAnnotations_colors parameter has a length smaller than 3 and thus will not be used. The default coloring will be used instead. ")
		donutplot_libAnnotations_colors = ['#ff8b00', '#6372b4', "#c6c6c6"]
	if donutplot_mzs_distr_colors is not None and  donutplot_mzs_distr_colors != [None] and len(donutplot_mzs_distr_colors) < 2:
		print(
			"WARNING: the provided list of colors in the donutplot_mzs_distr_colors parameter has a length smaller than 2 and thus will not be used. The default coloring will be used instead. ")
		donutplot_mzs_distr_colors = ["#0072c3", "#42be65"]
	
	# remove blanks if any
	sample_type_rm = []
	if rm_blanks and "BLANKS_TOTAL" in clean_counts_df.columns:
		clean_counts_df = clean_counts_df.loc[clean_counts_df.BLANKS_TOTAL == 0, :]
		sample_type_rm.append("blank")
	# remove beds and controls if requested
	if rm_beds and "BEDS_TOTAL" in clean_counts_df.columns:
		clean_counts_df = clean_counts_df.loc[clean_counts_df.BEDS_TOTAL == 0, :]
		sample_type_rm.append("bed")
	if rm_controls and "CONTROLS_TOTAL" in clean_counts_df.columns:
		clean_counts_df = clean_counts_df.loc[clean_counts_df.CONTROLS_TOTAL == 0, :]
		sample_type_rm.append("control")
	number_valid_mzs = clean_counts_df.shape[0]
	print("- Total number of m/z in the provided clean table =", n)
	print("- Number of valid m/z after filtering =", number_valid_mzs)
	if len(sample_type_rm) > 0:
		print("  - The m/z appearing in the following sample types were removed by the filters: ",
		      ','.join(sample_type_rm))
	if clean_counts_df.shape[0] == 0:
		sys.exit("ERROR: No remaining m/z after the filters. Post processing drug discovery analysis aborted.")
	
	mzs_selected = "" # all mzs
	if use_protonated:
		if "protonated_representative" in clean_counts_df.columns:
			print("- Filtering only the protonated representatives m/z")
			clean_counts_df = clean_counts_df.loc[clean_counts_df.protonated_representative == 1, :]
			number_valid_mzs = clean_counts_df.shape[0]
			print("  - Number of valid m/z after protonated filtering =", number_valid_mzs)
			mzs_selected = "[M+H]+ "  # protonated
			output_name = output_name+"_protonated"
			if clean_counts_df.shape[0] == 0:
				sys.exit("ERROR: No protonated m/z present. Post processing drug discovery analysis aborted.")
		else:
			print("- The protonated_representative column is not present in the provided clean table, not filtering the protonated m/z.")
	
	# get the samples that are present in the provided metadata as the set of samples to use in the plots
	# removing blanks and bed and control if requested
	samples_to_use = [sample_code+"_area"
	                  for sample_code in metadata_df.SAMPLE_CODE.values[~(metadata_df.SAMPLE_TYPE.isin(sample_type_rm))]]
	# get the name of the columns of all samples present in the clean table
	area_cols = [col for col in clean_counts_df.columns
	             if col.endswith('_area') and
	             not col in metadata_df.SAMPLE_CODE.values[(metadata_df.SAMPLE_TYPE.isin(sample_type_rm))]+"_area"]
	
	# create a column to store the m/z with an identification annotation
	clean_counts_df["annotated"] = False
	clean_counts_df.loc[clean_counts_df.curated_identification_best_origin.isin(["GNPS", "UNPD"]), "annotated"] = True
	annotation_lib = clean_counts_df['curated_identification_best_origin'].value_counts()
	annotation_lib["not_annotated"] = (~clean_counts_df.annotated).sum()
	if "GNPS" not in annotation_lib.index:
		annotation_lib["GNPS"] = 0
	if "UNPD" not in annotation_lib.index:
		annotation_lib["UNPD"] = 0
	# check consistency
	if annotation_lib.sum() != number_valid_mzs:
		sys.exit("The number of total m/z (library annotated and not annotated) does not match the number of valid mzs. Something went wrong in the processing.")
	print("- Creating the distribution of library annotation:", annotation_lib.to_dict())
	
	plot_donut_values(Path(output_path,output_name+"_mz_lib_annotation_dist.png"),
	                  values = annotation_lib[sorted(annotation_lib.index.values)].values,
	                  categories=sorted(annotation_lib.index.values),
	                  title = "Distribution of Library Annotations for the "+mzs_selected+"m/z", colors=donutplot_libAnnotations_colors,
	                  title_size=donutplots_title_size, text_size=donutplots_text_size)
	
	# filter only the samples columns, use all area columns here for a correct computation of the indicators
	clean_samples_df = clean_counts_df[area_cols]
	
	# count the number of m/z that appear in each sample (total), the number that are exclusive and/or redudant by sample
	# where the clean_samples_df is the clean count table with only the samples quantification columns
	# check the m/z that are present in each sample, which have a peak area > 0
	present_mzs_by_sample = clean_samples_df.fillna(0) > 0
	# count the number of samples that each m/z appear
	count_samples_by_mz = present_mzs_by_sample.sum(axis=1)
	# check the m/z that appear in only one sample, which are exclusive to a sample
	check_exclusive_mzs = (count_samples_by_mz == 1)
	# check the m/z that appear in more than one sample, which are redundant across the dataset
	check_redundant_mzs = (count_samples_by_mz > 1)
	# count the number of exclusive m/z by sample, which is present in the sample and is an exclusive mz
	exclusive_mzs_by_sample = present_mzs_by_sample.apply(lambda col: col & check_exclusive_mzs)
	
	# create a count table for the m/z with the total of m/z by sample, the number of exclusive m/z by sample and
	# # the number of redundant m/z by sample (which appear in more than one sample) for annotated and not annotated mzs
	mzs_count_by_sample_not_annotated = pd.concat({'total_not_annotated':present_mzs_by_sample[~clean_counts_df.annotated].sum(),
	                                               'exclusive_not_annotated':exclusive_mzs_by_sample[~clean_counts_df.annotated].sum()}, axis=1)
	mzs_count_by_sample_not_annotated['redundant_not_annotated'] = mzs_count_by_sample_not_annotated.total_not_annotated - mzs_count_by_sample_not_annotated.exclusive_not_annotated
	# for annotated mzs
	mzs_count_by_sample_annotated = pd.concat({'total_annotated':present_mzs_by_sample[clean_counts_df.annotated].sum(),
	                                           'exclusive_annotated':exclusive_mzs_by_sample[clean_counts_df.annotated].sum()}, axis=1)
	mzs_count_by_sample_annotated['redundant_annotated'] = mzs_count_by_sample_annotated.total_annotated - mzs_count_by_sample_annotated.exclusive_annotated
	
	print("- Computing the m/z distribution occurrence by sample and creating some plots:\n  "+
	      "- number of exclusive m/z (only appear in one sample);\n  "+
	      "- number of redundant m/z (appear in more than one sample)\n  "+
	      "- number of annotated m/z (received a curated library annotation from GNPS or UNPD).")
	# join the not annotated and the annotated m/z counts by sample
	# df_final_1
	mzs_count_by_sample = pd.concat([mzs_count_by_sample_annotated, mzs_count_by_sample_not_annotated], axis=1)
	mzs_count_by_sample['total'] = mzs_count_by_sample.loc[:,['total_annotated','total_not_annotated']].sum(axis=1)
	if samples_to_use != area_cols:
		mzs_count_by_sample = mzs_count_by_sample.loc[samples_to_use, :]
		print("- Filtering the information of only the selected samples for the output tables and plots. Number of selected samples = ", len(samples_to_use))
	# save table
	mzs_count_by_sample.index = mzs_count_by_sample.index.str.replace('_area', '')
	mzs_count_by_sample.to_csv(Path(output_path, output_name+"_mz_dist_count_by_sample.csv"), index_label='SAMPLE_CODE')
	
	# norm by selected cols, sort by exclusive not annotated - possible novel compounds - and filter topk
	# where the top novelty samples are the samples with more exclusive not annotated m/z
	cols_exclusive_redundant_ann = ['exclusive_not_annotated','exclusive_annotated', 'redundant_not_annotated','redundant_annotated']
	mzs_count_by_sample_norm = mzs_count_by_sample.loc[:,cols_exclusive_redundant_ann].div(mzs_count_by_sample.total, axis=0)*100
	mzs_count_by_sample_norm_topk = mzs_count_by_sample_norm.sort_values("exclusive_not_annotated", ascending=False).head(topk)
	selected_samples_topk = mzs_count_by_sample_norm_topk.index.values
	# plot exclusive and redundant m/z by sample
	plot_stacked_bar_pandas_df(Path(output_path,output_name+"_mz_distribution_top_"+str(topk)+"_samples_novelty.png"),
	                           mzs_count_by_sample_norm_topk,
	                           title='Exclusive and Redundant '+mzs_selected+'m/z Distribution with Library Annotation for the Top '+str(topk)+' Novelty Samples',
	                           xlabel='Samples', ylabel='Percentage of detected '+mzs_selected+'m/z',
	                           figsize=mzs_barplot_figsize, label_size=mzs_barplot_label_size, title_size=mzs_barplot_title_size,
	                           legend_bbox_to_anchor=mzs_barplot_legend_bbox, legend_fontsize=mzs_barplot_legend_fontsize,
	                           legend_ncol=mzs_barplot_legend_ncol,colors=mzs_barplot_colors)
	
	plot_donut_values(Path(output_path,output_name+"_mz_presence_dist_all_donutplot.png"),
	                  values = [check_exclusive_mzs.sum(), check_redundant_mzs.sum()],
	                  categories=["Exclusive", "Redundant"],
	                  title = "Distribution of Exclusive and Redundant "+mzs_selected+"m/z", colors=donutplot_mzs_distr_colors,
	                  title_size=donutplots_title_size, text_size=donutplots_text_size)
	plot_donut_values(Path(output_path,output_name+"_mz_presence_dist_annotated_donutplot.png"),
	                  values = [check_exclusive_mzs[clean_counts_df.annotated].sum(),
	                            check_redundant_mzs[clean_counts_df.annotated].sum()],
	                  categories=["Exclusive", "Redundant"],
	                  title = "Distribution of Exclusive and Redundant "+mzs_selected+"m/z with Library Annotation",
	                  colors=donutplot_mzs_distr_colors,
	                  title_size=donutplots_title_size, text_size=donutplots_text_size)
	
	print("- Computing the superclass grouping distribution of the m/z by sample and creating a plot.")
	# sum the occurrence of each superclass grouping by sample, transpose samples to rows and superclasses to columns
	clean_samples_by_superclass_grouping = present_mzs_by_sample.groupby(clean_counts_df.best_origin_curated_superclass_grouping).sum().T
	if samples_to_use != area_cols:
		clean_samples_by_superclass_grouping = clean_samples_by_superclass_grouping.loc[samples_to_use, :]
	# save table
	clean_samples_by_superclass_grouping.index = clean_samples_by_superclass_grouping.index.str.replace('_area', '')
	clean_samples_by_superclass_grouping.to_csv(Path(output_path, output_name+"_samples_composition_superclass_grouping_count_by_sample.csv"),
	                                            index_label='SAMPLE_CODE')
	
	# filter and order by topk and selected class, removing not annotated
	selected_samples_topk=[sample_code.replace("_area", "") for sample_code in selected_samples_topk]
	clean_samples_by_superclass_grouping_topk = clean_samples_by_superclass_grouping.loc[selected_samples_topk,list(superclass_colors.keys())]
	clean_samples_by_superclass_grouping_topk_norm = clean_samples_by_superclass_grouping_topk.div(clean_samples_by_superclass_grouping_topk.sum(axis=1), axis=0)*100
	plot_stacked_bar_pandas_df(Path(output_path,output_name+"_samples_composition_superclass_grouping_dist_top_"+str(topk)+"_samples_novelty.png"),
	                           clean_samples_by_superclass_grouping_topk_norm,
	                           title='Samples Composition by Superclass Grouping of the Top '+str(topk)+' Novelty Samples (normalized by sample)',
	                           xlabel='Samples', ylabel='Percentage of detected '+mzs_selected+'m/z',
	                           colors=list(superclass_colors.values()),
	                           figsize=superclass_barplot_figsize, label_size=superclass_barplot_label_size,
	                           title_size=superclass_barplot_title_size, legend_title="Superclass Grouping",
	                           legend_bbox_to_anchor=superclass_barplot_legend_bbox,
	                           legend_fontsize=superclass_barplot_legend_fontsize,
	                           legend_ncol=superclass_barplot_legend_ncol)

# convert string to list of text string, parse None values as well
def str2tlist(l):
	return [str(i) if i != "None" else None for i in l.split(',')]
# parse string to list of integers
def str2ilist(l):
	return [int(i) for i in l.split(',')]
# parse string to list of floats
def str2flist(l):
	return [float(i) for i in l.split(',')]
# parse string to boolean
def str2bool(v):
	return v.lower() in ('true', '1', 't')

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description='Post analysis routine for drug discovery research using the NP3 results. Create plots with the m/z occurrence and '+
		            'superclass grouping distributions by sample and store them in the output path together with the respective '+
		            'metrics table. ',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--clean_counts_path", type=str, required=True,
	                    help="Path to the clean counts table file of the NP3 job to be post analysed, the same of the metadata. The prefix of this filename will be used to name the output files.")
	parser.add_argument("--metadata_path", type=str, required=True,
	                    help="Path to the metadata file of the same NP3 job to be post analysed, it may contain the complete metadata or just a selection of the original samples, which will be used to filter the final plots and tables by sample. This filtering does not affect the metrics computations, its only for visualization purposes. The prefix of this filename will be used to name the output files.")
	
	parser.add_argument("--output_path", type=str, required=True,
	                    help="Path to the output directory where the plots and tables of the post analysis will be stored. If the output_path does not exists, it will be created. If it already exists, the created plots and tables will be overwritten.")
	# major filtering parms
	parser.add_argument("--topk", default=None, type=int,
	                    help="The number of samples to be selected and filtered in the plots with the top values of exclusive and not annotated m/z - top novelty samples. If None, disable filtering. This filtering does not affect the metrics computations, its only for visualization purposes.")
	parser.add_argument("--rm_blanks", default=True, type=str2bool,
	                    help="True or False to allow removing blank samples and m/z from the metrics computation.")
	parser.add_argument("--rm_beds", default=True, type=str2bool,
	                    help="True or False to allow removing culture media samples and m/z from the metrics computation.")
	parser.add_argument("--rm_controls", default=False, type=str2bool,
	                    help="True or False to allow removing control samples and m/z from the metrics computation.")
	parser.add_argument("--use_protonated", default=False, type=str2bool,
	                    help="True of False defining if only the putative [M+H] m/z should be used in the output tables and plots (filter the table with protonated_representative == 1). This will affect the metrics computation.")
	# plots parms
	parser.add_argument("--donutplots_title_size", default=16, type=int,
	                    help="The title size of the donut plots.")
	parser.add_argument("--donutplots_text_size", default=14, type=int,
	                    help="The axis and legend text sizes of the donut plots.")
	parser.add_argument("--donutplot_libAnnotations_colors", default='#ff8b00,#6372b4,#c6c6c6', type=str2tlist,
	                    help="The list of colors separated by comma for the library annotation distribution donut plot. Three colors are expected for 'GNPS', 'UNPD' and 'not_annotated' categories. Or None to use the default coloring of matplotlib.")
	parser.add_argument("--donutplot_mzs_distr_colors", default='#0072c3,#42be65', type=str2tlist,
	                    help="The list of colors separated by comma for the m/z occurrence distribution donut plot. Two colors are expected for 'Exclusive' and 'Redundant' categories. Or None to use the default coloring of matplotlib.")
	parser.add_argument("--mzs_barplot_figsize", default='18,8', type=str2ilist,
	                    help="The x,y figure size of the m/z occurrence distribution in a stacked bar plot by sample. Two integer values separated by comma.")
	parser.add_argument("--mzs_barplot_label_size", default=13, type=int,
	                    help="The axis label size of the m/z occurrence distribution in a stacked bar plot by sample.")
	parser.add_argument("--mzs_barplot_title_size", default=20, type=int,
	                    help="The title size of the m/z occurrence distribution in a stacked bar plot by sample.")
	parser.add_argument("--mzs_barplot_legend_bbox", default='0.5,-0.2', type=str2flist,
	                    help="The x,y anchoring coordinates relative to the plot area of the m/z occurrence distribution in a stacked bar plot by sample. Two float values separated by comma. (0, 0) is the bottom-left corner of the plot. (1, 1) is the top-right corner of the plot. Values smaller than 0 or greater than 1 will place the legend completely outside the plot area.")
	parser.add_argument("--mzs_barplot_legend_fontsize", default=17, type=int,
	                    help="The legend font size of the m/z occurrence distribution in a stacked bar plot by sample.")
	parser.add_argument("--mzs_barplot_legend_ncol", default=4, type=int,
	                    help="The number of columns to display the legends of the m/z occurrence distribution in a stacked bar plot by sample.")
	parser.add_argument("--mzs_barplot_colors", default='#0072c3,#FF8C00,#42be65,#FDDA0D', type=str2tlist,
	                    help="The list of colors of the m/z occurrence distribution in a stacked bar plot by sample. Four colors are expected for 'exclusive not annotated', 'exclusive annotated', 'redundant not annotated' and 'redundant annotated' m/z categories. If None, use the default matplotlib coloring.")
	parser.add_argument("--superclass_barplot_figsize", default='20,10', type=str2ilist,
	                    help="The x,y figure size of the superclass distribution in a stacked bar plot by sample. Two integer values separated by comma.")
	parser.add_argument("--superclass_barplot_label_size", default=16, type=int,
	                    help="The axis label size of the superclass distribution in a stacked bar plot by sample.")
	parser.add_argument("--superclass_barplot_title_size", default=20, type=int,
	                    help="The title size of the superclass distribution in a stacked bar plot by sample.")
	parser.add_argument("--superclass_barplot_legend_bbox", default='0.5,-0.15', type=str2flist,
	                    help="The x,y anchoring coordinates relative to the plot area of the superclass distribution in a stacked bar plot by sample. Two float values separated by comma. (0, 0) is the bottom-left corner of the plot. (1, 1) is the top-right corner of the plot. Values smaller than 0 or greater than 1 will place the legend completely outside the plot area.")
	parser.add_argument("--superclass_barplot_legend_fontsize", default=15, type=int,
	                    help="The legend font size and title size of the superclass distribution in a stacked bar plot by sample.")
	parser.add_argument("--superclass_barplot_legend_ncol", default=4, type=int,
	                    help="The number of columns to display the legends of the superclass distribution in a stacked bar plot by sample.")
		
	args = parser.parse_args()
	
	post_dd_analysis_plots(metadata_path=args.metadata_path, clean_counts_path=args.clean_counts_path,
	                       output_path=args.output_path, topk=args.topk,
	                       rm_blanks=args.rm_blanks,
	                       rm_beds=args.rm_beds, rm_controls=args.rm_controls,
	                       use_protonated=args.use_protonated,
	                       donutplots_title_size=args.donutplots_title_size,
	                       donutplots_text_size=args.donutplots_text_size,
                           donutplot_libAnnotations_colors=args.donutplot_libAnnotations_colors,
                           donutplot_mzs_distr_colors=args.donutplot_mzs_distr_colors,
                           mzs_barplot_figsize=args.mzs_barplot_figsize,
                           mzs_barplot_label_size=args.mzs_barplot_label_size,
                           mzs_barplot_title_size=args.mzs_barplot_title_size,
                           mzs_barplot_legend_bbox=args.mzs_barplot_legend_bbox,
                           mzs_barplot_legend_fontsize=args.mzs_barplot_legend_fontsize,
                           mzs_barplot_legend_ncol=args.mzs_barplot_legend_ncol,
                           mzs_barplot_colors=args.mzs_barplot_colors,
                           superclass_barplot_figsize=args.superclass_barplot_figsize,
                           superclass_barplot_label_size=args.superclass_barplot_label_size,
                           superclass_barplot_title_size=args.superclass_barplot_title_size,
                           superclass_barplot_legend_bbox=args.superclass_barplot_legend_bbox,
                           superclass_barplot_legend_fontsize=args.superclass_barplot_legend_fontsize,
                           superclass_barplot_legend_ncol=args.superclass_barplot_legend_ncol)
	