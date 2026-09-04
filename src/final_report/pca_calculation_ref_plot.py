# diff PCA from R and from python: https://stackoverflow.com/questions/77275858/discrepancy-in-pca-results-between-pythons-sklearn-and-rs-prcomp
# code PCA from https://stackoverflow.com/questions/39216897/plot-pca-loadings-and-loading-in-biplot-in-sklearn-like-rs-autoplot
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.colors as mcolors
import pandas as pd
from adjustText import adjust_text
import sys
from pathlib import Path

# modified from https://www.google.com/search?q=python+plot+pca+principal+components+unit+circle&rlz=1C1FCXM_pt-PTBR998BR998&oq=python+plot+pca+principal+components+unit+circle&gs_lcrp=EgZjaHJvbWUyBggAEEUYOTIHCAEQIRigAdIBCTE0MTU2ajBqN6gCALACAA&sourceid=chrome&ie=UTF-8
# Plotting the correlation circle
# visualize the quality of representation (cos2) of rows/columns from the results of Principal Component Analysis (PCA)
def plot_arrows_correlation_circle(output_path,scaled_loading, feature_names, pvars, color_cos2_cmap):
    cos2_loadings = (scaled_loading ** 2).sum(1)
    plt.figure(figsize=(10, 8))
    plt.set_cmap(color_cos2_cmap)
    sm = plt.cm.ScalarMappable(cmap=color_cos2_cmap,
                               norm=mcolors.Normalize(0, 1))  # Create a dummy mappable object for the colorbar
    sm.set_array(cos2_loadings)  # Required for ScalarMappable to work without actual plot data
    plt.rc('axes', edgecolor='#CCCCCC',linewidth=0.8) # change the axes color
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.grid(False)
    plt.xlabel(f"PC1 ({pvars[0]:.2f} %)")
    plt.ylabel(f"PC2 ({pvars[1]:.2f} %)")
    plt.title('PCA Correlation Circle - quality of representation (cos2)')
    circle = plt.Circle((0, 0), 1, color='grey', fill=False, linestyle='-', alpha=0.8) # Draw the unit circle
    plt.gca().add_patch(circle)
    texts = []
    for i, feature in enumerate(feature_names): # Plot the loading vectors for each original variable
        plt.arrow(0, 0, scaled_loading[i, 0], scaled_loading[i, 1],
                  head_width=0.05, head_length=0.05, color=sm.cmap(cos2_loadings[i]), length_includes_head=True)
        texts.append(plt.text(scaled_loading[i, 0], scaled_loading[i, 1], feature, color=sm.cmap(cos2_loadings[i])))
    adjust_text(texts, expand=(1.5, 3), arrowprops=dict(arrowstyle="-", color='lightgray', lw=0.5))
    plt.axhline(0, color='black', linewidth=1.5, linestyle='--', alpha=0.9)
    plt.axvline(0, color='black', linewidth=1.5, linestyle='--', alpha=0.9)
    plt.grid(which='major', color='#CCCCCC', alpha=0.7, linewidth=0.8, linestyle=':')
    plt.colorbar(sm, label="cos2", shrink=0.3, aspect=7, ticks=np.arange(0.5, 1.05, 0.1),
                 boundaries=np.linspace(0.5, 1, 256))
    plt.savefig(output_path / "pca_quality_representation_cos2_NP3_reference.png", dpi=300, bbox_inches='tight')
    plt.close()

# n equals the number of components, components equals the x, y values of each n arrow
# and feature_labels the names of the components
def plot_components_arrows(components, feature_labels):
    n = components.shape[0]
    width = -0.0025 * np.min([np.subtract(*plt.xlim()), np.subtract(*plt.ylim())])
    texts = []  # stores the components labels - to avoid overlapping
    for i in range(n):
        plt.arrow(0, 0, components[i, 0], components[i, 1], width=width, color='k', alpha=0.6, ec='none',
                  length_includes_head=True)
        if feature_labels is None:
            texts.append(plt.text(components[i, 0], components[i, 1], "Var" + str(i + 1), color='0.1', ha='center',
                                  va='center', size='x-large', fontweight='bold'))
        else:
            texts.append(plt.text(components[i, 0], components[i, 1], feature_labels[i], color='0.1', ha='center',
                                  va='center', size='x-large',
                                  fontweight='bold'))  # adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    adjust_text(texts, expand=(2.5, 1.2), arrowprops=dict(arrowstyle="-", color='0.1', lw=0.5))

# function to plot a PCA result with 2-cols stored in scores, one point by row, with the point labels listed in
# the point_labels list (must contain the UNPD label, which goes in the background) and formatted using the data_types_style dictionary containing the unique labels as keys and a 2 sized list
# in the values containing: (0) the name of the label to be used in the plot legend and (1) the color to be used in the points with
# the respective label; If components is informed as a 2-cols, the arrows of the informed components are plotted and labeled using the feature_labels
# if present with the descriptors names, instead they receive the default naming Var1, Var2, ...; The pvars may contain the explained variance of each
# axis in percentage (%) to be used in the axis labeling
def biplot_scatter_arrows(scores,point_labels, data_types_style,components=None,feature_labels=None,pvars=None):
    point_labels = np.asarray(point_labels)
    types_label = ['UNPD'] + list(np.unique(point_labels)[np.unique(point_labels) != 'UNPD']) # get unique point labels, put UNPD first to be in the background
    for name in types_label:
        if np.any(point_labels == name): # if no point of this label, skip it
            plt.scatter(*zip(*scores[point_labels == name]), label=data_types_style[name][0],
                        alpha = 0.8, c=data_types_style[name][1])
    plt.legend(title='Data Type', loc = 'upper right', fontsize='x-large')
    if components is not None:
        plot_components_arrows(components, feature_labels)
    if pvars is None:
        plt.xlabel("PC{}".format(1))
        plt.ylabel("PC{}".format(2))
    else:
        plt.xlabel(f"PC1 ({pvars[0]:.2f} %)")
        plt.ylabel(f"PC2 ({pvars[1]:.2f} %)")
    plt.grid(which='major', color='#CCCCCC', alpha=0.7, linewidth=0.8, linestyle=':')
    plt.axhline(0, color='black', linestyle='--', alpha=0.9)
    plt.axvline(0, color='black', linestyle='--', alpha=0.9)
    if np.any(point_labels == 'UNPD') or np.any(point_labels == 'GNPS'):
        plt.xlim(-30, 10) # with scale plt.xlim(-0.8, 0.25)
        plt.ylim(-7, 33) # with scale plt.ylim(-0.2, 0.9)
    plt.tight_layout()
    #plt.savefig(output_img_path)
    #plt.show()
    #plt.close()

# colors from https://matplotlib.org/stable/gallery/color/color_cycle_default.html#sphx-glr-gallery-color-color-cycle-default-py
data_types_style_NP3_reference = {'UNPD': ['UNPD 2018', 'tab:grey'],
                        'DrugBank': ['DrugBank 2019', 'tab:blue'],
                        'allosteric_review': ['Allosteric natural - PubMed 2019', 'tab:green'],
                        'mzs' : ['m/zs', 'tab:blue']}

def biplot_scatter_reference_PCA_wComponents(scores, point_labels, pvars=None, components=None, feature_labels=None):
    plt.subplots(1)
    # uses the refence data types style
    data_types_style = data_types_style_NP3_reference
    # plot PCA
    biplot_scatter_arrows(scores, point_labels, data_types_style, components, feature_labels, pvars)

# function to plot and save the reference PCA of NP3
def biplot_save_reference_PCA_wComponents(output_path, output_name, scores, point_labels, pvars=None, components=None,
                                          feature_labels=None):
    # plot PCA points
    biplot_scatter_reference_PCA_wComponents(scores, point_labels, pvars, components, feature_labels)
    # save the plot, show and close connection
    if components is not None:
        plt.savefig(output_path / (output_name+"chemical_space_NP3_reference_biplot_components.png"), dpi=300, bbox_inches='tight')
    else:
        plt.savefig(output_path / (output_name+"chemical_space_NP3_reference_biplot.png"), dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()

data_types_style_new = {'UNPD': ['NP3-UNPD', 'tab:pink'],
                        'GNPS': ['NP3-GNPS', 'tab:orange']}
# list of other colors to use for external table input, max to 8 classes
data_types_style_new_external = ['tab:pink','tab:orange','tab:purple','tab:brown','tab:red','tab:olive','tab:cyan','#ff9896']

# save PCA biplot of a new data - personalize name using data_type and output_name
def biplot_save_new_PCA_wComponents(output_path, output_name, data_type, scores, point_labels, pvars=None,
                                    components=None, feature_labels=None):
    #  create data_types_style for different data_types from external table
    if data_type not in ['UNPD', 'GNPS', 'UNPDxGNPS']:
        data_types_style = {label: [label,data_types_style_new_external[i]]
                            for i,label in enumerate(np.unique(point_labels))}
    else:
        data_types_style = data_types_style_new
    # plot PCA
    biplot_scatter_arrows(scores, point_labels, data_types_style, components=None, feature_labels=None, pvars=pvars)
    if components is not None:
        plot_components_arrows(components, feature_labels)
        plt.savefig(output_path / (output_name+"_chemical_space_NP3_"+data_type+"_PCA_biplot_components.png"), dpi=300, bbox_inches='tight')
    else:
        plt.savefig(output_path / (output_name+"_chemical_space_NP3_"+data_type+"_PCA_scores.png"), dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()


#set_descriptors_rcdk = ['nSmallRings','nAromRings','nRingBlocks','nAromBlocks','nRings3','nRings4','nRings5','nRings6','nRings7','nRings8','nRings9','tpsaEfficiency','Zagreb','WPATH','WPOL','WTPT.1','WTPT.2','WTPT.3','WTPT.4','WTPT.5','VAdjMat','VABC','TopoPSA','topoShape','geomShape','PetitjeanNumber','MDEC.11','MDEC.12','MDEC.13','MDEC.14','MDEC.22','MDEC.23','MDEC.24','MDEC.33','MDEC.34','MDEC.44','MDEO.11','MDEO.12','MDEO.22','MDEN.11','MDEN.12','MDEN.13','MDEN.22','MDEN.23','MDEN.33','khs.sLi','khs.ssBe','khs.ssssBe','khs.ssBH','khs.sssB','khs.ssssB','khs.sCH3','khs.dCH2','khs.ssCH2','khs.tCH','khs.dsCH','khs.aaCH','khs.sssCH','khs.ddC','khs.tsC','khs.dssC','khs.aasC','khs.aaaC','khs.ssssC','khs.sNH3','khs.sNH2','khs.ssNH2','khs.dNH','khs.ssNH','khs.aaNH','khs.tN','khs.sssNH','khs.dsN','khs.aaN','khs.sssN','khs.ddsN','khs.aasN','khs.ssssN','khs.sOH','khs.dO','khs.ssO','khs.aaO','khs.sF','khs.sSiH3','khs.ssSiH2','khs.sssSiH','khs.ssssSi','khs.sPH2','khs.ssPH','khs.sssP','khs.dsssP','khs.sssssP','khs.sSH','khs.dS','khs.ssS','khs.aaS','khs.dssS','khs.ddssS','khs.sCl','khs.sGeH3','khs.ssGeH2','khs.sssGeH','khs.ssssGe','khs.sAsH2','khs.ssAsH','khs.sssAs','khs.sssdAs','khs.sssssAs','khs.sSeH','khs.dSe','khs.ssSe','khs.aaSe','khs.dssSe','khs.ddssSe','khs.sBr','khs.sSnH3','khs.ssSnH2','khs.sssSnH','khs.ssssSn','khs.sI','khs.sPbH3','khs.ssPbH2','khs.sssPbH','khs.ssssPb','Kier1','Kier2','Kier3','HybRatio','fragC','FMF','ECCEN','SP.0','SP.1','SP.2','SP.3','SP.4','SP.5','SP.6','SP.7','VP.0','VP.1','VP.2','VP.3','VP.4','VP.5','VP.6','VP.7','SPC.4','SPC.5','SPC.6','VPC.4','VPC.5','VPC.6','SC.3','SC.4','SC.5','SC.6','VC.3','VC.4','VC.5','VC.6','SCH.3','SCH.4','SCH.5','SCH.6','SCH.7','VCH.3','VCH.4','VCH.5','VCH.6','VCH.7','C1SP1','C2SP1','C1SP2','C2SP2','C3SP2','C1SP3','C2SP3','C3SP3','C4SP3','ATSp1','ATSp2','ATSp3','ATSp4','ATSp5','ATSm1','ATSm2','ATSm3','ATSm4','ATSm5','ATSc1','ATSc2','ATSc3','ATSc4','ATSc5','topoShape1','geomShape1','MOMI.X','MOMI.Y','MOMI.Z','MOMI.XY','MOMI.XZ','MOMI.YZ','MOMI.R','LOBMAX','LOBMIN','GRAV.1','GRAV.2','GRAV.3','GRAVH.1','GRAVH.2','GRAVH.3','GRAV.4','GRAV.5','GRAV.6','PPSA.1','PPSA.2','PPSA.3','PNSA.1','PNSA.2','PNSA.3','DPSA.1','DPSA.2','DPSA.3','FPSA.1','FPSA.2','FPSA.3','FNSA.1','FNSA.2','FNSA.3','WPSA.1','WPSA.2','WPSA.3','WNSA.1','WNSA.2','WNSA.3','RPCG','RNCG','RPCS','RNCS','THSA','TPSA','RHSA','RPSA','Fsp3','XLogP','MW','LipinskiFailures','nRotB','MLogP','nAtomP','nAtomLC','nB','nBase','nAtom','nAromBond','naAromAtom','ALogP','ALogp2','AMR','nAcid']
set_reference_descriptors_bestCos_top24 = ['WPOL','ATSp1','ATSp2','SP.3','Zagreb','VP.3','naAromAtom','SP.2','ATSp3','VP.2','ATSm3','nB','ATSm2','nAromRings','nAromBlocks','C2SP2','khs.aasC','nAtomP','khs.aaCH','LipinskiFailures','MW','nAtom','VP.1','WTPT.1']


# Compute PCA from reference data and apply to transform new data
# This process effectively reuses the eigenvectors learned from the original data to transform subsequent data,
# ensuring consistency in the dimensionality reduction.
# the reference data is expected to have unique mostly smiles (duplicates are not removed) and the top 24 descriptors columns,
# also the Type, EntryID and SubType columns
# data_reference_path must contain the list of top 24 descriptors; new_data_path is the clean data for UNPD and the descriptors list for GNPS or an external dataset
# If data_type is UNPD: output_path is the folder to store the plots, inside the final_reports/chemical_report/chemical_space_identifications
# else if data_type is GNPS (also creates for UNPDxGNPS): output_path is the final clustering result, inside the outs dir
# otherwise if data_type is an external data, plot directly to the output_path
# for data_type == UNPD or GNPS, no blanks or beds are plotted
def pca_calculation_smiles_rcdk_ref_plot(data_reference_path, new_data_path, output_path, output_name, data_type="UNPD"):
    # check if files exists
    data_reference_path = Path(data_reference_path)
    new_data_path = Path(new_data_path)
    output_path = Path(output_path)
    if not data_reference_path.exists() or not data_reference_path.is_file():
        sys.exit("The provided path to the PCA data reference file does not exists: "+ data_reference_path.as_posix()+
                 ". PCA plotting aborted.")
    if not new_data_path.exists() or not new_data_path.is_file():
        sys.exit("The provided path to the new data file, to be used to calculate its reference PCA, does not exists: "+
                 new_data_path.as_posix()+
                 ". PCA plotting aborted.")
    if not output_path.exists() or not output_path.is_dir():
        sys.exit("The provided output path does not exists: "+output_path.as_posix() +". PCA plotting aborted.")
    
    # Original data
    # X_original is read here
    data_reference_descriptors = pd.read_csv(data_reference_path,
                                             low_memory=False)
    print("* Creating the NP3 reference PCA from SMILES for plotting",data_type,"*")
    # not removing duplicates here - the reference contain the minimum of duplicates
    #data_reference_descriptors = data_reference_descriptors.loc[~data_reference_descriptors.SMILES.duplicated(),:]
    # filter the top 24 descriptors selected
    X_original = data_reference_descriptors.loc[:,set_reference_descriptors_bestCos_top24]
    
    #Get valid Calculations cols
    valid_descriptors_cols = []
    for col in X_original.columns:
        try:
            if (X_original[col].var() > 0):
                valid_descriptors_cols.append(col)
        except: # rm NAs descriptors
            pass
    
    # remove zero variance cols
    X_filtered = X_original.loc[:, valid_descriptors_cols]
    # remove rows with NaNs
    valid_smiles_rows = (np.isnan(X_filtered).sum(1) == 0)
    X_filtered = X_filtered.loc[valid_smiles_rows, :]
    
    # filter label in Type
    point_labels = data_reference_descriptors.loc[valid_smiles_rows,"Type"]
    
    # extract final features
    features = X_filtered.columns.values
    
    # standardization
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_filtered)
    
    # PCA
    pca = PCA(n_components=2).fit(X_scaled)
    X_reduced = pca.transform(X_scaled)
    
    # Now pca.components_ holds the calculated eigenvectors
    # and pca.explained_variance_ holds the eigenvalues
    # proportions of variance explained by axes
    pvars = pca.explained_variance_ratio_ * 100
    print(f"  - Explained variance: PC1 {pvars[0]:.2f}% and PC2 {pvars[1]:.2f}%")
    print(f"  - Cumulative: {np.cumsum(pvars)[1]:.2f}%")
    
    # coordinates of samples (i.e., scores; let's take the first two axes)
    scores = X_reduced
    # mirror the first axis by convention
    scores[:,0] = scores[:,0]*-1
    # these could be scaled with the following code
    #scaled_scores = scores
    #scalex = 1.0/(scaled_scores[:,0].max() - scaled_scores[:,0].min())
    #scaley = 1.0/(scaled_scores[:,1].max() - scaled_scores[:,1].min())
    #scaled_scores[:,0] = scaled_scores[:,0] * scalex
    #scaled_scores[:,1] = scaled_scores[:,1] * scaley
    
    # coordinates of features - descriptors/components (i.e., loadings; note the transpose)
    loadings = pca.components_.T
    loadings[:,0] = loadings[:,0]*-1 # mirror data in the y-axis
    
    # Scale the features (arrows) properly to match the samples (points).
    # The following code scales by the maximum absolute value of samples on each axis.
    arrows = loadings * np.abs(scores).max(axis=0)
    
    ## could plot only top k features
    #k=20
    # Method 1: Find top k arrows that appear the longest (i.e., furthest from the origin) in the visible plot:
    #tops = (loadings ** 2).sum(axis=1).argsort()[-k:]
    #arrows = loadings[tops]
    # Method 2: Find top k features that drive most variance in the visible PCs:
    #tops = (loadings * pvars).sum(axis=1).argsort()[-k:]
    #arrows = loadings[tops]
    # scale arrows for features:  For all features, the sum of square loadings is always 1 per PC. With a small portion of features, we should bring them up such that the sum of square loadings of them is also 1.
    #arrows /= np.sqrt((arrows ** 2).sum(axis=0))
    
    # create default parameters for plotting
    params = {'legend.fontsize': 'x-large',
                'legend.title_fontsize': 'x-large',
                'figure.figsize': (9, 9),
                'axes.labelsize': 'x-large', 'axes.labelweight': 'bold', 'axes.titlesize': 'x-large',
                'xtick.labelsize': 'x-large', 'ytick.labelsize': 'x-large',
                'font.weight': 'bold'}
    pylab.rcParams.update(params) # set styling
    
    # use the following code to plot the default PCA with components and without components
    # biplot_save_reference_PCA_wComponents(output_path, scores, point_labels, pvars arrows, features)
    # biplot_save_reference_PCA_wComponents(output_path, scores, point_labels, pvars, None, None)
    
    print("  - Transforming the new", data_type,"identified data to the reference PCA chemical space and plotting")
    # And then scale, fit and transform the new data from NP3 to the
    # reference PCA and plot the result on top of it
    ###
    ##
    # New data identified that you want to transform using the same eigenvectors
    # read the new data here and filter the correct IDs
    new_data = pd.read_csv(new_data_path, low_memory=False)
    # check if the protonated_representative columns exists, if not create it with fake 0's
    if 'protonated_representative' not in new_data.columns:
        new_data.loc[:,'protonated_representative'] = 0
    # for UNPD extraxt the X_new by matching the best identified SMILES against the data_reference_descriptors SMILES,
    # # and return the corresponding not NA descriptors
    # for GNPS a table with the descriptors must be informed
    if data_type == "UNPD":
        # the new_data here is the clean table with identification curated for UNPD
        # the best identified UNPD ID is in column tremolo_UNPD_IDs_best
        # the new_data  must contain the tremolo_UNPD_IDs_best column with the UNPD ID of the best identification
        # and must contain a valid value in it (not NA)
        # the list of UNPD IDs will be matched against the reference table EntryID and its information will be filtered out
        # test this merge, merge to retrieve the descriptors instead and set to the x_filtered
        # check if curated result is present, if not skipt PCA
        # add the curated tag to the output name
        output_name = output_name + "_curated"
        if "tremolo_UNPD_score_best" not in new_data.columns:
            print("Invalid data from UNPD informed! The provided clean table does not have the tremolo-UNPD "+
                     "curated library annotation result, column 'tremolo_UNPD_score_best' is missing. Skipping PCA for library annotations.")
            return None
        if "tremolo_SMILES_best" not in new_data.columns:
            print("Invalid data from UNPD informed! The provided clean table does not have the tremolo-UNPD "+
                     "best library annotation result, column 'tremolo_SMILES_best' is missing. Skipping PCA for library annotations.")
            return None
        # filter only the not blank results
        if "BLANKS_TOTAL" in new_data.columns:
            new_data = new_data.loc[new_data.BLANKS_TOTAL == 0, :]
        if new_data.shape[0] == 0 or not (~new_data.tremolo_SMILES_best.isna()).any():
            print("No not blank and valid tremolo-UNPD curated library annotation is present in the provided clean table. "+
                  "Could not create the chemical space. Skipping PCA for library annotations.")
            return None
        # filter only the not bed results
        if "BEDS_TOTAL" in new_data.columns:
            new_data = new_data.loc[new_data.BEDS_TOTAL == 0, :]
        if new_data.shape[0] == 0 or not (~new_data.tremolo_SMILES_best.isna()).any():
            print("No not blank and not bed and valid tremolo-UNPD curated library annotation is present in the provided clean table. "+
                  "Could not create the chemical space. Skipping PCA for library annotations.")
        # filter only the curated result with tremolo_UNPD_score_best > 0
        new_data = new_data.loc[new_data.tremolo_UNPD_score_best > 0,:]
        if new_data.shape[0] == 0 or not (~new_data.tremolo_SMILES_best.isna()).any():
            print("No valid tremolo-UNPD curated library annotation is present in the provided clean table. "+
                  "Could not create the chemical space. Skipping PCA for library annotations.")
            return None
        # merge the new data with the reference descriptors using the best identified SMILES
        new_data_desc = pd.merge(new_data.loc[~new_data.tremolo_SMILES_best.isna(), ["tremolo_SMILES_best","protonated_representative"]],
                                 data_reference_descriptors.loc[:, ["SMILES"] + set_reference_descriptors_bestCos_top24],
                                 left_on="tremolo_SMILES_best", right_on="SMILES", how="left")
        for protonated in [0,1]:
            if protonated == 1:
                output_name = output_name + "_protonated"
                # filter only the protonated nodes
                new_data_desc = new_data_desc.loc[new_data_desc.protonated_representative == 1,:]
            # remove not matched entries
            new_data_desc_valid = new_data_desc.loc[~new_data_desc.SMILES.isna(), :]
            # check if new_data_desc_valid have valid descriptors and transform the data to the reference PCA
            if new_data_desc_valid.shape[0] > 0:
                # filter the descriptors columns and scale
                X_new = new_data_desc_valid.iloc[:, 3:]
                X_new_scaled = scaler.fit_transform(X_new)
                
                # Transform the new data using the previously fitted PCA object
                X_transformed = pca.transform(X_new_scaled)
                
                # now plot reference PCA without components as background for the new data
                biplot_scatter_reference_PCA_wComponents(scores, point_labels, pvars=pvars, components=None,
                                                         feature_labels=None)
                # plot the result new data in the PCA reference and save the result without components
                biplot_save_new_PCA_wComponents(output_path, output_name, data_type, scores=X_transformed,
                                                point_labels=[data_type] * new_data_desc_valid.shape[0],
                                                pvars=pvars)
                # then save it with components
                biplot_scatter_reference_PCA_wComponents(scores, point_labels, pvars=pvars, components=None,
                                                         feature_labels=None)
                biplot_save_new_PCA_wComponents(output_path, output_name, data_type, scores=X_transformed,
                                                point_labels=[data_type] * new_data_desc_valid.shape[0],
                                                pvars=pvars, components=arrows, feature_labels=features)
            else:
                print(f"No valid tremolo-UNPD curated {'and protonated ' if protonated == 1 else ''}library annotation after merging with reference SMILES. Skipping PCA for library annotations.")
                return None
    elif data_type == "GNPS":
        # add the curated tag to the output name
        output_name = output_name + "_curated"
        # make pca with GNPS best separated and together with the unpd best - UNPDxGNPS.
        # make pca with GNPS best result, read the calculated descriptors from output_path/identifications/gnps_results_smiles_descriptorsCDK.csv
        # here output_path is the final clustering result, inside the outs dir
        # the output path to store the results is output_path/final_reports/chemical_report/chemical_space_identifications
        # filter only the descriptors
        gnps_descriptors_file = output_path / "identifications" / "gnps_results_smiles_descriptorsCDK.csv"
        if not gnps_descriptors_file.exists() or not gnps_descriptors_file.is_file():
            print("The GNPS descriptors file is not present in: "+ gnps_descriptors_file.as_posix()+
                     ". Provide the correct output_path or check for errors in the RCDK descriptors calculation. "+
                     "PCA plotting for GNPS and UNPDxGNPS aborted.")
            return None
        output_path = output_path / "final_reports" / "chemical_report" / "chemical_space_identifications"
        output_path.mkdir(exist_ok=True)
        if not output_path.exists() or not output_path.is_dir():
            print("The final output directory of the chemical report to store the chemical space PCA result could not"+
                  " be created. PCA plotting for GNPS and UNPDxGNPS aborted.")
            return None
        # filter only the not blank results
        if "BLANKS_TOTAL" in new_data.columns:
            new_data = new_data.loc[new_data.BLANKS_TOTAL == 0, :]
        if new_data.shape[0] == 0 or not (~new_data.gnps_Smiles.isna()).any():
            print(
                "No not blank and valid GNPS curated library annotation is present in the provided clean table. " +
                "Could not create the chemical space. Skipping PCA for library annotations.")
            return None
        # filter only the not bed results
        if "BEDS_TOTAL" in new_data.columns:
            new_data = new_data.loc[new_data.BEDS_TOTAL == 0, :]
        if new_data.shape[0] == 0 or not (~new_data.gnps_Smiles.isna()).any():
            print(
                "No not blank and not bed and valid GNPS curated library annotation is present in the provided clean table. " +
                "Could not create the chemical space. Skipping PCA for library annotations.")
        # filter the valid gnps smiles with gnps_score > 0 and get their descriptors
        new_data_desc = new_data.loc[new_data.gnps_score > 0 ,["gnps_Smiles","protonated_representative"]]
        if new_data_desc.shape[0] == 0:
            print("  - No valid GNPS curated library annotation with score > 0. PCA plotting for GNPS and UNPDxGNPS aborted.")
            return None
        # read the descriptors calculation result and filter the valid ones
        gnps_descriptors_result = pd.read_csv(gnps_descriptors_file, low_memory=False)
        for protonated in [0,1]:
            if protonated == 1:
                output_name_plot = output_name + "_protonated"
                # filter only the protonated nodes
                new_data_desc = new_data_desc.loc[new_data_desc.protonated_representative == 1,:]
                if new_data_desc.shape[0] == 0:
                    print(
                        "  - No valid GNPS curated and protonated library annotation with score > 0. PCA plotting for GNPS and UNPDxGNPS aborted.")
                    break
            else:
                output_name_plot = output_name
            new_data_desc_valid = gnps_descriptors_result.loc[gnps_descriptors_result.gnps_Smiles.isin(new_data_desc.gnps_Smiles),:]
            # filter only the descriptors columns
            X_new = new_data_desc_valid.iloc[:, 4:]
            # remove rows with NaNs
            valid_smiles_rows = (np.isnan(X_new).sum(1) == 0)
            X_new_filtered = X_new.loc[valid_smiles_rows, :]
            if X_new_filtered.shape[0] == 0:
                if protonated == 0:
                    print(
                        "  - No valid GNPS curated library annotation with score > 0. PCA plotting for GNPS and UNPDxGNPS aborted.")
                    return None
                else:
                    print(
                        "  - No valid GNPS curated and protonated library annotation with score > 0. PCA plotting for GNPS protonated aborted.")
                    break # if no valid only for protonated, then go for UNPDxGNPS
            X_new_scaled = scaler.fit_transform(X_new_filtered)
            X_transformed = pca.transform(X_new_scaled)
            # now plot reference PCA without components as background for the new data
            biplot_scatter_reference_PCA_wComponents(scores, point_labels, pvars=pvars, components=None,
                                                     feature_labels=None)
            # plot the result new data in the PCA reference and save the result without components
            biplot_save_new_PCA_wComponents(output_path, output_name_plot, data_type, scores=X_transformed,
                                            point_labels=[data_type] * X_new_filtered.shape[0],
                                            pvars=pvars)
            # then save it with components
            biplot_scatter_reference_PCA_wComponents(scores, point_labels, pvars=pvars, components=None,
                                                     feature_labels=None)
            biplot_save_new_PCA_wComponents(output_path, output_name_plot, data_type, scores=X_transformed,
                                            point_labels=[data_type] * X_new_filtered.shape[0],
                                            pvars=pvars, components=arrows, feature_labels=features)
            
        ## Now get best origin identification data_type == "UNPDxGNPS", if both results are present
        # and plot PCA of UNPD and GNPS together using the curated_lib_annotation_origin
        data_type = "UNPDxGNPS"
        print("  - Transforming the new", data_type, "identified data to the reference PCA chemical space and plotting")
        if new_data.columns.isin(["curated_lib_annotation_origin", "gnps_Smiles", "tremolo_SMILES_best"]).sum() < 3:
            print("  - The tremolo-UNPD best result is not present, skipping UNPDxGNPS PCA plotting.")
            # skipt to the components circle plot
        else:
            new_data_desc = new_data.loc[~new_data.curated_lib_annotation_origin.isna(),
                                         ["curated_lib_annotation_origin", "gnps_Smiles", "tremolo_SMILES_best",
                                          "protonated_representative"]]
            if new_data_desc.curated_lib_annotation_origin.unique().size < 2:
                print("  - The curated_lib_annotation_origin only contains results from one source ('",
                      new_data_desc.curated_lib_annotation_origin.unique()[0],"'), skipping UNPDxGNPS PCA plotting.")
            else:
                for protonated in [0, 1]:
                    if protonated == 1:
                        output_name_plot = output_name + "_protonated"
                        # filter only the protonated nodes
                        new_data_desc = new_data_desc.loc[new_data_desc.protonated_representative == 1, :]
                        if new_data_desc.shape[0] == 0:
                            print(
                                "  - No valid UNPDxGNPS curated and protonated library annotation with score > 0. PCA plotting for UNPDxGNPS protonated aborted.")
                            break
                    else:
                        output_name_plot = output_name
                    new_data_unpd_desc = pd.merge(new_data_desc.loc[new_data_desc.curated_lib_annotation_origin == 'UNPD'],
                             data_reference_descriptors.loc[:, ["SMILES"] + set_reference_descriptors_bestCos_top24],
                             left_on="tremolo_SMILES_best", right_on="SMILES", how="left")
                    new_data_gnps_desc = gnps_descriptors_result.loc[gnps_descriptors_result.gnps_Smiles.isin(
                        new_data_desc.gnps_Smiles[new_data_desc.curated_lib_annotation_origin == 'GNPS']),:]
                    
                    X_new = pd.concat([new_data_unpd_desc.iloc[:, 5:], new_data_gnps_desc.iloc[:, 4:]], axis=0)
                    X_new.loc[:,"point_labels"] = ['UNPD'] * new_data_unpd_desc.shape[0] + ['GNPS'] * new_data_gnps_desc.shape[0]
                    # remove rows with NaNs
                    valid_smiles_rows = (np.isnan(X_new.iloc[:,:-1]).sum(1) == 0)
                    X_new_filtered = X_new.loc[valid_smiles_rows, :].copy()
                    point_labels_unpd_gnps = X_new_filtered.loc[:, "point_labels"]
                    X_new_filtered.drop(["point_labels"], axis=1, inplace=True)
                    X_new_scaled = scaler.fit_transform(X_new_filtered)
                    X_transformed = pca.transform(X_new_scaled)
                    # now plot reference PCA without components as background for the new data
                    biplot_scatter_reference_PCA_wComponents(scores, point_labels, pvars=pvars, components=None,
                                                             feature_labels=None)
                    # plot the result new data in the PCA reference and save the result without components
                    biplot_save_new_PCA_wComponents(output_path, output_name_plot, data_type, scores=X_transformed,
                                                    point_labels=point_labels_unpd_gnps,
                                                    pvars=pvars)
                    # then save it with components
                    biplot_scatter_reference_PCA_wComponents(scores, point_labels, pvars=pvars, components=None,
                                                             feature_labels=None)
                    biplot_save_new_PCA_wComponents(output_path, output_name_plot, data_type, scores=X_transformed,
                                                    point_labels=point_labels_unpd_gnps,
                                                    pvars=pvars, components=arrows, feature_labels=features)
    else: # external table
        # make PCA with external data containing the descriptors columns in the provided table
        # here output_path is the directory to save the outputs
        # filter only the descriptors and data_type column, which is equal to the provided type
        output_path.mkdir(exist_ok=True)
        if not output_path.exists() or not output_path.is_dir():
            sys.exit("The final output directory to store the PCA plots could not",
                  "be created. PCA plotting for external data type (",data_type, ") aborted.")
        output_name_plot = output_name
        # check if new data contains the top 24 descriptors columns and the data_type column
        if data_type not in new_data.columns:
            print("WARNING: The provided data_type is not present in the columns of the provided smiles_descriptors_table_path.",
                  "All the points will be colored with the same color and receive the same name equals to", data_type, ".")
            points_type = np.asarray([data_type] * new_data.shape[0])
        else:
            points_type = new_data.loc[:, data_type].values
            # convert any NA value to string
            if pd.isna(points_type).any():
                points_type[pd.isna(points_type)] = "NA"
            if np.unique(points_type).size > 8:
                print("WARNING: The provided data_type contains more than 8 different classes,",
                      "which is the maximum allowed number of classes. The points will be colored with the same color",
                      "instead and receive the same name equals to ", data_type, ".")
                points_type = np.asarray([data_type] * new_data.shape[0])
        if not pd.Series(set_reference_descriptors_bestCos_top24).isin(new_data.columns).all():
            sys.exit("The provided smiles_descriptors_table_path do not contain all the top 24 descriptors columns and thus",
                  "can not be transformed to the NP3 chemical space. PCA plotting for external data type (", data_type,
                  ") aborted.")
        # filter only the descriptors columns
        X_new = new_data.loc[:,set_reference_descriptors_bestCos_top24]
        # remove rows with NaNs
        valid_smiles_rows = (np.isnan(X_new).sum(1) == 0)
        X_new_filtered = X_new.loc[valid_smiles_rows, :]
        if X_new_filtered.shape[0] == 0:
            sys.exit("  - No valid descriptors in the provided external table, there are NAN in all rows.",
                "PCA plotting for external data type (", data_type, ") aborted.")
        if (~valid_smiles_rows).any():
            points_type = points_type[valid_smiles_rows]
            invalid_smiles = (~valid_smiles_rows).sum()
            print("    - There were a total of", invalid_smiles, "invalid entry with at least one NaN descriptor value.",
                  "The following rows (skipping the header) were removed from PCA:", np.where(~valid_smiles_rows)[0]+1)
        # scale and transform the new data
        X_new_scaled = scaler.fit_transform(X_new_filtered)
        X_transformed = pca.transform(X_new_scaled)
        # now plot reference PCA without components as background for the new data
        biplot_scatter_reference_PCA_wComponents(scores, point_labels, pvars=pvars, components=None,
                                                 feature_labels=None)
        # plot the result new data in the PCA reference and save the result without components
        biplot_save_new_PCA_wComponents(output_path, output_name_plot, data_type, scores=X_transformed,
                                        point_labels=points_type,
                                        pvars=pvars)
        # then save it with components
        biplot_scatter_reference_PCA_wComponents(scores, point_labels, pvars=pvars, components=None,
                                                 feature_labels=None)
        biplot_save_new_PCA_wComponents(output_path, output_name_plot, data_type, scores=X_transformed,
                                        point_labels=points_type,
                                        pvars=pvars, components=arrows, feature_labels=features)
    #
    print("  - Creating the PCA quality of representation circle plot with the reference components")
    # plot the PCA quality of representation (cos2) - circle plot with components
    # scale components loadings
    loadings_scaled = loadings * np.sqrt(pca.explained_variance_)
    gradient_color_cos2 = [(0.0, "#00AFBB"), (0.5, "#00AFBB"), (0.75, "#E7B800"), (1.0, "#FC4E07")]
    color_cos2_cmap = mcolors.LinearSegmentedColormap.from_list("OrRd_r", gradient_color_cos2)
    plot_arrows_correlation_circle(output_path, loadings_scaled, features, pvars, color_cos2_cmap)


def pca_calculation_mz_ref_plot(clean_data_path, output_path, output_name):
    # check if files exists
    clean_data_path = Path(clean_data_path)
    output_path = Path(output_path)

    if not clean_data_path.exists() or not clean_data_path.is_file():
        sys.exit("The provided path to the clean data file, to be used to calculate its reference mz PCA, does not exists." +
                 " PCA plotting aborted.")
    if not output_path.exists() or not output_path.is_dir():
        sys.exit("The provided output path does not exists. PCA plotting aborted.")
    
    # Original data
    # X_original is read here
    # TODO read the complete table and them filter the descriptor cols and the type cols - no filter for now - beta
    clean_data_mz = pd.read_csv(clean_data_path,
                                             low_memory=False,
                                             usecols=["mzConsensus","rtMean","basePeakInt","sumInts","maxArea"])
    print("* Creating the NP3 reference PCA from m/zs for plotting *")
    # not removing duplicates here - the clean data contains the minimum number of fragmented clusters
    X_original = clean_data_mz
    
    # Get valid Calculations cols
    valid_descriptors_cols = []
    for col in X_original.columns:
        try:
            if (X_original[col].var() > 0):
                valid_descriptors_cols.append(col)
        except:  # rm NAs descriptors
            pass
    
    # remove zero variance cols
    X_filtered = X_original.loc[:, valid_descriptors_cols]
    # remove rows with NaNs
    valid_mzs_rows = (np.isnan(X_filtered).sum(1) == 0)
    X_filtered = X_filtered.loc[valid_mzs_rows, :]
    
    # filter label in Type
    # TODO define Type to use here, blanks or beds or none - no filter for now - beta testing
    #point_labels = data_reference_descriptors.loc[valid_smiles_rows, "Type"]
    point_labels = ["mzs"]*X_filtered.shape[0]
    
    # extract final features
    features = X_filtered.columns.values
    
    # standardization
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_filtered)
    
    # PCA
    pca = PCA(n_components=2).fit(X_scaled)
    X_reduced = pca.transform(X_scaled)
    
    # Now pca.components_ holds the calculated eigenvectors
    # and pca.explained_variance_ holds the eigenvalues
    # proportions of variance explained by axes
    pvars = pca.explained_variance_ratio_ * 100
    print(f"  - Explained variance: PC1 {pvars[0]:.2f}% and PC2 {pvars[1]:.2f}%")
    print(f"  - Cumulative: {np.cumsum(pvars)[1]:.2f}%")
    
    # coordinates of samples (i.e., scores; let's take the first two axes)
    scores = X_reduced
    
    # coordinates of features - descriptors/components (i.e., loadings; note the transpose)
    loadings = pca.components_.T
    
    # Scale the features (arrows) properly to match the samples (points).
    # The following code scales by the maximum absolute value of samples on each axis.
    arrows = loadings * np.abs(scores).max(axis=0)
    
    # create default parameters for plotting
    params = {'legend.fontsize': 'x-large',
              'legend.title_fontsize': 'x-large',
              'figure.figsize': (9, 9),
              'axes.labelsize': 'x-large', 'axes.labelweight': 'bold', 'axes.titlesize': 'x-large',
              'xtick.labelsize': 'x-large', 'ytick.labelsize': 'x-large',
              'font.weight': 'bold'}
    pylab.rcParams.update(params)  # set styling

    # use the following code to plot the default PCA with components
    biplot_save_reference_PCA_wComponents(output_path, output_name+"_mz_quantification_", scores, point_labels, pvars, arrows, features)
    #
    print("  - Creating the PCA quality of representation circle plot with the reference components")
    # plot the PCA quality of representation (cos2) - circle plot with components
    # scale components loadings
    loadings_scaled = loadings * np.sqrt(pca.explained_variance_)
    gradient_color_cos2 = [(0.0, "#00AFBB"), (0.5, "#00AFBB"), (0.75, "#E7B800"), (1.0, "#FC4E07")]
    color_cos2_cmap = mcolors.LinearSegmentedColormap.from_list("OrRd_r", gradient_color_cos2)
    plot_arrows_correlation_circle(output_path, loadings_scaled, features, pvars, color_cos2_cmap)

if __name__ == "__main__":
    np3_chemical_space_reference_path = Path(__file__).resolve().parent / "Chemical_space_data" / "descriptors_reference_unpd_drugbank_allo_rev_natural_pubmedID_clean_top24.csv"
    if len(sys.argv) > 4:
        smiles_descriptors_table_path = sys.argv[1]
        type_column = sys.argv[2]
        output_path = sys.argv[3]
        output_name = sys.argv[4]
    else:
        print("Error: Four arguments must be supplied to created a PCA plot from an external table containing data_types and the top 24 descriptors using the NP3 reference datasets to create the chemical space:\n"
            "  1 - smiles_descriptors_table_path: Path to the external table containing the top 24 descriptors columns "
              "and one optional column for labeling the points (.csv);\n"
            "  2 - data_type: a name defining the type of the data to be used to label the points or the name of a "
              "column from smiles_descriptors_table_path containing the types/categories names of each entry "
              "to be used to label the points in the final PCA plot. This name must be different from 'UNPD' or 'GNPS' and, "
              "in the case of a column, it must contain at most 8 different classes names, otherwise only the provided "
              "data_type name will be used to label the points. These column values will go to the legend of the PCA;\n"
            "  3 - output_path: the path to the directory to store the PCA plots;\n"
            "  4 - output_name: the output name to be used for naming the output plots.\n")
        sys.exit(1)
    # call PCA creation
    pca_calculation_smiles_rcdk_ref_plot(np3_chemical_space_reference_path, smiles_descriptors_table_path, output_path, output_name, data_type=type_column)
    