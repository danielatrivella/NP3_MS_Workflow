import networkx as nx
import pandas as pd
import numpy as np
from itertools import chain
import sys, os
from mn_annotations_assign_protonated_representative import output_mn_annotations
from molecular_network_filtering_library import loading_network, print_net_info, filter_component, filter_top_k, filter_min_matched_peaks, add_selfloops, output_graph
"""
This script creates de SSMN [M+H]+ filtered and the IVAMN [M+H]+ networks without blanks. 
It uses as input the complete SSMN and the IVAMN networks, the IVAMN attribute table (protonated information) and the clean table. 
It works as follow: 
1. In IVAMN the blank nodes are removed together with its neighbours if argument blank_expansion != 0, else only remove the blank nodes;
2. Then the [M+H]+ are selected in the IVAMN with no blanks, resulting in the IVAMN [M+H]+ network;
3. The nodes from IVAMN [M+H]+ are used to select the final nodes from the SSMN, resulting in the SSMN [M+H]+ network;
4. Finally, the SSMN [M+H]+ is filtered using the min_matched_peaks, top_k and max_component size arguments, resulting in the SSMN [M+H]+ filtered.
"""

def load_direct_network(edges_file):
    ann_edges_table = pd.read_csv(edges_file)
    #
    # create undirected graph from edges list and remove selfloops
    G = nx.DiGraph()
    # source and target nodes are in the first two columns
    G.add_edges_from([(x[0], x[1], {'cosine': x[2], 'annotation': x[3],
                                              'mzError': x[4], 'rtError': x[5], "numCommonSamples": x[6]})
                      for x in ann_edges_table.itertuples(index=False)])
    G.remove_edges_from(nx.selfloop_edges(G))
    del ann_edges_table
    return G

def flatten_list_of_lists(lsts):
    return list(chain.from_iterable(lsts))

def create_protonated_SSMN_IVAMN(ivamn_file, ssmn_file, clean_table_file, blank_expansion, top_k,
                                 max_component_size, min_matched_peaks):
    if not os.path.isfile(ssmn_file):
        sys.exit("ERROR. The provided molecular network of similarities (SSMN) edges file '" + ssmn_file +
                 "' does not exists.")
    if not os.path.isfile(clean_table_file):
        sys.exit("ERROR. The provided clean table file '" + clean_table_file + "' does not exists.")

    if not os.path.isfile(ivamn_file):
        sys.exit("ERROR. The provided molecular network of annotations (IVAMN) edges file '" + ivamn_file +
              "' does not exists. Thus, all nodes are considered protonated nodes.")


    ## Create the IVAMN [M+H]+
    print("  * Creating the IVAMN [M+H]+ network *")
    # load ivamn
    ivamn = load_direct_network(ivamn_file)

    # 1. remove blanks from IVAMN using the blank_expansion,
    # remove the blank nodes ancestors without considering the edges direction
    print("\t1. Remove blanks from IVAMN using the blank_expansion")
    if 'BLANKS_TOTAL' in pd.read_csv(clean_table_file, nrows=1).columns:
        # read the clean table
        clean_table = pd.read_csv(clean_table_file,
                                  usecols=["msclusterID", "BLANKS_TOTAL", "protonated_representative"])
        blank_nodes_id = clean_table.loc[clean_table["BLANKS_TOTAL"] > 0, "msclusterID"].values
        uivamn = ivamn.to_undirected()
        nodes_to_remove = list(blank_nodes_id)
        if blank_expansion == -1:
            # remove all ancestors of the blank nodes
            for blank in blank_nodes_id:
                nodes_to_remove += list(nx.ancestors(uivamn, blank))
        elif blank_expansion > 0:
            # remove the ancestors of the blank nodes within a distance equal the blank_expansion
            for blank in blank_nodes_id:
                nodes_to_remove += flatten_list_of_lists([a[1] for a in nx.bfs_successors(uivamn, blank, blank_expansion)])
        # else blank_expansion == 0: only remove the blank nodes
        # remove duplicates
        nodes_to_remove = np.unique(nodes_to_remove)
        ivamn.remove_nodes_from(nodes_to_remove)
        del(uivamn)
        # filter in the clean table the not blank nodes from the IVAMN
        clean_table = clean_table.loc[np.isin(clean_table.msclusterID.values, np.array(ivamn.nodes())), :]
    else:
        # read the clean table
        clean_table = pd.read_csv(clean_table_file,
                                  usecols=["msclusterID", "protonated_representative"])

    # 2. Select the [M+H]+ ions in the remaining IVAMN
    print("\t2. Select the [M+H]+ ions in the remaining IVAMN")
    # filter the [M+H]+ nodes in ivamn
    nodes_to_remove = clean_table.loc[clean_table.protonated_representative == 0, "msclusterID"].values
    ivamn.remove_nodes_from(nodes_to_remove)
    # save the IVAMN [M+H]+ network
    output_mn_annotations(ivamn, ivamn_file.replace(".selfloop", "_protonated.selfloop"))
    print("\n    Done for IVAMN [M+H]+!")
    print_net_info(ivamn, "IVAMN [M+H]+")
    del ivamn

    # filter the [M+H]+ nodes in clean table
    clean_table = clean_table.loc[clean_table.protonated_representative == 1, :]

    ## Create the filtered SSMN [M+H]+
    print("  * Creating the SSMN [M+H]+ network *")
    # 3. Use the nodes from IVAMN [M+H]+ to filter the SSMN and obtain the SSMN [M+H]+ network;
    print("\t3. Use the nodes from IVAMN [M+H]+ to filter the SSMN")
    protonated_nodes = clean_table.msclusterID.values
    ssmn = loading_network(ssmn_file, remove_selfloops=True)
    ssmn.remove_nodes_from([n for n in ssmn if n not in protonated_nodes])

    # 4. Now, filter the SSMN [M+H]+ using the top_k, max_component size and min_matched_peaks arguments,
    # resulting in the SSMN [M+H]+ filtered.
    # filter min number of matched peaks
    print("\t4. Filter the SSMN [M+H]+ using the top_k, max_component size and min_matched_peaks")
    filter_min_matched_peaks(ssmn, min_matched_peaks)
    # filter top k and component size
    filter_top_k(ssmn, top_k)
    filter_component(ssmn, max_component_size)
    # Get single nodes and create selfloops
    add_selfloops(ssmn)
    # save result
    ssmn_protonated_file = str(ssmn_file.replace('.selfloop', '_protonated') +
                  '_mmp_' + str(min_matched_peaks) +
                  '_k_' + str(top_k) + '_x_' + str(max_component_size) +
                  '.selfloop')
    # Export the final SSMN [M+H]+ filtered
    output_graph(ssmn, ssmn_protonated_file)
    print("\n    Done for SSMN [M+H]+ filtered!")
    print_net_info(ssmn, "SSMN [M+H]+ filtered")
    return


if __name__ == "__main__":
    import sys, os
    if len(sys.argv) > 7:
        # print(sys.argv)
        ivamn_file = sys.argv[1]
        ssmn_file = sys.argv[2]
        clean_table = sys.argv[3]
        blank_expansion = int(sys.argv[4])
        top_k = int(sys.argv[5])
        max_component_size = int(sys.argv[6])
        min_matched_peaks = int(sys.argv[7])
    else:
        print("Error: Eight arguments must be supplied to created the SSMN [M+H]+ filtered and the IVAMN [M+H]+ networks:\n",
            " 1 - ivamn_file: Path to the molecular networking of annotations network IVAMN file (.selfloop);\n",
            " 2 - ssmn_file: Path to the complete molecular networking of similarity SSMN file (.selfloop), before "
            "filters - if the filtered SSMN is informed the resulting [M+H]+ network may be more fragmented;\n",
            " 3 - clean_table: Path to the clean table with the final list of consensus spectra - the nodes information "
            "of the networks (.csv);\n"
            " 4 - blank_expansion: an int with the distance of neighborhood nodes from the blanks in IVAMN to be "
            "selected for removal in the final protonated networks. (0) to only remove blanks nodes,"
            "(1) to remove nodes directly connected to a blank node, "
            "(2 or greater) to remove nodes in a distance equal to 2 or greater from a blank node, or "
            "(-1) to remove all possible neighbours and ancestors of a blank node - remove blank clusters - from IVAMN;\n"
            " 5 - net_top_k: the maximum number of neighbor nodes for one single node in SSMN [M+H]+. "
            "The edge between two nodes are kept only if both nodes are within each other's TopK most similar nodes. "
            "For example, if this value is set at 20, then a single node may be connected to up to 20 other nodes. "
            "Keeping this value low makes very large networks (many nodes) much easier to visualize. ;\n",
            " 6 - maximum_component_size: the maximum number of nodes that all components of the SSMN [M+H]+ must have. "
            "The edges of the network will be removed using an increasing cosine threshold until each component has at "
            "most maximum_component_size nodes;\n",
            " 7 - min_matched_peaks: The minimum number of common fragment ions that two separate consensus MS/MS "
            "spectra must share in order to be connected by an edge in the SSMN [M+H]+")
        sys.exit(1)

    create_protonated_SSMN_IVAMN(ivamn_file, ssmn_file, clean_table, blank_expansion, top_k,
                                 max_component_size, min_matched_peaks)
