import networkx as nx
import pandas as pd
import numpy as np
import sys, os

def output_mn_annotations(G, filename):
    # add selfloops
    G.add_edges_from([(isolated_node,isolated_node,{"cosine": 1.0, "annotation": "", "mzError": 0.0, "rtError": 0.0,
                                   "numCommonSamples": ""}) for isolated_node in nx.isolates(G)])
    #
    output_file = open(filename, "w")
    #Outputting the graph
    component_index = 0
    #write header
    output_file.write(",".join(["msclusterID_source","msclusterID_target","cosine","annotation",
                                "mzError","rtError","numCommonSamples", "componentIndex"]) + "\n")
    #write edges
    for component in nx.weakly_connected_components(G):
        component_index += 1
        for edge in G.subgraph(component).edges(data=True):
            output_list = []
            output_list.append(str(edge[0]))
            output_list.append(str(edge[1]))
            output_list.append("{:.3f}".format(edge[2]["cosine"]))
            output_list.append(str(edge[2]["annotation"]))
            output_list.append("{:.3f}".format(edge[2]["mzError"]))
            output_list.append("{:.3f}".format(edge[2]["rtError"]))
            output_list.append(str(edge[2]["numCommonSamples"]))
            output_list.append(str(component_index))
            output_file.write(",".join(output_list) + "\n")

def mn_annotation_find_protonated(mn_annotation_edges_file, nodes_filename):
    if not os.path.isfile(mn_annotation_edges_file):
        sys.exit("ERROR. The provided molecular network of annotations edges file '"+mn_annotation_edges_file+"' does not exists.")
    if not os.path.isfile(nodes_filename):
        sys.exit("ERROR. The provided clean table file '"+nodes_filename+"' does not exists.")

    # read the ivamn and fix the msclusterID columns types to int
    ann_edges_table = pd.read_csv(mn_annotation_edges_file,
                                  dtype={'msclusterID_source': np.int64, 'msclusterID_target':np.int64})

    # create undirected graph from edges list and remove selfloops
    G = nx.DiGraph()
    G.add_edges_from([(str(x[0]), str(x[1]), {'cosine': x[2], 'annotation': x[3],
                                              'mzError': x[4], 'rtError': x[5], "numCommonSamples": x[6]})
                      for x in ann_edges_table.itertuples(index=False)])
    G.remove_edges_from(nx.selfloop_edges(G))
    del ann_edges_table
    # execute hits algorithm
    # #The HITS algorithm computes two numbers for a node.
    # #Authorities estimates the node value based on the incoming links.
    # #Hubs estimates the node value based on outgoing links.
    # hits = nx.hits(G, normalized=False, max_iter=G.number_of_nodes())
    # Return all nodes having a path to \(source\) in G.
    number_ancestors = {str(node_name): len(nx.ancestors(G, str(node_name))) for node_name in list(G.nodes)}
    # execute page rank algorithm, use the number of ancestors as personalization if there is at least one node if an ancestor
    if (np.asarray(list(number_ancestors.values())) > 0).any():
        pagerank = nx.pagerank(G, max_iter=G.number_of_nodes()*2, tol=1e-07, personalization=number_ancestors)
    else:
        pagerank = nx.pagerank(G, max_iter=G.number_of_nodes()*2, tol=1e-07)


    G_nodes_info = pd.read_csv(nodes_filename, usecols=["msclusterID", "mzConsensus", "rtMin", "rtMax", "rtMean",
                                                        "multicharge_ion", "isotope_ion"])

    G_nodes_info['in_degree'] = -1
    G_nodes_info['total_degree'] = -1
    # G_nodes_info['hub'] = -1
    # G_nodes_info['authority'] = -1
    G_nodes_info['pagerank'] = -1
    G_nodes_info['number_ancestors'] = -1
    G_nodes_info['protonated_representative'] = -1
    # assign the link analysis algorithm scores to the graph nodes attributes and set protonated_representative as zero
    for i in range(G_nodes_info.shape[0]):
        #print(i)
        node_name = str(G_nodes_info.loc[i, 'msclusterID'])
        G_nodes_info.loc[i, 'in_degree'] = G.in_degree(node_name)
        G_nodes_info.loc[i, 'total_degree'] = G.degree(node_name)
        # G_nodes_info.loc[i, 'hub'] = hits[0][node_name]
        # G_nodes_info.loc[i, 'authority'] = hits[1][node_name]
        G_nodes_info.loc[i, 'pagerank'] = pagerank[node_name]
        G_nodes_info.loc[i, 'number_ancestors'] = number_ancestors[node_name]
        if not G.nodes.get(node_name) is None:
            G.nodes[node_name]['pagerank'] = G_nodes_info.loc[i, 'pagerank']
            G.nodes[node_name]['protonated_representative'] = 0
            G.nodes[node_name]['multicharge_ion'] = G_nodes_info.loc[i, 'multicharge_ion']
            G.nodes[node_name]['isotope_ion'] = G_nodes_info.loc[i, 'isotope_ion']
        else:
            sys.exit("ERROR. A inconsistency was found between the MN of annotations and the clean table. Node "+
                     node_name+ " is not present in the edges list")
            # G.add_nodes_from([(node_name, {'pagerank': G_nodes_info.loc[i, 'pagerank'],
            #                                'protonated_representative': 0})])

    # for each component assign a node as protonated_representative of the possible protonated ion
    # based on the maximum score of the pagerank algorithm. If more than one node have a maximum score,
    # deals with possible ambiguities due to cycles if no cycles is present leave them all as representatives.
    # Remove all the nodes ancestors to the protonated and repeat the process until there is no
    # remaining node in the component
    neutral_losses = ['[M+H-NH3]+', '[M+H-H2O]+', '[M+H-NH3-H2O]+']
    for i,component in enumerate(nx.weakly_connected_components(G)):
        #print(i)
        Gk = G.subgraph(component).copy()
        while Gk.number_of_nodes() > 0:
            nodes_pagerank = {node: Gk.nodes[node]['pagerank'] for node in list(Gk.nodes())}
            max_pagerank = max(nodes_pagerank.values())  # get the maximum score
            # node_star = max(list(Gk.nodes()), key=lambda i: Gk.nodes[i]['pagerank'])
            # return all nodes that have the maximum score
            node_star = [node for node,score in nodes_pagerank.items() if score == max_pagerank]
            if len(node_star) > 1:
                # get the cycles between the node_star[0] and the other node_stars
                node_star_neighbors = [neighbor for neighbor in Gk[node_star[0]] if (neighbor in node_star and node_star[0] in Gk[neighbor])]
                # if no cycles assign node_star[0] as the node_star
                # else assign the target of the minimum mzError edge in the cycles as the node_star
                if len(node_star_neighbors) == 0:
                    node_star = node_star[0]
                else:
                    Gk_multiple_star = Gk.subgraph(node_star_neighbors+[node_star[0]])
                    min_mzError_edge = min(Gk_multiple_star.edges(data=True), key=lambda edge: edge[2]['mzError'])
                    # assign node_star as the target of the edge with the minimum mzError
                    node_star = min_mzError_edge[1]
                    del Gk_multiple_star
            else:
                node_star = node_star[0]
            # check if node_star is a multicharge ion or and isotope ion or
            # have a neutral loss predecessor that is a multicharge ion or and isotope ion,
            # if yes remove it and its ancestors without assign it as a representative
            node_star_ancestors = nx.ancestors(G, node_star)
            node_star_ancestors.add(node_star)
            predecessors_neutral_losss_multicharge_isotope_ion = [predecessor for predecessor in G.predecessors(node_star)
                                                          if ((neutral_losses[0] in G[predecessor][node_star]['annotation'] or
                                                               neutral_losses[1] in G[predecessor][node_star]['annotation'] or
                                                               neutral_losses[2] in G[predecessor][node_star]['annotation']) and
                                                              ((G.nodes[predecessor]['multicharge_ion'] > 0) |
                                                               (G.nodes[predecessor]['isotope_ion'] > 0)))]
            if (G.nodes[node_star]['multicharge_ion'] > 0 or G.nodes[node_star]['isotope_ion'] > 0 or
                len(predecessors_neutral_losss_multicharge_isotope_ion) > 0):
                # also propagate the multicharge_ion and the isotope ion flags to its successors
                # that are neutral losses and not the monocharge or monoisotopic ion
                successors_neutral_loss = [neighbor for neighbor in Gk[node_star]
                                           if (neutral_losses[0] in Gk[node_star][neighbor]['annotation'] or
                                               neutral_losses[1] in Gk[node_star][neighbor]['annotation'] or
                                               neutral_losses[2] in Gk[node_star][neighbor]['annotation'])]
                if len(successors_neutral_loss) > 0:
                    # set both variant ions as 1 to remove the successors in both cases
                    for successor in successors_neutral_loss:
                        Gk.nodes[successor]['multicharge_ion'] = 1
                        Gk.nodes[successor]['isotope_ion'] = 1
                # remove the ancestors nodes of the multicharge ion and itself, to not assign them as representatives
                Gk.remove_nodes_from(node_star_ancestors)
            else:
                # assign the node_star, the node with the highest pagerank score in the current component,
                # as representative of the protonated ion
                G.nodes[node_star]['protonated_representative'] = 1
                # remove the ancestors nodes of the chosen star
                Gk.remove_nodes_from(node_star_ancestors)


    G_nodes_info['protonated_representative'] = [G.nodes[str(node)]['protonated_representative'] for node in G_nodes_info.msclusterID]

    print("\n* Molecular Network of Annotations Statistics *\n")
    print(nx.info(G))
    print("Number of components:", nx.number_weakly_connected_components(G))
    print("Number of isolated nodes: ", (G_nodes_info.total_degree == 0).sum())
    print("Number of protonated_representative nodes (with isolated nodes):", (G_nodes_info.protonated_representative == 1).sum())

    # for each m/z assigned as a protonated representative that is not an isolated node, compute the mz and rt error
    # associated with all its ancestors paths
    G_nodes_info['protonated_mzError_sum'] = 0.0
    G_nodes_info['protonated_rtError_sum'] = 0.0
    G_nodes_info['protonated_num_ancestors_edges'] = 0
    for idx, protonated_id in G_nodes_info[(G_nodes_info.protonated_representative == 1) &
                                           (G_nodes_info.total_degree > 0)].msclusterID.iteritems():
        protonated_id = str(protonated_id)
        protonated_ancestors = nx.ancestors(G, protonated_id)
        protonated_ancestors.add(protonated_id)
        Gk = G.subgraph(protonated_ancestors)
        protonated_mzerror = [edge[2]['mzError'] for edge in Gk.edges(data=True)]
        # pd.Series(protonated_mzerror).describe()
        protonated_rterror = [edge[2]['rtError'] for edge in Gk.edges(data=True)]
        G_nodes_info.loc[idx, 'protonated_mzError_sum'] = np.round(np.sum(protonated_mzerror), 3)
        G_nodes_info.loc[idx, 'protonated_rtError_sum'] = np.round(np.sum(protonated_rterror), 2)
        G_nodes_info.loc[idx, 'protonated_num_ancestors_edges'] = len(protonated_rterror)

    # save the mn of annotations attributes in a new table file
    G_nodes_info.to_csv(os.path.splitext(mn_annotation_edges_file)[0]+"_attributes.csv",
                        index=False)
    # overwrite the mn annotation edges file
    output_mn_annotations(G, mn_annotation_edges_file)

    # add the protonated_representative columns to the clean tables
    # G_nodes_info = G_nodes_info[['msclusterID', 'protonated_representative', 'protonated_mzError_sum', 'protonated_rtError_sum']]
    clean_table = pd.read_csv(nodes_filename, low_memory=False)
    clean_table[['protonated_representative', 'protonated_mzError_sum', 'protonated_rtError_sum']] = G_nodes_info[['protonated_representative', 'protonated_mzError_sum', 'protonated_rtError_sum']]
    clean_table.to_csv(nodes_filename, index=False)
    if 'peak_area' in nodes_filename:
        nodes_filename = nodes_filename.replace('peak_area', 'spectra')
        clean_table = pd.read_csv(nodes_filename, low_memory=False)
        clean_table[['protonated_representative', 'protonated_mzError_sum', 'protonated_rtError_sum']] = G_nodes_info[['protonated_representative', 'protonated_mzError_sum', 'protonated_rtError_sum']]
        clean_table.to_csv(nodes_filename, index=False)
    else:
        nodes_filename = nodes_filename.replace('spectra', 'peak_area')
        clean_table = pd.read_csv(nodes_filename, low_memory=False)
        clean_table[['protonated_representative', 'protonated_mzError_sum', 'protonated_rtError_sum']] = G_nodes_info[
            ['protonated_representative', 'protonated_mzError_sum', 'protonated_rtError_sum']]
        clean_table.to_csv(nodes_filename, index=False)


if __name__ == "__main__":
    if len(sys.argv) > 2:
        # print(sys.argv)
        mn_annotation_edges_file = sys.argv[1]
        nodes_filename = sys.argv[2]
    else:
        print("Error: Two arguments must be supplied to assign the protonated ion candidates in the molecular networking of annotations:\n",
        " 1 - Path to the molecular networking of annotations file ('<job_name>_molecular_networking_annotations.selfloop');\n",
        " 2 - Path to the clean table count with the final list of consensus spectra.")
        sys.exit(1)

    mn_annotation_find_protonated(mn_annotation_edges_file, nodes_filename)
