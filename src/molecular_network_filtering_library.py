#!/usr/bin/python

# code modified from https://github.com/CCMS-UCSD/GNPS_Workflows/blob/master/shared_code/molecular_network_library.py
# adapted for the NP3 graph file format, and implemented a damped approach in the edges prunning

import networkx as nx
import pandas as pd

def print_net_info(G, mn_name="network"):
    print("\n  * Network "+mn_name+" info *\n")
    G.name = mn_name
    print("    "+nx.info(G).replace("\n","\n    "))
    if G.is_directed():
        print("    Number of components:", nx.number_weakly_connected_components(G))
    else: # undirected
        print("    Number of components:", nx.number_connected_components(G))
    print("    Number of isolated nodes: ", nx.number_of_isolates(G)+nx.number_of_selfloops(G),"\n")

def loading_network(filename, remove_selfloops=False):
    G_edges_list = pd.read_csv(filename)
    G_edges_list = G_edges_list.fillna('')
    G = nx.Graph()
    # read edges file and information used in the filtering
    if 'annotation' in G_edges_list.columns:
        G.add_edges_from([(x[0], x[1], {'cosine': round(x[2], 3), 'num_matched_peaks': x[3], 'annotation': x[4],
                                        'num_peaks_source': x[5], 'num_peaks_target': x[6]})
                          for x in G_edges_list.itertuples(index=False)])
    else:
        G.add_edges_from([(x[0], x[1], {'cosine': round(x[2], 3), 'num_matched_peaks': x[3], 'annotation': '',
                                        'num_peaks_source': x[4], 'num_peaks_target': x[5]})
                          for x in G_edges_list.itertuples(index=False)])

    nx.set_node_attributes(G, G_edges_list[['msclusterID_source', 'num_peaks_source']].drop_duplicates().set_index(
        'msclusterID_source').to_dict()['num_peaks_source'], "num_peaks")
    nx.set_node_attributes(G, G_edges_list[['msclusterID_target', 'num_peaks_target']].drop_duplicates().set_index(
        'msclusterID_target').to_dict()['num_peaks_target'], "num_peaks")
    #remove selfloops
    if remove_selfloops:
        G.remove_edges_from(nx.selfloop_edges(G))

    return G

def nodes_topk_cutoff(G, top_k):
    node_cutoff_score = {}
    for node in G.nodes():
        node_edges = G.edges((node), data=True)
        node_edges = sorted(node_edges, key=lambda edge: edge[2]["cosine"], reverse=True)
        # take one more to use equality in the cosine score
        # guarantees top_k max neighbours, but can remove more than top_k nodes in case of a tied score
        edges_to_keep = node_edges[:(top_k + 1)]
        if len(edges_to_keep) == 0:
            continue
        # if the node has more than top_k neighbors, set the cutoff equals the cosine of the top_k + 1 neighbor
        # otherwise set the cutoff equals the last neighbor cosine - 0.001 to keep all neighbors
        if len(edges_to_keep) == (top_k + 1):
            node_cutoff_score[node] = edges_to_keep[-1][2]["cosine"]
        else:
            node_cutoff_score[node] = edges_to_keep[-1][2]["cosine"] - 0.001
    return node_cutoff_score

def filter_top_k(G, top_k):
    print("  Filter Top_K", top_k,"neighbors")
    #Keeping only the top K scoring edges per node
    print("  Starting Number of Edges", len(G.edges()))
    #Doing this for each pair, makes sure they are in each other's top K
    # remove the edges that are not in at least one top k
    node_cutoff_score = nodes_topk_cutoff(G, top_k)
    edge_to_remove = []
    for edge in G.edges(data=True):
        if edge[0] == edge[1]:  # do not remove selfloops
            continue
        cosine_score = edge[2]["cosine"]
        threshold1 = node_cutoff_score[edge[0]]
        threshold2 = node_cutoff_score[edge[1]]
        if cosine_score <= threshold1 or cosine_score <= threshold2:
            edge_to_remove.append(edge)
    for edge in edge_to_remove:
        G.remove_edge(edge[0], edge[1])
    print("  After Top K Mutual", len(G.edges()))

def filter_min_matched_peaks(G, min_matched):
    print("  Filter min_matched_peaks", min_matched, "fragment ions")
    # Keeping only the edges that have at least min_matched peaks in common
    print("  Starting Number of Edges", len(G.edges()))
    # Doing this for each pair,
    # remove the edges that do have at least min_matched peaks
    edge_to_remove = []
    for edge in G.edges(data=True):
        if edge[0] == edge[1]:  # do not remove selfloops
            continue
        num_matched_peaks = edge[2]["num_matched_peaks"]
        min_num_peaks = min([G.nodes[edge[0]]["num_peaks"], G.nodes[edge[1]]["num_peaks"]])
        # use the minimum number of peaks when applying this filter to prevent applying it for spectra with few peaks <= min_matched
        # in this case only remove the edge if it has only one match, otherwise keep it
        if num_matched_peaks < min_matched:
            if min_num_peaks > min_matched:
                edge_to_remove.append(edge)
            elif num_matched_peaks < 2:
                edge_to_remove.append(edge)
    for edge in edge_to_remove:
        G.remove_edge(edge[0], edge[1])
    print("  After Min Matched Peak ", len(G.edges()))

def filter_component(G, max_component_size):
    if max_component_size == 0:
        return

    big_components_present = True

    i = 1
    while big_components_present == True:
        big_components_present = False
        components = nx.connected_components(G)
        for component in components:
            if len(component) > max_component_size:
                # call prunning with decreasing cutoff proportional to the big component size
                if len(component) > 3*max_component_size:
                    prune_component(G, component, 0.01)
                elif len(component) > 2*max_component_size:
                    prune_component(G, component, 0.005)
                elif len(component) > 1.5*max_component_size:
                    prune_component(G, component, 0.002)
                else:
                    prune_component(G, component, 0.001)
                big_components_present = True
        i = i + 1
    print("  After", str(i),"rounds of Component Pruning", len(G.edges()))

def get_edges_of_component(G, component):
    component_edges = {}
    for node in component:
        node_edges = G.edges((node), data=True)
        for edge in node_edges:
            if edge[0] < edge[1]:
                key = str(edge[0]) + "-" + str(edge[1])
            else:
                key = str(edge[1]) + "-" + str(edge[0])
            component_edges[key] = edge

    component_edges = component_edges.values()
    return component_edges

def prune_component(G, component, cosine_delta=0.001):
    component_edges = get_edges_of_component(G, component)

    min_score = 1000
    for edge in component_edges:
        min_score = min(min_score, edge[2]["cosine"])

    cosine_threshold = cosine_delta + min_score
    for edge in component_edges:
        if edge[2]["cosine"] < cosine_threshold:
            #print(edge)
            G.remove_edge(edge[0], edge[1])

def add_selfloops(G):
    # Get single nodes and create selfloops
    container_edge = []
    for iso in nx.isolates(G):
        if G.degree(iso) == 0:
            container_edge.append((iso, iso, {'cosine': 1, "num_matched_peaks": -1, "annotation": "",
                                              "num_peaks_source": G.nodes[iso]["num_peaks"],
                                              "num_peaks_target": G.nodes[iso]["num_peaks"]}))

    print("  Added {} selfloop edges".format(len(container_edge)))

    # Insert all new edges
    G.add_edges_from(container_edge)

def output_graph(G, filename):
    output_file = open(filename, "w")
    #Outputting the graph
    component_index = 0
    #write header
    output_file.write(",".join(["msclusterID_source","msclusterID_target","cosine",
                                "num_matched_peaks", "annotation", "num_peaks_source", "num_peaks_target",
                                "componentIndex"]) + "\n")
    #write edges
    for component in nx.connected_components(G):
        component_index += 1
        for edge in get_edges_of_component(G, component):
            output_list = []
            if int(edge[0]) < int(edge[1]):
                output_list.append(str(edge[0]))
                output_list.append(str(edge[1]))
            else:
                output_list.append(str(edge[1]))
                output_list.append(str(edge[0]))
            output_list.append(str(edge[2]["cosine"]))
            output_list.append(str(edge[2]["num_matched_peaks"]))
            output_list.append(str(edge[2]["annotation"]))
            output_list.append(str(edge[2]["num_peaks_source"]))
            output_list.append(str(edge[2]["num_peaks_target"]))
            output_list.append(str(component_index))
            output_file.write(",".join(output_list) + "\n")

if __name__ == "__main__":
    import sys, os
    if len(sys.argv) > 4:
        # print(sys.argv)
        graph_file = sys.argv[1]
        top_k = int(sys.argv[2])
        max_component_size = int(sys.argv[3])
        min_matched_peaks = int(sys.argv[4])
    else:
      print("Error: Four arguments must be supplied to filter the NP3 spectral similarity molecular networking (SSMN):\n",
       " 1 - Path to the molecular networking of similarity file (.selfloop);\n",
       " 2 - net_top_k: the maximum number of neighbor nodes for one single node. The edge between two nodes are "
       "kept only if both nodes are within each other's TopK most similar nodes. For example, if this value is set "
       "at 20, then a single node may be connected to up to 20 other nodes. Keeping this value low makes very large "
       "networks (many nodes) much easier to visualize. ;\n",
       " 3 - maximum_component_size: the maximum number of nodes that all component of the network must have. The edges"
       "of the network will be removed using an increasing cosine threshold until each component has at most "
       "maximum_component_size nodes;\n",
       " 4 - minimum_matched_peaks: The minimum number of common fragment ions that two separate consensus MS/MS spectra "
       "must share in order to be connected by an edge in the SSMN")
      sys.exit(1)

    G = loading_network(graph_file)
    total_nodes = set(G.nodes) # Keep original nodes list
    # filter min number of matched peaks
    filter_min_matched_peaks(G, min_matched_peaks)
    # filter top k and component size
    filter_top_k(G, top_k)
    filter_component(G, max_component_size)
    outName = str(graph_file.replace('.selfloop', '')+
                  '_mmp_'+str(min_matched_peaks)+
                  '_k_'+str(top_k)+'_x_'+str(max_component_size)+
                  '.selfloop')

    # Get single nodes and create selfloops
    add_selfloops(G)
    #print net info
    print_net_info(G, "SSMN filtered")

    # Export
    output_graph(G, outName)
