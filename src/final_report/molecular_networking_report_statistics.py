#!/usr/bin/python

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import networkx as nx
import os, sys

# Add the parent directory to sys.path
current_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(parent_dir)

# Now you can import from parent_script.py
from molecular_network_filtering_library import loading_network
from mn_annotations_assign_protonated_representative import loading_direct_network

# compute some networking statistics for all the molecular networks present in the provided output_mn_path
# and save all report data in the mn_report_path
# output_mn_path points to the molecular_networking folder inside the final result folder
# of job named output_name
def compute_mn_report_statistics(output_mn_path, mn_report_path):
	if not output_mn_path.exists() or not output_mn_path.is_dir():
		print(output_mn_path)
		sys.exit("The provided molecular networking output result folder do not exists. Molecular Networking report aborted.")
	
	print("  - Computing the molecular networks statistics\n")
	# for each MN present in the output result folder, compute its statistics and save then in separated tables
	# inside the mn_report_path
	#  list of possible MN and compute statistics
	for mn_file in output_mn_path.glob("*.selfloop"):
		# check if this is a ssmn (undirect) or an ivamn network (direct), and read them properly
		if mn_file.name.find("ssmn") >= 0:
			G = loading_network(mn_file)
		elif mn_file.name.find("ivamn") >= 0:
			G = loading_direct_network(mn_file)
		else:
			print("Invalid molecular networking type:", mn_file.name, ". Skipping it.")
			pass
		G.name = mn_file.name.rstrip(".selfloop")
		# compute network statistics return dict information
		# Number of nodes, edges, degree, components and isolated, and also some clustering coefficients
		mn_statistics = {'Statistics': ['Network Name'],
		                 'Value': [G.name],
		                 'Description': ['The network file name']}
		# number of nodes of G
		n_nodes = G.number_of_nodes()
		mn_statistics['Statistics'].append("Number of nodes")
		mn_statistics['Value'].append(n_nodes)
		mn_statistics['Description'].append("The number of nodes in the network, the network size.")
		n_edges = G.number_of_edges()
		mn_statistics['Statistics'].append("Number of edges")
		mn_statistics['Value'].append(n_edges)
		mn_statistics['Description'].append("The number of edges in the network.")
		# number of isolated nodes
		n_isolated = nx.number_of_isolates(G) + nx.number_of_selfloops(G)
		mn_statistics['Statistics'].append("Number of isolated nodes")
		mn_statistics['Value'].append(f"{n_isolated} ({n_isolated / n_nodes * 100:.1f}%)")
		mn_statistics['Description'].append(
			"The number of isolated nodes (or selfloops) present in the network and its percentage overall. ")
		# compute average degree of nodes
		if G.is_directed():
			# also compute in-degree and out-degree
			avg_in_degree = sum(dict(G.in_degree()).values()) / G.number_of_nodes()
			mn_statistics['Statistics'].append("Average in-degree")
			mn_statistics['Value'].append(f"{avg_in_degree:.2f}")
			mn_statistics['Description'].append("The average in-degree of the nodes in the network, directed. The in-degree of a node is equal to the number of adjacent incoming connections (incoming edges) it have.")
			avg_out_degree = sum(dict(G.out_degree()).values()) / G.number_of_nodes()
			mn_statistics['Statistics'].append("Average out-degree")
			mn_statistics['Value'].append(f"{avg_out_degree:.2f}")
			mn_statistics['Description'].append("The average out-degree of the nodes in the network, directed. The out-degree of a node is equal to the number of adjacent outgoing connections (outgoing edges) it have.")
		avg_degree = sum(dict(G.degree()).values()) / G.number_of_nodes()
		mn_statistics['Statistics'].append("Average degree")
		mn_statistics['Value'].append(f"{avg_degree:.2f}")
		mn_statistics['Description'].append("The average degree of the nodes in the network, undirected. The degree of a node is equal to the number of adjacent connections (edges) it have.  Nodes with high degrees often serve as central points or 'hubs' in a network.")
		# compute number of components
		if G.is_directed():
			n_components = nx.number_weakly_connected_components(G)
			mn_statistics['Statistics'].append("Number of connected components")
			mn_statistics['Value'].append(f"{n_components} ({n_components-n_isolated})")
			mn_statistics['Description'].append(
				"The number of components present in the network, and the number of not isolated components between parenthesis (components with size >= 2 nodes). Isolated nodes count as one isolated component. For direct network this is the weakly connected components.")
		else:  # undirected
			n_components = nx.number_connected_components(G)
			mn_statistics['Statistics'].append("Number of connected components")
			mn_statistics['Value'].append(f"{n_components} ({n_components-n_isolated})")
			mn_statistics['Description'].append(
				"The number of components present in the network, and the number of not isolated components between parenthesis (components with size >= 2 nodes). Isolated nodes count as one isolated component. ")
		# Calculate global clustering coefficient
		global_clustering = nx.transitivity(G)
		mn_statistics['Statistics'].append("Global clustering coefficient ")
		mn_statistics['Value'].append(f"{global_clustering:.2f}")
		mn_statistics['Description'].append(
			"The clustering coefficient measures the degree to which nodes in a graph tend to cluster together, it measures the tendency of a node’s neighbors to be interconnected. The global coefficient is measured for the entire network, using the transitivity computation - the fraction of all possible triangles (nodes A, B and C are interconnected) present in the network.")
		# Calculate average clustering coefficient
		avg_clustering = nx.average_clustering(G)
		mn_statistics['Statistics'].append("Average clustering coefficient ")
		mn_statistics['Value'].append(f"{avg_clustering:.2f}")
		mn_statistics['Description'].append(
			"The network average clustering coefficient measures how likely nodes in a network are to form clusters by " +
			"assessing the density of triangles. A high average clustering coefficient indicates that a network's nodes " +
			"have a strong tendency to form tightly-knit local groups, triangles. The network average clustering coefficient" +
			" has a value range from 0 to 1. A value of 1 indicates that all possible triangles are present, meaning a " +
			"highly clustered or 'clique-like' network, while a value of 0 suggests a very sparse or unclustered network " +
			"with no triangles. For a single node, the local clustering coefficient measures the density of connections " +
			"among its neighbors. It is calculated as the ratio of the actual number of edges between the node's neighbors " +
			"to the total number of possible edges between them.")
		
		# save statistics of each network
		mn_statistics_table = pd.DataFrame(mn_statistics)
		mn_statistics_table.to_csv(mn_report_path / ("molecular_networking_statistics_"+G.name+".csv"),
		                           index=False)
