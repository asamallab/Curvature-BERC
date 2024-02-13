#!/usr/bin/env python
# coding: utf-8

#--------------------------------------------------------------------
# Calculates Robustness of a given network G
#--------------------------------------------------------------------

import random
import pandas as pd
import networkx as nx


def Efficiency(G):
    '''
    Input: G: a networkx graph;
    Output: efficiency of that graph 
    '''
    N = len(G)
    if N < 2:
        return 0
    inv_lengths = []
    for node in G:
        lengths = nx.single_source_dijkstra_path_length(G, node)
        inv = [1.0/x for x in lengths.values() if x !=  0]
        inv_lengths.extend(inv)
    return sum(inv_lengths)/(N*(N-1))


def Robustness_node(Dict, Graph, sort_order):
    '''
    Inputs: Dict: dictionary of vertices and its corresponding vertex measure, Graph: a networkx graph, sort_order (True/False): True for graph measures and False for curvature measures
    Outputs: node_removed: list of fraction of removed nodes at each step, efficiency: list of efficiency of the graph after removing each node
    '''
    efficiency=[]
    node_removed=[]
    # Calculating total number of nodes in the graph
    total_no_of_node = Graph.number_of_nodes()
    # Sort the dictionary by key and iterating
    no=1.0
    for c1 in sorted(Dict.items(),key=lambda t:t[1], reverse = sort_order):
        # Removing the node
        Graph.remove_node(c1[0])
        # Calculating number of nodes of the new graph
        no_of_node = Graph.number_of_nodes()
        # Ratio of edges removed and efficiency
        eff = Efficiency(Graph)
        efficiency.append(eff)
        frac_of_node_removed = float(total_no_of_node-no_of_node)/total_no_of_node
        node_removed.append(frac_of_node_removed)
        print(c1[0], ' is removed.', Graph, eff)
        no+=1
    return node_removed, efficiency


def Robustness_random(Graph):
    '''
    Inputs: Graph: a networkx graph,
    Outputs: node_removed: list of fraction of removed nodes at each step, mean_efficiency: list of mean efficiency 
             (repeted 10 times) of the graph after removing each node
    
    '''
    # Calculating total number of nodes in the graph
    total_no_of_node = Graph.number_of_nodes()
    allNodes = list(Graph.nodes())
    efficiencies = []
    for k in range(10):
        efficiency=[]
        node_removed=[]
        no=1.0
        random.shuffle(allNodes)
        nodes = allNodes
        G = Graph.copy()
        for i in nodes:
            # Removing the node
            G.remove_node(i)
            # Calculating number of nodes of the new graph
            no_of_node = G.number_of_nodes()
            # Ratio of edges removed and efficiency
            eff = Efficiency(G)
            efficiency.append(eff)
            frac_of_node_removed = float(total_no_of_node-no_of_node)/total_no_of_node
            node_removed.append(frac_of_node_removed)
            print(i, ' is removed.', G,eff)
            no+=1
        efficiencies.append(efficiency)
        print('done for random sample ', k)
    allefficiency = pd.DataFrame(efficiencies)   
    allefficiency = allefficiency.T
    mean_efficiency = list(allefficiency.mean(axis = 1))
    return node_removed, mean_efficiency

