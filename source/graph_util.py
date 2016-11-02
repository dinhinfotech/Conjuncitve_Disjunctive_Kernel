# -*- coding: utf-8 -*-
import networkx as nx
import numpy as np
import util
from multiprocessing import Pool


def load_graph(file_path=None):
    """
    Parameters:
    file_path: path to adjacency matrix of the graph
    
    Return: An undirected, unweighted graph
    """
    
    AM = np.loadtxt(file_path)
    N = AM.shape[0]  
    g = nx.Graph()        
    g.add_nodes_from(range(N))
    list_edges = []
    for u in range(N-1):
        for v in range(u+1,N):
            w = AM[u,v]
            if w != 0.0:
                list_edges.append((u,v))
    g.add_edges_from(list_edges,{'label':""})
    return g

def load_graphs(folder_path=None):
    """
    Parameters:
    folder_path: path to folder containing adjacency matrices to of corresponding graphs
    
    Return: A list of undirected, unweighted graphs
    """
    file_paths = [folder_path + file_name for file_name in util.list_files_in_folder(folder_path)]
    graphs = [load_graph(file_path) for file_path in file_paths]
    return graphs

def create_graph(A=None):
    """ Create a undirected, unweighted graph given adjacency matrix
    
    Parameters:
    - A: adjacency matrix
    
    Return: a graph
    """
    # Load adjacency matrix    
    N = A.shape[0]  
    G=nx.Graph()        
    G.add_nodes_from(range(N))
    list_edges = []
    for u in range(N-1):
        for v in range(u+1,N):
            w = A[u,v]
            if w != 0.0:
                list_edges.append((u,v))
    G.add_edges_from(list_edges,{'label':""})
    return G
    
def create_graphs(matrices=None):
    """ Create list of undirected, unweighted graphs given adjacency matrices
    
    Parameters:
    - matrices: list of adjacency matrices
    
    Return: List of graphs
    """
    # Load adjacency matrix    
    pool = Pool(processes=len(matrices))
    graphs = pool.map(create_graph, matrices)
    pool.close()
    pool.join()
    return graphs   

def union_graphs(graphs):
    """ Union a list of graphs
    
    Parameters:
    - list_graphs: List of graphs
    
    Return: A graph
    """    
    or_edges0 = set(nx.get_edge_attributes(graphs[0],'nesting').keys())
    and_edges0 = set(graphs[0].edges()) - or_edges0    
    node0 = set(graphs[0].nodes())
  

    for g in graphs[1:]:
        or_edges1 = set(nx.get_edge_attributes(g,'nesting').keys())
        and_edges1 = set(g.edges()) - or_edges1    
        node1 = set(g.nodes())
        
        or_edges0 = or_edges0.union(or_edges1)
        and_edges0 = and_edges0.union(and_edges1)
        node0 = node0.union(node1)            
    
    G = nx.Graph()
    G.add_nodes_from(list(node0))
    G.add_edges_from(list(and_edges0),{'label':""})
    G.add_edges_from(list(or_edges0),{'label':"",'nesting':True})
    
    return G    