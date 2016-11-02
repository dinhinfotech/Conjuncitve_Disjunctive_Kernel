# -*- coding: utf-8 -*-
import networkx as nx
from eden.graph import Vectorizer
from sklearn.metrics import pairwise
import ontology_processing as op

class Graph():
    def __init__(self,
                 deg_threshold=None,
                 cli_threshold=None,
                 max_node_ID=None,
                 r=1,
                 d=2,
                 nbits=20,
                 discrete=True,
                 n_jobs=1):
        self.deg_threshold = deg_threshold
        self.cli_threshold = cli_threshold
        self.max_node_ID = max_node_ID
        self.r = r
        self.d = d
        self.nbits = nbits
        self.discrete = discrete
        self.n_jobs = n_jobs  
        self.is_decomposed = False
    
    def kernel_matrix(self,G):
        
        if self.is_decomposed == False:
            [M, dict_newnode_clinodes] = self.vectorize(G)
            K = pairwise.linear_kernel(M,M)
            return [K, dict_newnode_clinodes]
        else:
            M = self.vectorize(G)
            K = pairwise.linear_kernel(M,M)
            return K

    def vectorize(self, G):
        # Initialize Vectorizer here and make new version of vectorize
        if self.is_decomposed == False:
            [G, dict_newnode_clinodes]  = self.decompose(G)
            vec = Vectorizer(nbits=self.nbits, 
                 discrete=self.discrete, 
                 n_jobs=self.n_jobs, 
                 r=self.r, 
                 d=self.d)
                             
            M = vec.vertex_transform([G])[0] 
            return [M, dict_newnode_clinodes]
        else:
                
            vec = Vectorizer(nbits=self.nbits, 
                             discrete=self.discrete, 
                             n_jobs=self.n_jobs, 
                             r=self.r, 
                             d=self.d)
                             
            M = vec.vertex_transform([G])[0]                

            return M

    def decompose(self, graph):
        [G, dict_newnode_clinodes] = self.kcore_clique_decompose(graph,
                                                                self.deg_threshold,
                                                                self.cli_threshold,
                                                                self.max_node_ID)
        self.is_decomposed = True
        return [G, dict_newnode_clinodes] 
    
    def labeling(self, G, n_clusters, dict_newnode_clinodes): 
        
        list_labels = op.node_clustering(n_clusters)
        dict_node_attvalues = {}
        for n in G.nodes_iter():
            if n not in dict_newnode_clinodes:
                dict_node_attvalues[n]= list_labels[n]
        
        for n, cli_nodes in dict_newnode_clinodes.iteritems():
            labels = [list_labels[v] for v in cli_nodes]
            label = max(set(labels), key=labels.count)
            dict_node_attvalues[n] = label
        nx.set_node_attributes(G,'label',dict_node_attvalues)         
   
    def clique_decompose(self, 
                         G, 
                         max_node_ID, 
                         dict_newnode_clinodes, 
                         remove_edges, 
                         cli_threshold):
        
        list_cliques = [cli for cli in nx.find_cliques(G) if len(cli) >=cli_threshold]

        new_node_IDs = []
        nesting_edges = []
        new_edges = []    
        
        for cli in list_cliques:        
            # Updage dict_nodeID_clinodes
            dict_newnode_clinodes[max_node_ID[0]] = cli
            # Add a new node
            new_node_IDs.append(max_node_ID[0])        
            # Add nesting edges connecting nodes in clique with the new node and edges
            # connecting all neighbors of clique's nodes to the new node
            for n in cli:
                nesting_edges.append((n,max_node_ID[0]))
                for v in nx.neighbors(G,n):
                    if v not in cli:
                        new_edges.append((v,max_node_ID[0]))
                        remove_edges.append((n,v))
            max_node_ID[0] = max_node_ID[0]+1
                    
                
        G.add_nodes_from(new_node_IDs)
        G.add_edges_from(nesting_edges, label="", nesting ='True')
        G.add_edges_from(new_edges,label="")        
                     
    def kcore_clique_decompose(self, 
                               G, 
                               deg_threshold, 
                               cli_threshold, 
                               max_node_ID):
        max_node_ID = [max_node_ID]
        dict_newnode_clinodes = {}
        remove_edges = []
        
        high_degree_nodes = []
        low_degree_nodes = []
        m = min(G.degree().values())
        t = max(deg_threshold,m)
        
        for n in G.nodes():
            if len(G.neighbors(n)) > t:
                high_degree_nodes.append(n)
            else:
                low_degree_nodes.append(n)
        G_high_degree0 = G.subgraph(high_degree_nodes)
        n_nodes = len(high_degree_nodes)
        
        G_low_degree0 = G.subgraph(low_degree_nodes)
        self.clique_decompose(G_low_degree0, max_node_ID, dict_newnode_clinodes,
                              remove_edges, cli_threshold)        
        G_union = G_low_degree0.copy()
    
        while n_nodes > 0:        
            high_degree_nodes = []
            low_degree_nodes = []
            m = min(G_high_degree0.degree().values())
            t = max(deg_threshold,m)
            
            for n in G_high_degree0.nodes():
                if len(G_high_degree0.neighbors(n)) > t:      
                    high_degree_nodes.append(n)
                else:
                    low_degree_nodes.append(n)
            G_high_degree = G_high_degree0.subgraph(high_degree_nodes)
            n_nodes = len(high_degree_nodes)
            
            G_low_degree = G_high_degree0.subgraph(low_degree_nodes)
            self.clique_decompose(G_low_degree, max_node_ID, dict_newnode_clinodes, 
                                  remove_edges, cli_threshold)
            
            G_union = nx.union(G_union,G_low_degree)
            
            G_high_degree0 = G_high_degree
            
        edges = G.edges()
        edges_union = G_union.edges()
        edges_nesting = list((set(edges) - set(edges_union))-set(remove_edges))
        G_union.add_edges_from(edges_nesting, label="", nesting ='True')
        # Remove edges inside cliques
        clique_edges = []
        for nodes in dict_newnode_clinodes.values():
            for idx1 in range(len(nodes)-1):
                for idx2 in range(idx1+1,len(nodes)):
                    if idx1 != idx2:
                        clique_edges.append((nodes[idx1],nodes[idx2]))
        G_union.remove_edges_from(clique_edges)
        G_union.remove_edges_from(remove_edges)
        
        return [G_union, dict_newnode_clinodes]           
    


        
    