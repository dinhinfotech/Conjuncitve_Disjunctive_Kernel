import util
import numpy as np
from sklearn.cluster import MiniBatchKMeans

def get_dict_goterm_index(goterm_ancestor_file_path):
    listlines = util.load_list_from_textfile(goterm_ancestor_file_path)    
    dict_goterm_index = {}   
    list_all_goterms = []
    idx = 0
    for line in listlines:
        strs = line.split()
        for st in strs:
            if st not in dict_goterm_index:
                dict_goterm_index[st] = idx
                list_all_goterms.append(st)
                idx+=1
    return [dict_goterm_index,list_all_goterms]

def get_dict_gene_goterms(gene_goterm_file_path):
    listlines = util.load_list_from_textfile(gene_goterm_file_path)    
    dict_gene_goterms = {}    
    for line in listlines:
        strs = line.split()
        dict_gene_goterms[strs[0]] = strs[1:]
    return dict_gene_goterms

def get_vectorizing_gene_goterms(list_all_genes_new_file_path, list_adjacency_matrices,
                                 goterm_ancestor_file_path, gene_goterm_file_path):
    # Using adjacency matrices in order to replace a gene does not appear on GO ontology by 
    # one of its neighbors that shows in the Go ontology
    [dict_goterm_index,list_all_goterms] = get_dict_goterm_index(goterm_ancestor_file_path) 
    dict_gene_goterms = get_dict_gene_goterms(gene_goterm_file_path) 
    dict_goterm_goterms = get_dict_gene_goterms(goterm_ancestor_file_path)
    list_all_genes_new = util.load_list_from_textfile(list_all_genes_new_file_path)
    N = list_adjacency_matrices[0].shape[0]
    list_all_genes = []
    for idx1, gene1 in enumerate(list_all_genes_new):
        if gene1 in dict_gene_goterms:
            list_all_genes.append(gene1)
        else:
            Flag = False
            for idx2, gene2 in enumerate(list_all_genes_new):
                if list_adjacency_matrices[0][idx1,idx2]!=0 and gene2 in dict_gene_goterms:
                    list_all_genes.append(gene2)
                    Flag = True
                    break
            if Flag==False:
                for idx2, gene2 in enumerate(list_all_genes_new):
                    if list_adjacency_matrices[1][idx1,idx2]!=0 and gene2 in dict_gene_goterms:
                        list_all_genes.append(gene2)
                        Flag = True
                        break
            if Flag==False:
                for idx2, gene2 in enumerate(list_all_genes_new):
                    if list_adjacency_matrices[2][idx1,idx2]!=0 and gene2 in dict_gene_goterms:
                        list_all_genes.append(gene2)
                        Flag = True
                        break                             
    M = len(list_all_goterms)
    G = np.zeros((N,M))
    
    for idx, gene in enumerate(list_all_genes):
        list_leaf_goterms = dict_gene_goterms[gene]
        for leaf_goterm in list_leaf_goterms:
            leaf_idx = dict_goterm_index[leaf_goterm]
            G[idx,leaf_idx] = 1
            list_goterms = dict_goterm_goterms[leaf_goterm]
            for goterm in list_goterms:
                goterm_idx = dict_goterm_index[goterm]
                G[idx,goterm_idx] = 1
    return G

def node_clustering(n_clusters):
    
    list_all_genes_new_file_path = "/lustre/guido.zampieri/Dinh/AOLGK/data/all_genes"
    list_adjacency_folder = "/lustre/guido.zampieri/Dinh/AOLGK/data/adjacency_matrix/"
    goterm_ancestor_file_path = "/lustre/guido.zampieri/GO_analysis/GObioProc/GObioProc_ancestor.txt"
    gene_goterm_file_path = "/lustre/guido.zampieri/GO_analysis/GObioProc/GObioProc_newList.txt"
    
    list_adjacency_matrices = util.load_matrices(list_adjacency_folder)
    
    G = get_vectorizing_gene_goterms(list_all_genes_new_file_path, list_adjacency_matrices,
                                 goterm_ancestor_file_path, gene_goterm_file_path)
                                     
    model = MiniBatchKMeans(n_clusters=n_clusters,init='k-means++',
                                       max_iter=10,n_init=10,
                                          random_state=n_clusters)
    model.fit(G)
    label_list = [str(element) for element in  model.predict(G)]
    return label_list    