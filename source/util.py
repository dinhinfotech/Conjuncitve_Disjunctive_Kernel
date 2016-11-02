# -*- coding: utf-8 -*-
"""
Util file includes utility functions
"""
from os import listdir
from os.path import isfile, join
import numpy as np
from cvxopt import matrix

def list_files_in_folder(folder_path):    
    """
    Return: A list of the file names in the folder
    """
          
    list = listdir(folder_path)
    onlyfiles = [ f for f in list  if isfile(join(folder_path,f)) ]
    return onlyfiles 

def load_list_from_file(file_path):
    """
    Return: A list saved in a file
    """
    
    f = open(file_path,'r')
    listlines = [line.rstrip() for line in f.readlines()]
    f.close()
    return listlines

def load_matrices(folder_path):
    """
    Return: A list of matrices saved in the folder
    """
    
    file_names = list_files_in_folder(folder_path)   
    matrices = []
    
    for file_name in file_names:
        matrices.append(np.loadtxt(folder_path + file_name))
        
    return matrices  

def filter_disease_genes(folder_in, folder_out, genes_path):
    """
    Parameters:
    - folder_in: Folder containing disease gene lists
    - folder_out: Saveing new disease gene lists
    - genes_path: list of genes corresponding to an adjacency matrix
    """
    
    names = list_files_in_folder(folder_in)
    genes = load_list_from_file(genes_path)
    
    for name in names:
        disease_genes = load_list_from_file(folder_in + name)
        
        new_disease_genes = []
        for gene in disease_genes:
            if gene in genes:
                new_disease_genes.append(gene)
        
        f = open(folder_out + name,'w')
        f.writelines([line + "\n" for line in new_disease_genes])
        f.close()

def extract_submatrix(row_indices, col_indices, A):
    """ Extracting a submatrix from  matrix A
    
    Parameter:
    row_indices: row index list that we want to extract
    col_indices: Column index list that we want to extract
    A: Matrix
    
    Return:
    submatrix of A
    """

    len_row = len(row_indices)
    len_col = len(col_indices)
    
    M = matrix(0.0,(len_row,len_col))
    for order1, idx_row in enumerate(row_indices):
        for order2, idx_col in enumerate(col_indices):
            M[order1,order2] = A[idx_row,idx_col]
    
    return M

def extract_test_adjacency_matrix(positive_genes_path, adjacency_path_in, adjacency_path_out, 
                                  genes_path):
    """
    Parameters:
    - folder_in: Folder containing disease gene lists
    - adjacency_path_in: full adjacency matrix
    - adjacency_path_out: file to save extracted adjacency matrix
    - genes_path: list of genes corresponding to an adjacency matrix
    """        
                             
    A = np.loadtxt(adjacency_path_in)        
    genes = load_list_from_file(genes_path)
    
    all_disease_genes = load_list_from_file(positive_genes_path)

    dict_gene_idx = {}
    for idx, gene in enumerate(genes):
        dict_gene_idx[gene] = idx
    
    list_idx = []
    for gene in all_disease_genes:
        list_idx.append(dict_gene_idx[gene])    
    test_matrix = extract_submatrix(list_idx, list_idx, A)    
    np.savetxt(adjacency_path_out,test_matrix)    