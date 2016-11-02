from graph_processing import Graph
import util
import graph_util as gu
from sklearn import cross_validation
from collections import defaultdict
from sklearn import svm
from sklearn import metrics
import networkx as nx

class Validation:
    def __init__(self,
                 adjacency_matrix_path=None,
                 list_training_genes_path=None,
                 list_training_labels_path=None,
                 list_all_genes_path=None,
                 list_deg_thresholds=None,
                 list_cli_thresholds=None,
                 n_clusters=None,
                 list_r=None,
                 list_d=None, 
                 list_c=None,
                 n_folds=3
                 ):
        self.adjacency_matrix_path = adjacency_matrix_path
        self.list_training_genes_path = list_training_genes_path
        self.list_training_labels_path = list_training_labels_path
        self.list_all_genes_path = list_all_genes_path
        self.list_deg_thresholds = list_deg_thresholds
        self.list_cli_thresholds = list_cli_thresholds
        self.n_clusters = n_clusters
        self.list_r = list_r
        self.list_d = list_d
        self.list_c = list_c
        self.n_folds = n_folds
        self.paras = None
        
    def para_selection(self):
        
        G = gu.load_graph(self.adjacency_matrix_path)
        training_genes = util.load_list_from_file(self.list_training_genes_path)
        training_labels = [int(e) for e in util.load_list_from_file(self.list_training_labels_path)]
        list_all_genes = util.load_list_from_file(self.list_all_genes_path)
        dict_gene_idx = {}
        for idx, gene in enumerate(list_all_genes):
            dict_gene_idx[gene]=idx
        
        max_node_ID = max(G.nodes())+1
        n_training_genes = len(training_genes)
        
        dec_paras = []
        for idx1 in self.list_deg_thresholds:
            for idx2 in self.list_cli_thresholds:
                for idx3 in self.n_clusters:
                    dec_paras.append((idx1,idx2,idx3))
                
        vec_paras = []
        for idx1 in self.list_r:
            for idx2 in self.list_d:
                vec_paras.append((idx1,idx2))
        
        
        dict_para_auc = defaultdict(lambda: 0)
        
        kf = cross_validation.KFold(n_training_genes, n_folds = self.n_folds)
        
        for train_index, test_index in kf:
            validation_genes = [training_genes[idx] for idx in train_index]
            validation_idx = [dict_gene_idx[gene] for gene in validation_genes]
            validation_labels = [training_labels[idx] for idx in train_index]
            test_genes = [training_genes[idx] for idx in test_index]
            test_labels = [training_labels[idx] for idx in test_index]
            
            list_unknown_genes = []
            list_unknown_genes.extend(test_genes)
            for gene in list_all_genes:
                if gene not in training_genes:
                    list_unknown_genes.append(gene)                    
            list_unknown_idx = [dict_gene_idx[gene] for gene in list_unknown_genes]                                
        
        
            for dec_para in dec_paras:
                GP = Graph(deg_threshold=dec_para[0], cli_threshold=dec_para[1],
                           max_node_ID = max_node_ID)
                [G_dec, dict_newnode_clinodes] = GP.decompose(G)
                
                dict_node_attvalues = {}
                for n in G_dec.nodes_iter():
                        dict_node_attvalues[n]= "A"
        
                nx.set_node_attributes(G_dec,'label',dict_node_attvalues)    
                #GP.labeling(G_dec, dec_para[2], dict_newnode_clinodes)
                
                for vec_para in vec_paras:
                    M = GP.vectorize(G_dec)
                    
                    M_val = M[validation_idx,:]
                    M_unknown = M[list_unknown_idx,:]
                    
                    for c in self.list_c:                        
                        clf = svm.LinearSVC(C=c)
                        clf.fit(M_val, validation_labels)
                        scores = clf.decision_function(M_unknown)
                        
                        qscores = []
                        
                        for s in scores[:len(test_genes)]:
                            qscore = float(sum([int(s >= value) for value in scores]))/len(scores)
                            qscores.append(qscore)
                        fpr, tpr, thresholds = metrics.roc_curve(test_labels, qscores, pos_label= 1)
                        auc = metrics.auc(fpr, tpr)
                        dict_para_auc[(dec_para,vec_para,c)]+=auc
                        
                        print auc
                    print "--------------"
        optimal_paras = max(dict_para_auc.iterkeys(), key=lambda k: dict_para_auc[k])
        self.paras = optimal_paras
        return optimal_paras

    def evaluation(self, optimal_paras=None):
        
        G = gu.load_graph(self.adjacency_matrix_path)
        training_genes = util.load_list_from_file(self.list_training_genes_path)
        training_labels = [int(e) for e in util.load_list_from_file(self.list_training_labels_path)]
        list_all_genes = util.load_list_from_file(self.list_all_genes_path)
        dict_gene_idx = {}
        for idx, gene in enumerate(list_all_genes):
            dict_gene_idx[gene]=idx
        
        max_node_ID = max(G.nodes())+1              
        
        GP = Graph(deg_threshold= optimal_paras[0][0],
                   cli_threshold= optimal_paras[0][1],
                   max_node_ID = max_node_ID,
                   r = optimal_paras[1][0],
                   d = optimal_paras[1][1])
        
        [G_dec, dict_newnode_clinodes] = GP.decompose(G)
        
        dict_node_attvalues = {}
        for n in G_dec.nodes_iter():
                dict_node_attvalues[n]= "A"

        nx.set_node_attributes(G_dec,'label',dict_node_attvalues)          
        
        M = GP.vectorize(G_dec)
        qscores = []
        
        for idx, gene in enumerate(training_genes):
            validation_genes = training_genes[:]
            validation_labels = training_labels[:]
            del validation_genes[idx]
            del validation_labels[idx]
            validation_idx = [dict_gene_idx[ge] for ge in validation_genes]
            M_val = M[validation_idx,:]            
            
            unknown_genes = [gene]
            for ge in list_all_genes:
                if ge not in training_genes:
                    unknown_genes.append(ge)
            unknown_idx = [dict_gene_idx[ge] for ge in unknown_genes]
            
            M_unknown = M[unknown_idx,:]
            
            
            clf = svm.LinearSVC(C=optimal_paras[2])
            clf.fit(M_val, validation_labels)
            scores = clf.decision_function(M_unknown)
            
            qscore = float(sum([int(scores[0] >= val) for val in scores]))/len(scores)
            
            qscores.append(qscore)
            
            print "Done gene ", idx
            
        fpr, tpr, thresholds = metrics.roc_curve(training_labels, qscores, pos_label= 1)
        auc = metrics.auc(fpr, tpr)
        
        return auc