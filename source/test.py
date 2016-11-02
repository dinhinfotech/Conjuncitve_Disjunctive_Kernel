import graph_util as gu
from graph_processing import Graph
from validation import Validation


adjacency_matrix_path = "/media/dinh/DATA/BIOINFOGEN/f3pc_data/test/adjacency_matrix_test/BioGPS.txt"
list_training_genes_path = "/media/dinh/DATA/BIOINFOGEN/f3pc_data/test/training_genes/0"
list_training_labels_path = "/media/dinh/DATA/BIOINFOGEN/f3pc_data/test/training_labels/0"
list_all_genes_path= "/media/dinh/DATA/BIOINFOGEN/f3pc_data/test/all_genes"
list_deg_thresholds = [10, 20]
list_cli_thresholds = [4,5] 
n_clusters = [20, 25]
list_r = [1,2]
list_d = [2,3]
list_c = [10**(-5),10**(-4),10**(-3),10**(-2),10**(-1),10**(0),10**(1),10**(2),10**(3),10**(4)]

Val = Validation(adjacency_matrix_path=adjacency_matrix_path,
                 list_training_genes_path=list_training_genes_path,
                 list_training_labels_path=list_training_labels_path,
                 list_all_genes_path=list_all_genes_path,
                 list_deg_thresholds=list_deg_thresholds,
                 list_cli_thresholds=list_cli_thresholds,
                 n_clusters=n_clusters,
                 list_r=list_r,
                 list_d=list_d,
                 list_c=list_c
                 )

optimal_paras = ((20, 4, 25), (1, 3), 100)
auc = Val.evaluation(optimal_paras)                
print "AUC: ", auc
