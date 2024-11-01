
import numpy as np
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, \
    homogeneity_score, completeness_score, v_measure_score
from sklearn.cluster import KMeans

from balanced_clustering import balanced_adjusted_rand_index, \
    balanced_adjusted_mutual_info, balanced_completeness, \
    balanced_homogeneity, balanced_v_measure

import pandas as pd
import os

import os

# if os.name == 'nt':  # Windows
#     # Correct the format of the Windows drive letter
#     current_dir = os.path.join('Y:\\', 'qiu-lab', 'Michael Recursive Clustering')
# else:
#     current_dir = os.path.join('/qiu-lab', 'Michael Recursive Clustering')
current_dir = os.getcwd()

def set_dir_read():
    setwd_path = os.path.join(
        current_dir,
        'Default',
        reference
    )
    os.chdir(setwd_path)

def set_dir_write():
    setwd_path = os.path.join(current_dir, 'Default')
    os.chdir(setwd_path)


results_df = pd.DataFrame({
    'Reference': [],
    'Method': [],
    'Annotation_Level': [],
    'Resolution': [],
    'Metric': [],
    'Value': [],
    'Number_of_Clusters': []
})

clusterings = ['recursive', 'seurat_equivalent']

for reference in ["PBMC", "adipose", "tonsil", "fetus"]:
    if reference in ["PBMC", "adipose", "tonsil"]:
        res_range = np.concatenate([np.arange(0.005, 0.05, 0.002), np.arange(0.05, 0.15, 0.005)])
    else:
        res_range = np.concatenate([np.arange(0.0025, 0.01, 0.0002), np.arange(0.015, 0.15, 0.005)])
    res_range = np.round(res_range, 4)

    ref_name = f'human_{reference}_integrated_'

    set_dir_read()

    for annotation_level in [1, 2]:
        csv_file = f'{reference} celltype.l{annotation_level}.csv'
        if os.path.exists(csv_file):
            # Read the CSV file for cell types
            cell_celltype_df = pd.read_csv(csv_file)
            cell_celltype = np.array(cell_celltype_df['celltype'].tolist())
        else:
            raise FileNotFoundError(f"CSV file {csv_file} not found")

        for c_resolution in res_range:
            print(f'{reference} {c_resolution} {annotation_level}')
            for clustering in clusterings:
                csv_file = f'{clustering}_cluster_assignments_{ref_name}{c_resolution}_0.csv'
                if os.path.exists(csv_file):
                    # Read the CSV file for cluster assignments
                    cell_cluster_assignments_df = pd.read_csv(csv_file)
                    cell_cluster_assignments = np.array(cell_cluster_assignments_df['cluster_assignment'].tolist())
                    num_clusters = len(set(cell_cluster_assignments))
                else:
                    raise FileNotFoundError(f"CSV file {csv_file} not found")

                # Compute metrics
                balanced_ari = balanced_adjusted_rand_index(labels_true=cell_celltype, labels_pred=cell_cluster_assignments)
                balanced_ami = balanced_adjusted_mutual_info(labels_true=cell_celltype, labels_pred=cell_cluster_assignments)
                balanced_compl = balanced_completeness(labels_true=cell_celltype, labels_pred=cell_cluster_assignments)
                balanced_homog = balanced_homogeneity(labels_true=cell_celltype, labels_pred=cell_cluster_assignments)
                balanced_vmeas = balanced_v_measure(labels_true=cell_celltype, labels_pred=cell_cluster_assignments)

                metrics = {
                    'Balanced ARI': balanced_ari,
                    'Balanced AMI': balanced_ami,
                    'Balanced Completeness': balanced_compl,
                    'Balanced Homogeneity': balanced_homog,
                    'Balanced V-Measure': balanced_vmeas,
                }

                method = 'Recursive' if clustering == 'recursive' else 'Single-Pass'

                # Add metrics to the results DataFrame
                for metric_type, value in metrics.items():
                    new_row = {
                        'Reference': reference if reference == "PBMC" else reference.capitalize(),
                        'Method': method,
                        'Annotation_Level': f'Annotation Level {annotation_level}',
                        'Resolution': c_resolution,
                        'Metric': metric_type,
                        'Value': value,
                        'Number_of_Clusters': num_clusters
                    }
                    results_df = pd.concat([results_df, pd.DataFrame([new_row])], ignore_index=True)

results_df.columns = results_df.columns.str.replace('_', ' ')

set_dir_write()

results_df.to_csv('balanced_metrics_df.csv', index=False)

print('done')