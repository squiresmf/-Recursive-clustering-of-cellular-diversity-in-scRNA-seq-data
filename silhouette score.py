import numpy as np
import pandas as pd
import time
from sklearn.metrics import silhouette_samples
import random
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import psutil


# Set the cluster sample size as needed
cluster_sample_size = None

# Function to determine optimal number of workers
def determine_optimal_workers():
    cpu_count = os.cpu_count()
    if cpu_count is None:
        return 4
    load = psutil.cpu_percent(interval=1)
    if load > 75:
        return max(1, cpu_count // 2)
    return cpu_count



# if os.name == 'nt':  # Windows
#     # Correct the format of the Windows drive letter
#     current_dir = os.path.join('Y:\\', 'qiu-lab', 'Michael Recursive Clustering')
# else:
#     current_dir = os.path.join('/qiu-lab', 'Michael Recursive Clustering')
current_dir = os.getcwd()

# Global dictionary to hold preloaded data (to be shared with workers)
preloaded_data = {}


# Initializer function for each worker process. Sets the global `preloaded_data` dictionary to the shared data.
def init_worker(shared_data):
    global preloaded_data
    preloaded_data = shared_data


def silhouette_score_func(X, labels, cluster_sample_size=None):
    start_time = time.time()

    unique_labels = np.unique(labels)

    if cluster_sample_size is not None:
        sampled_indices = []
        for label in unique_labels:
            cluster_indices = np.where(labels == label)[0]
            cluster_size = len(cluster_indices)

            if cluster_size <= cluster_sample_size:
                # Take all points if cluster size is small
                sampled_indices.extend(cluster_indices.tolist())
            else:
                # Randomly sample without replacement
                sampled = np.random.choice(cluster_indices, size=cluster_sample_size, replace=False)
                sampled_indices.extend(sampled.tolist())

        # Extract sampled data and labels
        X = X[sampled_indices]
        labels = labels[sampled_indices]

    # Compute silhouette scores for data
    silhouette_values = silhouette_samples(X, labels, metric='euclidean')

    # Compute average silhouette score per cluster
    cluster_scores = {}
    for label in unique_labels:
        cluster_silhouette = silhouette_values[labels == label]
        cluster_avg = np.mean(cluster_silhouette)
        cluster_scores[label] = cluster_avg

    end_time = time.time()
    duration = end_time - start_time
    print(f"Silhouette computation time: {duration:.2f} seconds")
    return cluster_scores


def process_reference_resolution(args):
    """
    Worker function to process a single (reference, resolution) task.
    Utilizes data_matrix from the global `preloaded_data`.
    """
    reference, c_resolution = args
    clusterings = ['recursive', 'seurat_equivalent']
    features_types = ['HVGs']
    leaf_clusters_scores = []

    ref_name = f'human_{reference}_integrated_'

    # Define paths
    default_dir = os.path.join(current_dir, 'Default', reference)

    # Access preloaded data from the global dictionary
    data_matrix = preloaded_data[reference]['data_matrix']

    for clustering in clusterings:
        method = 'Recursive' if clustering == 'recursive' else 'Single-Pass'

        cluster_assignments_path = os.path.join(
            default_dir,
            f'{clustering}_cluster_assignments_{ref_name}{c_resolution}_0.csv'
        )
        # Read cluster assignments
        cell_cluster_assignments = pd.read_csv(cluster_assignments_path)['cluster_assignment'].values

        for features_type in features_types:
            metric = f'Silhouette Score {features_type}'
            print(f"{reference} {c_resolution} {method} {features_type}")

            X = data_matrix

            # Compute silhouette scores
            cluster_scores = silhouette_score_func(
                X,
                cell_cluster_assignments,
                cluster_sample_size=cluster_sample_size
            )

            new_row = {
                'Reference': reference.capitalize() if reference != "PBMC" else reference,
                'Method': method,
                'Resolution': c_resolution,
                'Metric': metric,
                'Value': np.mean(list(cluster_scores.values())),
                'Number of Clusters': len(np.unique(cell_cluster_assignments))
            }

            leaf_clusters_scores.append(new_row)

    return leaf_clusters_scores


def main():
    # Start timing the main function
    main_start_time = time.time()

    leaf_clusters_scores_list = []
    references = ["PBMC", "adipose", "tonsil", "fetus"]
    res_range = np.arange(0.015, 0.155, 0.005)
    res_range = np.round(res_range, 4)

    # Preload all data matrices
    shared_data = {}
    for reference in references:
        print(f"Loading data for reference: {reference}")

        # Define paths
        default_dir = os.path.join(current_dir, 'Default', reference)
        matrix_path = os.path.join(default_dir, f"matrix {reference}.feather")

        # Load data
        data_matrix = pd.read_feather(matrix_path).values

        # Store in shared_data dictionary
        shared_data[reference] = {
            'data_matrix': data_matrix,
        }

    # Prepare arguments for parallel processing
    tasks = []
    for reference in references:
        for c_resolution in res_range:
            tasks.append((reference, c_resolution))

    # Determine the number of workers automatically
    try:
        num_workers = min(determine_optimal_workers(), 12)
    except Exception as e:
        print(f"Error determining optimal workers: {e}. Using default.")
        num_workers = os.cpu_count() or 4

    print(f"Using {num_workers} worker(s) for parallel processing.")

    # Use ProcessPoolExecutor for parallelism with initializer
    with ProcessPoolExecutor(max_workers=num_workers, initializer=init_worker, initargs=(shared_data,)) as executor:
        future_to_task = {executor.submit(process_reference_resolution, task): task for task in tasks}
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                result = future.result()
                leaf_clusters_scores_list.extend(result)
            except Exception as exc:
                print(f"Task {task} generated an exception: {exc}")

    leaf_clusters_scores_df = pd.DataFrame(leaf_clusters_scores_list)

    write_dir = os.path.join(current_dir, 'Default')
    title = 'leaf_clusters_scores.csv' if cluster_sample_size is None else f'leaf_clusters_scores_sample_{cluster_sample_size}.csv'
    write_path = os.path.join(write_dir, title)
    print(f"Saving results to {write_path}")
    leaf_clusters_scores_df.to_csv(write_path, index=False)

    main_end_time = time.time()
    total_duration = main_end_time - main_start_time
    print(f"Total execution time: {total_duration:.2f} seconds")


if __name__ == "__main__":
    main()
