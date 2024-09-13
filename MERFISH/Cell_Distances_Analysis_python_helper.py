import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm

def calculate_median_nearest_neighbor_distance(df):
    """
    Calculate the median nearest neighbor distance for a set of points.
    
    Parameters:
    df (pd.DataFrame): DataFrame containing 'center_x' and 'center_y' columns.
    
    Returns:
    float: The median nearest neighbor distance.
    """
    # Create a KDTree from the points
    tree = cKDTree(df[['center_x', 'center_y']].values)
    # Query the two nearest neighbors since the closest one is the point itself
    distances, _ = tree.query(df[['center_x', 'center_y']].values, k=2)
    # Take the second closest (which is the first nearest neighbor)
    nearest_neighbor_distances = distances[:, 1]
    # Compute the median of the nearest neighbor distances
    median_distance = np.median(nearest_neighbor_distances)
    return median_distance

def get_nearest_neighbors(df):
    """
    Compute the nearest neighbors for each cell type in the dataset and normalize distances.
    
    Parameters:
    df (pd.DataFrame): DataFrame containing 'center_x', 'center_y', 'cell_type', and 'cell_ID' columns.
    
    Returns:
    pd.DataFrame: DataFrame with additional columns for nearest neighbor distances and IDs for each cell type.
    """
    results = df.copy()
    cell_types = df['cell_type'].unique()
    # Calculate the median distance using all cells
    median_distance = calculate_median_nearest_neighbor_distance(df)

    # Precompute the KDTree for each cell type and subset dataframes
    trees = {ctype: cKDTree(df[df['cell_type'] == ctype][['center_x', 'center_y']])
             for ctype in cell_types}
    subsets = {ctype: df[df['cell_type'] == ctype].reset_index(drop=True) 
               for ctype in cell_types}

    # Initialize columns for distances and IDs
    for ctype in cell_types:
        results[f'{ctype}_dist'] = np.nan
        results[f'{ctype}_id'] = None

    # Iterate over each cell with a progress bar
    for index in tqdm(df.index, desc='Computing nearest neighbors'):
        row = df.loc[index]
        x, y = row['center_x'], row['center_y']

        # Search nearest neighbor for each cell type
        for ctype in cell_types:
            subset = subsets[ctype]
            tree = trees[ctype]

            if len(subset) == 1 and subset['cell_ID'].values[0] == row['cell_ID']:
                # Set the normalized distance to zero
                results.at[index, f'{ctype}_dist'] = 0
                results.at[index, f'{ctype}_id'] = row['cell_ID']
                continue

            if len(subset) > 1:
                dists, idxs = tree.query([x, y], k=2)
                # Select the second closest if the closest is the point itself
                if subset.iloc[idxs[0]]['cell_ID'] == row['cell_ID']:
                    nearest_dist = dists[1]
                    nearest_idx = idxs[1]
                else:
                    nearest_dist = dists[0]
                    nearest_idx = idxs[0]

                nearest_cell_id = subset.iloc[nearest_idx]['cell_ID']

                # Normalize the distance by the median_distance
                normalized_dist = nearest_dist / median_distance

                results.at[index, f'{ctype}_dist'] = normalized_dist
                results.at[index, f'{ctype}_id'] = nearest_cell_id

    return results

def calculate_distance_statistics(df):
    """
    Calculate distance statistics, including mean and median nearest neighbor distances.
    
    Parameters:
    df (pd.DataFrame): DataFrame containing 'center_x' and 'center_y' columns.
    
    Returns:
    dict: Dictionary containing the distance DataFrame and calculated statistics.
    """
    # Create a KDTree from the points
    tree = cKDTree(df[['center_x', 'center_y']].values)
    # Query the two nearest neighbors since the closest one is the point itself
    distances, _ = tree.query(df[['center_x', 'center_y']].values, k=2)
    # Take the second closest (which is the first nearest neighbor)
    nearest_neighbor_distances = distances[:, 1]
    # Compute the median and mean of the nearest neighbor distances
    median_distance = np.median(nearest_neighbor_distances)
    mean_distance = nearest_neighbor_distances.mean()
    nearest_neighbor_distances = nearest_neighbor_distances.tolist()  # Convert to list

    # Return a dictionary with the DataFrame and statistics
    results = {
        'distances_df': nearest_neighbor_distances,
        'mean_distance': mean_distance,
        'median_distance': median_distance
    }

    return results
