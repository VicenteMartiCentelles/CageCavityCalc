import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from sklearn.neighbors import KernelDensity
from sklearn.neighbors import NearestNeighbors
from kneed import KneeLocator
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA

def dbscan_clustering_windows(data_points, cage_atoms, atom_names):

    # Automated Parameter Determination for DBSCAN Clustering
    
    data_points = np.array(data_points)

    # Determine min_samples (k) as standard heuristic for 3D data: 2 * dimension
    min_samples_value = 2 * data_points.shape[1] 
    print(f"Set min_samples automatically to: {min_samples_value}")

    # Determine eps (epsilon) automatically using the k-distance method
    k = min_samples_value # k is min_samples for the k-distance graph
    k_neighbors = k + 1   # We look for the distance to the k-th nearest neighbor (k+1 index)

    # Compute the distance to the k-th nearest neighbor for every point
    neigh = NearestNeighbors(n_neighbors=k_neighbors)
    nbrs = neigh.fit(data_points)
    distances, indices = nbrs.kneighbors(data_points)

    # Sort the distances to the k-th neighbor
    k_distances = np.sort(distances[:, k], axis=0)

    # Use KneeLocator to find the "elbow" point (the optimal eps)
    # S=1.0 is the default sensitivity parameter. 'concave' and 'increasing' describe the curve shape.
    kn = KneeLocator(
        x=range(len(k_distances)), 
        y=k_distances, 
        S=1.0, 
        curve='convex', # The k-distance plot is convex, not concave.
        direction='increasing'
    )

    # The optimal epsilon is the value at the knee point
    eps_value = kn.elbow_y if kn.elbow_y is not None else 0.5 
    # Fallback to a small default value if no knee is found
    print(f"Optimal eps determined by KneeLocator: {eps_value:.4f}")

    if eps_value == 0 or np.isclose(eps_value, 0):
        # Avoid zero or near-zero eps which can crash or yield poor results
        print("\nWARNING: Auto-determined eps is too close to zero. Using fallback eps=0.5")
        eps_value = 0.5        


    # DBSCAN Clustering
    dbscan = DBSCAN(eps=eps_value, min_samples=min_samples_value)

    # Fit the model and predict the cluster for each point
    cluster_labels = dbscan.fit_predict(data_points)

    # Cluster labels: -1 indicates a Noise point, 0, 1, 2, ... indicate the cluster ID
    n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
    n_noise = list(cluster_labels).count(-1)

    print(f"DBSCAN Parameters: eps={eps_value}, min_samples={min_samples_value}")
    print(f"Number of estimated window clusters: {n_clusters}")
    print(f"Number of noise points: {n_noise}")
    #print("First 10 cluster labels:", cluster_labels[:10])





    # --- Ellipsoid Fitting and Plotting ---
    
    # Prepare the 3D plot
    #fig = plt.figure(figsize=(12, 10))
    #ax = fig.add_subplot(111, projection='3d')

    unique_labels = set(cluster_labels)
    colors = plt.cm.get_cmap('tab20', len(unique_labels)) # Use a distinct colormap
    all_fitted_ellipsoids = [] # To store data for ellipsoids to plot later

    print("\n--- Ellipsoid Fitting Results ---")
    print("-" * 35)

    for i, label in enumerate(unique_labels):
        if label == -1: # Noise points
            color_for_label = 'k' # Black for noise
            label_text = 'Noise'
            zorder_val = 1 # Plot noise behind clusters
        else:
            color_for_label = colors(label % colors.N) # Assign a color from the colormap
            label_text = f'Window {label}'
            zorder_val = 2 # Plot clusters in front

        cluster_points = data_points[cluster_labels == label]
        
        # Plot the cluster points
        ax.scatter(cluster_points[:, 0], cluster_points[:, 1], cluster_points[:, 2], 
                   color=color_for_label, marker='o', s=20, alpha=0.6, label=label_text, zorder=zorder_val)

        if label != -1 and len(cluster_points) >= 3: # Only fit ellipsoids to non-noise clusters with enough points
            # 1. Calculate the Center (Centroid)
            center = np.mean(cluster_points, axis=0)
            
            # 2. Fit PCA to find Axis Orientation and Lengths
            pca = PCA(n_components=3)
            pca.fit(cluster_points)
            
            N_sigma = 1
            std_devs = np.sqrt(pca.explained_variance_)
            axis_lengths = 2 * N_sigma * std_devs
            
            a_axis_len = axis_lengths[0]  # Longest axis
            b_axis_len = axis_lengths[1]  # Second longest axis
            c_axis_len = axis_lengths[2]  # Shortest axis (flatness)

            # The eigenvectors are the principal axes of the ellipsoid
            rotation_matrix = pca.components_.T # Columns are the eigenvectors

            print(f"Cluster {label} (Size: {len(cluster_points)}):")
            print(f"  Center (x, y, z): {center}")
            print(f"  Longest Axis (a): {a_axis_len:.4f}")
            print(f"  Second Axis (b): {b_axis_len:.4f}")
            print(f"  Flatness Axis (c): {c_axis_len:.4f}")
            print("-" * 35)
            
            # Store ellipsoid data for plotting
            all_fitted_ellipsoids.append({
                'center': center,
                'axes': axis_lengths,
                'rotation': rotation_matrix,
                'color': color_for_label
            })

    '''
    # --- Plotting the Ellipsoids ---
    for ellipsoid_data in all_fitted_ellipsoids:
        center = ellipsoid_data['center']
        axes = ellipsoid_data['axes']
        rotation = ellipsoid_data['rotation']
        color = ellipsoid_data['color']

        # Generate points on a unit sphere (50x50 resolution)
        u = np.linspace(0, 2 * np.pi, 50)
        v = np.linspace(0, np.pi, 50)
        x_sphere = np.outer(np.cos(u), np.sin(v))
        y_sphere = np.outer(np.sin(u), np.sin(v))
        z_sphere = np.outer(np.ones_like(u), np.cos(v))
        
        # Create an array of (x, y, z) coordinates for the sphere (flattened)
        points_sphere = np.vstack([x_sphere.flatten(), 
                                   y_sphere.flatten(), 
                                   z_sphere.flatten()])

        # 1. Scale by half the axis lengths (radii)
        # 2. Rotate by the PCA components (rotation matrix)
        # 3. Translate by the center
        
        # Ellipsoid axes are radii: axes/2. The scaling operation is a diagonal matrix S.
        # The transformation is: T = Rotation @ Scaling @ Sphere_Points + Center
        
        # Scale the sphere points by the radii (axes/2)
        # Note: numpy arrays allow element-wise multiplication with a 3xN array here.
        points_scaled = points_sphere * (axes / 2)[:, np.newaxis] 

        # Rotate the scaled points using the rotation matrix (PCA components)
        points_rotated = np.dot(rotation, points_scaled)

        # Translate the rotated points by the center and reshape
        points_final = points_rotated + center[:, np.newaxis]

        # Reshape the final points back into the 50x50 grid for plotting
        x_ellipsoid = points_final[0, :].reshape(x_sphere.shape)
        y_ellipsoid = points_final[1, :].reshape(y_sphere.shape)
        z_ellipsoid = points_final[2, :].reshape(z_sphere.shape)

        # Plot the ellipsoid surface
        ax.plot_surface(x_ellipsoid, y_ellipsoid, z_ellipsoid, 
                        color=color, alpha=0.15, rstride=4, cstride=4, 
                        linewidth=0.5, edgecolor=color, zorder=3)
        
        # Plot the center of the ellipsoid
        #ax.scatter(center[0], center[1], center[2], color=color, 
        #           marker='X', s=100, label=f'Center {label}', zorder=4)
    
    #Plot the atoms of the cage
    xyz = cage_atoms
    ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=50, color='k', label='Atoms') 

    ax.set_title(f'DBSCAN Clusters with Fitted Ellipsoids (eps={eps_value:.2f}, min_samples={min_samples_value})')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.legend(loc='lower left', bbox_to_anchor=(-0.1, 0.1)) # Adjust legend position for better visibility
    plt.tight_layout()
    plt.show()'''
    
    return all_fitted_ellipsoids



# Constants for PDB formatting
PDB_FORMAT_HETATM = (
    "HETATM{:5d} {:<4s} {:3s} {:1s}{:4d}    "
    "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:<2s}  \n"
)
# ANISOU record format: serial, name, altLoc, resName, chainID, resSeq, U11..U23 (scaled by 10000), element, charge
PDB_FORMAT_ANISOU = (
    "ANISOU{:5d} {:<4s}{:1s}{:3s} {:1s}{:4d} " 
    "{:7d}{:7d}{:7d}{:7d}{:7d}{:7d}       {:<2s}  \n" 
)

def pca_to_anisou_components(axes, rotation):
    """
    Converts PCA results (principal axes lengths, 2*sigma, and rotation matrix) 
    into the six anisotropic displacement parameters (Uij) for the PDB ANISOU record.
    
    The semi-axes of the ellipsoid are defined as the standard deviation (sigma). 
    The PDB U matrix is the variance-covariance matrix (U = sigma^2).
    
    Args:
        axes (np.array): Principal axis lengths (2*sigma1, 2*sigma2, 2*sigma3).
        rotation (np.array): Rotation matrix (eigenvectors).
        
    Returns:
        list: [U11, U22, U33, U12, U13, U23] scaled by 10000 and rounded to integer.
    """
    # 1. Convert the full axis length (2*sigma) to standard deviation (sigma)
    sigmas = axes / 2.0
    
    # 2. Calculate the variance matrix (U_diag) in the principal axis frame: U = sigma^2
    variances = (sigmas**2)*2
    U_diag = np.diag(variances)
    
    # 3. Rotate U_diag back to the original coordinate system: U = R * U_diag * R_transpose
    U_full = rotation @ U_diag @ rotation.T
    
    # 4. Extract the six unique components: U11, U22, U33, U12, U13, U23
    U11 = U_full[0, 0]
    U22 = U_full[1, 1]
    U33 = U_full[2, 2]
    U12 = U_full[0, 1]
    U13 = U_full[0, 2]
    U23 = U_full[1, 2]
    
    # 5. Scale by 10000 and convert to integer as required by PDB ANISOU format
    # PDB uses 7-digit integer fields for Uij * 10^4
    scaled_U = [
        int(round(U11 * 10000)),
        int(round(U22 * 10000)),
        int(round(U33 * 10000)),
        int(round(U12 * 10000)),
        int(round(U13 * 10000)),
        int(round(U23 * 10000)),
    ]
    
    return scaled_U


def save_ellipsoids_as_pdb(all_fitted_ellipsoids, cage_atoms, atom_names, filename="ellipsoids_dbscan_output.pdb"):
    """
    Generates a PDB file containing the cage atoms and the fitted ellipsoid centroids.
    
    The ellipsoid dimensions (a, b, c axes) and orientation are specified using the
    ANISOU record, allowing visualization software to draw the shape.
    Each ellipsoid is represented by a single HETATM record at its centroid.
    """
    print(f"\n--- Saving results to PDB file: {filename} ---")
    
    pdb_lines = []
    atom_serial = 1

    # PDB Header Information
    pdb_lines.append(f"HEADER    DBSCAN ELLIPSOID ANALYSIS\n")
    pdb_lines.append(f"REMARK 1  This file contains the host 'cage_atoms' and the centroids of the\n")
    pdb_lines.append(f"REMARK 1  detected 'window' clusters, represented as HETATM records.\n")
    pdb_lines.append(f"REMARK 1  Ellipsoid shape is defined by the ANISOU records following each HETATM.\n")
    pdb_lines.append(f"REMARK 1  Axis Lengths (2-sigma) are included in REMARK 2 for easy reading.\n")

    # 1. Write the original cage atoms (for context)
    # Use ATOM records for standard structure visualization
    for i, (x, y, z) in enumerate(cage_atoms):
        # Format: ATOM, serial, atom name (e.g., C), residue name (e.g., CAG), chain, res_seq, x, y, z, occupancy, temp_factor, element
        line = "ATOM  {:5d} {:<4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:<2s}  \n".format(
            atom_serial, atom_names[i], 'CAG', 'A', 1, x, y, z, 1.00, 0.00, atom_names[i]
        )
        pdb_lines.append(line)
        atom_serial += 1
    
    # 2. Write a TER record to mark the end of the first structure (atoms)
    pdb_lines.append("TER\n")
    
    # 3. Write the fitted ellipsoid centroids and ANISOU data
    for cluster_id, ellipsoid_data in enumerate(all_fitted_ellipsoids):
        center = ellipsoid_data['center']
        axes = ellipsoid_data['axes']
        rotation = ellipsoid_data['rotation']
        
        # axes are [a_axis_len, b_axis_len, c_axis_len] from PCA
        a_len, b_len, c_len = axes[0], axes[1], axes[2]
        
        res_seq = cluster_id + 1
        residue_name = f'EL{res_seq:02d}' # e.g., EL01, EL02
        chain_id = 'E' # Ellipsoid chain
        
        # --- 3a. Write Ellipsoid Dimensions as REMARKs (for human readability) ---
        pdb_lines.append(f"REMARK 2  Ellipsoid for Window Cluster ID: {cluster_id} (Residue {residue_name}, Chain {chain_id})\n")
        pdb_lines.append(f"REMARK 2  Axis Lengths (2-sigma): A={a_len:.4f}, B={b_len:.4f}, C={c_len:.4f}\n")
        pdb_lines.append(f"REMARK 2  Center: X={center[0]:.3f}, Y={center[1]:.3f}, Z={center[2]:.3f}\n")

        # --- 3b. Write Centroid as a single HETATM record ---
        atom_name = 'CEN'
        element = 'O'
        alt_loc = ''
        
        line = PDB_FORMAT_HETATM.format(
            atom_serial, atom_name, residue_name[:3], chain_id, res_seq,
            center[0], center[1], center[2], 1.00, 50.00, element # B-factor 50.00 for highlight
        )
        pdb_lines.append(line)
        
        # --- 3c. Calculate and Write ANISOU record ---
        anisou_components = pca_to_anisou_components(axes, rotation)
        U11, U22, U33, U12, U13, U23 = anisou_components
        
        anisou_line = PDB_FORMAT_ANISOU.format(
            atom_serial, atom_name, alt_loc, residue_name[:3], chain_id, res_seq,
            U11, U22, U33, U12, U13, U23, element 
        )
        pdb_lines.append(anisou_line)

        atom_serial += 1
        
    # 4. Write the final END record
    pdb_lines.append("END\n")

    # Write all lines to the file
    try:
        with open(filename, 'w') as f:
            f.writelines(pdb_lines)
        print(f"Successfully saved {atom_serial - 1} records (plus ANISOU) to {filename}.")
        print("Centroids saved as HETATM. Ellipsoid shape defined by ANISOU records.")
        print("You can now open this file in visualization software like PyMOL or VMD.")
    except Exception as e:
        print(f"Error saving PDB file: {e}")


