import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import matplotlib.pyplot as plt
from rdkit.Chem import Draw

# Load Excel file with a specific range (e.g., rows 10-30)
df = pd.read_excel("order_cmps_120_redo.xlsx", engine="openpyxl", skiprows=0, nrows=120)  # Adjust indices accordingly

# Extract SMILES from the first column
smiles_list = df.iloc[:, 0].dropna().tolist()
scores = df.iloc[:, 1].dropna().tolist()

smiles_list = smiles_list[::-1]
scores = scores[::-1]

# Convert SMILES to RDKit molecule objects
molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list if smiles]

# Print the molecules
for i, mol in enumerate(molecules):
    if mol:
        print(f"Molecule {i+1}: {Chem.MolToSmiles(mol)}")
    else:
        print(f"Molecule {i+1}: Invalid SMILES")

fpgen = AllChem.GetRDKitFPGenerator()
fps = [fpgen.GetFingerprint(x) for x in molecules]

from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

# Function to compute the full pairwise distance matrix
def fingerprint_distance_matrix(fps):
    num_fps = len(fps)
    dists = []
    for i in range(num_fps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])  # Compute similarities
        dists.extend([1 - sim for sim in sims])  # Convert to distances
    return dists

# Compute full distance matrix
dists = fingerprint_distance_matrix(fps)

# Perform clustering (0.2 threshold means 80% similarity cutoff)
clusters = Butina.ClusterData(dists, len(fps), 0.25, isDistData=True)

# Print cluster sizes
for i, cluster in enumerate(clusters):
    print(f"Cluster {i+1}: {len(cluster)} molecules")


# Assign cluster labels to molecules
cluster_labels = np.zeros(len(fps)) - 1  # Default: -1 for unassigned
for i, cluster in enumerate(clusters):
    for idx in cluster:
        cluster_labels[idx] = i

df.columns = ["ID", "scores"]
df["Cluster"] = cluster_labels

# Count the number of compounds in each cluster
cluster_sizes = df["Cluster"].value_counts()

# Filter for clusters with more than 1 compounds
valid_clusters = cluster_sizes[cluster_sizes > 1].index

# Filter the DataFrame to only include these clusters
df_filtered = df[df["Cluster"].isin(valid_clusters)]

import seaborn as sns

# Plot boxplot with all points as individual dots
plt.figure(figsize=(10, 6))
sns.boxplot(x="Cluster", y="scores", data=df_filtered, showfliers=False)  # Boxplot without outliers
sns.stripplot(x="Cluster", y="scores", data=df_filtered, color='black', size=4, jitter=True)  # Add individual points
plt.title("Docking Score Distribution by Cluster (Clusters > 1 Compounds)")
plt.xlabel("Cluster")
plt.ylabel("Docking Score")
plt.xticks(rotation=45)
plt.show()


# Calculate the mean of 'scores' for each cluster
mean_scores = df_filtered.groupby('Cluster')['scores'].median().sort_values()
# Create a list of clusters ordered by the median score
ordered_clusters = mean_scores.index
# Plot boxplot with all points as individual dots, ordered by average scores
plt.figure(figsize=(10, 6))
sns.boxplot(x="Cluster", y="scores", data=df_filtered, order=ordered_clusters, showfliers=False)  # Boxplot without outliers
sns.stripplot(x="Cluster", y="scores", data=df_filtered, color='black', size=4, jitter=True, order=ordered_clusters)  # Add individual points
plt.title("Docking Score Distribution by Cluster (Clusters > 1 Compounds)")
plt.xlabel("Cluster")
plt.ylabel("Docking Score")
plt.xticks(rotation=45)
plt.savefig('boxplot.png',dpi=300)
plt.show()


# Iterate through all clusters
for i, cluster in enumerate(clusters):
    # Only process clusters with at least 2 molecules
    if len(cluster) > 1:
        # Extract the SMILES strings for the current cluster
        cluster_smiles = [smiles_list[j] for j in cluster]
        
        # Convert SMILES strings to molecules
        cluster_molecules = [Chem.MolFromSmiles(smiles) for smiles in cluster_smiles if smiles]
        
        # Generate the image of the molecules in the current cluster
        img = Draw.MolsToGridImage(cluster_molecules, molsPerRow=5, subImgSize=(300, 300))

        # Optional: Save the image to a file
        img.save(f'cluster_{i}_image.png')

# all 120
img = Draw.MolsToGridImage(molecules, molsPerRow=20, subImgSize=(300, 300))
img.save(f'all.png')

import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from sklearn.decomposition import PCA
from rdkit.Chem import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

# Convert fingerprints to a NumPy array
def fingerprints_to_array(fps):
    arr = np.zeros((len(fps), len(fps[0])), dtype=np.int8)
    for i, fp in enumerate(fps):
        DataStructs.ConvertToNumpyArray(fp, arr[i, :])
    return arr

# Convert fingerprints to a matrix
fp_matrix = fingerprints_to_array(fps)

# Perform PCA (reduce to 2D for visualization)
pca = PCA(n_components=2)
pca_result = pca.fit_transform(fp_matrix)


plt.figure(figsize=(10, 7))
scatter = plt.scatter(pca_result[:, 0], pca_result[:, 1], c=scores, cmap="viridis", alpha=0.7)
plt.xlabel("PCA Component 1")
plt.ylabel("PCA Component 2")
plt.title("PCA of Molecules Colored by Scores")
plt.colorbar(label="Score")
plt.savefig('pca.png',dpi=300)
plt.show()

