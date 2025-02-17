from Bio.PDB import PDBParser
import numpy as np

# Parse the PDB file
parser = PDBParser(QUIET=True)
structure = parser.get_structure('NTD', '5mdi_ntrd22_redo_hit.pdb')

# Iterate over all residues and print small molecules
for model in structure:
    for chain in model:
        for residue in chain:
            # Heteroatoms have non-blank IDs
            if residue.id[0] != ' ':
                print(f'Small molecule: {residue.resname}, Chain: {chain.id}, Residue ID: {residue.id}')

from Bio.PDB import PDBIO

# Save the structure back to a PDB file
io = PDBIO()
io.set_structure(structure)
io.save('mod_structure.pdb')

import prody
from prody import *

# Load the modified structure
structure = parsePDB('5mdi_ntrd22_redo_hit.pdb')

# Print all atoms to verify the small molecule is included
print(structure)
print(f'Total number of atoms selected: {structure.numAtoms()}')

# Select all atoms, including the small molecule
calphas = structure.select('(protein and name CA) or (chain U)')
print(calphas)

# Initialize ANM
anm = ANM('ANM analysis')

# Build Hessian matrix using all atoms
anm.buildHessian(calphas)
hessian = anm.getHessian()
np.savetxt('hessian_matrix.txt', hessian, fmt='%.6f', header='Hessian Matrix')

anm.calcModes()

# Get fluctuations
fluctuations = calcSqFlucts(anm[:])

# Get residue numbers
residues = calphas.getResnums()
chains = calphas.getChids()

# Save the residue numbers and fluctuations to a text file
with open('residue_fluctuations.txt', 'w') as file:
    file.write('Residue\tChain\tFluctuation\n')
    for res, chain, fluct in zip(residues, chains, fluctuations):
        file.write(f'{res}\t{chain}\t{fluct}\n')


correlation_matrix = calcCrossCorr(anm)

# Save the correlation matrix to a text file
np.savetxt('correlation_matrix.txt', correlation_matrix[0:156,0:156], delimiter='\t')
unique_chains = list(set(chains))
# Extract and save correlations within each chain
for chain in unique_chains:
    chain_indices = [i for i, ch in enumerate(chains) if ch == chain]
    chain_corr_matrix = correlation_matrix[np.ix_(chain_indices, chain_indices)]
    np.savetxt(f'{chain}_corr.txt', chain_corr_matrix, delimiter='\t')

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Step 1: Load the correlation matrix from the file
correlation_matrix = np.loadtxt('correlation_matrix.txt', delimiter='\t')

# Step 2: Create residue labels (2 to 79 for the first set of residues and 1 to 79 for the second part)
residues_2_to_79 = np.arange(2, 80)  # Residues from 2 to 79 (total 78)
residues_1_to_79 = np.arange(1, 80)  # Residues from 1 to 79 (total 79)

# Step 3: Create a total label set for x and y axes (2-79 + 1-79)
# Combine the two sets of residue labels into one array for the axes
all_residues = np.concatenate([residues_2_to_79, residues_1_to_79])

# Step 4: Plot the heatmap with the custom labels
plt.figure(figsize=(12, 10))
sns.heatmap(correlation_matrix, xticklabels=all_residues, yticklabels=all_residues, cmap='coolwarm', annot=False, fmt='.2f',vmin=-0.8, vmax=1)

# Set the labels and title
plt.xlabel('Residue')
plt.ylabel('Residue')
plt.title('Correlation Matrix Heatmap')

# Adjust the font size of the tick labels (both x and y axes)
plt.xticks(fontsize=6)  # Adjust the font size for x-axis labels
plt.yticks(fontsize=6)  # Adjust the font size for y-axis labels
ticks_to_show = np.arange(0, len(all_residues), 6)  # Show every 5th residue
plt.xticks(ticks_to_show, all_residues[ticks_to_show], fontsize=8)  # Set x-axis ticks
plt.yticks(ticks_to_show, all_residues[ticks_to_show], fontsize=8)  # Set y-axis ticks
# Display the heatmap with demarcated regions
plt.savefig('heatmap.png',dpi=300)
plt.show()

