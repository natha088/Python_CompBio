# Python_CompBio
Python Scripts for Various Computational Biology Tasks

# hbl_all_res
This function calculates the water-protein hydrogen bond autocorrelation function on a per residue basis for a molecular dynamics (MD) trajectory using MDAnalysis.
Modify "structure.psf" and "structure.dcd" for your MD simulation. In addition modify: 

r = list(u.residues.resids[0:100]) --> modify 0:100 for sepcific number of residues in your protein
HBL_analysis = HBL(u, selection1, selection2, 0, 6250, 30) --> modify 6250 with number of frames in your trajectory; and 30 for time sliding window (ps)

Adapted from: https://docs.mdanalysis.org/1.0.1/documentation_pages/analysis/waterdynamics.html


# clustering
This script uses RDKit to perform chemical clustering using the Tanimoto similarity on a vector of compounds inputted as SMILES IDs, outputting the chemical structures in each cluster and plotting a PCA plot to visualize chemical diversity. The code also allows to enter features for each compound, in this case docking scores. This feature can be changed to any other computational or experimental metric derived for each specific compound. The script plots a boxplot of the inputted feature grouped by chemical cluster.


# ANM_flucs
This script uses Prody to perform Anisotropic Network Model (ANM) analysis on a protein inputted as PDB. The protein can be monomeric or multimeric and can contain het atoms such as small molecule ligands, the only requirement is that each molecular entity is labelled with a unique chain identifier. The script outputs the square fluctuations, motional correlation matrix and Hessian matrix as text files, as well as a heatmap plot of the motional correlations.
