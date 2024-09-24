# Python_CompBio
Python Scripts for Various Computational Biology Tasks

# hbl_all_res
This function calculates the water-protein hydrogen bond autocorrelation function on a per residue basis for a molecular dynamics (MD) trajectory using MDAnalysis.
Modify "structure.psf" and "structure.dcd" for your MD simulation. In addition modify: 

r = list(u.residues.resids[0:100]) --> modify 0:100 for sepcific number of residues in your protein
HBL_analysis = HBL(u, selection1, selection2, 0, 6250, 30) --> modify 6250 with number of frames in your trajectory; and 30 for time sliding window (ps)

Adapted from: https://docs.mdanalysis.org/1.0.1/documentation_pages/analysis/waterdynamics.html
