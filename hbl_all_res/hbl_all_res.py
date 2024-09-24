import MDAnalysis
import matplotlib.pyplot as plt
import matplotlib
import numpy
import math
from MDAnalysis.analysis.waterdynamics import HydrogenBondLifetimes as HBL

u = MDAnalysis.Universe('6pes_wb.psf', '6pes_100.dcd')
r = list(u.residues.resids[0:64])
time = 0
for x in r:
        selection1 = "byres name OH2 and sphzone 5.0 protein and resid " + str(x)
        selection2 = "resid " + str(x)
        HBL_analysis = HBL(u, selection1, selection2, 0, 6250, 30)
        HBL_analysis.run()
        for HBLc, HBLi in HBL_analysis.timeseries:
    	    print("{time} {HBLc}".format(time=time, HBLc=HBLc))
    	    time += 1

