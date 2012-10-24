Production of simulation and accuracy curve data for the SPME paper (Fig.4 in paper):

1) SIMULATION: run the NaCl.tcl script with espresso
~/Travail/Espresso/espresso-3.0.1/withFT/Espresso NaCl.tcl 

This produces 10 konfig_x files

2) CALCULATION OF REFERENCE FORCES
Run program EwaldSum (XCode project)
This program is in (MacBookPro)/~/Travail/Simulations/P3MAccuracy/EwaldSum
This creates Konfig_X.dat and Konfig_X.exact files that are recognized by Markus SPME program.

3) RUN SPME FORCE CALCULATION
Run program ~/ParticleMesh/AD_withoutSF/SPME-AD-SF.c

