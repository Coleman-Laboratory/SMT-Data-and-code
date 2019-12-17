# SMT-Data-and-code
Single Molecule Tracking (SMT) of Halo-PCNA and SNAP-PolD WT/DEAD

These parameters are used for slimFAST_2 generated localizations and binding trajectories(tracking) followed by processing from evalSPT_3.
slimFast2 and evalSPT_3 are compatible with Matlab R2014b.

parameters for slimFAST_2 localizations:
Error Rate: 10^-6
Detection: Box 11 pixels
Deflation Loops: 15
Optimization Settings
max # iterations: 50
term. tol: 10^-2
r0 tolerance: 50%
max. post refinement: 1.5

parameters for slimFAST_2 Acquisition:
Pixel Size: 0.084 um
Emission: 567nm
N.A. : 1.49
PSF Scaling: 1.2
PSF Std: 1.03 pixels
Counts/Photon: 20.2
Lag Time= 1340 ms

parameters for slimFast_2 Tracking:
Max. expected Diffusion Coefficient: 0.05 um^2/s
Searchexp. Factor: 1.2
Statistics Win.: 10 frames
Max. # Competitors: 3
Max. OFF-Time: 3 frames
Int. Fluc. Weight: 0.9
Local vs. Max. expected Diffusion Weight: 0.5


The output from evalSPT_3 contains 8 columns with the first four columns containing the relevant information required for downstream data processing.

Column 1 X-position of molecule
Column 2 Y-position of molecule
Column 3 Frame ID
Column 4 Binding Trajectory/Molecule ID

The output from evalSPT_3 is then used as input for Matlab script spt_process_v15CR.m file. The resulting data structure (sptana) contains metadata, nuclear segmented binding trajectories, and nuclear segmented dwell times for binding trajectories (ndwet). Each row is a separate cell. Sptana data structures should from 2 channel SMT should be organized to have subsequent rows containing cells from channel 1 data followed by subsequent rows of matched cells from channel 2. For example, 10 cells of 2 channel SMT should have rows 1-10 channel 1 information from Cells A-J followed by rows 11-20 channel 2 information from Cells A-J.

The sptana data structure containing matched cells with 2 channels of tracking is used by the livecellColocV1CR.m file to isolate spatially and temporally colocalized molecules. The output This program saves the following values, percentages of molecules colocalized in each cell (perColoc), percentages of molecules temporally colocalized in each cell (perClcT),aggregated arrival/departure times of channel/color 2 molecules with respect to channel/color 1 (ADBag), aggregated colocalization times (ClcTag), and a cell structure containing a table with all molecular IDs, arrival departure frames, XY position, arrival/departure times of channel/color 2 with respect to channel/color 1 (FinalMoldata)). 

The output from FinalMoldata file is:
column 1: average X position of colocalized binding event in channel 1
column 2: average Y position of colocalized binding event in channel 1
column 3: Molecular ID of colocalized binding event in channel 1
column 4: start frame of colocalized binding event in channel 1
column 5: end frame of colocalized binding event in channel 1
column 6: average X position of colocalized binding event in channel 2
column 7: average Y position of colocalized binding event in channel 2
column 8: Molecular ID of colocalized binding event in channel 2
column 9: start frame of colocalized binding event in channel 2
column 10: end frame of colocalized binding event in channel 2
column 11: arrival time (seconds) of molecule in channel 2 with respect to arrival time of molecule in channel 1
column 12: departure time (seconds) of molecule in channel 2 with respect to departure time of molecule in channel 1
column 13: time (seconds) of colocalization
column 14: residence time of colocalized binding event in channel 1
column 15: residence time of colocalized binding event in channel 2

The sptana data structure is used by sptclst_v24CR.m file to generate maps of hubs(clusters) of binding events. The resulting data structure (sptclst) contains metadata, cluster numbers, locations and Molecular IDs of binding events in clusters.
