This repository contains the script for creating a bidisperse grain pack and running a uniaxial compaction test on it using LIGGGHTS. Scripts for post-processing (calculating, plotting and visualizing coordination numbers and force chains) are also provided.

Created by Abhishek Bihani (Supervisor: Hugh Daigle).
May 2020

----------------------------------------------------------------------------------------------------------------------------------

Pre-requisities:

i) LIGGGHTS
ii) MATLAB

Recommended:

i) ParaView 

Please install LIGGGHTS on your system. 
Link- https://www.cfdem.com/media/DEM/docu/Section_start.html

If working with Windows, one option is to install the ubuntu app and then install LIGGGHTS as a library with pre-compiled codes. This is likely to reduce the possible issues compared to the other methods.

MATLAB is required for the post-processing the results and plotting graphs.

ParaView will allow visualization of the grain pack and force chain behavior during the compaction process.

----------------------------------------------------------------------------------------------------------------------------------

Procedure & Details-

I) LIGGGHTS script-

Navigate to the src folder using the linux command line, and run the following code (adjust the number of CPU processors to be used)-

mpirun -np 6 liggghts < in.uniaxialTest

The in.uniaxialTest script consists of three main stages-

a) Initializing input variables (Lines 3 to 76) 
b) Creating & modifying the grain pack (Lines 77 to 101)
c) Uniaxial compaction (Lines 102 to 141) 
  
The radius of the grains can be changed on lines 46 (small) and 47 (large)
The concentration by grain volume can be changed on lines 48 (small) and 49 (large)

On running the script, walls are imported to form the domain, and then grains are inserted upto step 10000, allowed to settle and grains outside set domain removed at step 16000. Compaction steps are from 16000 to 36000. The code returns an output every 2000 steps in form of force chain and grain properties. Both the output files are exported in two file types (without extension in tabular form and .vtp/vtr for opening in ParaView).
    

II) Post-processing-

Run forcechain_ana.m script in MATLAB for analyzing the strong/weak force chain networks and their relation with Large-Large, Large-Small, Small-Small grain contacts. It calculates the different properties as a function of compaction and saves it in post-processing.csv and also creates a .vtk file which can be used to study strong/weak force chains and grain positions in ParaView. 

The plot_creator.m script can be used to plot the change in calculated variables as a function of axial strain. 

The ParaView visualization works with version 4.0.1 and requires importing custom filter NetworkFlow.cpd. The grains can be observed using the NetworkPores filter and the force-chains using the NetworkThroats filter. 

The exported grain center coordinates and radii (xyzr.csv files) can be used as inputs for running the invasion percolation code.      
