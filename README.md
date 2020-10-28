# Uniaxial Compaction and Force-chain Analysis of Bidisperse Grain packs
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4021433.svg)](https://doi.org/10.5281/zenodo.4021433)

This repository contains the script for creating a bidisperse (two radii) grain pack under gravity and running a uniaxial compaction test on it using LIGGGHTS library (C++) on high-performance computing (HPC) resources. Scripts for post-processing like calculating, plotting and visualizing coordination numbers and force chains to study the grain behavior when subjected to compaction are also provided.

----------------------------------------------------------------------------------------------------------------------------------

## Pre-requisities

*i) LIGGGHTS*

Please install LIGGGHTS on your system. You can find downloading instructions [here](https://www.cfdem.com/media/DEM/docu/Section_start.html).

Note - If working with Windows, one option is to install the [ubuntu terminal](https://ubuntu.com/wsl) and then install LIGGGHTS as a library with pre-compiled codes. This is likely to reduce the possible issues compared to the other methods.

*ii) MATLAB/OCTAVE*

MATLAB or OCTAVE is required for the post-processing the results and plotting graphs.

## Recommended

*i) ParaView* 

[ParaView](https://www.paraview.org/download/) will allow visualization of the grain pack and force chain behavior during the compaction process. The ParaView visualization works with version 4.0.1 and requires importing custom filter [NetworkFlow.cpd](https://github.com/abhishekdbihani/compaction_LIGGGHTS/blob/master/post-processing/NetworkFlow.cpd).

----------------------------------------------------------------------------------------------------------------------------------

## Procedure & Details

### I) LIGGGHTS script

Navigate to the src folder using the linux command line, and run the following code (adjust the number of CPU processors to be used)-

 ``` mpirun -np 6 liggghts < in.uniaxialTest  ```

The in.uniaxialTest script consists of three main stages-

*a) Initializing input variables (Lines 3 to 76)* 

*b) Creating & modifying the grain pack (Lines 77 to 101)*

*c) Uniaxial compaction (Lines 102 to 141)* 
  
The radius of the grains can be changed on lines 46 (small) and 47 (large)
The concentration by grain volume can be changed on lines 48 (small) and 49 (large)

On running the script, walls are imported to form the domain, and then grains are inserted upto step 10000, allowed to settle and grains outside set domain removed at step 16000. Compaction steps are from 16000 to 36000. The code returns an output every 2000 steps in form of force chain and grain properties. Both the output files are exported in two file types (without extension in tabular form and .vtp/vtr for opening in ParaView).
    
### II) Post-processing

Run forcechain_ana.m script in MATLAB for analyzing the strong/weak force chain networks and their relation with Large-Large, Large-Small, Small-Small grain contacts. It calculates the different properties as a function of compaction and saves it in post-processing.csv and also creates a .vtk file which can be used to study strong/weak force chains and grain positions in ParaView. 

The plot_creator.m script can be used to plot the change in calculated variables as a function of axial strain. 

The ParaView visualization works with version 4.0.1 and requires importing custom filter NetworkFlow.cpd. The grains can be observed using the NetworkPores filter and the force-chains using the NetworkThroats filter. 

The exported grain center coordinates and radii (xyzr.csv files) can be used as inputs for running an invasion percolation code. 

----------------------------------------------------------------------------------------------------------------------------------

## Example

The figure below shows the visualization in PARAVIEW of a grain pack subjected to uniaxial testing in three stages- 1) No compaction, 2) Limited Compaction, 3) Final Compaction with a) 3D structure, b) 2D cross-section, and c) force chains for each stage.

<img src="https://github.com/abhishekdbihani/compaction_LIGGGHTS/blob/master/example%20compaction%20picture.png" align="middle" width="900" height="550" alt="compaction visualization" >

Datasets like [this](https://www.digitalrocksportal.org/projects/204) one can also be prepared using this code.

----------------------------------------------------------------------------------------------------------------------------------

## Citation

If you use this repository, please cite as-

Bihani, A., & Daigle, H. (2020, September 9). Uniaxial Compaction and Force-chain Analysis of Bidisperse Grain Packs (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.4021433

## Related Publications

Bihani A., Daigle H. (in-review). Seal Capacity, Force Chains, and Percolation in Silt-Clay Mixtures. Journal of Geophysical Research- Solid Earth. https://doi.org/10.1002/essoar.10504349.1

