# Earthquake iso-nuisance and iso-damage simulation

This respository holds code for MATLAB (version R2019a) to simulate the potential for nuisance and damage by earthquakes in Alberta. This code is associated with a paper titled *Earthquake iso-nuisance and iso-damage mapping for Alberta: applications for choosing magnitude thresholds to manage induced seismicity*, by Mauricio Reyes Canales, Elwyn Galloway, Steven Pawley, Javad Yusifbayov, and Greg Hartman. All authors are part of the Alberta Geological Survey.

By default, the script *AB_Simulation.m* performs simulations for all of Alberta. This script relies on some other scripts in this repository. The simulations are saved in the file *simulations.map*, though this file name can be configured.

After simulations are complete, *AB_Mapping.m* retrieves the simulations saved in *simulations.mat* and creates maps and plots of the results. The nuisance and damage tolerance thresholds are configured in this script.

## Performing a test simulation

Simulation time strongly depends on the simulation grid resolution, the size of the area being simulated, the number of realizations, and the number of earthquake magnitudes being simulated. Final simulations for our study took 24 to 48 hours to complete. It is recommended that a test simulation is performed to ensure the scripts are working correctly before pursuing a full-resolution simulation.

Simulation parameters all appear in the first ~60 lines of *AB_Simulation.m*. By changing the following parameters to the suggested values, a test simulation will take approximately 5 minutes:
* `source_dx = 0.4`
* `source_dy = 0.4`
* `Mag_rng = (2 : 1.5 : 5)`
* `Nr = 3`

Note: it may take a moment to read the large *Vs30_AB_Ext_3km_lit.mat* file.

## Vs30 Simulation

A time-averaged shear wave velocity for the upper 30 m of the land surface (Vs30, Holzer et al., 2005) was created for Alberta based on compiled shear wave velocity values for different bedrock and surficial lithological units. One hundred realizations of stochastic simulation were generated to estimate the Vs30 value on the impact grid. These realizations are included in this repository as *Vs30_AB_Ext_3km_lit.mat*. See the publication's Supplemental Information for more details.




