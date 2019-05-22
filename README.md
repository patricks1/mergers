# Galaxy Merger Uncertainties Project
This repository contains my scripts and Jupyter notebooks that 
- SHAM galaxies to the halos in Andrew Wetzel's treepm package 
- Identifies halo, subhalo, and galaxy mergers
- Calculates merger statistics
- Randomly varies the SHAM inputs a number of times and measures the resulting variance in merger statistics

The .py scripts do the heavy lifting, namely ''merg_rates.py''. The Jupyter notebooks serve to import the classes and methods written in the .py scripts, generate instances of the former, and execute the latter. If you want to understand what's going on, it's best to find a plot of interest in one of the .ipynb's, see what method generates the data, trace that back to the .py script, and read the code there.

## Most Important Scripts
### merg_rates.py
This is the heart of the project. It imports ''subhal_io_hack.py'', which it basically uses to generate an instance of the treepm simulation, and imports ''my_sham_hack.py'', which it uses to SHAMgalaxies to halos in the simulation. The script contains most of the methods I use to analyze galaxy merger statistics. Some of those methods are part of the shamedTreepmClass. Some are written as functions outside of the class.

Within this script, the most important methods are ''merg_tree(...)'' and ''mp_tree(...)''. ''merg_tree(...)'' generates merger trees for each galaxy. ''mp_tree(...)'' generates a dictionary for each galaxy, specifying what its index (i.e. location in the simulation data) was at each redshift, re 

### simdatbuilder.py
This is the script that I ultimately import into an ipython terminal and loop through any number of times to generate statistics for multiple runs of the simulation.

### my_sham_hack.py
This is my modified version of ChangHoon Hahn's ''sham_hack.py'', which he built from Dr. Wetzel's original ''sham.py''. I implemented a few error workarounds, added more stellar mass functions, and implemented functionality that keeps each galaxy's deviation from the mean SHMR constant over the duration of the simulation.

### subhalo_io_hack.py
The TreepmClass here basically generates an instance of the treepm simulation. The instance contains catalogs, broken up my redshift index, containing all the halo information for each redshift. Andrew Wetzel wrote the original subhalo_io and sham scripts. ''subhal_io_hack.py'' is my slightly modified version of that original. 

## Most Important Jupyter Notebooks
Be warned, the notebooks are somewhat messy. They are my ~on-the-fly implementations of the .py scripts. I use(d) the notebooks to generate plots for papers, presentations, and meetings.

### presentation_plots.ipynb
This notebook generates the main plots I use for my uncertainties papers and presentation.

### testing_hash_method.ipynb
At one point, identifying mergers and calculating statistics took ~12 hours for each run of the simulation. On Kilian Walsh's advice, I rewrote the heart of the process in ''merg_rates.py'' to generate merger "hash tables" (i.e. dictionaries) that have essentially zero lookup time. This sped up each run to ~6 minutes! This notebook tests merger statistics using that hash method to make sure it works. It does!  

## Note
This repository works for me and is a good representation of my work, but it will not run without the treepm backend data files, a few utilities, some comparison data files that I pulled from other papers, and all the correct environment paths specified so the environment knows where to find each package being imported.
