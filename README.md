# behavior
Code and data used in the paper "Game-theoretic modeling of collective decision-making during epidemics"

INSTRUCTIONS AND SYSTEM REQUIREMENTS:

The code for the simulation is written in MATLAB, using version 2021a and requires no add-ons. All the functions needed are included in the corresponding folders. Hence, to run the code it is sufficient to open the corresponding folder with MATLAB.

The code has been tested on a PC with 16GB RAM and CPU 1.9 GHz with OS Windows 10 (64-bit). The run time of a simulation with n=20,000 individuals and T=600 time-steps is approximately 2400 and 1400 seconds for the SIS ans SIR models, respectively.

ORGANIZATION:

The repository contains four folders: one for each of the two epidemic model (SIS and SIR), and one that contains the raw data of all the figures, and one that contains a demonstration of the output of a similations and the script used to generate it. To run it, it is sufficient to open the folder with MATLAB and run "plot_graph_dynamics_example.m". The run tume is of few second on a standard PC.

The folders sis and sir contain the code to simulate the SIS and SIR epidemic process coupled with the coevolutionary model, respectively. They also contain the scripts used to generate the corresponding figures in the paper, that is, Figure % (or the panel * of figure %) can be generated by running the file "code_for_fig#.m" or "code_for_fig%_a.m."

The folder "raw_data_figures" contains all the .tex files used to generate the figures in the paper

DEMO

The folder "Demo" contains a sample output of the simulation in the file "demo_output" and the code use to generate it. To run it, it is sufficient to open the folder with MATLAB and run "code_for_demo.m". The sample is a SIS model with n=10,000 individuals (similar to Figure 3a). The run tume is of approximately 1000 second on a standard PC.

LICENCE
Our code has no custom licence, just requires a standard MATLAB licence to run.

