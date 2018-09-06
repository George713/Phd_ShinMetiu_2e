# PhD_ShinMetiu_2e

This is the source code for simulations of a Shin-Metiu system extended by an additional electron. If you want to quickly use this code outside of the ARA-cluster of the Friedrich-Schiller-University Jena, Germany, you have to modify Makefile and src/modules to link to the correct libraries.

Source files are found in src/, while bin/ will contain the code's binary files after building. out/ will contain the program's output files.

File explanation:

input:
	Contains most system variables. Most noteworth: potential & grid properties including screening parameters (R_c, R_f & R_e), properties for the wave packet initializing system dynamics, properties for pump & probe pulses, length of propagation time step, spatial symmetry of electronic eigenstates (used for initializing the system with a certain spin configuration)

sc_vv:
	SCript File used for handling job attributes (as of this writing slurm is used). 'vv' is used for differentiating between configurations. For starting a job on ARA simply use 'sbatch sc_vv'.

TBC