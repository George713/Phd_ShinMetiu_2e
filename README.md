# PhD_ShinMetiu_2e

This is the source code for simulations of a Shin-Metiu system extended by an additional electron. If you want to quickly use this code outside of the ARA-cluster of the Friedrich-Schiller-University Jena, Germany, you have to modify Makefile and src/modules to link to the correct libraries.

Source files are found in src/, while bin/ will contain the code's binary files after building. out/ will contain the program's output files.

File explanation:

input:
	Contains most system variables. Most noteworth: potential & grid properties including screening parameters (R_c, R_f & R_e), properties for the wave packet initializing system dynamics, properties for pump & probe pulses, length of propagation time step, spatial symmetry of electronic eigenstates (used for initializing the system with a certain spin configuration)

sc_vv:
	SCript file used for handling job attributes (as of this writing slurm is being used). 'vv' is used for differentiating between configurations. For starting a job on ARA simply use 'sbatch sc_01' with 01 replacing vv.

src/main.f90:
	Master file bringing together the different elements of the calculation.

src/prop3d.f90:
	The heart of the program, as it solves the 3d TDSE. This piece of the program contributes the most to overall computational time. First, the initial wave packet is defined, followed by solving the TDSE via split-operator methode. By applying a cutoff-function to the wave function, an outer part is separated from the bound wave function. This outer part is considered ionized and propagated with an ionic potential. During propagation several observables are computed for later analysis. After solving the TDSE, the ionization spectrum and its asymmetry is computed.

src/potential.f90:
	Code calculating the potential arising from the 2e-Shin-Metiu model.

src/adiabatic.f90:
	Code for calculating electronic eigenstates and eigenenergies via relaxation methode.

src/input.f90:
	Code used for reading input variables from the input file.

src/modules.f90
	File declaring all variables used throughout the program. Most importantly, grid sizes are defined here. Also, conversion factors are set here.

src/subroutines.f90:
	Contains several subroutines used in one or more of the larger code parts.

src/nuclear_wf.f90:
	Calculates vibrational eigenstates via relaxation methode. Often not used.

src/prop1d.f90:
	Propagation in the Born-Oppenheimer approximation. Estimation of ionized part via perturbation theory. Often not used.