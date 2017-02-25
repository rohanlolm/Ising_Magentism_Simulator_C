# Ising_Magentism_Simulator_C
Simulates ferro-magnetism using the 2D Ising model implemented in C 

## Rohan Mehra: 
PHYM004 
Date Submitted: 20/01/15

## Files included in submission: 
-Ising_sim.c 					//Source code for simulation 
-Ising Simulator Documentation RM.pdf 		//Documentation and report
-default_results.txt				//Output datafile of simulation with default parameters
-S_ASCII_start.txt				//An 'ASCII art' representation of the initial spin configuration - only outputted if the -a option is applied
-S_ASCII_end.txt				//An 'ASCII art' representation of the final spin configuration - only outputted if the -a option is applied
-S_list_start.txt				//A list of the spin state at every coordinate initially - only outputted if the -a option is applied
-S_list_end.txt					//A list of the spin state at every coordinate finally- only outputted if the -a option is applied

This program simulates ferromagnetism in two spatial dimensions using the statistical mechanics approach of the Ising Model. 
An array of randomly distributed spins is generated, and the Metropolis algorithm is used to evolve the system in time until 
equilibrium is reached (macroscopic properties of the system are constant). At equilibrium the macroscopic properties of 
System Energy, Magnetisation, Specific Heat Capacity and Susceptibility are calculated. 
The system is evolved through a range of temperatures and the macroscopic properties at equilibrium for each temperature are printed to file. 
The source code must be compiled according to C99 standards. 

The program must be called with one command line argument, the output data filename. 
`` ./Ising_sim [-opt] <output_filename>``
The program can also be called with any combination of 7 options: 

## Options: 
-m: Set the maximum number of iterations for equilibrium to be reached at a given temperature. Default: 1,000,000,000.   
-t: Set equilibrium threshold. Equilibrium is found once the relative change in energy is no more than this value ever N time steps, 
	where N is the number of spins in the lattice (N = DxD). Default 1e-4   
-i: Set initial temperature. The maximum and initial temperature of the system in units of J/kb. Default: 3.5    
-f: Set final temperature, not inclusive. Final and minimum temperature will be final_temp + temp_step_size Default 0     
-s: Set temperature step size as the system is cooled. Default 0.1   
-d: Set the data record ratio. This number determines for how long energy and magnetisation values are recorded and averaged to   calculate    
	the macroscopic properties. A value of 2 means macroscopic quantities are calculated by averaging values for the same amount of time 
	it took to get to equilibrium. A value of 3 would mean values are averaged for twice the time it took to get to equilibrium. Default: 2    
-a: Enable state_writer() with ASCII art. This option turns on statewriter() before and after the system is temperature evolved. 
	This prints out the spin state of each particle and represents this in an ‘ASCII art’ format to visualise domains.    
