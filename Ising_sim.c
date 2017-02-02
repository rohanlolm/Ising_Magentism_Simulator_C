/* Project 2 - Ising Model
Rohan Mehra: 640052491
PHYM004 
Date Submitted: 20/01/15

Files included in submission: 
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

./Ising_sim [-opt] <output_filename>

The program can also be called with any combination of 7 options: 

Options: 
-m: Set the maximum number of iterations for equilibrium to be reached at a given temperature. Default: 1,000,000,000.
-t: Set equilibrium threshold. Equilibrium is found once the relative change in energy is no more than this value ever N time steps, 
	where N is the number of spins in the lattice (N = DxD). Default 1e-4
-i: Set initial temperature. The maximum and initial temperature of the system in units of J/kb. Default: 3.5
-f: Set final temperature, not inclusive. Final and minimum temperature will be final_temp + temp_step_size Default 0
-s: Set temperature step size as the system is cooled. Default 0.1
-d: Set the data record ratio. This number determines for how long energy and magnetisation values are recorded and averaged to calculate 
	the macroscopic properties. A value of 2 means macroscopic quantities are calculated by averaging values for the same amount of time 
	it took to get to equilibrium. A value of 3 would mean values are averaged for twice the time it took to get to equilibrium. Default: 2
-a: Enable state_writer() with ASCII art. This option turns on statewriter() before and after the system is temperature evolved. 
	This prints out the spin state of each particle and represents this in an ‘ASCII art’ format to visualise domains. 
*/

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <assert.h> 
#include <string.h> 
#include <errno.h> 
#include <time.h>
#include <stdbool.h>
#include <getopt.h>

#define muB 0			//External applied field
#define k_B 1			//1.38065e-23
#define D 100			//Size of the spin lattice - must be int. Change to 50 for a significantly faster execution time. 
#define J 1			//Exchange energy
#define SEED time(NULL)		//Random Seed
#define S(x,y) S[(x)%D][(y)%D]	//S access macro to simulate infinite lattice with periodic boundary conditions
#define FP_OFFSET 4		//flip_probs array offset 
#define JUMPSIZE (D*D)		//Increments to check for equilibrium 
#define DOWN -1.0		//Spin down 
#define UP 1.0			//Spin up

//Defaults - can be changed with options
#define EQ_THRESHOLD 1e-4	//Threshold value to find equilibrium 
#define INIT_TEMP 3.5		//Initial temperatures in units of J/k_B (maximum)
#define FINAL_TEMP 0 		//Final temperature in units of J/k_b (minimum)
#define MAX_STEPS 1000000000	//Maximum number of steps to find equilibrium
#define TEMP_STEP 0.1		//Temperature increments (cools)
#define DATA_RECORD_MULTIPLE 2.0	//Calculates physical properties using the same amount of time it took to get to equilibrium 

//**************************************************************************************************************************************************
//Structures
struct EM_sums{ 					//Stores information returned by Metrop_evolve()
	double sum_E_sq, sum_E, sum_M_sq, sum_M;	//sum_E_sq is the sum of E squared
	int counter, equilibrium_time; 			//counter is the amount of time the energies and magnetisations are averaged over 
};							//equilibrium_time is the time to equilibrium

struct properties{ 					//Stores physical properties calculated by properties_calc()
	double temperature, heat_capacity, susceptability, energy, magnetisation;
	int equilibrium_time; 
}; 

struct options{ 					//Stores options given by user
	int max_steps;
	double eq_threshold, init_temp, final_temp, temp_step, data_record_multiple; 
	bool write_state;			
};

//Prototypes *********************************************************************************************************************************
//Explanations given at each function at bottom of file. 
//xfunctions - Implementations of common functions with error checking which exit in the case of a failure. 
static FILE* xfopen(char *filename, char *option);

//File Writers - Write various properties and information to file. 
static int state_writer(char *data_filename, char *ASCII_filename, double S[D][D]);
static int properties_writer( struct properties *props, FILE* outfile);

//Populators - Populate arrays as required.
static int rand_pop_S(double S[D][D]);
static int pop_flip_probs(double *flip_probs, double temp);

//Property Calculators - Calculate specific physical properties of the system.
static int properties_calc( struct EM_sums *pEM_sums, struct properties *props, double S[D][D], double temp );
static double heat_capacity_calc( struct EM_sums *pEM_sums, double temp );
static double susceptability_calc( struct EM_sums *pEM_sums, double temp );
static double total_magnetisation_per_particle( double S[D][D] );
static double total_energy_per_particle( double S[D][D] );

//Other
static double get_rand_spin(void);
static int find_equilibrium(double current_energy, double old_energy, struct options *p_options );
static int Metrop_evolve( double S[D][D], double *flip_probs, struct EM_sums *pEM_sums, struct options *p_options, double temp );
static int temperature_evolve( double S[D][D], char *filename, struct options *p_options );
inline static int get_rand(void);

//*******************************************************************************************************************************
//Processes user inputs and initialises S as well as calls temperature_evolve. 
int main(int argc, char **argv){
	srand(SEED);					//Seed the generator 
	double S[D][D];					//The array of spins
	rand_pop_S( S );				//populate S with random spins
	char opt;
	struct options opts = { MAX_STEPS, EQ_THRESHOLD, INIT_TEMP, FINAL_TEMP, TEMP_STEP, DATA_RECORD_MULTIPLE, false };  //Struct with options

	while( ( opt = getopt(argc,argv,"m:t:i:s:d:f:a") ) !=-1 ) {			//Use getopt to loop through arguments
		if(opt == 'm') opts.max_steps = atoi(optarg); 
		else if(opt == 't') opts.eq_threshold = atof(optarg); 
		else if(opt == 'i') opts.init_temp = atof(optarg); 			//For arguments that need arguments, convert from string to double/integer
		else if(opt == 'f') opts.final_temp = atof(optarg);
		else if(opt == 's') opts.temp_step = atof(optarg);
		else if(opt == 'd') opts.data_record_multiple = atof(optarg);
		else if(opt == 'a') opts.write_state = true; 
	}
	if( (argc-1) != optind ){						//Check for incorrect number of 
			fprintf(stderr, "\nIncorrect number of arguments.\n"); 	//There should be one argument provided by user
			exit(EXIT_FAILURE);
	}

	if( opts.write_state ) state_writer( "S_list_start.txt", "S_ASCII_start.txt", S );	//Write initial and final states of S

	temperature_evolve( S, argv[optind], &opts );		//Evolve S through temps! 

	if( opts.write_state ) state_writer( "S_list_end.txt", "S_ASCII_end.txt", S );	
	return 0;
} 

//Takes user options, filename and S and evolves S through desired temperature range. Writes results to file.
static int temperature_evolve( double S[D][D], char *filename, struct options *p_options ){ 
	double temp; 
	double flip_probs[2*FP_OFFSET + 1];
	struct EM_sums EM_sums1 = {0};			//Initialise structs out of loop so memory is not reallocated  
	struct properties props;
	FILE *outfile = xfopen(filename, "w");		//Open file for writing
	fprintf(outfile, "#Ising Model Simulation. Lattice Size: %d, Eq. Threshold: %f, Data Record Ratio: %f\n", D, p_options->eq_threshold, p_options->data_record_multiple);
	fprintf(outfile, "#%s \t %20s \t %20s \t %20s \t %20s %20s\n\n", "Temperature", "Heat Capacity", "Susceptibility", "Final Energy", "Final Magnetisation", "Equilibrium Time");
	
	for(temp = p_options->init_temp; temp > p_options->final_temp; temp -= p_options->temp_step){ 	//Loop through temps
		pop_flip_probs( flip_probs, temp ); 					//populate Boltzmann probabilities at current temperature
		Metrop_evolve( S, flip_probs, &EM_sums1, p_options, temp );		//Metropolis algorithm at each temperature
		properties_calc( &EM_sums1, &props, S, temp );				//Calculate properties for each temperature
		props.temperature = temp;						
		properties_writer( &props, outfile ); 					//Write properties to file
		EM_sums1 = (struct EM_sums){0};						//Reset value of summing struct to 0 at each iteration
	}

	fclose(outfile);						
	return 0; 
}

//Applies the Metropolis algorithm to evolve the system in time until equilibrium is reached at each temeperature. 
static int Metrop_evolve( double S[D][D], double *flip_probs, struct EM_sums *pEM_sums, struct options *p_options, double temp ){ 	//pass pEM_sum as argument rather than return so we aren't allocating in a loop
	int x = 0;		//Lattice positions
	int y = 0;
	int time_step = 0; 
	bool isEquilibrium = false; 			//Are set to true if equilibrium found
	bool confirm_equilibrium = false;
	double uni_rand; 			//Uniformly distributed rand
	double delta_E, delta_M; 		//Change in Energy and Mag due to spin flip 
	double sum_adj_spin;			//Sum of adjacent spins 
	double flip_prob;			//Boltzmann probability to flip 
	double energy_shift, mag_shift; 	//Shift value used to calculate variance to prevent 'catastrophic cancellation' in floating point arithmetic. 
	double energy = total_energy_per_particle( S );		//total system energy
	double mag = total_magnetisation_per_particle( S );	//total system magnetisation
	double old_energy = energy;				//used in finding equilibrium
	const double inverse_N = 1.0 / (double)(D * D);			//Floating point multiplication faster than division. 

	while( time_step < p_options->max_steps ){ 
		time_step++;  
		x++;
		if(x%D == 0) y++;								//Go through sequentially equilibralise much faster
		
		sum_adj_spin = S(x-1,y) + S(x,y-1) + S(x+1,y) + S(x,y+1);
		delta_E = ( S(x,y) * ( 2.0 * J * sum_adj_spin + muB ) * inverse_N );		
		delta_M = -2*S(x,y) * inverse_N; 
		
		if( S(x,y) == DOWN ) flip_prob = 1 / flip_probs[(int)sum_adj_spin + FP_OFFSET]; 	//Exponential gets negated. Need to flip probabilities if spin down
		else flip_prob = flip_probs[(int)sum_adj_spin + FP_OFFSET]; 
		
		uni_rand = get_rand() / (double)RAND_MAX; 			
		if( uni_rand < flip_prob ){ 							//don't need to check delta_E < 0 as this case captures that as well 
			S(x,y) = -S(x,y); 

			energy 	+= delta_E;							//Modify total energy and magnetisation if spin flipped
			mag 	+= delta_M; 	
		}
	
		if( time_step % JUMPSIZE == 0 && !isEquilibrium ){ 				//Every JUMPSIZE and if not equilibrium check for equilibrium
			isEquilibrium = find_equilibrium(energy, old_energy, p_options); 			 
			old_energy = energy;
			if( isEquilibrium && !confirm_equilibrium ) {				//To ensure that noise didn't accidentally trigger equilibrium 
				confirm_equilibrium = true; 					//we check equilibrium is found twice in a row
				isEquilibrium = false;
			}
			else if( !isEquilibrium && confirm_equilibrium )	{		
				confirm_equilibrium = false;
			}
			else if( isEquilibrium && confirm_equilibrium )	{			//Equilibrium is found twice in a row we store equilibrium time
				pEM_sums->equilibrium_time = time_step;				//Energy and mag is stored to shift the variance calculation
				energy_shift = energy;
				mag_shift = mag;
				//state_writer( "S_Eqdata.txt", "S_Eqascii.txt", S );		//Uncomment to get state at equilibrium 
			}
		}

		if( isEquilibrium && confirm_equilibrium ){ 					//If equilibrium store energy and magnetisation data						
			pEM_sums->sum_E_sq	+= ( (energy - energy_shift)*(energy - energy_shift) ); 
			pEM_sums->sum_E	 	+= (energy - energy_shift);
			pEM_sums->sum_M_sq 	+= ( (mag - mag_shift)*(mag - mag_shift) ); 			 
			pEM_sums->sum_M 	+= (mag - mag_shift); 
			pEM_sums->counter++;
		}
		if( isEquilibrium && confirm_equilibrium && time_step > p_options->data_record_multiple * pEM_sums->equilibrium_time ) { //Break after a certain multiple of equilibrium time 
			printf("Equilibrium Found at Temp: %g Steps: %d\n", temp, pEM_sums->equilibrium_time );		//Only line temp variable is used in function. Temp argument not necessary 
			break; 		
		}	
	}
	if( !(isEquilibrium && confirm_equilibrium) ){ 
		printf("Failed to find Equilibrium. In file this is designated as Equilibrium Time = 0 \n");
	}
	return 0; 
} 

//xfunctions **********************************************************************************************************************
//Opens file, exits if file opening fails
static FILE* xfopen(char *filename, char *option){ 
	FILE *retp = fopen(filename, option); 
	if( retp == NULL ){ 
		fprintf(stderr,"File %s failed to open. \n",filename); 
		exit(EXIT_FAILURE); 
	}
	return retp;
}

//File Writers *******************************************************************************************************************
//Writes the spins of S at a particular time. Also creates an 'ASCII art' representation of the spin array S. 
static int state_writer(char *data_filename, char *ASCII_filename, double S[D][D]){ 	
	FILE *datafile = fopen(data_filename, "w"); 
	FILE *ASCII_file = fopen(ASCII_filename, "w");
	
	fprintf(datafile, "# x\ty\tSpin State \n");
	fprintf(datafile, "#_____________________\n");
	for(int x = 0; x < D; x++){ 
		for( int y = 0; y < D; y++){
			fprintf( datafile, "%d\t%d\t%g\n", x, y, S(x,y) );		//Write coordinate and spin
			if( S(x,y) == UP ) fprintf( ASCII_file,"+" );			//Spin UP is drawn as + 
			if( S(x,y) == DOWN ) fprintf( ASCII_file," " );			//Spin DOWN is blank
		}
		fprintf(ASCII_file, "\n");
	}
	fclose(datafile); 
	fclose(ASCII_file); 
	return 0; 
}

//writes properties as tab separated values
static int properties_writer( struct properties *props, FILE* outfile){ 	
	fprintf(outfile, "%f \t %20e \t %20e \t %20f \t %20f \t %20d \n", props->temperature, props->heat_capacity, props->susceptability, props->energy, props->magnetisation, props->equilibrium_time); 
	return 0;
} 

//Populators *******************************************************************************************************************
//Populate S with randomly orientated spins
static int rand_pop_S(double S[D][D]){
	for(int x = 0; x < D; x++) { 
		for(int y = 0; y < D; y++){ 
			S[x][y] = get_rand_spin();
		}
	}
	return 0;
}

//Populate the flip_probs array with Boltzmann probabilities. This is to avoid having to constantly calculate the exponential for every time step 
static int pop_flip_probs(double *flip_probs, double temp){
	double delta_E;
	
	for(int i = -FP_OFFSET; i < 5; i++){ 				//The sum of the adjacent spins can go from -4 to 4 hence the FP_OFFSET
		delta_E = 2.0*J*i + muB;
		flip_probs[i+FP_OFFSET] = exp( -delta_E / (k_B * temp) ); 
	}
	return 0;
}

//Property Calculators *********************************************************************************************************
//Calculates all the properties in one place and modifies the props struct.
static int properties_calc( struct EM_sums *pEM_sums, struct properties *props, double S[D][D], double temp ){

	props->heat_capacity	= heat_capacity_calc( pEM_sums, temp ); 
	props->susceptability 	= susceptability_calc( pEM_sums, temp );
	props->energy 		= total_energy_per_particle( S );
	props->magnetisation	= total_magnetisation_per_particle( S ); 
	props->equilibrium_time = pEM_sums->equilibrium_time;

	return 0;
}

//I think the following functions are self explanatory

static double heat_capacity_calc( struct EM_sums *pEM_sums, double temp ){

	double time_ave_E_sq = pEM_sums->sum_E_sq / (double)pEM_sums->counter; 			//typecast safety as counter integer 
	double E_time_ave_sq = (pEM_sums->sum_E / (double)pEM_sums->counter)*(pEM_sums->sum_E / (double)pEM_sums->counter); 
	return 1/(k_B * temp*temp)*(time_ave_E_sq - E_time_ave_sq); 
}

static double susceptability_calc( struct EM_sums *pEM_sums, double temp ){ 
	
	double time_ave_M_sq = pEM_sums->sum_M_sq / (double)pEM_sums->counter; 
	double M_time_ave_sq = (pEM_sums->sum_M / (double)pEM_sums->counter)*(pEM_sums->sum_M / (double)pEM_sums->counter); 
	return 1/(k_B * temp)*(time_ave_M_sq - M_time_ave_sq); 
}

static double total_energy_per_particle( double S[D][D] ){			//Sum energy due to each particle 
	double E_total = 0; 
	for(int x = 0; x < D; x++){ 
		for( int y = 0; y < D; y++){
			E_total -= S(x,y)*( J*(S(x+1, y) + S(x, y+1)) + muB ); 	//Doesn't double count
		}
	} 
	return E_total / (double)(D*D); 				//typecast safety as D is integer
} 

static double total_magnetisation_per_particle( double S[D][D] ){  		//Sum magnetisations from each particle 
	double mag = 0;
	for(int x = 0; x < D; x++){ 
		for( int y = 0; y < D; y++){
			mag += S(x,y); 
		} 
	}
	return mag / (double)(D*D); 
} 

//Other ******************************************************************************************************************************************
//Returns a random spin value
static double get_rand_spin(void){			//This won't be perfectly uniform but since we're modding by a very small number, 2, the skew will 
	int r = get_rand(); 				// be of the order 1/RAND_MAX - which is very small. 
	if(r%2 == 0) return UP; 
	else	return DOWN; 
}

//Takes two values for energy and options and decides if equilibrium is reached
static int find_equilibrium(double current_energy, double old_energy, struct options *p_options ){
	double change = fabs( (old_energy - current_energy) / current_energy ); 	//Divide by current_energy as old_energy could be asymptotically large	
	if( change < p_options->eq_threshold ) return true; 				//resulting in a small change incorrectly
	else return false; 
}

//This function exists encase you want to change the random number generator to something else later without having to greatly change the code 
inline static int get_rand(void){ 	
	return rand();			 
}
