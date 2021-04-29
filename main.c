/***************************************************+***************************
Main code:
- reads in parameter file
- initializes lattice, i.e., geomerty
  -> maybe initialze a pointer in the main that can be given to the geometry and is from there on the N^d array with all structs inside
- initialze hot or cold start
- loop through iterations of MCMC
  - write out energy and magnetization after each iteration
  - wirte out correlation function after maybe each 10th iteration
  - write out spin configuration after maybe each 100th iteration
- free array
***************************************************+***************************/


// lets get started with this


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "spin.h"
#include "observables.h"
#include "geometry.h"
#include "utilities.h"
#include "MCMC.h"

spinstruct_t *spinstruct_arr;
double *spinval_values;
double *spin_corr_vec;

int main(int argc, char const *argv[])
{
	if (argc != 11 ) {
		printf(" Hello, this is not the right amount of input parameters. \n"
		" The program requires the following 10 parameters (You had: {%d})\n  "
		"Remember: [filename] [N] [D] [boundary condition] [ordering] [spinmodel] [M_number_orient] [double B] [seed] [double T] [start_choice]\n", argc);
		exit(-1);
	}

	set_params(argc, argv); // setting up the geometry
	set_physics_params(argc,argv);  // and the physical parameters

	int N = get_N() ;
	int D = get_D() ;
	double B = get_B();
	int boundary_condition = get_boundary_condition();
	int ordering = get_ordering();
	int M_number_orient = get_M_number_orient();
	int spinmodel = get_spinmodel();
	int	seed = get_seed();
	double T = get_T();
	printf("\n \nYour chosen system: \n \n Size in one direction = %d\t Dimension = %d \t boundary_condition (0 for Dirichlet and 1 for periodic) = %d \t ordering (0 for lexo and 1 for blackwhite) = %d \n \n spinmodel (0 for Ising and 1 for Clock) = %d \t M_number_orient = %d \n\n"
	 " magnetic field B = %f \t seed = %d \t temperature T = %f     \n ",N, D, boundary_condition, ordering, spinmodel, M_number_orient, B, seed, T );


	spinval_values = set_spinval_values();  // creating the vector which has the real spinvalues for their orientation
	//spin_corr_vec = set_spincorrelation_vector();

	spinstruct_arr = set_spinarray();

	double magnetization = 0.0; // variables to store the magnetization and energy density in.
	double energy = 0.0;

// opening textfile for saving the values
FILE *fp, *fp2, *fp3;
int namesize = 60;
for (int i=1; i<=9; i++) { namesize += strlen(argv[9]); }
char filename[namesize];
char filename2[namesize];
char filename3[namesize];
snprintf( filename, sizeof(filename),
  "data/monte_carlo_em_%s_%s_%s_%s_%s_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9]);

snprintf( filename2, sizeof(filename),
	"data/monte_carlo_corr_%s_%s_%s_%s_%s_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9]);
snprintf( filename3, sizeof(filename),
	"data/monte_carlo_config_%s_%s_%s_%s_%s_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9]);

// file for the energy density and the magnetization density.
fp = fopen(filename, "w");
fprintf(fp, "\n energy density \tmagnetization_density\n");
//fclose(fp);
// file for ths spin-correlations.
fp2 = fopen(filename2, "w");
fprintf(fp2, "\n r \tcorrelation\n");
//fclose(fp2);
// file for the spin-configuration.
fp3 = fopen(filename3, "w");
fprintf(fp3, "\n idx \t spinvalue\n");
//fclose(fp3);
//fp = fopen(filename, "w");
// fp2 = fopen(filename2, "w");
// fp3 = fopen(filename3, "w");
/* Monte Carlo loop */
hot_start();
for (int i = 0 ; i< 100000; i++)
{
	MCMC_step(); // does one Monte Carlo step, thereby given the spins a new configuration
	energy = energy_density();
	magnetization = mag_density();
	// prints out the energy and magnetization (densities) every run.

	fprintf(fp,"%.16e\t%.16e\n",energy, magnetization);



	if ((i%10) == 0 ) // prints out the spin correlation every 10th run.
	{
		spin_corr_vec = set_spincorrelation_vector();

		// fprintf(fp2,"Monte Carlo round:%d \n",i );
		for (int j = 0; j< N/2; j++ )
		{
			fprintf(fp2,"%d\t%.16e\n",j, spin_corr_vec[j]);
		}

	}
	if ((i%1000) == 0 ) // print out the spin configuration every 100th run.
	{

		// fprintf(fp3,"Monte Carlo round:%d \n",i );
		for (int f = 0; f< int_pow(N,D); f++ )
		{
			if (spinmodel == 0)
			{
				fprintf(fp3,"%d\t%d \n",f, spinstruct_arr[f].spinval );
			}
			else if (spinmodel == 1)
			{
				fprintf(fp3,"%d\t%.16e\n",f, spinval_values[spinstruct_arr[f].spinval] );
			}
		}

	}
	//fclose(fp);
	// fclose(fp2);
	// fclose(fp3);

}
fclose(fp);
fclose(fp2);
fclose(fp3);




/* freeing the no longer needed vectors */
	free_spin_corr_vec();
	free_spinval_values();
	spinarray_free(spinstruct_arr);
	return 0 ;
}
