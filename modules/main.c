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

	set_params(argc, argv); // setting up the geometry
	set_physics_params(argc,argv);  // and the physical parameters
	if (argc != 10 ) {
		printf(" Hello, this is not the right amount of input parameters. \n"
    " The program requires the following 10 parameters (You had: {%d})\n  "
		"Remember: [filename] [N] [D] [boundary condition] [ordering] [spinmodel] [M_number_orient] [double B] [seed] [double T] \n", argc);
		exit(-1);
	}

	int N = get_N() ;
	int D = get_D() ;
	double B = get_B();
	int boundary_condition = get_boundary_condition();
	int ordering = get_ordering();
	int M_number_orient = get_M_number_orient();
	int spinmodel = get_spinmodel();
	int	seed = get_seed();
	double T = get_T();
	spinval_values = set_spinval_values();  // creating the vector which has the real spinvalues for their orientation
	spin_corr_vec = set_spincorrelation_vector();

	double magnetization = 0.0; // variables to store the magnetization and energy density in.
	double energy = 0.0;

// opening textfile for saving the values
FILE *fp *fp2, *fp3;
int namesize = 60;
for (int i=1; i<=9; i++) { namesize += strlen(argv[9]); }
char filename[namesize];
snprintf( filename, sizeof(filename),
  "data/monte_carlo_em_%s_%s_%s_%s_%s_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9]);

snprintf( filename2, sizeof(filename),
	"data/monte_carlo_corr_%s_%s_%s_%s_%s_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9]);
snprintf( filename3, sizeof(filename),
	"data/monte_carlo_config_%s_%s_%s_%s_%s_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9]);

// file for the energy density and the magnetization density.
fp = fopen(filename, "w");
fprintf(fp, "\n energy density \tmagnetization_density\n");
fclose(fp);
// file for ths spin-correlations.
fp2 = fopen(filename2, "w");
fprintf(fp2, "\n r \tcorrelation\n");
fclose(fp2);
// file for the spin-configuration.
fp3 = fopen(filename3, "w");
fprintf(fp3, "\n idx \t spinvalue\n");
fclose(fp3);

/* Monte Carlo loop */
for (int i = 0 ; i< 10000; i++)
{
	void MCMC_step(void)(); // does one Monte Carlo step, thereby given the spins a new configuration
	energy = energy_density();
	magnetization = mag_density();
	// prints out the energy and magnetization (densities) every run.
	fp = fopen(filename, "w");
	fprintf(fp,"%.16e\t%.16e\n",energy, magnetization);
	fclose(fp);


	if ((i%10) == 0 ) // prints out the spin correlation every 10th run.
	{
		spin_corr_vec = set_spincorrelation_vector();
		fp2 = fopen(filename2, "w");
		fprintf(fp2,"Monte Carlo round:%d \n",i );
		for (int j = 0; j< N/2; j++ )
		{
			fprintf(fp2,"%d\t%.16e\n",j, spin_corr_vec[j]);
		}
		fclose(fp2);
	}
	if ((i%100) == 0 ) // print out the spin configuration every 100th run.
	{
		fp3 = fopen(filename3, "w");
		fprintf(fp3,"Monte Carlo round:%d \n",i );
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
		fclose(fp3);
	}

}






/* freeing the no longer needed vectors */
	free_spin_corr_vec();
	free_spinval_values();
	spinarray_free(spinstruct_arr);
	return 0 ;
}
