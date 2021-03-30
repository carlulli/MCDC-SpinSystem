// need some tests for the observables

// the function magnetization density should be tested for the aligned spins and for the antialigned spin_spin_correlation_r

// for the aligned case we just excpect N^d times the spin   and for the antialigned case we excpect zero magnetization.

// first the aligned because

/******************************************************************************
Test of magnetization density for all spins aligned:
----------------------------------------------------
Procedure: 1. align all spins with function -> void aligned_start(spinstruct *spinstruct_arr, int direction)
					 2. calculate magnetization density with double magentization_density(spinstruct *spinstruct_arr)
					 3.


********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spin.h"
#include "observables.h"
#include "geometry.h"


double magnetization;
double energy_density;
double spin_corr_vec[N/2];
double spinval_values[M_number_orient] ; // !
double B ;
int N, D;



int main(int argc, char const *argv[])
{

	set_params(argc, argv);

	spinstruct_t *spinstruct_arr = set_spinarray();

	int N, D ;
	N = get_N() ;
	D = get_D() ;



	/******************************************************************************
	Test of magnetization density for all spins aligned:
	----------------------------------------------------
	Procedure: 1. align all spins with function -> void aligned_start(spinstruct *spinstruct_arr, int direction)
						 2. calculate magnetization density with double magentization_density(spinstruct *spinstruct_arr)
						 3.


	********************************************************************************/

// !! need to initialize seed somewhere


	printf("Your chosen system: \n Size in one direction = %d\t Dimension = %d \t boundary_condition = %s \t ordering = %s \n",
	    N, D, boundary_condition, ordering);  // model missing ! and optionally direction


	cold_start();  // for the Ising model it puts all the spins in 1 and for the Clock model m= 0, which results in spin 1 ;

	magnetization = magentization_density();    // calculate magnetization density

	printf("\n The magnetization density for the cold (aligned) start is: %f \n ",magnetization) ;

	printf("\n The expected magnetization density for the cold (aligned) start is exactly 1. \n ");


/*******************************************************************************
Test of magnetization density for the hot start
----------------------------------------------------
Expectation: We excpect the magnetization density to have random values, but with a mean of 0.
Procedure: 1. align spins randomly with hot_start
					 2. calculate magnetization density with double magentization_density(spinstruct *spinstruct_arr)

! at the moment it only prints out the magnetizations for 20 different hot starts.
********************************************************************************/
	printf("The magnetization density is expected to be zero in average. For different hot starts, we obtain the following magnetization_densities. \n", );
	for(i=0, i<20, i++)
	{
		hot_start();
		magnetization = magnetization_density() ;
		printf(" round %d: %f \t ",i, magnetization) ;
	}



/*******************************************************************************
test of the two-point spin-spin correlation function:
-----------------------------------------------------
Goal: Testing it for 3 different types of configuration
			1. for cold starts
			2. for a hot start
			3. for an anti-aligned start



********************************************************************************/
printf("\n \n Testing of the two-point spin-spin correlation function \n \n");

//cold start

printf("1. test for a cold start configuration. \n The resulting spin correlation vector is: \n")
cold_start();
double spin_corr_vec[N/2];  // initializing the vector, which we fill with the different correlations for the different r.
spincorrelation_vector(spin_corr_vec);
for(i=0, i<20,i++)
{
	printf(" r= %d corr = %f ",i,spin_corr_vec[i]);     // printing the correlations from r = 0 to r = 20
}

// hot start

printf("1. test for a hot start configuration. \n The resulting spin correlation vector is: \n")
hot_start;
spincorrelation_vector(spin_corr_vec);
for(i=0, i<20,i++)
{
	printf(" r= %d corr = %f ",i,spin_corr_vec[i]);     // printing the correlations from r = 0 to r = 20
}



/********************************************************************************
Test of the energy density:
---------------------------
Goal: See if the energy density for the cold and hot start are what we excpect
*********************************************************************************/

/* cold start: for each spin we have H = -B - 2D , so for a system of N length and Dimensions D, we expect H = N^D*(-B - 2D) and energy density = -B -2D
*/
printf("\n \n Test of the energy density: \n \n");
printf("1. test: The cold start \n");
cold_start();
energy_density = energy_density();
double expected_energy_density = - B - 2*D ;
printf("The energy density for the cold start is: %f. Whereas we expected: %f. \n",energy_density, expected_energy_density );


/* hot start: The expected energy density should be zero in average.
*/
printf("2. test: The hot start \n");
printf("Here we want to try a couple of different hot starts, to see whether the average energy density could be zero, which is what we expect.\n");
for(i=0,i<20,i++)
{
	hot_start();
	energy_density = energy_density();
	printf(" round:%d has energy_density = %f",i, energy_density);
}



}  //end of main
