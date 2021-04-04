#include "observables.h"
#include "utilities.h"
#include "geometry.h"
#include "spin.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
// initialize the static variables
// static int N;
// static int D;
// static int B;
double * spinval_values;
extern spinstruct_t *spinstruct_arr;

/*************************************************************************
magnetization density

***************************************************************************/
double mag_density(void)
{
  int N = get_N();
	int D = get_D();
  double total_magnetization = 0.0;
	int spinmodel = get_spinmodel();


    if (spinmodel == 0){ // Ising !
      for(int x=0; x< int_pow(N, D); x++)
      {
        total_magnetization  += spinstruct_arr[x].spinval;
      }
    }
    else if (spinmodel == 1){ // Clock !
      for(int x=0; x< int_pow(N, D); x++)
      {
        total_magnetization += spinval_values[spinstruct_arr[x].spinval]   ;  // takes the computed value of the exponential
          //from the precalculated vector   ! check if this is ok
      }
    }
  return total_magnetization/(int_pow(N, D));
}

// energy function   aka the hamiltonian


double energy_density(void)
{
	int N = get_N();
	int D = get_D();
	int B = get_B();
	int spinmodel = get_spinmodel();
  double total_energy = 0.0;
  double sum_next_neighbour_interaction = 0.0;

  int contrib_neighbor_index = 0;

	int prod_neighbor_orientation = 0 ;


  // loop over all the spins -> for each spin position the contribution of all the neighbours gets calculated

	if (spinmodel == 0) // Ising
  {
    for(int x = 0; x < int_pow(N, D) ; x++)
    {
      sum_next_neighbour_interaction = 0;            //for each spin the spins of each numbers get summed up
      for(int y = 0; y < 2*D ; y++)              //loop for the getting the sum over next neighbour spins
      {
        contrib_neighbor_index = spinstruct_arr[x].nnidx[y];
				if (contrib_neighbor_index == -1){
					continue;
				}
				// small m for the spin product of two neighbors
        prod_neighbor_orientation = spinmultiplication(spinstruct_arr[x].spinval, spinstruct_arr[contrib_neighbor_index].spinval);
				sum_next_neighbour_interaction += prod_neighbor_orientation; //
      }
      //total energy. Ising  in each iteration the energy contribution of one spin site is added
      total_energy -=  spinstruct_arr[x].spinval * B/*magnetic field*/ ; // !
			total_energy -= sum_next_neighbour_interaction/2 ;

    }
  }

  else if (spinmodel ==1) // clock
  {
    for(int x = 0; x < int_pow(N, D) ; x++)     // iterating over all the lattice sites
    {
      sum_next_neighbour_interaction = 0;            //for each spin the spins of all the neighbors get summed up
      for(int y = 0; y < 2*D ; y++)              //loop for the getting the sum over (2D) next neighbour spins
      {
        contrib_neighbor_index = spinstruct_arr[x].nnidx[y];
				if (contrib_neighbor_index == -1){
					continue;
				}   // gets the index of one contributing neighbor
				prod_neighbor_orientation = spinmultiplication(spinstruct_arr[x].spinval, spinstruct_arr[contrib_neighbor_index].spinval);  // orientation of product
				sum_next_neighbour_interaction += spinval_values[prod_neighbor_orientation]; // ! adds the real part
      }
      //total energy. Clock.In each iteration the energy contribution of one spin site is added
      total_energy -=  spinval_values[spinstruct_arr[x].spinval] *  B/*magnetic field*/;  // check  !
      total_energy -= sum_next_neighbour_interaction/2 ;   // changed formula
		}
  }

	//total_energy -= sum_next_neighbour_interaction ;   // the contribution of the next neighbor interaction gets added

  return total_energy/(int_pow(N, D));      //divide by size to get density
}
