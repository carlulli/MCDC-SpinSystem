#include "observables.h"
#include <math.h>
// initialize the static variables


/*************************************************************************
magnetization density

***************************************************************************/
double magentization_density(spinstruct *spinstruct_arr)
{
  // get size
  double total_magnetization = 0.0;

    if (spinmodule == 0){ // Ising
      for(x=0, x< pow(N, D), x++)
      {
        total_magnetization  += spinstruct_arr[x].spinval;
      }
    }
    else if (spinmodule == 1){ // Clock
      for(x=0, x< pow(N, D), x++)
      {
        total_magnetization += spinval_values[spinstruct_arr[x].spinval]     // takes the computed value of the exponential
          //from the precalculated vector#
      }
    }
  }
  return total_magnetization/(pow(N, D));
}











// energy function   aka the hamiltonian


double energy_density(spinstruct *spinstruct_arr)
{
  // !still need!-> functions to get the size, the magnetic field and the dimension of the problem

  double total_energy = 0.0;
  double sum_next_neighbour_spins = 0.0;
  int contrib_neighbor_index = 0;
  double real_part_for_m_orientation = 0.0 ;
  // loop over all the spins -> for each spin position the contribution of all the neighbours gets calculated
  for (spinmodule == 0) // Ising
  {
    for(x = 0, x < pow(N, D) , x++)
    {
      sum_next_neighbour_spins = 0;            //for each spin the spins of each numbers get summed up


      for(y = 0, y < 2*D , y++)              //loop for the getting the sum over next neighbour spins
      {
        contrib_neighbor_index = spinstruct_arr[x].nnidx[y]
        sum_next_neighbour_spins += spinstruct_arr[contrib_neighbor_index].spinval
      }

      //total energy. Ising  in each iteration the energy contribution of one spin site is added
      total_energy -=  spinstruct_arr[x].spinval * ( sum_next_neighbour_spins + B/*magnetic field*/);

    }
  }

  for (spinmodule ==1) // clock
  {
    for(x = 0, x < pow(N, D) , x++)     // iterating over all the lattice sites
    {
      sum_next_neighbour_spins = 0;            //for each spin the spins of all the neighbors get summed up


      for(y = 0, y < 2*D , y++)              //loop for the getting the sum over (2D) next neighbour spins
      {
        contrib_neighbor_index = spinstruct_arr[x].nnidx[y]   // gets the index of one contributing neighbor
        real_part_for_m_orientation = spinval_values( spinstruct_arr[contrib_neighbor_index].spinval ) ;  // gets the real part of its spin !
        sum_next_neighbour_spins += real_part_for_m_orientation ; // sums that real part to the sum of next neighbors
      }

      //total energy. Clock.In each iteration the energy contribution of one spin site is added
      total_energy -=  spinval_values( spinstruct_arr[x].spinval ) * ( sum_next_neighbour_spins + B/*magnetic field*/);  // check  !
      // ? need to check wether its actually just the real parts that need to be multiplied 
  }

  return total_energy/(pow(N, D);      //divide by size to get density


}





// two point spin correlation functions
