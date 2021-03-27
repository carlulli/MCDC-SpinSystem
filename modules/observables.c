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
  for(x=0, x< pow(N, D), x++)
  {
    total_magnetization  += spinstruct_arr[x].spinval;
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
  // loop over all the spins -> for each spin position the contribution of all the neighbours gets calculated
  for(x = 0, x < (pow(N, D) , x++)
  {
    sum_next_neighbour_spins = 0;            //for each spin the spins of each numbers get summed up


    for(y = 0, y < 2*D , y++)              //loop for the getting the sum over next neighbour spins
    {
      contrib_neighbor_index = spinstruct_arr[x].nnidx[y]
      sum_next_neighbour_spins += spinstruct_arr[contrib_neighbor_index].spinval
    }

    //total energy. in each iteration the energy contribution of one spin is added
    total_energy -=  spinstruct_arr[x].spinval * ( sum_next_neighbour_spins + B/*magnetic field*/);

  }

  return total_energy/(pow(N, D);      //divide by size to get density


}





// two point spin correlation functions
