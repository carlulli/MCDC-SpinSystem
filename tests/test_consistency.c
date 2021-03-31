/******************************************************************************
Test to check whether hamiltonian and therefore, spin products produce the same
results for different spin modules but the same configuration
********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "geometry.h"
// #include "MCMC.h"
#include "spin.h"
#include "observables.h"


int main(int argc, char const *argv[]) {
  /*
    to do:
    - set parameters (geo and phys)
    - initialize hot or cold start for one spinmodule
    - calculate and safe state energy
    - initialize hot or cold start for other spinmodule
    - calculate and safe state energy
  */
  // Set inital parameters
  set_params();
  set_physics_params();
  spinstruct_t *spinstruct_arr = set_spinarray();
  srand(get_seed());
  int hotstart = atoi(argv[9]); // if hotstart=1 -> hotstart, else cold start
  int N=get_N(), D=get_D();
  double energy_dense1, energy_dense2;

  // set initial spin orientations
  if (hotstart==1) { hot_start();// function that sets spins to special spin configuration }
  else { cold_start(); }

  energy_dense1 = energy_density();

  spinmodel = (1-spinmodel);
  energy_dense2 = energy_density();

  printf("\nThis test calculates the energy density twice with the same spin configuration.\n"
  "but with opposite spinmodels. (1=clock model and 0=ising model)\n");
  printf("Energy density of spinmodel %d = %f\n" "Energy density of spinmodel %d = %f\n",
    1-spinmodel,
    energy_dense1,
    spinmodel,
    energy_dense2);

    if (energy_dense1==energy_dense2) { printf("GREAT SUCCESS!\n"); }
    else { printf("Consistancy FAILED.\n");  }


  return 0;
}
