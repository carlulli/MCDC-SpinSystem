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
#include "utilities.h"

extern spinstruct_t *spinstruct_arr;
extern double *spinval_values;

int main(int argc, char const *argv[]) {
  /*
    to do:
    - set parameters (geo and phys)
    - initialize hot or cold start for one spinmodule
    - calculate and safe state energy
    - initialize hot or cold start for other spinmodule
    - calculate and safe state energy
  */
  printf("\nThis test calculates the energy density twice with the same spin configuration.\n"
  "but with opposite spinmodels. (1=clock model and 0=ising model)\n"
  "Therefore, the choice of spin orientation values of the clock model is restricted to even values,"
  "that are valid for the Ising model as well.\n");

  // Set inital parameters
  set_params(argc, argv);
  set_physics_params(argc, argv);

  // NEW STUFF BEGINNING
  int number_orient = get_M_number_orient();
  if (number_orient%2 != 0) {
    printf("ERROR Test only works for even number of orientations (M)\n"
    "Change the 5th input value to an even number (not counting the file name)\n");
    exit(-1);
  }
  //for now i couldnt find out a better way then essentially hardcoding 2
  // without changing stuff in other modules
  number_orient = number_orient/number_orient*2;
  // NEW STUFF FINISHED

  spinstruct_arr = set_spinarray();

  int spinmodel = get_spinmodel();
  int seed = get_seed();
  int start_choice = get_start_choice();
  double energy_dense1, energy_dense2;

  spinval_values = set_spinval_values();

  // set initial spin orientations
  if (start_choice==1) {
    srand(seed);
    hot_start(); } // function that sets spins to special spin configuration
  else { cold_start(); }

  energy_dense1 = energy_density();

  spinmodel = (1-spinmodel);
  energy_dense2 = energy_density();


  printf("Energy density of spinmodel %d = %f\n" "Energy density of spinmodel %d = %f\n",
    1-spinmodel,
    energy_dense1,
    spinmodel,
    energy_dense2);

    if (energy_dense1==energy_dense2) { printf("GREAT SUCCESS!\n"); }
    else { printf("Consistancy FAILED.\n");  }

    free_spinval_values();
    spinarray_free(spinstruct_arr);

  return 0;
}
