/******************************************************************************
Test to check if full magnetization can be reached for when temperature is
close to zero -> accepting only better energy spin flips
Input values:
[N] [D] [boundary_condition] [ordering] [spinmodel] [M_number_orient] [B] [seed] [T]

this needs way more iterations than proposed in the report!!!
and should probably only work correctly with the statistical analysis part
********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "geometry.h"
#include "MCMC.h"
#include "spin.h"
#include "observables.h"

extern spinstruct_t *spinstruct_arr;
extern double *spinval_values;
// extern double *spin_corr_vec;

int main(int argc, char const *argv[]) {
  /*
  program starts a 4x4 simulation with T=0 and writes out magetization after 100 runs
  for boundary condition? -> PERIODIC
  T=0, D=2, N=4, ordering=1 have to be input values
  */
  set_params(argc, argv);
  set_physics_params(argc, argv);
  // extern spinstruct_t *spinstruct_arr = set_spinarray();
  spinstruct_arr = set_spinarray();

  // void init_MCMC(void);

  int T, ordering, boundary_condition, spinmodel, seed;
  T = get_T();
  ordering = get_ordering();
  boundary_condition = get_boundary_condition();
  spinmodel = get_spinmodel();
  seed = get_seed();

  // spin_corr_vec = set_spincorrelation_vector();
  spinval_values = set_spinval_values();

  if (T>1.e-7 || ordering!=1 || boundary_condition!=1 || spinmodel!=1) {
    printf("Magentization test error. Remember:\n"
    "- Temperatur small < 1.e-7\n"
    "- Ordering black white = 1\n"
    "- Boundary Condition periodic = 1\n");
  }
  printf("This test runs 100 MCMC spin flip iterations with T->0 and for a small 2D lattice.\n"
  "This should take long enough for all spins to flip in one direction.\n"
  "For a 4x4 lattice with Ising model 16 spins should point in the same direction.\n");

  hot_start();
  srand(seed);

  int iterations=1000000;
  for (int i=0; i<iterations; i++) {
    MCMC_step();
    // printf("[testfullmag.c] DEBUGGING \n");
    // exit(-1);
    printf("Spin orientation values:\n"
    " %d, %d, %d, %d\n"
    " %d, %d, %d, %d\n"
    " %d, %d, %d, %d\n"
    " %d, %d, %d, %d\n",
    spinstruct_arr[14].spinval, spinstruct_arr[6].spinval, spinstruct_arr[15].spinval, spinstruct_arr[7].spinval,
    spinstruct_arr[4].spinval, spinstruct_arr[12].spinval, spinstruct_arr[5].spinval, spinstruct_arr[13].spinval,
    spinstruct_arr[10].spinval, spinstruct_arr[2].spinval, spinstruct_arr[11].spinval, spinstruct_arr[3].spinval,
    spinstruct_arr[0].spinval, spinstruct_arr[8].spinval, spinstruct_arr[1].spinval, spinstruct_arr[9].spinval );
  }


  printf("Simulation test result: Magnetisation density = %f\n", mag_density());

  // free_spin_corr_vec(); // free up the no longer needed get_arrays
  free_spinval_values();
  spinarray_free(spinstruct_arr);

  return 0;
}
