/***************************************************+***************************
Method:
we use the 1. proposed updating method (refering to desrciption) -> MARS step for each spin value proposal
for one iteration: We update first all even then all odd
Structure:
1. loop over all even spins:
    - update all even spins
2. MARS (accept/reject) step
3. loop over all odd spins:
    - update all odd spins
    - background can be calculated with already changed spins (but only if MARS was performed on them)
4. MARS step
***************************************************+***************************/
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "geometry.h"
#include "spin.h"

extern spinstruct_t *spinstruct_arr;
extern double *spinval_values;
// extern double *spin_corr_vec;

static inline void propose_spinvals(bool even);
static double one_spin_hamiltonian(int x, bool bb_spin);
static inline double min(double a, double b); // logically better in utilities, but since used here -> simpler
static void MARS(bool even);

void init_MCMC(void) {
  spinval_values = set_spinval_values();  // creating the vector which has the real spinvalues for their orientation
	// spin_corr_vec = set_spincorrelation_vector();
  spinstruct_arr = set_spinarray();
}

static inline void propose_spinvals(bool even) {
  /*
    functin that sets the proposed spinvalues to a new random value
      for all even or all odd spins in spinstruct_arr
  */
  int startidx, maxidx, D = get_D(), N = get_N();

  if (even==true) { startidx = 0; maxidx = (pow(N,D))/2; }
  else { startidx = (pow(N,D))/2; maxidx = pow(N,D); }

  for (int idx=startidx; idx<maxidx; idx++) {
    spinstruct_arr[idx].ppspinval = random_spin_orientation();
  }
}

// rewrite this function so it works for proposed spin and current spin
// static double one_spin_hamiltonian(int x) {
static double one_spin_hamiltonian(int x, bool pp_spin ) { // x is the index of the spin for which you are calculating the hamiltonian
/* function which returns the hamiltonian for a single spin with index x.*/
  int D;
  double B;
  D = get_D();
  B = get_B();
  int prod_neighbor_orientation, contrib_neighbor_index, xspinval, nnspinval, spinmodel;
  spinmodel=get_spinmodel();
  double one_spin_hamiltonian = 0.0;
  if (pp_spin==true) { xspinval = spinstruct_arr[x].ppspinval; }
  else { xspinval = spinstruct_arr[x].spinval; }

  for(int i=0; i<2*D; i++)	{
	      // getting the index of the i-th neighbor.
	     contrib_neighbor_index = spinstruct_arr[x].nnidx[i];
       if (contrib_neighbor_index == NOSITE) continue; //skip all next neighbors that are out of bounds (dbc)
	     // multiplying with said neighbor and saving the orientation that characterizes the product.
      // if (pp_spin==true) { nnspinval = spinstruct_arr[contrib_neighbor_index].ppspinval; }
      // else { nnspinval = spinstruct_arr[contrib_neighbor_index].spinval; }
      nnspinval = spinstruct_arr[contrib_neighbor_index].spinval; // current spinvalues of next neighbors is always the background!

	    prod_neighbor_orientation = spinmultiplication(xspinval, nnspinval);

      if (spinmodel==0) { one_spin_hamiltonian -= prod_neighbor_orientation;} //for Ising the orientation is the same as the real value
	    else if (spinmodel==1) {
        one_spin_hamiltonian -= spinval_values[prod_neighbor_orientation];
      }
  }
  one_spin_hamiltonian *= 0.5; //normalization
  // magnetic field contribution
  // magnetic field contribution. For the Clock value the real value of the spin for the orientation has to be taken from spinval_values.
  if (spinmodel==0) { one_spin_hamiltonian -=  xspinval * B; } // s*B which for Ising is simply the orientation*B
  else if (spinmodel==1) {
    one_spin_hamiltonian -=  spinval_values[xspinval] * B;
  }
  else { printf("[MCMC.c | one_spin_hamiltonian()] Invalid spinmodel chosen. (Only 0 and 1 valid).\n"); exit(-1); }

  return one_spin_hamiltonian;
}

static inline double min(double a, double b) {
  return ((a)<(b) ? (a):(b));
}

static void MARS(bool even) {
  /*
    Function makes accept/reject step for either all even or all odd spin values:
    takes decision about even or odd values and may change spinstruct_arr[even/odd].spinval or keep old ones
    1. decides on spinarray index range: 0-len(spinarry)/2 = even , len(spinarry)/2-len(spinarray) = odd
    2. calculates the sum of hamiltonians for this range for current (n or s) and proposed (p or pp) spinconfiguration
      -  for all x in range: h_n,p(x) = -2Re[s_n,p(x)*b(x)]-B*Re[s_n,p(x)]
    3. calculates accept value for ergodicity
    4. update decision:
      - if h_n(x)>=h_x(x): s_p(x)
      - else if random(x)<val_accept(x): s_p(x)
      - else: s_n(x)
  */
  // IF IT SHOULD BE WRITTEN FOR BOTH ORDERINGS -> loop has to be different (but idea with using all even and all odd updating is only good when this ordering is used)
  // double hsum_pp=0.0, hsum_s=0.0, val_accept; // sum of hamiltonian values for proposed and current spins value, accept value (q_s,n^a in desc.)
  double pp=0.0, s=0.0, val_accept;
  double T = get_T(); // how does random(x) work?
  int startidx, maxidx, D = get_D(), N = get_N();

  printf("MCMC running...\n");

  if (even==true) { startidx = 0; maxidx = (pow(N,D))/2; }
  else { startidx = (pow(N,D))/2; maxidx = pow(N,D); }

  for (int idx=startidx; idx<maxidx; idx++) {
    // hsum_pp += one_spin_hamiltonian(pp_arr[idx]);
    pp = one_spin_hamiltonian(spinstruct_arr[idx].idx, true); // if ppspinval in spinstruct_t
    s = one_spin_hamiltonian(spinstruct_arr[idx].idx, false); // 'spinstruct_arr[idx].idx' should be the same as just 'idx'

    val_accept = min( 1, exp((-1)*(double)(pp-s)/(double)T) );

    // if (s >= pp) { spinstruct_arr[idx].spinval = spinstruct_arr[idx].ppspinval; } // else if is same as (s >= pp || random[idx] < val_accept)
    // else if (random[idx] < val_accept) { spinstruct_arr[idx].spinval = spinstruct_arr[idx].ppspinval; }
    if (s >= pp || ((double) rand())/RAND_MAX < val_accept) { // random(x) -> make new random number in each loop step
      spinstruct_arr[idx].spinval = spinstruct_arr[idx].ppspinval;
      printf("Spin flip accepted...\n");
    }
    else { printf("Spin flip denied...\n"); }
  }
}

// for (int idx=startidx; idx<maxidx; idx++) {
//   // hsum_pp += one_spin_hamiltonian(pp_arr[idx]);
//   hsum_pp += one_spin_hamiltonian(spinstruct_arr[idx].idx, true); // if ppspinval in spinstruct_t
//   hsum_s += one_spin_hamiltonian(spinstruct_arr[idx].idx, false); // 'spinstruct_arr[idx].idx' should be the same as just 'idx'
// }
//
// val_accept = min( 1, exp((-1)*(double)(hsum_pp-hsum_s)/(double)T) );
//
// for (int idx=startidx; idx<maxidx; idx++) {
//   // if (hsum_s >= usum_pp) { spinstruct_arr[idx].spinval = spinstruct_arr[idx].ppspinval; } // else if is same as (hsum_s >= usum_pp || random[idx] < val_accept)
//   // else if (random[idx] < val_accept) { spinstruct_arr[idx].spinval = spinstruct_arr[idx].ppspinval; }
//   if (hsum_s >= hsum_pp || ((double) rand())/RAND_MAX < val_accept) { // random(x) -> make new random number in each loop step
//     spinstruct_arr[idx].spinval = spinstruct_arr[idx].ppspinval;
//   }
// }

void MCMC_step(void) {
  /*
    function that performs one Markov Chain Monte Carlo updating step
    updates values inside spinstruct_arr -> doesnt need input nor returns output
    1. loop over all even spins:
        - update all even spins
    2. MARS (accept/reject) step
    3. loop over all odd spins:
        - update all odd spins
        - background can be calculated with already changed spins (but only if MARS was performed on them)
    4. MARS step
  */
  int ordering = get_ordering();
  if (ordering != 1) {
    printf("[MCMC.c | MCMC_step()] ERROR. Actual MCMC updating only works with black white (odd/even) ordering.\n"
    "Choose 1 as 4th input variable (file name counting as 0)\n");
    exit(-1);
  }
  propose_spinvals(true);
  MARS(true);
  propose_spinvals(false);
  MARS(false);
}
