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

/***************************************************+***************************
Functions:
x -> index because of ordering: spin value at this lattice site = spinstruct_arr[idx]->spinval

+ BACKGROUND(x) (calculates background for one spin x )
  - b(x)=sum_(mu<D) s(x+e_mu)+(s(x-e_mu))
  - we need background for all updated spins
    -> should be able to calc through loop, as background of even spin does only depend on odd (and vice versa)
  - just loop over all next neighbors?
  - b(x)b(idx)=for (i<D) b += spinstruct_arr[spinstruct_arr[idx].nnidx[i]].spinval

+ HAMILTONIAN(x)_p,n (for one x)
  - -2Re[s_n,p(x)*b(x)]-B*Re[s_n,p(x)]
  - is this the same hamiltonian as in the hamiltonian module?
    -> should it be?
  - write in hamiltonian module?
  -

+ MARS accept/reject step for a number of flipped spins, but works same as with one:
    - input values:
            * proposed spin configuration (pointer to array with proposed values)
            * old spin value configuration (existing global variable spinstruct_arr)
            * sum of proposed spin value hamiltonians (for hamiltonians, background needs to be calculated)
            * sum of old spin value hamiltonians
    - is it model independent in general or does it need a wrapper function?
    - has if & if else & else condition with sum of spins and sum of hamiltonians
    - calculates value for which to accept: q(x) = min(1, e^{-1/kT(h_p(x)- h_n(x)))
    - what does it do: change or doesnt change the spinvalue (spinstruct_arr[i].spinval) for all even/odd spins

+ UPDATE loop through one set of indeces (even or odd)
    - for all even or all odd loop:
      - update per spin is performed by spin module

+ MCMC_step function that does one MCMC update step:
    1. loop over all even spins:
        - update all even spins
    2. MARS (accept/reject) step
    3. loop over all odd spins:
        - update all odd spins
        - background can be calculated with already changed spins (but only if MARS was performed on them)
    4. MARS step
***************************************************+***************************/
#include <math.h>
#include "geometry.h"
#include "spin.h"

static inline void propose_spinvals(bool even=TRUE);
static double one_spin_hamiltonian(int x);
static inline double min(a, b); // logically better in utilities, but since used here -> simpler
static void MARS(int prop_arr, bool even=TRUE);



static inline void propose_spinvals(bool even=TRUE) {
  /*
    functin that sets the proposed spinvalues to a new random value
      for all even or all odd spins in spinstruct_arr
  */
  int idx, maxidx, D = get_D(), N = get_N();

  if (even==TRUE) { idx = 0; maxidx = (pow(N,D)/2; }
  else { idx = (pow(N,D)/2; maxidx = pow(N,D); }

  for (idx; idx<maxidx; idx++) {
    spinstruct_arr[idx].ppspinval = random_spinval();
  }
}

// rewrtie this function so it works for proposed spin and curent spin
// static double one_spin_hamiltonian(int x) {
static double one_spin_hamiltonian(int x, bool pp_spin=TRUE ) { // x is the index of the spin for which you are calculating the hamiltonian
/* function which returns the hamiltonian for a single spin with index x.*/
    D = get_D();
    B = get_B();
    int prod_neighbor_orientation, contrib_neighbor_index, xspinval, nnspinval;
    double one_spin_hamiltonian = 0.0;
    if (pp_spin==TRUE) { xspinval = spinstruct_arr[x].ppspinval; }
    else { xspinval = spinstruct_arr[x].spinval; }

    for(i=0, i<2*D,i++)	{
	       // getting the index of the i-th neighbor.
	      contrib_neighbor_index = spinstruct_arr[x].nnidx[i];
	      // multiplying with said neighbor and saving the orientation that characterizes the product.
        if (pp_spin==TRUE) { nnspinval = spinstruct_arr[contrib_neighbor_index].ppspinval); }
        else { nnspinval = spinstruct_arr[contrib_neighbor_index].spinval); }

	      prod_neighbor_orientation = spinmultiplication(xspinval, nnspinval);
        if (spinmodel==0) { one_spin_hamiltonian -= 2* prod_neighbor_orientation;} //for Ising the orientation is the same as the real value
	      else if (spinmodel==1) { one_spin_hamiltonian -= 2* spinval_values[prod_neighbor_orientation];} // multiplied times two, to account for interaction going both ways.
	      // multiplied times two, to account for interaction going both ways
    }
    // magnetic field contribution
    // magnetic field contribution. For the Clock value the real value of the spin for the orientation has to be taken from spinval_values.
    if (spinmodel==0) { one_spin_hamiltonian -=  xspinval * B; } // s*B which for Ising is simply the orientation*B
    else if (spinmodel==1) { one_spin_hamiltonian -=  spinval_values[xspinval] * B; }

    return one_spin_hamiltonian;
}

static inline double min(a, b) {
  return ((a)<(b) ? (a):(b));
}

static void MARS(bool even=TRUE) {
  /*
    Function makes accept/reject step for either all even or all odd spin values:
    takes decision about even or odd values and may change spinstruct_arr[even/odd].spinval or keep old ones
    1. decides on spinarray index range (0-len(spinarry)/2 = even & len(spinarry)/2-len(spinarray) = odd)
    2. calculates the sum of hamiltonians for this range for current (n or s) and proposed (p or pp) spinconfiguration
      -  for all x in range: h_n,p(x) = -2Re[s_n,p(x)*b(x)]-B*Re[s_n,p(x)]
    3. calculates accept value for ergodicity
    4. update decision:
      - if h_n(x)>=h_x(x): s_p(x)
      - else if random(x)<val_accept(x): s_p(x)
      - else: s_n(x)
  */
  // IF IT SHOULD BE WRITTEN FOR BOTH ORDERINGS -> loop has to be different (but idea with using all even and all odd updating is only good when this ordering is used)
  double hsum_pp, hsum_s, val_accept; // sum of hamiltonian values for proposed and current spins value, accept value (q_s,n^a in desc.)
  double T = get_T(), random; // how does random(x) work?
  int idx, maxidx, D = get_D(), N = get_N();

  if (even==TRUE) { idx = 0; maxidx = (pow(N,D)/2; }
  else { idx = (pow(N,D)/2; maxidx = pow(N,D); }

  for (idx; idx<maxidx; idx++) {
    // hsum_pp += one_spin_hamiltonian(pp_arr[idx]);
    hsum_pp += one_spin_hamiltonian(spinstruct_arr[idx].idx, TRUE); // if ppspinval in spinstruct_t
    hsum_s += one_spin_hamiltonian(spinstruct_arr[idx].idx, FALSE);
  }

  val_accept = min( 1, exp((-1)*(double)(hsum_pp-hsum_s)/(double)T) );

  for (idx; idx<maxidx; idx++) {
    // if (hsum_s >= usum_pp) { spinstruct_arr[idx].spinval = spinstruct_arr[idx].ppspinval; } // else if is same as (hsum_s >= usum_pp || random[idx] < val_accept)
    // else if (random[idx] < val_accept) { spinstruct_arr[idx].spinval = spinstruct_arr[idx].ppspinval; }
    if (hsum_s >= usum_pp || ((double) rand())/RAND_MAX < val_accept) { // random(x) -> make new random number in each loop step
      spinstruct_arr[idx].spinval = spinstruct_arr[idx].ppspinval;
    }
  }
}

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
  propose_spinvals(TRUE);
  MARS(TRUE);
  propose_spinvals(FALSE);
  MARS(FALSE);
}
