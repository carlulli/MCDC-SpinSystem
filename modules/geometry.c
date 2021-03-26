/*******************************************************************************
This module handles the transformation of the real lattice to the discretized
version used for computing.
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "spin.h"

static int D=1; // D = Dimension
static int N=1; // N = Number of lattice sites in one direction (per dimension)
static int lattice_spacing=1; // lattice spacing fixed to 1 in project!
static int boundary_condition;
static int *bc_ptr=NULL;


void set_params(int argc, char const *argv[]) {
  // if argc not some value -> error message and exit
  if (N%2 != 0) {
    printf("[geomerty.c | set_params()] ARRR. Number of lattice points in one direction needs to be EVEN!\n");
    exit(-1);
  }
  if (D < 1) {
    printf("[geometry.c | set_params()] You need at least 1 dimension, smh.\n");
    exit(-1);
  }
  if (boundary_condition > 1) {
    printf("[geometry.c | set_params()] Boundary Condition out of bounds. Should be 0 for Dirichlet and 1 for periodic. \n");
    exit(-1);
  }
  N = atoi(argv[1]);
  D = atoi(argv[2]);
  if ()
  boundary_condition = atoi(argv[3]);
  bc_ptr = &boundary_condition;
}


/*******************************************************************************
Spin Struct (is this something for the header?)
*******************************************************************************/
struct spinstruct {
  /*
  structure that contains all values belonging to one lattive site:
  own index, spin value at this lattice site, own coordinates (array/list?), nn index array which is only 2D
  */
  int idx;
  double spinval; // maybe just int
  int *coord;
  int *nnidx; // doesnt need to be 2d as i am using array of structs
}

inline void spinstruct_alloc(spinstruct spstrct) {
  //allocate only for arrays?
  // safety feature to not allocate twice (initialzed pointers in struct are NULL )
  spstrct.coord = (int *) malloc(sizeof(int)*D);
  // spstrct.nnidx = (int *) malloc(sizeof(int)*pow(N,D)*2*D); // allocate 2d array!!! see above
  spstrct.nnidx = (int *) malloc(sizeof(int)*2*D);
}

inline void spinstruct_free(spinstruct spnstrct) {
  /* function to deallocate (free) memory used for spinstruct spnstrct */
  free(spstrct.coord);
  free(spstrct.nnidx);
}


/*******************************************************************************
Functions for distance calculation
*******************************************************************************/
inline double coord_distance_dirichlet(spinstruct *spnstrct1, spinstruct *spnstrct2) {
  /*
    function that caluates distance btw 2 spins for dirichlet bc
    input are spinstructs -> coordinates already known (using vector substr.)
  */
  double sum;
  for (int i=0; i<D; i++) {
    sum += (spnstrct1.coord[i] - spnstrct2.coord[i])*(spnstrct1.coord[i] - spnstrct2.coord[i]);
  }
  return sqrt(sum);
}

inline double coord_distance_periodic(spinstruct *spnstrct1, spinstruct *spnstrct2) {
  /*
    funtion that calclulates distacnce btw 2 spins for periodic bc
    input are spinstructs -> coordinates already known
    algo:
      1. dx = |x1_a - x1_b| for all xi in array[D]
      2. if dx > N/2: reduce dx for N -> dx = dx - N
      3. sum =  dx1^2 + dx2^2 + ... ->  sum can be calculated in same loop
      4. return sqrt(sum) (pythagoras)
  */
  double sum;
  int deltaxi[D];
  for (int i=0; i<D; i++) {
    deltaxi[i] = abs(spnstrct1.coord[i] - spnstrct2.coord[i]);
    if (deltaxi[i] > N/2) { deltaxi -= N; }
    sum += dx*dx;
  }
  return sqrt(sum);
}

double coord_distance(spinstruct *spnstrct1, spinstruct *spnstrct2) {
  /*
    wrapper function to calculate the distance between two indeces
    needs some way to respect boundary conditions
    returns the distance as double
  */
  double dist;
  if (bc_ptr=NULL) {
    printf("[geometry.c | coord_distance_blackwhite()] ERROR. Boundary Condition was not set!\n");
    exit(-1);
  }
  if (boundary_condition==0) { dist = coord_distance_dirichlet_blackwhite(spnstrct1, spnstrct2); }
  else if (boundary_condition==1) { dist = coord_distance_periodic_blackwhite(spnstrct1, spnstrct2); }
  else {
    printf("[geometry.c | coord_distance_blackwhite()] Upsi something else went wrong with the bc but it got set with the right value.\n");
    exit(-1);
  }
  return dist;
}

/*******************************************************************************
BLACK/WHITE ORDERING
*******************************************************************************/
void set_spinarray_blackwhite() {
  // GOAL: assigns elements (spin values) to an array in a way that it contains all spin values for the whole lattice
  //        this array can be called with the correct index and returns the element at this index
  //          that is: the spin value of spin at this lattice site where index points to
  // (assigns spin values to)? an array with all even lattice points first (black)
  //  and all white lattice points second (white)
  // above is the coordinate ordering, ie. first half indeces are even lattice points (black)
  // second half are odd indeces (white) (even and odd refer to the parity of coordinate tuple)
  // array is N^d long
  // array elements are spin orientation spin at this lattice point
  // -> function needs acces to spin module
  // does it even need acces to the function how the lattice points are created
  // partiy p(X) = (sum_mu=0^d-1 x_mu) % 2
  // an array has to exist where the needed and elements can be taken by using the index -> STRUCT idea?
  // the value at one index is not just the spin, but a struct
  // the struct contains: spin value, coordinates, array with nn indexes, own index
  /*
  1. set_params to get N and d
  2. allocate arrray of structs with size N^d
  3. loop through array and assign each point:
    -  the corresponding coord.
    -  a spin value
    - its index
    - next neighbor indeces (maybe by calling the corresponding function)
  */
  spinstruct *spinstruct_arr;
  set_params();
  spinstruct_arr = (spinstruct) malloc(sizeof(spinstruct)*pow(N, D));

  for (int i=0; i<pow(N, D), i++){
    spinstruct_alloc(spinstruct_arr[i]);
    spinstruct_arr[i].idx = i;
    spinstruct_arr[i].spinval = function_that_assigns_spinvalue();
    set_coord_blackwhite(i, spinstruct_arr[i].coord);
    set_nnarray_blackwhite(spinstruct_arr[i]);
  }
}

inline int parity(int *x) { // static inline if only needed in this file -> make declaration of all function at the top of file (with small comment and maybe even line)
  int sum;
  for (xi=0; xi<D; xi++) {
    sum += x[xi];
  }
  return sum%2;
}

inline int np_parity(int np) {
  return ((int) 2*np) / ((int) pow(N,D));
}

inline int iszero(int x) { // (x==0)
  int outcome;
  if (x==0) { outcome=1; }
  else { outcome=0; }
  return outcome;
}

int np_of_x(int *x) { // also inline?
  /* x = coordinate vector and np is index in black white ordering*/
  int np;
  np = (int)(parity(x)*pow(N,D)+n_of_x(x))/2;
  return np;
}

// maybe not set but "idx2coord" or something, as it returns something?
//inline void set_coord_blackwhite(); ??
void set_coord_blackwhite(int np, int *x) {
  /* np = current index, x = coordinate vector */
  int xprime[D]; // dynamic allocation better?
  for (int mu=0; mu<D; mu++) {
    xprime[mu] = ( 2*np-(2*np/pow(N,D))*pow(N,D)/pow(N,mu) )%N;
  }
  for (int xi=0; xi<D; xi++) {
    if (xi==0) { x[xi] = ( xprime[xi]+(np_parity(np)+parity(xprime))%2 )%N; }
    else { x[xi] = ( xprime[xi]+(np_parity(np)+parity(xprime))%2 * iszero(xprime[mu-1]) )%N; }
  }
}

//inline void set_nnarray_blackwhite(); ??
void set_nnarray_blackwhite(spinstruct *spnstrct) { // für lexo und blackwhite nutzen, dann mit wrapper funtion für np_of_x und n_of_x
    // WHAT ABOUT BOUNDARY CONDITION?
  /*
    function that sets the next neigbor indeces into a given spinstruct
    index array ordering:
      nnidx = [0, 1, 2, 3, ..., 2D-1, 2D] = [x0+1, x0-1, x1+1, x1-1,..., 2D+1, 2D-1]
    algo:
      1. copy coordinates to new array
      2. create new corrd array (x0+1, x0) & get index to assign to nnidx
      2. create new corrd array (x0-1, x0) & get index to assign to nnidx
      3. repeat in steps of 2 for i < D
  */
  int nnx[d];
  for (int j=0; j<D; j++) { nnx[j] = spnstrct.coord[j]; }

  for (int i=0; i<D; i+=2) {
    nnx[i] = spnstrct.coord[i]+1; // could also be nnx[i]+1 (would be the same)
    spnstrct.nnidx[i] = np_of_x(nnx);
    nnx[i] = spnstrct.coord[i]-1;
    spnstrct.nnidx[i+1] = np_of_x(nnx);
    nnx[i] = spnstrct.coord[i];
    // falls nnx[i] = N,-1 -> fallunterscheidung ob dirichlet / periodic
    // if dirichlet: einen festen wert, der nicht im gitter für spätere Fallunterscheidung (Variable: nosite = zb -1) -> wenn nn später benutzt mache if condition für nosite, im fall von nosite: kein nachbar element
    // if peridic: (coord+N)%N to reappear on coordinate of other side of lattice
  }
}

/*******************************************************************************
LEXOGRAPHIC ORDERING
*******************************************************************************/
void set_spinarray_lexo() {
  // GOAL: assigns elements (spin values) to an array in a way that it contains all spin values for the whole lattice
  //        this array can be called with the correct index and returns the element at this index
  //          that is: the spin value of spin at this lattice site where index points to
  // ordering of lattice sites to index is lexographic -> counted through in one direction
  // array is N^d long
  // array elements are spinstructs with spinvalue, own index, next neighbor index & own coordinates
  /*
  1. set_params to get N and d
  2. allocate arrray of structs with size N^d
  3. loop through array and assign each point:
    -  the corresponding coord.
    -  a spin value
    - its index
    - next neighbor indeces (maybe by calling the corresponding function)
  */
  spinstruct *spinstruct_arr;
  set_params();
  spinstruct_arr = (spinstruct) malloc(sizeof(spinstruct)*pow(N, D));

  for (int i=0; i<pow(N, D), i++){
    spinstruct_alloc(spinstruct_arr[i]);
    spinstruct_arr[i].idx = i;
    spinstruct_arr[i].spinval = function_that_assigns_spinvalue();
    set_coord_lexo(i, spinstruct_arr[i].coord);
    set_nnarray_lexo(spinstruct_arr[i]);
  }
}

void set_coord_lexo(int n, int *x) {
  /* n = current index, x = coordinate vector */
  for (int mu=0; mu<D; mu++) {
    x[mu] = (n/pow(N,mu))%N;
  }
}

int n_of_x(int *x) {
  /* x is coordinate vector and n is index in lexo ordering  */
  int sum;
  for (mu=0; muy<D; mu++) { sum = (x[mu]*pow(N, mu))/lattice_spacing; }
  return sum;
}

void set_nnarray_lexo(spinstruct *spnstrct) {
  // WHAT ABOUT BOUNDARY CONDITION
  /*
    function that sets the next neigbor indeces into a given spinstruct
    index array ordering:
      nnidx = [0, 1, 2, 3, ..., 2D-1, 2D] = [x0+1, x0-1, x1+1, x1-1,..., 2D+1, 2D-1]
      algo:
      1. nnidx[0] = nnidx+N^0 & nnidx[1] = nnidx[1] = nnidx-N^0
      2. keep going for power: mu<D and indexing in steps of 2 and smaller 2D/2=D
  */
  for (int mu=0; m<D; mu++) {
    for (int i=0; i<D; i+=2) {
      spnstrct.nnidx[i] += pow(N,mu);
      spnstrct.nnidx[i+1] -= pow(N,mu);
    }
  }
}
