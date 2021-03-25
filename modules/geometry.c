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


void set_params(int argc, char const *argv[]) {
  // if argc not some value -> error message and exit
  if (N%2 != 0) {
    printf("[geomerty.c | set_params()] ARRR. Number of lattice points in one direction needs to be EVEN!\n");
    exit(-1);
  }
  N = atoi(argv[1]);
  d = atoi(argv[2]);
}

struct spinstruct {
  /*
  structure that contains all values belonging to one lattive site:
  own index, spin value at this lattice site, own coordinates (array/list?), nn index array which is only 2D
  */
  int idx;
  double spinval;
  double coord;
  double nnidx; // doesnt need to be 2d as i am using array of structs
}

inline void spinstruct_alloc(spinstruct spstrct) {
  //allocate only for arrays?
  spstrct.coord = (double *) malloc(sizeof(double)*D);
  // spstrct.nnidx = (int *) malloc(sizeof(int)*pow(N,D)*2*D); // allocate 2d array!!! see above
  spstrct.nnidx = (int *) malloc(2*D);
}

inline void spinstruct_free(spinstruct spnstrct) {
  /* function to deallocate (free) memory used for spinstruct spnstrct */
  free(spstrct.coord);
  free(spstrct.nnidx);
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

  for (i=0; i<pow(N, D), i++){
    spinstruct_alloc(spinstruct_arr[i]);
    spinstruct_arr[i].idx = i;
    spinstruct_arr[i].spinval = function_that_assigns_spinvalue();
    set_coord_blackwhite(i, spinstruct_arr[i].coord);
    set_nnarray_blackwhite(spinstruct_arr[i]);
  }
}

inline int parity(int *x) {
  int sum;
  for (xi=0; xi<D; xi++) {
    sum += x[xi];
  }
  return sum%2;
}

inline int np_parity(int np) {
  return ((int) 2*np) / ((int) pow(N,D));
}

inline int isequal(int x) {
  int outcome;
  if (x==0) { outcome=1; }
  else { outcome=0; }
  return outcome;
}

// maybe not set but "idx2coord" or something, as it returns something?
//inline void set_coord_blackwhite(); ??
void set_coord_blackwhite(int np, int *x) {
  /* np = current index, x = coordinate vector */
  int xprime[D]; // dynamic allocation better?
  for (mu=0; mu<D; mu++) {
    xprime[mu] = ( 2*np-(2*np/pow(N,D))*pow(N,D)/pow(N,mu) )%N;
  }
  for (xi=0; xi<D; xi++) {
    if (xi==0) { x[xi] = ( xprime[xi]+(np_parity(np)+parity(xprime))%2 )%N; }
    else { x[xi] = ( xprime[xi]+(np_parity(np)+parity(xprime))%2 * isequal(xprime[mu-1]) )%N; }
  }
}

//inline void set_nnarray_blackwhite(); ??
void set_nnarray_blackwhite(spinstruct *spnstrct) {
  // GOAL: creates a 2d array: first dimension are lattice point indeces and second dimension are the next neighbors indeces
  //        when summing over nn, one can call this array to recover the lattice points next neighbor indeces (of the spin array)
  // WHY do i need 2 dimensions? As i wrote it above, the indeces would already work??? -> this basically are the two dimensions
  // nn[N^d][2d] two indexed array e.g. nn[4][xforward] = a value that i have to fill
  // the value should be the next neighbor index (depending on the indexing method)
  // PARAMTERS: rnng_idx is current index the spinstruct that is also a paramter;
  //            spnstrct is current lattice index' spinstruct & spnstrct.nnidx is 2*d array
  int nnx[d];
  for (j=0; j<D; j++) { nnx[j] = spnstrct.idx[j]; }

  for (i=0; i<D; i++) {
    nnx[i] = spnstrct.idx[i]+1; // could also be nnx[i]+1 (would be the same)
    spnstrct.nnidx[i] = np_of_x(nnx);
    nnx[i] = spnstrct.idx[i]-1;
    spnstrct.nnidx[i+1] = np_of_x(nnx);
    nnx[i] = spnstrct.idx[i];
  }
}

int np_of_x(int *x) {
  /* x = coordinate vector and np is index in black white ordering*/
  int np;
  np = (int)(parity(x)*pow(N,D)+n_of_x(x))/2;
  return np;
}

double coord_distance_blackwhite() {
  /*
    function to calculate the distance between two indeces
    needs some way to respect boundary conditions
  */
}

void get_coord_blackwhite(int n, int *coord) {
  // function returns specific coordinate for asked mu? and np
  // returns coordinate tuple (list, or array,...?)
  // tuple is d dimensional
  for (i=0; i<D; i++) {
    coord[i] =
  }
}

int get_nncoord_blackwhite(int n) {
  // funtion returns index of next neighbor of input index
  // returns array of the neighbor indeces
}


/*******************************************************************************
LEXOGRAPHIC ORDERING
*******************************************************************************/
void set_lexo_index() { // DOES IT?
  //creates an array with the lattice point indices in lexographical order
  // n(X) = sum_mu=0^d-1 x_mu * N^mu * 1/a
  double complex *n, **nn;
  n = (double complex*) malloc(dtype(double complex)*pow(N,d)); // double complex surely wrong
  nn = (double complex*) malloc(dtype(double complex)*pow(N,d)*(2*d));

  for (int i=0; i<N; i++) {
    n[i] = i;

    /*
      specific for 2 dim
      [0,1,2,3] = [xnext, xprevious, ynext/up, yprev/down]
      what about the boundary conditions??
    */
    n[i][0] = i+1;
    n[i][1] = i-1;
    n[i][2] = i+N;
    n[i][3] = i-N;

    /* for all dimensions d
      [0,1,2,3,4,5,...] = [xnext, xprevious, ynext/up, yprev/down, znext/forward, zprev/behind,...]
      nn[i][j] = i + (-1)^j * N^(j/2) with integer devision
    */
    for (int j=0; j<2*d; j++) {
      nn[i][j] = i + pow(-1, j) * pow(N, (j/2));
    }
  }
}

int n_of_x(int *x) {
  /* x is coordinate vector and n is index in lexo ordering  */
  int sum;
  for (mu=0; muy<D; mu++) { sum = (x[mu]*pow(N, mu))/lattice_spacing; }
  return sum;
}

int get_lexo_coordinate(int n) {
  // function returns specific coordinate for asked mu and n
  // x_mu(n) = (n/N^mu) % N with integer devision!
}

void coord_from_lexo(int n, int x[]) {
  // loop that fills in coordinate array
  // x array is d dimensional
}

int set_lexo_nextneighbor(int n, decision value for +1 or -1 neighbor ) {
  // function to return the next neighbor of a given index
  // should be able to return next (+1) or previous (-1) neighbor
  // should this be a seperate function? should this take pointers because of array
}


// something to create necessary arrays?
