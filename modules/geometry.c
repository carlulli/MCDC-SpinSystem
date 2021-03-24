/*******************************************************************************
This module handles the transformation of the real lattice to the discretized
version used for computing.
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "spin.h"

static int d=1;
static int N=1;


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
  own index, spin value at this lattice site, own coordinates (array/list?), nn index array
  */
  int idx;
  double spinval;
  double coord;
  double nnidx;
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
}

void set_nnarray_blackwhite() {
  // GOAL: creates a 2d array: first dimension are lattice point indeces and second dimension are the next neighbors indeces
  //        when summing over nn, one can call this array to recover the lattice points next neighbor indeces (of the spin array)
  // WHY do i need 2 dimensions? As i wrote it above, the indeces would already work??? -> this basically are the two dimensions
  // nn[N^d][2d] two indexed array e.g. nn[4][xforward] = a value that i have to fill
  // the value should be the next neighbor index (depending on the indexing method)
}

int get_coord_blackwhite(int n) {
  // function returns specific coordinate for asked mu? and np
  // returns coordinate tuple (list, or array,...?)
  // tuple is d dimensional
}

int get_nncoord_blackwhite(int n) {
  // funtion returns index of next neighbor of input index
  // returns array of the neighbor indeces
}

int set_blackwhite_nextneighbor(int np, decision value for +1 or -1 neighbor ) {
  // function to return the next neighbor of a given index
  // should be able to return next (+1) or previous (-1) neighbor
  // should this be a seperate function? should this take pointers because of array
}


/*******************************************************************************
LEXOGRAPHIC ORDERING
*******************************************************************************/
void set_lexo_index() { // DOES IT?
  //creates an array with the lattice point indices in lexographical order
  // n(X) = sum_mu=0^d-1 x_mu * N^mu * 1/a
  double complex *n, **nn;
  n = (double complex*) malloc(dtype(double complex)*pow(N,d));
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
