/*******************************************************************************
This module handles the transformation of the real lattice to the discretized
version used for computing.
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>


void set_geometry() {
  // sets all geometry values needed for the other functions
}

void set_lexo_index() {
  //creates an array with the lattice point indices in lexographical order
  // n(X) = sum_mu=0^d-1 x_mu * N^mu * 1/a
}

void set_blackwhite_index() {
  // creates an array with all even lattice points first (black)
  //  and all white lattice points second (white)
  // np(X) = 1/2 * ( p(X) * N^d + (sum_mu=0^d-1 x_mu * N^mu * 1/a) )
  // partiy p(X) = (sum_mu=0^d-1 x_mu) % 2
}

int get_lexo_coordinate(int n) {
  // function returns specific coordinate for asked mu and n
  // x_mu(n) = (n/N^mu) % N with integer devision!
}

int get_blackwhite_coordinate(int np) {
  // function returns specific coordinate for asked mu and np
  // x_mu(np) = ( (2*np + (2*np/N^d))/N^mu ) % N
}

int set_lexo_nextneighbor(int n, decision value for +1 or -1 neighbor ) {
  // function to return the next neighbor of a given index
  // should be able to return next (+1) or previous (-1) neighbor
  // should this be a seperate function? should this take pointers because of array
}

int set_blackwhite_nextneighbor(int np, decision value for +1 or -1 neighbor ) {
  // function to return the next neighbor of a given index
  // should be able to return next (+1) or previous (-1) neighbor
  // should this be a seperate function? should this take pointers because of array
}

// something to create necessary arrays? 
