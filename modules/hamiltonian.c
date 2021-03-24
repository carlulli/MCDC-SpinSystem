/*******************************************************************************
This module handles the hamiltonian for a ferromagnetic spin system in multiple
dimensions
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>


void H(double complex *in, double complex *out) {
/*
function to calculate the hamiltonian applied to a wavefunction
H(s) = - sum_[x,y]nextneighbors spinproduct(x,y) - B * sum_x s(x)
s(x) is spin at point x
x = (x1, x2, ... , xd) with d = Dimension
spinproduct = s(x)s(y) depends on indexing and boundary conditions

sum over next neighbors: nn should be an array that iterates through all next
  neighbors (size depends on d, e.g. 2D: 4 nn)
  this array should be given by geometry module

s(x) is a function given by updating method ?

B is given to simulation as parameter

 */
}
