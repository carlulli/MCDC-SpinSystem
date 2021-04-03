/*******************************************************************************
Header file of geometry module.
Here only global functions are listed.
At module's beginning, a detailed list of functions is written.
*******************************************************************************/
#ifndef GEOMETRY_H
#define GEOMETRY_H

typedef struct spinstruct { // can also have typedef in front
  /*
  structure that contains all values belonging to one lattive site:
  own index, spin value at this lattice site, some proposed spin value, own coordinates (array/list?), nn index array which is only 2D
  */
  int idx;
  int spinval; // maybe just int
  int ppspinval;
  int *coord;
  int *nnidx; // doesnt need to be 2d as i am using array of structs
} spinstruct_t; // if typedef at beginning -> use struct_tag here: spinstruct_t

// extern spinstruct_t *spinstruct_arr;

void set_params(int argc, char const *argv[]);
/* sets all nevessary parameters to create geometry */
int get_N();
/* returns N */
int get_D();
/* returns D */
int get_ordering();
/* returns ordering value (so far without safety function) */
int get_boundary_condition();
/* returns boundary condition value (so far without safety function) */
void spinarray_free(spinstruct_t *spnstrct_arr);
/* frees memory for array of struct (of size N^D) (wrapper of 3.) */
double coord_distance(spinstruct_t *spnstrct1, spinstruct_t *spnstrct2);
/* calculates distance btw 2 spinindeces (wrapper for 5. & 6.) */
// void set_spinarray(spinstruct_t *spinstruct_arr);
// void set_spinarray(void);
spinstruct_t* set_spinarray(void);
/* sets spinstructs for given array and ordering -> wrapper function for 8. & 13. */
int get_arraysize();
/* return size of allocated systems total spin array */
int n_of_x(int *x);


#endif /* GEOMETRY_H */
