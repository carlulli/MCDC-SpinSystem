/*******************************************************************************
This module handles the transformation of the real lattice to the discretized
version used for computing.
*******************************************************************************/

/*******************************************************************************
List of functions. (some declared here (S), others declared in header (H))
1. void set_params(argc, argv[]) (H)
    - sets all nevessary parameters to create geometry
2. int get_N()
    - returns N
3. int get_D()
    - returns D
4. static void spinstruct_alloc(spinstruct spstrct) (S)
    - allocates memory for struct arrays coord and nextneighbor_index
5. static inline void spinstruct_free(spinstruct spnstrct) (S)
    - frees for struct allocated memory
6. void spinarray_free(spinstruct *spnstrct_arr) (H)
    - frees memory for array of struct (of size N^D) (wrapper of 3.)
7. static inline double coord_distance_dirichlet(spinstruct *spnstrct1, spinstruct *spnstrct2) (S)
    - calculates distance btw 2 spinindeces for Dirichlet bc
8. static inline double coord_distance_periodic(spinstruct *spnstrct1, spinstruct *spnstrct2) (S)
    - calculates distance btw 2 spinindeces for periodic bc
9. double coord_distance(spinstruct *spnstrct1, spinstruct *spnstrct2) (H)
    - calculates distance btw 2 spinindeces (wrapper for 5. & 6.)
10. static inline int n_of_x_lexo(int *x) (S)
    - calculates index of coordinate for lexo ordering
11. static inline int n_of_x_blackwhite(int *x) (S)
    - calculates index of coordinate for blackwhite ordering
12. static int n_of_x(int *x) (S)
    - returns index of coordinate for both orderings (wrapper of 10. and 11.)
13. static void set_nnarray(spinstruct *spnstrct) (S)
    - sets a next neighbor index array for given struct
14. static void set_spinarray_blackwhite(spinstruct *spinstruct_arr) (S)
    - creates structs and assigns values for given array in blackwhite ordering
15. static inline int parity(int *x) (S)
    - returns parity of coordinate array
16. static inline int np_parity(int np) (S)
    - returns special parity 2np/n^D
17. static void set_coord_blackwhite(int np, int *x) (S)
    - sets the coordinates in coordinate vector for given index for blackwhite ordering
18. static void void set_spinarray_lexo(spinstruct *spinstruct_arr) (S)
    - creates structs and assigns values for given array in lexographic ordering
19. static inline void set_coord_lexo(int n, int *x) (S)
    - sets the coordinates in coordinate vector for given index for lexographic ordering
20. void set_spinarray(spinstruct *spinstruct_arr) (H)
    - sets spinstructs for given array and ordering -> wrapper function for 8. & 13.
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "spin.h"

static int D=1; // D = Dimension
static int N=0; // N = Number of lattice sites in one direction (per dimension)
static int lattice_spacing=1; // lattice spacing fixed to 1 in this project!
static int boundary_condition;
static int *bc_ptr=NULL;
static int ordering=1; // 0 for lexo and 1 for black white
int nosite = -1; // should be unchangeable but open for other modules

/* static function initialition inside module as not needed outside */
static void spinstruct_alloc(spinstruct spstrct);
static inline void spinstruct_free(spinstruct spnstrct);
static inline double coord_distance_dirichlet(spinstruct *spnstrct1, spinstruct *spnstrct2);
static inline double coord_distance_periodic(spinstruct *spnstrct1, spinstruct *spnstrct2);
static int n_of_x_lexo(int *x);
static inline int n_of_x_blackwhite(int *x);
static int n_of_x(int *x);
static void set_spinarray_blackwhite(spinstruct *spinstruct_arr);
static inline int parity(int *x);
static inline int np_parity(int np);
static void set_coord_blackwhite(int np, int *x);
static void set_spinarray_lexo(spinstruct *spinstruct_arr);
static inline void set_coord_lexo(int n, int *x);

/*******************************************************************************
Parameter handling
*******************************************************************************/
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
  ordering = atoi(argv[4]);
}

int get_N() {
    static int sN=0;
    if(N==0) { printf("[ geometry.c| get_N ] Error! N not yet set!\n"); exit(-1); }
    if(sN==0) { sN=N; }
    else {
      if((N!=sN))  {
         printf("[ geometry.c| get_N ] Error! (N) has changed: (%d) -> (%d)\n",sN,N);
         exit(-1);
      }
    }
    return N;
}

int get_D() {
    static int sD=0;
    if(D==0) { printf("[ geometry.c| get_D() ] Error! D not yet set!\n"); exit(-1); }
    if(sD==0) { sN=D; }
    else {
      if((N!=sD))  {
         printf("[ geometry.c| get_N ] Error! (N) has changed: (%d) -> (%d)\n",sD,D);
         exit(-1);
      }
    }
    return D;
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

static void spinstruct_alloc(spinstruct spstrct) {
  //allocate only for arrays?
  // safety feature to not allocate twice (initialzed pointers in struct are NULL )
  if (spstrct.coord==NULL | spstrct.nnidx==NULL) {
    spstrct.coord = (int *) malloc(sizeof(int)*D);
    spstrct.nnidx = (int *) malloc(sizeof(int)*2*D);
  }
  else {
    printf("[geomtry.c | spinstruct_alloc()] ERROR: apparently struct values are already allocated. \n");
    exit(-1);
  }
}

static inline void spinstruct_free(spinstruct spnstrct) {
  /* function to deallocate (free) memory used for spinstruct spnstrct */
  free(spstrct.coord);
  free(spstrct.nnidx);
}

void spinarray_free(spinstruct *spnstrct_arr) {
  for (int i=0; i<pow(N,D); i++) { spinstruct_free(spinstruct_arr[i]); }
}

/*******************************************************************************
Functions for distance calculation
*******************************************************************************/
static inline double coord_distance_dirichlet(spinstruct *spnstrct1, spinstruct *spnstrct2) {
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

static inline double coord_distance_periodic(spinstruct *spnstrct1, spinstruct *spnstrct2) {
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
Next neighbor part (for both orderings)
*******************************************************************************/

static inline int n_of_x_lexo(int *x) {
  /* x is coordinate vector and n is index in lexo ordering  */
  int sum;
  for (mu=0; muy<D; mu++) { sum = (x[mu]*pow(N, mu))/lattice_spacing; }
  return sum;
}

static inline int n_of_x_blackwhite(int *x) { // also inline?
  /* x = coordinate vector and np is index in black white ordering*/
  int np;
  np = (int)(parity(x)*pow(N,D)+n_of_x(x))/2;
  return np;
}

static int n_of_x(int *x) {
  if (ordering==0) { return n_of_x_lexo(x); }
  if (ordering==1) { return n_of_x_blackwhite(x); }
  else {
    printf("[geometry.c | n_of_x()] ERROR. No crrect odering value assined. 0 for lexo and 1 for blackwhite \n");
    exit(-1);
  }
}

static void set_nnarray(spinstruct *spnstrct) {
  /*
    function that sets the next neighbor indeces into a given spinstruct
    works for lexo and black white ordering and accounts for boundary condition
    index array ordering:
      nnidx = [0, 1, 2, 3, ..., 2D-1, 2D] = [x0+1, x0-1, x1+1, x1-1,..., 2D+1, 2D-1]
    algo:
      1. copy coordinates to new array
      2. create new coord array (x0+1, x0) & get index to assign to nnidx
      2. create new coord array (x0-1, x0) & get index to assign to nnidx
      3. repeat in steps of 2 for i < D
  */
  int nnx[D];
  for (int j=0; j<D; j++) { nnx[j] = spnstrct.coord[j]; }

  for (int i=0; i<2*D; i++) {
    nnx[i/2] = spnstrct.coord[i/2]+pow(-1,i); // because of integer devision in index: i=0 -> [0/2]=0, i=1 -> [1/2]=0, i=0 -> [2/2]=1, ...
    if (nnx[i/2]==N or nnx[i/2]==-1) {
      if (boundary_condition==0) { spnstrct.nnidx[i]=nosite; }
      else if (boundary_condition==1) {
        nnx[i/2]=(nnx[i/2]+N)%N;
        spnstrct.nnidx[i] = n_of_x(nnx);
      }
    }
    else { spnstrct.nnidx[i] = n_of_x(nnx);}
    nnx[i/2] = spnstrct.coord[i/2];
  }
}


/*******************************************************************************
BLACK/WHITE ORDERING
*******************************************************************************/
static void set_spinarray_blackwhite(spinstruct *spinstruct_arr) {
  // GOAL: assigns elements (spin values) to an array in a way that it contains all spin values for the whole lattice
  //        this array can be called with the correct index and returns the element at this index
  //          that is: the spin value of spin at this lattice site where index points to
  // (assigns spin values to)? an array with all even lattice points first (black)
  //  and all white lattice points second (white)
  // above is the coordinate ordering, ie. first half indeces are even lattice points (black)
  // second half are odd indeces (white) (even and odd refer to the parity of coordinate tuple)
  // array is N^d long
  // array elements are structs with spinvalue, own index, own coordinate and next neighbor indeces
  // -> function needs acces to spin module
  // an array has to exist where the needed and elements can be taken by using the index -> STRUCT idea?
  /*
  1. set_params to get N and d
  2. allocate arrray of structs with size N^d
  3. loop through array and assign each point:
    -  the corresponding coord.
    -  a spin value
    - its index
    - next neighbor indeces (maybe by calling the corresponding function)
  */
  // spinstruct *spinstruct_arr;
  set_params();
  spinstruct_arr = (spinstruct*) malloc(sizeof(spinstruct)*pow(N, D));

  for (int i=0; i<pow(N, D), i++){
    spinstruct_alloc(spinstruct_arr[i]);
    spinstruct_arr[i].idx = i;
    // spinstruct_arr[i].spinval = function_that_assigns_spinvalue(); maybe do this when initalizing the hot/cold start
    set_coord_blackwhite(i, spinstruct_arr[i].coord);
    set_nnarray_blackwhite(spinstruct_arr[i]);
  }
}

static inline int parity(int *x) { // static inline if only needed in this file -> make declaration of all function at the top of file (with small comment and maybe even line)
  int sum;
  for (xi=0; xi<D; xi++) {
    sum += x[xi];
  }
  return sum%2;
}

static inline int np_parity(int np) {
  return ((int) 2*np) / ((int) pow(N,D));
}

// maybe not set but "idx2coord" or something, as it returns something?
//inline void set_coord_blackwhite(); ??
static void set_coord_blackwhite(int np, int *x) {
  /* np = current index, x = coordinate vector */
  int xprime[D]; // dynamic allocation better?
  for (int mu=0; mu<D; mu++) {
    xprime[mu] = ( 2*np-(2*np/pow(N,D))*pow(N,D)/pow(N,mu) )%N;
  }
  for (int xi=0; xi<D; xi++) {
    if (xi==0) { x[xi] = ( xprime[xi]+(np_parity(np)+parity(xprime))%2 )%N; }
    else { x[xi] = ( xprime[xi]+(np_parity(np)+parity(xprime))%2 * (xprime[mu-1]==0) )%N; }
  }
}



/*******************************************************************************
LEXOGRAPHIC ORDERING
*******************************************************************************/
static void set_spinarray_lexo(spinstruct *spinstruct_arr) {
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

static inline void set_coord_lexo(int n, int *x) {
  /* n = current index, x = coordinate vector */
  for (int mu=0; mu<D; mu++) {
    x[mu] = (n/pow(N,mu))%N;
  }
}

/*******************************************************************************
Other public functions
*******************************************************************************/
void set_spinarray(spinstruct *spinstruct_arr) {
  if (boundary_condition==0) { set_spinarray_lexo(spinstruct_arr); }
  else if (boundary_condition==0) { set_spinarray_blackwhite(spinstruct_arr); }
  else {
     printf("[geometry.c | set_spinarray() ] ERROR. Upsi somethinf woring with the boundary condition. Remember 0=Dirichlet & 1=periodic.\n");
     exit(-1);
   }
}
