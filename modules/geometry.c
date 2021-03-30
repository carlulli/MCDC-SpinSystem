/*******************************************************************************
This module handles the transformation of the real lattice to the discretized
version used for computing.
*******************************************************************************/

/*******************************************************************************
List of functions. (some declared here (S), others declared in header (H))
1. void set_params(argc, argv[]) (H)
    - sets all nevessary parameters to create geometry
    - input params: argc=5; argv= [N] [D] [boundary_condition] [ordering]
2. int get_N() (H)
    - returns N
3. int get_D() (H)
    - returns D
4. int get_ordering() (H)
    - returns ordering (so far without any double checking)
5. int get_boundary_condition() (H)
    - returns boundary condition (so far without any double checking)
6. static void spinstruct_alloc(spinstruct_t *spstrct) (S)
    - allocates memory for struct arrays coord and nextneighbor_index
7. static inline void spinstruct_free(spinstruct_t *spnstrct) (S)
    - frees for struct allocated memory
8. void spinarray_free(spinstruct_t *spnstrct_arr) (H)
    - frees memory for array of struct (of size N^D) (wrapper of 3.)
9. static double coord_distance_dirichlet(spinstruct_t *spnstrct1, spinstruct_t *spnstrct2) (S)
    - calculates distance btw 2 spinindeces for Dirichlet bc
10. static double coord_distance_periodic(spinstruct_t *spnstrct1, spinstruct_t *spnstrct2) (S)
    - calculates distance btw 2 spinindeces for periodic bc
11. double coord_distance(spinstruct_t *spnstrct1, spinstruct_t *spnstrct2) (H)
    - calculates distance btw 2 spinindeces (wrapper for 5. & 6.)
12. static inline int parity(int *x) (S)
    - returns parity of coordinate array
13. static inline int np_parity(int np) (S)
    - returns special parity 2np/n^D
14. static inline int n_of_x_lexo(int *x) (S)
    - calculates index of coordinate for lexo ordering
15. static inline int n_of_x_blackwhite(int *x) (S)
    - calculates index of coordinate for blackwhite ordering
16. static int n_of_x(int *x) (S)
    - returns index of coordinate for both orderings (wrapper of 10. and 11.)
17. static void set_nnarray(spinstruct_t *spnstrct) (S)
    - sets a next neighbor index array for given struct
18. static void set_coord_blackwhite(int np, int *x) (S)
    - sets the coordinates in coordinate vector for given index for blackwhite ordering
19. static spinstruct_t* set_spinarray_blackwhite(spinstruct_t *spinstruct_arr) (S)
    - creates structs and assigns values for given array in blackwhite ordering
20. static inline void set_coord_lexo(int n, int *x) (S)
    - sets the coordinates in coordinate vector for given index for lexographic ordering
21. static spinstruct_t* set_spinarray_lexo(spinstruct_t *spinstruct_arr) (S)
    - creates structs and assigns values for given array in lexographic ordering
22. spinstruct_t* set_spinarray(spinstruct_t *spinstruct_arr) (H)
    - sets spinstructs for given array and ordering -> wrapper function for 8. & 13.
23. int get_arraysize() (H)
    - returns size of allocated systems total spin array
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "geometry.h"
#include "utilities.h"


static int D=1; // D = Dimension
static int N=0; // N = Number of lattice sites in one direction (per dimension)
static int lattice_spacing=1; // lattice spacing fixed to 1 in this project! if double -> problem in n_of_x_lexo
static int boundary_condition;
static int *bc_ptr=NULL;
static int ordering=1; // 0 for lexo and 1 for black white
int nosite = -1; // should be unchangeable but open for other modules


/* static function initialition inside module as not needed outside */
static void spinstruct_alloc(spinstruct_t *spstrct);
static inline void spinstruct_free(spinstruct_t *spnstrct);
static double coord_distance_dirichlet(spinstruct_t *spnstrct1, spinstruct_t *spnstrct2);
static double coord_distance_periodic(spinstruct_t *spnstrct1, spinstruct_t *spnstrct2);
static inline int parity(int *x);
static inline int np_parity(int np);
static int n_of_x_lexo(int *x);
static inline int n_of_x_blackwhite(int *x);
static int n_of_x(int *x);
static void set_coord_blackwhite(int np, int *x);
static spinstruct_t* set_spinarray_blackwhite(void);
static void set_coord_lexo(int n, int *x);
static spinstruct_t* set_spinarray_lexo(void);

/*******************************************************************************
Parameter handling
*******************************************************************************/
void set_params(int argc, char const *argv[]) {
  // if argc not some value -> error message and exit
  if (argc < 5) {
    printf("[goemetry.c | set_params()] ERROR: Not enough input parameters. \n"
    "\t\t You need min 5 and had: {%d}\n \t\t Remember: [filename] [N] [D] [boundary condition] [ordering]\n", argc);
  }
  N = atoi(argv[1]);
  D = atoi(argv[2]);
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
    if(sD==0) { sD=D; }
    else {
      if((N!=sD))  {
         printf("[ geometry.c| get_N ] Error! (N) has changed: (%d) -> (%d)\n",sD,D);
         exit(-1);
      }
    }
    return D;
}

int get_ordering() {
  return ordering;
}

int get_boundary_condition() {
  return boundary_condition;
}

/*******************************************************************************
Spin Struct (is this something for the header?) and spinstruct_arr constructor
*******************************************************************************/
static void spinstruct_alloc(spinstruct_t *spnstrct) {
  //allocate only for arrays?
  // safety feature to not allocate twice (initialzed pointers in struct are NULL )
  if (spnstrct->coord==NULL | spnstrct->nnidx==NULL) {
    spnstrct->coord = (int *) malloc(sizeof(int)*D);
    spnstrct->nnidx = (int *) malloc(sizeof(int)*2*D);
  }
  else {
    printf("[geometry.c | spinstruct_alloc()] ERROR: apparently struct values are already allocated. \n");
    exit(-1);
  }
}

static inline void spinstruct_free(spinstruct_t *spnstrct) {
  /* function to deallocate (free) memory used for spinstruct spnstrct */
  free(spnstrct->coord);
  free(spnstrct->nnidx);
  spnstrct->coord = NULL;
  spnstrct->nnidx = NULL;
}

void spinarray_free(spinstruct_t *spnstrct_arr) {
  for (int i=0; i<int_pow(N,D); i++) {
    spinstruct_free(&spnstrct_arr[i]);
  }
  free(spnstrct_arr);
  spnstrct_arr = NULL;
}

/*******************************************************************************
Functions for distance calculation
*******************************************************************************/
static double coord_distance_dirichlet(spinstruct_t *spnstrct1, spinstruct_t *spnstrct2) {
  /*
    function that caluates distance btw 2 spins for dirichlet bc
    input are spinstructs -> coordinates already known (using vector substr.)
  */
  int sum=0; // or double sum=0.0?
  for (int i=0; i<D; i++) {
    sum += (spnstrct1->coord[i] - spnstrct2->coord[i])*(spnstrct1->coord[i] - spnstrct2->coord[i]);
  }
  return sqrt(sum);
}

static double coord_distance_periodic(spinstruct_t *spnstrct1, spinstruct_t *spnstrct2) {
  /*
    funtion that calclulates distacnce btw 2 spins for periodic bc
    input are spinstructs -> coordinates already known
    algo:
      1. dx = |x1_a - x1_b| for all xi in array[D]
      2. if dx > N/2: reduce dx for N -> dx = dx - N
      3. sum =  dx1^2 + dx2^2 + ... ->  sum can be calculated in same loop
      4. return sqrt(sum) (pythagoras)
  */
  int sum=0;
  int deltaxi[D];
  for (int i=0; i<D; i++) {
    deltaxi[i] = abs(spnstrct1->coord[i] - spnstrct2->coord[i]);
    if (deltaxi[i] > N/2) { deltaxi[i] -= N; }
    sum += deltaxi[i]*deltaxi[i];
  }
  return sqrt(sum);
}

double coord_distance(spinstruct_t *spnstrct1, spinstruct_t *spnstrct2) {
  /*
    wrapper function to calculate the distance between two indeces
    needs some way to respect boundary conditions
    returns the distance as double
  */
  double dist=0;
  if (bc_ptr==NULL) {
    printf("[geometry.c | coord_distance_blackwhite()] ERROR. Boundary Condition was not set!\n");
    exit(-1);
  }
  if (boundary_condition==0) { dist = coord_distance_dirichlet(spnstrct1, spnstrct2); }
  else if (boundary_condition==1) { dist = coord_distance_periodic(spnstrct1, spnstrct2); }
  else {
    printf("[geometry.c | coord_distance_blackwhite()] Upsi something else went wrong with the bc but it got set with the right value.\n");
    exit(-1);
  }
  return dist;
}

/*******************************************************************************
Next neighbor part (for both orderings)
*******************************************************************************/
static inline int parity(int *x) {
  for (int xi=0; xi<D; xi++) {
    sum += x[xi];
  }
  return sum%2;
}

static inline int np_parity(int np) {
  return ( 2*np) / ( int_pow(N,D));
}

static inline int n_of_x_lexo(int *x) {
  /* x is coordinate vector and n is index in lexo ordering  */
  int sum=0;
  for (int mu=0; mu<D; mu++) { sum += ((x[mu]/lattice_spacing)*int_pow(N, mu)); }
  return sum;
}

static inline int n_of_x_blackwhite(int *x) { // also inline?
  /* x = coordinate vector and np = index in black white ordering*/
  int np=0;
  np = (int)(parity(x)*int_pow(N,D)+n_of_x_lexo(x))/2;
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

static void set_nnarray(spinstruct_t *spnstrct) {
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
  for (int j=0; j<D; j++) { nnx[j] = spnstrct->coord[j]; }
  for (int i=0; i<2*D; i++) {
    nnx[i/2] = spnstrct->coord[i/2]+int_pow(-1,i); // because of integer devision in index: i=0 -> [0/2]=0, i=1 -> [1/2]=0, i=0 -> [2/2]=1, ...

    if (nnx[i/2]==N || nnx[i/2]==-1) {
      if (boundary_condition==0) { spnstrct->nnidx[i]=nosite; }
      else if (boundary_condition==1) {
        nnx[i/2]=(nnx[i/2]+N)%N;
        spnstrct->nnidx[i] = n_of_x(nnx);
      }
    }
    else { spnstrct->nnidx[i] = n_of_x(nnx); }
    nnx[i/2] = spnstrct->coord[i/2];
  }
}


/*******************************************************************************
BLACK/WHITE ORDERING
*******************************************************************************/
// maybe not set but "idx2coord" or something, as it returns something?
static void set_coord_blackwhite(int np, int *x) {
  /* np = current index, x = coordinate vector */
  int xprime[D]; // dynamic allocation better?
  for (int mu=0; mu<D; mu++) {
    xprime[mu] = (( 2*np-(2*np/int_pow(N,D))*int_pow(N,D))/int_pow(N,mu) )%N;
  }
  for (int xi=0; xi<D; xi++) {
    if (xi==0) { x[xi] = ( xprime[xi]+(np_parity(np)+parity(xprime))%2 )%N; }
    else {
      x[xi] = (1-parity(xprime))*( (1-np_parity(np))*xprime[xi] + np_parity(np)*(x[xi-1]==0?xprime[xi]+1:xprime[xi]) )
              + (parity(xprime))*( (np_parity(np))*xprime[xi] + (1-np_parity(np))*(x[xi-1]==0?xprime[xi]+1:xprime[xi]) );
    }
  }
}

static spinstruct_t* set_spinarray_blackwhite(void) {
  // with constructor
  /*
  1. set_params to get N and d
  2. allocate temporary arrray of structs with size N^d
  3. loop through array and assign each point:
    -  the corresponding coord.
    -  a spin value
    - its index
    - next neighbor indeces (maybe by calling the corresponding function)
  4. return temporary array
  */
  spinstruct_t *this = (spinstruct_t*) malloc(sizeof(spinstruct_t)*int_pow(N, D));

  for (int i=0; i<int_pow(N, D); i++){
    spinstruct_alloc(&this[i]);
    this[i].idx = i;
    // spinstruct_arr[i].spinval = function_that_assigns_spinvalue(); maybe do this when initalizing the hot/cold start
    set_coord_blackwhite(i, this[i].coord);
    set_nnarray(&this[i]);
  }
  return this;
}

/*******************************************************************************
LEXOGRAPHIC ORDERING
*******************************************************************************/
static void set_coord_lexo(int n, int *x) {
  /* n = current index, x = coordinate vector */
  for (int mu=0; mu<D; mu++) {
    x[mu] = (n/int_pow(N,mu))%N;
  }
}

// static void set_spinarray_lexo(spinstruct_t *spinstruct_arr) {
static spinstruct_t* set_spinarray_lexo(void) {
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
  spinstruct_t *this = (spinstruct_t*) malloc(sizeof(spinstruct_t)*int_pow(N, D));

  for (int i=0; i<int_pow(N, D); i++){
    spinstruct_alloc(&this[i]);
    this[i].idx = i;
    // spinstruct_arr[i].spinval = function_that_assigns_spinvalue();
    set_coord_lexo(i, this[i].coord);
    set_nnarray(&this[i]);
  }
  return this;
}


/*******************************************************************************
Other public functions
*******************************************************************************/
spinstruct_t* set_spinarray(void) {
  spinstruct_t *spnstrctarr_ptr=NULL;
  if (bc_ptr==NULL) { printf("[geometry.c | set_spinarray()] ERROR! Params not set!\n"); exit(-1); }

  if (ordering==0) { spnstrctarr_ptr = set_spinarray_lexo(); }
  else if (ordering==1) { spnstrctarr_ptr = set_spinarray_blackwhite(); }

  else {
     printf("[geometry.c | set_spinarray() ] ERROR. Upsi something went woring with the boundary condition. Remember 0=Dirichlet & 1=periodic.\n");
     exit(-1);
   }
  return spnstrctarr_ptr;
}

int get_arraysize() {
  return int_pow(N,D);
}
