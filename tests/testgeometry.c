/*******************************************************************************
Test code for geometry module
tests if eg. coord for periodic and lexo ordering are same
idea: create array with index loop and with coordinate loop
-> passed if same combinations arrise
ONLY FOR 2D at the moment
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "geometry.h"


int main(int argc, char const *argv[]) {
  /* input at runtime: file name [N] [D] [boundary_condition] [ordering]  */
  set_params(argc, argv);

  // spinstruct_t *spinarray=NULL;
  // set_spinarray(spinarray);
  spinstruct_t *spinstruct_arr = set_spinarray();
  // spinarray = spinstruct_arr;

  int arrsize, dim;
  arrsize=get_arraysize();
  dim=get_D();
  printf("\nThis program prints multiple index and coordinate values of a 4x4 spin system to check if ordreing works corretly.\n\n");
  printf("Chosen ordering (lexo=0, blackwhite=1) and boundary condition (Dirichlet=0, Periodic=1) as recovered (lexo=0, blackwhite=1):\n"
          "\tordering=%d, boundary condition=%d\n",
          get_ordering(),
          get_boundary_condition());
  printf("\nGOAL ORDERING for N=4 and D=2 (a=1):\n");
  if (atoi(argv[4])==1) {
    printf("Blackwhite Ordering:\n");
    printf("\ty | \n"
          "\t3-| 14  6 15  7  \n"
          "\t2-|  4 12  5 13  \n"
          "\t1-| 10  2 11  3  \n"
          "\t0-|  0  8  1  9  \n"
          "\t  ---------------\n"
          "\t     0  1  2  3 x\n");
  }
  else if (atoi(argv[4])==0) {
    printf("Lexographic Ordering:\n");
    printf("\ty | \n"
          "\t3-| 12 13 14 15  \n"
          "\t2-|  8  9 10 11  \n"
          "\t1-|  4  5  6  7  \n"
          "\t0-|  0  1  2  3  \n"
          "\t  ---------------\n"
          "\t     0  1  2  3 x\n");
  }

  printf(
    "Your chosen system: \n Size in one direction = %d\t Dimension = %d \t boundary_condition = %s \t ordering = %s \n",
    get_N(), dim, argv[3], argv[4]
  );
  printf("LOOPING THROUGH INDECES.\n");
  printf("\tSpin index\t| Coordinates\n");
  printf("\t--------------------------------\n");
  for (int i=0; i<arrsize; i++) {
      printf("\t%d\t|\t(%d, %d)\n", spinstruct_arr[i].idx, spinstruct_arr[i].coord[0], spinstruct_arr[i].coord[1]);
  }

  printf("Some example coordinate distances for periodic ordering: \n"
         "\tdist(0,7) = %f\n" "\tdist(12,11) = %f\n" "\tdist(3,14) = %f\n"
         "\tdist(4,10) = %f\n" "\tdist(0,15) = %f\n" "\tdist(0,14) = %f\n",
         coord_distance(&spinstruct_arr[0], &spinstruct_arr[7]),
         coord_distance(&spinstruct_arr[12], &spinstruct_arr[11]),
         coord_distance(&spinstruct_arr[3], &spinstruct_arr[14]),
         coord_distance(&spinstruct_arr[4], &spinstruct_arr[10]),
         coord_distance(&spinstruct_arr[0], &spinstruct_arr[15]),
         coord_distance(&spinstruct_arr[0], &spinstruct_arr[14]) );

  printf("To compare with: \n");
  printf("Hypothenuse of 1x an 0y = %f\n", sqrt(1));
  printf("Hypothenuse of 1x an 1y = %f\n", sqrt(2));
  printf("Hypothenuse of 2x an 1y = %f\n", sqrt(5));
  printf("Hypothenuse of 3x an 1y = %f\n", sqrt(10));
  printf("Hypothenuse of 3x an 2y = %f\n", sqrt(13));
  printf("Hypothenuse of 3x an 3y = %f\n", sqrt(18));

  printf("\nWhat are the NEXT NEIGHBORS?\n");
  printf("Next neihbors for: \n");
  if (atoi(argv[3])==0) { printf("Dirichlet Boundary Conditions.\n");}
  else if (atoi(argv[3])==1) { printf("Periodic Boundary Conditions.\n"); }
  int idx[5] = {0, 2, 11, 4, 15};
  for (int i=0; i<5; i++) {
    printf("\tNNarray(%d)=[%d, %d, %d, %d]\n",
            spinstruct_arr[idx[i]].idx,
            spinstruct_arr[idx[i]].nnidx[0],
            spinstruct_arr[idx[i]].nnidx[1],
            spinstruct_arr[idx[i]].nnidx[2],
            spinstruct_arr[idx[i]].nnidx[3] );
  }

  printf("\nTest program finished.\n\n");
  spinarray_free(spinstruct_arr);
  return 0;
}
