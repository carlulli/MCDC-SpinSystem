/*******************************************************************************
Test code for geometry module
tests if eg. coord for periodic and lexo ordering are same
idea: create array with index loop and with coordinate loop
-> passed if same combinations arrise
ONLY FOR 2D at the moment
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>

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

  // printf("DEBUGG PRINT TYPE error trick %s\n", spinstruct_arr[0]->coord);
  // printf("DEBUGGING\n");
  // exit(-1);
  // int index=0;
  // index = spinstruct_arr[0].idx;


  printf(
    "Your chosen system: \n Size in one direction = %d\t Dimension = %d \t boundary_condition = %s \t ordering = %s \n",
    get_N(), dim, argv[3], argv[4]
  );
  printf("LOOPING THROUGH INDECES.\n");
  printf("Spin index\t| Coordinates\n");
  printf("--------------------------------\n");
  for (int i=0; i<arrsize; i++) {
      printf("%d\t|\t(%d, %d)\n", spinstruct_arr[i].idx, spinstruct_arr[i].coord[0], spinstruct_arr[i].coord[1]);
  }
  spinarray_free(spinstruct_arr);
  spinstruct_arr = set_spinarray_coord2d();

  printf(
    "Your chosen system: \n Size in one direction = %d\t Dimension = %d \t boundary_condition = %s \t ordering = %s \n",
    get_N(), dim, argv[3], argv[4]
  );
  printf("LOOPING THROUGH COORDINATES.\n");
  printf("Spin index\t| Coordinates\n");
  printf("--------------------------------\n");
  for (int i=0; i<arrsize; i++) {
      printf("%d\t|\t(%d, %d)\n", spinstruct_arr[i].idx, spinstruct_arr[i].coord[0], spinstruct_arr[i].coord[1]);
  }

  spinarray_free(spinstruct_arr);
  return 0;
}
