/*******************************************************************************
Header file of geometry module.
Here only global functions are listed.
At module's beginning, a detailed list of functions is written.
*******************************************************************************/

struct spinstruct;

void set_params(argc, argv[]) (H)
/* sets all nevessary parameters to create geometry */
int get_N() (H)
/* returns N */
int get_D() (H)
/* returns D */
void spinarray_free(spinstruct *spnstrct_arr);
/* frees memory for array of struct (of size N^D) (wrapper of 3.) */
double coord_distance(spinstruct *spnstrct1, spinstruct *spnstrct2)
/* calculates distance btw 2 spinindeces (wrapper for 5. & 6.) */
void set_spinarray(spinstruct *spinstruct_arr)
/* sets spinstructs for given array and ordering -> wrapper function for 8. & 13. */
