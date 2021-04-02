#ifndef SPIN_H
#define SPIN_H

extern double *spinval_values ;
extern double *spin_corr_vec ;

void set_physics_params (int argc, const char *argv[]);
void aligned_start_dir( int direction);
/*  gives every spin the orientation specified by direction */
void cold_start();
/* puts all the spins upwards  1 for Ising and orientation 0 for Clock model */
void hot_start();
/* gives each spin a random orientation */
int spinmultiplication(int spin1, int spin2);
/* gives the new orientation for the spin product, taking as input the two orientations of the two spins.*/
double * set_spinval_values(void); // ! uses spinval_values vector
/* calculates the real values of each spin, that is defined by the direction m */
double spin_spin_correlation_r(int r);
/*two-point spin-spin correlation function for the distance r */
double * set_spincorrelation_vector(void); // !
/* creates a vector filled with the spin_correlation functions for different r going from 0 to N/2 */
void free_spin_corr_vec();
 /* freeing up the spin correlation vector */
void free_spinval_values();
/*frees spinval values*/
int random_spin_orientation() ;

double get_B(void);

int get_seed(void) ;

int get_M_number_orient(void);

int get_spinmodel(void);

double get_T(void);


// void anti_aligned_start_black_white( int direction);    might not be needed







#endif
