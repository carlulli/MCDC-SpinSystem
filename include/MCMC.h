/*******************************************************************************
Header file of Markokv Chain Monte Carlo module.
List of all functions in .c file
Here are all public functions
*******************************************************************************/

#ifndef MCMC_H
#define MCMC_H

void init_MCMC(void);
/*
  allocats memory for spinval_values, spin_corr_vec and spinstruct_arr
  by calling corresponding functions
 */

void MCMC_step(void);
/*
  performs one whole updating step
   first all even, then all
*/

#endif /* MCMC_H */
