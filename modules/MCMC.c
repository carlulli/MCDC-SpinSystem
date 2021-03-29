/***************************************************+***************************
Method:
we use the 1. proposed updating method (refering to desrciption) -> MARS step for each spin value proposal
for one iteration: We update first all even then all odd
Structure:
1. loop over all even spins:
    - update all even spins
2. MARS (accept/reject) step
3. loop over all odd spins:
    - update all odd spins
    - background can be calculated with already changed spins (but only if MARS was performed on them)
4. MARS step
***************************************************+***************************/

/***************************************************+***************************
Functions:
x -> index because of ordering: spin value at this lattice site = spinstruct_arr[idx]->spinval

+ BACKGROUND(x) (calculates background for one spin x )
  - b(x)=sum_(mu<D) s(x+e_mu)+(s(x-e_mu))
  - we need background for all updated spins
    -> should be able to calc through loop, as background of even spin does only depend on odd (and vice versa)
  - just loop over all next neighbors?
  - b(x)b(idx)=for (i<D) b += spinstruct_arr[spinstruct_arr[idx].nnidx[i]].spinval

+ HAMILTONIAN(x)_p,n (for one x)
  - -2Re[s_n,p(x)*b(x)]-B*Re[s_n,p(x)]
  - is this the same hamiltonian as in the hamiltonian module?
    -> should it be?

+ MARS accept/reject step for a number of flipped spins, but works same as with one:
    - input values:
            * proposed spin configuration (pointer to array with proposed values)
            * old spin value configuration (existing global variable spinstruct_arr)
            * sum of proposed spin value hamiltonians (for hamiltonians, background needs to be calculated)
            * sum of old spin value hamiltonians
    - is it model independent in general or does it need a wrapper function?
    - has if & if else & else condition with sum of spins and sum of hamiltonians
    - calculates value for which to accept: q(x) = min(1, e^{-1/kT(h_p(x)- h_n(x)))
    - what does it do: change or doesnt change the spinvalue (spinstruct_arr[i].spinval) for all even/odd spins

+ UPDATE loop through one set of indeces (even or odd)
    - for all even or all odd loop:
      - update per spin is performed by spin module

+ MCMC_step function that does one MCMC update step:
    1. loop over all even spins:
        - update all even spins
    2. MARS (accept/reject) step
    3. loop over all odd spins:
        - update all odd spins
        - background can be calculated with already changed spins (but only if MARS was performed on them)
    4. MARS step
***************************************************+***************************/
