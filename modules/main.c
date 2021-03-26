/***************************************************+***************************
Main code:
- reads in parameter file
- initializes lattice, i.e., geomerty
  -> maybe initialze a pointer in the main that can be given to the geometry and is from there on the N^d array with all structs inside
- initialze hot or cold start
- loop through iterations of MCMC
  - write out energy and magnetization after each iteration
  - wirte out correlation function after maybe each 10th iteration
  - write out spin configuration after maybe each 100th iteration
- free array 
***************************************************+***************************/
