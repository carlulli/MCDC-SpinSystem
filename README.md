# MCMC-SpinSystem
This code simulates a ferromagnetic spin system using Markov-Chain Monte Carlo simulation for updating and a proper statistical analysis.

### Compilation and Running the code
- to **compile** run `make` in the working directory
- adjust your parameters in the *params.txt* file (the structure of the file is essential to work for the python script, so don't change it)
- to **run** run the python script *runSimulation.py* with python3 by running `python3 ./runSimulation.py` in the working directory
  - Note: currently *pandas* is required, because that was the easiest way to implement in short time. Obviously an implementation with the least amount of packages would be optimal 
  - For now, I have not figuered out how to properly kill a running process though (strg + c, strg + z)
