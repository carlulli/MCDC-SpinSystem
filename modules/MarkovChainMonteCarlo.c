/***************************************************+***************************
Method:
we use the 1. proposed updating method (refering to desrciption) -> MARS step for each spin value proposal
for one iteration: We update first all even then all odd 
Structure:
1. loop over all even spins:
    - update all even spins
2. MARS (accept/reject) step uh<yfghjlö
3. loop over all odd spins:
    - update all odd spins
    - background can be calculated with already changed spins (but only if MARS was performed on them)
4. MARS step
***************************************************+***************************/
