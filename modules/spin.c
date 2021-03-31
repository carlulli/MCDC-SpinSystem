/***************************************************+***************************
This module should contain:
- the three spin models
- a wrapper function for the spin models
- a function that proposes new spin values
  * the propsed spin values also depend on the spin model
  * here the random seed counter needs to me remembered. (im not sure to what extend this is relevant)
- probably a function to intialize spin values for a hot start (all spins in opposite directions)
- in that case also a function for a cold start (all spins aligned in one direction)
***************************************************+***************************/

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>






void aligned_start_dir( int direction)
{
  /* Function: It gives every spin the orientation specified by direction. Works for both the Ising, and the Clock model.
		 					 The direction can be defined in the parameter module.
		 Input: array of the spinstructs, direction of spin
  */
  for(x=0, x< pow(N, D), x++)
  {
    spinstruct_arr[x].spinval = direction;
  }


}

void cold_start()  // in the cold start all spins point the same way.
{
	if (spinmodel == 0) //Ising
	{
		for(x=0, x< pow(N, D), x++)
	  {
	    spinstruct_arr[x].spinval = 1;
	  }
	}

	if (spinmodel == 1) //Clock
	{
		for(x=0, x< pow(N, D), x++)
	  {
	    spinstruct_arr[x].spinval = 0;
	  }
	}

}



// void anti_aligned_start_black_white( int direction)
// {
// 	/* Function: Gives each spin the opposite direction than its neighboring spins.
// 	 						 In the black and white ordering the first half of the spins (all white) are  not in contact with each other,
// 							 but only in contact with spins out of the second half (all black), which also are not in contact with each other
// 							 (think of chessboard).
// 		 Input: array of the spinstructs, direction of spin
// 	*/
//   for( x=0 , x<pow(N,D)/2 , x++ )
//   {
//     // each black spin is aligned one way
//     spinstruct_arr[x].spinval = direction;
//   }
//   for( x=pow(N,D)/2 , x<pow(N,D) , x++ ){
//     // each white spin is aligned the opposite way
//     spinstruct_arr[x].spinval = - direction;
//   }
//
// }





void hot_start()
{
	/* Function: 	Gives each spin a random value, making it maximally unordered, therefore called hot.
		 Input: 		spinstruct array
	*/
  int placeholder;
  if (spinmodel == 0) //random spinset for Ising model, with spin -1 or 1.
	{
    for( x=0 , x<pow(N,D) , x++ ){
      placeholder = rand() % 2  ;	// randomly generates either a 0 or a 1.
      if (placeholder == 0){placeholder = -1 ;} // the zero than is switched to -1, thereby making it the spinvalue we actually want.
      spinstruct_arr[x].spinval = placeholder ;
    }
  }
  else if (spinmodel == 1){ // random spinset for the Clock model
    for( x=0 , x<pow(N,D) , x++ ){
      spinstruct_arr[x].spinval =  rand() % M_number_orient;  //  assigning from 0 to M-1 randomly  ... ! the M_number_orient
                                                              // needs to be defined in the parameters file

    }

  }

}

int random_spin_orientation()  // returns an integer which is a randomly selected possible orientation of the spin.
{

	if (spinmodel == 0) //Ising
	{
		random_spin_orientation = rand() % 2  ;
		if (random_spin_orientation == 0){random_spin_orientation = -1 ;}
	}
	else if (spinmodel == 1) //Clock
	{
		random_spin_orientation = rand() % M_number_orient;
	}
	return random_spin_orientation;
}









int spinmultiplication(int spin1, int spin2)
{
	/*	Function:	spin multiplication depending on which spin model is used.
		 	Input:		the two integers defining the spins.

	*/
  // for Ising spin1 and spin2 are the actual spin values, whereas for the Clock model they are the m values which determine the spin.
  int spin_product ;
  if (spinmodel == 0){  // Ising
    spinproduct = spin1 * spin2 ;
  }
  else if (spinmodel == 1){ // Clock
    spinproduct = (spin1 - spin2 + M_number_orient) % M_number_orient ;
  }
	return spinproduct;
}


// how does the spin addition work ?      in case of the  Isingmodel we just sum the spinval up, and for the clock model we sum up the the values
// given in an array with the spinval as indices.
// therefore we create a vector . then when adding spins under this model
//  ! is it really ok to only take the real part here?


double * set_spinval_values(void)
{ /*	Function: calculates the real values of each spin, that is defined by the direction m.*/
	double *spinval_values = (double*) malloc(sizeof(double)* M_number_orient) ;
  for (i = 0 , i < M_number_orient){
    spinval_values[i] = cos(2 * M_PI * i / M_number_orient ) ;     //
  }
	return spinval_values ;
}

void free_spinval_values(){
	free(spinval_values);
	spinval_values = NULL ;
}



/************************************************************************************************
two-point spin-spin correlation function: based on formula 2.21 in the description.pdf
**************************************************************************************************/
double spin_spin_correlation_r(int r)
{/*****************************************************************************
	 Function: two-point spin-spin correlation function
		Input: Distance r for which one wants to know the correlation_sum
		Return: spin-spin correlation for distance r
	*****************************************************************************/
	int coordvector_spin1[D]; //!
	int coordvector_spin2[D];
	int index_spin1 ;
	int index_spin2 ;
	int orient_spin1;
	int orient_spin2;
  int index_product ;
	double correlation_sum ;
	// different sums for the different dimensions
	// 1 Dimension
	if (D == 1)
	{
		for( x1=0 , x1< N-1 , x1++ )
		{
			// first the coordinate vectors get created using the indices of the sum.
			coordvector_spin1[0] = x1 ;
			coordvector_spin2[0] = x1 + r ;
			// here the coordinates need to be translated to the indices
			index_spin1 = n_of_x(coordvector_spin1);
			index_spin2 = n_of_x(coordvector_spin2);
			// then the spinorientations at those indices are looked up
			orient_spin1 = spinstruct_arr[index_spin1].spinval;
			orient_spin2 = spinstruct_arr[index_spin2].spinval;
			// here the two spins are multiplied  ! again the question arises whether taking the real part is the right thing to do

			index_product = spinmultiplication(orient_spin1,orient_spin2); // index of the spinproduct for each summand
			// ! and here the spinvalue corresponding to that
			correlation_sum += spinval_values[index_product];
		}
	}

	// 2 dimensions
	else if (D == 2)
	{
		for(x1=0, x1< N-1, x1++)
		{
			for(y1=0, y1< N-1, y1++)
			{
				for(y2=0, y2< N-1, y2++)
				{
					// first the coordinate vectors get created using the indices of the sum.
					coordvector_spin1[0] = x1 ;
					coordvector_spin2[0] = x1 + r ;
					coordvector_spin1[1] = y1 ;
					coordvector_spin2[1] = y2 ;
					// here the coordinates need to be translated to the indices
					index_spin1 = n_of_x(coordvector_spin1);
					index_spin2 = n_of_x(coordvector_spin2);
					// then the spinorientations at those indices are looked up
					orient_spin1 = spinstruct_arr[index_spin1].spinval;
					orient_spin2 = spinstruct_arr[index_spin2].spinval;
					// here the two spins are multiplied  ! again the question arises whether taking the real part is the right thing to do

					index_product = spinmultiplication(orient_spin1,orient_spin2); // index of the spinproduct for each summand
					// ! and here the spinvalue corresponding to that
					correlation_sum += spinval_values[index_product];
				}
			}
		}
	}

	// 3 dimensions
	else if (D == 2)
	{
		for(x1=0, x1< N-1, x1++)    // the huge sum over x1,y1,y2,z1,z2 starts. with regular coordinates over the whole space
		{
			for(y1=0, y1< N-1, y1++)
			{
				for(y2=0, y2< N-1, y2++)
				{
					for(z1=0, z1< N-1, z1++)
					{
						for(z2=0, z2< N-1, z2++)
						{
							// first the coordinate vectors get created using the indices of the sum.
							coordvector_spin1[0] = x1 ;
							coordvector_spin2[0] = x1 + r ;
							coordvector_spin1[1] = y1 ;
							coordvector_spin2[1] = y2 ;
							coordvector_spin1[2] = z1 ;
							coordvector_spin2[2] = z2 ;
							// here the coordinates need to be translated to the indices
							index_spin1 = n_of_x(coordvector_spin1);
							index_spin2 = n_of_x(coordvector_spin2);
							// then the spinorientations at those indices are looked up
							orient_spin1 = spinstruct_arr[index_spin1].spinval;
							orient_spin2 = spinstruct_arr[index_spin2].spinval;
							// here the two spins are multiplied  ! again the question arises whether taking the real part is the right thing to do

							index_product = spinmultiplication(orient_spin1,orient_spin2); // index of the spinproduct for each summand
							// ! and here the spinvalue corresponding to that
							correlation_sum += spinval_values[index_product];


						}
					}
				}
			}
		}
	}
	return correlation_sum / int_pow(N,D) ;   // ! here the correlation sum is returned, after being divided by the volume.
}

double * set_spincorrelation_vector(void)
{	/********************************************************************************************************************
		Function: creates a vector filled with the two-point spin_correlation functions for different r going from 0 to N/2
	**********************************************************************************************************************/
N = get_N();	 // length of one dimension
double *spin_corr_vec = (double*) malloc(sizeof(double)* N/2) ;  // ! check
/* memory is allocated for the vector spin_corr_vec, which is going to be filled with the spin correlation values for the N/2 different r (distance) values.  */
	for (r=0, r<N/2, r++)
	{
		spin_corr_vec[r] = spin_spin_correlation_r(r) ;
	}
	return spin_corr_vec ;
}

void free_spin_corr_vec(){ /* freeing up the spin correlation vector */
	free(spin_corr_vec);
	spin_corr_vec = NULL ;
}


/********************************************************************************
Setting the physical parameters
------------------------------
we need   the M_number_orient   which is number of orientations of the Clock model
we need the magnetic field  B
we need to select the spinmodel Ising 0 or clock 1




*********************************************************************************/
