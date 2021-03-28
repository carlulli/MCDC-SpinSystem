/*******************************************************************************
Module with small usefull extra shit
*******************************************************************************/

/* "int_pow" calculates the power of a given base to the exponent */
int int_pow(int base, int exponent) {
    int result = 1;
    for( int i = 1; i <= exponent; i++)
        result *= base;

    return result;
}
