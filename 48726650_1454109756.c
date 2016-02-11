/* mymath.c v. 0.57 */

#include "48726651_1454109756.h"

/************************  float SquareIt(float it)  ***********************/
/*  */
/* FUNCTION: calculates square */
/* ARGUMENTS:  number to square */
/* RETURN: result of calculation */
/* PROTOTYPE IN: mymath.h */
/* OTHER DEPENDENCIES: none */
/* NOTE: obviates some type casting problems with pow */
/*  */
/*****************************************************************************/

float SquareIt(float it)
{
return (it * it);
} /* end SquareIt */

/************************  int RoundUp(float id)  ***********************/
/* */
/* FUNCTION: finds nearest integer value of float, rounding up */
/* ARGUMENTS:  float */
/* RETURN: result as int type */
/* PROTOTYPE IN: mymath.h */
/* OTHER DEPENDENCIES: none */
/* */
/*****************************************************************************/

int RoundUp(float id)
{
int rnd;

if (id - (int) id >= 0.5) rnd = (int) id + 1.00;
else rnd = (int) id;
return(rnd);
}

/************************  double Factorial (int num)  ***********************/
/* */
/* FUNCTION: calculates factorial */
/* ARGUMENTS:  number of which to take factorial */
/* RETURN: result */
/* PROTOTYPE IN: mymath.h */
/* OTHER DEPENDENCIES: none */
/* NOTE: recursive, stops calling itself with input-1 when input = 1 */
/* */
/*****************************************************************************/

double Factorial (int num)
{
double fac;

fac = (double) num;

if (fac > 1) fac *= Factorial(fac - 1);

return (fac);
}

