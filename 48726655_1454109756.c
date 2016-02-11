/* randum.c v. 0.58 */

/************ initialization of seed for random number generator *************/
/* */
/* NOTE: initialize C's rand() generator by current time, then get a random */
/* number to be used as a seed.  This initialization is performed by calling */
/* InitRanSeed.  Call this at startup, then do not call it again. */
/* */
/*****************************************************************************/

#include <stddef.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "48726656_1454109756.h"

extern long gSeed; /* initialize at startup with InitRanSeed */

/************************** long InitRanSeed(void)  **************************/
/* */
/* FUNCTION: initializes rand and gets seed for NRuran2 */
/* ARGUMENTS: none */
/* RETURN: seed for NRuran2 */
/* PROTOTYPE IN: randum.h */
/* OTHER DEPENDENCIES: */
/* NOTE: */
/* */
/* syntax is globalSeed = InitRanSeed(); */
/* */
/*****************************************************************************/
long InitRanSeed(void)
{
srand((unsigned) time(NULL));
return(-rand());
}


/***************************  float NRuran2(void)  ***************************/
/* */
/* FUNCTION: high-falutin' pseudo-random number generator */
/* ARGUMENTS:  none */
/* RETURN: random number from 0-1, exclusive of endpoints */
/* PROTOTYPE IN: randum.h */
/* OTHER DEPENDENCIES:  randum.h has the necessary #defines */
/* NOTE: code is from p. 282 of _Numerical Recipes in C_: 'long period' */
/* random number generator of L'Ecuyer, with Bays-Durham shffle and added */
/* safeguards.  Returns a uniform random deviate between 0.0 and 1.0 */
/* (exclusive of the endpoint values).  Call with seed a negative integer to */
/* initialize; therafter, do not alter it between successive deviates in a */
/* sequence.  kRNMX should approximate the largest floating value that is */
/* less than 1. */
/* */
/*****************************************************************************/

float NRuran2(void)
{
int j;
long k;
static long seed2 = 123456789L;
static long iy = 0;
static long iv[kNTAB];
float temp, temp2;

if (gSeed <= 0)
	{
	if (-(gSeed) < 1) gSeed = 1;
	else gSeed = -(gSeed);
	seed2 = gSeed;
	for (j = kNTAB + 7; j >= 0; j--)
		{
		k = gSeed/kIQ1;
		gSeed = kIA1 * (gSeed - k * kIQ1) - kIR1 * k;
		if (gSeed < 0) gSeed += kIM1;
		if (j < kNTAB) iv[j] = gSeed;
		}
	iy = iv[0];
	}
k =  (gSeed)/kIQ1;
gSeed =  kIA1 * (gSeed - k * kIQ1) - kIR1 * k;
if (gSeed < 0) gSeed += kIM1;

k =  seed2/kIQ2;
seed2 =  kIA2 * (seed2 - k * kIQ2) - kIR2 * k;
if (seed2 < 0) seed2 += kIM2;

temp2 = iy / kNDIV;
j = (int) temp2;  /* trying to get rid of compiler warning */
iy = iv[j] - seed2;
iv[j] = gSeed;

if (iy < 1) iy += kIMM1;
if ((temp = kAM * iy) > kRNMX) return kRNMX;
else return temp;

} /* end of nruran2 generator */

/*****************  float NormRand(float mean, float stdDev)  ****************/
/*  */
/* FUNCTION: generates random numbers with a normal frequency distribution */
/* ARGUMENTS:  parameters mu and std dev, here taken from the observed set */
/* RETURN:  a normal random number */
/* PROTOTYPE IN: randum.h */
/* OTHER DEPENDENCIES: */
/* NOTE: Based on Box and Muller's method of generating normal random */
/* numbers from uniform random numbers.  The return of NRuran2 is converted */
/* into a normally distributed number with mean = 0 and stdDev = 1; this is */
/* then modified with the desired mu and stdDev */
/*  */
/*****************************************************************************/

float NormRand(float mean, float stdDev)
{
float fac, rsq, v1, v2;
static int iset = 0;
static float gset;

if (iset == 0)
	{
	do
		{
		v1 = 2.0 * NRuran2() - 1.0;
		v2 = 2.0 * NRuran2() - 1.0;
		rsq = v1 * v1 + v2 * v2;
		}
	while (rsq >= 1.0 || rsq == 0.0);

	fac = sqrt(-2.0 * log(rsq)/rsq);
	gset = mean + v1 * fac * stdDev;
	iset = 1;
	return mean + v2 * fac * stdDev;
	}
else
	{
	iset = 0;
	return gset;
	}

}/***  end of NormRand generator **/

