/* noisgene.c v. 0.58 */

#include <stdio.h>
#include <math.h>
#include <stddef.h>

#include "48726637_1454109756.h"
#include "48726648_1454109756.h"
#include "48726649_1454109756.h"
#include "48726656_1454109756.h"
#include "48726658_1454109756.h"

extern GeneData *gGeneDataPtr[];
extern ObsGeneInfo *gCurObsInfo[];

extern int GetNumValues(GeneData*);
extern void SortValues(GeneData*);
extern GeneData* GetGeneStruct(void);
extern void EmptyGeneList(GeneData*);

/* private functions: */

void EvolveASet(GeneData *intronSet, float delFreq, float slideFreq, float slideSD);
void CombineIntronList(GeneData *ptr);
void DivergeIntronPattern(GeneData *ptr, int gen, float delFreq, float slideFreq, float slideSD);


/*************  void EvolveASet(GeneData*, float, float, float)  *************/
/* */
/* FUNCTION: impose 1 generation of divergence on 1 gene structure */
/* ARGUMENTS: ptr to intron positions, and params: deletion frequency, slide frequency */
/* and standard deviation */
/* RETURN: none (the input set is altered) */
/* PROTOTYPE IN: noisgene.h */
/* OTHER DEPENDENCIES: randum.h, genetype.h */
/* NOTE: */
/* */
/*****************************************************************************/

void EvolveASet(GeneData *intronSet, float delFreq, float slideFreq, float slideSD)
{
int i, j, numValues, temp;
float test;

numValues = GetNumValues(intronSet);

for (i = 0; i < numValues; i++)
	{
	test = NRuran2();
	if (test < delFreq)
		{
		fprintf(stderr, " d");
		j = i;
		while (j < numValues - 1)
			{
			intronSet->number[j] = intronSet->number[j + 1];
			j++;
			}
		intronSet->number[j] = 0;
		numValues--;
		i--;
		}
	else
		{
		if ((test - delFreq) < slideFreq)
			{  /* Do . . . while is to make sure it doesn't slide out of the gene! */
			fprintf(stderr, " s");
			do
				{
				temp = intronSet->number[i];
				temp += (int) NormRand(0.00, slideSD);
				}
			while ((temp < 1) || (temp >= gCurObsInfo[INTRON]->size));
			intronSet->number[i] = temp;
			}
		else fprintf(stderr, " .");
		}
	}
SortValues(intronSet);

} /* end of EvolveASet */

/*****************   void CombineIntronList(GeneData *ptr)  *****************/
/*  */
/* FUNCTION: extracts a single set of distinct intron positions from a list */
/* of sets containing redundancies */
/* ARGUMENTS: ptr to first member of list */
/* RETURN: none */
/* PROTOTYPE IN: noisgene.h */
/* OTHER DEPENDENCIES: genetype.h */
/*  */
/*****************************************************************************/

void CombineIntronList(GeneData *ptr)
{
int i, numValues, nextNumValues;
boolean tooMany = false;

fprintf(stderr, "\ncombining . . . ");
while (ptr->nextSet != NULL)
	{
	numValues = GetNumValues(ptr);
	nextNumValues = GetNumValues(ptr->nextSet);
	if (nextNumValues)
		{
		for (i = 0; i < nextNumValues; i++)
			{
			if (numValues + i < kMaxNumValues)
				{
				ptr->number[numValues + i] = ptr->nextSet->number[i];
				}
			else tooMany = true;
			}
		}
	ptr->nextSet = ptr->nextSet->nextSet;  /* replace the next set with the next next set */
	SortValues(ptr);
	if (tooMany) break;
	}
fprintf(stderr, "done");
if (tooMany) fprintf(stderr, "\nNOTE: More positions were generated than were combined");
} /* end of CombineIntronList */

/******  void DivergeIntronPattern(GeneData*, int, float, float, float)  *****/
/* */
/* FUNCTION: diverges from progenitor pattern for up to 5 generations, then */
/* combines the results and puts them back in the progenitor's struct */
/* ARGUMENTS: ptr to progenitor, and model parameters */
/* RETURN: none */
/* PROTOTYPE IN: noisgene.h */
/* OTHER DEPENDENCIES: */
/* NOTE: effective use of this function requires that kMaxNumValues be set several */
/* times higher than the number of values in the progenitor set */
/* */
/*****************************************************************************/

void DivergeIntronPattern(GeneData *ptr, int gen, float delFreq, float slideFreq, float slideSD)
{
int i, j;
float temp;
GeneData *intronSet[32], *tempPtr;

fprintf(stderr, "\nDiverging . . . ");

for (i = 0; i < 32; i++) 
	{
	intronSet[i] = GetGeneStruct();
	if (i > 0) intronSet[i - 1]->nextSet = intronSet[i];
	}

*intronSet[0] = *ptr;
intronSet[0]->nextSet = intronSet[1];
/* note that last ptr->nextSet is already set to NULL */

for (i = 0; i < gen; i++)
	{
	temp = pow(2.0, (i + 1.0));
	j = (int) temp - 1;
	while (j >= 0)
		{
		fprintf(stderr, "\ngen %d child %d:", i + 1, j + 1);
		/* arggh!! it took me SOOO long to figure out that I had to reset */
		/* the pointers in the consecutive list */
		tempPtr = intronSet[j]->nextSet;
		*intronSet[j] = *intronSet[(int) j / 2];
		intronSet[j]->nextSet = tempPtr;
		EvolveASet(intronSet[j], delFreq, slideFreq, slideSD);
		j--;
		}
	HoldIt();
	}
	
/* that is, the parent of set 4 is set is (int) 4/2, which is 2. The parent of */
/* set 3 is set 1.  Each time throught the loop, each child set is evolved from its */
/* parent set, starting with the top-most child.   For instance, in gen 5, one starts */
/* with set 31, whose parent is set 15, rather than changing 15 first (which would  */
/* result in set 31 going through two rounds of evolution in one generation */

CombineIntronList(intronSet[0]);
*ptr = *intronSet[0];
ptr->nextSet = NULL;
EmptyGeneList(intronSet[0]);
} /* end of DivergeIntronPattern */

/***********************   void MakeNoisyGenes(void)   ***********************/
/* */
/* FUNCTION: gets parameters-- its just a front end for DivergeIntronPatterns */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: noisgene.h */
/* OTHER DEPENDENCIES: */
/* NOTE: */
/* */
/*****************************************************************************/

void MakeNoisyGenes(void)
{
int gen = 0;
float delFreq=0, slideFreq=0, slideSD=0;

fprintf(stderr, "\nThis routine evolves an intron-containing gene through a dichotomously");
fprintf(stderr, "\nbranching tree with g tiers, where g is the number of generations entered");
fprintf(stderr, "\nby the user.  In each generation, a gene from the previous generation is");
fprintf(stderr, "\nduplicated.  Each intron in each copy has a probability d of being deleted");
fprintf(stderr, "\nand a probability s of sliding to a new position (the standard deviation of");
fprintf(stderr, "\nsliding must be entered by the user).  In the end, there will be 2^g intron-");
fprintf(stderr, "\ncontaining genes.  The introns in these genes will be combined into a single");
fprintf(stderr, "\ngene.  Enter the parameters for sliding and loss below.");
do
	{
	fprintf(stderr, "\n\nNumber of generations (1 <= g <= 5): ");
	scanf( "%d", &gen);
	ClearLine();
	}
while ((gen < 1) || (gen > 5));
do
	{
	fprintf(stderr, "\nProbability of loss per intron per generation (0 <= d < 1): ");
	scanf( "%f", &delFreq);
	ClearLine();
	}
while ((delFreq < 0) || (delFreq >= 1));
do
	{
	fprintf(stderr, "\nProbability of slide per intron per generation (0 <= s <= 1): ");
	scanf( "%f", &slideFreq);
	ClearLine();
	}
while ((slideFreq < 0) || (slideFreq > 1));
if (slideFreq != 0)
	do
		{
		fprintf(stderr, "\nStandard deviation of sliding in bp (1 <= SD <= 50): ");
		scanf( "%f", &slideSD);
		ClearLine();
		}
	while ((slideSD < 1) || (slideSD > 50));

fprintf(stderr, "\ngen %d, del %f, slide %f with SD %f", gen, delFreq, slideFreq, slideSD);
HoldIt();
DivergeIntronPattern(gGeneDataPtr[0], gen, delFreq, slideFreq, slideSD);


} /* end of MakeNoisyGenes */
