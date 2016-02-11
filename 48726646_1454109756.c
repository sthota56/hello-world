/* genes.c v. 0.63 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "48726637_1454109756.h"
#include "48726648_1454109756.h"
#include "48726649_1454109756.h"
#include "48726651_1454109756.h"
#include "48726658_1454109756.h"
#include "48726645_1454109756.h"
#include "48726647_1454109756.h"

extern GeneData *gGeneDataPtr[];
extern ObsGeneInfo *gCurObsInfo[];
extern RefGeneInfo *gCurRefInfo[];
extern boolean gOkToTestSet[];
extern int gSlideLimit;

/* private function: */

void ApplySlideRule(GeneData *ptr);
boolean IntToExn(GeneData *ptr);

static char sGeneID[4][10] =
	{
	"OBSINTRN ",
	"REFINTRN ",
	"OBSEXONS ",
	"REFEXONS "
	};

static char vGeneID[4][30] =
	{
	"observed intron positions",
	"reference intron positions",
	"inferred ancestral exons",
	"reference exons"
	};

/**********************  GeneData* GetGeneStruct(void)  *********************/
/* */
/* FUNCTION: allocates and initializes new geneData struct */
/* ARGUMENTS:  none */
/* RETURN: pointer to the new struct */
/* PROTOTYPE IN: genes.h */
/* OTHER DEPENDENCIES: genetype.h */
/* NOTE: especially important for two reasons: first, much */
/* code is based on the fact that each list ends with a NULL pointer in the */
/* nextSet field; second, code also makes use of the fact that exon sizes and */
/* intron positions are always non-zero numbers */
/* */
/* syntax is "myNewDataPtr = GetGeneStruct();" */
/* */
/*****************************************************************************/

GeneData* GetGeneStruct(void)
{
int i = 0;
GeneData *newPtr;

newPtr = (GeneData*) Malloc( sizeof( GeneData ));
for (i = 0; i < kMaxNumValues; i++)
   	{
   	newPtr->number[i] = 0;
   	newPtr->score[i] = 0.00;
   	}

newPtr->geneScore = 0.00;
newPtr->stdDev = 0.00;
newPtr->nextSet = NULL;

return ( newPtr );
} /* end of GetGeneStruct */

/****************    void EmptyGeneList(GeneData *freePtr)     ****************/
/* */
/* FUNCTION: deallocates memory of struct or list or structs */
/* ARGUMENTS:  pointer to struct or first member of list of structs */
/* RETURN: none */
/* PROTOTYPE IN: genes.h */
/* OTHER DEPENDENCIES: genetype.h */
/* */
/*****************************************************************************/

void EmptyGeneList(GeneData *firstPtr)
{
GeneData *temp;

while (firstPtr->nextSet != NULL)
	{
	temp = firstPtr->nextSet->nextSet;
	free(firstPtr->nextSet);
	firstPtr->nextSet = temp;
	}
temp = GetGeneStruct();
*firstPtr = *temp;
}

/*****************   int GetNumValues (GeneData *getNumPtr)    ***************/
/* */
/* FUNCTION : counts the number of introns or exons in a set */
/* ARGUMENTS:  pointer to set of introns or exons to be counted */
/* RETURN: answer */
/* PROTOTYPE IN: genes.h */
/* OTHER DEPENDENCIES: genetype.h */
/* NOTE: none; could easily be modified to count the number of sets in */
/* a list */
/* */
/*****************************************************************************/

int GetNumValues (GeneData *getNumPtr)
{
int i = 0;

if (getNumPtr != NULL )
	{
	while (getNumPtr->number[i] != 0)
		{
		i++;
		}
	}
return(i);
}

/*****************  int GetNumCodons (GeneData *getNumPtr)   *****************/
/* */
/* FUNCTION : counts total codons in an exon set */
/* ARGUMENTS: pointer to set of exons to be counted */
/* RETURN: answer */
/* PROTOTYPE IN: genes.h */
/* OTHER DEPENDENCIES: genetype.h */
/* */
/*****************************************************************************/

int GetNumCodons (GeneData *getNumPtr)
{
int i, numCodons = 0;

if (getNumPtr == NULL )
	{
	fprintf( stderr,  "\nNo gene structure data have been entered yet..." );
	}
else
	{
	for (i = 0; i < kMaxNumValues; i++)
		{
		if (getNumPtr->number[i] != 0)
			{
			numCodons += getNumPtr->number[i];
			}
		}
	}
return(numCodons);
}/* end of getnumExons */

/**************** SortList and SortValues(GeneData *sortPtr)   ***************/
/* */
/* FUNCTION : sorts positions lowest to highest; removes duplicates */
/* ARGUMENTS:  pointer to the set to be ordered */
/* RETURN: none */
/* PROTOTYPE IN: genes.h */
/* OTHER DEPENDENCIES: genetype.h */
/* NOTE: Ordered intron positions are necessary for inter-converting intron */
/* positions and exon sizes; also necessary for centrality scoring; also */
/* makes the presentation more sensible */
/* */
/*****************************************************************************/

void SortList(GeneData *listPtr)
{

while (listPtr != NULL)
	{
	SortValues(listPtr);
	listPtr = listPtr->nextSet;
	}
}

void SortValues(GeneData *sortPtr)
{
int i, j, lowest, temp, numValues;

numValues = GetNumValues(sortPtr);

for (i = 0; i < numValues; i++)
	{
	j = i + 1;
	lowest = sortPtr->number[i];
	while (j < numValues)
		{
		if (sortPtr->number[j] <= lowest)
			{
			lowest = sortPtr->number[j];
			temp = sortPtr->number[i];
			if (temp == lowest)
				{
				while (j < numValues - 1)
					{
					sortPtr->number[j] = sortPtr->number[j + 1];
					j++;
					}
				sortPtr->number[j] = 0; /* sets last to zero */
				numValues--;
				j = i + 1;
				}
			else
				{
				sortPtr->number[i] = sortPtr->number[j];
				sortPtr->number[j] = temp;
				}
			}
		j++;
		}
	}
} /* end of SortValues */


/*******************  void EnterObsGeneData(void)  **********************/
/*  */
/* FUNCTION : converts Dibb-Newman positions into ABaCUS positions in memory */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: genes.h */
/* OTHER DEPENDENCIES: util.h, genetype.h */
/*  */
/*****************************************************************************/

void EnterObsGeneData(void)
{
int	i, k, m, codon, phase, numIntrons, numCodons = 0;

EmptyGeneList(gGeneDataPtr[OBSINTRON]);

fprintf( stderr,  "\nUse this routine to enter intron positions for the first time.  The data will");
fprintf( stderr,  "\nbe converted from intron positions in Dibb & Newman form to intron positions");
fprintf( stderr,  "\non a nucleotide scale.  These data will reside in memory and can be saved to");
fprintf( stderr,  "\ndisk for later use.");

fprintf( stderr,  "\n\nEnter the total number of codons in the gene: ");
scanf( "%d", &numCodons );
ClearLine();

fprintf( stderr,  "\nEnter the number of introns in the coding region: ");
scanf( "%d", &numIntrons);
ClearLine();

fprintf( stderr,  "\nEnter the codon and phase of each intron position (in 5' to 3' order)");
fprintf( stderr,  "\nseparating the numbers by a few spaces.\n\n");

for (i = 0; i < numIntrons; i++)
	{
	for (k = 0; k < i; k++) fprintf(stderr, " ");
	fprintf(stderr, "Intron %2d, codon and phase: ", 1 + i);
	scanf ( "%d %d", &codon, &phase);
	ClearLine();

	gGeneDataPtr[0]->number[i] = 3 * (codon - 1) + phase;
	}

/* updating current intron info */

gCurObsInfo[0]->size = 3 * numCodons;
gOkToTestSet[0] = false;
strcpy(gCurObsInfo[0]->sourceFile, "none saved");
gCurObsInfo[0]->numValues = GetNumValues(gGeneDataPtr[0]);
for (m = 0; m < 3; m++)
	{
	gCurObsInfo[0]->relativePhaseFreqs[m] = 0.0;
	}

/* updating info */

fprintf( stderr, "\nThe new intron positions are as follows: ");
WriteGeneData(0, 0, 0);
fprintf(stderr, "\n");
WriteGeneData(0, 1, 0);

} /* end of EnterObsGeneData */

/*****************   boolean IntToExn(GeneData *ptr)   **********************/
/*  */
/* FUNCTION: converts the intron positions in memory into exon sizes in  */
/* integral numbers of codons  */
/* ARGUMENTS: ptr to source of intron positions */
/* RETURN: whether conversion is successful or not */
/* PROTOTYPE IN: genes.h */
/* OTHER DEPENDENCIES: util.h, genetype.h */
/* NOTE: return(false) if two introns are less than 1 codon apart */
/*  */
/*****************************************************************************/

boolean IntToExn(GeneData *ptr)
{
int 	i, num, temp, leftBnd = 0;
GeneData *tempPtr;
boolean isOK = true;

tempPtr = GetGeneStruct();

num = GetNumValues(ptr);

for (i = 0; i < num; i++)
	{
	if ( fmod(ptr->number[i], 3) < 2) temp = (int) (ptr->number[i]/3) - leftBnd;
	else temp = (int) (ptr->number[i]/3) + 1 - leftBnd;
	if (temp > 0)
		{
		tempPtr->number[i] = temp;
		leftBnd += temp;
		}
	else (isOK = false);
	}

temp = (int) gCurObsInfo[0]->size/3 - leftBnd;
if (temp > 0) tempPtr->number[i] = temp;
else (isOK = false);

if (isOK)
	{
	*gGeneDataPtr[2] = *tempPtr;
	gCurObsInfo[1]->size = (int) gCurObsInfo[0]->size/3;
	gOkToTestSet[1] = false;
	}
return(isOK);
}

/*********************  void ApplySlideRule(GeneData *ptr)  ******************/
/* */
/* FUNCTION: amalgamates adjacent introns <= gSlideLimit apart */
/* ARGUMENTS: ptr to source of intron positions */
/* RETURN: none */
/* PROTOTYPE IN: genes.h */
/* OTHER DEPENDENCIES: util.h, genetype.h */
/* NOTE: change gSlideLimit in settings menu */
/* NOTE: the sliding limit begins with l = 1 bp and the function checks each */
/* neighboring pair of intron positions and combines them if they are <= l bp */
/* apart.  The limit l is increased by 1 and the process is repeated; and so */
/* on until gSlideLimit has been reached.  All positions begin with a weight */
/* 1, but if two positions are combined, the new position has a weight of 2, */
/* and so on.   */
/* */
/*****************************************************************************/

void ApplySlideRule(GeneData *ptr)
{
int i, j, l, num;
float temp;

num = GetNumValues(ptr);

/* pre-weight the positions with w = 1 */

for (i = 0; i < num; i++)
	{
	ptr->score[i] = 1.0;
	}

for (l = 1; l <= gSlideLimit; l++)
	{
	for (i = 1; i < num; i++)
		{
		if ((ptr->number[i] - ptr->number[i - 1]) <= l)
			{
			temp = (((ptr->number[i] * ptr->score[i]) + (ptr->number[i - 1] * ptr->score[i - 1]))/ (ptr->score[i] + ptr->score[i - 1]));
			ptr->number[i - 1] = RoundUp(temp);
			ptr->score[i - 1] += ptr->score[i]; /* i.e., combine weights */
			for (j = i; j < num - 1; j++)
				{
				ptr->number[j] = ptr->number[j + 1];
				ptr->score[j] = ptr->score[j + 1];
				}
			ptr->number[j] = 0; /* deleting last position */
			ptr->score[j] = 0;

			/* prepare for next round: */
			num--; /* decrementing number of positions */
			i--; /* decrementing i, in case successive pairs are <= limit */
			}
		}
	}

} /* end of function ApplySlideRule */

/************************  void InferExons(void)  ***************************/
/* */
/* FUNCTION: infers a set of exons from the observed set of intron positions */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: genes.h */
/* OTHER DEPENDENCIES: util.h, genetype.h */
/* */
/*****************************************************************************/

void InferExons(void)
{
int i,j;
char doSlides=0;
GeneData *tempIntrons;
boolean isOK=0;

tempIntrons = GetGeneStruct();
for (i = 0; i < kMaxNumValues; i++)
	{
	tempIntrons->number[i] = gGeneDataPtr[0]->number[i];
	}


fprintf(stderr, "\nThis routine infers an ancestral set of exons from the observed set of");
fprintf(stderr, "\nintron positions.  These exons become the default 'observed' set in ");
fprintf(stderr, "\nmemory and can be viewed, saved, etc.  A set of inferred ancestral introns");
fprintf(stderr, "\nreplaces the observed set of introns in memory.");
fprintf(stderr, "\nDo you wish to implement 'sliding' [y/n]?\n");
doSlides = GetCommandChar("yn");
if (doSlides == 'y') ApplySlideRule(tempIntrons);
isOK = IntToExn(tempIntrons);

if (isOK)
	{
	*gGeneDataPtr[0] = *tempIntrons;
	for (j = 0; j < 2; j++)
		{
		gCurObsInfo[j]->numValues = GetNumValues(gGeneDataPtr[j]);
		}
	fprintf(stderr, "\nThe inferred set of exon sizes is as follows:");
	WriteGeneData(2, 0, 0);
	fprintf(stderr, "\n");
	WriteGeneData(2, 1, 0);
	fprintf(stderr, "\n(NOTE: if sliding was invoked, the intron set may have changed also--");
	fprintf(stderr, "\nsave these data separately, if desired).");
	}
else
	{
	fprintf(stderr, "\nError!  No ancestral gene could be inferred.  ");
	if (doSlides == 'n')
		{
		fprintf(stderr, "Probably, some neighboring");
		fprintf(stderr, "\nintrons are within 1 codon of each other, and sliding was not invoked.");
		fprintf(stderr, "\nPlease try again, and invoke sliding this time (the original observed");
		fprintf(stderr, "\nset of intron positions has been preserved).  ");
		}
	else
		{
		fprintf(stderr, "It may help to look at the");
		fprintf(stderr, "\nfaulty set of exon sizes now in memory.");
		}
	HoldIt();
	}
} /* end of InferExons */

/***************  void WriteGeneData(int, boolean, boolean)  *****************/
/* */
/* FUNCTION: writes obs or ref introns or exons from memory into file */
/* ARGUMENTS: value of 0-3 for which set to write */
/* RETURN: none */
/* PROTOTYPE IN: genes.h */
/* OTHER DEPENDENCIES: */
/* */
/*****************************************************************************/

void WriteGeneData(int whichSet, boolean toDisk, boolean includeScores)
{
FILE *destPtr;
int i, numValues, j = 0;
char outGeneFile[kMaxNameLength] = "";
GeneData *writeDataPtr = NULL;

writeDataPtr = gGeneDataPtr[whichSet];
numValues = GetNumValues(writeDataPtr);

if (toDisk) 
	{
	destPtr = OpenFile("in which to store the gene data", "a", outGeneFile);
	PrintHeader(destPtr, outGeneFile);
	if (fmod(whichSet, 2)) /* i.e., if its a reference set */
		{
		fprintf(destPtr, "\nNULLTYPE  %s", gCurRefInfo[(int) (whichSet/2)]->refModel);
		}
	if (whichSet == 0)
		{
		fprintf(destPtr, "\nNUMCODON  %d", gCurObsInfo[0]->size/3);
		}
	}
else destPtr = stderr;

while ( writeDataPtr != NULL)
	{
	fprintf(destPtr, "\n%s ", sGeneID[whichSet]);
	j++;
	for (i = 0; i < numValues; i++)
		{
		fprintf(destPtr, "%5d", (writeDataPtr->number[i]));
		}
	if (includeScores)
		{
		fprintf(destPtr, "\nSCORES    ");
		for (i = 0; i < numValues; i++)
			{
			fprintf(destPtr, "%5.1f", (writeDataPtr->score[i]));
			}
		fprintf(destPtr, "\nMEANSCOR  %.4f", writeDataPtr->geneScore);
		fprintf(destPtr, "\tSTDDEV    %.4f", writeDataPtr->stdDev);
		}
	writeDataPtr = writeDataPtr->nextSet;
	}

if (toDisk)
	{
	fprintf(destPtr, "\n\n--- end of this entry ---\n");
	fclose(destPtr);
	if (!fmod(whichSet, 2)) 
		{
		strcpy(gCurObsInfo[(int) (whichSet/2)]->sourceFile, outGeneFile);
		}
	}
	
if (j > 1)
	{
	fprintf( stderr,  "\n\nData from %d sets with %d values each were ", j, numValues);
	if (toDisk) fprintf(stderr, "saved.");
	else fprintf(stderr, "displayed.");
	HoldIt();
	}

}/* end of WriteGeneData */

/******************  void LoadGeneData(int whichSet)  **********************/
/* */
/* FUNCTION: reads obs or ref introns or exons from file into memory */
/* ARGUMENTS: value of 'i' or 'e' for which set to load */
/* RETURN: none */
/* PROTOTYPE IN: genes.h */
/* OTHER DEPENDENCIES: util.h, fileutil.h, genetype.h */
/* NOTE: */
/* */
/*****************************************************************************/

void LoadGeneData(int whichSet)
{
char lineID[10], discard[100], inGeneFile[kMaxNameLength];
int j=0, k, temp, numCodons = 0;
boolean newData = false;
FILE *inGeneFilePtr;

inGeneFilePtr = OpenFile("containing the data to be loaded", "r", inGeneFile);

while (!feof(inGeneFilePtr))
	{
	fscanf(inGeneFilePtr, "%s", lineID);
	if (strstr(lineID, "NUMCODON" ))
		{
		fscanf(inGeneFilePtr, "%d", &numCodons);
		}
	if (!strncmp(lineID, sGeneID[whichSet * 2], 8 ))
		{
		newData = true;
		EmptyGeneList(gGeneDataPtr[whichSet * 2]);
		temp = 0;
		for (j = 0; j < kMaxNumValues; j++)
			{
			fscanf(inGeneFilePtr, " %d", &temp);
			if (temp != 0)
				{
				gGeneDataPtr[whichSet * 2]->number[j] = temp;
				temp = 0;
				}
			else
				{
				break;
				}
			}
		}
	else fgets(discard, 100, inGeneFilePtr);
	strcpy(lineID, "        ");
	}
fclose( inGeneFilePtr );

fprintf( stderr,  "\n\nfile %s has been read.", inGeneFile);

if (newData)
	{
	strcpy(gCurObsInfo[whichSet]->sourceFile, inGeneFile);
	if (numCodons == 0) 
		{
		if (whichSet) numCodons = GetNumCodons(gGeneDataPtr[whichSet * 2]);
		else 
			{		
			fprintf( stderr, "\nEnter the number of codons in the gene: ");
			scanf("%d", &numCodons);
			ClearLine();
			}
		}
	if (whichSet) gCurObsInfo[whichSet]->size = numCodons;
	else gCurObsInfo[whichSet]->size = 3*numCodons;
	gCurObsInfo[whichSet]->numValues = j;
	for (k = 0; k < 3; k++)
		{
		gCurObsInfo[whichSet]->relativePhaseFreqs[k] = 0.0;
		}
	gOkToTestSet[whichSet] = false;
	}
else
	{
	fprintf( stderr,  "\nOops! No %s were entered (check the input file)", vGeneID[whichSet * 2]);
	HoldIt();
	}

}/* end of LoadGeneData */
