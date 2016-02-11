/* dsttobnd.c v. 0.58 */

#include <stdio.h>
#include <string.h>

#include "48726637_1454109756.h"
#include "48726658_1454109756.h"
#include "48726645_1454109756.h"
#include "48726654_1454109756.h"
#include "48726648_1454109756.h"
#include "48726649_1454109756.h"
#include "48726643_1454109756.h"

extern void GetScoreStats(int type); /* this is in abacus.c */
extern int GetNumValues(GeneData*); /* this is in genes.c */

extern ScoreArray *gScoreArray[2];
extern ProtInfo *gCurProtInfo[3];
extern ScoringInfo *gCurScoreInfo[2];
extern boolean gPenalizeEnds;
extern GeneData *gGeneDataPtr[4];

char iERDist[] = "distance score, base pairs to nearest inter-element region";


/************************ void EnterElements(void)  **************************/
/* */
/* FUNCTION: creates a 0,1 array representing user-defined elements */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: dsttobnd.h */
/* OTHER DEPENDENCIES: util.h, prottype.h */
/* */
/* NOTE: the relation between codon position, intron position, and scoring */
/* matrix index.  An intron position is numbered with the nucleotide 5' to it, */
/* as explained in the function for entering observed gene data: */
/* */
/* intron position = 3(codon - 1) + phase; */
/* */
/* The number of possible intron positions is the gene length - 1, since */
/* an intron can follow any bp in the gene except the last one.  The scoring */
/* matrix indices start at 0 and go to gene length - 2.  Thus, the matrix */
/* element index for position 1 is 0, and, in general, */
/* */
/* matrix index for intron position: index i = p - 1; */
/* */
/* Finally, the protein structure elements are given as residue numbers, which */
/* correspond to codon numbers.  The first intron position implicated in the */
/* structural element is phase 1 of codon leftRes, which by the formulae above is */
/* */
/* left intron position = 3(leftRes - 1) + 1; */
/* left matrix element = left intron position - 1 = 3 (leftRes - 1); */
/* */
/* and the last intron in the structural element is in the second phase of */
/* codon rightRes, thus: */
/* */
/* right intron position = 3(rightRes - 1) + 2; */
/* right matrix element = right intron position - 1 = 3(rightRes - 1) + 1; */
/* */
/*****************************************************************************/

void EnterElements(void)
{
int i, left, right, numCodons;
char	response=0;

/* re-initialize array with 0's */

for (i = 0; i < kMaxArraySize; i++)
	{
	gScoreArray[0]->iScore[i] = 0;
	}
	
fprintf( stderr,  "\nEnter the number of residues in the protein: ");
scanf( "%d", &numCodons);
ClearLine();

fprintf( stderr, "\nNow enter the boundaries of each helix, sheet, turn, or other element,");
fprintf( stderr, "\none element at a time. \n\n");

while (response != 'd')
	{
	fprintf( stderr, "First and last residues of element (separated by spaces): ");
	scanf("%d %d", &left, &right );
	ClearLine();
	left = 3 * (left - 1);
	right = 3 * (right - 1) + 1;
	for (i = left; i <= right; i++)
		{
		gScoreArray[0]->iScore[i] = 1;
		}
	fprintf( stderr, "\nEnter another? (c=continue; d=done): ");
	response = getchar();
	ClearLine();
	}

/* converting the array to graduated form: */

ConvertArray(999);

/* updating current info */

gCurProtInfo[0]->size = (3 * numCodons) - 1;
strcpy(gCurProtInfo[0]->sourceFile, "none saved");

fprintf(stderr, "\nEnter a descriptive comment to be stored with these scores (<80 characters).\nCOMMENT: ");
gets(gCurProtInfo[0]->comment);

} /* end of EnterElements  */

/************************ void WriteDistArray(boolean toDisk) *************************/
/* */
/* FUNCTION: writes array to disk or screen */
/* ARGUMENTS: disk=1, screen=0 */
/* RETURN: none */
/* PROTOTYPE IN: dsttobnd.h */
/* OTHER DEPENDENCIES: util.h, prottype.h */
/* */

void WriteDistArray(boolean toDisk)
{
char outFile[kMaxNameLength];
int i;
FILE	*filePtr;

if (toDisk)
	{
	filePtr = OpenFile("to contain the distance scoring array", "w", outFile);
	PrintHeader(filePtr, outFile);
	}
else
	{
	filePtr = stderr;
	PrintHeader(filePtr, NULL);
	}

fprintf( filePtr, "\nCOMMENT  %s", gCurProtInfo[0]->comment);
fprintf( filePtr, "\nPROTFILE  %s", gCurProtInfo[0]->sourceFile);
fprintf( filePtr, "\nNUMCODON  %d", (int) (gCurProtInfo[0]->size + 1)/3);
fprintf( filePtr, "\nMAXSCORE  %d", gCurProtInfo[0]->param1);
fprintf( filePtr, "\nAVESCORE  %f", gCurProtInfo[0]->ave);
fprintf( filePtr, "\nSTDDEV    %f", gCurProtInfo[0]->stdDev);
fprintf( filePtr, "\nSCOARRAY   ");
if (!toDisk) HoldIt();
for (i = 0; i < gCurProtInfo[0]->size; i++)
	{
	fprintf( filePtr, "%d ", gScoreArray[0]->iScore[i]);
	}
if (toDisk)
	{
	fprintf( filePtr, "\n\nend of file\n\n");
	fclose(filePtr);
	fprintf( stderr, "\n%d scores written to disk in file %s.", gCurProtInfo[0]->size, outFile);
	strcpy(gCurProtInfo[0]->sourceFile, outFile);
	}
HoldIt();

} /* end WriteDistArray */

/*******************  void ConvertArray(int newMax)  ************************/
/* */
/* FUNCTION: converts global array to new maxValue */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: dsttobnd.h */
/* OTHER DEPENDENCIES: util.h, prottype.h */
/* NOTE: */
/* */
/*****************************************************************************/

void	ConvertArray(int newMax)
{
int i, j;

if (newMax == 0)
	{
	fprintf( stderr, "\nCONVERT ARRAY (converts the current array to a new maximum score)");
	fprintf( stderr, "\n\nEnter an integer for the new maximum possible score: ");
	scanf("%d", &newMax);
	ClearLine();
	}

/* converting the array into an array of 1's and 0's: */

for (j = 0; j < gCurProtInfo[0]->size; j++)
	{
	if (gScoreArray[0]->iScore[j] > 0)
		{
		gScoreArray[0]->iScore[j] = 1;
		}
	}

/* converting the array of 1's and 0's into a graded array:                */
/* the converter is very simple.  If the value of an element and both      */
/* flanking values are greater than j, then the value will be incremented; */
/* otherwise it will stay the same.  The values at the extreme ends of the */
/* array are left unchanged.  The result is a graded array, where the      */
/* original 0's are still there, but the 1's have been incremented if they */
/* are not close to a 0.  */

for (j = 0; j < (newMax - 1); j++)
	{
	if (gPenalizeEnds)
		{
		if ((gScoreArray[0]->iScore[0] > j) && (gScoreArray[0]->iScore[1] > j))
			{
			gScoreArray[0]->iScore[0] += 1;
			}
		}
	for (i = 1; i < (gCurProtInfo[0]->size - 1); i++)
		{
		if ((gScoreArray[0]->iScore[i] > j) && (gScoreArray[0]->iScore[i-1] > j) && (gScoreArray[0]->iScore[i+1] > j))
			{
			gScoreArray[0]->iScore[i] += 1;
			}
		}
	if (gPenalizeEnds)
		{
		if ((gScoreArray[0]->iScore[gCurProtInfo[0]->size - 1] > j) && (gScoreArray[0]->iScore[gCurProtInfo[0]->size - 2] > j))
			{
			gScoreArray[0]->iScore[gCurProtInfo[0]->size - 1] += 1;
			}
		}
	}

/* updating current info */

GetScoreStats(0);
gCurProtInfo[0]->param1 = newMax;

} /* end ConvertArray */

/*********************  void LoadArray(void)  ****************************/
/* */
/* FUNCTION: reads file contents into global array */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: dsttobnd.h */
/* OTHER DEPENDENCIES: util.h, prottype.h */
/* NOTE: */
/* */
/*****************************************************************************/

void LoadArray(void)
{
char	inArrayFile[kMaxNameLength], lineID[10], discard[200];
int i, temp, numCodons = 0;
FILE *arrayFilePtr;

/* initialize global array */

for (i = 0; i < kMaxArraySize; i++)
	{
	gScoreArray[0]->iScore[i] = 0;
	}

fprintf( stderr, "\nLOAD ARRAY (loads an array representing structural elements from disk)\n\n");
arrayFilePtr = OpenFile("containing the array to be loaded", "r", inArrayFile);

while (!feof(arrayFilePtr))
	{
	fscanf(arrayFilePtr, "%s", lineID);
	fprintf( stderr,  " %s . . .", lineID);

	if (strstr(lineID, "NUMCODON") != NULL)
		{
		fscanf(arrayFilePtr, " %d", &numCodons);
		gCurProtInfo[0]->size = (int) ( 3 * numCodons ) - 1;
		}
	if (strstr(lineID, "MAXSCORE") != NULL)
		{
		fscanf(arrayFilePtr, " %d", &gCurProtInfo[0]->param1);
		}
	if (strstr(lineID, "AVESCORE") != NULL)
		{
		fscanf(arrayFilePtr, " %f", &gCurProtInfo[0]->ave);
		}
	if (strstr(lineID, "STDDEV") != NULL)
		{
		fscanf(arrayFilePtr, " %f", &gCurProtInfo[0]->stdDev);
		}
	if (strstr(lineID, "SCOARRAY") != NULL)
		{
		for (i = 0; i < gCurProtInfo[0]->size; i++)
			{
			fscanf(arrayFilePtr, "%d", &temp);
			gScoreArray[0]->iScore[i] = temp;
			fprintf( stderr, "%d ", temp);
			temp = 0;
			}
		}
	fgets(discard, 200, arrayFilePtr);
	}
fclose(arrayFilePtr);

fprintf( stderr, "\n\nThe array from '%s' has been loaded into memory.  ", inArrayFile);

if (numCodons > 0)
	{
	fprintf( stderr, "The array represents\na gene of %d codons.", numCodons);
	}
else
	{
	fprintf( stderr,  "The number of codons\nrepresented by the array could not be determined.  Please enter the number");
	fprintf( stderr, "\nof codons now: ");
	scanf( "%d", &numCodons);
	ClearLine();
	}

/* updating current info */

strcpy(gCurProtInfo[0]->sourceFile, inArrayFile);
	
HoldIt();

} /* end of LoadArray() */

/**********************  void AssignArrayScore(void)  ****************************/
/* */
/* FUNCTION: assigns scores to intron positions using a scoring array */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: dsttobnd.h */
/* OTHER DEPENDENCIES: util.h, genes.h, genetype.h, prottype.h */
/* NOTE: */
/* */
/*****************************************************************************/

void AssignArrayScore(void)
{
int i, position, numIntrons;
GeneData *intronPtr;

numIntrons = GetNumValues(gGeneDataPtr[0]);

fprintf( stderr, "\nASSIGN SCORES (assigns scores to introns using the current array in memory)\n\n");

/* note: the array element position - 1 refers to the position */

intronPtr = gGeneDataPtr[0];
while (intronPtr != NULL)
	{
	for ( i = 0; i < numIntrons; i++)
		{
		position = intronPtr->number[i];
		intronPtr->score[i] = (float) gScoreArray[0]->iScore[position - 1];
		} 
	if (intronPtr != gGeneDataPtr[0]) intronPtr = intronPtr->nextSet;
	else intronPtr = gGeneDataPtr[1];
		
	fprintf( stderr, ". ");
	}

/* updating current info */

gCurScoreInfo[0]->whichType = 'a';
strcpy(gCurScoreInfo[0]->scoringRule, iERDist);
gCurScoreInfo[0]->param1 = 0;
gCurScoreInfo[0]->param2 = 0.00;
gCurScoreInfo[0]->weightAverage = false;

fprintf( stderr, "\nThe observed and any null introns have now been scored using the array.");
HoldIt();

} /* end of AssignArrayScore */
