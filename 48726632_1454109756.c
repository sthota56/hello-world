/* abacus.c main code block v. 0.62f, 26 Sept 1995 */

#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "48726637_1454109756.h"
#include "48726648_1454109756.h"
#include "48726654_1454109756.h"
#include "48726649_1454109756.h"
#include "48726658_1454109756.h"
#include "48726645_1454109756.h"
#include "48726647_1454109756.h"
#include "48726656_1454109756.h"
#include "48726651_1454109756.h"
#include "48726641_1454109756.h"
#include "48726653_1454109756.h"
#include "48726643_1454109756.h"
#include "48726635_1454109756.h"


#ifdef DOS
extern unsigned _stklen = 8000U;   /* Doubled size as fix for stack overflow*/
#endif

/******************************    globals      ******************************/
/* */
/* -- gene data are stored in 4 singly linked lists of GeneData structs */
/* -- CA coordinates are stored in a singly linked list of XYZCoord structs */
/* -- residue scores are stored in a 1-dimensional float array */
/* -- distance-to-boundary scores are stored in a 1-dimensional int array */
/*    (the two types of arrays are unionized to allow some common functions) */
/* */
/* -- for each type of data, there is an 'info' struct for current info */
/* */
/* -- a singly linked list of ExperimentInfo structs is maintained, each with */
/*    pointers to the pertinent gene and protein info, as well as other information */
/* -- the OkToTestSet booleans are used to determine whether the observed and reference */
/*    sets are comparable, so that a test may be performed */
/* */
/* -- several 'settings' globals are used */
/* */
/*****************************************************************************/


GeneData *gGeneDataPtr[4]={NULL,NULL,NULL,NULL}; /* for obs & ref introns & exons */

XYZCoord *gFirstResPtr=NULL;
ScoreArray *gScoreArray[2]={NULL,NULL};

ObsGeneInfo *gCurObsInfo[2]={NULL,NULL};  /* 0=intron, 1=exon */
RefGeneInfo *gCurRefInfo[2]={NULL,NULL};   /* ditto */
ScoringInfo *gCurScoreInfo[2]={NULL,NULL};  /* ditto */
ProtInfo *gCurProtInfo[3]={NULL,NULL,NULL};
ExptInfo *gFirstExptInfoPtr=NULL, *gCurExptInfoPtr=NULL;
boolean gOkToTestSet[2] = {false, false};


long gSeed; /* see randum.c for more information */

/* settings for optional saving of information: */

boolean gOutputAASeq = false;
int gOutputRefGenes = 0, gOutputRefScores = 0;

/* settings that affect hypothesis-testing (explained in settings menu): */

int gSlideLimit = 9, gWeightAverage = 0;
boolean gPenalizeEnds = false;
char gPhaseZero='r';



/******************  ExptInfo* GetExptInfoStruct(void)  **********************/
/* */
/* FUNCTION: allocates and initializes new struct for the experiment list */
/* ARGUMENTS:  none */
/* RETURN: pointer to the new struct */
/* PROTOTYPE IN: abacus.h */
/* DEPENDENCIES: Malloc in util.c; ExptInfo struct defined in infotype.h */
/* NOTE:  */
/* */
/*****************************************************************************/

ExptInfo* GetExptInfoStruct(void)
{
ExptInfo *newPtr;

newPtr = (ExptInfo*) Malloc( sizeof( ExptInfo ));

newPtr->exptNum = 0;
strcpy(newPtr->comment, "not available");
newPtr->whichType = 0;
newPtr->obsPtr = NULL;
newPtr->refPtr = NULL;
newPtr->protPtr = NULL;
newPtr->scorPtr = NULL;
newPtr->obsGeneScore = 0.00;
newPtr->refMean = 0.00;
newPtr->meanRefSD = 0.00;
newPtr->sDRefMean = 0.00;
newPtr->pValue = 0.00;
newPtr->next = NULL;

return ( newPtr );
} /* end of GetExptInfoStruct */

/******************  ScoreArray *InitArray(char utype)  ***********************/
/* */
/* FUNCTION: intializes one-dimensional array of ints or floats with 0's */
/* ARGUMENTS:  char 'i' or 'f' to indicate ints or floats */
/* RETURN: pointer to array */
/* OTHER DEPENDENCIES: Malloc in util.c; union of float or int arrays defined in prottype.h; */
/* */
/*****************************************************************************/

ScoreArray *InitArray(char utype)
{
ScoreArray *newPtr;
int i;

newPtr = (ScoreArray*) Malloc(sizeof(ScoreArray));
if (utype == 'i')
	{
	for (i = 0; i < kMaxArraySize; i++)
		{
		newPtr->iScore[i] = 0;
		}
	}
else
	{
	for (i = 0; i < ((int) kMaxArraySize/3); i++)
		{
		newPtr->fScore[i] = 0.00;
		}
	}
return ( newPtr );

} /* end InitArray */


/**********************  void InitGlobals(void)  *****************************/
/* */
/* FUNCTION: initializes current info to startup values */
/* ARGUMENTS:  none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h  */
/* OTHER DEPENDENCIES: genetype.h, infotype.h, prottype.h, util.h, genes.h */
/* */
/*****************************************************************************/

void InitGlobals(void)
{
int i,j;

gSeed = InitRanSeed();

/** intializing data & expt pointers: **/

for (i = 0; i < 4; i++) gGeneDataPtr[i] = GetGeneStruct();

gScoreArray[0] = InitArray('i');
gScoreArray[1] = InitArray('f');

gFirstResPtr = GetXYZCoordStruct();


for (i = 0; i < 2; i++)
	{
	gCurObsInfo[i] = (ObsGeneInfo*) Malloc(sizeof(ObsGeneInfo));
	strcpy(gCurObsInfo[i]->sourceFile, "none saved");
	gCurObsInfo[i]->size = 0;
	gCurObsInfo[i]->numValues = 0;
	for (j = 0; j < 3; j++)
		{
		gCurObsInfo[i]->relativePhaseFreqs[j] = (float) 0.0;
		}

	gCurRefInfo[i] = (RefGeneInfo*) Malloc(sizeof(RefGeneInfo));
	strcpy(gCurRefInfo[i]->refModel, "not available");
	strcpy(gCurRefInfo[i]->sourceFile, "none saved");
	gCurRefInfo[i]->numSets = 0;
	gCurRefInfo[i]->lowerLimit = 0;
	gCurRefInfo[i]->hasPhaseBias = false;
	gCurRefInfo[i]->hasPolarity = false;

	gCurScoreInfo[i] = (ScoringInfo*) Malloc(sizeof(ScoringInfo));
	gCurScoreInfo[i]->whichType = 0;
	strcpy(gCurScoreInfo[i]->scoringRule, "[none available]");
	gCurScoreInfo[i]->param1 = 0;
	gCurScoreInfo[i]->param2 = 0.00;
	gCurScoreInfo[i]->weightAverage = false;
	}

for (i = 0; i < 3; i++)
	{
	gCurProtInfo[i] = (ProtInfo*) Malloc(sizeof(ProtInfo));
	strcpy(gCurProtInfo[i]->sourceFile,"none saved");
	strcpy(gCurProtInfo[i]->comment, "none");
	gCurProtInfo[i]->size = 0;
	gCurProtInfo[i]->param1 = 0;
	gCurProtInfo[i]->ave = 0.00;
	gCurProtInfo[i]->stdDev = 0.00;
	}

/* separate initializers: */

gFirstExptInfoPtr = GetExptInfoStruct();
gFirstExptInfoPtr->exptNum = 1;
gCurExptInfoPtr = gFirstExptInfoPtr;


} /* end of InitGlobals */

/********************** ExptInfo* NextExptInfoPtr(void) **********************/
/* */
/* FUNCTION: gets new expt info struct and adds to list */
/* ARGUMENTS:  none */
/* RETURN: pointer to the new struct */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: struct defined in infotype.h */
/* NOTE: syntax is gCurExptInfoPtr = NextExptInfoPtr(); */
/* */
/*****************************************************************************/

ExptInfo* NextExptInfoPtr(void)
{
ExptInfo *nextExptPtr;

/* getting struct for the next experiment and filling in data: */

nextExptPtr = GetExptInfoStruct();
nextExptPtr->exptNum = 1 + gCurExptInfoPtr->exptNum;
gCurExptInfoPtr->next = nextExptPtr;

return (nextExptPtr);
} /* end of NextExperimentInfoPtr */


/********************   void ShowCurrentInfo(void)    ************************/
/* */
/* FUNCTION: gives info on data currently in memory */
/* ARGUMENTS:  none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* CALLS: PrintHeader */
/* OTHER DEPENDENCIES: util.h, infotype.h, fileutil.h, genes.h */
/* */
/*****************************************************************************/

void ShowCurrentInfo(void)
{
int i;
char geneID[2][20] = { "INTRON POSITIONS", "EXON SIZES" };

PrintHeader(stderr, NULL);

for (i = 0; i < 2; i++)
	{
	fprintf( stderr, "\n\n%s: \n  Observed: ", geneID[i]);
	if (GetNumValues(gGeneDataPtr[2*i]))
		{
		if (i == 0)
			{
			fprintf(stderr, "%d positions in gene of %d nt, from file \"%s\"", gCurObsInfo[i]->numValues, gCurObsInfo[i]->size, gCurObsInfo[i]->sourceFile);
			}
		else fprintf(stderr, "%d sizes in gene of %d codons, from file \"%s\"", gCurObsInfo[i]->numValues, gCurObsInfo[i]->size, gCurObsInfo[i]->sourceFile);

		fprintf(stderr, "\n  Reference: ");
		if (GetNumValues(gGeneDataPtr[1 + 2*i]) != 0)
			{
			fprintf(stderr, "%d sets %s", gCurRefInfo[i]->numSets, gCurRefInfo[i]->refModel );
			}
		else fprintf(stderr, "[none have been generated]");

		fprintf(stderr, "\n  Correspondence scores: %s", gCurScoreInfo[i]->scoringRule);
		if (gOkToTestSet[i]) fprintf(stderr, "\n    ==> ready to test reference hypothesis");
		}
	else
		{
		fprintf(stderr, " [none available]");
		}
	}

fprintf(stderr, "\n\nPROTEIN DATA FOR STRUCTURE-BASED SCORING: ");
fprintf(stderr, "\n  Distance-to-boundary scores: ");
if (gCurProtInfo[0]->size)
	{
	fprintf( stderr, "%d inter-bp sites (ave = %.1f; SD = %.1f) ",
	gCurProtInfo[0]->size, gCurProtInfo[0]->ave, gCurProtInfo[0]->stdDev );
	fprintf( stderr, "from file \"%s\"", gCurProtInfo[0]->sourceFile);
	}
else fprintf(stderr, "[none available] ");

fprintf(stderr, "\n\n  Residue scores: ");
if (gCurProtInfo[1]->size)
	{
	fprintf(stderr, "%d scores (ave = %.1f; SD = %.1f) from file \"%s\"", gCurProtInfo[1]->size,
	gCurProtInfo[1]->ave, gCurProtInfo[1]->stdDev, gCurProtInfo[1]->sourceFile);
	}
else fprintf(stderr, "[none available]");

fprintf(stderr, "\n\n  CA coordinates: ");
if (gCurProtInfo[2]->size)
	{
	fprintf(stderr, "%d CA coordinates from file \"%s\"", gCurProtInfo[2]->size, gCurProtInfo[2]->sourceFile);
	}
else fprintf(stderr, "[none available]");

fprintf(stderr, "\n\nEXPERIMENT LIST: ");
if (gCurExptInfoPtr->exptNum > 1)
	{
	fprintf(stderr, "results of %d test(s) in memory", gCurExptInfoPtr->exptNum - 1);
	}
else fprintf(stderr, "[no experiments have been performed yet]");

HoldIt();

} /* end of ShowCurrentInfo */


/********  void WriteExptInfo(ExptInfo *exptPtr, FILE *destPtr) **************/
/* */
/* FUNCTION: formats expt info for writing to screen or file */
/* ARGUMENTS:  pointers to info and destination */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: fileutil.h, infotype.h */
/* NOTE: to write to console, call with WriteExptInfo(myPtr, stderr) */
/* */
/*****************************************************************************/

void WriteExptInfo(ExptInfo *exptPtr, FILE *destPtr)
{
char geneType[20];

if (exptPtr->whichType == 'e') strcpy(geneType, "exon size");
else	strcpy(geneType, "intron position");

fprintf( destPtr, "\nEXPERIMENT %d: %s", exptPtr->exptNum, exptPtr->comment);
fprintf( destPtr, "\nGENE DATA: %d %ss from file \"%s\"", exptPtr->obsPtr->numValues, geneType, exptPtr->obsPtr->sourceFile);
fprintf( destPtr, "\nSIMULATIONS: %d sets %s", exptPtr->refPtr->numSets, exptPtr->refPtr->refModel);
fprintf( destPtr, "\nSCORING: ");
fprintf( destPtr, "%s", exptPtr->scorPtr->scoringRule);
if (exptPtr->whichType == 'a' || exptPtr->whichType == 'c')
	{
	fprintf( destPtr, "\n  from the set of scores (ave = %.2f; SD = %.2f) in \"%s\" ", exptPtr->protPtr->ave, exptPtr->protPtr->stdDev, exptPtr->protPtr->sourceFile);
	}
fprintf( destPtr, "\nMEAN SCORES:");
fprintf( destPtr, "\n   Observed: %f", exptPtr->obsGeneScore);
fprintf( destPtr, "\n   Reference: %f (SD = %f)", exptPtr->refMean, exptPtr->sDRefMean);
fprintf( destPtr, "\nRANK-ORDER PROBABILITY: %f", exptPtr->pValue);

} /* end of WriteExptInfo */

/***************   void WriteExptList(boolean toDisk)  ***********************/
/* */
/* FUNCTION: writes entire list of all complete expts to file or console */
/* ARGUMENTS:  whether to write to file or not */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: infotype.h, util.h, fileutil.h */
/* NOTE: */
/* */
/*****************************************************************************/

void WriteExptList(boolean toDisk)
{
char sumFileName[kMaxNameLength] = "";
FILE *destPtr = stderr;
ExptInfo *exptPtr;

exptPtr = gFirstExptInfoPtr;

if (toDisk)
	{
	destPtr = OpenFile("to contain summaries of all experiments\npresently in memory", "w", sumFileName);
	PrintHeader(destPtr, sumFileName);
	}
else PrintHeader(destPtr, NULL);


while (exptPtr->next != NULL) /* only writes tests, not untested info */
	{
	WriteExptInfo(exptPtr, destPtr);
	if (toDisk) fprintf(destPtr, "\n\n\n");
	else HoldIt();
	exptPtr = exptPtr->next;
	}
} /* end WriteExptList */

/*********************   void ShowStartupScreen(void)  ***********************/
/* */
/* FUNCTION: print startup info to console */
/* ARGUMENTS:  none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: util.h */
/* NOTE: maximum allowable number of introns/exons is maxNumValues-1 */
/* */
/*****************************************************************************/

void ShowStartupScreen(void)
{

fprintf( stderr,  "\n\n\n\t\t\t\t%s", kVersion );
fprintf( stderr,  "\n\t   Analysis of Blake's Conjecture Using Simulations");
fprintf( stderr,  "\n\n\t\t   Arlin Stoltzfus and David Spencer");
fprintf( stderr,  "\n\tDalhousie University, Halifax, Nova Scotia, CANADA B3H 4H7");

/**  reminder about #defined size limits and memory overhead **/

fprintf( stderr, "\n\n\n\t\t     Some reminders about limits:");

fprintf( stderr, "\n\n\t-maximum allowable number of exons or introns per set is %d", kMaxNumValues-1);
fprintf( stderr, "\n\t-maximum size of gene (or protein) is %d codons (or residues)", (int) kMaxArraySize/3);
fprintf( stderr, "\n\t-simulated genes will use %d Kbytes of memory per 1000", (int) sizeof(GeneData));

HoldIt();

} /* end ShowStartupScreen */


/*****************   long GetNumPairs( int numResidues )    ******************/
/* */
/* FUNCTION: number of pairwise comparisons of input number of residues */
/* ARGUMENTS: number of residues (viz., in an exon-encoded peptide) */
/* RETURN: answer */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES:  none */
/* NOTE: the number of distinct pairs of residues encoded by an exon is */
/* used in calculating the average inter-residue distance used in distance- */
/* based scoring */
/* */
/*****************************************************************************/

long GetNumPairs( int numResidues )
{
int i;
long numPairs = 0;

for (i = 1; i < numResidues; i++)
	{
	numPairs += i;
	}
return(numPairs);
}/** end getnumpairs **/

/***********  int ChooseExons(void), int ChooseRefs(void) ********************/
/* */
/* FUNCTION: prompts user to choose */
/* ARGUMENTS:  none */
/* RETURN: answer = int 0 or 1 */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: util.h */
/* */
/*****************************************************************************/

int ChooseExons(void)
{
int isExon=0;

fprintf( stderr,  "i=intron positions or e=exon sizes?");
if (GetCommandChar("ie") == 'e') isExon = 1;

return(isExon);
}

int ChooseRefs(void)
{
int isRef=0;

fprintf( stderr,  "o=observed data or r=reference data?");
if (GetCommandChar("or") == 'r') isRef = 1;

return(isRef);
}

/*******************  void GetStats(GeneData *statsPtr)  *********************/
/* */
/* FUNCTION: calculates means and stdDevs for use by EvaluateRef */
/* ARGUMENTS: pointer to GeneData struct */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: genetype.h, genes.h */
/* NOTE: implements option of weigting exon scores by exon lengths; this is */
/* useful for evaluating the effects of the size distribution on extensity */
/* scores */
/* */
/*****************************************************************************/

void GetStats(GeneData *statsPtr)
{
int i, numCodons, numValues;
float sumSqrdDev, sumScores, curMean, curStdDev;

numValues = GetNumValues(statsPtr);

if ((gWeightAverage) && ((statsPtr == gGeneDataPtr[2]) || (statsPtr == gGeneDataPtr[3])))
	{
	numCodons = GetNumCodons(statsPtr);
	while (statsPtr != NULL)
		{
		sumScores = 0.00;
		sumSqrdDev = 0.00;
		for (i = 0; i < numValues; i++)
			{
			sumScores += statsPtr->score[i] * statsPtr->number[i];
			}
		curMean = (float) sumScores / numCodons;
		statsPtr->geneScore = curMean;
		for (i = 0; i < numValues; i++)
			{
			sumSqrdDev += SquareIt( statsPtr->score[i] - curMean);
			}
		curStdDev = sqrt(sumSqrdDev / (numValues - 1));
		statsPtr->stdDev = curStdDev;
		statsPtr = statsPtr->nextSet;
		}
	}
else
	{
	while (statsPtr != NULL)
		{
		sumScores = 0.00;
		sumSqrdDev = 0.00;
		for (i = 0; i < numValues; i++)
			{
			sumScores += statsPtr->score[i];
			}
		curMean = (float) sumScores / numValues;
		statsPtr->geneScore = curMean;
		for (i = 0; i < numValues; i++)
			{
			sumSqrdDev += SquareIt( statsPtr->score[i] - curMean);
			}
		if (numValues > 1)
			{
			curStdDev = sqrt(sumSqrdDev / (numValues - 1));
			statsPtr->stdDev = curStdDev;
			}
		else statsPtr->stdDev = 0.00;
		statsPtr = statsPtr->nextSet;
		}
	}
} /* end GetStats */

/*******************    void GetScoreStats(int type)   *****************/
/* */
/* FUNCTION: calculates means and stdDevs for distance and residue scoring */
/* ARGUMENTS: which global array to use */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: genetype.h, genes.h, mymath.h, prottype.h */
/* NOTE: */
/* */
/*****************************************************************************/

void GetScoreStats(int type)   /* 0=dist-to-boundary array, 1=residues array */
{
int i;
float ave=0.00, stdDev=0.00;

if (type == 0)
	{
	for (i = 0; i < gCurProtInfo[type]->size; i++)
		{
		ave += (float) gScoreArray[type]->iScore[i]/gCurProtInfo[type]->size;
		}
	for (i = 0; i < gCurProtInfo[type]->size; i++)
		{
		stdDev += SquareIt(gScoreArray[type]->iScore[i] - ave)/gCurProtInfo[type]->size;
		}
	}
else
	{
	for  (i = 0; i < gCurProtInfo[type]->size; i++)
		{
		ave += (float) gScoreArray[type]->fScore[i]/gCurProtInfo[type]->size;
		}
	for  (i = 0; i < gCurProtInfo[type]->size; i++)
		{
		stdDev += SquareIt(gScoreArray[type]->fScore[i] - ave)/gCurProtInfo[type]->size;
		}
	}

stdDev = sqrt(stdDev);
gCurProtInfo[type]->stdDev = stdDev;
gCurProtInfo[type]->ave = ave;

}  /* end of GetScoreStats */

/*************   int PermuteValues(GeneData*, GeneData**)    *****************/
/* */
/* FUNCTION: creates a list of permuted sets of input set of values */
/* ARGUMENTS: input set of values, (address for) number of sets to generate */
/* RETURN: pointer to first member of list */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: randum.h, genes.h, genetype.h  */
/* NOTE: algorithm is explained in detail below */
/* */
/*****************************************************************************/

int PermuteValues(GeneData *sourcePtr, GeneData **destHandle)
{
GeneData *curPermPtr=NULL, *nextPermPtr=NULL, *tempSourcePtr=NULL;
int temp, k, p, f;
int numValues, numSets;
double j;

tempSourcePtr = GetGeneStruct();
*tempSourcePtr = *sourcePtr;

/* preparing to prompt user for number of sets: */

numValues = GetNumValues(tempSourcePtr);
j = Factorial(numValues);

if ( j <= 10000L)
	{
	fprintf( stderr, "\n\nWarning: there are only %.0f possible orders of %d non-identical sizes.", j, numValues);
	fprintf( stderr, "\nIf you choose to use a large proportion of this number, the reliability");
	fprintf( stderr, "\nof the resulting P value cannot be estimated easily.  Specifically, the P");
	fprintf( stderr, "\nvalue will have excessive variance because non-independence of scores");
	fprintf( stderr, "\nobtains when identical (and non-identical, but similar) exon orders are");
	fprintf( stderr, "\ngenerated.  ");
	if ( j <= 500L )
		{
		fprintf( stderr, "In this case, the number of orders is so low that a model of");
		fprintf( stderr, "this type is not recommended.");
		}
	}
else
	{
	fprintf( stderr, "\nThere are %.0f possible orders of %d non-identical sizes", j, numValues);
	}
fprintf( stderr, "\n\nEnter the number of sets to be generated by permutation: ");
scanf("%d", &numSets);
ClearLine();

/* setting pointer for new sets */

nextPermPtr = *destHandle;
k = 0;
while ((k < numSets) && (nextPermPtr != NULL))
	{
	curPermPtr = nextPermPtr;

	f = (numValues - 1);
	while (f >= 0)
		{
		p = (int) (f + 1) * NRuran2();
		temp = tempSourcePtr->number[p];
		curPermPtr->number[f] = temp;
		while (p < f)
			{
			tempSourcePtr->number[p] = tempSourcePtr->number[p + 1];
			p++;
			}
		tempSourcePtr->number[f] = temp;
		f--;
		}

/* p = picked from source, f=filled in destination. First the loop fills in  */
/* the last number in the reference array, number[f = numValues-1], giving it */
/* the value of one of the numbers in the observed set. The random pick      */
/* is made via the array index, an integer from 0 to f.  The while loop      */
/* bumps down the index of each subsequent value in the temp obs array,      */
/* then re-inserts the picked value held in temp at the end of the array.    */
/* When the next pick is made, the previous pick can't be picked again,      */
/* because the range of picked values decreases every time f is decremented  */
/* This repeats until f = 0, in                                              */
/* which case the random number function will always return p = 0, and the   */
/* the last unpicked obs value is picked and put into the remaining unfilled */
/* element of the reference array.  In the end, the temp array has the order */
/* of the reference array just created: thus, each shuffle is made from the  */
/* previous one, rather than being generated from the same order each time   */

	nextPermPtr = GetGeneStruct();
	curPermPtr->nextSet = nextPermPtr;
	k++;
	fprintf( stderr,  "\nset %d completed", k);
	}
curPermPtr->nextSet = NULL;
free(nextPermPtr);

return(numSets);

} /* end of PermuteValues */

/***********************  void UniformIntrons(void)  *************************/
/* */
/* FUNCTION: generates random introns with uniform probability per site */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: randum.h, genes.h, genetype.h, util.h */
/* NOTE: implements 'gOutputRefGenes' option to write sets of generated numbers */
/* to disk (useful for error checking) */
/* */
/*****************************************************************************/

/* identifier string for this model is initialized here: */
static char uniInt[] = "random intron positions, uniform probabilities per site";
static char cUniInt[] = "random intron positions, uniform probabilities per site,  with lower limit";

void UniformIntrons(void)
{
int i, j, k, m, n, temp, temp2, numIntrons, numSets, limit, geneSize;
boolean okToUse;
GeneData *intronPtr=NULL, *nextPtr=NULL;

geneSize = gCurObsInfo[0]->size;
numIntrons = GetNumValues(gGeneDataPtr[0]);

fprintf( stderr, "\nUNIFORM INTRONS (generates random introns with positions that follow");
fprintf( stderr, "\na uniform distribution)\n");
fprintf( stderr, "\nNumber of sets of %d introns to be generated: ", numIntrons);
scanf("%d", &numSets);
ClearLine();
fprintf( stderr, "\nEnter an integer >= 1 for the minimum bp between introns (entering 1 simply");
fprintf( stderr, "\nkeeps the same site from being picked twice in one set): ");
scanf("%d", &limit);
ClearLine();

nextPtr = gGeneDataPtr[1];
i = 0;
while ((i < numSets) && (nextPtr != NULL))
	{
	fprintf( stderr, "\nset %d ", i + 1);
	intronPtr = nextPtr;
	n = 0;
	okToUse = false;
	while ( okToUse == false )
		{
		fprintf( stderr, ". ");
		okToUse = true;
		for (j = 0; j < numIntrons; j++)
			{
			temp2 = 1 + (int) (geneSize - 1) * NRuran2();
			intronPtr->number[j] = temp2;
			}
		for (k = 0; k < numIntrons - 1; k++)
			{
			for (m = k + 1; m < numIntrons; m++)
				{
				temp = abs(intronPtr->number[k] - intronPtr->number[m]);
				if (( temp == 0 ) || ( temp < limit ))
					{
					okToUse = false;
					}
				}
			}
		n++;
		}
	nextPtr = GetGeneStruct();
	intronPtr->nextSet = nextPtr;
	fprintf( stderr, "completed in %d attempts", n);
	i++;
	}
intronPtr->nextSet = NULL;
free(nextPtr);

/* cleaning up */

SortList(gGeneDataPtr[1]);
if (gOutputRefGenes) WriteGeneData(1, true, false);

/* updating current info */

if (limit > 1) strcpy(gCurRefInfo[0]->refModel, cUniInt);
else strcpy(gCurRefInfo[0]->refModel, uniInt);

gCurRefInfo[0]->numSets = numSets;
gCurRefInfo[0]->lowerLimit = limit;
gCurRefInfo[0]->hasPhaseBias = false;
gCurRefInfo[0]->hasPolarity = false;

} /* end of uniform intron generator */

/***********************   void PermuteIIDs(void)  **************************/
/* */
/* FUNCTION: generates sets of intron positions */
/* ARGUMENTS: none */
/* RETURN: none */
/* CALLED BY: GenerateRefGenes submenu */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: genetype.h, infotype.h */
/* NOTE: */
/* */
/*****************************************************************************/

/* identifier string for this model is initialized here: */
static char pIIDInt[] = "permuted order of observed inter-intron distances (in bp)";

void PermuteIIDs(void)
{
int i, numSets=0, numIntrons, geneSize, current, previous, codons;
GeneData *obsIIDPtr, *tempPtr;

/* getting necessary data: */

numIntrons = GetNumValues(gGeneDataPtr[0]);

if (gCurObsInfo[0]->size == 0)
	{
	fprintf( stderr,  "\n\nEnter the total number of codons in the gene: ");
	scanf( "%d", &codons );
	gCurObsInfo[0]->size = 3 * codons;
	ClearLine();
	}
geneSize = gCurObsInfo[0]->size;

/* filling the obsIID struct with inter-intronic distances:  */

obsIIDPtr = GetGeneStruct();

previous = 0;
for (i = 0; i < numIntrons; i++)
	{
	current = gGeneDataPtr[0]->number[i];
	obsIIDPtr->number[i] = current - previous;
	previous = current;
	}
/* filling in the last inter-intron distance, using the gene length - last intron */
obsIIDPtr->number[numIntrons] = geneSize - previous;

/* getting the permuted IID sets */

numSets = PermuteValues(obsIIDPtr, &gGeneDataPtr[1]);

/* now we just need to take the inter-intron distance data, stored in the  */
/* list of reference intron structs, and convert it into intron position form */

tempPtr = gGeneDataPtr[1];

while (tempPtr != NULL)
	{
	previous = 0;
	for (i = 0; i < numIntrons; i++)
		{
		current = tempPtr->number[i];
		tempPtr->number[i] = current + previous;
		previous = current + previous;
		}
	tempPtr->number[numIntrons] = 0; /* i.e., last IID does not correspond to an intron position */
	tempPtr = tempPtr->nextSet;
	}

/* updating current info */

strcpy(gCurRefInfo[0]->refModel, pIIDInt);
gCurRefInfo[0]->numSets = numSets;
gCurRefInfo[0]->lowerLimit = 1;
gCurRefInfo[0]->hasPhaseBias = false;
gCurRefInfo[0]->hasPolarity = false;


}/* end of PermuteIIDs */

/***********************  void LogNormalExons(void)  *************************/
/* */
/* FUNCTION: generates sets of exon sizes with lognormal distribution */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: genetype.h, genes.h, randum.h */
/* NOTE: implements 'gOutputRefGenes' option to write sets of generated numbers */
/* to disk (useful for error checking) */
/* */
/*****************************************************************************/

/* identifier string for this model is initialized here: */
static char lognExn[] = "lognormally distributed random exon sizes (in whole codons)";

void LogNormalExons(void)
{
int	i, j, k, temp, numCodons, sumExonSizes, numExons, numSets;
float	logSize[kMaxNumValues], logNormMu=0, logNormStdDev, meanSqrdDev = 0;
GeneData	*curPtr=NULL, *nextPtr=NULL;

numCodons = GetNumCodons(gGeneDataPtr[2]);
numExons = GetNumValues(gGeneDataPtr[2]);

/* ShowProgress("LognormalExons", 1);*/

fprintf( stderr, "\nLOGNORMAL EXONS (generates random exons with a lognormal size distribution");
fprintf( stderr, "\nbased on the observed gene data)\n");

/* numExons should have already prevented any zeros from entering into */
/* the calculations below, since numExons stops counting when it hits  */
/* a zero in the observed set of exons */

/* first add up the logs of the exon sizes divided by the number of values to get the mean: */

for (i = 0; i < numExons; i++)
	{
	/* ShowLoopProgress("LognormalExons", 2, i);*/
	logSize[i] = (float) log((float) gGeneDataPtr[2]->number[i]);
	logNormMu += (float) (logSize[i]/ numExons);
	}

/* then get the squared deviations from the mean (now using 'logSize' for sqrd dev):*/

for (i = 0; i < numExons; i++)
	{
	/* ShowLoopProgress("LognormalExons", 3, i); */
	logSize[i] = SquareIt(logSize[i] - logNormMu);
	meanSqrdDev += (float) (logSize[i]/numExons);
	}

logNormStdDev = (float) sqrt((double) meanSqrdDev);

fprintf( stderr,  "\nThe lognormal parameters for the observed set are as follows: ");
fprintf( stderr,  "\n\n\tmean: %.4f\n\tstdDev: %.4f", logNormMu, logNormStdDev);

fprintf( stderr,  "\n\nEnter the number of reference sets of %d exons to be generated: ", numExons);
scanf( "%d", &numSets );
ClearLine();

nextPtr = gGeneDataPtr[3];
i = 0;
while ((i < numSets) && (nextPtr != NULL))
	{
	sumExonSizes = 0;
	k = 0;
	curPtr = nextPtr;
	while (sumExonSizes != numCodons)
		{
		sumExonSizes = 0;
		k++;
		for (j = 0; j < numExons; j++)
			{
			temp = 0;
			while (temp <= 0)
				{
				temp = (int) exp((double) NormRand(logNormMu, logNormStdDev));
				}
			curPtr->number[j] = temp;
			sumExonSizes += curPtr->number[j];
			}
		/* fprintf( stderr, "%d ", sumExonSizes);  */
		}
	nextPtr = GetGeneStruct();
	curPtr->nextSet = nextPtr;
	i++;

	fprintf( stderr, "\nset %d completed in %d attempts", i, k);
	}
curPtr->nextSet = NULL;
free(nextPtr);

if (gOutputRefGenes) WriteGeneData(3, true, false);

/* updating current info */

strcpy(gCurRefInfo[1]->refModel, lognExn);
gCurRefInfo[1]->numSets = numSets;
gCurRefInfo[1]->lowerLimit = 1;
gCurRefInfo[1]->hasPhaseBias = false;
gCurRefInfo[1]->hasPolarity = false;

}/** end of LognormalExons module **/

/***********************   void PermuteExons(void)  **************************/
/* */
/* FUNCTION: generates permuted sets of exon sizes */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: infotype.h, genetype.h */
/* NOTE: */
/* */
/*****************************************************************************/

/* identifier string for this model is initialized here: */
static char permExn[] = "permuted order of observed exon sizes (in whole codons)";

void PermuteExons(void)
{
int numSets;

numSets = PermuteValues(gGeneDataPtr[2], &gGeneDataPtr[3]);

/* updating current info */

strcpy(gCurRefInfo[1]->refModel, permExn);
gCurRefInfo[1]->numSets = numSets;
gCurRefInfo[1]->lowerLimit = 1;
gCurRefInfo[1]->hasPhaseBias = false;
gCurRefInfo[1]->hasPolarity = false;

}/** end of PermuteExons  ***/

/********************** void ExponentialExons(void)  *************************/
/*  */
/* FUNCTION: generates sets of random exons with an exponential size distribution */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: genetype.h, util.h, genes.h */
/* NOTE: this function isn't incorporated efficiently at present.  Much of the code  */
/* is the same for this and for UniformIntrons.  In a future revision, both will be  */
/* shells that repeatedly call a third function named IntronPicker, which takes as  */
/* arguments the number of introns to generate and returns a pointer to the  */
/* resulting genedata */
/* NOTE: implements 'gOutputRefGenes' option to write sets of generated numbers to  */
/* disk (useful for error checking) */
/*  */
/*****************************************************************************/

/* identifier string for this model is initialized here: */
static char expnExn[] = "exponentially disributed random exon sizes (in whole codons)";
/* identifier string for this model is initialized here: */
static char cExpnExn[] = "exponentially disributed random exon sizes (in whole codons) with lower limit";

void ExponentialExons(void)
{
int i, j, k, m, n, temp, temp2;
int nextLtBnd, curLtBnd, numCodons, numExons, numSets, limit;
boolean okToUse;
GeneData *exonPtr=NULL, *nextPtr=NULL;

numExons = GetNumValues(gGeneDataPtr[2]);
numCodons = GetNumCodons(gGeneDataPtr[2]);

fprintf( stderr, "\nEXPONENTIAL EXONS (generates random exons with an exponential size");
fprintf( stderr, "\ndistribution)\n");

fprintf( stderr, "\nNumber of sets of %d exons to be generated: ", numExons);
scanf("%d", &numSets);
ClearLine();

fprintf( stderr, "\nEnter an integer >= 1 for the minimum exon size (entering 1 simply");
fprintf( stderr, "\nprevents exons of length 0): ");
scanf("%d", &limit);
ClearLine();

nextPtr = gGeneDataPtr[3];
i = 0;
while ((i < numSets) && (nextPtr != NULL))
	{
	fprintf( stderr, "\nset %d ", i + 1);
	exonPtr = nextPtr;
	n = 0;
	okToUse = false;
	while ( okToUse == false )
		{
		fprintf( stderr, ". ");
		okToUse = true;
		for (j = 0; j < numExons - 1; j++)
			{
			temp2 = 1 + (int) ((numCodons - 1) * NRuran2());
			exonPtr->number[j] = temp2;
			}
		for (k = 0; k < numExons - 2; k++)
			{
			for (m = k + 1; m < numExons - 1; m++)
				{
				temp = abs(exonPtr->number[k] - exonPtr->number[m]);
				if (( temp == 0 ) || ( temp < limit ))
					{
					okToUse = false;
					}
				}
			}
		n++;
		}
	nextPtr = GetGeneStruct();
	exonPtr->nextSet = nextPtr;
	fprintf( stderr, "completed in %d attempts", n);
	i++;
	}
exonPtr->nextSet = NULL;
free(nextPtr);

SortList(gGeneDataPtr[3]);

exonPtr = gGeneDataPtr[3];
while (exonPtr != NULL)
	{
	curLtBnd = 0;
	for (j = 0; j < numExons - 1; j++)
		{
		nextLtBnd = exonPtr->number[j];
		exonPtr->number[j] -= curLtBnd;
		curLtBnd = nextLtBnd;
		}
	exonPtr->number[numExons - 1] = numCodons - curLtBnd;
	exonPtr = exonPtr->nextSet;
	}

if (gOutputRefGenes) WriteGeneData(3, true, false);

strcpy(gCurRefInfo[1]->refModel, expnExn);
if (limit > 1) strcpy(gCurRefInfo[1]->refModel, cExpnExn);
gCurRefInfo[1]->lowerLimit = limit;
gCurRefInfo[1]->numSets = numSets;
gCurRefInfo[1]->hasPhaseBias = false;
gCurRefInfo[1]->hasPolarity = false;

} /* end of ExponentialExons generator */

/* strings for scoring rules are initialized here */

static char extRule[8][80] =
	{
	"closeness score, number of CA-CA distances <= ",
	"extensity score, binary, 0 else 1 if diameter > ",
	"extensity score, number of CA-CA distances > ",
	"extensity score, percent of CA-CA distances > ",
	"extensity score, maximum CA-CA distance in Angstroms, from ",
	"extensity score, average CA-CA distance in Angstroms, from ",
	"extensity score, radius of gyration in Angstroms, from",
	"extensity score, (4/3)*pi*(radgyr^3)/numRes, Angstroms^3, from "
	};

/******************  void RadiusOfGyration(char whichRule)  *********************/
/* */
/* FUNCTION: calculates radius of gyration of exon encoded peptides */
/* ARGUMENTS: scoring rule ( int = 6 or 7) */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: prottype.h, genetype.h, crystal.h */
/* NOTE: rad. gyr: the root mean square deviation from the center of mass */
/* */
/* NOTE: rule number 7 is a crude approximation of the volume per residue; */
/* thus it is expected to be a property of state, like an inverted form of */
/* 'density', completely unlike the other metrics used here, which */
/* are all very dependent on the sizes of the exons.  My limited tests */
/* suggest that this expectation is roughly satisfied: random exon-encoded */
/* peptides obtain similar VPR scores regardless of their size. */
/* */
/*****************************************************************************/

void RadiusOfGyration(char whichRule)
{
int numExons, exonIndex, curLtBnd, curRtBnd;
float distSqrd, meanSqrdDev;
XYZCoord *center, *resPtr, *tempResPtr;
GeneData *exonPtr;

center = GetXYZCoordStruct();
exonPtr = gGeneDataPtr[2];
numExons = GetNumValues(exonPtr);

/* each of these below has to get re-initialized at the bottom of the loop: */
exonIndex = 0;  			/* NOTE: exonIndex 0 is the first exon */
resPtr = gFirstResPtr;

while (exonPtr != NULL)
	{
	fprintf( stderr, " .");
	while ( exonIndex < numExons )
		{
		/* ShowLoopProgress("RadiusOfGyration", 1, exonIndex); */
		center->xCoord = 0.00;
		center->yCoord = 0.00;
		center->zCoord = 0.00;
		tempResPtr = resPtr;
		curLtBnd = resPtr->resNum;
		curRtBnd = curLtBnd + (exonPtr->number[exonIndex] - 1);
		while ( (resPtr != NULL) && ( resPtr->resNum <= curRtBnd ))
			{
			/* ShowProgress("RadiusOfGyration", 2); */
			center->xCoord += (resPtr->xCoord)/(exonPtr->number[exonIndex]);
			center->yCoord += (resPtr->yCoord)/(exonPtr->number[exonIndex]);
			center->zCoord += (resPtr->zCoord)/(exonPtr->number[exonIndex]);
			resPtr = resPtr->nextRes;
			}
		resPtr = tempResPtr;
		meanSqrdDev = 0.0;
		exonPtr->score[exonIndex] = 0.0;
		while ( (resPtr != NULL) && ( resPtr->resNum <= curRtBnd ))
			{
			/* ShowProgress("RadiusOfGyration", 3); */
			distSqrd = SquareIt(InterAtomicDistance(resPtr, center));
			meanSqrdDev += distSqrd/(exonPtr->number[exonIndex]);
			resPtr = resPtr->nextRes;
			}
		if (whichRule == 'r')
			{
			exonPtr->score[exonIndex] = 4.18879 * pow(sqrt(meanSqrdDev), 3.0)/exonPtr->number[exonIndex];
			}
		else
			{
			exonPtr->score[exonIndex] = (float) sqrt(meanSqrdDev);
			}
		exonIndex++;
		}
	if (exonPtr == gGeneDataPtr[2]) exonPtr = gGeneDataPtr[3];
	else exonPtr = exonPtr->nextSet;
	exonIndex = 0;
	resPtr = gFirstResPtr;
	}
}/**  end of RadiusOfGyration ***/

/**********************  void ExtensityScores(void)  *************************/
/* */
/* FUNCTION: calculates extensity of exon-encoded peptides */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: genetype.h, prottype.h, infotype.h, crystal.h */
/* NOTE: rule N (#2) is similar to the rule of Gilbert & Glynias (1994) */
/* NOTE: rule VPR (#7) is still being tested for its applicability */
/* */
/*****************************************************************************/

void ExtensityScores(void)
{
char whichRule, rules[9] = "cbnpmarv";
int numExons, k, curLtBnd, curRtBnd, ruleNum;
long numPairs;
float dist, cutoff=0.0, percentFactor;
XYZCoord *xResPtr, *yResPtr;
GeneData *exonPtr;

exonPtr = gGeneDataPtr[2];
numExons = GetNumValues(exonPtr);

fprintf( stderr, "\nEXTENSITY SCORES (assigns scores to exons based on the extensity of exon-");
fprintf( stderr, "\nencoded peptides, using the atomic coordinates in memory)");

fprintf( stderr,  "\n\nChoose one of the metrics below to calculate extensity using alpha-carbon ");
fprintf( stderr,  "\ncoordinates.  The score assigned to an exon will be the extensity of the ");
fprintf( stderr,  "\nexon-encoded peptide; the score assigned to the entire gene is the average");
fprintf( stderr,  "\nof the exon scores.");
fprintf( stderr,  "\n\n\tRULE\t\tEXTENSITY DEFINED AS");
fprintf( stderr,  "\n\n\tb=binary \t1 if maximum distance > cutoff, else 0");
fprintf( stderr,  "\n\tn=number \tnumber of pairwise distances > cutoff");
fprintf( stderr,  "\n\tp=percent\tpercent of all pairwise distances > cutoff ");
fprintf( stderr,  "\n\tm=maximum\tmaximum of all pairwise distances");
fprintf( stderr,  "\n\ta=average\taverage of all pairwise distances");
fprintf( stderr,  "\n\tr=radius \tradius of gyration");
fprintf( stderr,  "\n\tv=volume \tvol per res, (4/3 * pi * radGyr^3)/numRes");
fprintf( stderr, "\n");

whichRule = GetCommandChar(rules); /* v is not fully tested and has been taken out of this version */
									   /* c does not really fit here, but I'm testing it right now */

ruleNum = (int) (strchr(rules, whichRule) - rules);  /* i.e., convert letter to int 0-7 */
if (ruleNum <= 3)
	{
	fprintf( stderr, "\nEnter a value for the distance cutoff in Angstroms: ");
	scanf("%f", &cutoff);
	ClearLine();
	}

fprintf( stderr, "\nCalculating distances");

if (strchr("rv", whichRule))
	{
	RadiusOfGyration(ruleNum);
	}
else
	{
	/* k & yResPtr also need to be reset at bottom of loop */
	k = 0;  			/* NOTE: exon index 0 is the first exon */
	yResPtr = gFirstResPtr;

	while (exonPtr != NULL)
		{
		fprintf( stderr, " .");
		while ( k < numExons )  /* before running out of exons in current set...*/
			{
			exonPtr->score[k] = 0; /* initialize score to 0 */
			xResPtr = yResPtr;
			curLtBnd = xResPtr->resNum;       /* set left boundary residue*/
			curRtBnd = curLtBnd + (exonPtr->number[k] - 1); /* = end of current exon  */
			numPairs = GetNumPairs(exonPtr->number[k]);
			if (numPairs == 0) numPairs = 1; /* prevents division by 0 in next line */
			percentFactor = (float) 100.00 / numPairs;
			while ( xResPtr->resNum < curRtBnd )
				{
			/* NOTE: an exon of length = 1 will bypass this loop, and will */
			/* receive a score of zero, as initialized above */
				yResPtr = xResPtr;
				while (yResPtr->resNum < curRtBnd)
					{
					yResPtr = yResPtr->nextRes;
					dist = InterAtomicDistance(xResPtr, yResPtr);

					switch(whichRule) /* rule-specific scores assigned here: */
						{
						case 'b':
							{
							if (dist >= cutoff) exonPtr->score[k] = 1;
							break;
							}
						case 'n':
							{
							if (dist >= cutoff) exonPtr->score[k] = exonPtr->score[k] + 1.0;
							break;
							}
						case 'c':
							{
							if (dist <= cutoff) exonPtr->score[k] = exonPtr->score[k] + 1.0;
							break;
							}
						case 'p':
							{
							if (dist >= cutoff) exonPtr->score[k] = exonPtr->score[k] + percentFactor;
							break;
							}
						case 'm':
							{
							if (dist > exonPtr->score[k]) exonPtr->score[k] = dist;
							break;
							}
						case 'a':
							{
							exonPtr->score[k] = exonPtr->score[k] + (float) (dist / numPairs);
							break;
							}
						}
					}
				xResPtr = xResPtr->nextRes;
				}
			yResPtr = yResPtr->nextRes;
			k++;  /* done with this exon, ready to start the next one */
			}
		if (exonPtr == gGeneDataPtr[2]) exonPtr = gGeneDataPtr[3];
		else exonPtr = exonPtr->nextSet;
		k = 0;
		yResPtr = gFirstResPtr;
		}  /* end of while loop */

	} /* end of else statement */

/* updating current info */

if (ruleNum <= 3)
	{
	sprintf(gCurScoreInfo[1]->scoringRule, "%s%.1f Angstroms, from \"%s\"", extRule[ruleNum], cutoff, gCurProtInfo[2]->sourceFile);
	}
sprintf(gCurScoreInfo[1]->scoringRule, "%s\"%s\"", extRule[ruleNum], gCurProtInfo[2]->sourceFile);
gCurScoreInfo[1]->whichType = 'e';
gCurScoreInfo[1]->param1 = 0;
gCurScoreInfo[1]->weightAverage = gWeightAverage;

}/**  end of ExtensityScores ***/


/**********************  void CentralityScores(void)  ************************/
/* */
/* FUNCTION: assigns centrality scores to crystal structure & to introns */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: util.h, genetype.h, prottype.h, crystal.h */
/* */
/*****************************************************************************/

char resRule[6][80] =
	{
	"centrality scores, proportion of inter-CA distances > ",
	"centrality scores, average inter-CA distance, from ",
	"centrality scores, maximum inter-CA distance, from ",
	"centrality scores, distance to center of mass, based on ",
	"centripetal profile scores of Go & Nosaka (1987), based on ",
	"surface accessibility (see details in source file)"
	};

void CentralityScores(void)
{
char response = 0;
int i, domNum, numSegments, rightBnd, whichRule, numResidues;
float dist, cutoff;
XYZCoord *resPtr, *otherResPtr, *domainCenter[10];  /* coordinates of domain center of mass */

for (i = 0; i < 10; i++)
	{
	domainCenter[i] = GetXYZCoordStruct();
	}

resPtr = gFirstResPtr;
numResidues = 0;
while (resPtr != NULL)
	{
	numResidues++;
	resPtr = resPtr->nextRes;
	}

fprintf( stderr, "\nCENTRALITY SCORES (assigns centrality scores to introns using the atomic");
fprintf( stderr, "\ncoordinates in memory)");
fprintf( stderr, "\n\nDoes the protein have more than one domain? (y = yes; n = no) ");
response = GetCommandChar("yn");
if (response == 'y')
	{
	fprintf( stderr, "\nThe contiguous segments composing each domain must be listed separately.");
	fprintf( stderr, "\nEnter the total number of separate contiguous segments (up to 10): ");
	scanf("%d", &numSegments);
	ClearLine();
	fprintf( stderr, "\nEnter the segments in the order in which they occur in the peptide chain.");
	fprintf( stderr, "\nFor each segment, you will be prompted for a number from 0 to 9 to designate");
	fprintf( stderr, "\nthe domain to which the segment belongs.\n");

	resPtr = gFirstResPtr;
	for (i = 0; i < numSegments; i++)
		{
		fprintf( stderr, "\nSegment %d extends from residue ", i + 1);
		fprintf( stderr, "%d to residue: ", resPtr->resNum);
		scanf("%d", &rightBnd);
		ClearLine();
		fprintf( stderr, "and belongs to domain number: ");
		scanf("%d", &domNum);
		ClearLine();
		while (resPtr->resNum <= rightBnd)
			{
			resPtr->domain = domNum;
			domainCenter[domNum]->resNum += 1; /* that is, use 'resNum' to hold domain size */
			resPtr = resPtr->nextRes;
			}
		}
	fprintf( stderr, "\nThe domains and their sizes are as follows:\n");
	for (i = 0; i < 10; i++)
		{
		if (domainCenter[i]->resNum != 0)
			{
			fprintf( stderr, "\n\tDomain %d:\t%d residues", i, domainCenter[i]->resNum);
			}
		}
	fprintf( stderr, "\n\nThis information is not saved, so you might want to print the screen here");
	HoldIt();
	}
else
	{
	domainCenter[0]->resNum = numResidues;
	} /* this allows for mono-domain proteins to be scored using the revised */
	  /* algorithm below, which is based on domSize */

fprintf( stderr,  "\nWhich scoring rule should be used to calculate the score for each intron (in");
fprintf( stderr,  "\neach case, the gene score is the average intron score)?");
fprintf( stderr,  "\n\n\t1. intron score = percent of intra-domain pairwise distances > cutoff ");
fprintf( stderr,  "\n\t2. intron score = average of intra-domain pairwise distances");
fprintf( stderr,  "\n\t3. intron score = maximum of intra-domain pairwise distances");
fprintf( stderr,  "\n\t4. intron score = distance from domain center of mass");
fprintf( stderr, "\n");

whichRule = (int) GetCommandChar("1234") - 48; /* convert to int 1-4 */

if (whichRule == 1)
	{
	fprintf( stderr, "\nEnter the value for the cutoff (Angstroms): ");
	scanf("%f", &cutoff);
	ClearLine();
	gCurScoreInfo[0]->param2 = cutoff;
	}
else cutoff = 0.00;

if (whichRule == 4)
	{
	for (i = 0; i < 10; i++)
		{  /* calculating center of mass for each domain */
		domainCenter[i]->xCoord = 0.00;
		domainCenter[i]->yCoord = 0.00;
		domainCenter[i]->zCoord = 0.00;
		resPtr = gFirstResPtr;
		while (resPtr != NULL)
			{
			if (resPtr->domain == i)
				{
				domainCenter[i]->xCoord += (resPtr->xCoord)/ domainCenter[i]->resNum;
				domainCenter[i]->yCoord += (resPtr->yCoord)/ domainCenter[i]->resNum;
				domainCenter[i]->zCoord += (resPtr->zCoord)/ domainCenter[i]->resNum;
				}
			resPtr = resPtr->nextRes;
			}
		}
	}

fprintf( stderr, "\nCalculating scores for all possible intron positions (a few seconds).");

resPtr = gFirstResPtr;
while (resPtr != NULL)
	{
	gScoreArray[1]->fScore[resPtr->resNum - 1] = 0.00;
	otherResPtr = gFirstResPtr;
	while (otherResPtr != NULL)
		{
		if (otherResPtr->domain == resPtr->domain)
			{
			dist = InterAtomicDistance(resPtr, otherResPtr);
			/* (debug Think C problem) fprintf( stderr, "x "); */
			switch (whichRule)
				{
				case 1:
					{
					if (dist > cutoff)
						{
						gScoreArray[1]->fScore[resPtr->resNum - 1] += (float) 100/(domainCenter[resPtr->domain]->resNum - 1);
						}
					break;
					}
				case 2:
					{
					gScoreArray[1]->fScore[resPtr->resNum - 1] += (float) dist/(domainCenter[resPtr->domain]->resNum - 1);
					break;
					}
				case 3:
					{
					if (dist > gScoreArray[1]->fScore[resPtr->resNum - 1])
						{
						gScoreArray[1]->fScore[resPtr->resNum - 1] = dist;
						}
					break;
					}
				case 4:
					{
					gScoreArray[1]->fScore[resPtr->resNum - 1] = InterAtomicDistance(resPtr, domainCenter[resPtr->domain]);
					break;
					}
				}
			}
		otherResPtr = otherResPtr->nextRes;
		}
	resPtr = resPtr->nextRes;
	}

/* updating current info */

gCurProtInfo[1]->size = numResidues;
GetScoreStats(1);
if (whichRule == 1)
	{
	sprintf(gCurProtInfo[1]->comment, "%s%.1f %s", resRule[whichRule-1], cutoff, gCurProtInfo[2]->sourceFile);
	}
else sprintf(gCurProtInfo[1]->comment, "%s%s", resRule[whichRule-1], gCurProtInfo[2]->sourceFile);


} /* end of CentralityScores */


/********************   void CentripetalProfile(void)    *********************/
/* */
/* FUNCTION: assigns scores to crystal structure & to introns */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: mymath.h, crystal.h, util.h */
/* */
/*****************************************************************************/

void CentripetalProfile(void)
{
XYZCoord *iPtr, *jPtr, *leftBndPtr;
int k=40, numDists, rtBnd, numResidues;
float sum;
boolean done=false;

/* set up */

fprintf( stderr, "\n\t\t\tCentripetal Profile Scores");
fprintf( stderr, "\n\nThis function assigns, to each residue i in the crystal structure, a");
fprintf( stderr, "\nscore equal to the mean squared distance from residue i to residue j,");
fprintf( stderr, "\nwhere j varies in a window from i - k to i + k (Go and Nosaka, 1987)");
fprintf( stderr, "\n\nEnter an integer value for k (40-90): ");
scanf( "%d", &k);
ClearLine();

numResidues = 0;
iPtr = gFirstResPtr;
while (iPtr != NULL)
	{
	numResidues++;
	iPtr = iPtr->nextRes;
	}

leftBndPtr = gFirstResPtr; /* the left window ptr is set at the first residue */
rtBnd = k; /* the right window bound starts at k + 1-- we'll add the 1 inside the loop */
iPtr = gFirstResPtr; /* we start with the first residue */

while (!done)
	{
	numDists = 0;
	sum = 0;

	/* setting left and right boundaries of the window.  As soon as resNum - k is greater
	than 0, then it is time to start moving the left window boundary.  The right bound
	is incremented until it reaches the last residue */

	if ((iPtr->resNum - k) > 0) leftBndPtr = leftBndPtr->nextRes;
	if ((iPtr->resNum + k) < numResidues) rtBnd++;
	else rtBnd = numResidues; /* necessary just in case the protein is < k residues */
	jPtr = leftBndPtr;

	while ((jPtr != NULL) && (jPtr->resNum <= rtBnd))
		{
		if (jPtr != iPtr)
			{
			sum += SquareIt(InterAtomicDistance(iPtr, jPtr));
			numDists++;
			}
		jPtr = jPtr->nextRes;
		}
	if (numDists != 0)
		{
		gScoreArray[1]->fScore[iPtr->resNum - 1] = (float) sum / numDists;
		}
	if (iPtr->resNum == numResidues) done = true;
	else iPtr = iPtr->nextRes;
	}

/* updating info on size, ave, stdDev, scoring rule */

gCurProtInfo[1]->size = numResidues;
GetScoreStats(1);
sprintf(gCurProtInfo[1]->comment, "%s%s", resRule[4], gCurProtInfo[2]->sourceFile);


} /* end of CentripetalProfile */

/*******************  void WriteResArray(boolean toDisk)   ******************/
/* */
/* FUNCTION: writes residue scores to console or file */
/* ARGUMENTS: whether or not to write to disk */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: util.h, fileutil.h, prottype.h, infotype.h */
/* */
/*****************************************************************************/

void WriteResArray(boolean toDisk)
{
char outFile[kMaxNameLength];
int i;
FILE	*filePtr;

if (toDisk)
	{
	filePtr = OpenFile("to contain the residue scores", "w", outFile);
	PrintHeader(filePtr, outFile);
	strcpy(gCurProtInfo[1]->sourceFile, outFile);
	}
else
	{
	filePtr = stderr;
	PrintHeader(filePtr, NULL);
	}
fprintf( filePtr, "\nCOMMENT:  %s", gCurProtInfo[1]->comment);
fprintf( filePtr, "\nPROTFILE  %s", gCurProtInfo[1]->sourceFile);
fprintf( filePtr, "\nNUMRESDU  %d", gCurProtInfo[1]->size);
fprintf( filePtr, "\nAVESCORE  %f", gCurProtInfo[1]->ave);
fprintf( filePtr, "\nSTDDEV    %f", gCurProtInfo[1]->stdDev);
fprintf( filePtr, "\nRESSCORE   ");

if (!toDisk) HoldIt();

for (i = 0; i < gCurProtInfo[1]->size; i++)
	{
	fprintf( filePtr, "%4.2f  ", gScoreArray[1]->fScore[i]);
	}
if (toDisk)
	{
	fclose(filePtr);
	fprintf( stderr, "\nScores for %d residues written to disk in file %s", gCurProtInfo[1]->size, outFile);
	}

} /* end of WriteResArray */

/***********************  void LoadResScores(void)  *************************/
/* */
/* FUNCTION: reads residues scores from disk */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: util.h, fileutil.h, prottype.h, infotype.h */
/* */
/****************************************************************************/

void LoadResScores(void)
{
FILE *filePtr;
char temp=0, lineID[10], discard[100];
int i=0, j=0;

filePtr = OpenFile("containing the scores", "r", gCurProtInfo[1]->sourceFile);

gCurProtInfo[1]->size = 0;
gCurProtInfo[1]->ave = 0.00;
gCurProtInfo[1]->stdDev = 0.00;
strcpy(gCurProtInfo[1]->comment, "none");

while(!feof(filePtr))
	{
	fscanf(filePtr, "%s", lineID);
	if(strstr(lineID, "COMMENT"))
		{
		do /* get rid of leading spaces */
			{
			temp = fgetc(filePtr);
			}
		while (temp == ' ');
		i=0;
		while ((temp != '\n') && (i < 80))
			{
			gCurProtInfo[1]->comment[i] = temp;
			i++;
			temp = fgetc(filePtr);
			}
		gCurProtInfo[1]->comment[i] = '\0';
		}
	if(strstr(lineID, "NUMRESDU"))
		{
		fscanf(filePtr, "%d", &gCurProtInfo[1]->size);
		}
	if (strstr(lineID, "AVESCORE") != NULL)
		{
		fscanf(filePtr, " %f", &gCurProtInfo[1]->ave);
		}
	if (strstr(lineID, "STDDEV") != NULL)
		{
		fscanf(filePtr, " %f", &gCurProtInfo[1]->stdDev);
		}
	if(strstr(lineID, "RESSCORE"))
		{
		if (gCurProtInfo[1]->size == 0)
			{
			fprintf( stderr, "\nfile does not have NUMRESDU line!  Pls. enter number of residues: ");
			scanf("%d", &gCurProtInfo[1]->size);
			ClearLine();
			}
		for (j = 0; j < gCurProtInfo[1]->size; j++)
			{
			fscanf(filePtr, "%f", &gScoreArray[1]->fScore[j]);
			fprintf(stderr, " .");
			}
		}
	fgets(discard, 100, filePtr);
	}
fclose(filePtr);
if (j != gCurProtInfo[1]->size) fprintf(stderr, "\nWarning: %d scores loaded, %d expected", j, gCurProtInfo[1]->size);
else fprintf( stderr, "\n\n%d scores have been loaded from %s.", gCurProtInfo[1]->size, gCurProtInfo[1]->sourceFile);
if (( gCurProtInfo[1]->ave <= 0.00) || (gCurProtInfo[1]->stdDev <= 0.00 ))
	{
	GetScoreStats(1);
	fprintf( stderr, "\nThe average and std dev of scores have been calculated, since the file");
	fprintf( stderr, "\ndid not contain complete information.  Choose s=save to save a complete");
	fprintf( stderr, "\nrecord to a new file.");
	HoldIt();
	}

} /* end of LoadResScores */

/*********************    void AssignResScores(void)   **********************/
/* */
/* FUNCTION: assigns residues scores to any introns in memory */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: prottype.h, genetype.h, genes.h, infotype.h */
/* */
/****************************************************************************/

void AssignResScores(void)
{
int j, numIntrons;
GeneData *intronPtr;
float rScore, lScore;

intronPtr = gGeneDataPtr[0];
numIntrons = GetNumValues(intronPtr);

fprintf( stderr, "\nAssigning scores for observed and any reference intron positions ");

while (intronPtr != NULL)
	{
	fprintf( stderr, ". ");
	j = 0; /* counter for introns in a set */
	while (j < numIntrons)
		{
		if (fmod(intronPtr->number[j], 3)) /* i.e., if its a phase-1 or phase-2 intron */
			{
			intronPtr->score[j] = gScoreArray[1]->fScore[(int) (intronPtr->number[j]/3)];
			}
		else  /* i.e., if its a phase-0 intron, then implement variable gPhaseZero interpretation */
			{
			rScore = gScoreArray[1]->fScore[(int) (intronPtr->number[j]/3)];
			lScore = gScoreArray[1]->fScore[(int) (intronPtr->number[j]/3) - 1 ];
			switch (gPhaseZero)
				{
				case 'r':  /* gets score of rightward (carboxy) residue */
					{
					intronPtr->score[j] = rScore;
					break;
					}
				case 'l':  /* gets score of leftward (amino) residue */
					{
					intronPtr->score[j] = lScore;
					break;
					}
				case 'a':  /* gets average of leftward and rightward residue */
					{
					intronPtr->score[j] = (lScore + rScore)/2;
					break;
					}
				case 'n': /* gets minimum of leftward and rightward scores */
					{
					if (rScore <= lScore) intronPtr->score[j] = rScore;
					else intronPtr->score[j]= lScore;
					break;
					}
				case 'x': /* gets maximum of leftward and rightward scores */
					{
					if (rScore >= lScore) intronPtr->score[j] = rScore;
					else intronPtr->score[j]= lScore;
					break;
					}
				}
			}
		j++;
		}
	if (intronPtr == gGeneDataPtr[0])
		{
		intronPtr = gGeneDataPtr[1];
		}
	else
		{
		intronPtr = intronPtr->nextSet;
		}
	}
strcpy(gCurScoreInfo[0]->scoringRule, gCurProtInfo[1]->comment);
gCurScoreInfo[0]->whichType = 'c';
gCurScoreInfo[0]->param1 = 0;
gCurScoreInfo[0]->weightAverage = false;

} /* end of AssignResScores */


/*******************  void RecordExptInfo(char whichTest)  *******************/
/*  */
/* FUNCTION: fills expt info struct with expt-specific info */
/* ARGUMENTS: which type of test is being done */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: infotype.h, util.h */
/* NOTE: The recorder should be called just before evaluating the data.  It fills */
/* the current expt info struct with information specific to the experiment */
/* that the user has designated for evaluation.  That is, when the user */
/* calls the evaluator to test introns vs. an array (for example) the */
/* recorder fills in the info on obs introns, ref introns, and the array. */
/* This is done by making copies of all necessary current data structs */
/* and giving  them to the current expt record.  Use of the recorder ensures */
/* that the expt info struct for an experiment will contain the data */
/* specific for that experiment that will not be modified by any future */
/* experiments. */
/*  */
/****************************************************************************/

void RecordExptInfo(char whichTest)
{
ObsGeneInfo *dupObsGeneInfo;
RefGeneInfo *dupRefGeneInfo;
ProtInfo *dupProtInfo;
ScoringInfo *dupScoreInfo;

gCurExptInfoPtr->whichType = whichTest;

/* gettings structs to duplicate info on obs, ref, prot & scoring: */

dupObsGeneInfo = (ObsGeneInfo*) Malloc (sizeof(ObsGeneInfo));
dupScoreInfo = (ScoringInfo*) Malloc (sizeof(ScoringInfo));
dupRefGeneInfo = (RefGeneInfo*) Malloc (sizeof(RefGeneInfo));
dupProtInfo = (ProtInfo*) Malloc (sizeof(ProtInfo));

/* picking which structs to duplicate and transferring info */

if ((whichTest == 'c') || (whichTest == 'a')) /* that is, if the test involved introns */
	{
	gOkToTestSet[0] = false;

	*dupObsGeneInfo = *gCurObsInfo[0];
	*dupRefGeneInfo = *gCurRefInfo[0];
	*dupScoreInfo = *gCurScoreInfo[0];
	}

if (whichTest == 'e')  /* test is for extensity */
	{
	gOkToTestSet[1] = false;

	*dupObsGeneInfo = *gCurObsInfo[1];
	*dupRefGeneInfo = *gCurRefInfo[1];
	*dupScoreInfo = *gCurScoreInfo[1];
	}

if (whichTest == 'e')  /* test uses atomic coordinates */
	{
	*dupProtInfo = *gCurProtInfo[2];
	}
if (whichTest == 'c')  /* test uses atomic coordinates */
	{
	*dupProtInfo = *gCurProtInfo[1];
	}
if (whichTest == 'a')  /* test uses array */
	{
	*dupProtInfo = *gCurProtInfo[0];
	}

gCurExptInfoPtr->scorPtr = dupScoreInfo;
gCurExptInfoPtr->obsPtr = dupObsGeneInfo;
gCurExptInfoPtr->refPtr = dupRefGeneInfo;
gCurExptInfoPtr->protPtr = dupProtInfo;

} /* end RecordExptInfo*/

/******************  void EvaluateRef(int whichSet)  *********************/
/* */
/* FUNCTION: reads scores, gets stats & rankings, fills in expt info */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES:  time.h, util.h, prottype.h, infotype.h, genetype.h */
/* NOTE: implements the 'gOutputScores' option to save the scores to disk */
/* NOTE: advances the experiment list */
/* */
/*****************************************************************************/

void EvaluateRef(int whichSet)
{
int numRefSets, lower, higher;
float  sumStdDev, meanStdDev, sumMeans, refMean, sumSqrdDev, sDRefMean;
float percentRank;
GeneData *whichObsPtr = NULL, *whichRefPtr = NULL, *tempPtr = NULL;
FILE *outFilePtr;
char	outFile[] = "refscor.out", geneType[8];
char whichTest = 0;

/* initialize counters */

sumMeans = 0.00;
lower = 0;
higher = 0;
numRefSets = 0;
percentRank = 0.00;

if (whichSet) strcpy(geneType, "exon" );
else strcpy(geneType, "intron");

whichObsPtr = gGeneDataPtr[whichSet * 2];
whichRefPtr = gGeneDataPtr[1 + whichSet * 2];
tempPtr = gGeneDataPtr[1 + whichSet * 2];
whichTest = gCurScoreInfo[whichSet]->whichType;

/* filling in mean and std dev for each set of scores: */

GetStats(whichObsPtr);
gCurExptInfoPtr->obsGeneScore = whichObsPtr->geneScore;
gCurExptInfoPtr->obsSD = whichObsPtr->stdDev;

GetStats(whichRefPtr);

/* ranking scores and calculating mean of all reference sets */

sumMeans = 0;
while (whichRefPtr != NULL)
	{
	if (whichRefPtr->geneScore <= whichObsPtr->geneScore) lower++;
	else higher++;
	sumMeans += whichRefPtr->geneScore;
	whichRefPtr = whichRefPtr->nextSet;
	numRefSets++;
	}

percentRank = (float) (100.00 * lower / numRefSets);
gCurExptInfoPtr->pValue = percentRank / 100.00;
refMean = sumMeans / numRefSets;
gCurExptInfoPtr->refMean = refMean;

/* calculating ave std dev per set, and std deviation of reference mean: */

sumSqrdDev = 0.00;
sumStdDev = 0.00;
whichRefPtr = tempPtr;
while (whichRefPtr != NULL)
	{
	sumStdDev += whichRefPtr->stdDev;
	sumSqrdDev += SquareIt(whichRefPtr->geneScore - refMean);
	whichRefPtr = whichRefPtr->nextSet;
	}
sDRefMean = sqrt(sumSqrdDev / numRefSets);
gCurExptInfoPtr->sDRefMean = sDRefMean;
meanStdDev = sumStdDev / numRefSets;
gCurExptInfoPtr->meanRefSD = meanStdDev;

/* this fills in everything else that has not just been filled in: */

RecordExptInfo(whichTest);

fprintf( stderr, "\n\t\t\t\tEVALUATION");
fprintf( stderr,  "\n\nOBSERVED GENE:\n\tgene score: %.4f", gCurExptInfoPtr->obsGeneScore);
fprintf( stderr,  "\n\tSD of %s scores: %.4f", geneType, whichObsPtr->stdDev);
fprintf( stderr,  "\n\nREFERENCE GENES, MEAN VALUES:\n\tgene score: %.4f (SD = %.4f)", refMean, sDRefMean);
fprintf( stderr,  "\n\tSD of %s scores: %.4f", geneType, meanStdDev);
fprintf( stderr,  "\n\nOf the %d reference genes", numRefSets);
fprintf( stderr,  "\n    %d had a score lower than or equal to that of the observed set", lower);
fprintf( stderr,  "\n    %d had a score higher than that of the observed set", higher);
fprintf( stderr,  "\n\nThus, the percentile ranking of the observed score is ====>  %.2f %%", percentRank );

if (percentRank <= 5)
	{
	fprintf( stderr,  "\n\nThe result is significant!  ");
	}
else
	{
	fprintf( stderr,  "\n\nThe result is not significant.  ");
	}

fprintf( stderr, "\n\nA descriptive comment (up to 80 characters) can be added to the permanent");
fprintf( stderr, "\nrecord of this experiment.");
fprintf( stderr, "\n\nCOMMENT: ");
gets(gCurExptInfoPtr->comment);

WriteExptInfo(gCurExptInfoPtr, stderr);
HoldIt();

if (gOutputRefScores)
	{
	if ((outFilePtr = fopen(outFile, "a")) != NULL)
		{
		PrintHeader(outFilePtr, outFile);
		fprintf( outFilePtr, "\nEXPRMENT  %d: %s", gCurExptInfoPtr->exptNum, gCurExptInfoPtr->comment);
		fprintf( outFilePtr, "\n\nOBSSCOR   %f", whichObsPtr->geneScore);
		whichRefPtr = tempPtr;
		fprintf( outFilePtr, "\n\nREFSCOR   ");
		while (whichRefPtr != NULL)
			{
			fprintf( outFilePtr, "%f  \n", whichRefPtr->geneScore);
			whichRefPtr = whichRefPtr->nextSet;
			}
		fprintf( outFilePtr, "\n------ end of this entry -----\n\n");
		fclose( outFilePtr );
		fprintf( stderr,  "\nThe observed and reference scores for this experiment have been written to disk");
		fprintf( stderr,  "\nin file %s", outFile);
		}
	else fprintf( stderr,  "Note: output file for reference scores could not be opened (check disk space)");
	HoldIt();
	}

}/***** end of EvaluateRef *****/

/******************  void GenerateRefGenes(void)  ****************************/
/* */
/* FUNCTION: submenu calls reference gene generators */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* */
/*****************************************************************************/

void GenerateRefGenes(void)
{
char whichRef=0;

while (whichRef != 'r')
	{
	fprintf( stderr, "\n\nSUBMENU to generate reference data");
	fprintf( stderr, "\n\nReference sets of intron positions");
	fprintf( stderr, "\n\n\tu = UNIFORM intron positions");
	fprintf( stderr, "\n\ti = INTER-INTRON distance permutation to generate intron positions");
/* 	fprintf( stderr, "\n\tb = BIASED phases (incorporate observed biases into reference model)"); */
	fprintf( stderr, "\n\nReference sets of ancestral exons");
	fprintf( stderr, "\n\n\tp = PERMUTED order of inferred ancestral exon sizes");
	fprintf( stderr, "\n\tl = LOGNORMALLY distributed exon sizes with observed mu and stdDev");
	fprintf( stderr, "\n\te = EXPONENTIALLY distributed exon sizes");
	fprintf( stderr, "\n\n\tr = RETURN to main");

	fprintf( stderr, "\n");
	whichRef = GetCommandChar("pleuibr");

	if (strchr("ple", whichRef))
		{
		if ( !GetNumValues(gGeneDataPtr[2]))
			{
			fprintf( stderr, "\nNo observed exons in memory; you must enter this first");
			whichRef = 'r';
			HoldIt();
			}
		else
			{
			strcpy(gCurRefInfo[1]->sourceFile, gCurObsInfo[1]->sourceFile);
			fprintf( stderr, "\nAny previous reference exons now in memory are being cleared");
			EmptyGeneList(gGeneDataPtr[3]);
			}
		}

	if (strchr("ui", whichRef))
		{
		if ( !GetNumValues(gGeneDataPtr[0]))
			{
			fprintf( stderr, "\nNo observed introns in memory; you must enter this first");
			whichRef = 'r';
			HoldIt();
			}
		else
			{
			strcpy(gCurRefInfo[0]->sourceFile, gCurObsInfo[0]->sourceFile);
			fprintf( stderr, "\nAny previous reference introns now in memory are being cleared");
			EmptyGeneList(gGeneDataPtr[1]);
			}
		}

	switch (whichRef)
		{
		case 'l':
			{
			LogNormalExons();
			gOkToTestSet[1] = false;
			break;
			}
		case 'p':
			{
			PermuteExons();
			gOkToTestSet[1] = false;
			break;
			}
		case 'e':
			{
			ExponentialExons();
			gOkToTestSet[1] = false;
			break;
			}
		case 'u':
			{
			UniformIntrons();
			gOkToTestSet[0] = false;
			break;
			}
		case 'i':
			{
			PermuteIIDs();
			gOkToTestSet[0] = false;
			break;
			}
		case 'b':
			{
			fprintf( stderr, "\nSorry. This option is not available in this version of ABaCUS.");
			HoldIt();
			break;
			}
		}
	}
}/** end  GenerateRefGenes()  **/


/**********************  void ElementArrays(void)  ***************************/
/* */
/* FUNCTION: submenu calls functions involving scoring arrays */
/* ARGUMENTS:  none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* */
/*****************************************************************************/

void ElementArrays(void)
{
char command = 0, option=0;

while (command != 'r')
	{
	fprintf( stderr, "\n\nSUBMENU to apply distance-to-boundary scores");
	fprintf( stderr, "\n\n\te = ENTER boundaries of elements & create a scoring array");
	fprintf( stderr, "\n\tl = LOAD a structural element scoring array from a file");
	fprintf( stderr, "\n\tv = VIEW the structural element scoring array currently in memory");
	fprintf( stderr, "\n\tc = CONVERT the scoring array in memory with a new maximum penalty");
	fprintf( stderr, "\n\ts = SAVE the scoring array in memory to a file");
	fprintf( stderr, "\n\ta = ASSIGN scores to introns using the current scoring array");
	fprintf( stderr, "\n\n\tr = RETURN to main menu");

	fprintf( stderr, "\n");
	command = GetCommandChar("elvcsar");

	if ((strchr("vcsa", command)) && (gCurProtInfo[0]->size == 0))
		{
		fprintf( stderr, "\nYou must enter or load an array before trying this!");
		HoldIt();
		command = 0;
		}
	switch (command)
		{
		case 'e':
			{
			option = 0;
			fprintf( stderr, "\nENTER DISCRETE ELEMENTS (creates a scoring array representing the elements)");
			if ((gCurProtInfo[0]->size > 0) && (strstr(gCurProtInfo[0]->sourceFile, "none saved")))
				{
				fprintf( stderr, "\n\nNOTE: there is an array in memory that has not been saved.  This array");
				fprintf( stderr, "\nwill be overwritten should you choose c=continue (s=stop).");
				option = GetCommandChar("cs");
				}
			if (option != 's') EnterElements();
			break;
			}
		case 'l':
			{
			LoadArray();
			break;
			}
		case 'v':
			{
			WriteDistArray(0);
			break;
			}
		case 'c':
			{
			ConvertArray(0);  /* 0=prompt user for limit */
			fprintf( stderr, "\nThe converted array looks like this: \n");
			WriteDistArray(0);
			break;
			}
		case 's':
			{
			WriteDistArray(1);
			break;
			}
		case 'a':
			{
			AssignArrayScore();
			if ((gGeneDataPtr[0] != NULL) && (gGeneDataPtr[1] != NULL))
				{
				gOkToTestSet[0] = true;
				}
			else gOkToTestSet[0] = false;
			break;
			}
		}
	}
}/* end of ElementArrays submenu */


/********************     void GeneDataMenu(void)     ***********************/
/* */
/* FUNCTION: submenu calls functions involving 'observed' gene data */
/* ARGUMENTS:  none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* */
/****************************************************************************/
void GeneDataMenu(void)
{
char option, command=0;
int whichSet;

while (command != 'r')
	{
	fprintf( stderr, "\n\nSUBMENU to manipulate gene data");
	fprintf( stderr, "\n\n\te = ENTER observed intron positions for the first time");
	fprintf( stderr, "\n\tv = VIEW gene data presently in memory");
	fprintf( stderr, "\n\ti = INFER exons from introns (with option to impose sliding rule)");
	fprintf( stderr, "\n\ts = SAVE gene data in memory to a file");
	fprintf( stderr, "\n\tl = LOAD gene data from a file into memory");
	fprintf( stderr, "\n\tn = generate NOISY data by sliding and deletion");
	fprintf( stderr, "\n\n\tr = RETURN to main menu");
	fprintf( stderr, "\n");

	command = GetCommandChar("evsilrn");
	switch (command)
		{
		case 'e':
			{
			option = 'c';
			if ((gGeneDataPtr[0] != NULL) || (gGeneDataPtr[2] != NULL))
				{
				fprintf( stderr,  "\nAny observed data previously in memory will be cleared to make room for");
				fprintf( stderr,  "\nthe new data.  Choose c=continue to proceed with entering new data (s=stop).");
				option = GetCommandChar("cs");
				}
			if (option == 'c') EnterObsGeneData();
			break;
			}
		case 'v':
			{
			fprintf( stderr,  "\nChoose which type of gene data to view, ");
			whichSet = (2 * ChooseExons()) + ChooseRefs();
			if (gGeneDataPtr[whichSet] != NULL) WriteGeneData(whichSet, false, false);
			else ErrNoDataTo("viewing gene data");
			break;
			}
		case 's':
			{
			fprintf( stderr,  "\nThis writes in the append mode.  If the file doesn't exist, it will be ");
			fprintf( stderr,  "\ncreated; if it does exist, new lines will be written at the end.");
			fprintf( stderr,  "\nThis allows you to add reference sets of exons to a previously-existing");
			fprintf( stderr,  "\nfile containing the observed exons.  You can also append hypothesis-");
			fprintf( stderr,  "\ntesting results to the same file.  ");
			fprintf( stderr,  "\n\nChoose which type of data to save, ");
			whichSet = (2 * ChooseExons()) + ChooseRefs();
			if (gGeneDataPtr[whichSet] != NULL) WriteGeneData(whichSet, true, true);
			else ErrNoDataTo("saving gene data");
			break;
			}
		case 'i':
			{
			if (GetNumValues(gGeneDataPtr[0]) > 0) InferExons();
			else ErrNoDataTo("inferring exons");
			break;
			}
		case 'l':
			{
			fprintf( stderr, "\nChoose which type of gene data to load, ");
			LoadGeneData(ChooseExons());
			break;
			}
		case 'n':
			{
			if (gGeneDataPtr[0] == NULL)
				{
				ErrNoDataTo("generate noisy data");
				}
			else
				{
				MakeNoisyGenes();
				fprintf(stderr, "\n\nThe noisy set of %d intron positions is now the 'observed' set in memory.", GetNumValues(gGeneDataPtr[0]));
				WriteGeneData(0, true, false);
				gCurObsInfo[0]->numValues = GetNumValues(gGeneDataPtr[0]);
				gOkToTestSet[0] = false;
				}
			break;
			}
		}
	}
} /* end of GeneDataMenu */


/*********************    void ResScoresMenu(void)     *********************/
/* */
/* FUNCTION: submenu calls functions involving residue scores */
/* ARGUMENTS:  none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* */
/****************************************************************************/
void ResScoresMenu(void)
{
char command=0;

while (command != 'r')
	{
	fprintf( stderr, "\n\nSUBMENU to apply centrality and centripetal profile scores");
	fprintf( stderr, "\n\n\tx = read crystal structure C-alpha coordinates from file");
	fprintf( stderr, "\n\tc = calculate CENTRALITY scores for residues");
	fprintf( stderr, "\n\tp = calculate CENTRIPETAL PROFILE scores for residues");
	fprintf( stderr, "\n\ts = SAVE residue scores in memory to a file");
	fprintf( stderr, "\n\tl = LOAD residue scores from file");
	fprintf( stderr, "\n\tv = VIEW residue scores in memory");
	fprintf( stderr, "\n\ta = ASSIGN residue scores to intron positions in memory");
	fprintf( stderr, "\n\n\tr = RETURN to main menu");
	fprintf( stderr, "\n");

	command = GetCommandChar("xcpslvar");

	if ((gFirstResPtr->nextRes == NULL) && strchr("cp", command))
		{
		fprintf( stderr,  "\nYou must load a crystal structure first!");
		command=0;
		HoldIt();
		}
	switch (command)
		{
		case 'x':
			{
			LoadCrystalStructure();
			break;
			}
		case 'c':
			{
			CentralityScores();
			break;
			}
		case 'p':
			{
			CentripetalProfile();
			break;
			}
		case 's':
			{
			if (gCurProtInfo[1]->size > 0) WriteResArray(1);
			else ErrNoDataTo("saving residue scores");
			break;
			}
		case 'l':
			{
			LoadResScores();
			break;
			}
		case 'a':
			{
			if ((gGeneDataPtr[0] != NULL) && (gCurProtInfo[1]->size > 0))
				{
				AssignResScores();
				if (gGeneDataPtr[1] != NULL)
					{
					gOkToTestSet[0] = true;
					}
				else gOkToTestSet[0] = false;
				}
			else ErrNoDataTo("assigning residue scores");
			break;
			}
		case 'v':
			{
			if (gCurProtInfo[1]->size > 0) WriteResArray(0);
			else ErrNoDataTo("viewing residue scores");
			break;
			}
		}
	}
} /* end of ResScoresMenu */


/************************  void ExtensityMenu(void)  ************************/
/* */
/* FUNCTION: submenu calls functions involving extensity scoring */
/* ARGUMENTS:  none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* */
/****************************************************************************/
void ExtensityMenu(void)
{
char command=0;

while (command != 'r')
	{
	command=0;
	fprintf( stderr, "\n\nSUBMENU to calculate and apply extensity scores ");
	fprintf( stderr, "\n\n\tx = read crystal structure CA coordinates from file");
	fprintf( stderr, "\n\te = calculate and assign EXTENSITY scores to exons");
	fprintf( stderr, "\n\n\tr = RETURN to main menu");
	fprintf( stderr, "\n");

	command = GetCommandChar("xer");

	if ((gFirstResPtr->nextRes == NULL) && (command == 'e'))
		{
		fprintf( stderr,  "\nTo do this, you must load a crystal structure first!");
		command=0;
		HoldIt();
		}
	switch (command)
		{
		case 'x':
			{
			LoadCrystalStructure();
			break;
			}
		case 'e':
			{
			if (gGeneDataPtr[2] != NULL)
				{
				ExtensityScores();
				if (gGeneDataPtr[3] != NULL) gOkToTestSet[1] = true;
				else gOkToTestSet[1] = false;
				}
			break;
			}
		}
	}
} /* end of ExtensityMenu */

/************************* void ChangeSettings(void)  ************************/
/* */
/* FUNCTION: does stuff */
/* ARGUMENTS:  none */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: util.h */
/* NOTE: these are run-time settings */
/* */
/*****************************************************************************/

void ChangeSettings(void)
{
char whichSetting = 0, response = 0;
unsigned int userSeed = 0;

while (whichSetting != 'r')
	{
	whichSetting = 0;
	fprintf( stderr, "\nSUBMENU to alter settings (current values in parentheses)");

	fprintf( stderr, "\n\t0 = initialize random number generator with user-defined seed");
	fprintf( stderr, "\n\t1 = toggle on/off option to write aa sequence of crystal to disk ");
	if (gOutputAASeq) fprintf(stderr, "(on)");
	else fprintf(stderr, "(off)");
	fprintf( stderr, "\n\t2 = change optional output from reference gene generators");
	if (gOutputRefGenes) fprintf( stderr, " (on).");
	else fprintf( stderr, " (off).");
	fprintf( stderr, "\n\t3 = change cryptic output of file with distribution of scores");
	if (gOutputRefScores) fprintf( stderr, " (on).");
	else fprintf( stderr, " (off).");
	fprintf( stderr, "\n\t4 = toggle between weighted and unweighted exon scores ");
	if (gWeightAverage) fprintf( stderr, "(weighted).");
	else fprintf( stderr, "(unweighted).");
	fprintf( stderr, "\n\t5 = treat gene edges as element edges when converting arrays");
	if (gPenalizeEnds) fprintf( stderr, " (off).");
	else fprintf( stderr, " (on) ");
	fprintf( stderr, "\n\t6 = change the sliding limit for intron positions (%d bp)", gSlideLimit);
	fprintf( stderr, "\n\t7 = change how phase-0 introns are interpreted (%c)", gPhaseZero);
	fprintf( 	stderr, 	"\n\n\ta = ACKNOWLEDGEMENTS");
	fprintf( 	stderr, 	"\n\th = HELP");
	fprintf( stderr, "\n\tr = RETURN to main menu");

	fprintf( stderr, "\n");
	whichSetting = GetCommandChar("1234567ahr");

	switch (whichSetting)
		{
		case '0':
			{
			fprintf( stderr, "\n\nThe random number generator is initialized at startup with computer");
			fprintf( stderr, "\nclock time (seconds elapsed since 0:00:00 GMT 1 Jan 1970).  In order");
			fprintf( stderr, "\nallow the possibility of exactly reproducing an experiment, the user");
			fprintf( stderr, "\nmay set the seed manually.  The ONLY reason to do this is if you are");
			fprintf( stderr, "\nan expert user who wishes to reproduce exactly the same set of 'random'");
			fprintf( stderr, "\nnumbers.  Choose y=yes or n=no to enter a manual seed now. ");
			response = GetCommandChar("yn");
			if (response == 'y')
				{
				fprintf( stderr, "\nEnter an unsigned 16-bit integer (i.e., 0 < number < 65,536): ");
				scanf( "%u", &userSeed);
				ClearLine();
				srand(userSeed);
				gSeed = -rand();
				fprintf( stderr, "\nThe number entered was: %u ", userSeed);
				HoldIt();
				}
			break;
			}
		case '1':
			{
			gOutputAASeq = !gOutputAASeq;
			break;
			}
		case '2':
			{
			fprintf( stderr, "\nThis output is normally off.  If changed to 1 = on, you will be prompted");
			fprintf( stderr, "\nto name output files for exponential exons (exon sizes in codons)");
			fprintf( stderr, "\nlognormal exons (logs of exon sizes in codons) and uniform introns");
			fprintf( stderr, "\n(intron position in bp).\n\nEnter 0 (off) or 1 (on): ");
			scanf("%d", &gOutputRefGenes);
			ClearLine();
			fprintf( stderr, "\nCryptic output from random generators is now ");
			if (gOutputRefGenes) fprintf( stderr, "on.");
			else fprintf( stderr, "off.");
			break;
			}
		case '3':
			{
			fprintf( stderr, "\nThe startup setting for this output is 0 = off.  If changed to 1 = on, ");
			fprintf( stderr, "\na file named 'refscor.out' with the means and SD for each reference");
			fprintf( stderr, "\nset will be written each time a hypothesis is evaluated.  Scores from ");
			fprintf( stderr, "\nsuccessive experiments will be appended to the same file, along with ");
			fprintf( stderr, "\ncomments that allow one to identify which scores go with which experiment.");
			fprintf( stderr, "\n\nEnter 0 (off) or 1 (on): ");
			scanf("%d", &gOutputRefScores);
			ClearLine();
			fprintf( stderr, "\nCryptic output of reference scores is now ");
			if (gOutputRefScores) fprintf( stderr, "on.");
			else fprintf( stderr, "off.");
			HoldIt();
			break;
			}
		case '4':
			{
			gWeightAverage = abs(gWeightAverage - 1);
			fprintf( stderr, "\n\nThe gene compactness score will now be the ");
			if (gWeightAverage)
				{
				fprintf( stderr, "weighted average of exon scores");
				fprintf( stderr, "\n(weighted by exon size, that is)");
				}
			else fprintf( stderr, "unweighted average of exon scores.");
			HoldIt();
			break;
			}
		case '5':
			{
			gPenalizeEnds = !gPenalizeEnds;
			fprintf( stderr, "\n\nFor the purposes of converting arrays, the edges of the gene will \n");
			if (gPenalizeEnds)
				{
				fprintf( stderr, "NOT ");
				}
			fprintf( stderr, "be treated as inter-element boundaries.");
			break;
			}
		case '6':
			{
			do
				{
				fprintf( stderr, "\n\nExtant introns as close or closer than the sliding limit are assumed");
				fprintf( stderr, "\nto have arisen from a single ancestral intron whose position is defined");
				fprintf( stderr, "\nas the average of the extant positions.   Allowable values range from");
				fprintf( stderr, "\n3 to 21 bp");
				fprintf( stderr, "\n\nEnter the sliding limit in bp (currently %d bp):  ", gSlideLimit);
				scanf("%d", &gSlideLimit);
				ClearLine();
				}
			while ((gSlideLimit < 3) || (gSlideLimit > 21));
			break;
			}
		case '7':
			{
			fprintf( stderr, "\nThe score assigned an intron between codon 'l' (encoding residue 'L') and");
			fprintf( stderr, "\nand the next codon downstream, codon 'r' (encoding residue R) is" );
			fprintf( stderr, "\ndetermined as follows:" );
			fprintf( stderr, "\n " );
			fprintf( stderr, "\n     when you pick:            the score of the intron is:" );
			fprintf( stderr, "\n" );
			fprintf( stderr, "\n    l=leftward (amino) --------> the score of residue L" );
			fprintf( stderr, "\n    r=rightward (carboxy) -----> the score of residue R" );
			fprintf( stderr, "\n    a=average -----------------> the average score of L and R" );
			fprintf( stderr, "\n    n=minimum -----------------> the lower score of L and R" );
			fprintf( stderr, "\n    x=maximum -----------------> the higher score of L and R" );
			fprintf( stderr, "\n" );
			fprintf( stderr, "\nThe arbitrary default is r=rightward.  To accomodate the equally valid" );
			fprintf( stderr, "\nmethods of Craik, et al., 1982 (Science 299, 180), the m=minimum" );
			fprintf( stderr, "\noption is implemented.  ");
			fprintf( stderr, "\n\nNote that the variable interpretation of phase-0 introns is only relevant");
			fprintf( stderr, "\nto residue-based scoring: centrality, centripetal profile, and surface");
			fprintf( stderr, "\naccessibility.  \n\n" );
			do
				{
				fprintf( stderr, "Enter l=left, r=right, a=ave, n=min, x=max (current setting is '%c'): ", gPhaseZero);
				gPhaseZero = getchar();
				ClearLine();
				}
			while (!strchr("lranx", gPhaseZero));
			break;
			}
		case 'a':
			{
			fprintf( stderr, "\n\n%s", kVersion);
			fprintf( stderr, "\n\nAnalysis of Blake's Conjecture Using Simulations");
			fprintf( stderr, "\nArlin Stoltzfus and David Spencer");
			fprintf( stderr, "\n\nABaCUS was initially developed using the Borland C++ integrated development");
			fprintf( stderr, "\nenvironment in MicroSoft DOS 6.0 on a CAST 486 PC.  ");
			fprintf( stderr, "\n\nThe code was designed entirely by A.S. and D.S., except for the random");
			fprintf( stderr, "\nnumber generator used in the algorithms for generating reference introns and");
			fprintf( stderr, "\nexons.  This pseudo-random number generator is a long-period generator with");
			fprintf( stderr, "\na shuffle, as described on page 282 of W.H. Press, S.A. Teukolsky, W.T. ");
			fprintf( stderr, "\nVetterling and B.P. Flannery, _Numerical Recipes in C_ (Cambridge, Univ. Press,");
			fprintf( stderr, "\n1992, 2nd edition).");
			fprintf( stderr, "\n\nThe radius of gyration metric for extensity was suggested by Michael Zuker.");
			fprintf( stderr, "\n\nThis work was funded indirectly by the Medical Research Council of Canada, the");
			fprintf( stderr, "\nNational Research Council of Canada, and the Canadian Institute for Advanced");
			fprintf( stderr, "\nResearch (CIAR) Program in Evolutionary Biology");
			HoldIt();
			}
		case 'h':
			{
			fprintf( stderr, "\t\t\t   HOW TO GET HELP");
			fprintf( stderr, "\n\nThere is no on-line help for ABaCUS.  Please read the help documents");
			fprintf( stderr, "\nprovided with the program, and check out 'Methods for Evaluating Exon-");
			fprintf( stderr, "\nProtein Correspondences' (Stoltzfus, Spencer and Doolittle, 1995, CABIOS");
			fprintf( stderr, "\n11(5), 509-515).  If these do not answer your questions, contact:");
			fprintf( stderr, "\n\n\t\t\t  Dr. Arlin Stoltzfus");
			fprintf( stderr, "\n\t\t      Department of Biochemistry");
			fprintf( stderr, "\n\t\t\t  Dalhousie University");
			fprintf( stderr, "\n\t\t   Halifax, Nova Scotia B0J 3J0 CANADA");
			fprintf( stderr, "\n\t\t       internet: arlin@is.dal.ca");
			fprintf( stderr, "\n\t\t\t  fax: (902) 494-1355");
			fprintf( stderr, "\n\nI will also be delighted to hear about any results obtained using ABaCUS.");
			fprintf( stderr, "\nRemember that both positive and negative results are important: they provide ");
			fprintf( stderr, "\nvalid information about the way the world is (or isn't, as the case may be)!");
			fprintf( stderr, "\n\nIf you experience a problem with ABaCUS, please indicate that this is");
			fprintf( stderr, "\n%s.  If ABaCUS refuses to read a PDB crystal", kVersion);
			fprintf( stderr, "\nstructure file, please tell me the name of the file and I'll suggest how to");
			fprintf( stderr, "\nedit the file, or send you a readable copy of it.");
			HoldIt();
			}
		}
	}
}
/* end of changesettings */

/**********************  char  GetMainCommand(void)  *************************/
/* */
/* FUNCTION: prints menu, prompts user, returns command to main */
/* ARGUMENTS:  none */
/* RETURN: single-letter command */
/* PROTOTYPE IN: abacus.h */
/* OTHER DEPENDENCIES: util.h */
/* */
/*****************************************************************************/

char	GetMainCommand(void)
{
char	command = 0;

fprintf( stderr, "\n\t\t\t\tMAIN MENU");
fprintf( stderr, "\n\nSubmenus");

fprintf( stderr, "\n\n\tg = Gene Data (to enter, load, save, view gene data)");
fprintf( stderr, "\n\tr = Reference Models (to generate reference data)");
fprintf( stderr, "\n\td = Distance-to-Boundary Scores (to score intron positions)");
fprintf( stderr, "\n\tc = Centrality & Centripetal Profile Scores (to score intron positions)");
fprintf( stderr, "\n\te = Extensity Scores  (to score exons)");

fprintf( stderr, "\n\nOther options");
fprintf( stderr, "\n\n\ts = Settings, Help, Info (miscellaneous run-time settings)");
fprintf( stderr, "\n\ti = Information (to view current data or experiment list)");
fprintf( stderr, "\n\tt = Test (to evaluate the reference hypothesis)");
fprintf( stderr, "\n\tq = Quit");
if (gCurExptInfoPtr->exptNum > 1)
	{
	fprintf( stderr, " (and save experiment summaries to disk)");
	}
else fprintf( stderr, " and get back to work");

fprintf( stderr, "\n");
command = GetCommandChar("grdcesitq");
return( command );
}

/*********************   void DoMainCommand(char)  ***************************/
/* */
/* FUNCTION: executes main command */
/* ARGUMENTS:  command supplied by GetMainCommand */
/* RETURN: none */
/* PROTOTYPE IN: abacus.h */
/* */
/*****************************************************************************/

void DoMainCommand(char command)
{
char whichType=0;
int whichSet;

switch(command)
	{
	case 'g':
		{
		GeneDataMenu();
		break;
		}
	case 'r':
		{
		GenerateRefGenes();
		break;
		}
	case 'd':
		{
		ElementArrays();
		break;
		}
	case 'c':
		{
		ResScoresMenu();
		break;
		}
	case 'e':
		{
		ExtensityMenu();
		break;
		}
	case 's':
		{
		ChangeSettings();
		break;
		}
	case 't':
		{
		if ((!gOkToTestSet[0]) && (!gOkToTestSet[1]))
			{
			fprintf( stderr, "\n\nNo test is possible at present, due to one of the following conditions:");
			fprintf( stderr, "\n\n -no complete (both observed & reference) set of gene data exists;");
			fprintf( stderr, "\n -a complete set of gene data exists, but has not been scored;");
			fprintf( stderr, "\n -a complete set of scored gene data exists, but has already been evaluated.");
			HoldIt();
			}
		else
			{
			if (gOkToTestSet[0])
				{
				if (gOkToTestSet[1])
					{
					fprintf( stderr, "Choose which type of data to evaluate, ");
					whichSet = ChooseExons();
					}
				else whichSet = 0;
				}
			else whichSet = 1;
			EvaluateRef(whichSet);
			gCurExptInfoPtr = NextExptInfoPtr();
			gOkToTestSet[whichSet] = false;
			}
		break;
		}
	 case 'i':
		{
		if (gCurExptInfoPtr->exptNum <= 1) whichType = 'c';
		else
			{
			fprintf( stderr, "\nInfo for p=past experiments or c=current data in memory?");
			whichType = GetCommandChar("pc");
			}
		if (whichType == 'c') ShowCurrentInfo();
		if (whichType == 'p') WriteExptList(false);
		break;
		}
	}
} /* end of v. 0.55a DoMainCommand */

void DoLastRites(void)
{
	char ch ;
if (gCurExptInfoPtr->exptNum > 1)
	{
	fprintf(stderr, "\n!\n!\n!ATTN: data in memory cannot be saved to disk; last chance is to record the");
	fprintf(stderr, "\ninfo that follows by hand or with screen dumps.  Press <enter> to proceed.");
	ch = gechar();
	while ( ch != '\n' );
	WriteExptList(false);
	}
fprintf(stderr, "\nexiting abacus due to error.  bye.");
exit(0);
}

/***************************    void Main(void)    ***************************/
/* */
/* FUNCTION: starts up, gets command, executes, writes summary file */
/* ARGUMENTS:  none */
/* RETURN: none */
/* */
/*****************************************************************************/

void main()
{
char  command=0;

InitGlobals();

ShowStartupScreen();

while (command != 'q')
	{
	DoMainCommand(command = GetMainCommand());
	}

if (gCurExptInfoPtr->exptNum > 1)
	{
	WriteExptList(true);
	}
} /* end of main */

/* end of main code block */
