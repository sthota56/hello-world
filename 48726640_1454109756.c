/* file crystal.c v. 0.58 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "48726637_1454109756.h"

#include "48726654_1454109756.h"
#include "48726649_1454109756.h"

#include "48726658_1454109756.h"
#include "48726645_1454109756.h"
#include "48726631_1454109756.h"
#include "48726651_1454109756.h"

#include "48726641_1454109756.h"

extern XYZCoord *gFirstResPtr;
extern ProtInfo *gCurProtInfo[3];
extern boolean gOutputAASeq;

void DumpOldCoordinates(XYZCoord*);  /* private function */

/********************** XYZCoord* GetXYZCoordStruct(void)  *******************/
/* */
/* FUNCTION: allocates and initializes new struct */
/* ARGUMENTS:  none */
/* RETURN: pointer to the new struct */
/* PROTOTYPE IN: crystal.h */
/* OTHER DEPENDENCIES: prottyp.h, util.h */
/* NOTE: syntax is "myNewDataPtr = GetXYZCoordStruct();" */
/* */
/*****************************************************************************/

XYZCoord* GetXYZCoordStruct(void)
{
XYZCoord *newPtr;

newPtr = (XYZCoord*) Malloc( sizeof( XYZCoord ));

newPtr->xCoord = 0.0;
newPtr->yCoord = 0.0;
newPtr->zCoord = 0.0;
newPtr->resNum = 0;
newPtr->nextRes = NULL;
newPtr->domain = 0;

return ( newPtr );
} /* end of GetXYZCoordStruct */

/****************   void DumpOldCoordinates(XYZCoord *freePtr)  **************/
/* */
/* FUNCTION: initializes first member; frees succeeding list of structs */
/* ARGUMENTS:  pointer to struct or first member of list of structs */
/* RETURN: none */
/* PROTOTYPE IN: crystal.h */
/* OTHER DEPENDENCIES: prottype.h */
/* CALLS: */
/* */
/*****************************************************************************/

void DumpOldCoordinates(XYZCoord *firstPtr)
{
XYZCoord *temp;

while (firstPtr->nextRes != NULL)
	{
	temp = firstPtr->nextRes->nextRes;
	free(firstPtr->nextRes);
	firstPtr->nextRes = temp;
	}
/* now clean out the first struct of its old values: */
temp = GetXYZCoordStruct();
*firstPtr = *temp;
}

/*********** float InterAtomicDistance(XYZCoord *ptr1, XYZCoord *ptr2)  *****/
/* */
/* FUNCTION: calculates distance between two atoms */
/* ARGUMENTS:  ptrs to coordinates */
/* RETURN: result of calculation */
/* PROTOTYPE IN: crystal.h */
/* OTHER DEPENDENCIES:  mymath.h, prottype.h */
/* NOTE: this function designed to obviate type casting problems with pow */
/* */
/*****************************************************************************/

float InterAtomicDistance(XYZCoord *ptr1, XYZCoord *ptr2)
{
float xSqrd, ySqrd, zSqrd, dist;
double sum;

xSqrd = SquareIt((ptr1->xCoord) - (ptr2->xCoord));
ySqrd = SquareIt((ptr1->yCoord) - (ptr2->yCoord));
zSqrd = SquareIt((ptr1->zCoord) - (ptr2->zCoord));

sum = xSqrd + ySqrd + zSqrd;
dist = (float) sqrt(sum);

return (dist);
} /* end InterAtomicDistance */

/********************  void LoadCrystalStructure(void)  **********************/
/* */
/* FUNCTION: reads PDB files (some have to be edited) */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: crystal.h */
/* OTHER DEPENDENCIES: util.h, fileutil.h, aaseq.h */
/* NOTE: reads user-defined file; writes automatic 'calpha.xyz' file; writes */
/* optional 'crystal.seq' file with single-letter code aa sequence */
/* */
/*****************************************************************************/

void LoadCrystalStructure(void)
{
char field1[10],field3[5],field4[4], discard[100];
char inFileName[kMaxNameLength], outFileName[kMaxNameLength] = "calpha.xyz";
int i, j, k=1, l=0, numOddballs = 0;
int field2, field5;
FILE	*CAlphaFPtr=NULL, *inPDBFilePtr=NULL, *AASeqFilePtr=NULL;
XYZCoord	*infoPtr=NULL, *lastPtr=NULL;
boolean consecutive = true;

char aaSeq[kMaxArraySize/3];
 
char aaSeq[1400];

fprintf( stderr, "\nLOAD ATOMIC COORDINATES\n");
fprintf( stderr, "\nThis function reads some PDB files. If the PDB file contains coordinates for");
fprintf( stderr, "\nmultiple subunits, you must first use a text editor to throw away data for");
fprintf( stderr, "\nall but one subunit, then remove the subunit designator field (e.g., search");
fprintf( stderr, "\nand replace ' A ' with '   ', if 'A' is the subunit designator).  The C-alpha");
fprintf( stderr, "\ncoordinates extracted from the PDB file are stored memory for use by other ");
fprintf( stderr, "\nfunctions, and are written to a file named 'calpha.xyz'");
fprintf( stderr, "\n");

inPDBFilePtr = OpenFile("containing the CA coordinates", "r", inFileName);

/* the following is to avoid conflict if the input name is calpha.xyz */
if (strstr(inFileName, outFileName)) strcpy(outFileName, "calpha2.xyz");

/* opening output file */
CAlphaFPtr = fopen(outFileName, "w"); /* This stores the stripped data */
PrintHeader(CAlphaFPtr, outFileName);
fprintf(CAlphaFPtr, "NOTE    contains CA lines extracted from source file %s\n", inFileName);

/* memory allocation for data; pointers to linked list of XYZ coordinates */

DumpOldCoordinates(gFirstResPtr); /* this deallocates the previous data set*/
infoPtr = gFirstResPtr;

j = 1;
while(!feof(inPDBFilePtr))
	{
	fscanf(inPDBFilePtr, "%s", field1);
	if(strstr(field1, "ATOM"))
		{
		fscanf(inPDBFilePtr, "%d %s", &field2, field3);
		if(strstr(field3, "CA"))
			{
			/* read remainder of the current line: */
			fscanf(inPDBFilePtr, "%s %d %f %f %f", field4, &field5, &infoPtr->xCoord, 
				&infoPtr->yCoord, &infoPtr->zCoord);
		
			/* print to output file: */
			fprintf(CAlphaFPtr, "ATOM  %d  %s  %s  %d  %f  %f  %f  \n", field2, field3, 
				field4, field5, infoPtr->xCoord, infoPtr->yCoord, infoPtr->zCoord);

			/* echo to console: */
			fprintf(stderr, "ATOM  %d  %s  %s  %d  %f  %f  %f  \n", field2, field3,
				field4, field5, infoPtr->xCoord, infoPtr->yCoord, infoPtr->zCoord);
			/* print to sequence file if necessary: */
			if (gOutputAASeq)
				{
				aaSeq[j-1] = SingleLetterCode(field4);
				if (aaSeq[j-1] == '?') numOddballs++;
				}

			/* assign resNum, note any discrepancy in numbering, and set pointers: */
			infoPtr->resNum = j++;
			if (infoPtr->resNum != field5) consecutive = false;
			infoPtr->nextRes = GetXYZCoordStruct();
			lastPtr = infoPtr;
			infoPtr = infoPtr->nextRes;
			}
		}
	fgets(discard, 100, inPDBFilePtr);
	}
fprintf(CAlphaFPtr, "END\n");
fclose(CAlphaFPtr);
fclose(inPDBFilePtr);
free(lastPtr->nextRes);
lastPtr->nextRes = NULL;

i = 0;

if (gOutputAASeq)
	{
	AASeqFilePtr = fopen("crystal.seq", "w");
	PrintHeader(AASeqFilePtr, "crystal.seq");
	fprintf(AASeqFilePtr, "NOTE  amino acid sequence from source file %s\n\n", inFileName);
	fprintf(stderr, "\nAmino acid sequence in single-letter code: \n\n");
	while (i < (j - 1))
		{
		fprintf(AASeqFilePtr, "%c", aaSeq[i]);
		fprintf(stderr, "%c", aaSeq[i]);
		if (k > 9)   /* space every 10 amino acids */
			{
			fprintf(AASeqFilePtr, " ");
			fprintf(stderr, " ");
			k = 0;
			l++;
			}
		if (l > 5) /* line break every 60 amino acids */
			{
			fprintf(AASeqFilePtr, "\n");
			fprintf(stderr, "\n");
			l = 0;
			}
		k++;
		i++;
		}
	fprintf(AASeqFilePtr, "\n\n-- end of sequence --");
	fclose(AASeqFilePtr);
	}

/* updating current info */

strcpy(gCurProtInfo[2]->sourceFile, inFileName);
gCurProtInfo[2]->size = j-1;
gCurProtInfo[2]->param1 = 0;
gCurProtInfo[2]->ave = 0.00;

fprintf( stderr, "\n\nThe file '%s' has been read, and %d C-alpha coordinates have been", inFileName, j - 1);
fprintf( stderr, "\nloaded.  ");
if (!consecutive)
	{
	fprintf( stderr, "The residues in the input file were numbered non-consecutively ");
	fprintf( stderr, "\n(for internal use only, ABaCUS assigns consecutive numbers to the");
	fprintf( stderr, "\nresidues).  ");
	}
if (gOutputAASeq)
	{
	fprintf(stderr, "The sequence has been recorded in the file 'crystal.seq'.");
	if (numOddballs > 0)
		{
		fprintf(stderr, "\n\nNOTE: %d residues could not be read.  These appear as 'x' in the sequence.", numOddballs);
		}
	}
HoldIt();

} /*end LoadCrystalStructure */
