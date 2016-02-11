/* fileutil.c v. 0.58 */

#include <stdio.h>
#include <time.h>

#include "48726637_1454109756.h"
#include "48726658_1454109756.h"
#include "48726645_1454109756.h"

/******************** void PrintHeader(FILE *, char[])  ******************/
/* */
/* FUNCTION: writes header for  output files and display */
/* ARGUMENTS:  pointer to destination, filename */
/* RETURN: nada */
/* PROTOTYPE IN: fileutil.h */
/* OTHER DEPENDENCIES: time.h */
/* NOTE: to write to screen, call with (stderr, NULL) */
/* */
/*************************************************************************/

void PrintHeader(FILE *destPtr, char fileName[kMaxNameLength])
{
struct tm *date;
time_t	now;
char timeString[80];

time ( &now );
date = localtime( &now );
strftime(timeString, 80, "%c", date);

if (fileName == NULL)
	{
	fprintf( destPtr, "\nData currently in memory ( %s )", timeString);
	}
else
	{
	fprintf( destPtr, "\nThis file originally named %s", fileName);
	fprintf( destPtr, "\nCreated by %s on %s\n", kVersion, timeString);
	}

} /* end PrintHeader */


/****** FILE *OpenFile(char whatFor[80], char mode, char &name[kMaxNameLength]) *******/
/* */
/* FUNCTION: opens a file for requested purpose, in requested mode */
/* ARGUMENTS: string description, mode, address of name */
/* RETURN: pointer to opened file */
/* PROTOTYPE IN: fileutil.h */
/* OTHER DEPENDENCIES: util.h */
/* NOTE: */
/* */
/**************************************************************************************/

FILE *OpenFile(char whatFor[80], char mode[5], char fileName[80])
{
FILE *fp;

fprintf( stderr, "\nEnter the name for the file %s: ", whatFor);
do
	{
	scanf("%s", fileName);
	ClearLine();
	fp = fopen ( fileName, mode );
	if (fp == NULL)
		{
		fprintf( stderr, "\nDidn't open. Please try again: ");
		}
	}
	while (fp == NULL);

return(fp);
}


