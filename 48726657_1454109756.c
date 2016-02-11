/* util.c v. 0.58 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>

#include "48726658_1454109756.h"

extern void DoLastRites(void); /* define this elsewhere as the exit function */
						 /* for the application */

/**************************  void *Malloc(size)    ***************************/
/* */
/* FUNCTION: calls malloc to allocate memory dynamically */
/* ARGUMENTS: size to be allocated */
/* RETURN: ptr to allocation */
/* PROTOTYPE IN: util.h */
/* OTHER DEPENDENCIES: */
/* NOTE: this function circumvents some non-portable type-casting problems; */
/* NOTE: the function DoLastRites is for a parting message or to save data, */
/* and must contain an exit. */
/* */
/*****************************************************************************/

void *Malloc(size_t size)
{
void *newPtr;

newPtr = (void *) malloc(size);
if (!newPtr)
	{
	fprintf(stderr, "\nOut of memory!");
	DoLastRites();
	}

return((void *) newPtr);

} /* end of Malloc */

/**************************  void ClearLine(void)  ***************************/
/* */
/* FUNCTION: clears the hanging carriage return after a user-supplied command */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: util.h */
/* OTHER DEPENDENCIES: none */
/* */
/*****************************************************************************/

void	ClearLine(void)
{
while ( getchar() != '\n' );
}

/***************************  void HoldIt(void)  *****************************/
/* */
/* FUNCTION: prompts for a <CR> in order to allow the user a chance to view what */
/* is being displayed before it scrolls out of view; or just to make a pause */
/* before going to the next step */
/* ARGUMENTS: none */
/* RETURN: none */
/* PROTOTYPE IN: util.h */
/* OTHER DEPENDENCIES: */
/* */
/*****************************************************************************/

void 	HoldIt()
{
fprintf( stderr, "\n\n\t\t\t(hit 'enter' to continue)");
ClearLine();
}

/*********************  ErrNoDataTo(char operation[30])  ************************/
/* */
/* FUNCTION: gives message to user when operation cannot be completed */
/* ARGUMENTS: description of operation */
/* RETURN: none */
/* PROTOTYPE IN: util.h */
/* OTHER DEPENDENCIES: */
/* */
/*****************************************************************************/

void ErrNoDataTo(char operation[30])
{
fprintf(stderr, "\nThe data necessary for %s do not yet exist.", operation);
HoldIt();
}

/*****************  char GetCommandChar(char options[20])  *******************/
/* */
/* FUNCTION: gets a single letter command from the set of options */
/* ARGUMENTS: string containing all letter options */
/* RETURN: single char */
/* PROTOTYPE IN: util.h */
/* OTHER DEPENDENCIES: */
/* NOTE: up to 20 commands possible (this can be changed); DON'T use 'z' as */
/* a command, since this letter was chosen arbitrarily to initialize */
/* 'response' with a non-matching value. */
/* */
/*****************************************************************************/

char GetCommandChar(char options[20])
{
char response='z';

fprintf( stderr, "\n");   /* puts an extra blank line before the first query */
while (!strchr(options, response))
	{
	fprintf( stderr, "Enter command: ");
	response = getchar();
	ClearLine();
	}
return (response);
} /* end of GetCommandChar */


/*****   debugging toys ShowProgress, HoldProgress and ShowLoopProgress  *****/
/* */
/* use theses to place markers in functions to aid debugging without a real */
/* debugger.  e.g., ShowProgress("DoStuff", 3) would print a message */
/* informing the user that the program is "in function DoStuff at */
/* progress marker 3" whenever the code is reached.  If there is no scrollable */
/* console or the program is blowing up at the bug, use HoldProgress.  To */
/* track within iterated loops, use ShowLoopProgress */
/* */
/*****************************************************************************/

void ShowProgress(char function[20], int marker)
{
fprintf(stderr, "\nin function %s, at progress marker %d;", function, marker);
}

void HoldProgress(char function[20], int marker)
{
fprintf(stderr, "\nin function %s, holding at progress marker %d;", function, marker);
HoldIt();
}

void ShowLoopProgress(char function[20], int marker, int index)
{
fprintf(stderr, "\nin function %s, at loop marker %d, at index = %d;", function, marker, index);
/* HoldIt(); */
}

