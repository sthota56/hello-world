/* aaseq.c v. 0.57 */

#include <string.h>
#include "48726631_1454109756.h"

/* residues in rough order of abundance; this represents an optimization for speed, */
/* given that the algorithm below moves through the list from left to right */

static triplet codeStrings[20] = {"ALA","GLY","LYS","LEU","VAL","THR","SER","ASP","GLU","PHE","ASN","PRO","ILE","HIS","ARG","GLN","TYR","CYS","MET","TRP"};
static char codeChars[21] = "AGKLVTSDEFNPIHRQYCMW";

/*************** char SingleLetterCode(char source[4])  **********************/
/* */
/* FUNCTION: gets single-letter code corresponding to triplet code string */
/* ARGUMENTS:  string with triplet to be interpreted */
/* RETURN: single-letter code */
/* INCLUDE:  aaseq.h */
/* NOTE: returns '?' for no match */
/* */
/*****************************************************************************/

char SingleLetterCode(char source[4])
{
int i=0;
char slc = '?'; /* initialize to '?': if no match, then function returns '?' */

while (i < 20)
	{
	if (strstr(source, codeStrings[i]))
		{
		slc = codeChars[i];
		i = 20;
		}
	i++;
	}

return(slc);
}

