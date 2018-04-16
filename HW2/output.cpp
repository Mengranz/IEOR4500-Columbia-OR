#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int output(char *solutionfilename, double *solution, int numsec, int status)
{
 	int retcode = 0;  
	int j;  
	FILE *out = NULL; 

	out = fopen(solutionfilename, "w");
    if (!out){
	  printf("cannot open final output file %s for writing\n", solutionfilename); retcode = 400; 
	  goto BACK;
    }
  /* read until finding Optimal */

	fprintf(out, "status %d\n", status);
	if (status != 0){
		fclose(out); goto BACK;
	}

	for (j = 0; j <= numsec; j++){
	  fprintf(out, "%g\n", solution[j]); 
	}

	fprintf(out, "\nEND\n");

	fclose(out);
   
BACK:
  return retcode;
}
 