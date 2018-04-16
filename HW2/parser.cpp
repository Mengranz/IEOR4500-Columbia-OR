#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int parser(char *sourcefilename, double *solution, int numsec, int *pstatus)
{
 	int retcode = 0; 
	char mybuffer[100];
	int j, k, numnonz; 
	double optimalvalue, xvalue;
	FILE *results = NULL; 

	results = fopen(sourcefilename, "r");
    if (!results){
	  printf("cannot open %s\n", sourcefilename); retcode = 300; 
	  goto BACK;
    }
  /* read until finding Optimal */

	*pstatus = 100;

    for (;;){
	  fscanf(results, "%s", mybuffer);
	  /* compare mybuffer to 'Optimal'*/
	  if (strcmp(mybuffer, "Optimal") == 0){
		  /* now read three more*/
		  fscanf(results, "%s", mybuffer);
		  fscanf(results, "%s", mybuffer);
		  fscanf(results, "%s", mybuffer);
		  optimalvalue = atof(mybuffer);
		  printf(" value = %g\n", optimalvalue);
		  if (optimalvalue < -.0001){
			  printf("type A arbitrage exists!\n");
			  *pstatus = 0;
			  /* read again to get the number of nonzeros*/
			  fscanf(results, "%s", mybuffer);
			  numnonz = atoi(mybuffer);
			  fscanf(results, "%s", mybuffer); fscanf(results, "%s", mybuffer); fscanf(results, "%s", mybuffer);
			  fscanf(results, "%s", mybuffer);
			  for (k = 0; k < numnonz; k++){
				  fscanf(results, "%s", mybuffer);
				  j = atoi(mybuffer + 1);
				  fscanf(results, "%s", mybuffer); fscanf(results, "%s", mybuffer);
				  xvalue = atof(mybuffer);
				  solution[j] = xvalue;
				  printf("%d -> %g\n", j, xvalue);
			  }
		  }
		  else{
			  *pstatus = 1;
			  printf("no type A\n"); break;
		  }
	  }
	  else if (strcmp(mybuffer, "bye.") == 0){
		  break;
	  }
	}

	fclose(results);
   
BACK:
  return retcode;
}
 