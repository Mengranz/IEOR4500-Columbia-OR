#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <sys/types.h>
#include <sys/stat.h>

#include <windows.h>
#include <float.h>
#include <cmath>

using namespace std;
double power(int numsec, double *covariance);

int __stdcall eigen(char *datafile, int numsec, int numday, int rolldays) {/*number of assets;number of days*/

	FILE *in = NULL, *outputFile = NULL,*outputFile2 = NULL, *outputFile3 = NULL;
	/*int k, n, j;*/
	int retcode = 0;
	char mybuffer[100];
	double *p, *ret, *sum, *mean, *covarsum, *covar;
	double *eigenvalue;

	in = fopen(datafile, "r"); /*price datafile*/


	if (in == NULL) {
		printf("could not open %s for reading\n", datafile);
		retcode = 200; goto BACK;
	}

	/*p is the price data read from the data file*/
	p = (double *)calloc(numsec*numday, sizeof(double));
	if (p == NULL) {
		printf("no memory\n"); retcode = 400; goto BACK;
	}
	for (int k = 0; k < numsec; k++) {
		for (int n = 0; n < numday; n++) {
			fscanf(in, "%s", mybuffer);
			p[k*numday + n] = atof(mybuffer);
		}
	}
	printf("Finish reading price data!\n");
	fclose(in);

	ret = (double *)calloc(numsec*(numday-1), sizeof(double));
	sum = (double *)calloc(numsec, sizeof(double));
	mean = (double *)calloc(numsec, sizeof(double));
	covarsum = (double *)calloc(numsec*numsec, sizeof(double));
	covar = (double *)calloc(numsec*numsec, sizeof(double));

	/* compute daily return */
	outputFile3 = fopen("return.csv", "w");
	for (int k = 0; k < numsec; k++) {
		for (int n = 0; n < numday-1; n++) {
			ret[k*(numday - 1) + n] = (p[k*numday + n + 1] - p[k*numday + n])/p[k*numday + n]*100;
			fprintf(outputFile3, "%g,", ret[k*(numday - 1) + n]);
		}
		fprintf(outputFile3, "\n");
	}



	/*int i = 0;*/
	for (int k = 0; k < numsec; k++) {
		for (int n = 0; n <= rolldays - 1; n++) {
			sum[k] += ret[k*(numday-1) + n];
		}
		mean[k] = sum[k] / rolldays;
	}


	for (int k = 0; k < numsec; k++) {
		for (int j = 0; j < numsec; j++) {
			for (int n = 0; n <= rolldays - 1; n++) {
				covarsum[numsec*k + j] += (ret[k*(numday-1) + n] - mean[k])*(ret[j*(numday - 1) + n] - mean[j]);

			}
			covar[numsec*k + j] = covarsum[numsec*k + j] / (rolldays - 1);
		}
	}

	eigenvalue = (double *)calloc(numday - rolldays + 1, sizeof(double));

	eigenvalue[0] = power(numsec, covar);
	printf("eigenvalue is: %g\n", eigenvalue[0]);

	for (int i = 1; i < numday - rolldays; i++) {
		/* compute mean for each asset over 'rolldays' days*/
		for (int k = 0; k < numsec; k++) {
			mean[k] = mean[k] + (ret[k*(numday - 1) + i - 1 + rolldays] - ret[k*(numday - 1) + i - 1]) / rolldays;
		}

		/* compute covariance over 'rolldays' days*/
		for (int k = 0; k < numsec; k++) {
			for (int j = 0; j < numsec; j++) {
				covar[numsec*k + j] =
					covar[numsec*k + j] + ((ret[k*(numday - 1) + i - 1 + rolldays]
						- mean[k])*(ret[j*(numday - 1) + i - 1 + rolldays] - mean[j])
						- (ret[k*(numday - 1) + i - 1] - mean[k])*(ret[j*(numday - 1) + i - 1]
							- mean[j])) / (rolldays - 1);
			}

		}
		eigenvalue[i] = power(numsec, covar);
		printf("eigenvalue is: %g\n", eigenvalue[i]);
	}


	outputFile = fopen("eigenvalue.txt", "w");
	for (int i = 0; i < numday - rolldays ; i++) {
		fprintf(outputFile, "%g\n", eigenvalue[i]);
	}
	fclose(outputFile);

	/*print out to check whether covariance is computed correctly or not; Checked correct*/
	/*no need to print;*/
	outputFile2 = fopen("covariance.csv", "w");
	for (int k = 0; k < numsec; k++) {
		for(int j=0;j<numsec;j++){
			fprintf(outputFile2, "%g,", covar[numsec*k + j]);
		}
		fprintf(outputFile2, "\n");
	}

	printf("Done!");
	return retcode;

BACK:
	return retcode;
}

double power(int numsec, double *covariance) {
	
	int returncode = 0, i, j;
	double *new_array = NULL, *diff_array = NULL, *rowmult = NULL, *rowmultnew = NULL, *random_array = NULL;
	double absw;
	double abswsum, diffsum; //abs value sum of w
	new_array = (double *)calloc(numsec, sizeof(double));
	diff_array = (double *)calloc(numsec, sizeof(double));
	rowmult = (double *)calloc(numsec, sizeof(double));
	rowmultnew = (double *)calloc(numsec, sizeof(double));

	abswsum = 0.00;
	absw = 0.00;
	diffsum = 0.00;

	/* create random w */
	random_array = (double *)calloc(numsec, sizeof(double));
	srand(1);
	for (i = 0; i < numsec; i++) {
		random_array[i] = rand() % 10 +1;  //random 1 to 10 
		//        random_array[i] = 1;
		abswsum += random_array[i] * random_array[i];

	} //normalize w//
	absw = sqrt(abswsum);
	for (i = 0; i < numsec; i++) {
		new_array[i] = random_array[i] / absw;
		//        new_array[i] = random_array[i];
	}

	int counter;
	double lambda;
	/* iteration to find eigenvalue & eigenvector */
	for (counter = 0; counter < 200; counter++) {
		//        rowmult[i] = 0;
		printf("iteration #: %i\n", counter);
		// printf("new_arrary: %f\n", new_array[0]);
		for (i = 0; i < numsec; i++) {
			for (j = 0; j < numsec; j++) {
				/*   Q*W   */
				rowmult[i] += covariance[i*numsec + j] * new_array[j];
			}
		}
		/*normalize*/
		abswsum = 0.0;
		absw = 0.0;
		for (i = 0; i < numsec; i++) {
			abswsum += rowmult[i] * rowmult[i];
		}
		absw = sqrt(abswsum);

		diffsum = 0.0;

		for (i = 0; i < numsec; i++) {

			rowmultnew[i] = rowmult[i] / absw;
			/* compare the new array with the old */
			diff_array[i] = new_array[i] - rowmultnew[i];
			// printf("new_array %f\n", new_array[i]);
			diffsum += diff_array[i] * diff_array[i];

		}

	
		/*pass value to new_array */
		for (i = 0; i < numsec; i++) {
			new_array[i] = rowmultnew[i];
		}

		/*stop iteration when similar*/
		if (sqrt(diffsum) <= 0.0001) {
			break;

		}
		//        printf("diffsum is %f\n", diffsum);

		/*calculate lambda*/
		lambda = 0.0;
		for (i = 0; i < numsec; i++) {
			for (j = 0; j < numsec; j++) {

				lambda += new_array[i] * covariance[i*numsec + j] * new_array[j];  //wt*Q*W

			}
		}

		printf("lambda is %g\n", lambda);


	}

	return lambda; 
}


	
