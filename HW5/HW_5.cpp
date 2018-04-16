
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <sys/types.h>
#include <sys/stat.h>

#include <float.h>
#include <cmath>

#include <gurobi_c++.h>

using namespace std;

double *p;
int read(int numsec, int numday, int colnum, int year, char *filename);

double *covar = NULL, *covar_origin = NULL, *mean = NULL;
int covariance(int numsec, int numday, double *p);

double *eigenvalue = NULL, *eigenvector = NULL;
double power(int numsec, int numprinc, double *covariance);

double *x;
int numname = 0;
int portfolio_optimize(int numsec, int numprinc, int type, int namelimit, double *qval);
int longshort(int numsec, int numprinc, int type, int namelimit, double *qval, double *xstart);


int main(int argc, char **argv)
{
	/* a) */
	FILE *in = NULL, *outputFile = NULL;
	/*int k, n, j;*/
	int numsec, numday, colnum,numprinc;
	int retcode = 0, year;
	char mybuffer[200];
	double *ret, *sum;
	FILE *risklevel = NULL;
	char *filename;
	numprinc = 10;

	if (argc != 5) {
		printf("Usage:  Eigen10.exe datafilename numsec numday colnum\n"); retcode = 100; goto BACK;
	}

	filename = argv[1];
	numsec = atoi(argv[2]);
	numday = atoi(argv[3]);
	colnum = atoi(argv[4]);

	/* read first year data */
	year = 1;
	read(numsec, numday, colnum, year, filename);
	
	/* compute covariance matrix */
	covariance(numsec, numday, p);

	/* implement the power method */
	power(numsec, numprinc, covar);

	outputFile = fopen("eigenvalue.txt", "w");
	for (int k = 0; k < numprinc; k++) {
		fprintf(outputFile, "%g\n", eigenvalue[k]);
	}
	fclose(outputFile);

	printf("Done with PCA computation!");

	/* start portfolio optimization */
	int type, namelimit, j, Nq;
	double risklambda, newrisklambda;
	double *qval;
	
	/* b) */
	type = 0;
	namelimit = numsec;
	
	Nq = numsec + numprinc;

	/* start the risk aversion level from 1.0 */
	risklambda = 1.0; newrisklambda = 1.0;

	/* the coefficients of the quardratic terms at initial risk aversion level */
	qval = (double *)calloc(Nq, sizeof(double));  /** values **/
	 
	for (j = 0; j < numsec; j++) qval[j] = covar_origin[j*numsec + j] * risklambda;
	for (j = numsec; j < Nq; j++) qval[j] = eigenvalue[j - numsec]* risklambda;
	
	/* try different risk aversion levels until the number of names is at least 200*/
	while (numname < 200) {
		for (j = 0; j < Nq; j++) {
			qval[j] *= newrisklambda / risklambda;
		}
		numname = portfolio_optimize(numsec, numprinc, type, namelimit, qval);
		risklambda = newrisklambda;
		newrisklambda += 0.1; /* change by 0.1 every time */
	}
	
	risklevel = fopen("risklevel.txt", "w");	
	fprintf(risklevel, "%g\n", risklambda);

	double *xstart;
	xstart = (double *)calloc(numsec, sizeof(double));
	for (int k = 0; k < numsec; k++) {
		xstart[k] = x[k];
	}


	/* c) */
	type = 1; /*add names constraints*/
	namelimit = 100; /* start from limit = 100 */

	while (namelimit >= 50) {
		numname = portfolio_optimize(numsec, numprinc, type, namelimit, qval);
		namelimit -= 10;
	}

	/* d) extra credit */

	/* read second year data */
	year = 2;
	read(numsec, numday, colnum, year, filename);

	/* compute covariance matrix */
	covariance(numsec, numday, p);

	/* implement the power method */
	power(numsec, numprinc, covar);

	/* output eigenvalue */
	outputFile = fopen("eigenvalue2.txt", "w");
	for (int k = 0; k < numprinc; k++) {
		fprintf(outputFile, "%g\n", eigenvalue[k]);
	}
	fclose(outputFile);
	
	/* compute new qval */
	for (j = 0; j < numsec; j++) qval[j] = covar_origin[j*numsec + j] * risklambda;
	for (j = numsec; j < Nq; j++) qval[j] = eigenvalue[j - numsec] * risklambda;

	type = 0;
	namelimit = numsec;

	/* longshort portfolio optimization */
	longshort(numsec, numprinc, type, namelimit, qval, xstart);

	return retcode;

BACK:
	return retcode;
}

int portfolio_optimize(int numsec,int numprinc,int type,int namelimit,double *qval) {
	int retcode = 0;
	GRBenv   *env = NULL;
	GRBmodel *model = NULL;
	int n, j;
	double *obj = NULL;
	double *lb = NULL;
	double *ub = NULL;
	int *qrow, *qcol, Nq;
	int *cind;
	double rhs;
	char sense;
	double *cval;
	int numnonz;

	char **names, *vartype;
	
	n = numsec + numprinc+ numsec*type; 
	/* numsec 'x' variables, numprinc 'y' variables,numsec*type 'z' binary variables */
    /* type: 1 if add names constraints */

	retcode = GRBloadenv(&env, "portfolio.log");
	if (retcode) goto BACK;

	/* Create initial model */
	retcode = GRBnewmodel(env, &model, "portfolio", n, NULL, NULL, NULL, NULL, NULL);
	if (retcode) goto BACK;

	names = (char **)calloc(n, sizeof(char *));

	/** next we create the remaining attributes for the n columns **/
	obj = (double *)calloc(n, sizeof(double));
	ub = (double *)calloc(n, sizeof(double));
	lb = (double *)calloc(n, sizeof(double));
	x = (double *)calloc(n, sizeof(double));
	vartype = (char *)calloc(n, sizeof(char));

	for (j = 0; j < numsec; j++) {
		names[j] = (char *)calloc(4, sizeof(char));
		sprintf(names[j], "x%d", j);
		lb[j] = 0.0;
		ub[j] = 0.02;
		obj[j] =-mean[j];
	}

	for (j = numsec; j < numsec+numprinc; j++) {
		names[j] = (char *)calloc(4, sizeof(char));
		if (names[j] == NULL) {
			retcode = 1; goto BACK;
		}
		sprintf(names[j], "y%d", j - numsec);
		lb[j] = -100.0;
		ub[j] = 100.0;
		obj[j] = 0.0;
	}

	/* Binary variables */
	if (type == 1) {
		for (j = numsec + numprinc; j < n; j++) {
			names[j] = (char *)calloc(4, sizeof(char));
			if (names[j] == NULL) {
				retcode = 1; goto BACK;
			}
			sprintf(names[j], "z%d", j - numsec - numprinc);
			lb[j] = 0.0; /* redundant because of calloc */
			ub[j] = 1.0;
			obj[j] = 0.0; /* redundant, again */
		}
	}
	
	/* initialize variables */
	for (j = 0; j < n; j++) {
		retcode = GRBsetstrattrelement(model, "VarName", j, names[j]);
		if (retcode) goto BACK;
		
		retcode = GRBsetdblattrelement(model, "Obj", j, obj[j]);
		if (retcode) goto BACK;
		
		retcode = GRBsetdblattrelement(model, "LB", j, lb[j]);
		if (retcode) goto BACK;
		
		retcode = GRBsetdblattrelement(model, "UB", j, ub[j]);
		if (retcode) goto BACK;
		
		if (j < numsec + numprinc) vartype[j] = GRB_CONTINUOUS;
		else vartype[j] = GRB_BINARY;
		
		retcode = GRBsetcharattrelement(model, "VTYPE", j, vartype[j]);
		if (retcode) goto BACK;

	}

	/** next, the quadratic -- there are numsec^2 terms **/

	Nq = numsec+numprinc;
	qrow = (int *)calloc(Nq, sizeof(int));  /** row indices **/
	qcol = (int *)calloc(Nq, sizeof(int));  /** column indices **/

	if ((qrow == NULL) || (qcol == NULL) || (qval == NULL)) {
		printf("could not create quadratic\n");
		retcode = 1; goto BACK;
	}

	for (j = 0; j < Nq; j++) {
			qrow[j] = j;
			qcol[j] = j;		
		}
	
	retcode = GRBaddqpterms(model, Nq, qrow, qcol, qval);
	if (retcode) goto BACK;


	/** now we will add one constraint at a time **/
	/** we need to have a couple of auxiliary arrays **/

	cind = (int *)calloc(n, sizeof(int));  /** n is over the top since no constraint is totally dense;
											but it's not too bad here **/
	cval = (double *)calloc(n, sizeof(double));

	/** sum of x variables = 1 **/
	for (j = 0; j < numsec; j++) {
		cval[j] = 1.0;
		cind[j] = j;
	}

	numnonz = numsec;
	rhs = 1.0;
	sense = GRB_EQUAL;

	retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, "convexity");
	printf("retcode2: %d", retcode);
	if (retcode) goto BACK;

	/** y equals V*x  **/
	for (int i = 0; i < numprinc; i++) {
		for (j = 0; j < numsec; j++) {
			cind[j] = j;
			cval[j] = -eigenvector[i*numsec + j];
		}
		cval[numsec] = 1.0;
		cind[numsec] = numsec + i;

		numnonz = numsec + 1;
		rhs = 0.0;
		sense = GRB_EQUAL;

		retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, NULL);
		if (retcode) goto BACK;
		}
	

	/** sum of z variables <= namelimit **/
	if (type == 1) {
		for (j = numsec + numprinc; j < n; j++) {
			cval[j - numsec - numprinc] = 1.0;
			cind[j - numsec - numprinc] = j;
		}
		numnonz = numsec;
		rhs = namelimit;
		sense = GRB_LESS_EQUAL;

		retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, "namelimit_constraint");
		if (retcode) goto BACK;

		/** now say xj = 0 unless zj = 1 **/
		for (j = 0; j < numsec; j++) {
			cval[0] = 1.0;  cind[0] = j;
			cval[1] = -0.02;  cind[1] = numsec + numprinc + j;

			numnonz = 2;
			rhs = 0.0;
			sense = GRB_LESS_EQUAL;

			/* let's reuse some space */
			sprintf(names[0], "control%d", j);

			retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, names[0]);
			if (retcode) goto BACK;
		}		
	}

	retcode = GRBupdatemodel(model);
	if (retcode) goto BACK;

	/** optional: write the problem **/

	retcode = GRBwrite(model, "portfolio.lp");
	if (retcode) goto BACK;

	retcode = GRBoptimize(model);
	if (retcode) goto BACK;

	/** get solution **/

	retcode = GRBgetdblattrarray(model,
		GRB_DBL_ATTR_X, 0, n, x);
	if (retcode) goto BACK;

	/** now let's see the values **/

	/*for (j = 0; j < numsec; j++) {
		printf("%s = %g\n", names[j], x[j]);
	}*/
	
	numname = 0;
	for (j = 0; j < numsec; j++) {
		if (x[j] > 0.000000005) {
			numname += 1;
		}
	}	

	FILE *outputFile5 = NULL;
	char buf[100];
	char name[100];
	sprintf(name, "%d", namelimit);
	sprintf(buf, "HW5_portfolio_%s.csv", name);
	
	outputFile5 = fopen(buf, "w");
	fprintf(outputFile5, "Limit on Number of Names, %d\n", namelimit);
	fprintf(outputFile5, "Number of Names, %d\n", numname);

	for (j = 0; j < numsec; j++) {
		fprintf(outputFile5, "%s, %g\n", names[j], x[j]);
	}
	return numname;

	GRBfreeenv(env);
	
BACK:
	printf("\nexiting with retcode %d\n", retcode);
	return retcode;
}

double power(int numsec, int numprinc, double *covariance) {

	int returncode = 0, i, j;
	double *new_array = NULL, *diff_array = NULL, *rowmult = NULL, *rowmultnew = NULL, *random_array = NULL, *initial_array = NULL;
	double absw;
	double abswsum, diffsum; //abs value sum of w
	abswsum = 0.00;
	absw = 0.00;
	diffsum = 0.00;
	//double A[10];  //store eigenvalues
	// double hold[numsec];

	new_array = (double *)calloc(numsec, sizeof(double));
	diff_array = (double *)calloc(numsec, sizeof(double));
	rowmult = (double *)calloc(numsec, sizeof(double));
	rowmultnew = (double *)calloc(numsec, sizeof(double));

	eigenvalue = (double *)calloc(numprinc, sizeof(double));
	eigenvector = (double *)calloc(numprinc * numsec, sizeof(double));

	/* create random w */
	random_array = (double *)calloc(numsec, sizeof(double));
	srand(1);
	for (i = 0; i < numsec; i++) {
		random_array[i] = rand() % numprinc + 1;  //random 1 to numprinc
											//        random_array[i] = 1;
		abswsum += random_array[i] * random_array[i];

	} //normalize w//
	absw = sqrt(abswsum);
	for (i = 0; i < numsec; i++) {
		new_array[i] = random_array[i] / absw;
		//        new_array[i] = random_array[i];
	}

	initial_array = (double *)calloc(numsec, sizeof(double));
	for (i = 0; i < numsec; i++) {
		initial_array[i] = new_array[i];
		//printf("initial_first: %g\n", initial_array[i]);
		//        new_array[i] = random_array[i];
	}


	int bigcounter;
	//double lambda;
	for (bigcounter = 0; bigcounter < numprinc; bigcounter++) {   //10 eigenvalues
		printf("bigiteration #: %i\n", bigcounter);


		int counter;
		double lambda;
		/* iteration to find eigenvalue & eigenvector */
		for (counter = 0; counter < 1000; counter++) {
			for (i = 0; i < numsec; i++) {
				rowmult[i] = 0;
			}
			//printf("iteration #: %i\n", counter);
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
				diff_array[i] = 0;
			}

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
			if (sqrt(diffsum) <= 0.000001) {
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

			//printf("lambda is %g\n", lambda);
		}
		eigenvalue[bigcounter] = lambda;
		printf("lambda is %g\n", lambda);

		for (int i = 1; i < numsec + 1; i++) {
			for (int j = 0; j < numsec; j++) {
				eigenvector[bigcounter*i + j] = new_array[j];
			}
		}
		//printf("newarray(1) is %g\n", new_array[0]);
		//printf("eigenvec(1) is %g\n", Eigenvec_Array[30]);

		for (int i = 0; i < numsec; i++) {
			for (int j = 0; j < numsec; j++) {
				covariance[i*numsec + j] -= lambda * new_array[i] * new_array[j];
				//covariance[i*numsec + j] =  fabs(covariance[i*numsec + j]);
			}
		}
		double wT_w0 = 0;

		for (int i = 0; i < numsec; i++)
		{
			wT_w0 += new_array[i] * initial_array[i];
		}

		for (int i = 0; i < numsec; i++)
		{
			new_array[i] = initial_array[i] - wT_w0 * new_array[i];
		}


		for (int i = 0; i < numsec; i++)
		{
			initial_array[i] = new_array[i];
		}
	}
	printf("first_eigenvalue is %g\n", eigenvalue[0]);
	return returncode;
}

int covariance(int numsec, int numday, double *p) {
	 
	int retcode = 0;
	double *ret, *sum, *covarsum;
	FILE *outputFile2, *outputFile3, *outputFile4;

	ret = (double *)calloc(numsec*(numday - 1), sizeof(double));
	sum = (double *)calloc(numsec, sizeof(double));
	mean = (double *)calloc(numsec, sizeof(double));
	covarsum = (double *)calloc(numsec*numsec, sizeof(double));
	covar = (double *)calloc(numsec*numsec, sizeof(double));
	covar_origin = (double *)calloc(numsec*numsec, sizeof(double));

	/* compute daily return */
	outputFile3 = fopen("return.csv", "w");
	for (int k = 0; k < numsec; k++) {
		for (int n = 0; n < numday - 1; n++) {
			ret[k*(numday - 1) + n] = (p[k*numday + n + 1] - p[k*numday + n]) / p[k*numday + n];
			fprintf(outputFile3, "%g,", ret[k*(numday - 1) + n]);
		}
		fprintf(outputFile3, "\n");
	}

	/*int i = 0;*/
	for (int k = 0; k < numsec; k++) {
		for (int n = 0; n < numday; n++) {
			sum[k] += ret[k*(numday - 1) + n];
		}
		mean[k] = sum[k] / numday;
	}

	outputFile4 = fopen("avg_return.csv", "w");
	for (int k = 0; k < numsec; k++) {
		fprintf(outputFile4, "%g,", mean[k]);
		fprintf(outputFile4, "\n");
	}

	for (int k = 0; k < numsec; k++) {
		for (int j = 0; j < numsec; j++) {
			for (int n = 0; n < numday; n++) {
				covarsum[numsec*k + j] += (ret[k*(numday - 1) + n] - mean[k])*(ret[j*(numday - 1) + n] - mean[j]);

			}
			covar[numsec*k + j] = covarsum[numsec*k + j] / (numday - 1);
		}
	}

	outputFile2 = fopen("covariance.csv", "w");
	for (int k = 0; k < numsec; k++) {
		for (int j = 0; j<numsec; j++) {
			fprintf(outputFile2, "%g,", covar[numsec*k + j]);
		}
		fprintf(outputFile2, "\n");
	}

	for (int i = 0; i < numsec*numsec; i++) {
		covar_origin[i] = covar[i];
	}

	return retcode;
}

int longshort(int numsec, int numprinc, int type, int namelimit, double *qval,double *xstart) {
	int retcode = 0;
	GRBenv   *env = NULL;
	GRBmodel *model = NULL;
	int n, j;
	double *obj = NULL;
	double *lb = NULL;
	double *ub = NULL;
	int *qrow, *qcol, Nq;
	int *cind;
	double rhs;
	char sense;
	double *cval;
	int numnonz;

	char **names, *vartype;

	n = numsec*5 + numprinc + numsec*type;
	/* numsec 'x' variables, numprinc 'y' variables,numsec*type 'z' binary variables */
	/* type: 1 if add names constraints */

	retcode = GRBloadenv(&env, "portfolio.log");
	if (retcode) goto BACK;

	/* Create initial model */
	retcode = GRBnewmodel(env, &model, "portfolio", n, NULL, NULL, NULL, NULL, NULL);
	if (retcode) goto BACK;

	names = (char **)calloc(n, sizeof(char *));

	/** next we create the remaining attributes for the n columns **/
	obj = (double *)calloc(n, sizeof(double));
	ub = (double *)calloc(n, sizeof(double));
	lb = (double *)calloc(n, sizeof(double));
	x = (double *)calloc(n, sizeof(double));
	vartype = (char *)calloc(n, sizeof(char));

	/* weight variables x */
	for (j = 0; j < numsec; j++) {
		names[j] = (char *)calloc(4, sizeof(char));
		sprintf(names[j], "x%d", j);
		lb[j] = -100.0;
		ub[j] = 0.02;
		obj[j] = -mean[j];
	}
	printf("1");

	/* principla variables y*/
	for (j = numsec; j < numsec + numprinc; j++) {
		names[j] = (char *)calloc(4, sizeof(char));
		if (names[j] == NULL) {
			retcode = 1; goto BACK;
		}
		sprintf(names[j], "y%d", j - numsec);
		lb[j] = -100.0;
		ub[j] = 100.0;
		obj[j] = 0.0;
	}
	printf("2");
	/* Binary variables */
	if (type == 1) {
		for (j = numsec + numprinc; j < numsec + numprinc + numsec*type; j++) {
			names[j] = (char *)calloc(4, sizeof(char));
			if (names[j] == NULL) {
				retcode = 1; goto BACK;
			}
			sprintf(names[j], "z%d", j - numsec - numprinc);
			lb[j] = 0.0; /* redundant because of calloc */
			ub[j] = 1.0;
			obj[j] = 0.0; /* redundant, again */
		}
	}
	printf("3");
	/* p,m,t variables*/
	for (j = numsec + numprinc + numsec*type; j < numsec*2 + numprinc + numsec*type; j++) {
		names[j] = (char *)calloc(4, sizeof(char));
		if (names[j] == NULL) {
			retcode = 1; goto BACK;
		}
		sprintf(names[j], "p%d", j - numsec - numprinc - numsec*type);
		lb[j] = 0.0;
		ub[j] = 100.0;
		obj[j] = 0.0;
	}
	printf("4");
	for (j = numsec*2 + numprinc + numsec*type; j < numsec * 3 + numprinc + numsec*type; j++) {
		names[j] = (char *)calloc(4, sizeof(char));
		if (names[j] == NULL) {
			retcode = 1; goto BACK;
		}
		sprintf(names[j], "m%d", j - numsec*2 - numprinc - numsec*type);
		lb[j] = 0.0;
		ub[j] = 0.6;
		obj[j] = 0.0;
	}
	printf("5");
	for (j = numsec * 3 + numprinc + numsec*type; j < numsec * 4 + numprinc + numsec*type; j++) {
		names[j] = (char *)calloc(4, sizeof(char));
		if (names[j] == NULL) {
			retcode = 1; goto BACK;
		}
		sprintf(names[j], "t%d", j - numsec * 3 - numprinc - numsec*type);
		lb[j] = 0.0;
		ub[j] = 100.0;
		obj[j] = 0.0;
	}
	printf("6");

	/* new binary variable */
	for (j = numsec * 4 + numprinc + numsec*type; j < numsec * 5 + numprinc + numsec*type; j++) {
		names[j] = (char *)calloc(4, sizeof(char));
		if (names[j] == NULL) {
			retcode = 1; goto BACK;
		}
		sprintf(names[j], "a%d", j - numsec * 3 - numprinc - numsec*type);
		lb[j] = 0.0;
		ub[j] = 1.0;
		obj[j] = 0.0;
	}
	printf("7");

	/* initialize variables */
	for (j = 0; j < n; j++) {
		retcode = GRBsetstrattrelement(model, "VarName", j, names[j]);
		if (retcode) goto BACK;

		retcode = GRBsetdblattrelement(model, "Obj", j, obj[j]);
		if (retcode) goto BACK;

		retcode = GRBsetdblattrelement(model, "LB", j, lb[j]);
		if (retcode) goto BACK;

		retcode = GRBsetdblattrelement(model, "UB", j, ub[j]);
		if (retcode) goto BACK;

		if (((j > numsec + numprinc) & (j<numsec + numprinc+numsec*type))|| (j>numsec*4 + numprinc + numsec*type)) vartype[j] = GRB_BINARY;
		else vartype[j] = GRB_CONTINUOUS;

		retcode = GRBsetcharattrelement(model, "VTYPE", j, vartype[j]);
		if (retcode) goto BACK;

	}
	printf("8");

	/** next, the quadratic -- there are numsec^2 terms **/

	Nq = numsec + numprinc;
	qrow = (int *)calloc(Nq, sizeof(int));  /** row indices **/
	qcol = (int *)calloc(Nq, sizeof(int));  /** column indices **/
	printf("9");
	if ((qrow == NULL) || (qcol == NULL) || (qval == NULL)) {
		printf("could not create quadratic\n");
		retcode = 1; goto BACK;
	}

	for (j = 0; j < Nq; j++) {
		qrow[j] = j;
		qcol[j] = j;
	}

	retcode = GRBaddqpterms(model, Nq, qrow, qcol, qval);
	if (retcode) goto BACK;
	printf("10");

	/** now we will add one constraint at a time **/
	/** we need to have a couple of auxiliary arrays **/

	cind = (int *)calloc(n, sizeof(int));  /** n is over the top since no constraint is totally dense;
										   but it's not too bad here **/
	cval = (double *)calloc(n, sizeof(double));

	/** sum of x variables = 1 **/
	for (j = 0; j < numsec; j++) {
		cval[j] = 1.0;
		cind[j] = j;
	}

	numnonz = numsec;
	rhs = 1.0;
	sense = GRB_EQUAL;

	retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, "convexity");
	if (retcode) goto BACK;
	printf("11");

	/** y equals V*x  **/
	for (int i = 0; i < numprinc; i++) {
		for (j = 0; j < numsec; j++) {
			cind[j] = j;
			cval[j] = -eigenvector[i*numsec + j];
		}
		cval[numsec] = 1.0;
		cind[numsec] = numsec + i;

		numnonz = numsec + 1;
		rhs = 0.0;
		sense = GRB_EQUAL;

		retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, NULL);
		if (retcode) goto BACK;
	}

	printf("12");

	/** sum of z variables <= namelimit **/
	if (type == 1) {
		for (j = numsec + numprinc; j < numsec + numprinc + numsec*type; j++) {
			cval[j - numsec - numprinc] = 1.0;
			cind[j - numsec - numprinc] = j;
		}
		numnonz = numsec;
		rhs = namelimit;
		sense = GRB_LESS_EQUAL;

		retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, "namelimit_constraint");
		if (retcode) goto BACK;

		/** now say xj = 0 unless zj = 1 **/
		for (j = 0; j < numsec; j++) {
			cval[0] = 1.0;  cind[0] = j;
			cval[1] = -0.02;  cind[1] = numsec + numprinc + j;

			numnonz = 2;
			rhs = 0.0;
			sense = GRB_LESS_EQUAL;

			/* let's reuse some space */
			sprintf(names[0], "control%d", j);

			retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, names[0]);
			if (retcode) goto BACK;
		}
	}
	printf("13");

	/* sum of long position p - t less than 1.6 - xstart */

	/*for (j = numsec + numprinc + numsec*type; j < numsec * 2 + numprinc + numsec*type; j++) {
		cval[j] = 1.0;
		cind[j] = j;
	}

	numnonz = numsec;
	rhs = 0.6;
	sense = GRB_LESS_EQUAL;

	retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, NULL);
	if (retcode) goto BACK;*/

	/* sum of short position mi greater than 0.4 and less than 0.6 */
	for (j = numsec * 2 + numprinc + numsec*type; j < numsec * 3 + numprinc + numsec*type; j++) {
		cval[j - (numsec * 2 + numprinc + numsec*type)] = 1.0;
		cind[j - (numsec * 2 + numprinc + numsec*type)] = j;
	}

	numnonz = numsec;
	rhs = 0.6;
	sense = GRB_LESS_EQUAL;

	retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs,NULL);
	if (retcode) goto BACK;
	printf("14");

	for (j = numsec * 2 + numprinc + numsec*type; j < numsec * 3 + numprinc + numsec*type; j++) {
		cval[j - (numsec * 2 + numprinc + numsec*type)] = 1.0;
		cind[j - (numsec * 2 + numprinc + numsec*type)] = j;
	}

	numnonz = numsec;
	rhs = 0.4;
	sense = GRB_GREATER_EQUAL;

	retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, NULL);
	if (retcode) goto BACK;
	printf("15");

	/* sum of ti less than 0.2 */
	for (j = numsec * 3 + numprinc + numsec*type; j < numsec * 4 + numprinc + numsec*type; j++) {
		cval[j - (numsec * 3 + numprinc + numsec*type)] = 1.0;
		cind[j - (numsec * 3 + numprinc + numsec*type)] = j;
	}

	numnonz = numsec;
	rhs = 0.2;
	sense = GRB_LESS_EQUAL;

	retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, NULL);
	if (retcode) goto BACK;
	printf("16");

	/* xi - pi + mi + ti = xstarti */
	for (int j = 0; j < numsec; j++) {
		cind[0] = j; cval[0] = 1.0;
		cind[1] = j + numsec + numprinc + numsec*type; cval[1] = -1.0;
		cind[2] = j + numsec * 2 + numprinc + numsec*type; cval[2] = 1.0;
		cind[3] = j + numsec * 3 + numprinc + numsec*type; cval[3] = 1.0;	

		numnonz = 4;
		rhs = xstart[j];
		sense = GRB_EQUAL;

		retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, NULL);
		if (retcode) goto BACK;
	}
	printf("17");

	double limit = 0.00000005;

	/* xstarti*ai - ti <=0 for xstarti>0; pi-0.02ai<=0 for xstart = 0 */
	for (int j = 0; j < numsec; j++) {
		if (xstart[j] > 0.0) {
			cind[0] = j + numsec * 3 + numprinc + numsec*type;//ti
			cind[1] = j + numsec * 4 + numprinc + numsec*type;//ai
			cval[0] = -1.0; //ti
			cval[1] = xstart[j]; //ai
		}
		else{
			cind[0] = j + numsec + numprinc + numsec*type;//pi
			cind[1] = j + numsec * 4 + numprinc + numsec*type;//ai
			cval[0] = 1.0; //pi
			cval[1] = -0.02; //ai
		}
		numnonz = 2;
		rhs = 0.0;
		sense = GRB_LESS_EQUAL;

		retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, NULL);
		if (retcode) goto BACK;	
	}
	printf("18");


	/* mi - 0.6*ai <=0 for xstarti>0; mi + 0.6ai<=0.6 for xstart = 0 */
	for (int j = 0; j < numsec; j++) {
		cind[0] = j + numsec * 2 + numprinc + numsec*type;
		cind[1] = j + numsec * 4 + numprinc + numsec*type;
		if (xstart[j] > 0.0) {
			cval[0] = 1.0; //mi
			cval[1] = -0.6; //ai
			rhs = 0.0;
		}
		else {
			cval[0] = 1.0; //mi
			cval[1] = 0.6; //ai
			rhs = 0.6;
		}
		numnonz = 2;
		sense = GRB_LESS_EQUAL;

		retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, NULL);
		if (retcode) goto BACK;
	}

	printf("19");

	retcode = GRBupdatemodel(model);
	if (retcode) goto BACK;

	/** optional: write the problem **/

	retcode = GRBwrite(model, "portfolio.lp");
	if (retcode) goto BACK;

	retcode = GRBoptimize(model);
	if (retcode) goto BACK;

	/** get solution **/

	retcode = GRBgetdblattrarray(model,
		GRB_DBL_ATTR_X, 0, n, x);
	if (retcode) goto BACK;

	/** now let's see the values **/

	/*for (j = 0; j < numsec; j++) {
	printf("%s = %g\n", names[j], x[j]);
	}*/

	numname = 0;
	for (j = 0; j < numsec; j++) {
		if (x[j] > 0.000000005) {
			numname += 1;
		}
	}

	FILE *outputFile6 = NULL;
    outputFile6 = fopen("HW5_portfolio_longshort.csv", "w");
	fprintf(outputFile6, "Limit on Number of Names, %d\n", namelimit);
	fprintf(outputFile6, "Number of Names, %d\n", numname);

	for (j = 0; j < numsec; j++) {
		fprintf(outputFile6, "%s, %g\n", names[j], x[j]);
	}
	return numname;

	GRBfreeenv(env);

BACK:
	printf("\nexiting with retcode %d\n", retcode);
	return retcode;
}

int read(int numsec, int numday, int colnum, int year,char *filename) {
	FILE *in;
	int retcode = 0;
	char mybuffer[200];

	in = fopen(filename, "r"); /*price datafile*/

	if (in == NULL) {
		printf("could not open %s for reading\n", filename);
		retcode = 200; goto BACK;
	}

	/*p is the price data read from the data file*/
	p = (double *)calloc(numsec*numday, sizeof(double));
	if (p == NULL) {
		printf("no memory\n"); retcode = 400; goto BACK;
	}

	//skip to read first sec's price
	for (int k = 0; k < colnum + 252*(year-1)+1; k++) {
		fscanf(in, "%s", mybuffer);
	}

	for (int k = 0; k < numsec; k++) {

		for (int n = 0; n < numday; n++) {
			fscanf(in, "%s", mybuffer);
			p[k*numday + n] = atof(mybuffer);
		}

		for (int k = 0; k < colnum - numday; k++) {  //skip data not in first year
			fscanf(in, "%s", mybuffer);
		}
	}

	printf("Finish reading price data!\n");
	fclose(in);

	return retcode;

BACK:
	return retcode;
}
