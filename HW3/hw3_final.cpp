#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <cmath>
// input: example.txt
// professor gives the first feasible solution algorithm: mainqp1
using namespace std;
typedef struct bag{
	int n;
	double *x;
	double *lb;
	double *ub;
	double *mu;
	double *covariance;
	double *gradient;
	double lambda;
	int num_iter;
	FILE *output;
	int iter_counter;
	int *sortindex;
	//double *copyofgradient;
	double *y;
	double *ymin;
	double s;
	double fmin;
};

int readit(char *nameoffile, int *addressofn, double **, double **, double **, double *, double **);

int algo(int n, double *x, double *lb, double *ub, double *mu, double *covariance, double lambda);

int improve(struct bag *pmybag);

int feasible(int n, double *x, double *lb, double *ub);

int gradient(struct bag* pmybag);

int stepdirection(struct bag* pmybag);

int stepsize(struct bag* pmybag);

int move(struct bag* pmybag);

struct gradient2darray;
struct gradient2darray {
	double grad;
	int index;
};

int partition(struct gradient2darray* A, int l, int r);

int quickSort(struct gradient2darray* A, int l, int r);

int main(int argc, char **argv)
{
	int retcode = 0;
	int n;
	double *lb, *ub, *covariance, *mu, lambda, *x;

	if (argc != 2) {
		printf("usage: qp1 filename\n");  retcode = 1;
		goto BACK;
	}

	//& means address
	retcode = readit(argv[1], &n, &lb, &ub, &mu, &lambda, &covariance);
	if (retcode) goto BACK;

	x = (double *)calloc(n, sizeof(double));
	if (x == NULL) {
		printf(" no memory for x\n"); retcode = 1; goto BACK;
	}

	retcode = algo(n, x, lb, ub, mu, covariance, lambda);
BACK:
	return retcode;
}

int readit(char *filename, int *address_of_n, double **plb, double **pub,
	double **pmu, double *plambda, double **pcovariance)
{
	int readcode = 0, fscancode;
	FILE *datafile = NULL;
	char buffer[100];
	int n, i, j;
	double *lb = NULL, *ub = NULL, *mu = NULL, *covariance = NULL;

	datafile = fopen(filename, "r");
	if (!datafile) {
		printf("cannot open file %s\n", filename);
		readcode = 2;  goto BACK;
	}

	printf("reading data file %s\n", filename);

	fscanf(datafile, "%s", buffer);
	fscancode = fscanf(datafile, "%s", buffer);
	if (fscancode == EOF) {
		printf("problem: premature file end at ...\n");
		readcode = 4; goto BACK;
	}
//n
	n = *address_of_n = atoi(buffer);

	printf("n = %d\n", n);

	lb = (double *)calloc(n, sizeof(double));
	*plb = lb;
	ub = (double *)calloc(n, sizeof(double));
	*pub = ub;
	mu = (double *)calloc(n, sizeof(double));
	*pmu = mu;
	covariance = (double *)calloc(n*n, sizeof(double));
	*pcovariance = covariance;

	if (!lb || !ub || !mu || !covariance) {
		printf("not enough memory for lb ub mu covariance\n"); readcode = 3; goto BACK;
	}
//skip j_lower_upper_mu
	fscanf(datafile, "%s", buffer);

	for (j = 0; j < n; j++) {
		fscanf(datafile, "%s", buffer);
		fscanf(datafile, "%s", buffer);
		lb[j] = atof(buffer);
		fscanf(datafile, "%s", buffer);
		ub[j] = atof(buffer);
		fscanf(datafile, "%s", buffer);
		mu[j] = atof(buffer);
		printf("j = %d lb = %g ub = %g mu = %g\n", j, lb[j], ub[j], mu[j]);
	}

//lambda
	fscanf(datafile, "%s", buffer);
	fscanf(datafile, "%s", buffer);
	*plambda = atof(buffer);

	fscanf(datafile, "%s", buffer); /* reading 'covariance'*/

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			fscanf(datafile, "%s", buffer);
			covariance[i*n + j] = atof(buffer);
		}
	}

//strcmp compares two strings and returns 0 if they are the same 
	fscanf(datafile, "%s", buffer);
	if (strcmp(buffer, "END") != 0) {
		printf("possible error in data file: 'END' missing\n");
	}
x

	fclose(datafile);

BACK:

	return readcode;
}


int algo(int n, double *x, double *lb, double *ub, double *mu, double *covariance, double lambda)
{
	int returncode = 0;
	double *gradient = NULL;// *copyofgradient;
	struct bag *pmybag = NULL;

	printf("\n running algorithm\n");

	returncode = feasible(n, x, lb, ub);
	// added on 1014
	if (returncode) {
		goto BACK;
	}

	gradient = (double*)calloc(n, sizeof(double));
	// if failed to allocate memory
	if (!gradient) {
		returncode = 100; // the reason to assign different code is for easy identificaition of error
		goto BACK;
	}
	//copyofgradient = gradient + n;

	pmybag = (struct bag*)calloc(1, sizeof(struct bag));

	if (!pmybag) {
		returncode = 100; goto BACK;
	}

	// put things into bag
	pmybag->n = n; ///*copy n on the rhs to the left; go to where pmybag store n, address*/
	pmybag->fmin = DBL_MAX;
	pmybag->covariance = covariance;
	pmybag->gradient = gradient;
	//pmybag->copyofgradient = copyofgradient;
	pmybag->lambda = lambda;
	pmybag->lb = lb;
	pmybag->mu = mu;
	pmybag->ub = ub;
	pmybag->x = x;

	pmybag->sortindex = (int*)calloc(n, sizeof(int));
	pmybag->y = (double*)calloc(n, sizeof(double));
	pmybag->ymin = (double*)calloc(n, sizeof(double));


	pmybag->num_iter = 2000;

	pmybag->output = fopen("theoutput.dat", "w");

	if (!pmybag->output) {
		printf("could not open the file");
		returncode = 200; goto BACK;

	}

	returncode = improve(pmybag);

BACK:
	return returncode;
}

int improve(struct bag *pmybag) /*find gradi; find direction ; first order algo; how far you move stepsize==>one iteration*/
{
	int counter;
	int returncode = 0;
	for (counter = 0; counter < pmybag->num_iter; counter++) {
		fprintf(pmybag->output, "iteration # %d\n", counter);
		printf("iteration # %d\n", counter);

		pmybag->iter_counter = counter;
		/* compute gradient */
		//(returncode) => returncode != 0 
		returncode = gradient(pmybag);
		if (returncode) goto BACK;

		/* compute step direction */
		returncode = stepdirection(pmybag);
		if (returncode) goto BACK;

		/* compute stepsize */
		returncode = stepsize(pmybag);
		if (returncode) goto BACK;

		/* move in the descent direction by stepsize */
		returncode = move(pmybag);
		if (returncode) goto BACK;

		fprintf(pmybag->output, "done with iteration # %d\n\n", counter);
		printf("done with iteration # %d\n\n", counter);

	}
BACK:
	FILE *outputFile = NULL;
	outputFile = fopen("HW3_bigdata.csv", "w");
	fprintf(outputFile, "Asset, Weight\n");
	for (int i = 0; i < pmybag->n; i++) {
		//printf("x%d: %g\n", i, pmybag->x[i]);
		fprintf(outputFile, "x%d,%g\n", i, pmybag->x[i]);
	}

	return returncode;
}

int gradient(struct bag* pmybag)
{
	int returncode = 0, j, n = pmybag->n, i;
	double lambda = pmybag->lambda, first, second, third;


	for (j = 0; j < n; j++) {
		first = 0;
		first = 2 * lambda*pmybag->x[j] * pmybag->covariance[j*(n + 1)];

		second = 0;
		for (i = 0; i < n; i++)if (i != j) {
			second += pmybag->covariance[j*n + i] * pmybag->x[i]; //might be faster than [i*n+j]

		}
		second *= 2 * lambda;


		third = -pmybag->mu[j];

		pmybag->gradient[j] = first + second + third;
		// printf("gradient: %g\n", pmybag->gradient[j]);

	}


	//fprintf(pmybag->output, "gradient iteration # %d\n", pmybag->iter_counter);
	//printf("gradient iteration # %d\n", pmybag->iter_counter);

	return returncode;
}

int stepdirection(struct bag* pmybag)
{
	int returncode = 0;
	/* sort the gradient */
	// quick sort // involves two recursive calls //search on the website
	// need to keep the index

	struct gradient2darray *A;
	A = (gradient2darray*)calloc(pmybag->n, sizeof(gradient2darray));
	for (int i = 0; i < pmybag->n; i++) {
		A[i].grad = pmybag->gradient[i];
		A[i].index = i;
	}
	quickSort(A, 0, pmybag->n - 1);

	/*for (int i = 0; i < pmybag->n; i++) {
	printf("grad: %g\n", A[i].grad);
	printf("index: %d\n", A[i].index);
	}
	*/

	/* perform enumeration */
	double tmin = DBL_MAX;

	double sum, t;
	// index m that the optimal solution y* is given by...//emuerate to get best(lowest) g*y 
	for (int m = 0; m < pmybag->n; m++) {
		for (int j = 0; j < pmybag->n; j++) {
			if (j < m) {
				pmybag->y[A[j].index] =  pmybag->ub[A[j].index]  - pmybag->x[A[j].index];
				//printf("j<m,X%d: %g, y%d: %g\n", A[j].index,pmybag->x[A[j].index], A[j].index, pmybag->y[A[j].index]);
			}
			else if (j > m) {
				pmybag->y[A[j].index] =  pmybag->lb[A[j].index] - pmybag->x[A[j].index];
				//printf("j>m,X%d: %g, y%d: %g\n", A[j].index, pmybag->x[A[j].index], A[j].index, pmybag->y[A[j].index]);
			}


		}
		//printf("ymoriginal: %g\n", pmybag->y[A[m].index]);
		sum = 0;
		for (int j = 0; j < pmybag->n; j++) if (j != m) {
			sum += pmybag->y[A[j].index];
			// printf("sum: %g\n", sum);
		}
		// the sum of y* is 0 --> cal the ym* 
		pmybag->y[A[m].index] = -sum;
	

		t = 0;
		//if satisfy the constrain
		if (pmybag->y[A[m].index] <= pmybag->ub[A[m].index] - pmybag->x[A[m].index] &&
			pmybag->y[A[m].index] >= pmybag->lb[A[m].index] - pmybag->x[A[m].index]) {
			// the objective funtion gk*y 	
			for (int i = 0; i < pmybag->n; i++) {
				t += pmybag->gradient[i] * pmybag->y[i];
				//printf("gradient: %g\n", pmybag->gradient[i]);
				//printf("y: %g\n", pmybag->y[i]);
				//printf("t: %g\n", t);
			}
			//replace tmin with t for the best g*y we got so far ; store ymin
			if (t < tmin) {
				tmin = t;
				//printf("tmin: %g\n", tmin);
				for (int i = 0; i < pmybag->n; i++) {
					pmybag->ymin[i] = pmybag->y[i];
				}


				//printf("ymin1: %g\n", pmybag->ymin[1]);
			}

		}

		//printf("%g\n", sum);
		/*printf("%g\n", pmybag->ub[A[m].index]);
		printf("%g\n", pmybag->x[A[m].index]);
		printf("%g\n", pmybag->y[A[m].index]);*/



	}
	//after getting the lowest g*y tmin ; if the change smaller than 0.001 stop iteration 
	printf("tminfinal: %g\n", tmin);
	if (tmin > -0.001) {
		returncode = 1;
	}

	/*for (int j = 0; j < pmybag->n; j++) {
	printf("ymin%d: %g\n",j, pmybag->ymin[j]);
	}*/

	return returncode;
}



int quickSort(struct gradient2darray *A, int l, int r)
{
	int j;
	int returncode = 0;
	if (l < r)
	{
		// divide and conquer
		j = partition(A, l, r);
		quickSort(A, l, j - 1);
		quickSort(A, j + 1, r);
	}
	return returncode;
}

int partition(struct gradient2darray *A, int l, int r) {
	int  i, j;
	gradient2darray t;
	double pivot;
	pivot = A[l].grad;
	i = l; j = r + 1;

	while (1)
	{
		do ++i; while (A[i].grad <= pivot && i <= r);
		do --j; while (A[j].grad > pivot);
		if (i >= j) break;
		t = A[i]; A[i] = A[j]; A[j] = t;
	}
	t = A[l]; A[l] = A[j]; A[j] = t;
	return j;
}

int stepsize(struct bag* pmybag)
{
	int returncode = 0;
	double lambda = pmybag->lambda;

	int i, j, n = pmybag->n;
	//derive the G'(s) = 0 s = ...
	double first = 0, second = 0, third = 0, fourth = 0, fifth = 0;

	for (j = 0; j < n; j++) {

		first += pmybag->ymin[j] * pmybag->ymin[j] * pmybag->covariance[j*(n + 1)];
		//printf("ymin: %g\n", pmybag->ymin[j]);

		second += pmybag->mu[j] * pmybag->ymin[j] / pmybag->lambda;

		third += pmybag->covariance[j*(n + 1)] * pmybag->x[j] * pmybag->ymin[j];

		for (i = 0; i <n; i++) if(i!=j){
			fourth += pmybag->covariance[j*n + i] * pmybag->ymin[i] * pmybag->ymin[j];
			fifth += pmybag->covariance[j*n + i] * (pmybag->ymin[j] * pmybag->x[i] + pmybag->ymin[i] * pmybag->x[j])/2;

		}

	}
	pmybag->s = (second - 2 * (third + fifth)) / (first + fourth)/2;
	printf("soriginal: %g\n", pmybag->s);
	//0<s<1
	//if s <0, 0; if s>0 & s>1, 1; else s 
	pmybag->s = pmybag->s < 0 ? 0 : pmybag->s>1 ? 1 : pmybag->s;

	/*printf("first: %g\n", first);
	printf("second: %g\n", second);
	printf("third: %g\n", third);
	printf("fourth: %g\n", fourth);
	printf("fifth: %g\n", fifth);*/

	//printf("s: %g\n", pmybag->s);
	return returncode;
}

int move(struct bag* pmybag)
{
	int returncode = 0;
	int n = pmybag->n;
	//change the weight with ymin  
	for (int i = 0; i < n; i++) {
		pmybag->x[i] = pmybag->x[i] + pmybag->s*pmybag->ymin[i];

	}

	// calculate objective function value //another way to stop iteration
	// int n = pmybag->n;
	// double lambda = pmybag->lambda, first, second, third, f_new, f_incr;
	// f_new = 0;
	// f_incr = 0;

	// for (int j = 0; j < n; j++) {
	// 	first = 0;
	// 	first = lambda* pow(pmybag->x[j], 2) * pmybag->covariance[j*(n + 1)];

	// 	second = 0;
	// 	for (int i = 0; i < n; i++)if (i > j) {
	// 		second += pmybag->covariance[j*n + i] * pmybag->x[j] * pmybag->x[i]; //might be faster than [i*n+j]

	// 	}
	// 	second *= 2 * lambda ;


	// 	third = -pmybag->mu[j] * pmybag->x[j];

	// 	f_new += first + second + third;
	// 	//printf("gradient: %g\n", pmybag->gradient[j]);

	// }
	// f_incr = pmybag->fmin - f_new;
	// printf("f_incr: %g; f_new: %g \n", f_incr, f_new);
	// pmybag->fmin = f_new;
	// if (f_incr < 0.00001) {
	// 	returncode = 1;
	// }
	return returncode;
}

int feasible(int n, double *x, double *lb, double *ub)
{
	int returncode = 0, j;
	double sum, increase;

	printf(" computing first feasible solution\n");

	sum = 0;
	/* set variables to lower bounds */
	for (j = 0; j < n; j++) {
		x[j] = lb[j];
		sum += lb[j];
	}
	printf("after first pass, sum = %g\n", sum);
	for (j = 0; j < n; j++) {
		increase = 1.0 - sum;
		if (increase > ub[j] - x[j])
			increase = ub[j] - x[j];
		x[j] += increase;
		sum += increase;
		printf("after j = %d, sum = %g\n", j, sum);
		if (sum >= 1.0)
			break;
	}
	/*for (j = 0; j < n; j++) {
	printf("x: %g\n", x[j]);
	}*/
	return returncode;
}

// in each iteration, compute the gradient; allocate memory once; 
