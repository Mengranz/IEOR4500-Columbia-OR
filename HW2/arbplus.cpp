#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include <ctime>

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <windows.h>
#include <fstream>
#include <string>

char does_it_exist(char *filename);
int parser(char *sourcefilename, double *solution, int numsec, int *pstatus);
int output(char *solutionfilename, double *solution, int numsec, int status);
int main(int argc, char* argv[])
{
	int retcode = 0;
	FILE *in = NULL, *out = NULL, *res = NULL, *outputFile = NULL;
	char mybuffer[100];
	int numsec, numscen, j, k, numnonz, solutionstatus, *score, numsimu, *allscore, i, *freq;
	double r;
	double *p, optimalvalue, xvalue, *portfolio, *v, *pv, num,x;
	FILE *results = NULL;

	if (argc != 3) {
		printf("Usage:  arb1.exe datafilename lpfilename\n"); retcode = 100; goto BACK;
	}

	in = fopen(argv[1], "r");
	if (in == NULL) {
		printf("could not open %s for reading\n", argv[1]);
		retcode = 200; goto BACK;
	}

	fscanf(in, "%s", mybuffer);
	fscanf(in, "%s", mybuffer);
	numsec = atoi(mybuffer);
	fscanf(in, "%s", mybuffer);
	fscanf(in, "%s", mybuffer);
	numscen = atoi(mybuffer);
	fscanf(in, "%s", mybuffer);
	fscanf(in, "%s", mybuffer);
	r = atof(mybuffer);

	printf("securities: %d, scenarios: %d;;  r = %g\n",
		numsec, numscen, r);

	p = (double *)calloc((1 + numscen)*(1 + numsec), sizeof(double));
	if (p == NULL) {
		printf("no memory\n"); retcode = 400; goto BACK;
	}
	for (k = 0; k <= numscen; k++) {
		fscanf(in, "%s", mybuffer);
		p[k*(1 + numsec)] = 1 + r*(k != 0);
		for (j = 1; j <= numsec; j++) {
			fscanf(in, "%s", mybuffer);
			p[k*(1 + numsec) + j] = atof(mybuffer);
		}
	}

	fscanf(in, "%s", mybuffer);

	fclose(in);

	out = fopen(argv[2], "w");
	if (out == NULL) {
		printf("can't open %s\n", argv[2]); retcode = 500; goto BACK;
	}
	printf("printing LP to file %s\n", argv[2]);

	fprintf(out, "Minimize ");
	for (j = 0; j <= numsec; j++) {
		if (p[j] >= 0) fprintf(out, "+ "); fprintf(out, "%g x%d ", p[j], j);
	}
	fprintf(out, "\n");
	fprintf(out, "Subject to\n");

	for (k = 1; k <= numscen; k++) {
		fprintf(out, "scen%d: ", k);

		for (j = 0; j <= numsec; j++) {
			if (p[k*(1 + numsec) + j] >= 0) fprintf(out, "+ ");
			fprintf(out, "%g x%d ", p[k*(1 + numsec) + j], j);
		}
		fprintf(out, " >= 0\n");
	}

	fprintf(out, "Bounds\n");
	for (j = 0; j <= numsec; j++) {
		fprintf(out, "-1 <= x%d <= 1\n", j);
	}
	fprintf(out, "End\n");

	fclose(out);

	//free(p);

	out = fopen("hidden.dat", "w");
	fclose(out);

	if (does_it_exist("script.py") == 0) {
		printf("need 'script.py'\n"); retcode = 1; goto BACK;
	}

	sprintf(mybuffer, "python script.py %s hidden.dat nothidden.dat", argv[2]);

	printf("mybuffer: %s\n", mybuffer);

	if (does_it_exist("nothidden.dat")) {
		remove("nothidden.dat");
	}

	system(mybuffer);

	/** sleep-wake cycle **/

	for (;;) {
		if (does_it_exist("nothidden.dat")) {
			printf("\ngurobi done!\n");
			Sleep(1000);
			break;
		}
		else {
			Sleep(100);
		}
	}

	portfolio = (double *)calloc(1 + numsec, sizeof(double));
	if (!portfolio) {
		printf("cannot allocate portfolio solution\n");
		retcode = 200; goto BACK;
	}


	/** next, read mygurobi.log **/
	solutionstatus = 700;
	retcode = parser("mygurobi.log", portfolio, numsec, &solutionstatus);
	if (retcode)
		goto BACK;
	retcode = output("solution.dat", portfolio, numsec, solutionstatus); 



    /* read portfolio weight from solution.dat */
	res = fopen("solution.dat", "r");
	v = (double *)calloc(numsec + 1, sizeof(double));

	/* escape the first two words: status 0 */
	fscanf(res, "%s", mybuffer);
	fscanf(res, "%s", mybuffer);

	/* store weights into v */
	for (j = 0; j < numsec + 1; j++) {
		fscanf(res, "%s", mybuffer);
		v[j] = atof(mybuffer);
	}
	
	/* set number of trials to 1,000,000 */
	numsimu = 1000000;

	/* initialize the portfolio value, score and allscore array */
	pv = (double *)calloc(numscen, sizeof(double));
	score = (int*)calloc(numscen, sizeof(double));
	allscore = (int*)calloc(numsimu, sizeof(double));

	

	/* main part */
	srand(static_cast <unsigned> (time(NULL)));
	for (i = 0; i < numsimu; i++) {		
		for (k = 1; k <= numscen; k++) {
			x = p[k*(1 + numsec)];
			pv[k - 1] = x * v[0];
			for (j = 1; j <= numsec; j++) {
				/* generate random float between -1 and 1 */
				num = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 2))-1;			
				//pnew[(k - 1)*(1 + numsec) + j] = x + num*0.05*x;

				/* compute portfolio value */
				x = p[k*(1 + numsec) + j];
				pv[k - 1] += (x + num*0.05*x) * v[j];
			}
			/* assign score 1 to portfolio with value greater than 0; assign 0 otherwise */
			if (pv[k - 1] > 0) {
				score[k - 1] = 1;
			}
			else {
				score[k - 1] = 0;
			}
		}

		/* compute the sum of scores of all scenarios in each simulation */
		int sum = 0;
		for (int a = 0; a < numscen; a++) {
			sum += score[a];
		}
		allscore[i] = sum;
	}
	

	/* compute the allscore frequency */

	freq = (int*)calloc(numsimu, sizeof(double));

	/* initilize every element to 0 */
	for (i = 0; i < numsimu; i++) {
		freq[i] = 0;
	}

	for (i = 0; i < numsimu; i++) {
		++freq[allscore[i]];
	}
    
	/* print number - freq */
	//using namespace std;
	//ofstream outputFile;
	//ofstream fs;
	//std::string filename = "HW2_smalldata.csv";
	//fs.open(filename);
	//outputFile << "Number of Scenarios > 0" << "," << "Frequency" << endl;
	outputFile = fopen("HW2_bigdata.csv", "w");
	fprintf(outputFile, "Number of Scenarios,Frequency\n");
	for (k = 0; k < numscen; k++) {
		if (freq[k] == 0) {
			; 
		}
		else {
			printf("%d : %d\n", k, freq[k]); retcode = 800;
			fprintf(outputFile, "%d,%d\n", k, freq[k]);
			//fs << k << "," << freq[k] << endl;
			/*for (int l = 0; l < freq[k]; l++) std::cout << '*';
				std::cout << std::endl;*/
			
		}
	}

BACK:
	return retcode;
	

}



char does_it_exist(char *filename)
{
	struct stat buf;

	// the function stat returns 0 if the file exists

	if (0 == stat(filename, &buf)){
		return 1;
	}
	else return 0;
}