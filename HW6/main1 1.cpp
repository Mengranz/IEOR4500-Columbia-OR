#include <windows.h> 
#include <process.h>
#include <stdio.h>
#include <stdlib.h> 
#include "baggie.h" 

/****
 This program obtains two parameters from the command-line:  N and W.  The program 
 then launches N worker threads.  Each worker will then attempt to run the computationally heavy
 section (two inner loops) of the "comp" function.   However, at most W threads are allowed
 to do so simultaneously.  This condition will be regulated with a mutex.  Let's use the terminology
 "busy" to describe a worker thread that is employed in the computationally heavy section.  Busy threads
 will either be done with the computation, or be told by the master to terminate according some iteration
 limit that is checked by the master.  As each "busy" worker terminates, any of the remaining "unbusy" 
 workers try to become "busy".  One of them will succeed.  Eventually all workers terminate and the program
 stops.
***/

unsigned _stdcall comp_wrapper(void *foo);



//void comp(double, double, double, double *); 

int main(int argc, char *argv[])
{
	HANDLE *pThread;
	unsigned *pthreadID;
	HANDLE consolemutex;
	HANDLE *mastermutexes;
	int retcode = 0;
	int N, W;
	int *NumShare_arr, *T_arr;
	double *alpha_arr, *pi1_arr, *pi2_arr, *p1_arr, *p2_arr, *p3_arr;	

	baggie **ppbaggies;

	if(argc != 4){
		printf("usage: heavy.exe N W filename\n");
		retcode = 1; goto BACK;
	}

	N = atoi(argv[1]);
	W = atoi(argv[2]);

	if( (N <= 0) || (W <= 0)){
		printf("bad value of N or W: %d %d\n", N, W);
		retcode = 1; goto BACK;
	}

	/** in this program we will have N worker threads simultaneously running; however at any time at most
	  W of them will be allowed to actually be computing in a "heavy" way.  This simulates a multi-agent
	  search where different search algorithms take turns.   Controlling W allows us to control how CPU time this process
	  effectively is given **/


	ppbaggies = (baggie **)calloc(N, sizeof(baggie *));
	/** ppbaggies is an array, each of whose members is the address of a baggie, and so the type of ppbaggies is baggie ** **/
	if(ppbaggies == NULL){
		cout << "cannot allocate" << N <<"baggies\n";
		retcode = 1; goto BACK;
	}
	pThread = (HANDLE *)calloc(N, sizeof(HANDLE));
	pthreadID= (unsigned *)calloc(N, sizeof(unsigned));
	mastermutexes = (HANDLE *)calloc(N, sizeof(HANDLE));
	if((pThread == NULL) || (pthreadID == NULL) || (mastermutexes == NULL)){
		cout << "cannot allocate" << N << "handles and threadids\n";
		retcode = 1; goto BACK;
	}

	FILE *in;
	char mybuffer[200];

	in = fopen(argv[3], "r"); /*datafile*/

	if (in == NULL) {
		printf("could not open %s for reading\n", argv[3]);
		retcode = 200; goto BACK;
	}

	/*read parameters from the data file*/
	NumShare_arr = (int *)calloc(N, sizeof(int));
	T_arr = (int *)calloc(N, sizeof(int));
	alpha_arr = (double *)calloc(N, sizeof(double));
	pi1_arr = (double *)calloc(N, sizeof(double));
	pi2_arr = (double *)calloc(N, sizeof(double));
	p1_arr = (double *)calloc(N, sizeof(double));
	p2_arr = (double *)calloc(N, sizeof(double));
	p3_arr = (double *)calloc(N, sizeof(double));

	//skip to read first line
	for (int k = 0; k < 8; k++) {
		fscanf(in, "%s", mybuffer);
	}

	for (int k = 0; k < N; k++) {

		fscanf(in, "%s", mybuffer);
		NumShare_arr[k] = atoi(mybuffer);

		fscanf(in, "%s", mybuffer);
		T_arr[k] = atoi(mybuffer);

		fscanf(in, "%s", mybuffer);
		alpha_arr[k] = atof(mybuffer);

		fscanf(in, "%s", mybuffer);
		pi1_arr[k] = atof(mybuffer); 

		fscanf(in, "%s", mybuffer);
		pi2_arr[k] = atof(mybuffer);

		fscanf(in, "%s", mybuffer);
		p1_arr[k] = atof(mybuffer);

		fscanf(in, "%s", mybuffer);
		p2_arr[k] = atof(mybuffer);

		fscanf(in, "%s", mybuffer);
		p3_arr[k] = atof(mybuffer);

	}

	printf("Finish reading data file!\n");
	fclose(in);

	for(int j = 0; j < N; j++){
		ppbaggies[j] = new baggie(NumShare_arr[j], T_arr[j],alpha_arr[j],pi1_arr[j],
									pi2_arr[j],p1_arr[j],p2_arr[j],p3_arr[j],j);  
	}


	consolemutex = CreateMutex(NULL, 0, NULL);

	for(int j = 0; j < N; j++){
		ppbaggies[j]->setconsolemutex(consolemutex); // consolemutex shared across workers plus master
	}

	HANDLE heavymutex;
	heavymutex = CreateMutex(NULL, 0, NULL);

	int nowinheavy = 0;

	for(int j = 0; j < N; j++){
		mastermutexes[j] = CreateMutex(NULL, 0, NULL);
		ppbaggies[j]->setmastermutex(mastermutexes[j]);

		ppbaggies[j]->setmaxworkersinheavysection(W);
		ppbaggies[j]->setheavysectionmutex(heavymutex); 

		ppbaggies[j]->setnowinheavyaddress( &nowinheavy );


	}

	for(int j = 0; j < N; j++){
		pThread[j] = (HANDLE)_beginthreadex( NULL, 0, &comp_wrapper, (void *) ppbaggies[j], 
			0, 		&pthreadID[j] );
	}
	
	

	int numberrunning = N;

	for (; numberrunning > 0;) {
		Sleep(5000);
		printf("master will now check on workers\n"); fflush(stdout);

		 for(int j = 0; j < N; j++){
			//double jiterations;
			char jstatus = RUNNING;

			WaitForSingleObject(mastermutexes[j], INFINITE);

			jstatus = ppbaggies[j]->getstatus();
			//fprintf(outputFile, "%g\n", ppbaggies[j]->optimalvalue);
			ReleaseMutex(mastermutexes[j]);

			if(jstatus == RUNNING){
		
				WaitForSingleObject(mastermutexes[j], INFINITE);

				//jiterations = ppbaggies[j]->getmeits();

				--numberrunning;

				ReleaseMutex(mastermutexes[j]);
			
				WaitForSingleObject(consolemutex, INFINITE);
				//printf("master: worker %d has done %g iterations\n", j, 
					//jiterations);
				ReleaseMutex(consolemutex);
			}
		}
	}

	
	

	for(int j = 0; j < N; j++){
		WaitForSingleObject(pThread[j], INFINITE);
		printf("--> thread %d done\n", j); 
		//delete ppbaggies[j]; // calls destructor
	}
	
	FILE *outputFile;
	outputFile = fopen("optimal.csv", "w");

	for (int i = 0; i < N; i++) {
		printf("opt: %g\n", ppbaggies[i]->optimalvalue);
		fprintf(outputFile, "%g,", ppbaggies[i]->optimalvalue);
		for (int t = 0; t < ppbaggies[i]->T; t++) {
			fprintf(outputFile, "%d,", ppbaggies[i]->salevalue[t]);
		}
		fprintf(outputFile, "\n");

		delete ppbaggies[i];
	}

	fclose(outputFile);
	
	free(ppbaggies);
BACK:
	return retcode;
}



unsigned _stdcall comp_wrapper(void *genericaddress)
{
	baggie *pbaggie = (baggie *) genericaddress;

	pbaggie->baggiecomp();

	return 0;
}

