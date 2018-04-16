
 
#include <windows.h> 
#include <process.h>
#include "baggie.h"

double mytimecheck(void);

// implementation file for class baggie

baggie :: baggie(int NumShare_in, int T_in, double alpha_in, double pi1_in, double pi2_in, double p1_in,
				double p2_in, double p3_in, int name_in)
{
	NumShare = NumShare_in; T = T_in; alpha = alpha_in; pi1 = pi1_in; pi2 = pi2_in; p1 = p1_in; p2 = p2_in;
	p2 = p2_in; p3 = p3_in;
	name = name_in;
	status = WAITING;
}

void baggie :: setconsolemutex(HANDLE consolemutexinput)
{
	consolemutex = consolemutexinput;
}
void baggie :: setmastermutex(HANDLE mastermutexinput)
{
	mastermutex = mastermutexinput;
}

void baggie :: letmein(void)
{
	char icangoin;
	int localinheavysection;
	
		icangoin = 0;
		while(icangoin == 0){
			Sleep(1000);
			WaitForSingleObject(heavysectionmutex, INFINITE);
			 
			if( (*address_of_nowinheavysection) < maxworkersinheavysection){
				/** key logic: it checks to see if the number of workers in the heavy section is less than the
				number we want to allow **/
				icangoin = 1;
				++*address_of_nowinheavysection; //** increase the count
				localinheavysection = *address_of_nowinheavysection;  
				// so localinheavysection will have the count of busy workers
			}

			ReleaseMutex(heavysectionmutex);
		}
		WaitForSingleObject(consolemutex, INFINITE);
		cout << "******worker" << name <<": I'm in and there are " << localinheavysection <<" total busy workers\n";
		// we can use localinheavysection without protecting it with a mutex, because it is a local variable to this function, i.e.
		// it is not shared with other mutexes
		ReleaseMutex(consolemutex);
}

void baggie :: seeya(void)
{
	
		WaitForSingleObject(heavysectionmutex, INFINITE);
		--*address_of_nowinheavysection;
		ReleaseMutex(heavysectionmutex); 

}

void baggie :: baggiecomp(void)
{
	int i, j, h, t, k;
	int othercounter, bestsell, cumsale;
	double newprice, candidate, bestone, *shift;
	double *optimal,*sale;

	int retcode = 0;
	othercounter = 0;
	status = RUNNING;

	shift = (double *)calloc(NumShare + 1, sizeof(double));
	optimal = (double *)calloc((NumShare + 1)*T, sizeof(double));
	sale = (double *)calloc((NumShare + 1)*T, sizeof(double));
	salevalue = (int *)calloc(T, sizeof(double));

	for (j = 0; j <= NumShare; j++)
		shift[j] = p1*(1 - alpha*pow((double)j, pi1)) + p2*(1 - alpha*pow((double)j, pi2));

	/** do last stage **/
	for (j = 0; j <= NumShare; j++) {
		bestone = 0;
		bestsell = 0;
		for (h = 0; h <= j; h++) {
			//newprice = 1 - alpha*pow((double)j, pi);
			newprice = shift[h];
			candidate = (1 - p3)*newprice*h;

			if (candidate > bestone) {
				bestone = candidate;
				bestsell = h;
			}

			optimal[(T - 1)*(NumShare + 1) + j] = bestone;
			sale[(T - 1)*(NumShare + 1) + j] = bestsell;
		}

		// V[k,t] stored at optimal[t*(NumShare+1) + k] 
	}

	for (t = T - 2; t >= 0; t--) {

		letmein(); // check to see if we can become busy

		double t1 = mytimecheck();  // see the comments below.  mytimecheck() returns the time of day in milliseconds
									// it is defined in mytimer.cpp

		for (j = 0; j <= NumShare; j++) {

			bestone = 0;
			bestsell = 0;
			/** enumerate possibilities **/
			for (h = 0; h <= j; h++) {
				newprice = shift[h];
				candidate = p3*newprice* optimal[(t + 1)*(NumShare + 1) + j]
					+ (1 - p3)*newprice*(h + optimal[(t + 1)*(NumShare + 1) + j - h]);

				if (candidate > bestone) {
					bestone = candidate;
					bestsell = h;
				}
			}

			optimal[t*(NumShare + 1) + j] = bestone;
			sale[t*(NumShare + 1) + j] = bestsell;

	}

		printf("done with stage t = %d\n", t);

		double t2 = mytimecheck();  // check out to see how this function works, it's in mytimer.cpp
									// mytimecheck simply returns the time of day in milliseconds
		double tdiff;

		tdiff = t2 - t1;  // t1 was set above 

		WaitForSingleObject(consolemutex, INFINITE);
		printf(" >> worker %d:  I have completed heavy loop in time %g\n", name, tdiff);
		ReleaseMutex(consolemutex);

		seeya();

		WaitForSingleObject(consolemutex, INFINITE);
		printf(" >> worker %d:  I am out\n", name);
		ReleaseMutex(consolemutex);
	}

	printf("optimal value for trade sequencing = %g\n", optimal[NumShare]);

	optimalvalue = optimal[NumShare];

	
	salevalue[0] = sale[NumShare];
	cumsale = 0;
	printf("optimal sale for trade sequencing = %g\n", salevalue[0]);
	for (int t = 1; t < T; t++) {
		cumsale += salevalue[t - 1];
		salevalue[t] = sale[t*(NumShare+1)+(NumShare - cumsale +1)];
		printf("optimal sale for trade sequencing = %d\n", salevalue[t]);
	}
	//printf("optimal value for trade sequencing = %g\n", optimalvalue);
	

BACK:
	printf("\nran with code %d\n", retcode);
	
}

