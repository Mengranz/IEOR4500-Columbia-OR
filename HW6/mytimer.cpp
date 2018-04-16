#include<sys\timeb.h>

double mytimecheck(void)
{
   double seconds, millis;
   struct timeb mytimeb;

   ftime(&mytimeb);

   seconds = (double) mytimeb.time;
   millis = ( (double) mytimeb.millitm)/1000.00;

   return seconds+millis;
}


 
