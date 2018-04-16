#include <iostream> 

#define RUNNING 1
#define WAITING 0
#define FINISHED 2 
 

class baggie{
public:
	baggie(int NumShare_in, int T_in, double alpha_in, double pi1_in, double pi2_in, double p1_in,
		double p2_in, double p3_in, int name_in);
	~baggie(){ printf("worker %d says goodbye\n", name); } 
  void setconsolemutex(HANDLE consolemutexinput);
  void setmastermutex(HANDLE consolemutexinput);
  void baggiecomp();
  int NumShare, T;
  double getmeits(void){return iterationsdone;}
  void setstatustofinished(void){status = FINISHED;}
  int getstatus(void){ return status; }
  double optimalvalue; 
  int *salevalue;
  void setheavysectionmutex(HANDLE heavysectioninput){heavysectionmutex = heavysectioninput;}
  void setmaxworkersinheavysection(int maxheavy){
	  maxworkersinheavysection = maxheavy;}
  void setnowinheavyaddress(int *paddress){address_of_nowinheavysection = paddress;}
 private:
  
  double alpha, pi1, pi2, p1, p2, p3;
  int name;
  double iterationsdone;
  int status;
  int maxworkersinheavysection;
  int *address_of_nowinheavysection;  /** this is the address of the integer keeping track of how many workers are busy **/
  HANDLE heavysectionmutex;
  HANDLE consolemutex;
  HANDLE mastermutex;
  double result;
  void letmein(void);
  void seeya(void);
};

using namespace std;

