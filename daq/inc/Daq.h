#ifndef _DAQ_H_
#define _DAQ_H_

//Control defines:
//#define SOR_EOB_MODE
//#define USE_DPP_PSD
//#define USE_INTERRUPTS
//#define USE_DPP_ZLE

#define CAEN_USE_DIGITIZERS

#define INT_N_EVENTS 100
#define N_READS   100
//#define TIME_TO_RUN   5 //Time to run in seconds

#ifdef SOR_EOB_MODE 
#define TIME_TO_RUN   50 //Time to run in seconds
#else
#define TIME_TO_RUN   120  //Time to run in seconds
#endif

// Or enable SOR_EOB_MODE

#define BUFSIZE 1000000000

extern int handle;
extern int inburst;
extern int BreakSignal;

void DAQ_exit();
int DAQ_init();

#endif
