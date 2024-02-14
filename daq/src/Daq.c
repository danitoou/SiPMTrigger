#include <stdlib.h>
#include <stdio.h>
#include "CAENDigitizer.h"
#include <time.h>
#include <sys/stat.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/time.h>
#include <signal.h>

#include "Daq.h"
#include "FADC.h"
#include "config.h"

int handle;
int inburst;
int BreakSignal = 0;

config_t *cfg = NULL;
CAEN_DGTZ_DPP_PSD_Params_t DPPParams;
DigitizerParams_t Params;
char *mybuf = NULL;
int *mybsize = NULL;
FILE *logfile;


void
termination_handler (int signum)
{

  if (inburst > 0)
    {
      //Cannot send signal to the board ... it's busy
      printf ("Still in burst mode ... \n");
      //sleep(5);
      BreakSignal = 1;
      signal (SIGINT, SIG_IGN);
      printf ("Received signal: %d  \n", signum);
    }
  else
    {
      signal (SIGINT, SIG_DFL);
      fflush (NULL);
      printf ("Received signal: %d \n", signum);
      printf ("Out of burst mode ... \n");

      DAQ_exit ();
    }
  /* //  CAEN_DGTZ_Reset (handle); */
  /* CAEN_DGTZ_Reset (handle); */
  /* CAEN_DGTZ_CloseDigitizer(handle); */


  /* //  DAQ_exit(); */
}



int
acq_run (int handle)
{
  int reg = 0x8104;
  int data;
  CAEN_DGTZ_ReadRegister (handle, reg, &data);
  return (data >> 2) & 0x1;
}




int
DAQ_init ()
{
  default_config ();
  if (cfg->usedpp)
    {
      //Configure the DPP
      default_config_dpp ();
    }

  BreakSignal = 0;
  inburst = 0;
  return 0;
}

void
DAQ_exit ()
{
  int ret = -1;
  //Relies only on the global variables
  //Should be called after closing all the files

  //free the alocated memories:
  //Close the connections to the board:
  //while (ret!=0) {

  printf ("Enterring DAQ EXIT\n");
  CAEN_DGTZ_Reset (handle);
  /* CAEN_DGTZ_CloseDigitizer(handle); */

  /* exit(0); */
  ret |= FADCStopAcquisition ();	//stop the acquisition, if running
  ret |= CAEN_DGTZ_Reset (handle);
  //}
  printf ("Board reset \n");

  fflush (NULL);

  if (mybuf)
    free (mybuf);
  if (mybsize)
    free (mybsize);
  if (cfg)
    free (cfg);

  printf ("Run stopped \n");
  //write all buffers:
  //Close the connection
  CAEN_DGTZ_CloseDigitizer (handle);
  printf ("Digitizer closed \n");

  exit (0);
}



double
time_diff (struct timeval *t1, struct timeval *t2)
{
  //  time_difference in seconds:
  return (((double) t2->tv_sec - (double) t1->tv_sec) * 1000000. +
	  ((double) t2->tv_usec - (double) t1->tv_usec)) / 1000000.;

}


int
main (int argc, char **argv)
{

  CAEN_DGTZ_ErrorCode ret;
  CAEN_DGTZ_EventInfo_t eventInfo;
  CAEN_DGTZ_DPP_PSD_Event_t *Events[8];
  CAEN_DGTZ_DPP_PSD_Waveforms_t *Waveform = NULL;


  void *Evt = NULL;
  char *buffer = NULL;

  int burst = 0;

  int i, b, k;
  int c = 0, count = 0;
  char *evtptr = NULL;
  uint32_t size, bsize;
  uint32_t numEvents;

  uint32_t NumEvents[8];

  uint32_t data;
  uint32_t reg;

  int filehandle;
  char ofname[256];

  uint32_t mybufsize;
  int pos = 0;

  double tot_size = 0.;


  time_t tim;
  struct tm *tt;
  time_t tbeg;
  time_t tend;
  struct tm *tt_beg;

  struct timeval tv_beg;
  struct timeval tv_end;
  double dt;

  int iprev_ev = -1;
  int nmissed_ev = 0;

  int bsize_limit;

  char logfilename[256];
  int nsofttrig = 0;


  //Signal processing
  signal (SIGINT, termination_handler);
  printf ("Trap signals: %d, ", SIGINT);
  signal (SIGHUP, termination_handler);
  printf ("%d, ", SIGHUP);
  signal (SIGTERM, termination_handler);
  printf ("%d, ", SIGTERM);
  signal (SIGKILL, termination_handler);
  printf ("%d, ", SIGKILL);
  signal (SIGUSR2, termination_handler);
  printf ("%d, ", SIGUSR2);
  signal (SIGFPE, termination_handler);
  printf ("%d\n", SIGFPE);

  memset (&Params, 0, sizeof (DigitizerParams_t));
  memset (&DPPParams, 0, sizeof (CAEN_DGTZ_DPP_PSD_Params_t));


  DAQ_init ();

  //Board opening and basic reset
  FADCBoardInit ();

  FADCPostConfig ();

  //Channels setup
  FADCChannelConfig ();

  //Record length and information
  FADCDigitizationSetup ();

  //Trigger setup and input signals levels
  FADCTriggerConfig ();

  //Transfer setup
  ret = FADCIOSetup ();

  FADCPostConfig ();

  if (cfg->usedpp)
    {
      FADCDPPConfig ();
    }



  if (ret != CAEN_DGTZ_Success)
    {
      printf ("Errors during Digitizer Configuration.\n");
      DAQ_exit ();
      return 0;
    }

  /*

     //Initialize the readout buffer - automatic size calculation:
     ret |= CAEN_DGTZ_MallocReadoutBuffer(handle, &buffer, &size);
     printf("Allocated data buffer with size: %d\n",size);
   */


  mybuf = malloc (BUFSIZE);	//~1200MB buffer
  if (mybuf == NULL)
    {
      return 0;
    }
  printf ("Allocated memory buffer with size %dMB\n", BUFSIZE / 1.e6);
  bsize_limit = (int) (0.8 * BUFSIZE);


  //  mybsize = malloc(1000000*sizeof(int)); //array for the buf sizes
  mybsize = malloc (100000 * sizeof (int));	//array for the buf sizes
  if (mybsize == NULL)
    {
      return 0;
    }
  printf ("Allocated buffer for transfer sizes storage\n");


  mybufsize = 0;
  tot_size = 0.;
  count = 0;
  numEvents = 0;
  tbeg = tim = time (0);
  tt = localtime (&tim);
  printf ("TIME START %d %.2d %.2d %.2d %d %d \n", tt->tm_year + 1900,
	  tt->tm_mon + 1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
  gettimeofday (&tv_beg, NULL);



  if (cfg->usedpp)
    {
      uint32_t AllocatedSize, BufferSize;
      /* Allocate memory for the readout buffer */
      ret =
	CAEN_DGTZ_MallocReadoutBuffer (handle, (char **) (&buffer),
				       &AllocatedSize);
      /* Allocate memory for the events */
      ret |=
	CAEN_DGTZ_MallocDPPEvents (handle, (void **) Events, &AllocatedSize);
      /* Allocate memory for the waveforms */
      ret |=
	CAEN_DGTZ_MallocDPPWaveforms (handle, (void **) &Waveform,
				      &AllocatedSize);
      if (ret)
	{
	  printf ("Can't allocate memory buffers\n");
	  DAQ_exit ();
	  return 0;
	}

      //CAEN_DGTZ_SWStartAcquisition(handle);
      printf ("DPP Acquisition Started\n");

      int nread = 0;
      mybufsize = 0;
      FADCStartAcquisition ();

      gettimeofday (&tv_beg, NULL);
      tim = time (0);
      tt_beg = localtime (&tim);
      numEvents = 0;
      while (nread < 1000)
	{
	  //  CAEN_DGTZ_SendSWtrigger(handle);
	  CAEN_DGTZ_ReadData (handle, CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,
			      &mybuf[mybufsize], &BufferSize);

	  if (BufferSize > 0)
	    {
	      ret |=
		CAEN_DGTZ_GetDPPEvents (handle, &mybuf[mybufsize], BufferSize,
					(void **) Events, NumEvents);
	      /*
	         printf("Read out buffer with size %d\n",BufferSize);
	       */
	      int j;
	      for (j = 0; j < 8; j++)
		{
		  printf ("Channel %d: events: %d\n", j, NumEvents[j]);
		}

	      numEvents += NumEvents[0];
	    }

	  nread++;
	  mybufsize += BufferSize;
	}
      CAEN_DGTZ_SWStopAcquisition (handle);
      gettimeofday (&tv_end, NULL);

      printf ("DPP Acquisition Stopped\n");

      dt = time_diff (&tv_beg, &tv_end);
      tot_size = mybufsize / (1024. * 1024.);
      printf ("Burst END    %d : %d \n", tt->tm_min, tt->tm_sec);
      printf ("Total number of events: %d, Total size read: %lf MB\n",
	      numEvents, tot_size);
      if (dt > 0)
	{
	  printf ("Time: %d s, time: %f s, Transfer rate: %7.2f MB/s\n",
		  (tend - tbeg), dt, tot_size / dt);
	  printf ("Total rate: %7.2f ev/s\n", numEvents / dt);
	}

      //Write data to file:
      //   tend =tim = time(0);
      //tt_beg  = localtime(&tim);   
      sprintf (ofname, "data/dpprawdata_%d_%.2d_%.2d_%.2d_%.2d_%.2d",
	       tt_beg->tm_year + 1900, tt_beg->tm_mon + 1, tt_beg->tm_mday,
	       tt_beg->tm_hour, tt_beg->tm_min, tt_beg->tm_sec);
      printf ("Writing %2.10f s to file %s \n", dt, ofname);
      filehandle = open (ofname, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);
      write (filehandle, mybuf, mybufsize);
      close (filehandle);






      DAQ_exit ();
      return 0;
    }





#ifdef SOR_EOB_MODE
  //Start of the acquisition
  FADCStartAcquisition ();
  inburst = 1;
  burst = 0;

  //Constant reading
  //  while(burst <= 10)  {
  gettimeofday (&tv_end, NULL);
  while (1)
    {

      if (daq_quit ())
	{
	  printf ("Taken %d number of bursts, EXITING\n", burst);
	  DAQ_exit ();
	}

      gettimeofday (&tv_beg, NULL);
      tim = time (0);
      tt_beg = localtime (&tim);
      while (acq_run (handle) && (time_diff (&tv_beg, &tv_end) < TIME_TO_RUN)
	     && mybufsize < bsize_limit && nsofttrig < 10)
	{
	  if (cfg->exttrigger == 0 && cfg->chtrigger == 0)
	    {
	      CAEN_DGTZ_SendSWtrigger (handle);
	      usleep (10000);
	      nsofttrig++;
	      //Send a software trigger
	    }
	  inburst = 2;
	  //Continuous reading of the data
	  CAEN_DGTZ_ReadData (handle, CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,
			      &mybuf[mybufsize], &bsize);
	  mybufsize += bsize;
	  gettimeofday (&tv_end, NULL);
	  //      inburst = 2;
	}

      //sleep(1);

      //Postprocessing the burst
      if (inburst == 2)
	{
	  //Stop the acquisition until we process the output data
	  FADCStopAcquisition ();
	  inburst = 0;

	  //inburst = 0;
	  gettimeofday (&tv_end, NULL);
	  //Doing a final read in case there are still some events stored in the FADC buffers

	  /*
	     do {
	     CAEN_DGTZ_ReadData(handle,CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,&mybuf[mybufsize],&bsize); 
	     mybufsize += bsize;
	     } while (bsize > 0);
	     if(daq_quit()){
	     CAEN_DGTZ_Reset(handle); 
	     }
	     //HACK - forget about the rest of the data, anyway we are resetting 
	   */

	  //Get timing
	  printf ("Postprocessing:\n");

	  dt = time_diff (&tv_beg, &tv_end);
	  tend = tim = time (0);
	  tt = localtime (&tim);

	  CAEN_DGTZ_GetNumEvents (handle, mybuf, mybufsize, &numEvents);
	  //printf("Number of events in buffer: %d\n",numEvents);

	  tot_size = mybufsize / (1024. * 1024.);
	  printf ("Burst END:   %d : %d : %d \n", tt->tm_hour, tt->tm_min,
		  tt->tm_sec);
	  printf ("Total number of events: %d, Total size read: %lf MB\n",
		  numEvents, tot_size);
	  if (dt > 0)
	    {
	      printf
		("Time working: %d s, time for burst: %f s, Transfer rate: %7.2f MB/s\n",
		 (tend - tbeg), dt, tot_size / dt);
	      printf ("Total rate: %7.2f ev/s\n", numEvents / dt);
	    }

	  //Write data to file:
	  sprintf (ofname, "data/rawdata_%d_%.2d_%.2d_%.2d_%.2d_%.2d",
		   tt_beg->tm_year + 1900, tt_beg->tm_mon + 1,
		   tt_beg->tm_mday, tt_beg->tm_hour, tt_beg->tm_min,
		   tt_beg->tm_sec);
	  printf ("Writing %2.10f s to file %s \n", dt, ofname);
	  filehandle = open (ofname, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);
	  write (filehandle, mybuf, mybufsize);
	  close (filehandle);

	  sprintf (ofname,
		   "status/daq-ready/rawdata_%d_%.2d_%.2d_%.2d_%.2d_%.2d",
		   tt_beg->tm_year + 1900, tt_beg->tm_mon + 1,
		   tt_beg->tm_mday, tt_beg->tm_hour, tt_beg->tm_min,
		   tt_beg->tm_sec);
	  filehandle = open (ofname, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);
	  close (filehandle);



	  //Start the buffer from the beginning
	  mybufsize = 0;
	  inburst = 0;
	  nsofttrig = 0;

	  burst++;
	  if (BreakSignal == 1 || daq_quit ())
	    {

	      printf ("Bursts on disk: %d; EXITING\n", burst);
	      DAQ_exit ();
	    }
	  else
	    {
	      //Clean all error buffers and start the acquisition again
	      /* FADCCleanErrors(); */
	      /* FADCCleanBuffers(); */
	      inburst = 0;
	      //ret = CAEN_DGTZ_Reset (handle);
	      if (ret)
		{
		  printf ("Cannot reset the board, exiting: %d\n", ret);
		  DAQ_exit ();
		  break;
		}
	      //CAEN_DGTZ_CloseDigitizer(handle);

	      //FADCBoardInit();

	      /* FADCChannelConfig (); */

	      /* FADCDigitizationSetup (); */

	      //FADCTriggerConfig ();

	      /* ret = FADCIOSetup (); */

	      /* FADCPostConfig (); */

	      FADCStartAcquisition ();
	      inburst = 1;
	    }

	}

      if (BreakSignal == 1)
	{			// check signals caught
	  printf ("\n");
	  break;
	}

    }

  DAQ_exit ();



#else
  //Define the burst as a fixed amount of time
  burst = 0;
  //  if(cfg->exttrigger == 0 && cfg->chtrigger == 0) {
  nsofttrig = 0;
  // }
  while (1)
    {

      if (daq_quit ())
	{
	  printf ("Taken %d number of bursts, EXITING\n", burst);
	  DAQ_exit ();
	}

      //===============================================================
      // Start a new burst
      FADCStartAcquisition ();
      inburst = 1;

      gettimeofday (&tv_beg, NULL);
      tim = time (0);
      tt_beg = localtime (&tim);

      do
	{
	  if (0 && cfg->exttrigger == 0 && cfg->chtrigger == 0)
	    {
	      //Send a software trigger
	      CAEN_DGTZ_SendSWtrigger (handle);
	      usleep (10000);	//Sleep for 10ms (RO window is ~2ms)
	      nsofttrig++;
	    }
	  //for( k=0;k<N_READS;k++) {
	  //Continuous reading of the data
	  CAEN_DGTZ_ReadData (handle, CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,
			      &mybuf[mybufsize], &bsize);
	  mybufsize += bsize;
	  gettimeofday (&tv_end, NULL);
	}
      while (time_diff (&tv_beg, &tv_end) < TIME_TO_RUN
	     && mybufsize < bsize_limit && nsofttrig < 100);

      FADCStopAcquisition ();
      // Finish the burst
      //===============================================================
      inburst = 0;

      //      sleep(1);

      //  gettimeofday(&tv_end,NULL); //Update the time
      dt = time_diff (&tv_beg, &tv_end);

      CAEN_DGTZ_GetNumEvents (handle, mybuf, mybufsize, &numEvents);
      //printf("Number of events in buffer: %d\n",numEvents);
      tend = tim = time (0);
      tt = localtime (&tim);


      tot_size = mybufsize / (1024. * 1024.);
      printf ("Burst END    %d : %d \n", tt->tm_min, tt->tm_sec);
      printf ("Total number of events: %d, Total size read: %lf MB\n",
	      numEvents, tot_size);
      if (dt > 0)
	{
	  printf ("Time: %d s, time: %f s, Transfer rate: %7.2f MB/s\n",
		  (tend - tbeg), dt, tot_size / dt);
	  printf ("Total rate: %7.2f ev/s\n", numEvents / dt);
	}

      //Write data to file:
      sprintf (ofname, "data/rawdata_%d_%.2d_%.2d_%.2d_%.2d_%.2d",
	       tt_beg->tm_year + 1900, tt_beg->tm_mon + 1, tt_beg->tm_mday,
	       tt_beg->tm_hour, tt_beg->tm_min, tt_beg->tm_sec);
      printf ("Writing %2.10f s to file %s \n", dt, ofname);
      filehandle = open (ofname, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);
      write (filehandle, mybuf, mybufsize);
      close (filehandle);

      sprintf (ofname, "status/daq-ready/rawdata_%d_%.2d_%.2d_%.2d_%.2d_%.2d",
	       tt_beg->tm_year + 1900, tt_beg->tm_mon + 1, tt_beg->tm_mday,
	       tt_beg->tm_hour, tt_beg->tm_min, tt_beg->tm_sec);
      filehandle = open (ofname, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);
      close (filehandle);

      if (BreakSignal == 1)
	{			// check signals caught
	  printf ("\n");
	  break;
	}

      //Start the buffer from the beginning
      CAEN_DGTZ_Reset (handle);
      //FADCBoardInit();

      FADCChannelConfig ();

      FADCDigitizationSetup ();

      FADCTriggerConfig ();

      ret = FADCIOSetup ();

      FADCPostConfig ();


      FADCCleanErrors();
      FADCCleanBuffers(); 
      mybufsize = 0;
      burst++;
      if (cfg->exttrigger == 0 && cfg->chtrigger == 0)
	{
	  nsofttrig = 0;
	}

      if (BreakSignal == 1)
	{			// check signals caught
	  printf ("\n");
	  break;
	}


    }

  DAQ_exit ();


#endif



  //Polling based reading
#ifndef USE_INTERRUPTS
  for (k = 0; k < N_READS; k++)
    {
      //ret |= CAEN_DGTZ_SendSWtrigger(handle); /* Send a SW Trigger */

      //    ret |= CAEN_DGTZ_ReadData(handle,CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,buffer,&bsize); 
      ret |=
	CAEN_DGTZ_ReadData (handle, CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,
			    &mybuf[mybufsize], &mybsize[k]);
      mybufsize += mybsize[k];

      continue;

      //printf("Buffer size: %d\n",bsize);

      ret |=
	CAEN_DGTZ_GetNumEvents (handle, &mybuf[mybufsize], bsize, &numEvents);
      //    ret |= CAEN_DGTZ_GetNumEvents(handle,buffer,bsize,&numEvents);
      tot_size += (bsize) / (1024. * 1024.);


      if (numEvents >= 1)
	{
	  // printf("Read %d events\n",numEvents);
	  count += numEvents;

	  if (0)
	    for (i = 0; i < numEvents; i++)
	      {
		/* Get the Infos and pointer to the event */
		ret |=
		  CAEN_DGTZ_GetEventInfo (handle, buffer, bsize, i,
					  &eventInfo, &evtptr);

		/* Decode the event to get the data */
		ret |= CAEN_DGTZ_DecodeEvent (handle, evtptr, &Evt);
		//*************************************
		// Event Elaboration
		//*************************************
		ret |= CAEN_DGTZ_FreeEvent (handle, &Evt);
	      }

	}
    }


#else


  for (k = 0; k < N_READS; k++)
    {
      ret = CAEN_DGTZ_IRQWait (handle, 10000);	//Wait for an interruption request, timeout 100s
      //printf("Interrupt recieved: %d\n",ret);
      ret |=
	CAEN_DGTZ_ReadData (handle, CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,
			    &mybuf[mybufsize], &mybsize[k]);
      mybufsize += mybsize[k];
      ret =
	CAEN_DGTZ_SetInterruptConfig (handle, CAEN_DGTZ_DISABLE, 1, 0,
				      INT_N_EVENTS, CAEN_DGTZ_IRQ_MODE_RORA);
      // printf("Interrupt acknoledged: %d\n",ret);
      CAEN_DGTZ_SetInterruptConfig (handle, CAEN_DGTZ_ENABLE, 1, 0,
				    INT_N_EVENTS, CAEN_DGTZ_IRQ_MODE_RORA);

      //    ret = CAEN_DGTZ_IACKCycle( handle,&board_id);


      //CAEN_DGTZ_RearmInterrupt (handle);
    }


#endif

  CAEN_DGTZ_WriteRegister (handle, reg, 0);

  //CAEN_DGTZ_SWStopAcquisition(handle);
  gettimeofday (&tv_end, NULL);


  dt = time_diff (&tv_beg, &tv_end);

  tend = tim = time (0);
  tt = localtime (&tim);


  //Process time:




  //checks on the data:
  printf ("Number of events transferred\n");
  for (i = 0; i < k; i++)
    {
      CAEN_DGTZ_GetNumEvents (handle, &mybuf[pos], mybsize[i], &numEvents);
      pos += mybsize[i];
      printf ("%d\t", numEvents);
      if (i % 8 == 7)
	printf ("\n");
      count += numEvents;
      //    tot_size += mybsize[i]/(1024.*1024.);
    }
  printf ("\n");

  printf ("mybufsize: %d, last position: %d\n", mybufsize, pos);

  tot_size = mybufsize / (1024. * 1024.);

  printf ("Postprocessing:\n");
  CAEN_DGTZ_GetNumEvents (handle, mybuf, mybufsize, &numEvents);
  printf ("Number of events sum: %d \t Number of events in buffer: %d\n",
	  count, numEvents);


  nmissed_ev = 0;
  iprev_ev = 0;
  if (0)
    for (i = 0; i < numEvents; i++)
      {
	/* Get the Infos and pointer to the event */
	ret |=
	  CAEN_DGTZ_GetEventInfo (handle, mybuf, mybufsize, i, &eventInfo,
				  &evtptr);
	/*
	   if(iprev_ev == -1) {
	   iprev_ev = 1;
	   } else 
	 */
	//printf("EventCounter: %d\n",eventInfo.EventCounter);
	if (eventInfo.EventCounter - iprev_ev != 1)
	  {
	    nmissed_ev++;
	    printf ("Event number mismatch: prev %d, curr %d\n", iprev_ev,
		    eventInfo.EventCounter);
	  }

	iprev_ev = eventInfo.EventCounter;
      }

  ret |=
    CAEN_DGTZ_GetEventInfo (handle, mybuf, mybufsize, numEvents - 1,
			    &eventInfo, &evtptr);

  printf ("numEvents: %d, nmissed_ev: %d, eventInfo.EventCounter: %d\n",
	  numEvents, nmissed_ev, eventInfo.EventCounter);

  nmissed_ev = eventInfo.EventCounter - numEvents;

  if (nmissed_ev > 0)
    {
      printf
	("========= Event number differs - numEvents: %d, nmissed_ev: %d, eventInfo.EventCounter: %d\n",
	 numEvents, nmissed_ev, eventInfo.EventCounter);
    }




  printf ("TIME END    %d : %d \n", tt->tm_min, tt->tm_sec);
  printf ("Total number of events: %d, Total size read: %lf MB\n", count,
	  tot_size);
  if (dt > 0)
    {
      printf ("Time: %d s, time: %f s, Transfer rate: %7.2f MB/s\n",
	      (tend - tbeg), dt, tot_size / dt);
      printf ("Total rate: %7.2f ev/s\n", count / dt);
    }

  //Write data to file:
  sprintf (ofname, "rawdata_%d_%d_%d_%d_%d_%d", tt->tm_year + 1900,
	   tt->tm_mon, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
  printf ("Writing to file %s\n", ofname);
  filehandle = open (ofname, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);
  write (filehandle, mybuf, mybufsize);
  close (filehandle);


  DAQ_exit ();


  return 0;
}
