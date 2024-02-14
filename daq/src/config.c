#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <unistd.h>
#include "config.h"



config_t * default_config(){
  // config_t *cfg;
  int i;

  cfg =  malloc(sizeof(config_t));
  memset(cfg,0,sizeof(config_t));

  strcpy(cfg->cfg_file,"cfg/daq.cfg");
  strcpy(cfg->quit_file,"cfg/quit");
  
  cfg->nchan = 8;//4;
  cfg->chmask = 0xff;//0xf;//0b11;//0xff;
  cfg->chtrigmask = 0xff;// 0xf;// 0b1;//0xff;
  cfg->dual = 0;

  //Boolean variables 
  cfg->chtrigger = 1;
  cfg->exttrigger = 0;

  cfg->trigthr = 20.; //Thrigger threshold in mV
  //cfg->trigthr = 5.; //Thrigger threshold in mV

  //cfg->useburst = 1;
  //cfg->useburst = 1;
  
  //cfg->iolevel = IOTTL;
  cfg->iolevel = IONIM;
  cfg->usedpp = 0;
  
  cfg->nev = 100;
  cfg->tr_overlap = 0;

  for(i=0;i<8;i++) {
    cfg->choffset[i] =  0x0;
  }
  for(i=0;i<8;i++) {
    cfg->chped[i] = 970;
  }

  cfg->chped[0] = 977;
  cfg->chped[1] = 970;
  /* cfg->chped[2] = 1016; */
  /* cfg->chped[3] = 1017; */
  /* cfg->chped[4] = 1016; */
  /* cfg->chped[5] = 1016; */
  /* cfg->chped[6] = 1016; */
  /* cfg->chped[7] = 1016; */
  /* cfg->chped[0] = 966.900024; */
  /* cfg->chped[1] = 963.842102; */
  /* cfg->chped[2] = 965.950012; */
  /* cfg->chped[3] = 729.950012; */
  /* cfg->chped[4] = 990.950012; */
  /* cfg->chped[5] = 984.950012; */
  /* cfg->chped[6] = 981.799988; */
  /* cfg->chped[7] = 974.333313; */

  
  cfg->chped[0] = 978;
  cfg->chped[1] = 975;
  cfg->chped[2] = 978;
  cfg->chped[3] = 742;
  cfg->chped[4] = 1002;
  cfg->chped[5] = 996;
  cfg->chped[6] = 993;
  cfg->chped[7] = 986;

  


  cfg->irq=0;

  cfg->nblt = 1024;

  return cfg;
}


int default_config_dpp(){
  int i;
  
  memset(&Params,0,sizeof(DigitizerParams_t));
  memset(&DPPParams,0,sizeof(CAEN_DGTZ_DPP_PSD_Params_t));

  //  Params.LinkType = CAEN_DGTZ_PCI_OpticalLink;
  Params.LinkType = CAEN_DGTZ_OpticalLink;
  Params.VMEBaseAddress = 0; 
  Params.IOlev = CAEN_DGTZ_IOLevel_NIM;
  Params.ChannelMask = cfg->chmask;
  Params.PulsePolarity = CAEN_DGTZ_PulsePolarityNegative;
  Params.RecordLength=48;
  //Params.AcqMode = CAEN_DGTZ_DPP_ACQ_MODE_Mixed;          // CAEN_DGTZ_DPP_ACQ_MODE_List or CAEN_DGTZ_DPP_ACQ_MODE_Oscilloscope
  Params.AcqMode = CAEN_DGTZ_DPP_ACQ_MODE_List;
  Params.EventAggr = 512;

  for(i=0;i<cfg->nchan;i++) {
    //Setting the defaults in the DPP structure:
    DPPParams.thr[i] = 10;  // Trigger Threshold
    /* The following parameter is used to specifiy the number of samples for the baseline averaging: 
       0 -> absolute Bl
       1 -> 4samp
       2 -> 8samp
       3 -> 16samp
       4 -> 32samp
       5 -> 64samp
       6 -> 128samp */
    DPPParams.nsbl[i] = 2;
    DPPParams.lgate[i] = 10;   // Long Gate Width (N*4ns)
    DPPParams.sgate[i] = 5;    // Short Gate Width (N*4ns)
    DPPParams.pgate[i] = 5;   // Pre Gate Width (N*4ns)
    /* Self Trigger Mode:
       0 -> Disabled
       1 -> Enabled */
    DPPParams.selft[i] = 1;
    /* Trigger configuration:
       CAEN_DGTZ_DPP_TriggerConfig_Peak -> trigger on peak
       CAEN_DGTZ_DPP_TriggerConfig_Threshold -> trigger on threshold */
    DPPParams.trgc[i] = CAEN_DGTZ_DPP_TriggerConfig_Threshold; // Ignored for x751
    /* Trigger Validation Acquisition Window */
    DPPParams.tvaw[i] = 50;
    /* Charge sensibility: 0->40fc/LSB; 1->160fc/LSB; 2->640fc/LSB; 3->2,5pc/LSB */
    DPPParams.csens[i] = 0;
  }
  /* Pile-Up rejection Mode
     CAEN_DGTZ_DPP_PSD_PUR_DetectOnly -> Only Detect Pile-Up
     CAEN_DGTZ_DPP_PSD_PUR_Enabled -> Reject Pile-Up */
  DPPParams.purh = CAEN_DGTZ_DPP_PSD_PUR_DetectOnly; // Ignored for x751
  /* Purity Gap */
  DPPParams.purgap = 100;  // Ignored for x751
  DPPParams.blthr = 3;     // Baseline Threshold
  DPPParams.bltmo = 100;   // Baseline Timeout
  
}




int read_config(config_t *cfg,char *fname){
  FILE *fin;
  if (cfg == NULL) {
    printf("Please read the config after initialization of the default configuration values\n");
    return 0;
  }
  if(fname!=NULL) {
    strcpy(cfg->cfg_file,fname);
  }
  if( (fin = fopen(cfg->cfg_file,"r")) == NULL) {
    printf("Cannot open configuration file %s for reading\n",cfg->cfg_file);
    return 0;
  }
  
  
  return 1;
}



int daq_quit(){
  if( access( cfg->quit_file, F_OK ) != -1 ) {
    // file exists
    return 1;
  } else {
    // file doesn't exist
    return 0;
  } 
}


int write_config(char *fname,config_t *cfg){
  return 0;
}
