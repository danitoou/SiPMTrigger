#ifndef _config_h_
#define _config_h_
#include "CAENDigitizer.h"

#define IONIM 1
#define IOTTL 0

typedef struct config_s {
  char cfg_file[256];
  char quit_file[256];
  
  int nchan;
  short int chmask;
  short int chtrigmask;
  int choffset[8];
  float chped[8];

  int iolevel;
  int chtrigger;
  int exttrigger;
  float   trigthr;
  int useburst;
  int tr_overlap;
  int dual;

  int usedpp;
  
  int nev;  
  int irq;
  int nblt;

  int inburts;
  int handle;
}  config_t;

typedef struct
{
  CAEN_DGTZ_ConnectionType LinkType;
  uint32_t VMEBaseAddress;
  uint32_t RecordLength;
  uint32_t ChannelMask;
  int EventAggr;
  CAEN_DGTZ_PulsePolarity_t PulsePolarity;
  CAEN_DGTZ_DPP_AcqMode_t AcqMode;
  CAEN_DGTZ_IOLevel_t IOlev;
} DigitizerParams_t;



int read_config(config_t *cfg,char *fname);
int write_config(char *fname,config_t *cfg);
config_t * default_config();
int daq_quit();

extern config_t *cfg;
extern CAEN_DGTZ_DPP_PSD_Params_t DPPParams;
extern DigitizerParams_t Params;

#endif
