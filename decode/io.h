#ifndef _MY_IO_H_
#define _MY_IO_H_

//#define USE_CASTOR


#ifdef USE_CASTOR

#include "shift/rfio_api.h"
#define FOPEN   rfio_fopen
#define FCLOSE  rfio_fclose
#define STAT    rfio_stat
#define FREAD   rfio_fread
#define FEOF    rfio_feof


# else 

#define FOPEN   fopen
#define FCLOSE  fclose
#define STAT    stat
#define FREAD   fread
#define FEOF    feof


#endif


#endif
