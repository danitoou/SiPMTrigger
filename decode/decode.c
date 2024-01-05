#include "stdio.h"
#include "time.h"
#include "stdlib.h"
#include <string.h>
#include "math.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//#include "rfio_api.h"


//histogramming - ROOT:
#include "TROOT.h"
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMinuit.h"

#include "TNtuple.h"
#include "TProfile.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TLorentzVector.h>
#include <TArrow.h>

//#define USE_CASTOR

#include "io.h"

TROOT root("FADC","FADC");


//All sizes are defined in terms of words !
#define HDR_SIZE 4      //4 fields of 4 bytes:
#define MAX_DATA 10000000   //Maximum of 30e6 samples, 3 per word
#define MAX_EVENT_SIZE ((HDR_SIZE) + (MAX_DATA))


#define MAX_N_SAMPLES   30000000
#define MAX_N_HITS      10
#define MAX_N_CHANNELS  8

//Header bit definitions and offsets
#define HDR_EVT_SIZE    0x0fffffff
#define HDR_CHECK       0xf0000000
#define HDR_BOARD_ID    0xf8000000
#define HDR_RES         0x06000000
#define HDR_ZERO        0x01000000
#define HDR_PATTERN     0x00ffff00
#define HDR_CH_MASK     0x000000ff
#define HDR_RESERVED    0xff000000
#define HDR_EVT_COUNTER 0x00ffffff
#define HDR_TRIG_TIME   0xffffffff

#define OFF_HDR_EVT_SIZE    0
#define OFF_HDR_CHECK       28
#define OFF_HDR_BOARD_ID    27
#define OFF_HDR_RES         25      
#define OFF_HDR_ZERO        24
#define OFF_HDR_PATTERN     8
#define OFF_HDR_CH_MASK     0
#define OFF_HDR_RESERVED    24
#define OFF_HDR_EVT_COUNTER 0
#define OFF_HDR_TRIG_TIME   0


//Channel aggregate header:


#define AGD_HDR_SIZE      0x004fffff
#define AGR_HDR_SMPL      0x00000fff
#define AGR_HDR_DP        0x00040000
#define AGR_HDR_DT        0x80000000
#define AGR_HDR_EQ        0x40000000
#define AGR_HDR_ET        0x20000000
#define AGR_HDR_EB        0x10000000
#define AGR_HDR_ES        0x08000000



#define OFF_AGD_HDR_SIZE  0
#define OFF_AGR_HDR_SMPL  0
#define OFF_AGR_HDR_DP    16
#define OFF_AGR_HDR_DT    31
#define OFF_AGR_HDR_EQ    30
#define OFF_AGR_HDR_ET    29
#define OFF_AGR_HDR_EB    28
#define OFF_AGR_HDR_ES    27









#define DATA_HEADER     0xc0000000
#define DATA_SAMPLE     0x000003ff
#define DATA_NSBITS     10      //Number of bits in a sample
#define DATA_NSAMPLES   3       //Number of samples per word
//#define DATA_CH_END     0x80000000
#define DATA_CH_CONT    0xc0000000

//Signal properties:
#define CH0_PED 550
#define CH1_PED 600
#define PEDESTAL 1010
#define SIGNAL_THRESHOLD 5
#define SIGNAL_WIDTH  40 //Total signal width = SIGNAL_WIDTH + PED_SAMPLE_OFFSET/2
#define SIGNAL_PERCENTAGE 0.1
#define PED_SAMPLE_OFFSET 30 //Offset to start the calculation of PED wrt maximal signal
#define PED_MAXSAMPLE  20
#define TH_RATIO_1   0.1
#define TH_RATIO_2   0.9

#define E_CHARGE 1.6e-19

#define TOT_THRESH 10. //Threshold for time over threshold
#define NTHRESHOLDS 10
#define THRESH_STEP 3.
#define THRESH_START 70.

#define N_LAV_TH 0


static char buf[MAX_EVENT_SIZE];
float hv=0;
float gain=0;
int nfiles = 0;
int n1e_sac = 0;



typedef struct {
  unsigned int check;
  unsigned int eventsize;
  unsigned int size;
  unsigned int bid;
  unsigned int res;
  unsigned int zero;
  unsigned int pattern;
  unsigned int mask;
  unsigned int reserved;
  unsigned int evcnt;
  unsigned int time;
} eventHeader;


typedef struct {
  float time;
  float rise;
  float fall;
  float width;
} timeprop;


typedef struct {
  float localtime;
  float time;
  float height;
  float integral;
  float energy;
  float charge;
  float ne;
  float npe;
  float ped;
  int nsamples;
  int presamples;   //
  int postsamples; //Including also the maximal sample
} myhit;


typedef struct {
  Double_t par[10];
  Double_t err[10];
  Double_t chi2;
  Double_t ndf;
  int success;
}fitresults;



typedef struct{
  int nsamples;
  int n;
  short int sample[MAX_N_SAMPLES];
  int nhit;
  myhit hit[MAX_N_HITS];
  float ped;
  float ev_ped;
  float integral;
  float hv;
  float gain;

  myhit hit_thresh;
  myhit hit_width;
  myhit hit_relthresh;
  float tot;

  fitresults fit;

  float energy;
  float ev_energy;
  float cal;
  float max;
  int imax;
  unsigned int qs;
  unsigned int ql;
  unsigned int time;
  unsigned int base;
  float t;
  float t_2;
} channel;


typedef struct {
  int size;
  int nwords;
  int dt;    // Dual trace
  int eq;   //Charge is enabled
  int et;  //Time Tag is enabled
  int eb;  //Baseline is enabled
  int es;  //Waveform (samples) is enabled
  int dp;  //Digital Probe Selection
} evtaggrheader;

typedef struct {
  float charge;
  float charge_thresh;
  float charge_relthresh;
  float energy;  
  float time;  
  float nhits;

} sacdata;


typedef struct  {
  eventHeader hdr;
  short int nch;
  channel ch[8];
  short int nhit;
  myhit hit[MAX_N_HITS];
  int dpp;
  int ql;
  int qs;
  evtaggrheader aggrhdr;
  sacdata sac;
} myevent;




typedef struct {
  TFile *hisfile;
  TH1F  *nhit;
  TH1F *dt;
  TH1F *beam_rate;
  TH1F *beam_1e_rate;
  TH1F *rate;
  TH1F *tot_energy;
  TH1F *ch_all_trig[MAX_N_CHANNELS];
  // TH1F *ch_four_seven[MAX_N_CHANNELS];
  TH1F *ch_int[MAX_N_CHANNELS];
  TH1F *ch_energy[MAX_N_CHANNELS];
  TH1F *ch_ev_ped[MAX_N_CHANNELS];
  TH1F *ch_ev_energy[MAX_N_CHANNELS];
  TH2F *ch_ch_energy[(MAX_N_CHANNELS*(MAX_N_CHANNELS-1)) / 2];
  TH2F *ch_ch_energy_ev_ped[(MAX_N_CHANNELS*(MAX_N_CHANNELS-1)) / 2];
  
  TH1F *ch_ch_time_diff[(MAX_N_CHANNELS*(MAX_N_CHANNELS-1)) / 2];
  TH1F *ch_ch_time_diff_ecut[(MAX_N_CHANNELS*(MAX_N_CHANNELS-1)) / 2];

  TH1F *ch_imax[MAX_N_CHANNELS];
  TH1F *ch_imax_cut[MAX_N_CHANNELS];
  TH1F *ch_max[MAX_N_CHANNELS];
  TH1F *ch_max_ampl[MAX_N_CHANNELS];
  TH1F *ch_max_ampl_outt[MAX_N_CHANNELS];
  TH1F *ch_tot[MAX_N_CHANNELS];



  TH1F *ch_width_int[MAX_N_CHANNELS];
  TH1F *ch_thresh_int[MAX_N_CHANNELS];
  TH1F *ch_relthresh_int[MAX_N_CHANNELS];

  TH1F *ch_width_q[MAX_N_CHANNELS];
  TH1F *ch_thresh_q[MAX_N_CHANNELS];
  TH1F *ch_relthresh_q[MAX_N_CHANNELS];

  TH1F *ch_width_q_outt[MAX_N_CHANNELS];
  TH1F *ch_thresh_q_outt[MAX_N_CHANNELS];
  TH1F *ch_relthresh_q_outt[MAX_N_CHANNELS];


  TH1F *ch_width_ne[MAX_N_CHANNELS];
  TH1F *ch_thresh_ne[MAX_N_CHANNELS];
  TH1F *ch_relthresh_ne[MAX_N_CHANNELS];

  TH1F *ch_width_npe[MAX_N_CHANNELS];
  TH1F *ch_thresh_npe[MAX_N_CHANNELS];
  TH1F *ch_relthresh_npe[MAX_N_CHANNELS];

  
  TH1F *ch_width_ns[MAX_N_CHANNELS];
  TH1F *ch_thresh_ns[MAX_N_CHANNELS];
  TH1F *ch_relthresh_ns[MAX_N_CHANNELS];
  
  TH1F *ch_width_npre[MAX_N_CHANNELS];
  TH1F *ch_thresh_npre[MAX_N_CHANNELS];
  TH1F *ch_relthresh_npre[MAX_N_CHANNELS];
  
  TH1F *ch_width_npost[MAX_N_CHANNELS];
  TH1F *ch_thresh_npost[MAX_N_CHANNELS];
  TH1F *ch_relthresh_npost[MAX_N_CHANNELS];

  TH1F *ch_width_q_he[MAX_N_CHANNELS];
  TH1F *ch_thresh_q_he[MAX_N_CHANNELS];
  TH1F *ch_relthresh_q_he[MAX_N_CHANNELS];


  //Histograms for the fit results:
  TH1F *ch_fit_norm[MAX_N_CHANNELS];
  TH1F *ch_fit_tscint[MAX_N_CHANNELS];
  TH1F *ch_fit_tpmt[MAX_N_CHANNELS];
  TH1F *ch_fit_a[MAX_N_CHANNELS];
  TH1F *ch_fit_chi2[MAX_N_CHANNELS];
  TH1F *ch_fit_time[MAX_N_CHANNELS];
  TH1F *ch_fit_dtime[MAX_N_CHANNELS];

  TH2F *ch_fit_norm_vs_integral[MAX_N_CHANNELS];
  
  TH2F *ch_charge_vs_ampl[MAX_N_CHANNELS];
  TH2F *ch_charge_vs_tot[MAX_N_CHANNELS];
  TH2F *ch_ampl_vs_tot[MAX_N_CHANNELS];

  TH1F *ch_ampl_trig_thresh[N_LAV_TH][MAX_N_CHANNELS];
  TH1F *ch_charge_trig_thresh[N_LAV_TH][MAX_N_CHANNELS];
  TH1F *ch_time_trig_thresh[N_LAV_TH][MAX_N_CHANNELS];
  

  TProfile *ch_fit_norm_vs_integral_prof[MAX_N_CHANNELS];
  TProfile *ch_fit_tscint_vs_integral[MAX_N_CHANNELS];
  TProfile *ch_fit_tpmt_vs_integral[MAX_N_CHANNELS];
  TProfile *ch_fit_a_vs_integral[MAX_N_CHANNELS];

  TProfile *ch_fit_norm_vs_charge_prof[MAX_N_CHANNELS];
  TProfile *ch_fit_tscint_vs_charge[MAX_N_CHANNELS];
  TProfile *ch_fit_tpmt_vs_charge[MAX_N_CHANNELS];
  TProfile *ch_fit_a_vs_charge[MAX_N_CHANNELS];


  TH1F *ch_ampl_scan[NTHRESHOLDS];
  TH1F *ch_eff_all_ev;
  TH1F *ch_eff_pass_ev;
  

  TH1F *ch_fit_norm_ov_int[MAX_N_CHANNELS];
  
  TH1F *ch_fit_dt;

  TH1F *ch_12_dt;
  TH1F *ch_47_dt;


  //BTF histos

  //SAC total charge collected
  TH1F *sac_charge;
  TH1F *sac_charge_thresh;
  TH1F *sac_charge_relthresh;
  TH1F *sac_energy;

  TH1F *sac_charge_trigg;
  TH1F *sac_charge_trigg_thresh;
  TH1F *sac_charge_trigg_relthresh;

  TH1F *sac_energy_trigg;
  TH2F *sac_tr1_tr2_imax;

   
  //Information filled for each file
  TH1F *tmp_charge[MAX_N_CHANNELS];
  TGraphErrors *gr_channel_charge[MAX_N_CHANNELS];
  



} myhistos;

myhistos *histo;


myevent *evt;



//Function declarations:
int soft_trigger_ok(myevent *evt, int ich,float th);
int trigger_single_electron(myevent *evt);
float get_charge(float uint);

Double_t my_exp(Double_t *x,Double_t *par){
  
  return   par[0] + par[1]*exp(-(x[0] - par[3]) / par[2]);
}



Double_t signal_func(Double_t *x,Double_t *par){
  Double_t f;
  Double_t y;
  Double_t c;
  //Parameter description:
  //par[0] - Normalization (~Amplitude)
  //par[1] - start of the signal (t0)
  //par[2] - scintillator decay time (tau)
  //par[3] - c = (1/tau - 1/b), b~ decay const for the PMT
  //   OR
  //par[3] - b: ~ decay const for the PMT
  //par[4] - a:  ~ signal width for PMT
  
  if(x[0] < par[1]) {
    f=0.;
  } else {
    y = x[0] - par[1];
    c = (1./par[2]) - (1./par[3]);
    //c = par[3] ;


    /* f = - (par[0]/((par[4]*par[4]) + c * c ) ) * (exp(-y/par[2])/par[2]) * */
    /*   (exp(c*y)* */
    /*    ( (c*sin(par[4]*y)) - (par[4]*cos(par[4]*y)) )  + par[4] ) ; */


    
    if(par[4]*y < M_PI) {
      /* f = - par[0] * (exp(-y/par[2])/par[2]) *  */
      /* 	(exp(c*y)* */
      /* 	 ( (c*sin(par[4]*y)) - (par[4]*cos(par[4]*y)) )  + par[4] ) ; */
      
      f = - (par[0]/((par[4]*par[4]) + c * c ) ) * (exp(-y/par[2])/par[2]) *
        (exp(c*y)*
         ( (c*sin(par[4]*y)) - (par[4]*cos(par[4]*y)) )  + par[4] ) ;
      
    } else { 
      /* f = - par[0] * (exp(-y/par[2])/par[2]) *  */
      /* 	( (exp(c* M_PI/par[4])*par[4])+ par[4] ) ; */

      f = - (par[0]/((par[4]*par[4]) + c * c ) ) * (exp(-y/par[2])/par[2]) *
	( (exp(c* M_PI/par[4])*par[4])+ par[4] ) ;
      

    }
    
    

  }
  return f;
}

Double_t pmt_gain(Double_t *x,Double_t *par){
  //par[0] - constant
  //par[1] - power
  return par[0]*pow(x[0],par[1]);

}

void swap(unsigned int *a){
  int b=0;
  int i;
  for(i=0;i<4;i++) {
    b |= (( (*a & (0xff<<(8*i))) >> 8*i) << 8*(3-i));
  }
  printf("a: 0x%8x b: 0x%8x\n",*a,b);
  *a=b;
}





void histo_init(char *hfile){
  int i;
  char name[64];
  char title[64];
  int ihis;

  histo = new myhistos;
  histo->hisfile = new TFile(hfile,"recreate");
  histo->hisfile->cd();

  //Histograms initialization
  histo->nhit       =  new TH1F("nhit","Number of TDC hit data",40,0.0,40.0);
  histo->dt         =  new TH1F("ev_dt","Event time difference",1000000,0.0,1000000);
  histo->beam_rate  =  new TH1F("beam_rate","Rate in the SAC",10000,0.5,10000.5);
  histo->rate       =  new TH1F("rate","Rate in the SAC",1000,0,100000);

  histo->beam_1e_rate  =  new TH1F("beam_1e_rate","1e Rate in the SAC",10000,0.5,10000.5);


  ihis=0;
  
  histo->ch_12_dt = new TH1F("ch03_time_diff","Maximal sample time difference",100,-50.,50.);
  histo->ch_47_dt = new TH1F("ch47_time_diff","Maximal sample time difference",100,-50.,50.);


  
  for(i=0;i<MAX_N_CHANNELS;i++) {
    sprintf(name,"ch%d_pe_all",i);
    sprintf(title,"Photoelectrons in channel %d when all 4 channels are triggered",i);
    histo->ch_all_trig[i] = new TH1F(name,title,300,0,300);
    // else histo->ch_four_seven[i] = new TH1F(name,title,2200,0,300);
    sprintf(name,"ch%d_int",i);
    sprintf(title,"Integral of the signal in channel %d",i);
    histo->ch_int[i]   =  new TH1F(name, title ,157696,0.0,1024.*1540);

    sprintf(name,"ch%d_energy",i);
    sprintf(title,"Energy of the signal in channel %d",i);
    histo->ch_energy[i]  =  new TH1F(name, title,159696,-2000.0,1024*1540);

    sprintf(name,"ch%d_ev_ped",i);
    sprintf(title,"Pedestal of channel %d calculated event by event",i);
    //    histo->ch_ev_ped[i]  =  new TH1F(name, title,62*7,900.0,1024);
    histo->ch_ev_ped[i]  =  new TH1F(name, title,2048,0.0,1024);

    sprintf(name,"ch%d_ev_energy",i);
    sprintf(title,"Energy of the signal in channel %d using event pedestal",i);
    //histo->ch_ev_energy[i]  =  new TH1F(name, title,159696,-2000.0,1024*1540);
    histo->ch_ev_energy[i]  =  new TH1F(name, title,1024,-2000.0,1024*1540);

    for(int j=i+1;j<MAX_N_CHANNELS;j++) {
      sprintf(name,"ch%d_vs_ch%d_energy",i,j);
      sprintf(title,"Energy ch %d vs energy ch %d",i,j);
      histo->ch_ch_energy[ihis] = new TH2F(name,title,1020,0,1020,1020,0,1020);
      
      sprintf(name,"ch%d_vs_ch%d_energy_ev_ped",i,j);
      sprintf(title,"Charge ch %d vs charge ch %d using event pedestals",i,j);
      histo->ch_ch_energy_ev_ped[ihis] = new TH2F(name,title,1010,-2.e-11,200.e-12,1010,-2.e-11,200.e-12);
      
      sprintf(name,"ch%d_vs_ch%d_time",i,j);
      sprintf(title,"Time ch %d minus time ch %d",i,j);
      histo->ch_ch_time_diff[ihis] = new TH1F(name,title,200,-100.,100.);

      sprintf(name,"ch%d_vs_ch%d_time_ecut",i,j);
      sprintf(title,"Time ch %d minus time ch %d with charge cut",i,j);
      histo->ch_ch_time_diff_ecut[ihis] = new TH1F(name,title,200,-100.,100.);

      ihis++;      
    }
    
    sprintf(name,"ch%d_max_sample",i);
    histo->ch_imax[i] = new TH1F(name,"Index of maximal sample",1024,0,1024);

    sprintf(name,"ch%d_max_sample_cut",i);
    histo->ch_imax_cut[i] = new TH1F(name,"Index of maximal sample if signal",1024,0,1024);

    sprintf(name,"ch%d_max_ampl",i);
    histo->ch_max[i] = new TH1F(name,"Maximal amplitude",1024,0.0,1024);
    sprintf(name,"ch%d_max_ampl_no_ped",i);
    histo->ch_max_ampl[i] = new TH1F(name,"Maximal amplitude with subtracted pedestal",1024,0.0,1024);
    sprintf(name,"ch%d_max_ampl_outt",i);
    histo->ch_max_ampl_outt[i] = new TH1F(name,"Maximal amplitude if out of time",1024,0.0,1024);

    sprintf(name,"ch%d_tot",i);
    histo->ch_tot[i] = new TH1F(name,"Time over threshold",1000,0.0,1000);

    sprintf(name,"ch%d_width_int",i);
    histo->ch_width_int[i] = new TH1F(name,"Signal integral using fixed width",120,-200.0,1000.);
    sprintf(name,"ch%d_width_q",i);
    histo->ch_width_q[i] = new TH1F(name,"Total charge using fixed width",2200,-1.e-11,2e-10);
    sprintf(name,"ch%d_width_q_he",i);
    histo->ch_width_q_he[i] = new TH1F(name,"Total charge using fixed width",2200,-1.e-11,1e-9);
    sprintf(name,"ch%d_width_q_outt",i);
    histo->ch_width_q_outt[i] = new TH1F(name,"Total charge using fixed width out of time",101,-2.e-12,2.e-10);
    sprintf(name,"ch%d_width_ne",i);
    histo->ch_width_ne[i] = new TH1F(name,"Number of electrons using fixed width",100,0.0,1.e8);
    sprintf(name,"ch%d_width_npe",i);
    histo->ch_width_npe[i] = new TH1F(name,"Number of photoelectrons using fixed width",500,0.0,5.e2);
    sprintf(name,"ch%d_width_ns",i);
    histo->ch_width_ns[i] = new TH1F(name,"Number of samples to get the signal using fixed width",150,0.0,150);    
    sprintf(name,"ch%d_width_npre",i);
    histo->ch_width_npre[i] = new TH1F(name,"Number of pre-samples to get the signal using fixed width",50,0.0,50);    
    sprintf(name,"ch%d_width_npost",i);
    histo->ch_width_npost[i] = new TH1F(name,"Number of post-samples to get the signal using fixed width",100,0.0,100);    

    sprintf(name,"ch%d_charge_ampl",i);
    histo->ch_charge_vs_ampl[i] = new TH2F(name,"Charge versus amplitude",101,-2.e-12,2.e-10,102,0.0,1024);
    sprintf(name,"ch%d_charge_vs_tot",i);
    histo->ch_charge_vs_tot[i] = new TH2F(name,"Charge versus ToT",101,-2.e-12,2.e-10,100,0.0,1000);
    sprintf(name,"ch%d_ampl_vs_tot",i);
    histo->ch_ampl_vs_tot[i] = new TH2F(name,"Amplitude versus ToT",102,0.0,1024,100,0.0,1000.);



    
    sprintf(name,"ch%d_thresh_int",i);
    histo->ch_thresh_int[i] = new TH1F(name,"Signal integral using fixed thresh",120,-200.0,1000.);
    sprintf(name,"ch%d_thresh_q",i);
    histo->ch_thresh_q[i] = new TH1F(name,"Total charge using fixed thresh",1010,-2.e-12,2.e-10);
    sprintf(name,"ch%d_thresh_q_he",i);
    histo->ch_thresh_q_he[i] = new TH1F(name,"Total charge using fixed thresh",2200,-1.e-11,1e-9);
    sprintf(name,"ch%d_thresh_q_outt",i);
    histo->ch_thresh_q_outt[i] = new TH1F(name,"Total charge using fixed thresh out of time",1010,-2.e-12,2.e-10);
    sprintf(name,"ch%d_thresh_ne",i);
    histo->ch_thresh_ne[i] = new TH1F(name,"Number of electrons using fixed thresh",100,0.0,1.e8);
    sprintf(name,"ch%d_thresh_npe",i);
    histo->ch_thresh_npe[i] = new TH1F(name,"Number of photoelectrons using fixed thresh",5000,0.0,5.e3);
    sprintf(name,"ch%d_thresh_ns",i);
    histo->ch_thresh_ns[i] = new TH1F(name,"Number of samples to get the signal using fixed thresh",150,0.0,150);    
    sprintf(name,"ch%d_thresh_npre",i);
    histo->ch_thresh_npre[i] = new TH1F(name,"Number of pre-samples to get the signal using fixed thresh",50,0.0,50);    
    sprintf(name,"ch%d_thresh_npost",i);
    histo->ch_thresh_npost[i] = new TH1F(name,"Number of post-samples to get the signal using fixed thresh",150,0.0,150);    

    
    sprintf(name,"ch%d_relthresh_int",i);
    histo->ch_relthresh_int[i] = new TH1F(name,"Signal integral using fixed relthresh",120,-200.0,1000.);
    sprintf(name,"ch%d_relthresh_q",i);
    histo->ch_relthresh_q[i] = new TH1F(name,"Total charge using fixed relthresh",1010,-2.e-12,.1e-10);
    sprintf(name,"ch%d_relthresh_q_he",i);
    histo->ch_relthresh_q_he[i] = new TH1F(name,"Total charge using fixed relthresh",2200,-1.e-11,1e-9);
    sprintf(name,"ch%d_relthresh_q_outt",i);
    histo->ch_relthresh_q_outt[i] = new TH1F(name,"Total charge using fixed relthresh out of time",1010,-2.e-12,.1e-10);
    sprintf(name,"ch%d_relthresh_ne",i);
    histo->ch_relthresh_ne[i] = new TH1F(name,"Number of electrons using fixed relthresh",100,0.0,1.e8);
    sprintf(name,"ch%d_relthresh_npe",i);
    histo->ch_relthresh_npe[i] = new TH1F(name,"Number of photoelectrons using fixed relthresh",5000,0.0,5.e3);

    sprintf(name,"ch%d_relthresh_ns",i);
    histo->ch_relthresh_ns[i] = new TH1F(name,"Number of samples to get the signal using fixed relthresh",150,0.0,150);    
    sprintf(name,"ch%d_relthresh_npre",i);
    histo->ch_relthresh_npre[i] = new TH1F(name,"Number of pre-samples to get the signal using fixed relthresh",50,0.0,50);    
    sprintf(name,"ch%d_relthresh_npost",i);
    histo->ch_relthresh_npost[i] = new TH1F(name,"Number of post-samples to get the signal using fixed relthresh",150,0.0,150);    

    //Histos for LAV trigger threshold scans

    for(int kk = 0; kk < N_LAV_TH;kk++) {
      sprintf(name,"ch%d_max_ampl_no_ped_th%d",i,kk+5);
      histo->ch_ampl_trig_thresh[kk][i] = new TH1F(name,"Maximal amplitude with subtracted pedestal",1024,0.0,1024);
      sprintf(name,"ch%d_width_charge_th%d",i,kk+5);
      histo->ch_charge_trig_thresh[kk][i] = new TH1F(name,"Maximal amplitude with subtracted pedestal",1010,-2.e-12,2.e-10);
      sprintf(name,"ch%d_time_th%d",i,kk+5);
      histo->ch_time_trig_thresh[kk][i] = new TH1F(name,"Maximal amplitude with subtracted pedestal",256,0.0,256);
    }



    sprintf(name,"ch%d_fit_norm",i);
    histo->ch_fit_norm[i] = new TH1F(name,"Normalization factor",5200,-200.0,10000);
    sprintf(name,"ch%d_fit_tscint",i);
    histo->ch_fit_tscint[i] = new TH1F(name,"Decay time of scintillator",100,0.0,100);    
    sprintf(name,"ch%d_fit_tpmt",i);
    histo->ch_fit_tpmt[i] = new TH1F(name,"Decay constant connected with PMT",100,0.0,100);
    sprintf(name,"ch%d_fit_a",i);
    histo->ch_fit_a[i] = new TH1F(name,"PMT scale parameter",100,0.0,2);
    sprintf(name,"ch%d_fit_chi2",i);
    histo->ch_fit_chi2[i] = new TH1F(name,"Chi2 of the fit",100,0.0,100);

    sprintf(name,"ch%d_fit_norm_vs_integral",i);
    histo->ch_fit_norm_vs_integral[i] = new TH2F(name,"Norm vs Integral",520,-200.0,5000,520,-200.0,5000.);

    sprintf(name,"ch%d_fit_norm_vs_integral_prof",i);
    histo->ch_fit_norm_vs_integral_prof[i] = new TProfile(name,"Norm vs Integral",520,-200.0,5000,-200.0,5000.);
    sprintf(name,"ch%d_fit_tscint_vs_integral",i);
    histo->ch_fit_tscint_vs_integral[i] = new TProfile(name,"Decay time of scintillator vs INT",520,-200.0,5000., 2.,20.);
    sprintf(name,"ch%d_fit_tpmt_vs_integral",i);
    histo->ch_fit_tpmt_vs_integral[i] = new TProfile(name,"Decay time of PMT vs INT",520,-200.0,5000., 1.,20.);
    sprintf(name,"ch%d_fit_a_vs_integral",i);
    histo->ch_fit_a_vs_integral[i] = new TProfile(name,"PMT scale parameter vs INT",520,-200.0,5000., 0,2.);
    
    sprintf(name,"ch%d_fit_norm_vs_charge_prof",i);
    histo->ch_fit_norm_vs_charge_prof[i] = new TProfile(name,"Norm vs Charge",520,get_charge(-200.0),get_charge(5000),-200.0,5000.);
    sprintf(name,"ch%d_fit_tscint_vs_charge",i);
    histo->ch_fit_tscint_vs_charge[i] = new TProfile(name,"Decay time of scintillator vs Charge",520,get_charge(-200.0),get_charge(5000), 2.,20.);
    sprintf(name,"ch%d_fit_tpmt_vs_charge",i);
    histo->ch_fit_tpmt_vs_charge[i] = new TProfile(name,"Decay time of PMT vs Charge",520,get_charge(-200.0),get_charge(5000), 1.,20.);
    sprintf(name,"ch%d_fit_a_vs_charge",i);
    histo->ch_fit_a_vs_charge[i] = new TProfile(name,"PMT scale parameter vs Charge",520,get_charge(-200.0),get_charge(5000), 0,2.);




    sprintf(name,"ch%d_fit_norm_ov_int",i);
    histo->ch_fit_norm_ov_int[i] = new TH1F(name,"Fit normalization over INT",100,0.0,2.0);
    
    sprintf(name,"ch%d_fit_time",i);
    histo->ch_fit_time[i] = new TH1F(name,"Time of the hit",512,0.0,512.0);
    sprintf(name,"ch%d_fit_dtime",i);
    histo->ch_fit_dtime[i] = new TH1F(name,"Delta-Time of the hit",100,-50.0,50.0);


    //Histos to be reset for each file
    sprintf(name,"tmp_%d_charge",i);
    histo->tmp_charge[i] = new TH1F(name,"Total charge using fixed width",1010,-2.e-11,2.e-9);

    sprintf(name,"gr_%d_charge",i);
    histo->gr_channel_charge[i] = new TGraphErrors();


  }
  histo->ch_fit_dt = new TH1F("ch_fit_dt","Time difference",100,-20.,20.);
  
  histo->tot_energy = new TH1F("tot_energy","Energy in the event",10200,-20,1000);
  
  for(int i=0;i<NTHRESHOLDS;i++) {
    sprintf(name,"ch_ampl_scan_%dmv",(int) (THRESH_START + (i*THRESH_STEP)) );
    sprintf(title,"LeadGlass amplitude for %d mV trigger threshold", (int) (THRESH_START + (i*THRESH_STEP)) );
    histo->ch_ampl_scan[i] = new TH1F(name,title,1024,0.0,1024);
  }
  histo->ch_eff_all_ev = new TH1F("ch_eff_all_ev","Number of triggered events",NTHRESHOLDS,THRESH_START - (0.5*THRESH_STEP) ,THRESH_START + ((NTHRESHOLDS-0.5)*THRESH_STEP) );
  histo->ch_eff_pass_ev = new TH1F("ch_eff_pass_ev","Number of triggered events above the threshold",NTHRESHOLDS,THRESH_START - (0.5*THRESH_STEP) ,THRESH_START + ((NTHRESHOLDS-0.5)*THRESH_STEP) );
  
  
  histo->sac_charge = new TH1F("sac_charge","Total charge in SAC",100,0.0,300e-12);
  histo->sac_charge_thresh = new TH1F("sac_charge_thresh","Total charge in SAC",100,0.0,300e-12);
  histo->sac_charge_relthresh = new TH1F("sac_chargerelthresh","Total charge in SAC",100,0.0,300e-12);
  histo->sac_charge_trigg = new TH1F("sac_charge_trigg","Total charge in SAC if trigger",100,0.0,300e-12);
  histo->sac_charge_trigg_thresh = new TH1F("sac_charge_trigg_thresh","Total charge in SAC if trigger",100,0.0,300e-12);
  histo->sac_charge_trigg_relthresh = new TH1F("sac_charge_trigg_relthresh","Total charge in SAC if trigger",100,0.0,300e-12);
  
  histo->sac_energy = new TH1F("sac_energy","Total energy in SAC",100,0.0,1000.0);
  histo->sac_energy_trigg = new TH1F("sac_energy_trigg","Total energy in SAC if trigger",100,0.0,1000.0);

  histo->sac_tr1_tr2_imax = new TH2F("sac_tr1_tr2_imax","Trigger1 vs trigger2 time",100,0.0,500.0,100,0.0,500.0);

}

void histo_exit(){
  //  histo->dt->Write();
  char name[256];
  /* for(int i=0;i<MAX_N_CHANNELS;i++){ */
  /*   histo->ch_width_npe[i]->Scale(1./nfiles); */
  /*   histo->ch_thresh_npe[i]->Scale(1./nfiles); */
  /*   histo->ch_relthresh_npe[i]->Scale(1./nfiles); */
  /* } */
  for(int i=0;i<MAX_N_CHANNELS;i++) {
    sprintf(name,"gr_%d_charge",i);
    histo->gr_channel_charge[i]->Write(name);
  }
  histo->hisfile->Write(0);
  histo->hisfile->Close();
  delete histo->hisfile;
}

void burst_init(){
  int i=0;
  for(i = 0; i < MAX_N_CHANNELS;i++) {
    histo->tmp_charge[i]->Reset();
  }
}

void burst_exit(){
  //Perform operations at the end of each burst

  float charge;
  float ch_rms;
  int ent;
  int i;
  int ipoint;

  for(i = 0; i<MAX_N_CHANNELS;i++) {
    charge =  histo->tmp_charge[i]->GetMean();
    ch_rms = histo->tmp_charge[i]->GetRMS();
    ent = histo->tmp_charge[i]->GetEntries();
    ipoint = histo->gr_channel_charge[i]->GetN();
    histo->gr_channel_charge[i]->SetPoint(ipoint,ipoint+1,charge );
    histo->gr_channel_charge[i]->SetPointError(ipoint,0,ch_rms/sqrt(ent) );

  }
  
}


void draw_results(){
  printf("Drawing the results \n");
  
  root.Reset();
  
  TCanvas *c = new TCanvas("c","",0,0,800,600);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetGrid();
  
  histo->dt->Draw();
  
  sleep(2);
  
  return;
  
}



void fit_signal_shape(channel *ch){
  
  ch->fit.success = 0;
  
  Double_t f_beg = ch->imax - 10.;
  //Double_t f_end = ch->imax + 100.;
  Double_t f_end = ch->imax + 40.;

  
  if ( f_beg<0 || f_end >= ch->nsamples ) {
    return;
  }
  
  if(fabs(ch->max - ch->ev_ped) < 5.) return;

  TF1 *func = new TF1("fit_func",signal_func,f_beg,f_end,5);
  //  Double_t f_pars[5] = {ch->hit_width.integral*0.5,ch->imax-4.,4.2,4.1,0.32};
  //Double_t f_pars[5] = {ch->hit_width.integral*0.2,ch->imax-20.,5.5,4.,0.25};
  //  Double_t f_pars[5] = {ch->hit_width.integral, ch->imax-20 ,5.5,9.,0.3};
  //  Double_t f_pars[5] = {ch->hit_width.integral, ch->imax-20 ,10,20.,0.3};
  Double_t f_pars[5] = {0.22*ch->hit_width.integral, ch->imax-10 ,3.,3.,0.4};
  func->SetParameters(f_pars);
  /* func->FixParameter(1,ch->imax-12); */
  /* func->FixParameter(4,0.17); */
  /* func->FixParameter(2,15.1); */
  //func->FixParameter(3,2.1);
  //  func->FixParameter(0,ch->hit_width.integral*0.15);
  
  

  
  double x[ch->nsamples];
  double y[ch->nsamples];
  double xe[ch->nsamples];
  double ye[ch->nsamples];
 
  
  for (int i=0;i<ch->nsamples;i++) {
    x[i]=i;
    y[i]= ch->sample[i] - ch->ev_ped;
    ye[i] = 1;
    xe[i] = 1;

  }
  
  //  TGraph *gr=new TGraph(ch->nsamples, x , y);
  TGraphErrors *gr=new TGraphErrors(ch->nsamples, x , y,xe,ye);
  //  gMinuit->Command("SET PRINT -1");
  gr->Fit(func,"QR");
  
  for(int i=0;i<func->GetNpar();i++) {
    ch->fit.par[i] = func->GetParameter(i);
    ch->fit.err[i] = func->GetParError(i);
  }
  ch->fit.chi2 = func->GetChisquare();
  ch->fit.ndf = func->GetNDF();
  
  
  ch->fit.success=1;

  if(func) delete func;
  if(gr) delete gr;
  
}




void fill_evt_header(eventHeader *hdr,void *vbuf){
  unsigned int *buf = ( unsigned int *) vbuf;

  int i;
  /*
  for (i=0;i<4;i++) {
    //  swap(&buf[i]);
    printf("buf[%d] =  0x%8x \n",i,buf[i]);


  }
  */
  
  
  hdr->check       = ( buf[0] & HDR_CHECK       ) >> OFF_HDR_CHECK      ;
  hdr->eventsize   = ( buf[0] & HDR_EVT_SIZE    ) >> OFF_HDR_EVT_SIZE   ;
  hdr->bid         = ( buf[1] & HDR_BOARD_ID    ) >> OFF_HDR_BOARD_ID   ;
  hdr->res         = ( buf[1] & HDR_RES         ) >> OFF_HDR_RES        ;
  hdr->zero        = ( buf[1] & HDR_ZERO        ) >> OFF_HDR_ZERO       ;
  hdr->pattern     = ( buf[1] & HDR_PATTERN     ) >> OFF_HDR_PATTERN    ;
  hdr->mask        = ( buf[1] & HDR_CH_MASK     ) >> OFF_HDR_CH_MASK    ;
  hdr->reserved    = ( buf[2] & HDR_RESERVED    ) >> OFF_HDR_RESERVED   ;
  hdr->evcnt       = ( buf[2] & HDR_EVT_COUNTER ) >> OFF_HDR_EVT_COUNTER;
  hdr->time        = ( buf[3] & HDR_TRIG_TIME   ) >> OFF_HDR_TRIG_TIME  ;

}


int active_channel(int n, eventHeader *hdr ) {
  return (hdr->mask>>n) & 0x1;
}





void print_evt_hdr(eventHeader *hdr){
  printf(" Event check: \t 0x%1x \n Event size:\t %d \n Board id \t %d \n Channel mask:\t 0x%02x \n Event count:\t %d\n Event time:\t %u \n",
	 hdr->check     ,
	 hdr->eventsize ,
	 hdr->bid       ,
	 //hdr->res      
	 //hdr->zero     
	 //hdr->pattern  
	 hdr->mask      ,
	 //hdr->reserved 
	 hdr->evcnt     ,
	 hdr->time      
	 );
}


void fill_evt(myevent * evt,void *vbuf){
  unsigned int *buf = ( unsigned int *) vbuf;
  int i;
  int isample = 0;
  int ich=0;
  
  int nsamples;
  
  int nchannels=0;

  int ch_size;
  
  //  printf("Filling the event \n");  
  
  for (i=0;i< MAX_N_CHANNELS;i++) {
    if ( (evt->hdr.mask >> i ) & 0x1) 
      nchannels++;
  }
  
  //printf("Number of channels enabled in the event: %d\n",nchannels);

  ch_size = (evt->hdr.eventsize - HDR_SIZE) / nchannels;

  for (i=0;i < evt->hdr.eventsize - HDR_SIZE;i++) {
    // printf("Processing word %d\n",i);
    
    nsamples = (buf[i] & DATA_HEADER)>>30;
    for(int j=0;j<nsamples;j++) {      
      evt->ch[ich].sample[isample++] = (buf[i] >> (j*DATA_NSBITS)) & DATA_SAMPLE;
    }
    
    //printf("DATA header: 0x%8x \n"   ,buf[i] & DATA_HEADER );
    if( (i+1)%ch_size == 0  )  {
      evt->ch[ich].nsamples = isample;
      //Start filling the next channel
      //  printf("Done with channel %d\n",ich);
      ich++; 
      isample = 0;
    }
    if(ich>8) {
      printf("====ERROR:  Number of channels more than 8!\n");
      exit(0);
    }					      
  }
  evt->nch = ich; 
  
}


float maxAmplitude(channel *ch) {
  int i;
  ch->max = 1024; //minimal amplitude in fact, signal is negative
 
  for(int i=0;i<ch->nsamples;i++){
    if(ch->max > ch->sample[i]) {
      ch->imax = i;
      ch->max = ch->sample[i];
    }
  }
  return ch->max;
}




int plot_event(myevent *evt){

  char fname[256];
  static int j = 0;
  // if(maxAmplitude(&(evt->ch[4])) > 800  ||  (evt->ch[4].imax >250 ) ||(evt->ch[4].imax < 50 ) ) return 0;
  //if(maxAmplitude(&(evt->ch[0])) > 590) return 0;

  //Trigger OK
  // if( !soft_trigger_ok(evt,2,THRESH_START) ) return 0;

  //
  //if(! trigger_single_electron(evt)) return 0;

  //  if (evt->ch[6].ped - evt->ch[6].max < 20) return 0;

  //Selection on the LeadGlass channel
  if(1
     /* evt->ch[2].ev_ped-evt->ch[2].max < 90  */
     /* || evt->ch[2].ev_ped-evt->ch[2].max > 160 */
     //     && evt->ch[2].imax > 125. 
     //     && evt->ch[2].imax < 170. 
     
     
     //     evt->ch[2].imax < 125. ||  evt->ch[2].imax > 170. || evt->ch[0]. 
     ) {
    //Good events, no need to study more :)
    //    return 0;
  }
  
  


  root.Reset();
  
  // printf("Drawing the event\n");

  TCanvas *c = new TCanvas("c","",0,0,600,450);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetGrid();


  TH1F *his = c->DrawFrame(0,-200,1022,20,"");
  his->GetXaxis()->SetTitle("Time [ns]");
  his->GetYaxis()->SetTitle("Amplitude [mV]");
  his->Draw();
  //TH1F *his = c->DrawFrame(200,570,300,620,"");
  //TH1F *his = c->DrawFrame(290,570,350,620,"");

  //TH1F *his = c->DrawFrame(0,900,200,1024,"");

  //TH1F *his = new TH1I("his","Sample values",evt->ch[0].nsamples,
  
  float x[evt->ch[0].nsamples];  
  float xe[evt->ch[0].nsamples];

  float y[evt->ch[0].nsamples];
  float ye[evt->ch[0].nsamples];

  float n=0;

  float y2[evt->ch[0].nsamples];
  float y3[evt->ch[0].nsamples];
  float y4[evt->ch[0].nsamples];
  float y5[evt->ch[0].nsamples];
  float y6[evt->ch[0].nsamples];
  float y7[evt->ch[0].nsamples];
  float y8[evt->ch[0].nsamples];

  printf("======PEDESTALS========\n");
  for(int i = 0;i<8;i++) {
    printf("Pedestal for channel %d - %f\n",i, evt->ch[i].ev_ped);
  } 

  

  for (int i=0;i<evt->ch[0].nsamples;i++) {
    x[i]=i;
    y[i]= evt->ch[0].sample[i] - evt->ch[0].ev_ped;
    ye[i] = 0.5;
    xe[i] = 0.5;
    // printf("Sample[%d]: %d\n",i,y[i]);
    y2[i]= evt->ch[1].sample[i] - evt->ch[1].ev_ped;   
    y3[i]= evt->ch[2].sample[i] - evt->ch[2].ev_ped;   
    y4[i]= evt->ch[3].sample[i] - evt->ch[3].ev_ped;   
    y5[i]= evt->ch[4].sample[i] - evt->ch[4].ev_ped;   
    y6[i]= evt->ch[5].sample[i] - evt->ch[5].ev_ped;   
    //  printf( "ADC: %d\n" ,y7[i]= evt->ch[6].sample[i]); // - evt->ch[6].ev_ped;   
    y7[i]= evt->ch[6].sample[i] - evt->ch[6].ev_ped;   
    y8[i]= evt->ch[7].sample[i] - evt->ch[7].ev_ped;   


    
    n++;
  }
  
  
  //  printf("Arrays filled\n");
  
  his->Draw();  
  //  his->GetXaxis()->SetTitle("FADC");

  Double_t f_beg = evt->ch[1].imax - 10;
  //Double_t f_end = evt->ch[0].imax + 100.;
  Double_t f_end = evt->ch[1].imax + 40.;

  TF1 *func = new TF1("fit_func",signal_func,f_beg,f_end,5);
  //Double_t f_pars[5] = {0.3*evt->ch[2].hit_width.integral, evt->ch[2].imax-20 ,15.,7.,0.16};
  Double_t f_pars[5] = {0.22*evt->ch[1].hit_width.integral, evt->ch[1].imax-5 ,13,3.,0.17};
  printf("Integral/2: %f \n", evt->ch[1].hit_width.integral/3);
  func->SetParameters(f_pars);
  /* func->FixParameter(2,14.9);  */
  /* func->FixParameter(3,8.); */
  /* func->FixParameter(4,0.17);   */

  /*
  TF1 *func1 = new TF1("my_exp",my_exp,f_beg,f_end,4);
  Double_t f_parse[4] = {0, -evt->ch[0].max, 20, evt->ch[0].imax };
  func1->SetParameters(f_parse);
  */

  //  TGraph *gr=new TGraph(n, x , y3); 
  TGraphErrors *gr=new TGraphErrors(n, x , y,xe,ye);

  gr->SetMarkerSize(1.2);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(1);
  gr->Draw("P");

  printf("Ch. int: %f\n",evt->ch[2].hit_width.integral);
  // gr->Fit(func,"R");gr->Fit(func,"R");
  
  // func->Draw("same");


  

  TGraph *gr2=new TGraph(n, x , y2);

  gr2->SetMarkerSize(1.2);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerColor(2);
  gr2->Draw("P same");
  //gr2->Fit(func,"R");
  //  func->Draw("same");

  
  TGraph *gr3=new TGraph(n, x , y3);

  gr3->SetMarkerSize(1.2);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerColor(3);
  gr3->Draw("P same");

  TGraph *gr4=new TGraph(n, x , y4);

  gr4->SetMarkerSize(1.2);
  gr4->SetMarkerStyle(20);
  gr4->SetMarkerColor(4);
   gr4->Draw("P same");


   TGraph *gr5=new TGraph(n, x , y5);

  gr5->SetMarkerSize(1.2);
  gr5->SetMarkerStyle(20);
  gr5->SetMarkerColor(5);
  gr5->Draw("P same");


   TGraph *gr6=new TGraph(n, x , y6);

  gr6->SetMarkerSize(1.2);
  gr6->SetMarkerStyle(20);
  gr6->SetMarkerColor(6);
  gr6->Draw("P same");
  

   TGraph *gr7=new TGraph(n, x , y7);

  gr7->SetMarkerSize(1.2);
  gr7->SetMarkerStyle(20);
  gr7->SetMarkerColor(1);
  gr7->Draw("P same");
  //gr7->Draw("P");



  TGraph *gr8=new TGraph(n, x , y8);
   
  gr8->SetMarkerSize(1.2);
  gr8->SetMarkerStyle(20);
  gr8->SetMarkerColor(8);
  gr8->Draw("P same");

  

  c->Update();
  
  //c->Update();
  c->Modified();
  sprintf(fname,"tmp-plot/plot%d-low.gif",j);
  //  sprintf(fname,"tmp-lavblock/plot%d-high.pdf",j);
  //sprintf(fname,"pics/plot.gif");
  c->Print(fname);
  sprintf(fname,"tmp-plot/plot%d-low.pdf",j);
  //sprintf(fname,"pics/plot.gif");
  c->Print(fname);

  j++;

  //getchar();
  sleep(2);
  
  
  //removing the unnecessary part
  /*
  delete gr;
  delete gr2;
  delete gr3;
  delete gr4;
  delete gr5;
  delete his;
  delete c;
  delete func;
  */

  return 0;  
}

void init(myevent *evt){
  int i;

  for(i=0;i<MAX_N_CHANNELS;i++){
    evt->ch[i].ped = PEDESTAL;
  }
  
  evt->ch[0].ped = 1009.59;
  evt->ch[1].ped = 1008.65;
  evt->ch[2].ped = 1008.27;
  evt->ch[3].ped = 1009.74;
  

  evt->ch[0].cal = 332.6;
  evt->ch[1].cal = 283.5;
  evt->ch[2].cal = 214.5;
  evt->ch[3].cal = 283.4;
  
}


void calc_ped(channel *ch){
  int i;
  int ns;
  /*
    Using only the first 200 samples ..... not always available
  for(int i=0;i<200;i++) {
    ch->ev_ped += ch->sample[i];
  }
  ch->ev_ped/=200.;
  */
  

  //Start for pedestal calculation - 5 samples from the maximum
  if(ch->imax < PED_SAMPLE_OFFSET ) {
    //Use the default pedestal from a database
  } else {
    ch->ev_ped = 0;
    ns=0;
    while( (ns < PED_MAXSAMPLE ) && ( ch->imax-ns - PED_SAMPLE_OFFSET >= 0)) {
      ch->ev_ped+=ch->sample[ch->imax-ns-PED_SAMPLE_OFFSET];
      ns++;
    } 
    ch->ev_ped/=ns;
  }
}

float get_charge(float uint) {
  float r=50.;
  float dt = 1.e-9;  //1 ns time bin
  return uint*1.e-3*dt/r; //uint is in milivolts
}


void calc_int_width(channel *ch){
  int i;
  ch->gain = 6.44e5;
  ch->hit_width.integral = 0;
  ch->hit_width.nsamples = 0;
  ch->hit_width.ped = ch->ev_ped;
  ch->hit_width.presamples = PED_SAMPLE_OFFSET/2;
  ch->hit_width.postsamples = SIGNAL_WIDTH;


  for(i = ch->imax - PED_SAMPLE_OFFSET; i <ch->imax+SIGNAL_WIDTH;i++) {
    if(i>=0 && i<ch->nsamples  ) {
      ch->hit_width.integral += (ch->sample[i] - ch->hit_width.ped);
      ch->hit_width.nsamples++;
    }
  }
  
  ch->hit_width.integral*=-1.;

  
  ch->hit_width.charge = get_charge(ch->hit_width.integral);
  ch->hit_width.ne = ch->hit_width.charge/E_CHARGE;
  ch->hit_width.npe = ch->hit_width.ne/ch->gain;


}

void calc_int_thresh(channel *ch){
  int i;
  
  ch->hit_thresh.integral = 0;
  ch->hit_thresh.nsamples = 0;
  ch->hit_thresh.ped = ch->ev_ped;
  ch->hit_thresh.presamples = 0;
  ch->hit_thresh.postsamples = 0;

  //Going back:
  i = ch->imax-1;
  while( i > 0 && ( ch->sample[i] - ch->hit_thresh.ped < -SIGNAL_THRESHOLD) ) {
    ch->hit_thresh.integral += ( ch->sample[i] - ch->hit_thresh.ped);
    ch->hit_thresh.nsamples++;
    ch->hit_thresh.presamples++;
    i--;
  }
  
  i = ch->imax;
  while( i < ch->nsamples && ( ch->sample[i] - ch->hit_thresh.ped < -SIGNAL_THRESHOLD) ) {
    ch->hit_thresh.integral += (ch->sample[i] - ch->hit_thresh.ped);
    ch->hit_thresh.nsamples++;
    ch->hit_thresh.postsamples++;
    i++;
  }
  ch->hit_thresh.integral*=-1.;

  ch->hit_thresh.charge = get_charge(ch->hit_thresh.integral);
  ch->hit_thresh.ne = ch->hit_thresh.charge/E_CHARGE;
  ch->hit_thresh.npe = ch->hit_thresh.ne/ch->gain;

}




void calc_int_relthresh(channel *ch){
  int i;
  
  ch->hit_relthresh.integral = 0;
  ch->hit_relthresh.nsamples = 0;
  ch->hit_relthresh.ped = ch->ev_ped;
  ch->hit_relthresh.presamples=0;
  ch->hit_relthresh.postsamples=0;

  //  printf("====Getting the amplitude using relative thresholds ======== \n");
  //Going back:
  i = ch->imax-1;
  while( i > 0 && ( (1.*ch->sample[i] - ch->hit_relthresh.ped) < (SIGNAL_PERCENTAGE*(1.*ch->sample[ch->imax] -ch->hit_relthresh.ped) ) )  ) {
    //printf("Thresh: %f, Offset: %f, Sample[%d]=%d \n", (SIGNAL_PERCENTAGE*(1.*ch->sample[ch->imax] -ch->hit_relthresh.ped) )  , (1.*ch->sample[i] - ch->hit_relthresh.ped) , i, ch->sample[i]);;

    ch->hit_relthresh.integral += (ch->sample[i] - ch->hit_relthresh.ped );
    ch->hit_relthresh.nsamples++;
    ch->hit_relthresh.presamples++;
    i--;
  }
  i = ch->imax;
  while( i < ch->nsamples && ( (1.*ch->sample[i] - ch->hit_relthresh.ped) < (SIGNAL_PERCENTAGE*(1.*ch->sample[ch->imax] - ch->hit_relthresh.ped )) ) ) {
    
    //printf("Thresh: %f, Offset: %f, Sample[%d]=%d \n", (SIGNAL_PERCENTAGE*(1.*ch->sample[ch->imax] -ch->hit_relthresh.ped) )  , (1.*ch->sample[i] - ch->hit_relthresh.ped) , i, ch->sample[i]);;
    ch->hit_relthresh.integral += (1.*ch->sample[i] - ch->hit_relthresh.ped);
    ch->hit_relthresh.nsamples++;
    ch->hit_relthresh.postsamples++;
    i++;
  }
  ch->hit_relthresh.integral*=-1.;
  
  ch->hit_relthresh.charge = get_charge(ch->hit_relthresh.integral);
  ch->hit_relthresh.ne = ch->hit_relthresh.charge/E_CHARGE;
  ch->hit_relthresh.npe = ch->hit_relthresh.ne/ch->gain;
  
}


void calc_hit_time(channel *ch){
  float t1;
  float t2;

  float t3;
  float t4;

  
  float val1;
  float val2;

  int t1_ok=0;
  int t2_ok=0;

  int t3_ok=0;
  int t4_ok=0;

  float max = ( ch->ev_ped - ch->max);
  
  
  for(int i=0;i< ch->nsamples-1;i++) {
    float val1 = - (1.*ch->sample[i] - ch->ev_ped);
    float val2 = - (1.*ch->sample[i+1] - ch->ev_ped);
    
    if( t1_ok == 0 && val1 < TH_RATIO_1*max && val2 > TH_RATIO_1*max) {
      t1 = i +  (TH_RATIO_1*max - val1)/(val2 - val1)  ;
      t1_ok = 1;
    }
    if( t1_ok = 1 && t2_ok == 0 && val1 < TH_RATIO_2*max && val2 > TH_RATIO_2*max) {
      t2 = i+ (TH_RATIO_2*max - val1)/(val2 - val1)   ;
      t2_ok = 1;
    }

    if( t3_ok == 0 && val1 <= 20. && val2 >20.) {
      t3 = i + (20. - val1)/(val2 - val1);
      t3_ok = 1;
    }
    if( t3_ok = 1 && t4_ok == 0 && val1 <= 40. && val2 > 40.) {
      t4 = i + (40. - val1)/(val2 - val1);
      t4_ok = 1;
    }



  }
  
  ch->t = (t1 + t2)/2.;
  ch->t_2 = t3 - (t4-t3);
  
}


void calc_tot(channel *ch){
  int i;
  ch->tot = 0;
  if(ch->ev_ped-ch->max < TOT_THRESH) {
    ch->tot = -1.;
    return;
  }
  
  i = ch->imax-1;
  while( i > 0 && (ch->ev_ped - ch->sample[i] > TOT_THRESH) ) {
    ch->tot+=1.;
    i--;
  } 
  i = ch->imax;
  while( i < ch->nsamples  && (ch->ev_ped - ch->sample[i] > TOT_THRESH) ) {
    ch->tot+=1.;
    i++;
  }
  
}

void reco_channel(channel *ch){
  int i;

  ch->hv = hv;
  /* TF1 *pmtgain = new TF1("pmt_gain", pmt_gain,500,2500,2); */
  /* Double_t gain_pars[2] = {1.4973334268703837e-15,6.784231145378457}; */
  /* pmtgain->SetParameters(gain_pars); */
  /* ch->gain = pmtgain->Eval(ch->hv); */
  

  ch->gain = 6.44e5;

  maxAmplitude(ch);
  calc_ped(ch);
  calc_tot(ch);
  


  calc_int_width(ch);
  calc_int_thresh(ch);
  calc_int_relthresh(ch);
  
  calc_hit_time(ch);
  
  ch->integral=0;
  //  for(i=0;i<ch->nsamples;i++) {
  
  for(i=ch->nsamples-154;i<ch->nsamples;i++) {
    ch->integral += ch->sample[i];
  }
  
  //Define a signal - maximum till 
  ch->fit.success = 0;
  //  fit_signal_shape(ch);
  
  
  // ch->energy = (ch->nsamples*ch->ped - ch->integral);  
  // ch->ev_energy = (ch->nsamples*ch->ev_ped - ch->integral);  
  ch->energy = (154*ch->ped - ch->integral);  
  ch->ev_energy = (154*ch->ev_ped - ch->integral);  

  //if(pmtgain) delete pmtgain;
}



void reco_event(myevent *evt){
  int i;

  for(i=0;i<evt->nch;i++) {
    reco_channel(&(evt->ch[i]));
  }


  evt->sac.charge=0;
  evt->sac.charge_thresh=0;
  evt->sac.charge_relthresh=0;
  for (i=0;i<4;i++) {
    evt->sac.charge += evt->ch[i].hit_width.charge;
    evt->sac.charge_thresh += evt->ch[i].hit_thresh.charge;
    evt->sac.charge_relthresh += evt->ch[i].hit_relthresh.charge;
    
    evt->sac.energy += evt->ch[i].hit_relthresh.charge*1.e13;
  }
  
  
}


int soft_trigger_ok(myevent *evt, int ich,float th) {
  if(ich == 0 || ich == 1) return 1;
  
  if(ich == 2) 
    {
    //Channel 2 is the one under study:
    // Check the amplitudes of the scintillators

    if( evt->ch[0].ev_ped - evt->ch[0].max > th && evt->ch[1].ev_ped-evt->ch[1].max > th //Event amplitudes
	//	&& evt->ch[0].ev_ped - evt->ch[0].max < 20. && evt->ch[1].ev_ped-evt->ch[1].max < 20. 
	&& evt->ch[0].imax > 110  && evt->ch[0].imax < 140  //Ch 0 hit on time
	&& evt->ch[1].imax > 110  && evt->ch[1].imax < 140  //Ch 1 hit on time
	) {
      histo->ch_12_dt->Fill(evt->ch[0].imax - evt->ch[1].imax);
      if( fabs(evt->ch[0].imax - evt->ch[1].imax ) < 1.) 
      return 1;
    }     

  }
  return 0;
}

int soft_trigger_ok_lav(myevent *evt, int ich,float th) {
  if(ich == 0 || ich == 3 || ich == 4 || ich == 7 ) return 1;
  
  if(ich < 3) 
    {
      //Channel 1 & 2 are the ones under study:
      if( evt->ch[0].ev_ped - evt->ch[0].max > th && evt->ch[3].ev_ped-evt->ch[3].max > th //Event amplitudes
	  && evt->ch[0].imax > 80  && evt->ch[0].imax < 150  //Ch 0 hit on time
	  && evt->ch[3].imax > 80  && evt->ch[3].imax < 150  //Ch 1 hit on time
	  ) {
	if( th > 7.) 
	  histo->ch_12_dt->Fill(evt->ch[0].imax - evt->ch[1].imax);
	if( fabs(evt->ch[0].imax - evt->ch[1].imax ) < 5.) 
	  return 1;
      } 
      return 0; //No trigger
    }
      
  if(ich > 4) 
    {
      //Channel 5 & 6 are the ones under study:
      if( evt->ch[4].ev_ped - evt->ch[4].max > th && evt->ch[7].ev_ped-evt->ch[7].max > th //Event amplitudes
	  && evt->ch[4].imax > 80  && evt->ch[4].imax < 150  //Ch 0 hit on time
	  && evt->ch[7].imax > 80  && evt->ch[7].imax < 150  //Ch 1 hit on time
	  ) {
	if( th > 7.) 
	  histo->ch_47_dt->Fill(evt->ch[4].imax - evt->ch[7].imax);
	if( fabs(evt->ch[0].imax - evt->ch[1].imax ) < 5.) 
	  return 1;
      }       
      return 0;      
    }
  return 0;
}





int trigger_single_electron(myevent *evt) {
  
  int n = 1;

  //timing cut
  if (evt->ch[4].imax < 325. || evt->ch[4].imax > 350) return 0;
  //total charge cut:
  if (evt->ch[4].hit_width.charge < n*25e-12 || evt->ch[4].hit_width.charge  > n*45e-12) return 0;
  
  
  //timing cut
  if (evt->ch[5].imax < 325. || evt->ch[5].imax > 350) return 0;
  //total charge cut:
  //  if (evt->ch[5].hit_width.charge < 40e-12 || evt->ch[5].hit_width.charge  > 60e-12) return 0;
  if (evt->ch[5].hit_width.charge < n*25e-12 || evt->ch[5].hit_width.charge  > n*45e-12) return 0;
 
  if(evt->ch[4].imax - evt->ch[5].imax < -0.5 || evt->ch[4].imax - evt->ch[5].imax > 2.5 ) return 0;

  

  return 1;
}


int ana_event(myevent *evt){
  int i,j;
  int ihis;
  static int prevt=-1;

  if(prevt > 0) {
    histo->dt->Fill(evt->hdr.time - prevt);
  } 
  prevt = evt->hdr.time;
  //  return 0;

  reco_event(evt);
  
  //  if(evt->ch[0].imax < 115 || evt->ch[0].imax>135) return 0;

  histo->sac_charge->Fill(evt->sac.charge);
  
  if ( evt->sac.charge  > 30.e-12 &&  evt->sac.charge  < 50.e-12  ) {
    n1e_sac++; 
  }
  histo->sac_charge_thresh->Fill(evt->sac.charge_thresh);
  histo->sac_charge_relthresh->Fill(evt->sac.charge_relthresh);
  histo->sac_energy->Fill(evt->sac.energy);

  histo->sac_tr1_tr2_imax->Fill(evt->ch[4].imax,evt->ch[5].imax);

  histo->ch_fit_dt->Fill(evt->ch[4].imax - evt->ch[5].imax);
  
  if(1 ||  trigger_single_electron(evt)      ) {
    //Fill in the sac histos
    histo->sac_energy_trigg->Fill(evt->sac.energy);
    /* histo->ch_max[0]->Fill(evt->ch[0].ev_ped - evt->ch[0].max); */
    /* histo->ch_max[6]->Fill(evt->ch[6].ev_ped - evt->ch[6].max); */
    /* histo->ch_max[7]->Fill(evt->ch[7].ev_ped - evt->ch[7].max); */

    for(i=0;i<evt->nch;i++){

      histo->ch_max[i]->Fill(evt->ch[i].ev_ped - evt->ch[i].max);

      histo->ch_relthresh_int[i]->Fill(evt->ch[i].hit_relthresh.integral);
      histo->ch_relthresh_q[i]->Fill(evt->ch[i].hit_relthresh.charge);
      histo->ch_relthresh_q_he[i]->Fill(evt->ch[i].hit_relthresh.charge);
      histo->ch_relthresh_ne[i]->Fill(evt->ch[i].hit_relthresh.ne);
      histo->ch_relthresh_npe[i]->Fill(evt->ch[i].hit_relthresh.npe);
      histo->ch_relthresh_ns[i]->Fill(evt->ch[i].hit_relthresh.nsamples);
      histo->ch_relthresh_npre[i]->Fill(evt->ch[i].hit_relthresh.presamples);
      histo->ch_relthresh_npost[i]->Fill(evt->ch[i].hit_relthresh.postsamples);
      
      
      
    }

    if (evt->ch[7].hit_relthresh.charge <  0.3e-12) {
      histo->sac_charge_trigg->Fill(evt->sac.charge);
      histo->sac_charge_trigg_thresh->Fill(evt->sac.charge_thresh);
      histo->sac_charge_trigg_relthresh->Fill(evt->sac.charge_relthresh);
    }
    



  }

  


  ihis=0;
  /* histo->ch_ch_time_diff[0]->Fill(evt->ch[0].t - evt->ch[1].t); */
  /* histo->ch_ch_time_diff[1]->Fill(evt->ch[1].t - evt->ch[2].t); */

  /* histo->ch_ch_time_diff[2]->Fill(evt->ch[0].t_2 - evt->ch[1].t_2); */
  /* histo->ch_ch_time_diff[3]->Fill(evt->ch[1].t_2 - evt->ch[2].t_2); */

  for(i=0;i<evt->nch;i++){
    
    //check if the event selection for that channel is fine:
    if(!soft_trigger_ok(evt,i,THRESH_START)) {
      // continue;
    }
    
    //printf("Number of samples in channel %d: %d\n",i,evt->ch[i].nsamples);
    histo->ch_int[i]->Fill(evt->ch[i].integral);
    //printf("channel %d: integral: %f\n",i,evt->ch[i].integral);
    histo->ch_energy[i]->Fill(evt->ch[i].energy);
    
    // printf("channel %d: event pedestal: %f\n",i,evt->ch[i].ev_ped);
    histo->ch_ev_ped[i]->Fill(evt->ch[i].ev_ped);
    histo->ch_ev_energy[i] -> Fill(evt->ch[i].ev_energy);
    
    for(j=i+1;j<evt->nch;j++) {

//      histo->ch_ch_energy[ihis]->Fill(evt->ch[i].energy,evt->ch[j].energy);
//      histo->ch_ch_energy_ev_ped[ihis]->Fill(evt->ch[i].ev_energy,evt->ch[j].ev_energy);
      histo->ch_ch_energy[ihis]->Fill(evt->ch[i].max,evt->ch[j].max);
      if(1 || (evt->ch[i].imax < 88 && evt->ch[j].imax < 117)) {
	histo->ch_ch_energy_ev_ped[ihis]->Fill(evt->ch[i].hit_width.charge, evt->ch[j].hit_width.charge);
      }
      //      histo->ch_ch_time_diff[ihis]->Fill(evt->ch[i].t_2 - evt->ch[j].t_2);
      histo->ch_ch_time_diff[ihis]->Fill(evt->ch[i].t - evt->ch[j].t);

      ihis++;      

    }
    /* if(i!=0 && i!=6 && i!=7 ) */
    /*   histo->ch_max[i]->Fill(evt->ch[i].max); */
    
    histo->ch_imax[i]->Fill(evt->ch[i].imax);
    histo->ch_max_ampl[i]->Fill(evt->ch[i].ev_ped-evt->ch[i].max);
    histo->ch_tot[i]->Fill(evt->ch[i].tot);

    if(evt->ch[i].max < 990.) {
      histo->ch_imax_cut[i]->Fill(evt->ch[i].imax);
    }


    //Out  of time (noise) information:

    //    if(evt->ch[i].imax > 20  && evt->ch[i].imax < 320) {
    if(evt->ch[i].imax > 20 & evt->ch[i].imax < 400) {
      histo->ch_max_ampl_outt[i]->Fill(evt->ch[i].ev_ped-evt->ch[i].max);
      histo->ch_width_q_outt[i]->Fill(evt->ch[i].hit_width.charge);
      histo->ch_thresh_q_outt[i]->Fill(evt->ch[i].hit_thresh.charge);
      histo->ch_relthresh_q_outt[i]->Fill(evt->ch[i].hit_relthresh.charge);
    }
    
    histo->ch_width_int[i]->Fill(evt->ch[i].hit_width.integral);
    histo->ch_width_q[i]->Fill(evt->ch[i].hit_width.charge);
    histo->ch_width_q_he[i]->Fill(evt->ch[i].hit_width.charge);
    histo->tmp_charge[i] ->Fill(evt->ch[i].hit_width.charge);

    histo->ch_width_ne[i]->Fill(evt->ch[i].hit_width.ne);
    histo->ch_width_npe[i]->Fill(evt->ch[i].hit_width.npe);
    histo->ch_width_ns[i]->Fill(evt->ch[i].hit_width.nsamples);
    histo->ch_width_npre[i]->Fill(evt->ch[i].hit_width.presamples);
    histo->ch_width_npost[i]->Fill(evt->ch[i].hit_width.postsamples);
    histo->ch_charge_vs_ampl[i]->Fill(evt->ch[i].hit_width.charge,evt->ch[i].ev_ped-evt->ch[i].max);
    histo->ch_charge_vs_tot[i]->Fill(evt->ch[i].hit_width.charge,evt->ch[i].tot);
    histo->ch_ampl_vs_tot[i]->Fill(evt->ch[i].ev_ped-evt->ch[i].max,evt->ch[i].tot);

    histo->ch_thresh_int[i]->Fill(evt->ch[i].hit_thresh.integral);
    histo->ch_thresh_q[i]->Fill(evt->ch[i].hit_thresh.charge);
    histo->ch_thresh_q_he[i]->Fill(evt->ch[i].hit_thresh.charge);
    histo->ch_thresh_ne[i]->Fill(evt->ch[i].hit_thresh.ne);
    histo->ch_thresh_npe[i]->Fill(evt->ch[i].hit_thresh.npe);
    histo->ch_thresh_ns[i]->Fill(evt->ch[i].hit_thresh.nsamples);
    histo->ch_thresh_npre[i]->Fill(evt->ch[i].hit_thresh.presamples);
    histo->ch_thresh_npost[i]->Fill(evt->ch[i].hit_thresh.postsamples);


    //Histos for LAV block efficiency studies:

    for (int kk = 0; kk < N_LAV_TH; kk++) {
      if(soft_trigger_ok_lav(evt,i,5. + kk)) {
	histo->ch_ampl_trig_thresh[kk][i]->Fill(evt->ch[i].ev_ped-evt->ch[i].max);
	histo->ch_charge_trig_thresh[kk][i]->Fill(evt->ch[i].hit_width.charge);
	histo->ch_time_trig_thresh[kk][i]-> Fill(evt->ch[i].imax);

	
      }      
    }
    


      
    //if (i == 1  && evt->ch[i].hit_width.npe > 5.) {
    //}
    //continue;
    if(i == 9) {
      //For the LeadGlass channel:
      for(int j=0;j<NTHRESHOLDS;j++) {
	if(soft_trigger_ok(evt, i, THRESH_START + j*THRESH_STEP )) {
	  
	  histo->ch_ampl_scan[j]->Fill(evt->ch[i].ev_ped-evt->ch[i].max);
	  histo->ch_eff_all_ev->Fill( THRESH_START + j*THRESH_STEP  );
	  
	  if(evt->ch[i].ev_ped-evt->ch[i].max > 5. //Channel above a threshold
	     && evt->ch[i].imax > 125. 
	     && evt->ch[i].imax < 170. ) {
	    histo->ch_eff_pass_ev->Fill(THRESH_START + j*THRESH_STEP );
	  }
	  
	}
	
      }
      
    }
    
    //    if( !soft_trigger_ok(evt,2,THRESH_START) ) return 0;
    //    fit_signal_shape(&(evt->ch[i]) );
    if( evt->ch[i].fit.success
	//&& (evt->ch[i].fit.chi2/evt->ch[i].fit.ndf < 2.)  
	//&& fabs(evt->ch[0].max) > 20 && fabs(evt->ch[1].max) > 20 ){
	//	&& evt->ch[0].hit_width.integral>100. && evt->ch[1].hit_width.integral>100. ) {
	//&& evt->ch[i].hit_width.npe > 800
	//	&& evt->ch[1].hit_width.npe > 5. 
	//	&& evt->ch[2].fit.chi2 < 40 
	&& evt->ch[i].ev_ped - evt->ch[i].max > 5.
	//	&& evt->ch[i].fit.par[0] < 200
	){
      //Successful fit of the signal
      histo->ch_fit_norm[i]->Fill(evt->ch[i].fit.par[0]);
      histo->ch_fit_norm_vs_integral[i]->Fill(evt->ch[i].fit.par[0],evt->ch[i].hit_width.integral);

      histo->ch_fit_norm_vs_integral_prof[i]->Fill(evt->ch[i].hit_width.integral,evt->ch[i].fit.par[0]);

      histo->ch_fit_norm_vs_charge_prof[i]->Fill(evt->ch[i].hit_width.charge,evt->ch[i].fit.par[0]);

      double ratio = evt->ch[i].fit.par[0]/evt->ch[i].hit_width.integral;
      histo->ch_fit_norm_ov_int[i]->Fill(ratio);
      histo->ch_fit_time[i]->Fill(evt->ch[i].fit.par[1]);
      histo->ch_fit_dtime[i]->Fill(evt->ch[i].fit.par[1] - evt->ch[i].imax);
      histo->ch_fit_chi2[i]->Fill(evt->ch[i].fit.chi2);

      //Consistency check: - central value is arround 0.5
      //      if(1){
      if( 1
	  //&& fabs (ratio-0.5) < 0.2  
	 ){
	 // && evt->ch[0].hit_width.integral>100. && evt->ch[1].hit_width.integral>100. ) {	
	histo->ch_fit_tscint[i]->Fill(evt->ch[i].fit.par[2]);
	histo->ch_fit_tpmt[i]->Fill(evt->ch[i].fit.par[3]);
	histo->ch_fit_a[i]->Fill(evt->ch[i].fit.par[4]);
	histo->ch_fit_tscint_vs_integral[i]->Fill(evt->ch[i].hit_width.integral,evt->ch[i].fit.par[2]);
	histo->ch_fit_tpmt_vs_integral[i]->Fill(evt->ch[i].hit_width.integral,evt->ch[i].fit.par[3]);
	histo->ch_fit_a_vs_integral[i]->Fill(evt->ch[i].hit_width.integral,evt->ch[i].fit.par[4]);

	histo->ch_fit_tscint_vs_charge[i]->Fill(evt->ch[i].hit_width.charge,evt->ch[i].fit.par[2]);
	histo->ch_fit_tpmt_vs_charge[i]->Fill(evt->ch[i].hit_width.charge,evt->ch[i].fit.par[3]);
	histo->ch_fit_a_vs_charge[i]->Fill(evt->ch[i].hit_width.charge,evt->ch[i].fit.par[4]);
      }
    }    
    

  }
  
  if(evt->ch[0].fit.success 
     // && (evt->ch[0].fit.chi2/evt->ch[0].fit.ndf < 1.5)
     &&evt->ch[1].fit.success 
     //     && (evt->ch[1].fit.chi2/evt->ch[1].fit.ndf < 1.5)  
     //     && evt->ch[0].hit_width.integral>100. && evt->ch[1].hit_width.integral>100. 
     // && evt->ch[0].hit_width.npe > 800. && evt->ch[1].hit_width.npe > 800.
     ) {
    
    // histo->ch_fit_dt->Fill(evt->ch[4].fit.par[1] - evt->ch[5].fit.par[1]);
  }
  

  evt->hit[0].energy = 0;
  for(i=0;i<evt->nch;i++){
    evt->hit[0].energy+= evt->ch[i].ev_energy/evt->ch[i].cal;
  }
  
  histo->tot_energy->Fill(evt->hit[0].energy);
  //histo->nhit->Fill(evt->nhit);

  //plot data:

  if(evt->ch[0].imax >500 && evt->ch[0].imax < 520 && evt->ch[0].ev_ped - evt->ch[0].max > 20  ) {
    //  plot_event( evt);
  }
  
  ihis=0;
  for(int i = 0; i < evt->nch; i++) {
    for(int j = i+1; j < evt->nch; j++) {
      if(evt->ch[i].hit_width.charge > 15e-12 && evt->ch[j].hit_width.charge > 15e-12) {
	      histo->ch_ch_time_diff_ecut[ihis]->Fill(evt->ch[i].t - evt->ch[j].t);
      }
      ihis++;
    }   
  }

  int voltage = 56;
  float gain[] = {
    736055.9375*voltage-39378268.28125,
    742631.875*voltage-39621693.125,
    696619.0625*voltage-37172320,
    691407.8125*voltage-36945108.59375,
    836545.3125*voltage-44894339.0625,
    651700.9375*voltage-34622918.59375,
    817228.5625*voltage-43756791.125,
    772854*voltage-41316392.40625
    };

  bool all_4_ch = true;
  for(int i = 0; i < 4; i++) {
    if(evt->ch[i].hit_width.charge < 15e-12) {
      all_4_ch = false;
    }
  }
  if(all_4_ch) {
    for(int i = 0; i < 4; i++) {
      histo->ch_all_trig[i]->Fill(evt->ch[i].hit_width.charge/((1.6e-19)*gain[i]));
    }
  }

  all_4_ch = true;
  for(int i = 4; i < evt->nch; i++) {
    if(evt->ch[i].hit_width.charge < 15e-12) {
      all_4_ch = false;
    }
  }
  if(all_4_ch) {
    for(int i = 4; i < 8; i++) {
      histo->ch_all_trig[i]->Fill(evt->ch[i].hit_width.charge/((1.6e-19)*gain[i]));
    }
  }

  
  
  
  



  return 0;

}

void print_evt(myevent *evt){
  int i,j;
  printf("\n");
  printf("Number of channels active in this event: %d\n", evt->nch);

  int ch_ped[8];
  ch_ped[0]=553;
  ch_ped[1]=604;
 


  for(i=0;i<evt->nch;i++){
    printf("-------- Channel %d content------------\n",i);
    int print = 1;
    for(j=0;j<evt->ch[i].nsamples;j++) {
      //if(fabs(evt->ch[i].sample[j] - ch_ped[i]) > 10) print = 20;
      if(print>0){
	//	printf("%10d\t %5d\n",j,evt->ch[i].sample[j]);
	printf("%5d",evt->ch[i].sample[j]);
	//	print--;
      }
    }
    printf("\n");
  }
}

int read_event(FILE *fptr, myevent *evt) {
  unsigned int eventsize;
  static int size=0;
  static int n=0;

  //Clearing the event:
  memset(evt,0,sizeof(myevent));
  //Get event header:
  if (fread(buf,4,HDR_SIZE,fptr) != HDR_SIZE ) return 0;
  //process the event header
  //  printf("================ New event ================= \n");
  fill_evt_header(&(evt->hdr),buf);
  // print_evt_hdr(&(evt->hdr));
  
  size = ( evt->hdr.eventsize - HDR_SIZE);
  
  if (fread(buf,4*size,1,fptr) != size ) return 0;
  fill_evt(evt,buf);
  // print_evt(evt);
 
  //return 0; //stop after the first event
  return evt->hdr.eventsize;
}



int print_ch_aggr_hdr(evtaggrheader *hdr){
  printf("=======Channel aggregate header========\n");
  printf("Agregate size:   %d\n",hdr->size);
  printf("Number of words: %d\n",hdr->nwords);
  printf("Dual trace:      %d\n",hdr->dt);
  printf("Digital Probe:   %d\n",hdr->dp);
  printf("Time tag:        %d\n",hdr->et);
  printf("Charge:          %d\n",hdr->eq);
  printf("Baseline:        %d\n",hdr->eb);
  printf("Waveform:        %d\n",hdr->es);
  return 0;
}

int fill_ch_aggr_hdr(evtaggrheader *hdr,int *buf){

  hdr->size      =   (  buf[0] & AGD_HDR_SIZE  )     >> OFF_AGD_HDR_SIZE   ;
  hdr->nwords    =   (  buf[1] & AGR_HDR_SMPL  )     >> OFF_AGR_HDR_SMPL   ;
  hdr->dp        =   (  buf[1] & AGR_HDR_DP    )     >> OFF_AGR_HDR_DP     ;
  hdr->dt        =   (  buf[1] & AGR_HDR_DT    )     >> OFF_AGR_HDR_DT     ;
  hdr->eq        =   (  buf[1] & AGR_HDR_EQ    )     >> OFF_AGR_HDR_EQ     ;
  hdr->et        =   (  buf[1] & AGR_HDR_ET    )     >> OFF_AGR_HDR_ET     ;
  hdr->eb        =   (  buf[1] & AGR_HDR_EB    )     >> OFF_AGR_HDR_EB     ;
  hdr->es        =   (  buf[1] & AGR_HDR_ES    )     >> OFF_AGR_HDR_ES     ;
  return hdr->size;
}



int fill_dpp_event(myevent *evt,int *buf,int ich){
  int size=0;
  int i=0;
  int isample = 0;

  //Time stamp
  if(evt->aggrhdr.et) {
    evt->ch[ich].time = buf[size];
    size++; //1 word for the trigger time tag
  }
  //Samples
  if(evt->aggrhdr.es) {
    evt->ch[ich].nsamples = evt->aggrhdr.nwords*4 * 3; //Number of samples = nwords * 3s/w
    //Fill in the samples information
    
    /* for(i=0;i<evt->aggrhdr.nwords*4;i++) { */
    /*   for(int j=0;j<3;j++) {       */
    /* 	evt->ch[ich].sample[isample++] = (buf[size+i] >> (j*DATA_NSBITS)) & DATA_SAMPLE; */
    /*   } */
    /* } */
    /* size += evt->aggrhdr.nwords*4 ; */

    while(i < evt->ch[ich].nsamples ) {
      evt->ch[ich].sample[i] =  (buf[size] >> ( (i%3) * DATA_NSBITS) ) & DATA_SAMPLE  ;
      i++;
      if(i%3 == 0) size++; //we have read 3 samples, go to the next
    }
    if(i%3!=0) size++;
  }

  printf("Recorded samples: %d\n",i);
  
  
  //Baseline
  if(evt->aggrhdr.eb) {
    evt->ch[ich].base = buf[size];
    size++;
  }

  //Qshort, Qlong
  if(evt->aggrhdr.eq) {
    evt->ch[ich].qs = buf[size] && 0x00007fff;
    evt->ch[ich].ql = (buf[size] >> 16 ) && 0xffff;
    size++;
  }
  printf("Total number of words processed: %d\n",size);

  return size;
}


int process_channel_aggregate(myevent *evt,int *buf,int ich){
  int size=0;
  fill_ch_aggr_hdr(&(evt->aggrhdr),buf);
  print_ch_aggr_hdr(&(evt->aggrhdr));
  size = 2; //2 words are reserved for the aggregate header
  
  while (size != evt->aggrhdr.size) {
    size+=fill_dpp_event(evt,&(buf[size]),ich);
    ana_event(evt);
  }
  
  return evt->aggrhdr.size;
}

int print_dpp_hdr(myevent *evt){
  printf("========DPP Event header ============\n");
  printf("DPP aggregate size in words: %d\n",evt->hdr.eventsize);
  printf("DPP event count: %d\n",evt->hdr.evcnt);
  printf("DPP Channel mask: %d\n",evt->hdr.mask);
  printf("DPP event time tag: %d\n",evt->hdr.time);
  printf("======================================\n");
}

int process_event(void *vbuf,myevent *evt){
  int *buf = (int *) vbuf;
  int i;
  int size=0;
  int ich=0;

  //  memset(evt,0,sizeof(myevent)); 
  fill_evt_header(&(evt->hdr),buf);  
  
  if(evt->dpp) { 
    print_dpp_hdr(evt);
    size += 4; //number of 32bit words read up to the moment - just header
    // process_evt_dpp();
    while( size != evt->hdr.eventsize) {
      size += process_channel_aggregate(evt,&(buf[size]),ich);
      ich++;
    }
    printf("Number of channels in this aggregate: %d\n",ich);
  } else {
    // print_evt_hdr(&evt->hdr);
    fill_evt(evt,buf+4);
    ana_event(evt);    
  }
  return 4*evt->hdr.eventsize;
}


int process_file(char *fname){
  FILE *fptr;
  unsigned int nevents=0;
  unsigned int word;
  //  char * filebuf;
  static char filebuf[1000000000];

  struct stat fstat;
  unsigned int currpos=0;
  static int nfile=0;
  
  char dpp[] = "dpp";

  /* evt = (myevent *) malloc(sizeof(myevent)); */
  memset(evt,0,sizeof(myevent)); 
  init(evt); 

  burst_init();

  //Getting info for the input file:
  STAT(fname,&fstat);
  printf("File: %s ; size: %d\n",fname,fstat.st_size);
  //  return 0;
  
  if(strstr(fname,dpp)!=NULL ) {
    evt->dpp=1;
    printf("File with DPP-PSD events \n");
  }
  
  //Prepare the buffers:
  //filebuf = (char *) malloc(fstat.st_size);
  if( filebuf == NULL) return 0;
  
  if((fptr=FOPEN(fname,"r"))==NULL) return -1;
  printf("File  \"%s\" opened successfully\n",fname);
  
  //reading the file:
  FREAD(filebuf,fstat.st_size,1,fptr);
  
  FCLOSE(fptr);

  printf("File read in a buffer\n");

  n1e_sac = 0;

  
  while (currpos != fstat.st_size ) {
    
    currpos += process_event(&filebuf[currpos],evt);
    nevents++;

    if(nevents < 10) {
      // print_evt(evt);
    }
    if( nevents % 10000 == 0 ) {
      printf("%d events read\n",nevents);
    }
  }
  

  float nsec=1.;
  
  histo->beam_1e_rate->Fill(nfile,(1.*n1e_sac )/nevents);
  
  histo->beam_rate->SetBinContent(nfile,nevents/nsec);
  histo->rate->Fill(nevents/nsec);

  printf("Processed %ld events in file: %s \t event counter: %d, SAC 1e events: %f\n",nevents,fname,evt->hdr.evcnt,(1.*n1e_sac )/nevents);

  // if(filebuf) free(filebuf);
  //if(evt) free(evt);
  
  burst_exit();
  
  nfile ++;

  return nevents;

  //return 0;
  // LOOP OVER THE EVENTS
  while (read_event(fptr,evt) > 0  ) {
    //ana_event(evt);
    nevents++;
    if( nevents% 100 == 0 ) {
      printf("%d events read\n",nevents);
    }
    /* if( nevents% 200 == 0 ) { */
    /*   exit(0); */
    /* } */
    
  }
  fclose(fptr);

  printf("Processed %ld events in file: %s\n",nevents,fname);
  return 0;
}



int main(int argc, char **argv) {
  int usage=0;
  int i;
  FILE *flist;
  char fname[256];

  int nopt=0;
  char *opts[20];
  TApplication theApp("App",&nopt,opts);

  long int nevents = 0;

  /*
  printf("Size of int - %d\n",sizeof(int));
  printf("Size of unsigned int - %d\n",sizeof(unsigned int));
  printf("Size of long int - %d\n",sizeof(long int));
  printf("Size of short int - %d\n",sizeof(short int));
  printf("Size of char - %d\n",sizeof(char));
  //printf("Size of void - %d\n",sizeof(void));

 // printf("EvtIdx_data = %d\n",EvtIdx_data);
 */

  evt = (myevent *) malloc(sizeof(myevent));
  memset(evt,0,sizeof(myevent)); 
  init(evt);

  long int totnevents=0;
  histo_init("histos.root");


  if((argc==3)&&(strcmp(argv[1],"-i")==0))
    {
      printf("Operating on a data file \n");
      usage++;
      totnevents = i = process_file(argv[2]);

      if(i < 0) {
	printf("It is not possible to open file %s for reading\n",argv[2]);
	exit(0);
      } else {
	nfiles++;
      }
    }

  if((argc==3)&&(strcmp(argv[1],"-l")==0)) {
    usage++;
    printf("Operating on a list file \n");

    if( (flist = fopen(argv[2],"r")) == NULL) {
        printf("Cannot open the list of input files\n");
        return 0;
    }
    
    while(fscanf(flist,"%s",fname)!=EOF) {
      if(fname[0] == '%') {
	fscanf(flist,"%f",&hv);
	printf("High voltage used for the list: %f\n",hv);
      } else {
	if(fname[0] != '#')  {
	  totnevents += process_file(fname);
	  nfiles++;
	}
      }
    }
    
  }
  
  printf("Processed %d files\n",nfiles);
  printf("Analised %ld events \n",totnevents);

  if(usage==0)
    {
      printf("Usage of the programme: \n \n");
      printf("%s -i file_name\n\n",argv[0]);
      printf("%s -l file_list\n\n",argv[0]);
    }

  // draw_results();

  if(evt) free(evt);

  histo_exit();

  return 0;
}
