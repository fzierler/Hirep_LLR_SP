/****************************************************************************
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Main HMC program
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "ranlux.h"
#include "geometry.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "dirac.h"
#include "logger.h"
#include "llr_hmc_utils.h"
#include "memory.h"
#include "communications.h"
#include "observables.h"
#include "utils.h"
#include "spectrum.h"
#include "cinfo.c"
#include "wilsonflow.h"


/* LLR parameters */
typedef struct _input_llr {
  char make[256];
  int nmc,nth,it;
  double starta,S0,dS;
  /* for the reading function */
  input_record_t read[8];
} input_llr;




#define init_input_llr(varname) \
  { \
  .read={\
    {"make llr iterations", "llr:make = %s", STRING_T, &((varname).make)}, \
    {"Number of MC steps per RM iteration ", "llr:nmc = %d", INT_T, &((varname).nmc)}, \
    {"Number of MC therm steps per RM iteration", "llr:nth = %d", INT_T, &((varname).nth)}, \
    {"Initial a", "llr:starta = %lf", DOUBLE_T, &((varname).starta)}, \
    {"Robbins Monro startint iteration", "llr:it = %d", INT_T, &((varname).it)}, \
    {"Cental action", "llr:S0 = %lf", DOUBLE_T, &((varname).S0)}, \
    {"Delta S", "llr:dS = %lf", DOUBLE_T, &((varname).dS)}, \
    {NULL, NULL, 0, NULL}				\
    }\
}
 
input_llr llr_var=init_input_llr(llr_var);

typedef struct _input_WF {
  double tmax;
  int nmeas;
  int nint;
  char make[20];

  /* for the reading function */
  input_record_t read[5];

} input_WF;

#define init_input_WF(varname) \
  { \
  .read={\
    {"make WF", "WF:make = %s", STRING_T, &((varname).make)}, \
    {"WF max integration time", "WF:tmax = %lf", DOUBLE_T, &((varname).tmax)},\
    {"WF number of measures", "WF:nmeas = %d", DOUBLE_T, &((varname).nmeas)},\
    {"WF number of integration steps between measures", "WF:nint = %d", INT_T, &((varname).nint)},\
    {NULL, NULL, 0, NULL}\
    }\
  }

input_WF WF_var = init_input_WF(WF_var);

hmc_flow flow=init_hmc_flow(flow);

char input_filename[256] = "input_file";
char output_filename[256] = "out_0";
char error_filename[256] = "err_0";

static void read_cmdline(int argc, char* argv[]) {
  int i, ai=0, ao=0, am=0, requested=1;

  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) {ai=i+1;requested+=2;}
    else if (strcmp(argv[i],"-o")==0) {ao=i+1;requested+=2;}
    else if (strcmp(argv[i],"-m")==0) {am=i;requested+=1;}
  }

  if (am != 0) {
    print_compiling_info();
    exit(0);
  }

  error(argc!=requested,1,"read_cmdline [hmc.c]",
      "Arguments: [-i <input file>] [-o <output file>] [-m]");

  if (ao!=0) strcpy(output_filename,argv[ao]);
  if (ai!=0) strcpy(input_filename,argv[ai]);
}



int main(int argc,char *argv[]) {
  struct timeval startmain, endmain, etimemain; /* //for trajectory timing */
  gettimeofday(&startmain,0);
  int i;
  char sbuf[128];

  read_cmdline(argc,argv);
  
  /* setup process communications */
  setup_process(&argc,&argv);
  
  /* read global variables file */
  read_input(glb_var.read,input_filename);
  
  setup_replicas();
  
  /* logger setup */
  read_input(logger_var.read,input_filename);
  logger_set_input(&logger_var);
  if (PID!=0) { logger_disable(); }   /* disable logger for MPI processes != 0 */
  else {
    FILE* stderrp;
    sprintf(sbuf,">>%s",output_filename);  logger_stdout(sbuf);
    stderrp=freopen(error_filename,"w",stderr);
    error(stderrp==NULL,1,"main [hmc.c]",
	  "Cannot redirect the stderr");
  }
  
  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS);
  lprintf("MAIN",0,"[RepID: %d][world_size: %d]\n[MPI_ID: %d][MPI_size: %d]\n",RID,WORLD_SIZE,MPI_PID,MPI_WORLD_SIZE);
  lprintf("MAIN",0,"SVN Revision: %d\n", CI_svnrevision);

  //  lprintf("MAIN",0,"Logger lelvel: %d\n",logger_getlevel(0));
  
  /* setup lattice geometry */
  if (geometry_init() == 1) { finalize_process(); return 0; }
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */ 

  /* setup random numbers */
  read_input(rlx_var.read,input_filename);
  lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */

  if(strcmp(rlx_var.rlxd_start,"continue")==0 && rlx_var.rlxd_state[0]!='\0')
  {
    /*load saved state*/
    lprintf("MAIN",0,"Loading rlxd state from file [%s]\n",rlx_var.rlxd_state);
    read_ranlxd_state(rlx_var.rlxd_state);
  }

#ifdef GAUGE_SUN
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
#elif GAUGE_SON
  lprintf("MAIN",0,"Gauge group: SO(%d)\n",NG);
#else
  lprintf("MAIN",0,"Default gauge group: SU(%d)\n",NG);
#endif
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

  /* read input for llr update */
  read_input(llr_var.read,input_filename);

  lprintf("MAIN",0,"LLR nunber of mc steps per RM: %d\n",llr_var.nmc);
  lprintf("MAIN",0,"LLR nunber of therm steps per RM %d\n",llr_var.nth);
  lprintf("MAIN",0,"LLR Initial a %f\n",llr_var.starta);
  lprintf("MAIN",0,"LLR RM start value iteration %d\n",llr_var.it);
  lprintf("MAIN",0,"LLR S0 Central action %f\n",llr_var.S0);
  lprintf("MAIN",0,"LLR Delta S %f\n",llr_var.dS);
 
  /* read input for WF generation and measure */
  read_input(WF_var.read,input_filename);

  lprintf("MAIN",0,"WF make %s\n",WF_var.make);
  lprintf("MAIN",0,"WF tmax: %e\n",WF_var.tmax);
  lprintf("MAIN",0,"WF number of measures: %d\n",WF_var.nmeas);
  lprintf("MAIN",0,"WF time lapse between measures: %e\n",WF_var.tmax/WF_var.nmeas);
  lprintf("MAIN",0,"WF number of integration intervals per measure: %d\n",WF_var.nint);
  lprintf("MAIN",0,"WF number of integration intervals: %d\n",WF_var.nint*WF_var.nmeas);
  lprintf("MAIN",0,"WF integration step: %e\n",WF_var.tmax/(WF_var.nmeas*WF_var.nint));


  /* Init Monte Carlo */

  init_mc(&flow, input_filename);
  lprintf("MAIN",0,"MVM during HMC initialzation: %ld\n",getMVM());
  lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());

  
  init_robbinsmonro(llr_var.nmc,llr_var.nth,llr_var.starta,llr_var.it,llr_var.dS,llr_var.S0);

  WF_initialize();

  int j=0;
  for(j=0;j<flow.rmrestart;++j) {
    
    restart_robbinsmonro();  
    
    for (i=0;i<flow.therm;++i){
      struct timeval start, end, etime; /* //for trajectory timing */
      
      lprintf("MAIN",0,"Starting thermalization, thermalization steps = %d\n",flow.therm);
      lprintf("MAIN",0,"Thermalization #%d...\n",i);
      gettimeofday(&start,0);
      thermrobbinsmonro();
      gettimeofday(&end,0);
      timeval_subtract(&etime,&end,&start);
      lprintf("MAIN",0,"Thermalization sequence #%d: generated in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);    
      
      
#ifdef MEASURE_FORCELLR
      if (force_ave==NULL){
	force_ave = (double*) malloc(num_mon()*sizeof(double));
	force_max = (double*) malloc(num_mon()*sizeof(double));
	n_inv_iter = (int*) malloc(num_mon()*sizeof(int));
      }
      for (int k=0;k<num_mon();k++){
	force_ave[k]=0.0;
	force_max[k]=0.0;
	n_inv_iter[k]=0;
      }
#endif
      
      
    }
    
    
    
    lprintf("MAIN",0,"Thermalization done.\n");
    

    for(i=flow.start;i<flow.end;++i) {
      struct timeval start, end, etime; /* //for trajectory timing */
      lprintf("MAIN",0,"Trajectory #%d...\n",i);
      
      gettimeofday(&start,0);
      
#ifdef MEASURE_FORCELLR
      if (force_ave==NULL){
	force_ave = (double*) malloc(num_mon()*sizeof(double));
	force_max = (double*) malloc(num_mon()*sizeof(double));
	n_inv_iter = (int*) malloc(num_mon()*sizeof(int));
      }
      for (int k=0;k<num_mon();k++){
	force_ave[k]=0.0;
	force_max[k]=0.0;
	n_inv_iter[k]=0;
      }
#endif
      
      robbinsmonro();
      gettimeofday(&end,0);
      timeval_subtract(&etime,&end,&start);
      lprintf("MAIN",0,"Robbins Monro sequence #%d: generated in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);

      lprintf("MAIN",0,"Plaq a fixed %lf \n",avr_plaquette());    
      lprintf("MAIN",0,"<a_rho(%d,%d,%lf)>= %f\n",j,i,getS0(),get_llr_a());
      
      
    }
    
    
    
    lprintf("MAIN",0,"Robins Monro update done.\n");
    int k, n;
    double epsilon=WF_var.tmax/(WF_var.nmeas*WF_var.nint),t=0.;                                                                                                                                      
    double wf_E,wf_Esym,wf_TC;
    
    for(i=0;i<flow.obsnmeas;++i) {
      struct timeval start, end, etime; /* //for trajectory timing */
      
      llr_fixed_a_update();
      
      lprintf("MAIN",0,"Obs measure for fixed E=%f dE=%f a=%f T\n",getS0(),getdS(),get_llr_a());
      
      if((i%flow.meas_freq)==0) {                                                               
	/* plaquette */                                                                                                                                 
   	lprintf("MAIN",0,"Plaquette: %1.8e\n",avr_plaquette());                                                                                                            
        if(strcmp(WF_var.make,"true")==0) {
	  t=0.;	
	  for(n=-1;n<WF_var.nmeas;n++) {
	    if(n>-1) {
	      for(k=0;k<WF_var.nint;k++) {
		WilsonFlow3(u_gauge,epsilon);
		t+=epsilon;
	      }
	    }
	    wf_E=WF_E(u_gauge);
	    wf_Esym=WF_Esym(u_gauge);
	    wf_TC=WF_topo(u_gauge);
	    lprintf("WILSONFLOW",0,"WF (ncnfg,t,E,t2*E,Esym,t2*Esym,TC) = %d %e %e %e %e %e %e\n",i,t,wf_E,t*t*wf_E,wf_Esym,t*t*wf_Esym,wf_TC);
	  }
	}
	
	timeval_subtract(&etime,&end,&start);
	lprintf("MAIN",0,"LLR Update sequence for Obs #%d: generated in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
      }
    }
  }
  /* save final configuration */
  save_conf(&flow, (flow.obsnmeas + flow.end-flow.start)*flow.rmrestart);
  /* Only save state if we have a file to save to */
  if(rlx_var.rlxd_state[0]!='\0') {
    lprintf("MAIN",0,"Saving rlxd state to file %s\n",rlx_var.rlxd_state);
    write_ranlxd_state(rlx_var.rlxd_state);
  }
  gettimeofday(&endmain,0);
  timeval_subtract(&etimemain,&endmain,&startmain);
  
  lprintf("MAIN",0,"Total simulation time =[%ld sec %ld usec]\n",etimemain.tv_sec,etimemain.tv_usec);
  /* finalize Monte Carlo */
  end_mc();
  
  /* close communications */
  finalize_process();
  
  return 0;
  
}
