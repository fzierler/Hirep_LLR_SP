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
#include <unistd.h>
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




#define init_input_llr(varname)                                         \
  {                                                                     \
    .read={                                                             \
      {"make llr iterations", "llr:make = %s", STRING_T, &((varname).make)}, \
      {"Number of MC steps per RM iteration ", "llr:nmc = %d", INT_T, &((varname).nmc)}, \
      {"Number of MC therm steps per RM iteration", "llr:nth = %d", INT_T, &((varname).nth)}, \
      {"Initial a", "llr:starta = %lf", DOUBLE_T, &((varname).starta)}, \
      {"Robbins Monro startint iteration", "llr:it = %d", INT_T, &((varname).it)}, \
      {"Cental action", "llr:S0 = %lf", DOUBLE_T, &((varname).S0)},     \
      {"Delta S", "llr:dS = %lf", DOUBLE_T, &((varname).dS)},           \
      {NULL, NULL, 0, NULL}                                             \
    }                                                                   \
  }
 
input_llr llr_var=init_input_llr(llr_var);

/*WF parameters*/
typedef struct _input_WF {
  char make[256];
  double tmax;
  int nmeas,nint;
  double eps;
  double delta;
  int def_glueball;

  /* for the reading function */
  input_record_t read[8];

} input_WF;

#define init_input_WF(varname)						\
  {									\
    .read={								\
      {"make WF", "WF:make = %s", STRING_T, (varname).make},		\
      {"WF max integration time", "WF:tmax = %lf", DOUBLE_T, &((varname).tmax)}, \
      {"WF number of measures", "WF:nmeas = %d", DOUBLE_T, &((varname).nmeas)}, \
      {"WF number of integration steps between measures", "WF:nint = %d", INT_T, &((varname).nint)}, \
      {"WF initial epsilon", "WF:eps = %lf", DOUBLE_T, &((varname).eps)}, \
      {"WF delta", "WF:delta = %lf", DOUBLE_T, &((varname).delta)},	\
      {"enable glueball", "WF:def_glueball = %d",INT_T, &(varname).def_glueball}, \
      {NULL, NULL,0,NULL}						\
    }									\
  }

input_WF WF_var = init_input_WF(WF_var);
suNg_field* wf_gauge=NULL;


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
  int i,j;
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
  gethostname(sbuf,128);
  lprintf("MAIN",0,"Hostname: %s\n",sbuf);
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
 


  /* Init Monte Carlo */

  init_mc(&flow, input_filename);
  
  lprintf("MAIN",0,"Initial Thermalization steps: %d\n",flow.therm);
  lprintf("MAIN",0,"Inter RM Thermalization steps: %d\n",flow.interrm_therm);
  lprintf("MAIN",0,"RM steps: %d\n",flow.end-flow.start);
  lprintf("MAIN",0,"RM start value: %d\n",flow.rmstart);
  lprintf("MAIN",0,"RM global restart: %d\n",flow.rmrestart);
  lprintf("MAIN",0,"Observable measurements start value: %d\n",flow.obsmeas_start);
  lprintf("MAIN",0,"Observable measurements after RM: %d\n",flow.obsnmeas);
  
  /* read input for llr update */
  read_input(llr_var.read,input_filename);

  lprintf("MAIN",0,"LLR number of mc steps per RM: %d\n",llr_var.nmc);
  lprintf("MAIN",0,"LLR number of therm steps per RM %d\n",llr_var.nth);
  lprintf("MAIN",0,"LLR Initial a %f\n",llr_var.starta);
  lprintf("MAIN",0,"LLR RM start value iteration %d\n",llr_var.it);
  lprintf("MAIN",0,"LLR S0 Central action %f\n",llr_var.S0);
  lprintf("MAIN",0,"LLR Delta S %f\n",llr_var.dS);
 
  /* read input for WF generation and measure */
  read_input(WF_var.read,input_filename);

  lprintf("MAIN",0,"WF make %s\n",WF_var.make);
  lprintf("MAIN",0,"WF tmax: %e\n",WF_var.tmax);
  lprintf("MAIN",0,"WF number of measurements: %d\n",WF_var.nmeas);
  lprintf("MAIN",0,"WF time lapse between measurements: %e\n",WF_var.tmax/WF_var.nmeas);
  lprintf("MAIN",0,"WF number of integration intervals per measurement: %d\n",WF_var.nint);
  lprintf("MAIN",0,"WF number of integration intervals: %d\n",WF_var.nint*WF_var.nmeas);
  lprintf("MAIN",0,"WF integration step: %e\n",WF_var.tmax/(WF_var.nmeas*WF_var.nint));


  WF_initialize();
  init_robbinsmonro(llr_var.nmc,llr_var.nth,llr_var.starta,llr_var.it,flow.umb_freq,llr_var.dS,llr_var.S0);

  lprintf("MAIN",0,"MVM during HMC initialzation: %ld\n",getMVM());
  lprintf("MAIN",0,"Initial Plaquette: %1.8e\n",avr_plaquette());
  lprintf("MAIN",0,"Start Initial Thermalization ------------------------------------------\n");


  for (i=1;i<=flow.therm;++i){
    struct timeval start, end, etime; /* //for trajectory timing */
    lprintf("MAIN",0,"Initial Thermalization step #%d/%d\n",i,flow.therm);
    gettimeofday(&start,0);
    thermrobbinsmonro();
    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Initial Thermalization Plaquette: %1.8e\n",avr_plaquette());
    lprintf("MAIN",0,"Initial Thermalization step #%d/%d: generated in [%ld sec %ld usec]\n",i,flow.therm,etime.tv_sec,etime.tv_usec);    
  }
  lprintf("MAIN",0,"End Initial Thermalization ------------------------------------------\n");

  
  for(j=flow.rmrestart;j<flow.rmrestart;++j) {
    
    lprintf("MAIN",0,"Start RM%d Inter Restart Thermalization ------------------------------------------\n",j);
    lprintf("MAIN",0,"-------------------------------------------------------------------\n");
    restart_robbinsmonro(llr_var.it);  
    
    for (i=1;i<=flow.interrm_therm;++i){
      struct timeval start, end, etime; /* //for trajectory timing */
      lprintf("MAIN",0,"RM%d Inter Restart Thermalization step #%d/%d\n",j,i,flow.interrm_therm);
      gettimeofday(&start,0);
      thermrobbinsmonro();
      gettimeofday(&end,0);
      timeval_subtract(&etime,&end,&start);
      lprintf("MAIN",0,"RM%d Inter Restart Thermalization Plaquette: %lf\n",j,avr_plaquette());
      lprintf("MAIN",0,"RM%d Inter Restart Thermalization step #%d/%d: generated in [%ld sec %ld usec]\n",j,i,flow.interrm_therm,etime.tv_sec,etime.tv_usec);    
      
      
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
    lprintf("MAIN",0,"End RM%d Inter Restart Thermalization ------------------------------------------\n",j);
    
   
    lprintf("MAIN",0,"Start RM%d Phase ------------------------------------------\n",j);
   
    lprintf("MAIN",0,"<a_rho(%d,%d,%lf)>= %f\n",j,flow.start-1,getS0(),get_llr_a()); // restart will have set a to starta

    for(i=flow.start;i<flow.end;++i) {
      struct timeval start, end, etime; /* //for trajectory timing */
      lprintf("MAIN",0,"RM%d sequence #%d/%d\n",j,i,flow.end-flow.start);
      
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
      lprintf("MAIN",0,"RM%d sequence %d Plaquette: %lf \n",j,i,avr_plaquette());    
      lprintf("MAIN",0,"<a_rho(%d,%d,%lf)>= %f\n",j,i,getS0(),get_llr_a());
      lprintf("MAIN",0,"RM%d sequence #%d/%d: generated in [%ld sec %ld usec]\n",j,i,flow.end-flow.start,etime.tv_sec,etime.tv_usec);
 
    }
    
    lprintf("MAIN",0,"RM%d Final S0= %lf a_llr= %f\n",j,getS0(),get_llr_a());
    lprintf("MAIN",0,"End RM%d Phase ------------------------------------------\n",j);
    lprintf("MAIN",0,"Start RM%d Measurement Phase ------------------------------------------\n",j);
    
    for(i=flow.obsmeas_start;i<flow.obsnmeas;++i) {
      struct timeval start, end, etime; /* //for trajectory timing */

      lprintf("MAIN",0,"RM%d Trajectory #%d of %d\n",j,i,flow.obsnmeas);
        
      gettimeofday(&start,0);
    
      llr_fixed_a_update();
          
      gettimeofday(&end,0);
      timeval_subtract(&etime,&end,&start);
      lprintf("MAIN",0,"RM%d Trajectory #%d: generated in [%ld sec %ld usec]\n",j,i,etime.tv_sec,etime.tv_usec);
      
      lprintf("MAIN",0,"RM%d Measurement %d for fixed S0= %f dS= %f a_llr= %f\n",j,i,getS0(),getdS(),get_llr_a());
      
      if((i%flow.meas_freq)==0) {                                                               
	    /* plaquette */                                                                                                                                 
   	    lprintf("MAIN",0,"RM%d Measurement %d Plaquette: %1.8e\n",j,i,avr_plaquette());
        
        if(strcmp(WF_var.make,"true")==0) {
          double E, Esym, TC;
          int k;
          double epsilon=WF_var.eps;
          double t=0.;
          double dt = (double)WF_var.tmax/(double)WF_var.nmeas;
          
          if(wf_gauge==NULL) wf_gauge=alloc_gfield(&glattice);
          
          
          E=WF_E(u_gauge);
          Esym=WF_Esym(u_gauge);
          TC=WF_topo(u_gauge);
          lprintf("WILSONFLOW",0,"WF (ncnfg,t,E,t2*E,Esym,t2*Esym,TC) = %d %e %e %e %e %e %e\n",i,t,E,t*t*E,Esym,t*t*Esym,TC);
          
          suNg_field_copy(wf_gauge,u_gauge);
          
          k=1;	
          double epsilon_new=0;
          while (t < WF_var.tmax)
            {	
              if (t+epsilon > (double)k*dt)
                epsilon = (double)k*dt - t; 
              
              epsilon_new=WilsonFlow3_adaptative(wf_gauge,epsilon,WF_var.delta);
              
              if ( fabs(epsilon_new+1.) > 1e-7) 
                t=t+epsilon;
              
              if ( fabs(t - (double)k*dt ) < 1e-7 ) {
                k=k+1;
                E=WF_E(wf_gauge);
                Esym=WF_Esym(wf_gauge);
                TC=WF_topo(wf_gauge);
                lprintf("WILSONFLOW",0,"WF (ncnfg,t,E,t2*E,Esym,t2*Esym,TC) = %d %e %e %e %e %e %e\n",i,t,E,t*t*E,Esym,t*t*Esym,TC);
              }
              if (fabs(epsilon_new + 1.) > 1e-7) epsilon=epsilon_new;	
              if (fabs(epsilon_new +1.) < 1e-7 ) epsilon=epsilon/2;	
              
            }
          
        }
        gettimeofday(&end,0);
        timeval_subtract(&etime,&end,&start);
        lprintf("MAIN",0,"RM%d Measurements on trajectory %d performed in [%ld sec %ld usec]\n",j,i,etime.tv_sec,etime.tv_usec);
  
        if((i%flow.save_freq)==0) {
          //save_conf(&flow, flow.obsnmeas*j + i);
          save_conf(&flow, j, i);
          /* Only save state if we have a file to save to */
          if(rlx_var.rlxd_state[0]!='\0') {
            lprintf("MAIN",0,"Saving rlxd state to file %s\n",rlx_var.rlxd_state);
          write_ranlxd_state(rlx_var.rlxd_state);
          } 
        }
      }
    }


    lprintf("MAIN",0,"End RM%d Measurement Phase ------------------------------------------\n",j);
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
