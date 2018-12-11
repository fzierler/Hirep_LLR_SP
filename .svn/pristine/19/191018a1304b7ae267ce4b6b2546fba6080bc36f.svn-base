/***************************************************************************\
 * Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
 * All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "memory.h"
#include "random.h"
#include "observables.h"
#include "dirac.h"
#include "representation.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"

/* State quantities for HMC */

static suNg_av_field *momenta=NULL;
static suNg_field *u_gauge_old=NULL;
static scalar_field *la=NULL; /* local action field for Metropolis test */
static scalar_field *la_llr=NULL; /* local action field for llr energies in Metropolis test */
static ghmc_par update_par;
static int init=0;
/*Adaptive thermalization for LLR*/
#ifdef ADAPTIVE

typedef struct{
  int stepsmeasure;
  int stepsadaptive;
  double deltaS;
  double deltaSold;
  double mulcoeff;
  double targetSmin; 
  double targetSmax;
 

} adaptive_par;

static adaptive_par adp;

static inline double abs_val(double x){
  if(x>=0)
    return x;
  else
    return -x;
}


void total_llr_action( double * S_llr){
  _MASTER_FOR(&glattice,i) {
    *_FIELD_AT(la_llr,i)=0.;
  }

  int nmon=num_mon();

  for (int i=0;i<nmon;++i) {
    const monomial *m=mon_n(i);
    if(is_llr[m->data.type])
      m->add_llr_local_action(m,la_llr);
  }
  *S_llr=0.0;

  _MASTER_FOR(&glattice,i) {
    *S_llr+=*_FIELD_AT(la_llr,i);
  }

  global_sum(S_llr,1);
}


void setstep(){
  if(adp.stepsadaptive>adp.stepsmeasure)
    update_par.integrator->nsteps=adp.stepsadaptive;
  else
    update_par.integrator->nsteps=adp.stepsmeasure;
}

static void adapt(int therm){

  if(therm==0){
    if(adp.stepsadaptive>adp.stepsmeasure)
      update_par.integrator->nsteps=adp.stepsadaptive;
    else
      update_par.integrator->nsteps=adp.stepsmeasure;
  }
  else{
    if (abs_val(adp.deltaS) < adp.targetSmin){
      
      if (abs_val(adp.deltaSold) < adp.targetSmin){
	if(adp.mulcoeff<10.)
	adp.mulcoeff*=2.;
      }else
	if(adp.mulcoeff>0.02)
	  adp.mulcoeff/=2.;
      
      adp.stepsadaptive/=(1+adp.mulcoeff);
      if (adp.stepsadaptive<4)
	adp.stepsadaptive=4;
      
      
    }else if(abs_val(adp.deltaS)>adp.targetSmax){
      
      if (abs_val(adp.deltaSold) > adp.targetSmin){
	if(adp.mulcoeff<10.)
	adp.mulcoeff*=2;
      }else
	  if(adp.mulcoeff>0.02)
	    adp.mulcoeff/=2;
      
      adp.stepsadaptive*=(1+adp.mulcoeff);
	
    }else{
      if(adp.mulcoeff>0.02)
	adp.mulcoeff/=2;
      
      
    }
    update_par.integrator->nsteps=adp.stepsadaptive;
  }
  
}



static void init_adaptive(void){
  adp.stepsmeasure=update_par.integrator->nsteps;
  adp.stepsadaptive=update_par.integrator->nsteps;
  adp.deltaSold=1.;
  adp.mulcoeff=1.;
  adp.targetSmin=0.1;
  adp.targetSmax=0.35;
}

#endif

void init_ghmc(ghmc_par *par){
  
  if (init) {
    /* already initialized */
    lprintf("GHMC",0,"WARNING: GHMC already initialized!\nWARNING: Ignoring call to init_ghmc.\n");
    return;
  }
	
  lprintf("GHMC",0,"Initializing...\n");
	  
  /* allocate space for the backup copy of gfield */
  if(u_gauge_old==NULL) u_gauge_old=alloc_gfield(&glattice);
  suNg_field_copy(u_gauge_old,u_gauge);
  
  /* allocate momenta */
  if(momenta==NULL) momenta = alloc_avfield(&glattice);
  
  /* allocate pseudofermions */
  /* we allocate one more pseudofermion for the computation
   * of the final action density
   */

  /* allocate memory for the local action */
  /* NOTE: should this be moved into local_action.c ? */
  if(la==NULL){
    la=alloc_sfield(2, &glattice);  
    /* la_llr=la+glattice.gsize_spinor; */
    la_llr=alloc_sfield(2, &glattice);
  }
  /* represent gauge field */
  represent_gauge_field();
  
  /* copy update parameters */
  update_par = *par;
#ifdef ADAPTIVE
  init_adaptive();
#endif

//#ifdef ROTATED_SF
//  hmc_action_par.SF_ct = _update_par.SF_ct;
//#endif
  init = 1;
  
  lprintf("HMC",0,"Initialization done.\n");
}

void free_ghmc()
{
	if (!init)
	{
		/* not initialized */
		lprintf("HMC",0,"WARNING: HMC not initialized!\nWARNNG: Ignoring call to free_hmc.\n");
		return;
	}

	/* free momenta */
	if(u_gauge_old!=NULL) free_gfield(u_gauge_old); u_gauge_old=NULL;
	if(momenta!=NULL) free_avfield(momenta); momenta=NULL;
	if(la!=NULL) free_sfield(la); la=NULL;
  
	/*Free integrator */
	integrator_par *ip = update_par.integrator;
	while(ip != NULL)
	{
		update_par.integrator=ip->next;
		free(ip->mon_list);
		free(ip);
		ip = update_par.integrator;
	}
	update_par.integrator = NULL;

	//free_force_hmc();
	init = 0;
	lprintf("HMC",0,"Memory deallocated.\n");
}

int update_ghmc()
{
  double deltaH;

  if(!init)
    {
      /* not initialized */
      lprintf("HMC",0,"WARNING: GHMC not initialized!\nWARNNG: Ignoring call to update_ghmc.\n");
      return -1;
    }

  /* init monomials */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->init_traj(m);
  }
  
  /* generate new momenta */
  lprintf("HMC",30,"Generating gaussian momenta and pseudofermions...\n");
  gaussian_momenta(momenta);

  /* generate new pseudofermions */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->gaussian_pf(m);
  }

  /* compute starting action */
  lprintf("HMC",30,"Computing action density...\n");
  local_hmc_action(NEW, la, momenta);

  /* correct pseudofermion distribution */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->correct_pf(m);
  }

  /* integrate molecular dynamics */
  lprintf("HMC",30,"MD integration...\n");

  update_par.integrator->integrator(momenta,update_par.tlen,update_par.integrator);

  /* project and represent gauge field */
  project_gauge_field();
  represent_gauge_field();

  /* compute new action */
  lprintf("HMC",30,"Computing new action density...\n");
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->correct_la_pf(m);
  }

  local_hmc_action(DELTA, la, momenta);



  /* Metropolis test */
  deltaH = 0.0;
  _MASTER_FOR_SUM(la->type,i,deltaH) {
    deltaH += *_FIELD_AT(la,i);
  }

  global_sum(&deltaH, 1);

  lprintf("HMC",10,"[DeltaS = %1.8e][exp(-DS) = %1.8e]\n",deltaH,exp(-deltaH));

  if(deltaH < 0) {
    suNg_field_copy(u_gauge_old,u_gauge);
  } else {
    double r;
    if(PID == 0) {
      ranlxd(&r,1);
      if(r < exp(-deltaH)) {
        r = 1.0;
      } else {
        r = -1.0;
      }
    }

    bcast(&r, 1);

    if(r > 0) {
      suNg_field_copy(u_gauge_old,u_gauge);
    } else {
      lprintf("HMC",10,"Configuration rejected.\n");
      suNg_field_copy(u_gauge,u_gauge_old);
      start_gf_sendrecv(u_gauge); /* this may not be needed if we always guarantee that we copy also the buffers */
      represent_gauge_field();
      return 0;
    }
  }

  lprintf("HMC",10,"Configuration accepted.\n");
  return 1;
}


#ifdef MEASURE_FORCEHMC
/*Functions to check forces */
void corret_pf_dist_hmc(){
    /* init monomials */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->init_traj(m);
  }
  
  /* generate new momenta */
  lprintf("HMC",30,"Generating gaussian momenta and pseudofermions...\n");
  gaussian_momenta(momenta);

  /* generate new pseudofermions */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->gaussian_pf(m);
  }

  /* compute starting action */
  lprintf("HMC",30,"Computing action density...\n");
  local_hmc_action(NEW, la, momenta);

  /* correct pseudofermion distribution */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->correct_pf(m);
  }
}

void calc_one_force(int n_force){
  integrator_par *ip=update_par.integrator;
  for (;;){
    error(ip==NULL,1,"calc_one_force","Error in force index\n");
    for(int n = 0; n < ip->nmon; n++){
      const monomial *m=ip->mon_list[n];
      if (m->data.id==n_force){
        m->force_f(1,momenta,m->force_par);
        return;
      }
    }
    ip=ip->next;
  }
}





#endif



int update_llr_ghmc(double  *ret_S_llr, double  *ret_S_non_llr,int therm)
{

  lprintf("LLR HMC",0," Starting update...\n");
  double deltaH;
  double tmp[6];
  double *Snew=tmp, *Sold=tmp+1, *momold=tmp+2, *momnew=tmp+3,*Snew_llr=tmp+4,*Sold_llr=tmp+5;

  if(!init)
    {
      /* not initialized */
      lprintf("LLR HMC",0,"WARNING: GHMC not initialized!\nWARNNG: Ignoring call to update_llr_ghmc.\n");
      return -1;
    }
  
  /* init monomials */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->init_traj(m);
  }
  
  /* generate new momenta */
  lprintf("LLR HMC",30,"Generating gaussian momenta and pseudofermions...\n");
  gaussian_momenta(momenta);

  /* generate new pseudofermions */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->gaussian_pf(m);
  }

  /* compute starting action */
  lprintf("LLR HMC",30,"Computing initial action density...\n");
  local_llr_hmc_action(Sold, Sold_llr, momold, la, la_llr, momenta);

  /* correct pseudofermion distribution */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->correct_pf(m);
  }


  /* integrate molecular dynamics */
  lprintf("LLR HMC",10,"MD integration...\n");
#ifdef ADAPTIVE
  lprintf("LLR HMC",10,"Adaptive stepsize definition enabled \n");
  adapt(therm);
#endif
  update_par.integrator->integrator(momenta,update_par.tlen,update_par.integrator);

  /* project and represent gauge field */
  project_gauge_field();
  represent_gauge_field();

  /* compute new action */
  lprintf("LLR HMC",30,"Computing new action density...\n");
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->correct_la_pf(m);
  }

  local_llr_hmc_action(Snew, Snew_llr, momnew, la, la_llr, momenta);


  /* Metropolis test */

  /*this is horrible!!! this is the actual global sum of Snew Snew_llr Sold Sold_llr momnew momold*/
  global_sum(Snew,6);

  deltaH=*Snew+(*Snew_llr - getS0())*(*Snew_llr - getS0())/(2.*getdS()*getdS()) + *momnew - get_llr_a()*(*Snew_llr) - *Sold -(*Sold_llr - getS0())*(*Sold_llr - getS0())/(2.*getdS()*getdS()) + get_llr_a()*(*Sold_llr) -*momold;


#ifdef ADAPTIVE
  adp.deltaSold=adp.deltaS;
  adp.deltaS=deltaH;
#endif
  lprintf("LLR HMC",30,"[Snew non llr = %1.8e ][Sold non llr = %1.8e ][Snew llr = %1.8e ][Sold llr = %1.8e ]\n",*Snew,*Sold,*Snew_llr,*Sold_llr);
  lprintf("LLR HMC",30,"[momnew = %1.8e ][momold = %1.8e ]\n",*momnew,*momold);
  lprintf("LLR HMC",30,"[DeltaS = %1.8e ][exp(-DS) = %1.8e ]\n",deltaH,exp(-deltaH));

 
  if(deltaH < 0) {
    suNg_field_copy(u_gauge_old,u_gauge);
  } else {
    double r;
    if(PID == 0) {
      ranlxd(&r,1);
      if(r < exp(-deltaH)) {
        r = 1.0;
      } else {
        r = -1.0;
      }
    }

    bcast(&r, 1);

    if(r > 0) {
      suNg_field_copy(u_gauge_old,u_gauge);
    } else {
      lprintf("LLR HMC",0,"Configuration rejected.\n");
      suNg_field_copy(u_gauge,u_gauge_old);
      start_gf_sendrecv(u_gauge); /* this may not be needed if we always guarantee that we copy also the buffers */
      represent_gauge_field();
      *ret_S_llr=*Sold_llr;
      *ret_S_non_llr=*Sold;
      return 0;
    }
  }

  lprintf("LLR HMC",10,"Configuration accepted.\n");
  *ret_S_llr=*Snew_llr;
  *ret_S_non_llr=*Snew;

  return 1;
}


