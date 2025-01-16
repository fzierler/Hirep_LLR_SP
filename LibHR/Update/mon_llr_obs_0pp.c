/***************************************************************************\
 * Copyright (c) 2015, Martin Hansen, Claudio Pica                        *
 * All rights reserved.                                                   *
 \***************************************************************************/

#include "global.h"
#include <math.h>
#include "update.h"
#include "logger.h"
#include "memory.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "observables.h"
#include <stdlib.h>

void llr_obs_0pp_init_traj(const struct _monomial *m) {
  /* empty */
}

void llr_obs_0pp_gaussian_pf(const struct _monomial *m) {
  /* empty */
}

void llr_obs_0pp_correct_pf(const struct _monomial *m) {
  /* empty */
}

void llr_obs_0pp_correct_la_pf(const struct _monomial *m) {
  /* empty */
}

const spinor_field* llr_obs_0pp_pseudofermion(const struct _monomial *m) {
  return NULL;
}

void llr_obs_0pp_add_local_action(const struct _monomial *m, scalar_field *loc_action) {
  
  int t0=((force_obs_0pp_par*)(m->force_par))->t0;
  int t1=((force_obs_0pp_par*)(m->force_par))->t1;
  double p0=spatial_wall_plaq(t0);
  double p1=spatial_wall_plaq(t1);
    
  if(PID==0){
    double shift=((force_obs_0pp_par*)(m->force_par))->shift;
    *_FIELD_AT(loc_action,ipt(1,1,1,1)) += log(p0*p1+shift); 

  }
}

void llr_obs_0pp_free(struct _monomial *m) {
  free(m->force_par);
  free(m);
}

struct _monomial* llr_obs_0pp_create(const monomial_data *data) {
  monomial *m = malloc(sizeof(*m));    
  // Copy data structure

  mon_obs_0pp_par *par = (mon_obs_0pp_par*)(data->par);

  m->data = *data;
  
  // Allocate memory for spinor field
  /* empty */
  
  // Setup force parameters
  /* empty */
  
  // Setup pointers to update functions
  m->free = &llr_obs_0pp_free;
  
  m->force_f = &force_obs_0pp;



  m->force_par = (void*)(&(par->fpar));
  par->fpar.t0=par->t0;  
  par->fpar.t1=par->t1;  
  par->fpar.shift=par->shift;
  

  m->pseudofermion = &llr_obs_0pp_pseudofermion;
  m->init_traj = &llr_obs_0pp_init_traj;
  m->gaussian_pf = &llr_obs_0pp_gaussian_pf;
  m->correct_pf = &llr_obs_0pp_correct_pf;
  m->correct_la_pf = &llr_obs_0pp_correct_la_pf;
  m->add_local_action = &llr_obs_0pp_add_local_action;
  m->add_llr_local_action = &llr_obs_0pp_add_local_action;
  
  return m;
}
