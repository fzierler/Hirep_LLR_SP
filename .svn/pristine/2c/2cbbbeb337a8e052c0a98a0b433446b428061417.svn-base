/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File plaquette.c
*
* Routines for the average plaquette
*
*******************************************************************************/

#include "global.h"
#include "geometry.h"
#include "suN_repr_func.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "memory.h"
#include "random.h"
#include "dirac.h"
#include "representation.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"
#include "observables.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

void clover_F(suNg_algebra_vector *F, int ix, int mu, int nu) {
  int iy,iz,iw;
  suNg *v1,*v2,*v3,*v4,w1,w2,w3;
 
  _suNg_unit(w3);
  _suNg_mul(w3,-4.,w3);
 
  iy=iup(ix,mu);
  iz=iup(ix,nu);
   
  v1=pu_gauge(ix,mu);
  v2=pu_gauge(iy,nu);
  v3=pu_gauge(iz,mu);
  v4=pu_gauge(ix,nu);
   
  _suNg_times_suNg(w1,(*v1),(*v2));
  _suNg_times_suNg_dagger(w2,w1,(*v3));
  _suNg_times_suNg_dagger(w1,w2,(*v4));
#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    _suNg_mul(w1,plaq_weight[ix*16+nu*4+mu],w1);
  }
#endif
  _suNg_add_assign(w3,w1);
   
  iy=idn(ix,mu);
  iz=iup(iy,nu);
   
  v1=pu_gauge(ix,nu);
  v2=pu_gauge(iz,mu);
  v3=pu_gauge(iy,nu);
  v4=pu_gauge(iy,mu);
   
  _suNg_times_suNg_dagger(w1,(*v1),(*v2));
  _suNg_times_suNg_dagger(w2,w1,(*v3));
  _suNg_times_suNg(w1,w2,(*v4));
#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    _suNg_mul(w1,plaq_weight[iy*16+nu*4+mu],w1);
  }
#endif
  _suNg_add_assign(w3,w1);
   
  iy=idn(ix,mu);
  iz=idn(iy,nu);
  iw=idn(ix,nu);
   
  v1=pu_gauge(iy,mu);
  v2=pu_gauge(iz,nu);
  v3=pu_gauge(iz,mu);
  v4=pu_gauge(iw,nu);
   
  _suNg_times_suNg(w1,(*v2),(*v1));
  _suNg_dagger_times_suNg(w2,w1,(*v3));
  _suNg_times_suNg(w1,w2,(*v4));
#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    _suNg_mul(w1,plaq_weight[iz*16+nu*4+mu],w1);
  }
#endif
  _suNg_add_assign(w3,w1);
 
  iy=idn(ix,nu);
  iz=iup(iy,mu);
   
  v1=pu_gauge(iy,nu);
  v2=pu_gauge(iy,mu);
  v3=pu_gauge(iz,nu);
  v4=pu_gauge(ix,mu);
   
  _suNg_dagger_times_suNg(w1,(*v1),(*v2));
  _suNg_times_suNg(w2,w1,(*v3));
  _suNg_times_suNg_dagger(w1,w2,(*v4));
#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    _suNg_mul(w1,plaq_weight[iy*16+nu*4+mu],w1);
  }
#endif  
  _suNg_add_assign(w3,w1);
   
  _fund_algebra_project(*F,w3);
   
  _algebra_vector_mul_g(*F,1/4.,*F);
}


double EMT_E(int ix) {
  int mu,nu,iy;
  double ret;
   
  ret=24.*NG;
 
  for(mu=0;mu<4;mu++) for(nu=mu+1;nu<4;nu++) {
      ret-=plaq(ix, mu, nu);
      iy=idn(ix,mu);
      ret-=plaq(iy, mu, nu);
      iy=idn(iy,nu);
      ret-=plaq(iy, mu, nu);
      iy=idn(ix,nu);
      ret-=plaq(iy, mu, nu);
    }
   
  return ret/2.;
}



double EMT_Es(int ix) {
  int mu,nu,iy;
  double ret;
   
  ret=12.*NG;
 
  for(mu=1;mu<4;mu++) for(nu=mu+1;nu<4;nu++) {
      ret-=plaq(ix, mu, nu);
      iy=idn(ix,mu);
      ret-=plaq(iy, mu, nu);
      iy=idn(iy,nu);
      ret-=plaq(iy, mu, nu);
      iy=idn(ix,nu);
      ret-=plaq(iy, mu, nu);
    }
   
  return ret/2.;
}
 
 
 
double EMT_Et(int ix) {
  int nu,iy;
  double ret;
   
  ret=12.*NG;
 
  for(nu=1;nu<4;nu++) {
    ret-=plaq(ix, 0, nu);
    iy=idn(ix,0);
    ret-=plaq(iy, 0, nu);
    iy=idn(iy,nu);
    ret-=plaq(iy, 0, nu);
    iy=idn(ix,nu);
    ret-=plaq(iy, 0, nu);
  }
   
  return ret*2.;
}


double EMT_Ec(int ix) {
  int mu,nu;
  suNg_algebra_vector F;
  double p;
  double ret;
 
  ret=0.;
 
  for(mu=0;mu<4;mu++) for(nu=mu+1;nu<4;nu++) {
      clover_F(&F,ix,mu,nu);
      _algebra_vector_sqnorm_g(p,F);
      ret += p;
    }
 
  return ret*_FUND_NORM2;
}
 
 
 
double EMT_Ecs(int ix) {
  int mu,nu;
  suNg_algebra_vector F;
  double p;
  double ret;
 
  ret=0.;
 
  for(mu=1;mu<4;mu++) for(nu=mu+1;nu<4;nu++) {
      clover_F(&F,ix,mu,nu);
      _algebra_vector_sqnorm_g(p,F);
      ret += p;
    }
 
  return ret*_FUND_NORM2;
}
 
 
 
double EMT_Ect(int ix) {
  int nu;
  suNg_algebra_vector F;
  double p;
  double ret;
 
  ret=0.;
 
  for(nu=1;nu<4;nu++) {
    clover_F(&F,ix,0,nu);
    _algebra_vector_sqnorm_g(p,F);
    ret += p;
  }
 
  return ret*_FUND_NORM2;
}

double EMT_Y(int ix, int mu, int nu) {
  int sigma;
  suNg_algebra_vector F1,F2;
  double p;
  double ret;
 
  ret=0.;
 
  for(sigma=0;sigma<4;sigma++) {
    if(sigma==mu || sigma==nu) continue;
    clover_F(&F1,ix,mu,sigma);
    clover_F(&F2,ix,nu,sigma);
    _algebra_vector_prod_g(p,F1,F2);
    ret += p;
  }
   
  ret*=2.*_FUND_NORM2;
 
  if(mu==nu) ret-=EMT_Ec(ix);
 
  return ret;
}

double EMT_Y00(int ix) {
  int sigma;
  suNg_algebra_vector F1,F2;
  double p;
  double ret;
 
  ret=0.;
 
  for(sigma=1;sigma<4;sigma++) {
    clover_F(&F1,ix,0,sigma);
    clover_F(&F2,ix,0,sigma);
    _algebra_vector_prod_g(p,F1,F2);
    ret += p;
  }
   
  ret*=2.*_FUND_NORM2;
 
  ret-=EMT_Ec(ix);
 
  return ret;
}


double avg_EMT_Y00()
{
  double pa=0.;

  _PIECE_FOR(&glattice,ixp) {
    if(ixp==glattice.inner_master_pieces) {
      _OMP_PRAGMA( master )
	/* wait for gauge field to be transfered */
	complete_gf_sendrecv(u_gauge);
      _OMP_PRAGMA( barrier )
	}
    _SITE_FOR_SUM(&glattice,ixp,ix,pa) {
      pa+=EMT_Y00(ix);
    }
  }

  global_sum(&pa, 1);

  return pa/(double)(NG)/GLB_VOLUME;

}

double EMT_Y0k(int ix,int k) {
  int sigma;
  suNg_algebra_vector F1,F2;
  double p;
  double ret;
 
  ret=0.;
 
  for(sigma=1;sigma<4;sigma++) {
    clover_F(&F1,ix,0,sigma);
    clover_F(&F2,ix,k,sigma);
    _algebra_vector_prod_g(p,F1,F2);
    ret += p;
  }
   
  ret*=2.*_FUND_NORM2;
 
  //ret-=EMT_Ec(ix);
 
  return ret;
}

double avg_EMT_Y0k(int k)
{
  double pa=0.;

  _PIECE_FOR(&glattice,ixp) {
    if(ixp==glattice.inner_master_pieces) {
      _OMP_PRAGMA( master )
	/* wait for gauge field to be transfered */
	complete_gf_sendrecv(u_gauge);
      _OMP_PRAGMA( barrier )
	}
    _SITE_FOR_SUM(&glattice,ixp,ix,pa,k) {
      pa+=EMT_Y0k(ix,k);
    }
  }

  global_sum(&pa, 1);

  return pa/GLB_VOLUME;

}




