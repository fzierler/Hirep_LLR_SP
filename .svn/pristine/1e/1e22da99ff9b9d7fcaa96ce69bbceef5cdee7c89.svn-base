/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "representation.h"
#include "logger.h"
#include "observables.h"
#include "communications.h"

#include <stdio.h>
#include <math.h>



void force_obs_0pp(double dt, suNg_av_field *force, void *vpar){

  /* check input types */
  _TWO_SPINORS_MATCHING(u_gauge,force);
  force_obs_0pp_par *fpar=(force_obs_0pp_par*)vpar;

  int t0=fpar->t0;
  int t1=fpar->t1;
  

  int coord_zero[4];

  origin_coord(coord_zero);
  

  double p0=.0,p1=.0,coeff=.0;


  p0=spatial_wall_plaq(t0);
  p1=spatial_wall_plaq(t1);

  suNg s1,s2;
  suNg_algebra_vector f;
  double shift=fpar->shift;
  int g[4],i;



  
  if(coord_zero[0] <= t0 && coord_zero[0]+T > t0){ 

    coeff=p1/(p0*p1+shift);

 

    g[0]= t0 - coord_zero[0];
    for(g[3]=0;g[3]<Z;g[3]++)
      for(g[2]=0;g[2]<Y;g[2]++)
	for(g[1]=0;g[1]<X;g[1]++){
	  i=ipt(g[0],g[1],g[2],g[3]); 

	  for (int mu=1; mu<4; ++mu) {
	    spatial_staples(i,mu,&s1);
	    _suNg_times_suNg_dagger(s2,*_4FIELD_AT(u_gauge,i,mu),s1);
      
	    /* the projection itself takes the TA: proj(M) = proj(TA(M)) */
	    _fund_algebra_project(f,s2);
	    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,i,mu), dt*coeff, f);
      
	  }
	}
  }




  if(coord_zero[0] <= t1 && coord_zero[0]+T > t1){ 
  
    coeff=p0/(p0*p1+shift);


    g[0]= t1 - coord_zero[0];
    for(g[3]=0;g[3]<Z;g[3]++)
      for(g[2]=0;g[2]<Y;g[2]++)
	for(g[1]=0;g[1]<X;g[1]++){
	  i=ipt(g[0],g[1],g[2],g[3]); 

	  for (int mu=1; mu<4; ++mu) {
	    spatial_staples(i,mu,&s1);
	    _suNg_times_suNg_dagger(s2,*_4FIELD_AT(u_gauge,i,mu),s1);
      
	    /* the projection itself takes the TA: proj(M) = proj(TA(M)) */
	    _fund_algebra_project(f,s2);
	    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,i,mu), dt*coeff, f);
      
	  }
	}
  }
  
  apply_BCs_on_momentum_field(force);  
  }

