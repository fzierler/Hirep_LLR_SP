/***************************************************************************\
* Copyright (c) 2008, Antonio Rago                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/


#include "global.h"
#include "update.h"
#include "suN.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void force_llr_hmc(double dt, suNg_av_field *force, void *vpar,double* action){
force_hmc(dt,force, vpar);
spinor_field *Xs=get_hmc_spinor_Xs();
  *action=spinor_field_sqnorm_f(Xs);
}

