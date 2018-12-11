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
#include "communications.h"
#include "observables.h"


#include <stdio.h>
#include <math.h>



void force_llr_0(double dt, suNg_av_field *force, void *vpar,double* action){
 force0( dt,force,  vpar);
 double beta = *((double*)vpar);
 *action = (1.0-avr_plaquette())*(beta*6.*GLB_VOLUME);
}

