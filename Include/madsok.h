/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef MADSOK_H
#define MADSOK_H

#include "suN.h"
#include "observables.h"
//#include "update_u1.h"

void initmadsok(int a,int b,int c,int d);
////static double corr(double *energies,int dis,int size);
////static double madsok(double *energies,int lenght,int *cont);
void termal(double beta,int nhb,int nor);

#endif
