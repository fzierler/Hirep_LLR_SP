/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
 *
 * File stopfile.c
 *
 * Allow for graceful stop if a "stopfile" is found
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "io.h"
#include "error.h"
#include "global.h"
#include "logger.h"
#include "observables.h"
#include "communications.h"
#include "utils.h"

/* One stopfile to be placed in directory Rep_0  */ 
int check_stopfile(char filename[])  {
    int foundfile=0;
    if(MPI_PID==0) {
        if(fopen(filename,"r")!=NULL) foundfile=1;
    }
#ifdef WITH_MPI
    MPI_Bcast(&foundfile,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
    return foundfile;
}
