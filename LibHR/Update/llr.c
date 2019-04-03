/*************************************************************************** \
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
 *
 * File update_llr.c
 *
 * Update programs
 *
 *******************************************************************************/
#include "suN.h"
#include "utils.h"
#include "global.h"
#include "logger.h"
#include "random.h"
#include "communications.h"
#include <math.h>
#include "update.h"
#include "observables.h"
#include <stdlib.h>

typedef struct
{
  int nrm, nth;
  int umb_RM_freq, umb_meas_freq, umb_therm_freq;
  int cfactor;
  int it;
  double starta;
  double a;
  double S0;
  double dS;
} llrparams;

static llrparams llrp;

//reset it to 0 and rhoa to initial value
void restart_robbinsmonro(int startit)
{
  llrp.it = startit;
  llrp.a = llrp.starta;
}

void init_robbinsmonro(int nrm, int nth, double starta, int it, double cfactor, int RMswap, int measswap, int thermswap, double dS, double S0)
{
  llrp.nrm = nrm;
  llrp.nth = nth;
  llrp.it = it;
  llrp.cfactor = cfactor;
  llrp.umb_RM_freq = RMswap;
  llrp.umb_meas_freq = measswap;
  llrp.umb_therm_freq = thermswap;
  llrp.starta = starta;
  llrp.dS = dS;
  llrp.S0 = S0;
}

double get_llr_a(void)
{
  return llrp.a;
}

double getS0(void)
{
  return llrp.S0;
}

double getdS(void)
{
  return llrp.dS;
}

void thermrobbinsmonro(void)
{
  static int step = 1;
  int no_metropolit_test = 0;
  if (step < 20)
    no_metropolit_test = 1;

  double S_llr, S_non_llr;
  update_llr_ghmc(&S_llr, &S_non_llr, no_metropolit_test);
#ifdef WITH_UMBRELLA
  if (step % (llrp.umb_therm_freq) == 0)
    umbrella_swap(&S_llr, &llrp.S0, &llrp.a, &llrp.dS, &llrp.starta);
#endif
  step++;
}

void llr_fixed_a_update(void)
{
  static int step = 1;

  double S_llr, S_non_llr;
  update_llr_ghmc(&S_llr, &S_non_llr, 0);

#ifdef WITH_UMBRELLA
  if (step % (llrp.umb_meas_freq) == 0)
    umbrella_swap(&S_llr, &llrp.S0, &llrp.a, &llrp.dS, &llrp.starta);
#endif
  step++;
}

void robbinsmonro(void)
{

  int rmstep;
  double S_llr;
  double S_non_llr;

  for (rmstep = 1; rmstep <= llrp.nth; rmstep++)
  {
    lprintf("ROBBINSMONRO", 30, "Inter Sequence Thermalization: %d\n", rmstep);
    update_llr_ghmc(&S_llr, &S_non_llr, 0);
#ifdef WITH_UMBRELLA
    if (rmstep % (llrp.umb_therm_freq) == 1)
      umbrella_swap(&S_llr, &llrp.S0, &llrp.a, &llrp.dS, &llrp.starta);
#endif
  }

  double S_llr_avr = 0., S_non_llr_avr = 0.;
  for (rmstep = 1; rmstep <= llrp.nrm; rmstep++)
  {
    update_llr_ghmc(&S_llr, &S_non_llr, 0);

    S_llr_avr += S_llr;
    S_non_llr_avr += S_non_llr;
#ifdef WITH_UMBRELLA
    if (rmstep % (llrp.umb_RM_freq) == 0)
      umbrella_swap(&S_llr, &llrp.S0, &llrp.a, &llrp.dS, &llrp.starta);
#endif
  }

  S_llr_avr /= (double)llrp.nrm;
  S_non_llr_avr /= (double)llrp.nrm;
  llrp.a += (S_llr_avr - llrp.S0) * llrp.cfactor / (llrp.dS * llrp.dS * llrp.it);
  lprintf("ROBBINSMONRO", 0, "RM Iteration: %d S_llr= %lf S_non_llr= %lf  a_llr= %lf \n", llrp.it, S_llr_avr, S_non_llr_avr, llrp.a);
  llrp.it++;
}

typedef struct
{
  double S_llr;
  double dS;
  double S0;
  double a;
  double starta;
  double deltaS;
  int rep;
  int repnext;

} reppar;

#ifdef WITH_UMBRELLA
static int compare_S0(const void *p, const void *q)
{
  //lprintf("llr:compare",0,"inside compare \n");
  reppar x = *(reppar *)p;
  reppar y = *(reppar *)q;

  if (x.S0 < y.S0)
    return -1; // Return -1 if you want ascending, 1 if you want descending order.
  else if (x.S0 > y.S0)
    return 1; // Return 1 if you want ascending, -1 if you want descending order.

  return 0;
}

static int compare_deltaS(const void *p, const void *q)
{
  //lprintf("llr:compare",0,"inside compare \n");
  reppar x = *(reppar *)p;
  reppar y = *(reppar *)q;

  if (x.deltaS < y.deltaS)
    return -1; // Return -1 if you want ascending, 1 if you want descending order.
  else if (x.deltaS > y.deltaS)
    return 1; // Return 1 if you want ascending, -1 if you want descending order.

  return 0;
}

void swap(double *data)
{
  reppar drep[N_REP];
  int toswap[N_REP];
  double rand, temp;
  int i, j;

  toswap[N_REP - 1] = 1;

  //sorting the replicas with respect of their energy  S0 so we can swap nearest neighbourgh

  for (i = 0; i < N_REP; i++)
  {

    drep[i].S_llr = data[5 * i];
    drep[i].S0 = data[5 * i + 1];
    drep[i].a = data[5 * i + 2];
    drep[i].dS = data[5 * i + 3];
    drep[i].starta = data[5 * i + 4];
    drep[i].rep = i;
  }

  qsort(drep, N_REP, sizeof(drep[0]), compare_S0);

  //sorting with respect of the difference in hamiltonian between r and r+1 where r is the replica

  for (i = 0; i < N_REP - 1; i++)
  {
    toswap[i] = 1;
    drep[i].repnext = drep[i + 1].rep;
    double S1 = drep[i + 1].a * drep[i].S_llr + drep[i].a * drep[i + 1].S_llr;
    double tmp1 = (drep[i].S_llr - drep[i + 1].S0) / drep[i + 1].dS;
    double tmp2 = (drep[i + 1].S_llr - drep[i].S0) / drep[i].dS;
    double S2 = (tmp1 * tmp1 + tmp2 * tmp2) / 2.;
    double S3 = drep[i + 1].a * drep[i + 1].S_llr + drep[i].a * drep[i].S_llr;
    tmp1 = (drep[i].S_llr - drep[i].S0) / drep[i].dS;
    tmp2 = (drep[i + 1].S_llr - drep[i + 1].S0) / drep[i + 1].dS;
    double S4 = (tmp1 * tmp1 + tmp2 * tmp2) / 2.;
    drep[i].deltaS = S1 + S2 - S3 - S4;  //djl
  }

  qsort(drep, N_REP - 1, sizeof(drep[0]), compare_deltaS);

  //swap of replicas

  for (i = 0; i < N_REP - 1; i++)
  {
    if (toswap[drep[i].rep] && toswap[drep[i].repnext])
    {
      toswap[drep[i].rep] = 0;
      toswap[drep[i].repnext] = 0;
      if (drep[i].deltaS < 0)
      {
        for (j = 0; j < 5; j++)
        {
          temp = data[5 * drep[i].rep + j];
          data[5 * drep[i].rep + j] = data[5 * drep[i].repnext + j];
          data[5 * drep[i].repnext + j] = temp;
        }
      }
      else
      {
        ranlxd(&rand, 1);
        if (rand < exp(-drep[i].deltaS))
        {
          for (j = 0; j < 5; j++)
          {
            temp = data[5 * drep[i].rep + j];
            data[5 * drep[i].rep + j] = data[5 * drep[i].repnext + j];
            data[5 * drep[i].repnext + j] = temp;
          }
        }
      }
    }
  }
}

void setreplica(double *data)
{
//  lprintf("llr:setreplica", 0, "Updating OLD LLR Param: S0 %lf,  a  %lf , dS %lf  \n", llrp.S0, llrp.a, llrp.dS);
  llrp.S0 = data[1];
  llrp.a = data[2];
  llrp.dS = data[3];
  llrp.starta = data[4];

//  lprintf("llr:setreplica", 0, "New LLR Param: S0 %lf,  a  %lf , dS %lf  \n", llrp.S0, llrp.a, llrp.dS);
}

#endif //WITH_UMBRELLA
