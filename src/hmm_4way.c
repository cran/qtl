/**********************************************************************
 * 
 * hmm_4way.c
 * 
 * copyright (c) 2001, 2002, Karl W Broman, Johns Hopkins University
 * 
 * last modified Mar, 2002
 * first written Feb, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 * 
 * C functions for the R/qtl package
 * 
 * Contains: init_4way, emit_4way, step_4way, nrec_4way, nrec_4way1,
 *           nrec_4way2, calc_genoprob_4way, sim_geno_4way, 
 *           est_map_4way, argmax_geno_4way, errorlod_4way, 
 *           calc_errorlod_4way, nrec2_4way, logprec_4way, est_rf_4way
 *           calc_pairprob_4way
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the "4-way" cross (autosomal data)
 *
 * Genotype codes:  0=AC, 1=AD, 2=BC, 3=BD
 * Phenotype codes: 0=missing, 1=AC, 2=BC, 3=AD, 4=BD,
 *                  5=(AC/AD=A), 6=(BC/BD=B),7=(AC/BC=C),8=(AD/BD=D)
 *                  9=(AC/BD), 10=(AD/BC)
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "hmm_main.h"

double init_4way(int true_gen)
{
  return(log(0.25));
}

double emit_4way(int obs_gen, int true_gen, double error_prob)
{
  switch(obs_gen) {
  case 0: return(0.0);
  case 1: 
    switch(true_gen) {
    case 1: return(log(1.0-error_prob));
    case 2: case 3: case 4: return(log(error_prob/3.0));
    }
  case 2: 
    switch(true_gen) {
    case 2: return(log(1.0-error_prob));
    case 1: case 3: case 4: return(log(error_prob/3.0));
    }
  case 3: 
    switch(true_gen) {
    case 3: return(log(1.0-error_prob));
    case 1: case 2: case 4: return(log(error_prob/3.0));
    }
  case 4: 
    switch(true_gen) {
    case 4: return(log(1.0-error_prob));
    case 1: case 2: case 3: return(log(error_prob/3.0));
    }
  case 5:
    switch(true_gen) {
    case 1: case 3: return(log(1.0-2.0*error_prob/3.0));
    case 2: case 4: return(log(2.0*error_prob/3.0));
    }
  case 6:
    switch(true_gen) {
    case 2: case 4: return(log(1.0-2.0*error_prob/3.0));
    case 1: case 3: return(log(2.0*error_prob/3.0));
    }
  case 7:
    switch(true_gen) {
    case 1: case 2: return(log(1.0-2.0*error_prob/3.0));
    case 3: case 4: return(log(2.0*error_prob/3.0));
    }
  case 8:
    switch(true_gen) {
    case 3: case 4: return(log(1.0-2.0*error_prob/3.0));
    case 1: case 2: return(log(2.0*error_prob/3.0));
    }
  case 9:
    switch(true_gen) {
    case 1: case 4: return(log(1.0-2.0*error_prob/3.0));
    case 2: case 3: return(log(2.0*error_prob/3.0));
    }
  case 10:
    switch(true_gen) {
    case 2: case 3: return(log(1.0-2.0*error_prob/3.0));
    case 1: case 4: return(log(2.0*error_prob/3.0));
    }
  }
  return(0.0); /* shouldn't get here */
}

double step_4way(int gen1, int gen2, double rf1, double rf2)
{
  switch(gen1) {
  case 1:
    switch(gen2) {
    case 1: return(log(1.0-rf1)+log(1.0-rf2));
    case 2: return(log(rf1)+log(1.0-rf2));
    case 3: return(log(1.0-rf1)+log(rf2));
    case 4: return(log(rf1)+log(rf2));
    }
  case 2:
    switch(gen2) {
    case 1: return(log(rf1)+log(1.0-rf2));
    case 2: return(log(1.0-rf1)+log(1.0-rf2));
    case 3: return(log(rf1)+log(rf2));
    case 4: return(log(1.0-rf1)+log(rf2));
    }
  case 3:
    switch(gen2) {
    case 1: return(log(1.0-rf1)+log(rf2));
    case 2: return(log(rf1)+log(rf2));
    case 3: return(log(1.0-rf1)+log(1.0-rf2));
    case 4: return(log(rf1)+log(1.0-rf2));
    }
  case 4:
    switch(gen2) {
    case 1: return(log(rf1)+log(rf2));
    case 2: return(log(1.0-rf1)+log(rf2));
    case 3: return(log(rf1)+log(1.0-rf2));
    case 4: return(log(1.0-rf1)+log(1.0-rf2));
    }
  }
  return(log(-1.0)); /* shouldn't get here */
}

double nrec_4way(int gen1, int gen2)
{
  switch(gen1) {
  case 1:
    switch(gen2) {
    case 1: return(0.0);
    case 2: case 3: return(0.5);
    case 4: return(1.0);
    }
  case 2:
    switch(gen2) {
    case 1: case 4: return(0.5);
    case 2: return(0.0);
    case 3: return(1.0);
    }
  case 3:
    switch(gen2) {
    case 1: case 4: return(0.5);
    case 3: return(0.0);
    case 2: return(1.0);
    }
  case 4:
    switch(gen2) {
    case 2: case 3: return(0.5);
    case 4: return(0.0);
    case 1: return(1.0);
    }
  }
  return(log(-1.0)); /* shouldn't get here */
}

double nrec_4way1(int gen1, int gen2)
{
  switch(gen1) {
  case 1: case 3:
    switch(gen2) {
    case 1: case 3: return(0.0);
    case 2: case 4: return(1.0);
    }
  case 2: case 4:
    switch(gen2) {
    case 1: case 3: return(1.0);
    case 2: case 4: return(0.0);
    }
  }
  return(log(-1.0)); /* shouldn't get here */
}

double nrec_4way2(int gen1, int gen2)
{
  switch(gen1) {
  case 1: case 2:
    switch(gen2) {
    case 1: case 2: return(0.0);
    case 3: case 4: return(1.0);
    }
  case 3: case 4:
    switch(gen2) {
    case 1: case 2: return(1.0);
    case 3: case 4: return(0.0);
    }
  }
  return(log(-1.0)); /* shouldn't get here */
}

void calc_genoprob_4way(int *n_ind, int *n_mar, int *geno, 
			double *rf1, double *rf2, double *error_prob, 
			double *genoprob) 
{
  calc_genoprob(*n_ind, *n_mar, 4, geno, rf1, rf2, *error_prob, genoprob,
		init_4way, emit_4way, step_4way);
}

void sim_geno_4way(int *n_ind, int *n_pos, int *n_draws, int *geno,
		   double *rf1, double *rf2, double *error_prob, int *draws)
{
  sim_geno(*n_ind, *n_pos, 4, *n_draws, geno, rf1, rf2, *error_prob,
	   draws, init_4way, emit_4way, step_4way);
}

void est_map_4way(int *n_ind, int *n_mar, int *geno, double *rf1, double *rf2,
		  double *error_prob, double *loglik, int *maxit, 
		  double *tol, int *sexsp, int *trace)
{
  est_map(*n_ind, *n_mar, 4, geno, rf1, rf2, *error_prob, 
	  init_4way, emit_4way, step_4way, nrec_4way1, nrec_4way2, 
	  loglik, *maxit, *tol, *sexsp, *trace);
}

void argmax_geno_4way(int *n_ind, int *n_pos, int *geno, 
		      double *rf1, double *rf2, 
		      double *error_prob, int *argmax)
{		    
  argmax_geno(*n_ind, *n_pos, 4, geno, rf1, rf2, *error_prob,
	      argmax, init_4way, emit_4way, step_4way);
}

double errorlod_4way(int obs, double *prob, double error_prob)
{
  double p=0.0;

  switch(obs) {
  case 0: p=1.0; break;
  case 1: case 2: case 3: case 4: p=prob[obs-1]; break;
  case 5: p=(prob[0]+prob[2]); break;
  case 6: p=(prob[1]+prob[3]); break;
  case 7: p=(prob[0]+prob[1]); break;
  case 8: p=(prob[2]+prob[3]); break;
  case 9: p=(prob[0]+prob[3]); break;
  case 10: p=(prob[1]+prob[2]); break;
  }
  
  p = (1.0-p)/p;

  if(obs>4) p *= (1.0-2.0*error_prob/3.0)/(2.0*error_prob/3.0);
  else p *= (1.0-error_prob)/error_prob;

  if(p < TOL) return(-12.0);
  else return(log10(p));
}

void calc_errorlod_4way(int *n_ind, int *n_mar, int *geno, 
		      double *error_prob, double *genoprob, 
		      double *errlod)
{
  calc_errorlod(*n_ind, *n_mar, 4, geno, *error_prob, genoprob,
		errlod, errorlod_4way);
}

double nrec2_4way(int obs1, int obs2, double rf)
{
  int temp;

  /* make obs1 <= obs2 */
  if(obs1 > obs2) {
    temp = obs2;
    obs2 = obs1;
    obs1 = temp;
  }

  switch(obs1) {
  case 1: 
    switch(obs2) {
    case 1: return(0.0);
    case 2: case 3: return(1.0);
    case 4: return(2.0);
    case 5: case 7: return(rf);
    case 6: case 8: return(1.0+rf);
    case 9: return(2*rf*rf/(rf*rf+(1.0-rf)*(1.0-rf)));
    case 10: return(1.0);
    }
  case 2:
    switch(obs2) {
    case 2: return(0.0);
    case 3: return(2.0);
    case 4: return(1.0);
    case 5: case 8: return(1.0+rf);
    case 6: case 7: return(rf);
    case 9: return(1.0);
    case 10: return(2*rf*rf/(rf*rf+(1.0-rf)*(1.0-rf)));
    }
  case 3:
    switch(obs2) {
    case 3: return(0.0);
    case 4: return(1.0);
    case 5: case 8: return(rf);
    case 6: case 7: return(1.0+rf);
    case 9: return(1.0);
    case 10: return(2*rf*rf/(rf*rf+(1.0-rf)*(1.0-rf)));
    }
  case 4: 
    switch(obs2) {
    case 4: return(0.0);
    case 5: case 7: return(1.0+rf);
    case 6: case 8: return(rf);
    case 9: return(2*rf*rf/(rf*rf+(1.0-rf)*(1.0-rf)));
    case 10: return(1.0);
    }
  case 5: 
    switch(obs2) {
    case 5: return(rf);
    case 6: return(1.0+rf);
    case 7: case 8: case 9: case 10: return(2.0*rf);
    }
  case 6: 
    switch(obs2) {
    case 6: return(rf);
    case 7: case 8: case 9: case 10: return(2.0*rf);
    }
  case 7:
    switch(obs2) {
    case 7: return(rf);
    case 8: return(1.0+rf);
    case 9: case 10: return(2.0*rf);
    }
  case 8:
    switch(obs2) {
    case 8: return(rf);
    case 9: case 10: return(2.0*rf);
    }
  case 9: case 10: 
    if(obs1==obs2) return(2*rf*rf/(rf*rf+(1.0-rf)*(1.0-rf)));
    else return(1.0);
  }
  return(log(-1.0)); /* shouldn't get here */
}

double logprec_4way(int obs1, int obs2, double rf)
{
  int temp;

  /* make obs1 <= obs2 */
  if(obs1 > obs2) {
    temp = obs2;
    obs2 = obs1;
    obs1 = temp;
  }

  switch(obs1) {
  case 1: 
    switch(obs2) {
    case 1: return(2.0*log(1.0-rf));
    case 2: case 3: return(log(rf)+log(1.0-rf));
    case 4: return(2.0*log(rf));
    case 5: case 7: return(log(1.0-rf));
    case 6: case 8: return(log(rf));
    case 9: return(log(rf*rf+(1.0-rf)*(1.0-rf)));
    case 10: return(log(2.0)+log(rf)+log(1.0-rf));
    }
  case 2:
    switch(obs2) {
    case 2: return(2.0*log(1.0-rf));
    case 3: return(2.0*log(rf));
    case 4: return(log(rf)+log(1.0-rf));
    case 5: case 8: return(log(rf));
    case 6: case 7: return(log(1.0-rf));
    case 9: return(log(2.0)+log(rf)+log(1.0-rf));
    case 10: return(log(rf*rf+(1.0-rf)*(1.0-rf)));
    }
  case 3:
    switch(obs2) {
    case 3: return(2.0*log(1.0-rf));
    case 4: return(log(rf)+log(1.0-rf));
    case 5: case 8: return(log(1.0-rf));
    case 6: case 7: return(log(rf));
    case 9: return(log(2.0)+log(rf)+log(1.0-rf));
    case 10: return(log(rf*rf+(1.0-rf)*(1.0-rf)));
    }
  case 4: 
    switch(obs2) {
    case 4: return(2.0*log(1.0-rf));
    case 5: case 7: return(log(rf));
    case 6: case 8: return(log(1.0-rf));
    case 9: return(log(rf*rf+(1.0-rf)*(1.0-rf)));
    case 10: return(log(2.0)+log(rf)+log(1.0-rf));
    }
  case 5: 
    switch(obs2) {
    case 5: return(log(1.0-rf));
    case 6: return(log(rf));
    case 7: case 8: case 9: case 10: return(0.0);
    }
  case 6: 
    switch(obs2) {
    case 6: return(log(1.0-rf));
    case 7: case 8: case 9: case 10: return(0.0);
    }
  case 7:
    switch(obs2) {
    case 7: return(log(1.0-rf));
    case 8: return(log(rf));
    case 9: case 10: return(0.0);
    }
  case 8:
    switch(obs2) {
    case 8: return(log(1.0-rf));
    case 9: case 10: return(0.0);
    }
  case 9: case 10: 
    if(obs1==obs2) return(log(rf*rf+(1.0-rf)*(1.0-rf)));
    else return(log(2.0)+log(rf)+log(1.0-rf));
  }
  return(log(-1.0)); /* shouldn't get here */
}

void est_rf_4way(int *n_ind, int *n_mar, int *geno, double *rf, 
	       int *maxit, double *tol)
{
  est_rf(*n_ind, *n_mar, geno, rf, nrec2_4way, logprec_4way, 
	 *maxit, *tol);
}

void calc_pairprob_4way(int *n_ind, int *n_mar, int *geno, 
			double *rf1, double *rf2, double *error_prob, 
			double *genoprob, double *pairprob) 
{
  calc_pairprob(*n_ind, *n_mar, 4, geno, rf1, rf2, *error_prob, genoprob,
		pairprob, init_4way, emit_4way, step_4way);
}

/* end of hmm_4way.c */
