/**********************************************************************
 * 
 * hmm_f2ss.c
 * 
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 * Sept, 2001
 * Licensed under the GNU General Public License version 2 (June, 1991)
 * 
 * C functions for the R/qtl package
 *
 * Contains: init_f2ss, emit_f2ss, step_f2ss, nrec_f2ss1,
 *           nrec_f2ss2, est_map_f2ss, argmax_geno_f2ss
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the F2 intercross with 
 * sex-specific maps.
 *
 * Genotype codes:  1=AA; 2=AB; 3=BA; 4=BB
 * Phenotype codes: 0=missing; 1=AA; 2=AB; 3=BB; 4=not BB; 5=not AA
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "hmm_main.h"
#include "hmm_f2ss.h"

double init_f2ss(int true_gen)
{
  return(log(0.25));
}

double emit_f2ss(int obs_gen, int true_gen, double error_prob)
{
  switch(obs_gen) {
  case 0: return(0.0);
  case 1: 
    if(true_gen==1) return(log(1.0-error_prob));
    else return(log(error_prob/2.0));
  case 2:
    if(true_gen==2 || true_gen==3) return(log(1.0-error_prob));
    else return(log(error_prob/2.0));
  case 3:
    if(true_gen==4) return(log(1.0-error_prob));
    else return(log(error_prob/2.0));
  case 4: /* AA or AB (not BB) */
    if(true_gen != 4) return(log(1.0-error_prob/2.0));
    else return(log(error_prob/2.0));
  case 5: /* AB or BB (not AA) */
    if(true_gen != 1) return(log(1.0-error_prob/2.0));
    else return(log(error_prob/2.0));
  }
  return(log(-1.0)); /* shouldn't get here */
}
    
  
double step_f2ss(int gen1, int gen2, double rf1, double rf2) 
{
  switch(gen1) {
  case 1:
    switch(gen2) {
    case 1: return(log(1.0-rf1)+log(1.0-rf2));
    case 2: return(log(rf2)+log(1.0-rf1));
    case 3: return(log(rf1)+log(1.0-rf2));
    case 4: return(log(rf1)+log(rf2));
    }
  case 2:
    switch(gen2) {
    case 1: return(log(rf2)+log(1.0-rf1));
    case 2: return(log(1.0-rf1)+log(1.0-rf2));
    case 3: return(log(rf1)+log(rf2));
    case 4: return(log(rf1)+log(1.0-rf2));
    }
  case 3:
    switch(gen2) {
    case 4: return(log(rf2)+log(1.0-rf1));
    case 3: return(log(1.0-rf1)+log(1.0-rf2));
    case 2: return(log(rf1)+log(rf2));
    case 1: return(log(rf1)+log(1.0-rf2));
    }
  case 4:
    switch(gen2) {
    case 4: return(log(1.0-rf1)+log(1.0-rf2));
    case 3: return(log(rf2)+log(1.0-rf1));
    case 2: return(log(rf1)+log(1.0-rf2));
    case 1: return(log(rf1)+log(rf2));
    }
  }
  return(log(-1.0)); /* shouldn't get here */
}

double nrec_f2ss1(int gen1, int gen2)
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


double nrec_f2ss2(int gen1, int gen2)
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


void est_map_f2ss(int *n_ind, int *n_mar, int *geno, double *rf1,
		  double *rf2, double *error_prob, double *loglik,
		  int *maxit, double *tol, int *junk, int *prnt)
{
  est_map(*n_ind, *n_mar, 4, geno, rf1, rf2, *error_prob, 
	  init_f2ss, emit_f2ss, step_f2ss, nrec_f2ss1, nrec_f2ss2,
	  loglik, *maxit, *tol, 1, *prnt);
}


void argmax_geno_f2ss(int *n_ind, int *n_pos, int *geno, 
		      double *rf1, double *rf2, 
		      double *error_prob, int *argmax)
{		    
  argmax_geno(*n_ind, *n_pos, 4, geno, rf1, rf2, *error_prob,
	      argmax, init_f2ss, emit_f2ss, step_f2ss);
}

/* end of hmm_f2ss.c */
