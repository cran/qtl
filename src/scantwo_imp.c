/**********************************************************************
 *
 * scantwo_imp.c
 *
 * copyright (c) 2001-2, Karl W Broman, Johns Hopkins University
 *                     and Hao Wu, The Jackson Lab
 *
 * This file was written by Hao Wu with modifications by 
 * Karl Broman.
 *
 * last modified Oct, 2002 
 * first written Nov, 2001 
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a 2-dimensional genome scan 
 * with a 2-QTL model by imputation.
 *
 * Contains: R_scantwo_imp, scantwo_imp, altRss2
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include "util.h"
#include "scantwo_imp.h"
#include "scanone_imp.h"

#define TOL 1.0e-10

/**********************************************************************
 *
 * R_scantwo_imp
 *
 * Wrapper for call from R; reorganizes genotype prob, additive and 
 * interactive covariates and result matrix. Then calls scantwo_imp.
 *
 **********************************************************************/

void R_scantwo_imp(int *n_ind, int *same_chr, int *n_pos1, int *n_pos2, 
		   int *n_gen1, int *n_gen2, int *n_draws, int *draws1, 
		   int *draws2, double *addcov, int *n_addcov, 
		   double *intcov, int *n_intcov, double *pheno, 
		   double *weights, double *result)
{
  int ***Draws1, ***Draws2;
  double **Addcov, **Intcov;

  /* reorganize draws */
  reorg_draws(*n_ind, *n_pos1, *n_draws, draws1, &Draws1);
  if(!(*same_chr)) reorg_draws(*n_ind, *n_pos2, *n_draws, draws2, &Draws2);

  /* reorganize addcov and intcov (if they are not empty) */
  /* currently reorg_geno function is used to reorganized the data */
  if(*n_addcov != 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov != 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  /* call the engine function scantwo_imp */
  scantwo_imp(*n_ind, *same_chr, *n_pos1, *n_pos2, *n_gen1, *n_gen2, 
	      *n_draws, Draws1, Draws2, Addcov, *n_addcov, 
	      Intcov, *n_intcov, pheno, weights, result);
}

/**********************************************************************
 * 
 * scantwo_imp
 *
 * Performs genotype pair scan using the pseudomarker algorithm 
 * (imputation) method of Sen and Churchill (2001).
 * 
 * n_ind        Number of individuals
 *
 * same_chr     If = 1, work only with Draws1 and do 2-QTL model with
 *              QTLs on the same chromosome.
 *
 * chr2         Chromesome id 2
 *
 * n_pos1       Number of marker positions in chromesome 1
 *
 * n_pos2       Number of marker positions in chromesome 2
 *
 * n_gen1       Number of different genotypes on chr 1
 *
 * n_gen2       Number of different genotypes on chr 2
 *
 * n_draws      Number of impiutations
 *
 * Draws1       Array of genotype imputations in chromesome 1, 
 *              indexed as Draws1[repl][mar][ind]
 * 
 * Draws2       Array of genotype imputations in chromesome 2, 
 *              indexed as Draws2[repl][mar][ind]
 *
 * addcov	Additive covariates matrix, addcov[mar][ind]
 *
 * n_addcov     Number of additive covariates
 *
 * intcov	Interacting covariates matrix, intcov[mar][ind]
 *
 * n_intcov     Number of interacting covariates
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * result       Result vector of length [n_pos1*n_pos2];
 *
 **********************************************************************/

void scantwo_imp(int n_ind, int same_chr, int n_pos1, int n_pos2, 
		 int n_gen1, int n_gen2, int n_draws, int ***Draws1, 
		 int ***Draws2, double **Addcov, int n_addcov, 
		 double **Intcov, int n_intcov, double *pheno, 
		 double *weights, double *result)
{

  /* create local variables */
  int i, i1, i2, j; /* loop variants */
  double lrss0, lrss_add, lrss_full, *LODfull, *LODint;
  int n_col_f, n_gen_sq, *iwork, idx; 
  double *dwork;

  /* constants */
  n_gen_sq = n_gen1*n_gen2;
  n_col_f = n_gen_sq + n_addcov + n_intcov*(n_gen_sq-1);

  /* allocate space */
  LODfull = (double *)R_alloc(n_draws, sizeof(double));
  LODint = (double *)R_alloc(n_draws, sizeof(double));
  iwork = (int *)R_alloc(n_col_f, sizeof(int));
  dwork = (double *)R_alloc(n_col_f*n_ind + 2*n_ind + 4*n_col_f,
			    sizeof(double));
			 
  /* adjust phenotypes and covariates using weights */
  /* Note: these are actually square-root of weights */
  for(i=0; i<n_ind; i++) {
    pheno[i] *= weights[i];
    for(j=0; j<n_addcov; j++) Addcov[j][i] *= weights[i];
    for(j=0; j<n_intcov; j++) Intcov[j][i] *= weights[i];
  }

  /* Note that lrss0 is the log10(RSS) for model E(Yi) = b0;
     lrss_add is the log10(RSS) for model 
                 E(Yi) = b0 + b1*q1 + b2*q2;
     lrss_full is the log10(RSS) for model 
                 E(Yi) = b0 + b1*q1 + b2*q2 + b3*(q1*q2); 
     Additive and interactive covariates are included (if any) */

  /* Call nullRss to calculate the RSS for the null model */
  lrss0 = log10(nullRss(pheno, weights, n_ind, Addcov, 
			n_addcov, dwork, iwork));

  /* calculate the LOD score for each pair of markers */
  if(same_chr) { /* if the pair is on the same chromesome */
    for(i1=0; i1<n_pos1-1; i1++) {
      for (i2=i1+1; i2<n_pos1; i2++) {
	for(j=0; j<n_draws; j++) { /* loop over imputations */

	  /* rss for alternative model */
	  altRss2(pheno, weights, n_ind, n_gen1, n_gen1, Draws1[j][i1], 
		  Draws1[j][i2], Addcov, n_addcov, Intcov, n_intcov, 
		  &lrss_add, &lrss_full, dwork, iwork);
	  
	  /* calculate 2 different LOD scores */
	  LODfull[j] = (double)n_ind/2.0*(lrss0-lrss_full);
	  LODint[j] = (double)n_ind/2.0*(lrss_add-lrss_full);

	}
	/* calculate the weight average on the two LOD score vector
	   and fill the result matrix */
	result[i1*n_pos1+i2] = wtaverage(LODfull, n_draws);
	result[i2*n_pos1+i1] = wtaverage(LODint, n_draws);

      } /* end loop over position 1 */
    } /* end loop over position 2 */
  }

  else { /* the pair is for different chromesome */
    idx = n_pos1*n_pos2;
    for(i1=0; i1<n_pos1; i1++) { /* loop over markers on chro 1 */
      for(i2=0; i2<n_pos2; i2++) { /* loop over markers on chr 2 */
	for(j=0; j<n_draws; j++) { /* loop over imputations */
	  /* rss for alternative model */
	  altRss2(pheno, weights, n_ind, n_gen1, n_gen2, Draws1[j][i1], 
		  Draws2[j][i2], Addcov, n_addcov, Intcov, n_intcov, 
		  &lrss_add, &lrss_full, dwork, iwork);
	  
	  /* calculate 2 different LOD scores */
	  LODfull[j] = (double)n_ind/2.0*(lrss0-lrss_full);
	  LODint[j] = (double)n_ind/2.0*(lrss_add-lrss_full);
	}
	/* calculate the weight average on the two LOD score vector
	   and fill the result matrix */
	result[i1 + n_pos1*i2] = wtaverage(LODint, n_draws);
	result[idx + i1 + n_pos1*i2] = wtaverage(LODfull, n_draws);

      } /* end loop over chromesome 2 */
    } /* end loop over chromesome 1 */
  } 

/* end of scantwo_imp() */
}


void altRss2(double *pheno, double *weights, 
	     int n_ind, int n_gen1, int n_gen2, 
	     int *Draws1, int *Draws2, double **Addcov, int n_addcov, 
	     double **Intcov, int n_intcov, double *lrss_add, 
	     double *lrss_full, double *dwork, int *iwork)
{
  int ny, *jpvt, i, k, s;
  int n_col_a, n_col_f, n_gen_sq;
  double *work, *x, *qty, *qraux, *coef, *resid, tol;
  
  /* constants */
  n_gen_sq = n_gen1*n_gen2;
  n_col_a = (n_gen1+n_gen2-1) + n_addcov + n_intcov*(n_gen1+n_gen2-2);
  n_col_f = n_gen_sq + n_addcov + n_intcov*(n_gen_sq-1);
  tol = TOL; /* tolerance for linear regression */
  ny = 1; /* number of phenotypes */

  /* dwork must be length n_ind*n_col_f + 2*n_ind + 4*n_col_f */
  x = dwork;
  coef = x + n_ind*n_col_f;
  resid = coef + n_col_f;
  qty = resid + n_ind;
  qraux = qty + n_ind;
  work = qraux + n_col_f;
  
  /* iwork must be length n_col_f */
  jpvt = iwork;

  /* ADDITIVE MODEL */
  /* zero out X matrix */
  for(i=0; i<n_ind*n_col_a; i++) x[i] = 0.0;

  /* fill up X matrix */
  for(i=0; i<n_ind; i++) {
    x[i+(Draws1[i]-1)*n_ind] = weights[i]; /* QTL 1 */
    s = n_gen1;
    if(Draws2[i] < n_gen2) /* QTL 2 */
      x[i+(Draws2[i]-1+s)*n_ind] = weights[i]; 
    s += (n_gen2-1);
    for(k=0; k<n_addcov; k++) /* add cov */
      x[i+(k+s)*n_ind] = Addcov[k][i];
    s += n_addcov;
    for(k=0; k<n_intcov; k++) {
      if(Draws1[i] < n_gen1) /* QTL1 x int cov */
	x[i+(Draws1[i]-1+s)*n_ind] = Intcov[k][i];
      s += (n_gen1-1);
      if(Draws2[i] < n_gen2) /* QTL 2 x int cov*/
	x[i+(Draws2[i]-1+s)*n_ind] = Intcov[k][i];
      s += (n_gen2-1);
    }
  } /* end loop over individuals */

  /* regression */
  dqrls_(x, &n_ind, &n_col_a, pheno, &ny, &tol, coef, resid,
	 qty, &k, jpvt, qraux, work);

  /* calculate RSS */
  *lrss_add = 0.0;
  for(i=0; i<n_ind; i++)
    (*lrss_add) += (resid[i]*resid[i]);
  (*lrss_add) = log10(*lrss_add);

  /* INTERACTIVE MODEL */
  /* zero out X matrix */
  for(i=0; i<n_ind*n_col_f; i++) x[i] = 0.0;

  /* fill up X matrix */
  for(i=0; i<n_ind; i++) {
    x[i+(Draws1[i]-1)*n_ind] = weights[i]; /* QTL 1 */
    s = n_gen1;
    if(Draws2[i] < n_gen2) /* QTL 2 */
      x[i+(Draws2[i]-1+s)*n_ind] = weights[i]; 
    s += (n_gen2-1);
    if(Draws1[i] < n_gen1 && Draws2[i] < n_gen2) /* QTL x QTL */
      x[i+((Draws1[i]-1)*(n_gen2-1)+Draws2[i]-1+s)*n_ind] = weights[i];
    s += ((n_gen1-1)*(n_gen2-1));
    for(k=0; k<n_addcov; k++) /* add cov */
      x[i+(k+s)*n_ind] = Addcov[k][i];
    s += n_addcov;
    for(k=0; k<n_intcov; k++) {
      if(Draws1[i] < n_gen1) /* QTL1 x int cov */
	x[i+(Draws1[i]-1+s)*n_ind] = Intcov[k][i];
      s += (n_gen1-1);
      if(Draws2[i] < n_gen2) /* QTL 2 x int cov */
	x[i+(Draws2[i]-1+s)*n_ind] = Intcov[k][i];
      s += (n_gen2-1);
      /* QTL x QTL x int cov */
      if(Draws1[i] < n_gen1 && Draws2[i] < n_gen2) 
	x[i+((Draws1[i]-1)*(n_gen2-1)+Draws2[i]-1+s)*n_ind] = 
	  Intcov[k][i];
      s += ((n_gen1-1)*(n_gen2-1));
    }
  } /* end loop over individuals */

  /* regression */
  dqrls_(x, &n_ind, &n_col_f, pheno, &ny, &tol, coef, resid,
	 qty, &k, jpvt, qraux, work);

  /* calculate RSS */
  *lrss_full = 0.0;
  for(i=0; i<n_ind; i++)
    (*lrss_full) += (resid[i]*resid[i]);
  (*lrss_full) = log10(*lrss_full);
}

/* end of scantwo_imp.c */
