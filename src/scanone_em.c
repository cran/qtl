/**********************************************************************
 * 
 * scanone_em.c
 *
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 *
 * last modified Nov, 2001
 * first written Nov, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a genome scan with a 
 * single QTL model by interval mapping (the EM algorithm).
 *
 * Contains: R_scanone_em, scanone_em
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
#include "scanone_em.h"
#include "scanone_em_covar.h"

/**********************************************************************
 * 
 * R_scanone_em
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_em.
 * 
 **********************************************************************/

void R_scanone_em(int *n_ind, int *n_pos, int *n_gen, 
		  double *genoprob, double *addcov, int *n_addcov,
		  double *intcov, int *n_intcov, double *pheno,
		  double *result, int *std_start, double *start,
		  int *maxit, double *tol, int *trace)
{
  double ***Genoprob, **Result, **work, **Addcov, **Intcov;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  reorg_errlod(*n_pos, *n_gen+2, result, &Result);
  allocate_dmatrix(4,*n_gen, &work);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  if(*n_addcov == 0 && *n_intcov == 0) { /* no covariates */
    /* Read R's random seed */
    GetRNGstate();

    scanone_em(*n_ind, *n_pos, *n_gen, Genoprob, pheno, Result, 
	       *std_start, start, *maxit, *tol, work);

    /* Write R's random seed */
    PutRNGstate();
  }
  else { /* interval mapping with covariates */
    scanone_em_covar(*n_ind, *n_pos, *n_gen, Genoprob, Addcov,
		     *n_addcov, Intcov, *n_intcov, pheno, result,
		     *maxit, *tol, *trace);
  }
}

/**********************************************************************
 * 
 * scanone_em
 *
 * Performs genotype scan using interval mapping.  (The multipoint
 * genotype probabilities have already been calculated in 
 * calc.genoprob)
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *
 * pheno        Phenotype data, as a vector
 *
 * Result       Result matrix of size [n_pos x (n_gen+2)]; upon return, 
 *              the first column contains the log10 likelihood, the 
 *              next set contain estimated genotype-specific means, and 
 *              the last column contains the estimated residual SD
 *
 * std_start    If 1, use the usual starting points [initial weights as 
 *                    Pr(QTL geno | marker genotypes)]
 *              If -1, use iid Bernoulli(1/2) for initial weights.
 *              If 0, use the specified values for the means and SD as
 *                    the starting point.
 *
 * start        If std_start = 0, use these as initial estimates of the
 *              genotype-specific means and the residual SD.
 * 
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * work         Workspace of dimension 4 x n_gen
 *
 **********************************************************************/

void scanone_em(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		double *pheno, double **Result, int std_start, 
                double *start,
		int maxit, double tol, double **work)
{
  int i, j, k, s, flag=0;
  double s1, s2, s3, oldsig, r;

  for(i=0; i<n_pos; i++) { /* loop over marker positions */

    /* initiate EM */
    s1 = 0.0;

    if(std_start==0) { /* specified starting point */
      for(k=0; k<n_gen; k++) work[1][k] = start[k];
      oldsig = start[n_gen];
    }
    else {
      if(std_start == 1) { /* the usual starting points */
	for(k=0; k<n_gen; k++) {
	  work[1][k] = s2 = s3 = 0.0;
	  for(j=0; j<n_ind; j++) {
	    s2 += Genoprob[k][i][j]; /* count up numbers */
	    work[1][k] += Genoprob[k][i][j]*pheno[j]; /* means */
	    s3 += Genoprob[k][i][j]*pheno[j]*pheno[j]; /* RSS */
	  }
	  s1 += (s3 - work[1][k]*work[1][k]/s2);
	  work[1][k] /= s2;
	}
	oldsig = sqrt(s1/(double)(n_ind));
      }
      else { /* start using random weights */
	for(k=0; k<n_gen; k++) {
	  work[1][k] = s2 = s3 = 0.0;
	  for(j=0; j<n_ind; j++) {
	    r = unif_rand()/(double)(n_gen); 
	    s2 += r; /* count up numbers */
	    work[1][k] += r*pheno[j]; /* means */
	    s3 += r*pheno[j]*pheno[j]; /* RSS */
	  }
	  s1 += (s3 - work[1][k]*work[1][k]/s2);
	  work[1][k] /= s2;
	}
	oldsig = sqrt(s1/(double)(n_ind));
      }
    }

    for(s=0; s < maxit; s++) { /* EM iterations */
    
      for(k=0; k<n_gen; k++) 
	Result[k+1][i] = work[2][k] = work[3][k] = 0.0;
      Result[n_gen+1][i] = 0.0;

      for(j=0; j<n_ind; j++) { /* loop over individuals */
	/* E-step */
	s1=0.0;
	for(k=0; k<n_gen; k++) 
	  s1 += (work[0][k] = Genoprob[k][i][j]*
		 dnorm(pheno[j],work[1][k],oldsig,0));
	for(k=0; k<n_gen; k++) 
	  work[0][k] /= s1;
	
	/* M-step */
	for(k=0; k<n_gen; k++) {
	  work[2][k] += work[0][k]; /* count up numbers */
	  Result[k+1][i] += work[0][k] * pheno[j]; /* means */
	  work[3][k] += work[0][k] * pheno[j] * pheno[j]; /* RSS */
	}
      }
      
      /* complete M-step */
      for(k=0; k<n_gen; k++) {
	Result[1+n_gen][i] += (work[3][k] - Result[k+1][i]*Result[k+1][i]/
			       work[2][k]);
	Result[k+1][i] /= work[2][k];
      }
      Result[1+n_gen][i] = sqrt(Result[1+n_gen][i]/(double)n_ind);

      /* check for convergence */
      flag = 0;
      for(k=0; k<n_gen; k++) {
	if(fabs(Result[k+1][i] - work[1][k]) > tol*(fabs(work[1][k])+tol*100.0)) {
	  flag = 1;
	  break;
	}
      }
      if(fabs(Result[n_gen+1][i] - oldsig) > tol*(oldsig+tol*100.0)) flag = 1;

      if(!flag) break;

      oldsig = Result[n_gen+1][i];
      for(k=0; k<n_gen; k++) work[1][k] = Result[1+k][i];

    } /* end of EM iterations */

    if(flag) warning("Didn't converge!\n");

    /* calculate negative log lik */
    Result[0][i] = 0.0;
    for(j=0; j<n_ind; j++) {
      s1 = 0.0;
      for(k=0; k<n_gen; k++) 
	s1 += Genoprob[k][i][j] * dnorm(pheno[j], Result[1+k][i], 
					Result[1+n_gen][i], 0);
      Result[0][i] -= log10(s1);
    }
    Result[1+n_gen][i] *= sqrt((double)n_ind/(double)(n_ind-n_gen));

  } /* end loop over marker positions */
}

/* end of scanone_em.c */

