/**********************************************************************
 * 
 * scanone_imp.c
 *
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 *                 and Hao Wu, The Jackson Laboratory
 *
 * This file is written by Hao Wu (hao@jax.org), 
 * with slight modifications by Karl Broman.
 *
 * last modified Nov, 2001
 * first written Nov, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a genome scan with a 
 * single QTL model by imputation.  
 *
 * Contains: R_scanone_imp, scanone_imp, nullRss, altRss
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
#include "scanone_imp.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scanone_imp
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_imp.
 * 
 **********************************************************************/

void R_scanone_imp(int *n_ind, int *n_pos, int *n_gen, int *n_draws, 
		   int *draws, double *addcov, int *n_addcov, 
		   double *intcov, int *n_intcov, double *pheno, 
		   double *result, int *trim, int *direct)
{
  /* reorganize draws */
  int ***Draws;
  double **Addcov, **Intcov;
  
  reorg_draws(*n_ind, *n_pos, *n_draws, draws, &Draws);

  /* reorganize addcov and intcov (if they are not empty) */
  /* currently reorg_errlod function is used to reorganize the data */
  if(*n_addcov != 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov != 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);
      
  scanone_imp(*n_ind, *n_pos, *n_gen, *n_draws, Draws, 
	      Addcov, *n_addcov, Intcov, *n_intcov, pheno, result,
	      *trim, *direct);
}


/**********************************************************************
 * 
 * scanone_imp
 *
 * Performs genotype scan using the pseudomarker algorithm (imputation) 
 * method of Sen and Churchill (2001).
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * n_draws      Number of impiutations
 *
 * Draws        Array of genotype imputations, indexed as 
 *              Draws[repl][mar][ind]
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
 * result       Result vector of length [n_pos]; upon return, contains
 *              the "LPD" (log posterior distribution of QTL location).
 * 
 * trim         If 1, trim off the top and bottom log2(n.draws) LODs 
 *
 * direct       If 1, return log10[mean(10^LOD)]; if 0, return
 *              mean(LOD) + var(LOD)/2
 *
 **********************************************************************/

void scanone_imp(int n_ind, int n_pos, int n_gen, int n_draws, 
		 int ***Draws, double **addcov, int n_addcov, 
		 double **intcov, int n_intcov, double *pheno, 
		 double *result, int trim, int direct)
{

  /* create local variables */
  int i, j, k; /* loop variants */
  int idx; /* the smallest and biggest idx LOD scores need to be thrown */
  double tol, rss, lrss0, sum, sums, meanLOD, varLOD, *LOD;
  double *dwork;
  int *iwork, sizefull;

  sizefull = n_gen + n_addcov + n_intcov*(n_gen-1);
  dwork = (double *)R_alloc(n_ind*(sizefull+2)+4*sizefull, sizeof(double));
  iwork = (int *)R_alloc(sizefull, sizeof(int));

  LOD = (double *)R_alloc(n_draws, sizeof(double));

  /* tolerance for linear regression */
  tol = TOL;

  /* calculate the number of LOD needs to be thrown */
  if(trim) idx = (int) floor( 0.5*log(n_draws)/log(2) );
  else idx=0;

  /* Call nullRss to calculate the RSS for the null model */
  lrss0 = log10(nullRss(pheno, n_ind, addcov, n_addcov, dwork, iwork));

  /* calculate the LOD score for each marker */
  for(i=0; i<n_pos; i++) { /* loop over positions */

    for(j=0; j<n_draws; j++) { /* loop over imputations */

      /* call altRss to calcualte the RSS for alternative model,
	 given marker and imputatin number */
      rss = altRss(pheno, n_ind, n_gen, Draws[j][i], addcov, n_addcov,
		   intcov, n_intcov, dwork, iwork);

      /* calculate the LOD score for this marker in this imputation */
      LOD[j] = (double)n_ind/2.0*(lrss0-log10(rss));

    } /* end loop over imputations */

    if(trim) /* sort the LOD scores in ascending order */
      R_rsort(LOD, n_draws);

    if(direct) { /* result = log10[mean(10^LOD)] */
      result[i] = 0.0;
      for(k=idx; k<n_draws-idx; k++) 
	result[i] += exp(LOD[k]*log(10.0));
      result[i] = log10(result[i]/(double)(n_draws-2*idx));
    }
    else { /* result = mean(LOD) + var(LOD)/2 */
      /* get a new list of LOD scores, throwing the biggest 
	 and smallest idx LOD scores. */
      for(k=idx, sum=0.0; k<n_draws-idx; k++) 
	sum += LOD[k]; /* calculate the sum of newLOD in the same loop */

      /* calculate the mean of newLOD */
      meanLOD = sum / (n_draws-2*idx); 

      /* calculate the variance of newLOD */
      for(k=idx,sums=0.0; k<n_draws-idx; k++) 
	sums += (LOD[k]-meanLOD) * (LOD[k]-meanLOD);
      varLOD = sums/(n_draws-2*idx-1);

      /* The return value */
      result[i] = meanLOD + log(10.0)*0.5*varLOD;
    }

  } /* end loop over positions */
}


/* function to calculate the null model RSS for scanone_imp */

double nullRss(double *pheno, int n_ind, double **addcov, int n_addcov, 
	       double *dwork, int *iwork)
{
  /* create local variables */
  int i, j, ny, k, *jpvt0, ncolx0;
  double tol, rss0=0.0, s=0.0, m; 
  double *work0, *x0, *qraux0, *coef0, *qty0, *resid0;

  tol = TOL;
  ny = 1; /* number of phenotypes */

  /*if there's no additive covariates, use two-pass 
    algorithm to calculate null RSS */
  if ( n_addcov == 0 ) {
    for(i=0; i<n_ind; i++) 
      s += pheno[i];
    m = s / (double)n_ind;
    
    for(i=0.0; i<n_ind; i++) 
      rss0 += (pheno[i]-m)*(pheno[i]-m);
  }
  else {
    /* if there are additive covariates, fit linear regression model */
    
    ncolx0 = 1 + n_addcov; /* number of columns in x0 matrix */

    /* point to areas of workspaces */
    /* dwork length n_ind*ncolx0 + 4*ncolx0 + 2*n_ind */
    x0 = dwork; /* length n_ind*ncolx0 */
    coef0 = x0 + n_ind*ncolx0; /* length ncolx0 */
    work0 = coef0 + ncolx0;  /* length 2*ncolx0 */
    qraux0 = work0 + 2*ncolx0; /* length ncolx0 */
    resid0 = qraux0 + ncolx0; /* length n_ind */
    qty0 = resid0 + n_ind; /* length n_ind */
    
    jpvt0 = iwork; /* length ncolx0 */
    
    /* redidual and qty vector will use the same pointers 
       as in the alternative model*/
    
    /* fill up x0 matrix */
    for (i=0; i<n_ind; i++) {
      x0[i] = 1.0; /* the first row (column in Fortran) are all 1s */
      for(j=0; j<n_addcov; j++) 
	/* note that addcov[mar][ind] */
    	x0[(j+1)*n_ind+i] = addcov[j][i]; 
    }

    k = 0;

    /* fit linear regression model */
    dqrls_(x0, &n_ind, &ncolx0, pheno, &ny, &tol, coef0, resid0,
	   qty0, &k, jpvt0, qraux0, work0);

    /* calculate the null RSS */
    for (i=0,rss0=0.0; i<n_ind; i++) 
      rss0 += resid0[i]*resid0[i];
  }

  /* return rss0 */
  return(rss0);

}


/* function to calculate the alternative model RSS. 
   This function is called by scanone_imp */

double altRss(double *pheno, int n_ind, int n_gen, int *Draws,
	      double **addcov, int n_addcov, double **intcov, int n_intcov,
	      double *dwork, int *iwork)
{
  /* create local variables */

  int ny, *jpvt, k, i, j, s, s2, ncolx;
  double *work, *x, *qty, *qraux, *coef, *resid;
  double tol, rss;

  /* tolerance for linear regression */
  tol = TOL;

  /* number of columns in design matrix X */
  ncolx = n_gen + n_addcov + n_intcov*(n_gen-1); 

  ny = 1; /* number of phenotypes */
  k = 0;

  /* fill up design matrix */
  x = dwork;
  for(i=0; i<n_ind; i++) {
    /* QTL genotypes */
    for(s=0; s<n_gen; s++) {
      if(Draws[i] == s+1) x[i+s*n_ind] = 1.0;
      else x[i+s*n_ind] = 0.0;
    }

    /* Additive covariates */
    for(s=0, s2=n_gen; s<n_addcov; s++, s2++) 
      x[i+s2*n_ind] = addcov[s][i];
      
    /* Interactive covariates */
    for(s=0; s<n_intcov; s++) {
      for(j=0; j<n_gen-1; j++, s2++) {
	if(Draws[i] == j+1) x[i+n_ind*s2] = intcov[s][i];
	else x[i+n_ind*s2] = 0.0;
      }
    }
    
  } /* end loop over individuals */
  /* Done filling up X matrix */

  /* point to rest of workspace */
  jpvt = iwork;
  work = dwork + n_ind*ncolx;
  qty = work + 2*ncolx;
  qraux = qty + n_ind;
  coef = qraux + ncolx;
  resid = qraux+ncolx;
  ny = 1;

  /* call Fortran function to fit the linear regression model */
  dqrls_(x, &n_ind, &ncolx, pheno, &ny, &tol, coef, resid,
	 qty, &k, jpvt, qraux, work);

  /* calculate RSS */
  for(i=0, rss=0.0; i<n_ind; i++)
    rss += resid[i]*resid[i];

  /* return rss */
  return(rss);

  /* end of function */
}

/* end of scanone_imp.c */
