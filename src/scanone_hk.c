/**********************************************************************
 * 
 * scanone_hk.c
 *
 * copyright (c) 2001-2, Karl W Broman, Johns Hopkins University
 *
 * last modified Oct, 2002
 * first written Nov, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a genome scan with a 
 * single QTL model by Haley-Knott regression
 *
 * Contains: R_scanone_hk, scanone_hk
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
#include "scanone_hk.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scanone_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_hk.
 * 
 **********************************************************************/

void R_scanone_hk(int *n_ind, int *n_pos, int *n_gen,
		  double *genoprob, double *addcov, int *n_addcov, 
                  double *intcov, int *n_intcov, double *pheno,
		  double *weights, double *result)
{
  double ***Genoprob, **Result, **Addcov, **Intcov;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  reorg_errlod(*n_pos, *n_gen+2, result, &Result);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scanone_hk(*n_ind, *n_pos, *n_gen, Genoprob, Addcov, *n_addcov,
	     Intcov, *n_intcov, pheno, weights, Result);
}

/**********************************************************************
 * 
 * scanone_hk
 *
 * Performs genotype scan using the Haley-Knott regression method
 * (regressing phenotypes on conditional genotype probabilities; the
 * multipoint genotype probabilities have already been calculated in
 * calc.genoprob)
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *              Indexed as Genoprob[gen][pos][ind]
 *
 * Addcov       Matrix of additive covariates: Addcov[cov][ind]
 * 
 * n_addcov     Number of columns of Addcov
 *
 * Intcov       Number of interactive covariates: Intcov[cov][ind]
 *
 * n_intcov     Number of columns of Intcov
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * Result       Result matrix of size [n_pos x (n_gen+2)]; upon return, 
 *              the first column contains the LOD, the next set contain
 *              estimated genotype-specific means, and the last column
 *              contains the estimated residual SD
 *
 **********************************************************************/

void scanone_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
                double **Addcov, int n_addcov, double **Intcov, 
		int n_intcov, double *pheno, double *weights, 
		double **Result)
{
  int ny, *jpvt, k, i, j, ncol, ncol0, k2, s;
  double *work, *x, *qty, *qraux, *coef, *resid, tol;

  /* tolerance for linear regression */
  tol = TOL;

  ncol = n_gen + (n_gen-1)*n_intcov+n_addcov;
  ncol0 = n_addcov+1;

  /* allocate space and set things up*/
  x = (double *)R_alloc(n_ind*ncol, sizeof(double));
  coef = (double *)R_alloc(ncol, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  jpvt = (int *)R_alloc(ncol, sizeof(int));
  qraux = (double *)R_alloc(ncol, sizeof(double));
  work = (double *)R_alloc(2 * ncol, sizeof(double));
  ny = 1;

  /* NULL model is now done in R ********************
     (only do it once!)
  for(j=0; j<n_ind; j++) {
    x[j] = 1.0;
    for(k=0; k<n_addcov; k++) 
      x[j+(k+1)*n_ind] = Addcov[k][j];
  }
  dqrls_(x, &n_ind, &ncol0, pheno, &ny, &tol, coef, resid,
	 qty, &k, jpvt, qraux, work);
  rss0 = 0.0;
  for(j=0; j<n_ind; j++)  rss0 += (resid[j]*resid[j]);
  Null model is now done in R ********************/

  for(j=0; j<n_ind; j++) 
    pheno[j] *= weights[j];
  /* note: weights are really square-root of weights */

  for(i=0; i<n_pos; i++) { /* loop over positions */
    for(k=0; k<n_gen; k++) jpvt[k] = k;

    /* fill up X matrix */
    for(j=0; j<n_ind; j++) {
      for(k=0; k<n_gen; k++)
	x[j+k*n_ind] = Genoprob[k][i][j]*weights[j]; 
      for(k=0; k<n_addcov; k++)
	x[j+(k+n_gen)*n_ind] = Addcov[k][j]*weights[j];
      for(k=0,s=0; k<n_gen-1; k++)
	for(k2=0; k2<n_intcov; k2++,s++) 
	  x[j+(n_gen+n_addcov+s)*n_ind] = Genoprob[k][i][j]*Intcov[k2][j]*weights[j];
    }

    /* linear regression of phenotype on QTL genotype probabilities */
    dqrls_(x, &n_ind, &ncol, pheno, &ny, &tol, coef, resid,
	   qty, &k, jpvt, qraux, work);

    /* RSS */
    Result[0][i] = 0.0;
    for(j=0; j<n_ind; j++) Result[0][i] += (resid[j]*resid[j]);

    if(n_addcov == 0 && n_intcov == 0) {
      /* save coefficients only if no covariates */
      /* re-scramble coefficients */
      for(j=0; j<n_gen; j++) Result[1+j][i] = coef[jpvt[j]];
    
      /* residual SD */
      Result[1+n_gen][i] = sqrt(Result[0][i]/(double)(n_ind-n_gen));
    }

    /* log10 likelihood */
    Result[0][i] = (double)n_ind/2.0*log10(Result[0][i]);
  } /* end loop over positions */

}

/* end of scanone_hk.c */
