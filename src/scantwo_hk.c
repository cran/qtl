/**********************************************************************
 * 
 * scantwo_hk.c
 *
 * copyright (c) 2001-3, Karl W Broman, Johns Hopkins University
 *
 * last modified Dec, 2003
 * first written Nov, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a two-dimensional genome scan with  
 * a two-QTL model by Haley-Knott regression
 *
 * Contains: R_scantwo_1chr_hk, scantwo_1chr_hk, 
 *           R_scantwo_2chr_hk, scantwo_2chr_hk
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
#include "scantwo_hk.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scantwo_1chr_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_1chr_hk.
 * 
 **********************************************************************/

void R_scantwo_1chr_hk(int *n_ind, int *n_pos, int *n_gen,
		       double *genoprob, double *pairprob, 
		       double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *weights, 
		       double *result)
{
  double ***Genoprob, **Result, **Addcov, **Intcov, *****Pairprob;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  reorg_pairprob(*n_ind, *n_pos, *n_gen, pairprob, &Pairprob);
  reorg_errlod(*n_pos, *n_pos, result, &Result);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scantwo_1chr_hk(*n_ind, *n_pos, *n_gen, Genoprob, Pairprob, 
		  Addcov, *n_addcov, Intcov, *n_intcov, 
		  pheno, weights, Result);
}

/**********************************************************************
 * 
 * scantwo_1chr_hk
 *
 * Performs a 2-dimensional genome scan using the Haley-Knott 
 * regression method (regressing phenotypes on conditional genotype 
 * probabilities) for a two-QTL model with the two QTL residing on
 * the same chromosome.
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
 * Pairprob     Array of joint genotype probabilities for QTL
 *              pairs; indexed as Pairprob[gen1][gen2][pos1][pos2][ind]
 *              where pos2 > pos1 (for pos2 <= pos1, points to nothing)
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
 * Result       Result matrix of size [n_pos x n_pos]; the lower
 *              triangle (row > col) contains the joint LODs while 
 *              the upper triangle (row < col) contains the LODs for 
 *              testing epistasis.
 *              Note: indexed as Result[col][row]
 *
 **********************************************************************/

void scantwo_1chr_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		     double *****Pairprob, double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     double *weights, double **Result)
{
  int ny, *jpvt, i, i2, j, k, k2, k3, s;
  int n_col_a, n_col_f, n_gen_sq;
  double *work, *x, *qty, *qraux, *coef, *resid, tol;

  /* tolerance for linear regression */
  tol = TOL;

  n_gen_sq = n_gen*n_gen;
  /* no. param in additive QTL model */
  n_col_a = (n_gen*2-1)+n_addcov+n_intcov*(n_gen-1)*2; 
  /* no. param full model */
  n_col_f = n_gen_sq+n_addcov+n_intcov*(n_gen_sq-1); 

  /* allocate space and set things up*/
  x = (double *)R_alloc(n_ind*n_col_f, sizeof(double));
  coef = (double *)R_alloc(n_col_f, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  jpvt = (int *)R_alloc(n_col_f, sizeof(int));
  qraux = (double *)R_alloc(n_col_f, sizeof(double));
  work = (double *)R_alloc(2 * n_col_f, sizeof(double));
  ny = 1;

  /* modify pheno, Addcov and Intcov with weights */
  for(j=0; j<n_ind; j++) {
    pheno[j] *= weights[j];
    for(k=0; k<n_addcov; k++) Addcov[k][j] *= weights[j];
    for(k=0; k<n_intcov; k++) Intcov[k][j] *= weights[j];
  }    

  /* NULL model is now done in R ********************
     (only do it once!)
  for(j=0; j<n_ind; j++) {
    x[j] = 1.0;
    for(k=0; k<n_addcov; k++) 
      x[j+(k+1)*n_ind] = Addcov[k][j];
  }
  F77_CALL(dqrls)(x, &n_ind, &n_col_0, pheno, &ny, &tol, coef, resid,
   	          qty, &k, jpvt, qraux, work);
  lrss0 = 0.0;
  for(j=0; j<n_ind; j++)  lrss0 += (resid[j]*resid[j]);
  lrss0 = log10(lrss0);
  Null model is now done in R ********************/

  for(i=0; i<n_pos-1; i++) { 
    for(i2=i+1; i2<n_pos; i2++) { /* loop over pairs of positions */

      /* ADDITIVE MODEL */
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob[k][i][j]*weights[j];  /* s keeps track of column */
	for(k=0; k<n_gen-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob[k][i2][j]*weights[j];
	for(k=0; k<n_addcov; k++, s++) /* additive covariates */
	  x[j+s*n_ind] = Addcov[k][j];
	for(k=0; k<n_gen-1; k++) /* interactive x QTL 1 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob[k][i][j]*Intcov[k2][j];
	for(k=0; k<n_gen-1; k++) /* interactive x QTL 2 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob[k][i2][j]*Intcov[k2][j];
      }
      /* linear regression of phenotype on QTL genotype probabilities */
      F77_CALL(dqrls)(x, &n_ind, &n_col_a, pheno, &ny, &tol, coef, resid,
		      qty, &k, jpvt, qraux, work);
      /* RSS */
      Result[i2][i] = 0.0;
      for(j=0; j<n_ind; j++) Result[i2][i] += (resid[j]*resid[j]);
      Result[i2][i] = log10(Result[i2][i]); /* take log base 10*/

      /* INTERACTIVE MODEL */
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob[k][i][j]*weights[j];  /* s keeps track of column */
	for(k=0; k<n_gen-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob[k][i2][j]*weights[j];
	for(k=0; k<n_gen-1; k++)
	  for(k2=0; k2<n_gen-1; k2++,s++) /* QTL 2 x QTL 2 */
	    x[j+s*n_ind] = Pairprob[k][k2][i][i2][j]*weights[j];
	for(k=0; k<n_addcov; k++, s++) /* additive covariates */
	  x[j+s*n_ind] = Addcov[k][j];
	for(k=0; k<n_gen-1; k++) /* interactive x QTL 1 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob[k][i][j]*Intcov[k2][j];
	for(k=0; k<n_gen-1; k++) /* interactive x QTL 2 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob[k][i2][j]*Intcov[k2][j];
	for(k=0; k<n_gen-1; k++) /* interactive x QTL 1 x QTL 2 */
	  for(k2=0; k2<n_gen-1; k2++) 
	    for(k3=0; k3<n_intcov; k3++,s++)
	      x[j+s*n_ind] = Pairprob[k][k2][i][i2][j]*Intcov[k3][j];
      }

      /* linear regression of phenotype on QTL genotype probabilities */
      F77_CALL(dqrls)(x, &n_ind, &n_col_f, pheno, &ny, &tol, coef, resid,
		      qty, &k, jpvt, qraux, work);

      /* RSS */
      Result[i][i2] = 0.0;
      for(j=0; j<n_ind; j++) Result[i][i2] += (resid[j]*resid[j]);
      Result[i][i2] = log10(Result[i][i2]); /* take log base 10*/

      /* convert to LODs */
      Result[i2][i] = (double)n_ind/2.0*(Result[i2][i]-Result[i][i2]);
      Result[i][i2] = (double)n_ind/2.0*Result[i][i2];

    } /* end loop over positions */
  } 
}

/**********************************************************************
 * 
 * R_scantwo_2chr_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_2chr_hk.
 * 
 **********************************************************************/

void R_scantwo_2chr_hk(int *n_ind, int *n_pos1, int *n_pos2, 
		       int *n_gen1, int *n_gen2,
		       double *genoprob1, double *genoprob2,
		       double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *weights,
		       double *result_full, double *result_int)
{
  double ***Genoprob1, ***Genoprob2, **Result_full, **Result_int;
  double **Addcov, **Intcov;

  reorg_genoprob(*n_ind, *n_pos1, *n_gen1, genoprob1, &Genoprob1);
  reorg_genoprob(*n_ind, *n_pos2, *n_gen2, genoprob2, &Genoprob2);
  reorg_errlod(*n_pos1, *n_pos2, result_full, &Result_full);
  reorg_errlod(*n_pos1, *n_pos2, result_int, &Result_int);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scantwo_2chr_hk(*n_ind, *n_pos1, *n_pos2, *n_gen1, *n_gen2, 
		  Genoprob1, Genoprob2, Addcov, *n_addcov, Intcov, 
		  *n_intcov, pheno, weights, Result_full, Result_int);
}

/**********************************************************************
 * 
 * scantwo_2chr_hk
 *
 * Performs a 2-dimensional genome scan using the Haley-Knott 
 * regression method (regressing phenotypes on conditional genotype 
 * probabilities) for a two-QTL model with the two QTL residing on
 * the different chromosomes.
 * 
 * n_ind        Number of individuals
 *
 * n_pos1       Number of marker positions on first chromosome
 *
 * n_pos2       Number of marker positions on second chromosome
 *
 * n_gen1       Number of different genotypes for first chromosome
 *
 * n_gen2       Number of different genotypes for second chromosome
 *
 * Genoprob1    Array of conditional genotype probs for 1st chr
 *              Indexed as Genoprob[gen][pos][ind]
 *
 * Genoprob2    Array of conditional genotype probs for 2nd chr
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
 * Result_full  Result matrix of size [n_pos1 x n_pos2]
 *              containing the joint LODs
 *              Note: indexed as Result[pos2][pos1]
 *
 * Result_int   Result matrix of size [n_pos2 x n_pos1] 
 *              containing the LODs testing interactions
 *              also indexed as Result[pos2][pos1]
 *
 **********************************************************************/

void scantwo_2chr_hk(int n_ind, int n_pos1, int n_pos2, int n_gen1, 
		     int n_gen2, double ***Genoprob1, 
		     double ***Genoprob2, 
		     double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     double *weights,
		     double **Result_full, double **Result_int)
{
  int ny, *jpvt, i, i2, j, k, k2, k3, s;
  int n_col_a, n_col_f, n_gen_sq;
  double *work, *x, *qty, *qraux, *coef, *resid, tol;

  /* tolerance for linear regression */
  tol = TOL;

  n_gen_sq = n_gen1*n_gen2;
  /* no. param in additive QTL model */
  n_col_a = (n_gen1+n_gen2-1)+n_addcov+n_intcov*(n_gen1+n_gen2-2); 
  /* no. param full model */
  n_col_f = n_gen_sq+n_addcov+n_intcov*(n_gen_sq-1); 

  /* allocate space and set things up*/
  x = (double *)R_alloc(n_ind*n_col_f, sizeof(double));
  coef = (double *)R_alloc(n_col_f, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  jpvt = (int *)R_alloc(n_col_f, sizeof(int));
  qraux = (double *)R_alloc(n_col_f, sizeof(double));
  work = (double *)R_alloc(2 * n_col_f, sizeof(double));
  ny = 1;

  /* modify pheno, Addcov and Intcov with weights */
  for(j=0; j<n_ind; j++) {
    pheno[j] *= weights[j];
    for(k=0; k<n_addcov; k++) Addcov[k][j] *= weights[j];
    for(k=0; k<n_intcov; k++) Intcov[k][j] *= weights[j];
  }    

  /* NULL model is now done in R ********************
     (only do it once!)
  for(j=0; j<n_ind; j++) {
    x[j] = 1.0;
    for(k=0; k<n_addcov; k++) 
      x[j+(k+1)*n_ind] = Addcov[k][j];
  }
  F77_CALL(dqrls)(x, &n_ind, &n_col_0, pheno, &ny, &tol, coef, resid,
  	          qty, &k, jpvt, qraux, work);
  lrss0 = 0.0;
  for(j=0; j<n_ind; j++)  lrss0 += (resid[j]*resid[j]);
  lrss0 = log10(lrss0);
  Null model is now done in R ********************/


  for(i=0; i<n_pos1; i++) { 
    for(i2=0; i2<n_pos2; i2++) { /* loop over pairs of positions */

      /* ADDITIVE MODEL */
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen1; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob1[k][i][j]*weights[j];  /* s keeps track of column */
	for(k=0; k<n_gen2-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob2[k][i2][j]*weights[j];
	for(k=0; k<n_addcov; k++, s++) /* additive covariates */
	  x[j+s*n_ind] = Addcov[k][j];
	for(k=0; k<n_gen1-1; k++) /* interactive x QTL 1 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob1[k][i][j]*Intcov[k2][j];
	for(k=0; k<n_gen2-1; k++) /* interactive x QTL 2 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob2[k][i2][j]*Intcov[k2][j];
      }
      /* linear regression of phenotype on QTL genotype probabilities */
      F77_CALL(dqrls)(x, &n_ind, &n_col_a, pheno, &ny, &tol, coef, resid,
		      qty, &k, jpvt, qraux, work);
      /* RSS */
      Result_int[i2][i] = 0.0;
      for(j=0; j<n_ind; j++) Result_int[i2][i] += (resid[j]*resid[j]);
      Result_int[i2][i] = log10(Result_int[i2][i]); /* take log base 10*/

      /* INTERACTIVE MODEL */
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen1; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob1[k][i][j]*weights[j];  /* s keeps track of column */
	for(k=0; k<n_gen2-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob2[k][i2][j]*weights[j];
	for(k=0; k<n_gen1-1; k++)
	  for(k2=0; k2<n_gen2-1; k2++,s++) /* QTL 2 x QTL 2 */
	    x[j+s*n_ind] = Genoprob1[k][i][j]*Genoprob2[k2][i2][j]*weights[j];
	for(k=0; k<n_addcov; k++, s++) /* additive covariates */
	  x[j+s*n_ind] = Addcov[k][j];
	for(k=0; k<n_gen1-1; k++) /* interactive x QTL 1 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob1[k][i][j]*Intcov[k2][j];
	for(k=0; k<n_gen2-1; k++) /* interactive x QTL 2 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob2[k][i2][j]*Intcov[k2][j];
	for(k=0; k<n_gen1-1; k++) /* interactive x QTL 1 x QTL 2 */
	  for(k2=0; k2<n_gen2-1; k2++) 
	    for(k3=0; k3<n_intcov; k3++,s++)
	      x[j+s*n_ind] = Genoprob1[k][i][j]*Genoprob2[k2][i2][j]*
		Intcov[k3][j];

      }

      /* linear regression of phenotype on QTL genotype probabilities */
      F77_CALL(dqrls)(x, &n_ind, &n_col_f, pheno, &ny, &tol, coef, resid,
		      qty, &k, jpvt, qraux, work);

      /* RSS */
      Result_full[i2][i] = 0.0;
      for(j=0; j<n_ind; j++) Result_full[i2][i] += (resid[j]*resid[j]);
      Result_full[i2][i] = log10(Result_full[i2][i]); /* take log base 10*/

      /* convert to LODs */
      Result_int[i2][i] = (double)n_ind/2.0*
	(Result_int[i2][i] - Result_full[i2][i]);
      Result_full[i2][i] = (double)n_ind/2.0*Result_full[i2][i];

    } /* end loop over positions */
  } 
}

/* end of scantwo_hk.c */
