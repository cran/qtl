/**********************************************************************
 * 
 * discan.c
 *
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 *
 * last modified Oct, 2001
 * first written Oct, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a genome scan with a 
 * single QTL model
 *
 * Contains: R_discan_mr, discan_mr,
 *           R_discan_im, discan_im
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "util.h"
#include "discan.h"

/**********************************************************************
 * 
 * R_discan_mr
 *
 * Wrapper for call from R; reorganizes genotype and result matrix
 * and calls discan_mr.
 * 
 **********************************************************************/

void R_discan_mr(int *n_ind, int *n_pos, int *n_gen,
		    int *geno, double *pheno, double *result)
{
  int **Geno;
  double **Result;

  reorg_geno(*n_ind, *n_pos, geno, &Geno);
  reorg_errlod(*n_pos, *n_gen+1, result, &Result);

  discan_mr(*n_ind, *n_pos, *n_gen, Geno, pheno, Result);
}

/**********************************************************************
 * 
 * R_discan_im
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls discan_im.
 * 
 **********************************************************************/

void R_discan_im(int *n_ind, int *n_pos, int *n_gen, 
		 double *genoprob, double *pheno, double *result, 
		 int *maxit, double *tol)
{
  double ***Genoprob, **Result, **work;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  reorg_errlod(*n_pos, *n_gen+1, result, &Result);
  allocate_dmatrix(3, *n_gen, &work);

  discan_im(*n_ind, *n_pos, *n_gen, Genoprob, pheno, Result, 
	    *maxit, *tol, work);

}

/**********************************************************************
 * 
 * discan_mr
 *
 * Performs genotype scan using marker regression for a dichotomous trait
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Geno         Genotype matrix
 *
 * pheno        Phenotype data, as a vector
 *
 * Result       Result matrix of size [n_pos x (n_gen+1)]; upon return, 
 *              the first column contains the RSS and the rest contain
 *              estimated genotype-specific probabilities
 *
 **********************************************************************/

void discan_mr(int n_ind, int n_pos, int n_gen, int **Geno, 
		  double *pheno, double **Result)
{
  int i, j, k, n, tp, *ng, *np;

  /* number of individuals in each genotype group */
  ng = (int *)R_alloc(n_gen, sizeof(int));

  /* number of individuals with phenotype=1 in each genotype group */
  np = (int *)R_alloc(n_gen, sizeof(int));

  for(i=0; i<n_pos; i++) {
    Result[0][i] = 0.0; n=tp=0; 
    for(j=0; j< n_gen; j++) {
      ng[j] = np[j] = 0;

      for(k=0; k<n_ind; k++) {
	if(Geno[i][k] == j+1) {
	  if(pheno[k]) {
	    np[j]++;
	    tp++;
	  }
	  n++; ng[j]++;
	}
      }
      
      if(ng[j] > 0)
	Result[j+1][i] = (double)np[j]/(double)ng[j];
      else /* no individuals with this genotype */
	Result[j+1][i] = NA_REAL;
    } /* loop over genotype groups */

    /* calculate LOD score */
    for(j=0; j<n_gen; j++) {
      if(np[j] > 0 && np[j] < ng[j]) 
	Result[0][i] += ((double)np[j]*log10(Result[j+1][i]) +
			 (double)(ng[j]-np[j])*log10(1.0-Result[j+1][i]));
    }
    if(tp > 0 && tp < n) 
      Result[0][i] -= ((double)tp*log10((double)tp/(double)n) +
		       (double)(n-tp)*log10((double)(n-tp)/(double)n));

  } /* loop over marker positions */
}

/**********************************************************************
 * 
 * discan_im
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
 * Result       Result matrix of size [n_pos x (n_gen+1)]; upon return, 
 *              the first column contains the log10 likelihood and the 
 *              rest contain estimated genotype-specific probabilities
 *
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * work         Workspace of dimension 3 x n_gen
 *
 **********************************************************************/

void discan_im(int n_ind, int n_pos, int n_gen, double ***Genoprob,
	       double *pheno, double **Result, 
	       int maxit, double tol, double **work)
{
  int i, j, k, s, flag=0;
  double sw;

  for(i=0; i<n_pos; i++) { /* loop over marker positions */

    /* initiate EM */
    for(k=0; k<n_gen; k++) {
      Result[k+1][i] = sw=0.0;
      for(j=0; j<n_ind; j++) {
	sw += Genoprob[k][i][j]; /* k=gen, i=pos, j=ind */
	if(pheno[j]) Result[k+1][i] += Genoprob[k][i][j];
      }
      Result[k+1][i] /= sw;
    }

    for(s=0; s < maxit; s++) { /* EM iterations */
    
      /* copy over current estimates */
      for(k=0; k<n_gen; k++) {
	work[0][k] = Result[k+1][i]; 
	Result[k+1][i] = work[1][k] = 0.0;
      }
      
      for(j=0; j<n_ind; j++) { /* loop over individuals */
	/* E-step */
	sw = 0.0;
	for(k=0; k<n_gen; k++) {
	  work[2][k] = Genoprob[k][i][j];
	  if(pheno[j]) work[2][k] *= work[0][k];
	  else work[2][k] *= (1.0-work[0][k]);
	  sw += work[2][k];
	}
	for(k=0; k<n_gen; k++) work[2][k] /= sw;

	/* M-step */
	for(k=0; k<n_gen; k++) {
	  work[1][k] += work[2][k];
	  if(pheno[j]) Result[k+1][i] += work[2][k];
	}
      }
      
      /* complete M-step */
      for(k=0; k<n_gen; k++) 
	Result[k+1][i] /= work[1][k];

      /* check for convergence */
      flag = 0;
      for(k=0; k<n_gen; k++) {
	if(fabs(Result[k+1][i] - work[0][k]) > tol*(fabs(work[0][k])+tol*100.0)) {
	  flag = 1;
	  break;
	}
      }

      if(!flag) break;
    } /* end of EM iterations */

    if(flag) warning("Didn't converge!\n");

    /* calculate log lik */
    Result[0][i] = 0.0;
    for(j=0; j<n_ind; j++) {
      sw = 0.0;
      if(pheno[j]) 
	for(k=0; k<n_gen; k++) sw += Genoprob[k][i][j] * Result[k+1][i];
      else
	for(k=0; k<n_gen; k++) sw += Genoprob[k][i][j] * (1.0-Result[k+1][i]);
      Result[0][i] += log10(sw);
    }

  } /* end loop over marker positions */
}



/* end of discan.c */

