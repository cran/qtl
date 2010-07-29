/**********************************************************************
 * 
 * fitqtl_imp_binary.c
 *
 * copyright (c) 2010, Karl W. Broman
 *
 * last modified Jun, 2010
 * first written Jun, 2010
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 * 
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 * C functions for the R/qtl package
 *
 * These functions are for fitting a fixed multiple-QTL model by 
 * imputation, for a binary trait.
 *
 * Contains: R_fitqtl_imp_binary, fitqtl_imp_binary, galtLODimpbin
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "fitqtl_imp_binary.h"
#include "fitqtl_hk_binary.h"
#define TOL 1e-12
#define IDXINTQ 15  /* maximum no. QTLs in an interaction */
#define IDXINTC 10  /* maximum no. covariates in an interaction */

void R_fitqtl_imp_binary(int *n_ind, int *n_qtl, int *n_gen, int *n_draws,
			 int *draws, int *n_cov, double *cov, int *model, 
			 int *n_int, double *pheno, int *get_ests,
			 /* return variables */
			 double *lod, int *df, double *ests, double *ests_covar,
			 double *design_mat,
			 /* convergence */
			 double *tol, int *maxit)
{
  int ***Draws;
  double **Cov;

  /* reorganize draws */
  reorg_draws(*n_ind, *n_qtl, *n_draws, draws, &Draws);

  /* reorganize cov (if they are not empty) */
  /* currently reorg_errlod function is used to reorganize the data */
  if(*n_cov != 0) reorg_errlod(*n_ind, *n_cov, cov, &Cov);

  fitqtl_imp_binary(*n_ind, *n_qtl, n_gen, *n_draws, Draws,
		    Cov, *n_cov, model, *n_int, pheno, *get_ests, lod, df,
		    ests, ests_covar, design_mat, *tol, *maxit);
}


/**********************************************************************
 * 
 * fitqtl_imp_binary
 *
 * Fits a fixed multiple-QTL model by multiple imputation.
 * 
 * n_ind        Number of individuals
 *
 * n_qtl        Number of QTLs in the model 
 *
 * n_gen        Number of different genotypes
 *
 * n_draws      Number of impiutations
 *
 * Draws        Array of genotype imputations, indexed as 
 *              Draws[draw][mar][ind]
 *
 * Cov          covariates matrix, Cov[mar][ind]
 *
 * n_cov        Number of covariates
 *
 * model        Model matrix
 *
 * n_int        Number of interactions in the model
 *
 * pheno        Phenotype data, as a vector
 *
 * get_ests     0/1: If 1, return estimated effects and their variances
 *
 * lod          Return LOD score
 *
 * df           Return degree of freedom
 *
 * ests         Return ests (vector of length sizefull)
 *
 * ests_covar   Return covariance matrix of ests (sizefull^2 matrix)
 *
 * tol          Tolerance for convergence
 * 
 * maxit        Maximum number of iterations in IRLS
 *
 **********************************************************************/

void fitqtl_imp_binary(int n_ind, int n_qtl, int *n_gen, int n_draws, 
		       int ***Draws, double **Cov, int n_cov, 
		       int *model, int n_int, double *pheno, int get_ests,
		       double *lod, int *df, double *ests, double *ests_covar,
		       double *design_mat, double tol, int maxit) 
{

  /* create local variables */
  int i, j, ii, jj, n_qc, itmp; /* loop variants and temp variables */
  double llik, llik0, *LOD_array;
  double *the_ests, *the_covar, **TheEsts, ***TheCovar;
  double *dwork, **Ests_covar, tot_wt=0.0, *wts;
  double **Covar_mean, **Mean_covar, *mean_ests; /* for ests and cov matrix */
  int *iwork, sizefull, n_trim, *index;

  /* number to trim from each end of the imputations */
  n_trim = (int) floor( 0.5*log(n_draws)/log(2.0) );

  /* initialization */
  sizefull = 1;

  /* calculate the dimension of the design matrix for full model */
  n_qc = n_qtl+n_cov; /* total number of QTLs and covariates */
  /* for additive QTLs and covariates*/
  for(i=0; i<n_qc; i++)
    sizefull += n_gen[i];
  /* for interactions, loop thru all interactions */
  for(i=0; i<n_int; i++) { 
    for(j=0, itmp=1; j<n_qc; j++) {
      if(model[i*n_qc+j])
	itmp *= n_gen[j];
    }
    sizefull += itmp; 
  }

  /* reorganize Ests_covar for easy use later */
  /* and make space for estimates and covariance matrix */
  if(get_ests) {
    reorg_errlod(sizefull, sizefull, ests_covar, &Ests_covar);

    allocate_double(sizefull*n_draws, &the_ests);
    allocate_double(sizefull*sizefull*n_draws, &the_covar);

    /* I need to save all of the estimates and covariance matrices */
    reorg_errlod(sizefull, n_draws, the_ests, &TheEsts);
    reorg_genoprob(sizefull, sizefull, n_draws, the_covar, &TheCovar);

    allocate_dmatrix(sizefull, sizefull, &Mean_covar);
    allocate_dmatrix(sizefull, sizefull, &Covar_mean);
    allocate_double(sizefull, &mean_ests);
    allocate_double(n_draws, &wts);
  }

  /* allocate memory for working arrays, total memory is
     sizefull*n_ind+6*n_ind+4*sizefull for double array, 
     and sizefull for integer array.
     All memory will be allocated one time and split later */
  dwork = (double *)R_alloc(sizefull*n_ind+6*n_ind+4*sizefull,
			    sizeof(double));
  iwork = (int *)R_alloc(sizefull, sizeof(int));
  index = (int *)R_alloc(n_draws, sizeof(int));
  LOD_array = (double *)R_alloc(n_draws, sizeof(double));


  /* calculate null model log10 likelihood */
  llik0 = nullLODbin(pheno, n_ind);

  /* loop over imputations */
  for(i=0; i<n_draws; i++) {

    R_CheckUserInterrupt(); /* check for ^C */

    /* calculate alternative model RSS */
    llik = galtLODimpbin(pheno, n_ind, n_gen, n_qtl, Draws[i], 
			 Cov, n_cov, model, n_int, dwork, iwork, sizefull,
			 get_ests, ests, Ests_covar, design_mat,
			 tol, maxit);

    /* calculate the LOD score in this imputation */
    LOD_array[i] = (llik - llik0);

    /* if getting estimates, calculate the weights */
    if(get_ests) { 
      wts[i] = LOD_array[i]*log(10.0);
      if(i==0) tot_wt = wts[i];
      else tot_wt = addlog(tot_wt, wts[i]);
      
      for(ii=0; ii<sizefull; ii++) {
	TheEsts[i][ii] = ests[ii];
	for(jj=ii; jj<sizefull; jj++) 
	  TheCovar[i][ii][jj] = Ests_covar[ii][jj];
      }
    }

  } /* end loop over imputations */

  /* sort the lod scores, and trim the weights */
  if(get_ests) {
    for(i=0; i<n_draws; i++) {
      index[i] = i;
      wts[i] = exp(wts[i]-tot_wt);
    }

    rsort_with_index(LOD_array, index, n_draws);

    for(i=0; i<n_trim; i++) 
      wts[index[i]] = wts[index[n_draws-i-1]] = 0.0;

    /* re-scale wts */
    tot_wt = 0.0;
    for(i=0; i<n_draws; i++) tot_wt += wts[i];
    for(i=0; i<n_draws; i++) wts[i] /= tot_wt;
  } 

  /* calculate the result LOD score */
  *lod = wtaverage(LOD_array, n_draws);

  /* degree of freedom equals to the number of columns of x minus 1 (mean) */
  *df = sizefull - 1;

  /* get means and variances and covariances of estimates */
  if(get_ests) { 
    for(i=0; i<n_draws; i++) {
      if(i==0) {
	for(ii=0; ii<sizefull; ii++) {
	  mean_ests[ii] = TheEsts[i][ii] * wts[i];
	  for(jj=ii; jj<sizefull; jj++) {
	    Mean_covar[ii][jj] = TheCovar[i][ii][jj] * wts[i];
	    Covar_mean[ii][jj] = 0.0;
	  }
	}
      }
      else {
	for(ii=0; ii<sizefull; ii++) {
	  mean_ests[ii] += TheEsts[i][ii]*wts[i];
	  for(jj=ii; jj<sizefull; jj++) {
	    Mean_covar[ii][jj] += TheCovar[i][ii][jj]*wts[i];
	    Covar_mean[ii][jj] += (TheEsts[i][ii]-TheEsts[0][ii])*
	      (TheEsts[i][jj]-TheEsts[0][jj])*wts[i];
	  }
	}
      }
    }
      
    for(i=0; i<sizefull; i++) {
      ests[i] = mean_ests[i];
      for(j=i; j<sizefull; j++) {
	Covar_mean[i][j] = (Covar_mean[i][j] - (mean_ests[i]-TheEsts[0][i])*
			    (mean_ests[j]-TheEsts[0][j]))*(double)n_draws/(double)(n_draws-1);
	Ests_covar[i][j] = Ests_covar[j][i] = Mean_covar[i][j] + Covar_mean[i][j];
      }
    }
  } /* done getting estimates */

}



/* galtRss - calculate RSS for full model by multiple imputation */
double galtLODimpbin(double *pheno, int n_ind, int *n_gen, int n_qtl, 
		     int **Draws, double **Cov, int n_cov, int *model, 
		     int n_int, double *dwork, int *iwork, int sizefull,
		     int get_ests, double *ests, double **Ests_covar,
		     double *designmat, double tol, int maxit) 
{
  /* local variables */
  int i, j, k, kk, s, *jpvt, ny, idx_col, n_qc, itmp1, itmp2, n_int_col, tmp_idx, job, flag;
  double *work, *x, *qty, *qraux, *coef, *resid, tol2;
  double *nu, *pi, *z, *wt;
  /* The following vars are used for parsing model. 
     the dimension of idx_int_q and idx_int_c are set to be arbitrary number
     for the ease of programming. But I think 15 and 10 are big enough. 
     Is there any body want to try a 16 way interaction? */
  int n_int_q, n_int_c, idx_int_q[IDXINTQ], idx_int_c[IDXINTC]; 
  /* return variable */
  double curllik, llik=0.0;

  /* initialization */
  ny = 1; 
  idx_col = 0;
  tol2 = TOL;
  n_qc = n_qtl + n_cov;

  /* split the memory block: 
     design matrix x will be (n_ind x sizefull), coef will be (1 x sizefull),
     resid will be (1 x n_ind), qty will be (1 x n_ind), 
     qraux will be (1 x sizefull), work will be (2 x sizefull) */
  x = dwork;
  coef = x + n_ind*sizefull;
  resid = coef + sizefull;
  qty = resid + n_ind;
  pi = qty + n_ind;
  z = pi + n_ind;
  nu = z + n_ind;
  wt = nu + n_ind;
  qraux = wt + n_ind;
  work = qraux + sizefull; 
  /* integer array */
  jpvt = iwork;

  /******************************************************
   The following part will construct the design matrix x 
   ******************************************************/
  /* fill first row with 1s. It's corresponding to the mean */
  for(i=0; i<n_ind; i++) x[i] = 1.0;
  idx_col += 1;  /* increment column index */

  /* zero out the rest of x */
  for(i=idx_col*n_ind; i<n_ind*sizefull; i++) x[i] = 0.0;

  /*****************
   * Additive terms 
   *****************/
  /* loop thru QTLs */
  /* if the geno type is one, do nothing (the effects go to the means).
     Otherwise, set proper entry in x to be 1. The idea is that for 
     backcross, genotype 1 -> 0; 2 -> 1. For F2, genotype 1  -> [0 0]; 
     2 -> [1 0]; 3 ->[0 1]. For 4-way, 1 -> [0 0 0], 2 -> [1 0 0],
     3 -> [0 1 0], 4 -> [0 0 1] and so on */
  for(i=0; i<n_qtl; i++) {
    for(j=0; j<n_ind; j++) { /*loop thru individuals */
      if(Draws[i][j] != 1) 
	x[(idx_col+Draws[i][j]-2)*n_ind + j] = 1.0;
    } /* finish the current QTL */
    /* increment idx_col */
    idx_col += n_gen[i];
  }

  /* loop thru covariates */
  for(i=0; i<n_cov; i++) {
    for(j=0; j<n_ind; j++)  /* loop individuals */
      x[idx_col*n_ind + j] = Cov[i][j];
    idx_col ++; /* increment idx_col by 1 */
  }

  /*******************
   * interactive terms 
   *******************/
  /* loop thru interactions */
  for(i=0; i<n_int; i++) {
    n_int_q = 0;
    n_int_c = 0;
   /* total number of columns in the design matrix for this interaction */
    n_int_col = 1; 
    /* parse the model matrix */
    for(j=0; j<n_qtl; j++) { 
      if(model[i*n_qc+j]) { /* this QTL is in the interaction */
	idx_int_q[n_int_q] = j;
	n_int_q ++;
	n_int_col *= n_gen[j]; 
      }
    }
    for(j=n_qtl; j<n_qc; j++) {
      if(model[i*n_qc+j]) { /* this covariate is in the interaction */
	idx_int_c[n_int_c] = j - n_qtl;
	n_int_c ++; 
	/* note that n_gen for covariate are always one (one column for a covariate),
	   so there's no n_int_col *= n_gen[j] here */
      }
    }

    /* construct the design matrix for this interaction */
    for(j=0; j<n_ind; j++) { /* loop thru the individuals */
      if( n_int_q ) { /* there IS QTL(s) in interaction */
	itmp2 = 1;
	for(k=0; k<n_int_q; k++) {
	  itmp1 = idx_int_q[k]; /* index for QTL in the interaction */
	  if(Draws[itmp1][j] == 1) {
	    /* if any genotype in this interaction is 1, 
	       all entries for this individual will be 0.
	       break the loop and go to next individual */
	    itmp2 = 0; /* a flag to indicate some genotype is 1 */
	    break;
	  }
	}
	/* if itmp2 is 1 (all genotypes are not 1), set corresponding 
	   entry in design matrix to be 1 */
	if( itmp2 ) {
	  /* find the proper column index */
	  itmp1 = idx_int_q[n_int_q-1];
	  tmp_idx = Draws[itmp1][j] - 2;
	  /* reuse itmp2 */
	  itmp2 = n_gen[idx_int_q[n_int_q-1]];
	  for(k=n_int_q-2; k>=0; k--) {
	    itmp1 = idx_int_q[k];
	    tmp_idx += (Draws[itmp1][j]-2)*itmp2;
	    itmp2 *= n_gen[idx_int_q[k]];
	  }
	  if(tmp_idx != 0) 
	    kk = 0;
	  x[(idx_col+tmp_idx)*n_ind+j] = 1;
	  /* interaction with covariates */
	  for(k=0; k<n_int_c; k++)
	    x[(idx_col+tmp_idx)*n_ind+j] *= Cov[idx_int_c[k]][j];
	} /* finish this interaction (with QTL case) */
      }
      else { 
	/* there's NO QTL in interaction (the interaction for covariates).
	   can this happen? I'll put it here anyway */
	x[idx_col*n_ind+j] = 1; /* this entry is the product of all convariates */
	for(k=0; k<n_int_c; k++)
	  x[idx_col*n_ind+j] *= Cov[idx_int_c[k]][j];
      } /* finish this interaction (without QTL case) */
    } /* finish the loop for individuals */

    idx_col += n_int_col;
  } /* finish the loop for interaction */
  /* finish design matrix construction */

  /* save design matrix */
  memcpy(designmat, x, n_ind*sizefull*sizeof(double));

  /* starting point for IRLS */
  curllik = 0.0;
  for(j=0; j<n_ind; j++) {
    pi[j] = (pheno[j] + 0.5)/2.0;
    wt[j] = sqrt(pi[j] * (1.0-pi[j]));
    nu[j] = log(pi[j]) - log(1.0-pi[j]);
    z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
    curllik += pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]);
  }

  /* multiply design matrix by current wts */
  for(i=0; i<sizefull; i++) 
    for(j=0; j<n_ind; j++)
      x[i*n_ind+j] *= wt[j];

  flag = 0;
  for(s=0; s<maxit; s++) { /* IRLS iterations */
  
    R_CheckUserInterrupt(); /* check for ^C */

    /* make jpvt = numbers 0, 1, ..., (sizefull-1) */
    /*      jpvt keeps track of any permutation of X columns */
    for(i=0; i<sizefull; i++) jpvt[i] = i;

    /* call dqrls to fit regression model */
    F77_CALL(dqrls)(x, &n_ind, &sizefull, z, &ny, &tol2, coef, resid,
  		  qty, &kk, jpvt, qraux, work);

    /* get ests; need to permute back */
    for(i=0; i<kk; i++) ests[jpvt[i]] = coef[i];
    for(i=kk; i<sizefull; i++) ests[jpvt[i]] = 0.0;
    
    /* re-form design matrix */
    memcpy(x, designmat, n_ind*sizefull*sizeof(double));

    /* calculate fitted values, probs, new wts, new z's */
    llik = 0.0;
    for(j=0; j<n_ind; j++) {
      nu[j] = 0.0;
      for(i=0; i<sizefull; i++) 
	nu[j] += x[i*n_ind+j] * ests[i];
      pi[j] = exp(nu[j]);
      pi[j] /= (1.0 + pi[j]);
      wt[j] = sqrt(pi[j] * (1.0-pi[j]));
      z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
      llik += (pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]));

      /* multiply design matrix by new weights */
      for(i=0; i<sizefull; i++) 
	x[i*n_ind+j] *= wt[j];
    }

    if(fabs(llik - curllik) < tol) { /* converged? */
      flag = 1;
      break;
    }
    curllik = llik;

  } /* end of IRLS iterations */

  if(!flag)
    warning("Didn't converge.");

  if(get_ests) { /* get the estimates and their covariance matrix */

    /* need to re-run the last regression */
    F77_CALL(dqrls)(x, &n_ind, &sizefull, z, &ny, &tol2, coef, resid,
		    qty, &kk, jpvt, qraux, work);

    /* get covariance matrix: dpodi to get (X'X)^-1; re-sort; multiply by sigma_hat^2 */
    job = 1; /* indicates to dpodi to get inverse and not determinant */
    F77_CALL(dpodi)(x, &n_ind, &sizefull, work, &job);

    for(i=0; i<kk; i++) 
      for(j=i; j<kk; j++) 
	Ests_covar[jpvt[i]][jpvt[j]] = Ests_covar[jpvt[j]][jpvt[i]] = 
	  x[j*n_ind+i];
    for(i=kk; i<sizefull; i++)
      for(j=0; j<kk; j++)
	Ests_covar[jpvt[i]][j] = Ests_covar[j][jpvt[i]] = 0.0;
  }

  return(llik);
}

/* end of fitqtl_imp_binary.c */
