/**********************************************************************
 * 
 * fitqtl_imp.c
 *
 * copyright (c) 2002, Hao Wu, The Jackson Laboratory
 *
 * last modified June, 2002
 * first written Apr, 2002
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a general genome scan 
 *
 * Contains: R_fitqtl_imp, fitqtl_imp, nullRss0, galtRss
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
#include "fitqtl_imp.h"
#define TOL 1e-12

void R_fitqtl_imp(int *n_ind, int *n_qtl, int *n_gen, int *n_draws,
		  int *draws, int *n_cov, double *cov, int *model, 
		  int *n_int, double *pheno,
		  /* return variables */
		  double *lod, int *df)
{
  /* reorganize draws */
  int ***Draws;
  double **Cov; 

  reorg_draws(*n_ind, *n_qtl, *n_draws, draws, &Draws);

  /* reorganize cov (if they are not empty) */
  /* currently reorg_errlod function is used to reorganize the data */
  if(*n_cov != 0) reorg_errlod(*n_ind, *n_cov, cov, &Cov);

  fitqtl_imp(*n_ind, *n_qtl, n_gen, *n_draws, Draws,
              Cov, *n_cov, model, *n_int, pheno, lod, df);
}


/**********************************************************************
 * 
 * fitqlt_imp
 *
 * Performs general genotype scan using the pseudomarker algorithm 
 * (imputation) method of Sen and Churchill (2001).
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
 * lod          Return LOD score
 *
 * df           Return degree of freedom
 *
 **********************************************************************/

void fitqtl_imp(int n_ind, int n_qtl, int *n_gen, int n_draws, 
		int ***Draws, double **Cov, int n_cov, 
		int *model, int n_int, double *pheno, double *lod, int *df) 
{

  /* create local variables */
  int i, j, n_qc, itmp; /* loop variants and temp variables */
  double tol, lrss, lrss0, *LOD_array;
  double *dwork;
  int *iwork, sizefull;

  /* initialization */
  sizefull = 1;
  tol = TOL;

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

  /* allocate memory for working arrays, total memory is
     sizefull*n_ind+2*n_ind+4*sizefull for double array, 
     and sizefull for integer array.
     All memory will be allocated one time and split later */
  dwork = (double *)R_alloc(sizefull*n_ind+2*n_ind+4*sizefull,
			    sizeof(double));
  iwork = (int *)R_alloc(sizefull, sizeof(int));
  LOD_array = (double *)R_alloc(n_draws, sizeof(double));


  /* calculate null model RSS */
  lrss0 = log10(nullRss0(pheno, n_ind));

  /* loop over imputations */
  for(i=0; i<n_draws; i++) {
    /* calculate alternative model RSS */
    lrss = log10( galtRss(pheno, n_ind, n_gen, n_qtl, Draws[i], 
			  Cov, n_cov, model, n_int, dwork, iwork, sizefull) );

    /* calculate the LOD score in this imputation */
    LOD_array[i] = (double)n_ind/2.0*(lrss0-lrss);
  } /* end loop over imputations */

  /* calculate the result LOD score */
  *lod = wtaverage(LOD_array, n_draws);

  /* degree of freedom equals to the number of columns of x minus 1 (mean) */
  *df = sizefull - 1;

}



/* function to calculate the null model RSS. This function is different
   from the function used in scanone_imp and scantwo_imp, which contain 
   the additive covariate in the null model */
double nullRss0(double *pheno, int n_ind)
{
  /* local variables */
  int i;
  double s, rss, tol, m;
  tol = TOL;

  /* calculate the mean of phenotype */
  s = 0.0;
  for(i=0; i<n_ind; i++)
    s += pheno[i];
  m = s / (double)n_ind;

  /* calculate RSS */
  rss = 0.0;
  for(i=0; i<n_ind; i++)
    rss += (pheno[i]-m)*(pheno[i]-m);

  return(rss);
}


/* galtRss - calculate RSS for full model in general scan */
double galtRss(double *pheno, int n_ind, int *n_gen, int n_qtl, 
	       int **Draws, double **Cov, int n_cov, int *model, 
	       int n_int, double *dwork, int *iwork, int sizefull) 
{
  /* local variables */
  int i, j, k, kk, *jpvt, ny, idx_col, n_qc, itmp1, itmp2, n_int_col, tmp_idx;
  double *work, *x, *qty, *qraux, *coef, *resid, tol;
  /* The following vars are used for parsing model. 
     the dimension of idx_int_q and idx_int_c are set to be arbitrary number
     for the ease of programming. But I think 15 and 10 are big enought. 
     Is there any body want to try a 16 way interaction? */
  int n_int_q, n_int_c, idx_int_q[15], idx_int_c[10]; 
  /* return variable */
  double rss_full;

  /* initialization */
  ny = 1; 
  idx_col = 0;
  rss_full = 0.0;
  tol = TOL;
  n_qc = n_qtl + n_cov;

  /* split the memory block: 
     design matrix x will be (n_ind x sizefull), coef will be (1 x sizefull),
     resid will be (1 x n_ind), qty will be (1 x n_ind), 
     qraux will be (1 x sizefull), work will be (2 x sizefull) */
  x = dwork;
  coef = x + n_ind*sizefull;
  resid = coef + sizefull;
  qty = resid + n_ind;
  qraux = qty + n_ind;
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
	       all entrys for this individual will be 0.
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


  /* call dqrls to fit regression model */
  dqrls_(x, &n_ind, &sizefull, pheno, &ny, &tol, coef, resid,
         qty, &k, jpvt, qraux, work);

  /* calculate RSS */
  for(i=0; i<n_ind; i++)
    rss_full += resid[i]*resid[i];

  return(rss_full);
}

/* end of fitqtl_imp.c */
