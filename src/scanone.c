/**********************************************************************
 * 
 * scanone.c
 *
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 * Oct, 2001; Aug, 2001; May, 2001
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a genome scan with a 
 * single QTL model
 *
 * Contains: R_scanone_anova, scanone_anova, scanone_anova_perm
 *           R_scanone_im, scanone_im, scanone_im_perm
 *           R_scanone_hk, scanone_hk, scanone_hk_perm
 *           R_scanone_imp, scanone_imp, scanone_imp_perm
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
#include "scanone.h"
#define TOL 1e-10

/**********************************************************************
 * 
 * R_scanone_anova
 *
 * Wrapper for call from R; reorganizes genotype and result matrix
 * and calls scanone_anova.
 * 
 **********************************************************************/

void R_scanone_anova(int *n_ind, int *n_pos, int *n_gen,
		     int *geno, double *pheno,
		     double *result)
{
  int **Geno;
  double **Result;

  reorg_geno(*n_ind, *n_pos, geno, &Geno);
  reorg_errlod(*n_pos, *n_gen+2, result, &Result);

  scanone_anova(*n_ind, *n_pos, *n_gen, Geno, pheno, Result);
}

/**********************************************************************
 * 
 * R_scanone_im
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_im.
 * 
 **********************************************************************/

void R_scanone_im(int *n_ind, int *n_pos, int *n_gen, 
		  double *genoprob, double *pheno,
		  double *result, int *std_start, double *start,
		  int *maxit, double *tol)
{
  double ***Genoprob, **Result, **work;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  reorg_errlod(*n_pos, *n_gen+2, result, &Result);
  allocate_dmatrix(4,*n_gen, &work);

  /* Read R's random seed */
  GetRNGstate();

  scanone_im(*n_ind, *n_pos, *n_gen, Genoprob, pheno, Result, 
	     *std_start, start, *maxit, *tol, work);

  /* Write R's random seed */
  PutRNGstate();
}

/**********************************************************************
 * 
 * R_scanone_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_hk.
 * 
 **********************************************************************/

void R_scanone_hk(int *n_ind, int *n_pos, int *n_gen,
		  double *genoprob, double *pheno,
		  double *result)
{
  double ***Genoprob, **Result;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  reorg_errlod(*n_pos, *n_gen+2, result, &Result);

  scanone_hk(*n_ind, *n_pos, *n_gen, Genoprob, pheno, Result);
}


/**********************************************************************
 * 
 * R_scanone_imp
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_imp.
 * 
 **********************************************************************/

void R_scanone_imp(int *n_ind, int *n_pos, int *n_gen, int *n_draws,
		   int *draws, double *pheno, double *result)
{
  int ***Draws;

  reorg_draws(*n_ind, *n_pos, *n_draws, draws, &Draws);

  scanone_imp(*n_ind, *n_pos, *n_gen, *n_draws, Draws, pheno, result);
}


/**********************************************************************
 * 
 * scanone_anova
 *
 * Performs genotype scan using ANOVA
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
 * Result       Result matrix of size [n_pos x (n_gen+2)]; upon return, 
 *              the first column contains the RSS, the next set contain
 *              estimated genotype-specific means, and the last column
 *              contains the estimated residual SD
 *
 **********************************************************************/

void scanone_anova(int n_ind, int n_pos, int n_gen, int **Geno, 
		   double *pheno, double **Result)
{
  int i, j, k, n, *ng;
  double s1, s2, ss0, s0;

  /* number of individuals in each genotype group */
  ng = (int *)R_alloc(n_gen, sizeof(int));

  for(i=0; i<n_pos; i++) {
    Result[0][i] = 0.0; n=0; ss0=s0=0.0;
    for(j=0; j< n_gen; j++) {

      Result[j+1][i] = s2 = 0.0;
      ng[j] = 0;

      for(k=0; k<n_ind; k++) {
	if(Geno[i][k] == j+1) {
	  Result[j+1][i] += (s1 = pheno[k]);
	  s2 += (s1*s1);
	  n++; ng[j]++;
	}
      }
      /* overall sum of squares and sum of phenotypes*/
      ss0 += s2;
      s0 += Result[j+1][i];
      
      if(ng[j] > 0) {
	/* RSS for genotype group */
	Result[0][i] += s2 - (Result[j+1][i]*Result[j+1][i])/(double)ng[j];
	
	/* mean phenotype for genotype group */
	Result[j+1][i] = Result[j+1][i]/(double)ng[j];
      }
      else /* no individuals with this genotype */
	Result[j+1][i] = NA_REAL;
	
    } /* loop over genotype groups */

    /* overall RSS */
    ss0 = ss0 - (s0*s0)/(double)n;

    Result[n_gen+1][i] = sqrt(Result[0][i]/(double)(n - n_gen));

    Result[0][i] = (double)n*(log10(ss0) - log10(Result[0][i]))/2.0;

  } /* loop over marker positions */

}

/**********************************************************************
 * 
 * scanone_im
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

void scanone_im(int n_ind, int n_pos, int n_gen, double ***Genoprob,
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

      for(j=0; j<n_ind; j++) { /* loop over families */
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
	if(fabs(Result[k+1][i] - work[1][k]) > tol) {
	  flag = 1;
	  break;
	}
      }
      if(fabs(Result[n_gen+1][i] - oldsig) > tol) flag = 1;

      if(!flag) break;

      oldsig = Result[n_gen+1][i];
      for(k=0; k<n_gen; k++) work[1][k] = Result[1+k][i];

    } /* end of EM iterations */

    if(flag) warning("Didn't converge!\n");

    /* calculate log lik */
    Result[0][i] = 0.0;
    for(j=0; j<n_ind; j++) {
      s1 = 0.0;
      for(k=0; k<n_gen; k++) 
	s1 += Genoprob[k][i][j] * dnorm(pheno[j], Result[1+k][i], 
					Result[1+n_gen][i], 0);
      Result[0][i] += log10(s1);
    }
    Result[1+n_gen][i] *= sqrt((double)n_ind/(double)(n_ind-n_gen));

  } /* end loop over marker positions */
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
 *
 * pheno        Phenotype data, as a vector
 *
 * Result       Result matrix of size [n_pos x (n_gen+2)]; upon return, 
 *              the first column contains the RSS, the next set contain
 *              estimated genotype-specific means, and the last column
 *              contains the estimated residual SD
 *
 **********************************************************************/

void scanone_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		double *pheno, double **Result)
{
  int ny, *jpvt, k, i, j;
  double *work, *x, *qty, *qraux, *coef, *resid, tol;

  /* tolerance for linear regression */
  tol = TOL;

  /* allocate space and set things up*/
  jpvt = (int *)R_alloc(n_gen, sizeof(int));
  work = (double *)R_alloc(2 * n_gen, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  qraux = (double *)R_alloc(n_gen, sizeof(double));
  x = (double *)R_alloc(n_ind*n_gen, sizeof(double));
  coef = (double *)R_alloc(n_gen, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  ny = 1;
  for(i=0; i<n_gen; i++) jpvt[i] = i;

  for(i=0; i<n_pos; i++) { /* loop over positions */
    /* fill up X matrix with genotype probabilities */
    for(j=0; j<n_ind; j++) 
      for(k=0; k<n_gen; k++)
	x[j+k*n_ind] = Genoprob[k][i][j]; 

    /* linear regression of phenotype on QTL genotype probabilities */
    k=0;
    dqrls_(x, &n_ind, &n_gen, pheno, &ny, &tol, coef, resid,
	   qty, &k, jpvt, qraux, work);

    /* re-scramble coefficients */
    for(j=0; j<n_gen; j++)
      Result[1+j][i] = coef[jpvt[j]];
    
    /* RSS */
    Result[0][i] = 0.0;
    for(j=0; j<n_ind; j++)
      Result[0][i] += (resid[j]*resid[j]);

    /* residual SD */
    Result[1+n_gen][i] = sqrt(Result[0][i]/(double)(n_ind-n_gen));

    /* log10 likelihood */
    Result[0][i] = log10(Result[0][i]);

  } /* end loop over positions */

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
 * pheno        Phenotype data, as a vector
 *
 * result       Result vector of length [n_pos]; upon return, contains
 *              the "LPD" (log posterior distribution of QTL location).
 *
 **********************************************************************/

void scanone_imp(int n_ind, int n_pos, int n_gen, int n_draws,
		 int ***Draws, double *pheno, double *result)
{
  int ny, *jpvt, k, i, j, s;
  double *work, *x, *qty, *qraux, *coef, *resid, tol, rss, rss0;
  double sums, sumss;

  /* tolerance for linear regression */
  tol = TOL;

  /* allocate space and set things up*/
  jpvt = (int *)R_alloc(n_gen, sizeof(int));
  work = (double *)R_alloc(2 * n_gen, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  qraux = (double *)R_alloc(n_gen, sizeof(double));
  x = (double *)R_alloc(n_ind*n_gen, sizeof(double));
  coef = (double *)R_alloc(n_gen, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  ny = 1;

  /* two-pass algorithm to calculate null RSS */
  for(i=0, rss=0.0; i<n_ind; i++) 
    rss += pheno[i];
  rss /= (double)n_ind;
  for(i=0, rss0=0.0; i<n_ind; i++) 
    rss0 += (pheno[i]-rss)*(pheno[i]-rss);

  for(i=0; i<n_pos; i++) { /* loop over positions */

    sums = sumss = 0.0;
    for(j=0; j<n_draws; j++) { /* loop over imputations */

      /* fill up X matrix with genotype data */
      for(k=0; k<n_ind; k++) {
	for(s=0; s<n_gen; s++)
	  x[k+s*n_ind] = 0.0;
	x[k+n_ind*(Draws[j][i][k]-1)] = 1.0;
      }

      /* linear regression of phenotype on QTL genotype probabilities */
      k=0;
      dqrls_(x, &n_ind, &n_gen, pheno, &ny, &tol, coef, resid,
	     qty, &k, jpvt, qraux, work);

      /* calculate RSS */
      for(k=0, rss=0.0; k<n_ind; k++)
	rss += (resid[k]*resid[k]);

      rss = (double)n_ind/2.0*(log10(rss0)-log10(rss));
      sums += rss;
      sumss += (rss*rss);

    } /* end loop over imputations */

    result[i] = sumss/(double)(n_draws) +
      0.5*(sumss-sums*sums/(double)(n_draws))/(double)(n_draws-1);

  } /* end loop over positions */
}


/**********************************************************************
 * 
 * scanone_anova_perm
 * 
 * Do permutation test with anova
 * 
 * arguments are like that for R_scanone_anova, except result, 
 * which is a vector of length n_perm (to contain the output)
 *
 **********************************************************************/

void scanone_anova_perm(int *n_ind, int *n_pos, int *n_gen, int *geno,
			double *pheno, int *n_perm, double *result)
{
  int **Geno, i, j;
  double **scanres;

  /* Read R's random seed */
  GetRNGstate();

  allocate_dmatrix(*n_gen+2, *n_pos, &scanres);
  reorg_geno(*n_ind, *n_pos, geno, &Geno);

  for(i=0; i < *n_perm; i++) {
    double_permute(pheno, *n_ind); 

    scanone_anova(*n_ind, *n_pos, *n_gen, Geno, pheno, scanres);
  
    result[i] = scanres[0][0];
    for(j=0; j< *n_pos; j++)
      if(result[i] > scanres[0][j]) 
	result[i] = scanres[0][j];
  }

  /* write R's random seed */
  PutRNGstate();
}


void scanone_im_perm(int *n_ind, int *n_pos, int *n_gen, 
		     double *genoprob, double *pheno, 
		     int *n_perm, double *result,
		     int *std_start, double *start,
		     int *maxit, double *tol)
{
  int i, j;
  double ***Genoprob, **scanres, **work;
  
  /* Read R's random seed */
  GetRNGstate();

  allocate_dmatrix(*n_gen+2, *n_pos, &scanres);
  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  allocate_dmatrix(4,*n_gen, &work);

  for(i=0; i < *n_perm; i++) {
    double_permute(pheno, *n_ind); 

    scanone_im(*n_ind, *n_pos, *n_gen, Genoprob, pheno,
	       scanres, *std_start, start, *maxit, *tol, work);
  
    result[i] = scanres[0][0];
    for(j=0; j< *n_pos; j++)
      if(result[i] < scanres[0][j]) 
	result[i] = scanres[0][j];
  }

  /* write R's random seed */
  PutRNGstate();
}

void scanone_hk_perm(int *n_ind, int *n_pos, int *n_gen, 
		     double *genoprob, double *pheno, 
		     int *n_perm, double *result)
{
  int i, j;
  double ***Genoprob, **scanres;
  
  /* Read R's random seed */
  GetRNGstate();

  allocate_dmatrix(*n_gen+2, *n_pos, &scanres);
  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);

  for(i=0; i < *n_perm; i++) {
    double_permute(pheno, *n_ind); 

    scanone_hk(*n_ind, *n_pos, *n_gen, Genoprob, pheno,
	       scanres);
  
    result[i] = scanres[0][0];
    for(j=0; j< *n_pos; j++)
      if(result[i] > scanres[0][j]) 
	result[i] = scanres[0][j];
  }

  /* write R's random seed */
  PutRNGstate();
}



/* end of scanone.c */

