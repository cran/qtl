/**********************************************************************
 * 
 * hmm_main.c
 *
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 * Aug, 2001; Feb, 2001
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for the main HMM engine
 *
 * Contains: calc_genoprob, sim_geno, est_map, argmax_geno
 *           calc_errorlod, est_rf
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "hmm_main.h"
#include "util.h"

/**********************************************************************
 * 
 * calc_genoprob
 *
 * This function uses the hidden Markov model technology to calculate 
 * the genotype probabilities at each of marker and (optionally) at 
 * points between markers, conditional on all marker data for a 
 * chromosome.  This assumes data on a single chromosome
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              calculate the genotype probabilities)
 *
 * n_gen        Number of different genotypes
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps (may be ignored)
 *
 * error_prob   Genotyping error probability
 *
 * genoprob     Genotype probabilities (the output); a single vector
 *              stored by columns (ind moves fastest, then mar, then
 *              genotype
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void calc_genoprob(int n_ind, int n_pos, int n_gen, int *geno, 
		   double *rf, double *rf2, 
		   double error_prob, double *genoprob, 
		   double initf(int), 
		   double emitf(int, int, double),
		   double stepf(int, int, double, double)) 
{
  int i, j, j2, v, v2;
  double s, **alpha, **beta;
  int **Geno;
  double ***Genoprob;
  
  /* allocate space for alpha and beta and 
     reorganize geno and genoprob */
  reorg_geno(n_ind, n_pos, geno, &Geno);
  reorg_genoprob(n_ind, n_pos, n_gen, genoprob, &Genoprob);
  allocate_alpha(n_pos, n_gen, &alpha);
  allocate_alpha(n_pos, n_gen, &beta);

  for(i=0; i<n_ind; i++) { /* i = individual */

    /* initialize alpha and beta */
    for(v=0; v<n_gen; v++) {
      alpha[v][0] = initf(v+1) + emitf(Geno[0][i], v+1, error_prob);
      beta[v][n_pos-1] = 0.0;
    }

    /* forward-backward equations */
    for(j=1,j2=n_pos-2; j<n_pos; j++, j2--) {
      
      for(v=0; v<n_gen; v++) {
	alpha[v][j] = alpha[0][j-1] + stepf(1, v+1, rf[j-1], rf2[j-1]);
	
	beta[v][j2] = beta[0][j2+1] + stepf(v+1,1,rf[j2], rf2[j2]) + 
	  emitf(Geno[j2+1][i],1,error_prob);

	for(v2=1; v2<n_gen; v2++) {
	  alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
			       stepf(v2+1,v+1,rf[j-1],rf2[j-1]));
	  beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
			       stepf(v+1,v2+1,rf[j2],rf2[j2]) +
			       emitf(Geno[j2+1][i],v2+1,error_prob));
	}

	alpha[v][j] += emitf(Geno[j][i],v+1,error_prob);
      }
    }

    /* calculate genotype probabilities */
    for(j=0; j<n_pos; j++) {
      s = 0.0;
      for(v=0; v<n_gen; v++) 
	s += (Genoprob[v][j][i] = exp(alpha[v][j] + beta[v][j]));

      for(v=0; v<n_gen; v++) 
	Genoprob[v][j][i] /= s;

    }


  } /* loop over individuals */
  

}



/**********************************************************************
 * 
 * sim_geno
 *
 * This function simulates from the joint distribution Pr(g | O)
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              simulate genotypes)
 *
 * n_gen        Number of different genotypes
 *
 * n_draws      Number of simulation replicates
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps
 *
 * error_prob   Genotyping error probability
 *
 * draws        Simulated genotypes (the output), a single vector
 *              stored by columns (ind moves fastest, then mar, then
 *              draws
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void sim_geno(int n_ind, int n_pos, int n_gen, int n_draws,
	      int *geno, double *rf, double *rf2, 
	      double error_prob, int *draws,
	      double initf(int), 
	      double emitf(int, int, double),
	      double stepf(int, int, double, double)) 
{
  int i, k, j, v, v2;
  double s, **beta, *probs;
  int **Geno, ***Draws, curstate;
  
  /* allocate space for beta and 
     reorganize geno and draws */
  reorg_geno(n_ind, n_pos, geno, &Geno);
  reorg_draws(n_ind, n_pos, n_draws, draws, &Draws);
  allocate_alpha(n_pos, n_gen, &beta);
  allocate_double(n_gen, &probs);

  /* Read R's random seed */
  GetRNGstate();

  for(i=0; i<n_ind; i++) { /* i = individual */

    /* do backward equations */
    /* initialize beta */
    for(v=0; v<n_gen; v++) 
      beta[v][n_pos-1] = 0.0;

    /* backward equations */
    for(j=n_pos-2; j>=0; j--) {
      
      for(v=0; v<n_gen; v++) {
	beta[v][j] = beta[0][j+1] + stepf(v+1,1,rf[j], rf2[j]) + 
	  emitf(Geno[j+1][i],1,error_prob);

	for(v2=1; v2<n_gen; v2++) 
	  beta[v][j] = addlog(beta[v][j], beta[v2][j+1] + 
			      stepf(v+1,v2+1,rf[j],rf2[j]) +
			      emitf(Geno[j+1][i],v2+1,error_prob));

      }
    }

    for(k=0; k<n_draws; k++) { /* k = simulate replicate */

      /* first draw */
      /* calculate probs */
      for(v=0, s=0; v<n_gen; v++) 
	s += (probs[v] = exp(initf(v+1) + emitf(Geno[0][i], v+1, error_prob) +
			     beta[v][0]));
      for(v=0; v<n_gen; v++)
	probs[v] /= s;

      /* make draw */
      curstate = Draws[k][0][i] = sample_int(n_gen, probs);
      
      /* move along chromosome */
      for(j=1; j<n_pos; j++) {
	/* calculate probs */
	for(v=0, s=0; v<n_gen; v++) 
	  probs[v] = exp(stepf(curstate,v+1,rf[j-1],rf2[j-1]) +
			 emitf(Geno[j][i],v+1,error_prob) +
			 beta[v][j] - beta[curstate-1][j-1]);
	/* make draw */
	curstate = Draws[k][j][i] = sample_int(n_gen, probs);

      }

    } /* loop over replicates */

  } /* loop over individuals */
  
  /* write R's random seed */
  PutRNGstate();

}



/**********************************************************************
 * 
 * est_map
 *
 * This function re-estimates the genetic map for a chromosome
 *
 * n_ind        Number of individuals
 *
 * n_mar        Number of markers 
 *
 * n_gen        Number of different genotypes
 *
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          Second set of recombination fractions (may not be needed)
 *
 * error_prob   Genotyping error probability
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 * nrecf1       Function returning number of recombinations associated
 *              with (g_1, g_2)
 *
 * nrecf2       Another such function, used only in the case of a sex-
 *              specific map
 *
 * loglik       Loglik at final estimates of recombination fractions
 *
 * maxit        Maximum number of iterations to perform
 * 
 * tol          Tolerance for determining convergence
 * 
 * sexsp        Indicates whether sex-specific maps should be estimated
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void est_map(int n_ind, int n_mar, int n_gen, int *geno, double *rf, 
	     double *rf2, double error_prob, double initf(int), 
	     double emitf(int, int, double),
	     double stepf(int, int, double, double), 
	     double nrecf1(int, int), double nrecf2(int, int), 
	     double *loglik, int maxit, double tol, int sexsp, 
	     int prnt)
{
  int i, j, j2, v, v2, it, flag=0, **Geno;
  double s, **alpha, **beta, **gamma, *cur_rf, *cur_rf2;
  double curloglik;
  
  /* allocate space for beta and reorganize geno */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  allocate_alpha(n_mar, n_gen, &alpha);
  allocate_alpha(n_mar, n_gen, &beta);
  allocate_dmatrix(n_gen, n_gen, &gamma);
  allocate_double(n_mar-1, &cur_rf);
  allocate_double(n_mar-1, &cur_rf2);

  if(prnt) {
    /* print initial estimates */
    Rprintf("      "); 
    for(j=0; j<n_mar-1; j++) Rprintf("%.3lf ", rf[j]);
    Rprintf("\n");
  }

  /* begin EM algorithm */
  for(it=0; it<maxit; it++) {

    for(j=0; j<n_mar-1; j++) {
      cur_rf[j] = cur_rf2[j] = rf[j];
      rf[j] = 0.0;
      if(sexsp) {
	cur_rf2[j] = rf2[j];
	rf2[j] = 0.0;
      }
    }

    for(i=0; i<n_ind; i++) { /* i = individual */

      /* initialize alpha and beta */
      for(v=0; v<n_gen; v++) {
	alpha[v][0] = initf(v+1) + emitf(Geno[0][i], v+1, error_prob);
	beta[v][n_mar-1] = 0.0;
      }

      /* forward-backward equations */
      for(j=1,j2=n_mar-2; j<n_mar; j++, j2--) {
	
	for(v=0; v<n_gen; v++) {
	  alpha[v][j] = alpha[0][j-1] + stepf(1, v+1, cur_rf[j-1], cur_rf2[j-1]);
	  
	  beta[v][j2] = beta[0][j2+1] + stepf(v+1,1,cur_rf[j2], cur_rf2[j2]) + 
	    emitf(Geno[j2+1][i],1,error_prob);
	  
	  for(v2=1; v2<n_gen; v2++) {
	    alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
				 stepf(v2+1,v+1,cur_rf[j-1],cur_rf2[j-1]));
	    beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
				 stepf(v+1,v2+1,cur_rf[j2],cur_rf2[j2]) +
				 emitf(Geno[j2+1][i],v2+1,error_prob));
	  }
	  
	  alpha[v][j] += emitf(Geno[j][i],v+1,error_prob);
		 
	}

      }

      for(j=0; j<n_mar-1; j++) {

	/* calculate gamma = log Pr(v1, v2, O) */
	for(v=0, s=0.0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    gamma[v][v2] = alpha[v][j] + beta[v2][j+1] + 
	      emitf(Geno[j+1][i], v2+1, error_prob) +
	      stepf(v+1, v2+1, cur_rf[j], cur_rf2[j]);

	    if(v==0 && v2==0) s = gamma[v][v2];
	    else s = addlog(s, gamma[v][v2]);
	  }
	}

	for(v=0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    rf[j] += nrecf1(v+1,v2+1) * exp(gamma[v][v2] - s);
	    if(sexsp) rf2[j] += nrecf2(v+1,v2+1) * exp(gamma[v][v2] - s);
	  }
	}
      }

    } /* loop over individuals */

    /* rescale */
    for(j=0; j<n_mar-1; j++) {
      rf[j] /= (double)n_ind;
      if(rf[j] < tol) rf[j] = tol;
      else if(rf[j] > 0.5-tol) rf[j] = 0.5-tol;
      
      if(sexsp) {
	rf2[j] /= (double)n_ind;
	if(rf2[j] < tol) rf2[j] = tol;
	else if(rf2[j] > 0.5-tol) rf2[j] = 0.5-tol;
      }
    }

    /* check convergence */
    for(j=0, flag=0; j<n_mar-1; j++) {
      if(fabs(rf[j] - cur_rf[j]) > tol || 
	 (sexsp && fabs(rf2[j] - cur_rf2[j]) > tol)) {
	flag = 1; 
	break;
      }
    }

    if(!flag) break;

  } /* end EM algorithm */
  
  if(flag) warning("Didn't converge!\n");

  /* calculate log likelihood */
  *loglik = 0.0;
  for(i=0; i<n_ind; i++) { /* i = individual */
    /* initialize alpha and beta */
    for(v=0; v<n_gen; v++) 
      alpha[v][0] = initf(v+1) + emitf(Geno[0][i], v+1, error_prob);
    /* forward-backward equations */
    for(j=1; j<n_mar; j++) {
      for(v=0; v<n_gen; v++) {
	alpha[v][j] = alpha[0][j-1] + 
	  stepf(1, v+1, rf[j-1], rf2[j-1]);
	for(v2=1; v2<n_gen; v2++) 
	  alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
			       stepf(v2+1,v+1,rf[j-1],rf2[j-1]));
	alpha[v][j] += emitf(Geno[j][i],v+1,error_prob);
      }
    }

    curloglik = alpha[0][n_mar-1];
    for(v=1; v<n_gen; v++) 
      curloglik = addlog(curloglik, alpha[v][n_mar-1]);
    *loglik += curloglik;
  }

  if(prnt) {
    /* print final estimates */
    Rprintf(" %4d ", it+1);
    for(j=0; j<n_mar-1; j++) Rprintf("%.3lf ", rf[j]);
    Rprintf("\n");
    
    Rprintf("loglik: %10.4lf\n\n", *loglik);
  }

}



/**********************************************************************
 * 
 * argmax_geno
 *
 * This function uses the Viterbi algorithm to calculate the most 
 * likely sequence of underlying genotypes, given the observed marker
 * data for a chromosome.
 * This assumes data on a single chromosome
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              find most likely genotypes)
 *
 * n_gen        Number of different genotypes
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps (may be ignored)
 *
 * error_prob   Genotyping error probability
 *
 * argmax       Matrix of most likely genotypes (the output); a single 
 *              vector stored by columns (ind moves fastest, then pos)
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void argmax_geno(int n_ind, int n_pos, int n_gen, int *geno, 
		 double *rf, double *rf2, 
		 double error_prob, int *argmax, 
		 double initf(int), 
		 double emitf(int, int, double),
		 double stepf(int, int, double, double)) 
{
  int i, j, v, v2;
  double s, t, *gamma, *tempgamma;
  int **Geno, **Argmax, **traceback;
  
  /* Read R's random seed */
  /* in the case of multiple "most likely" genotype sequences, 
     we pick from them at random */
  GetRNGstate();

  /* allocate space and 
     reorganize geno and argmax */
  reorg_geno(n_ind, n_pos, geno, &Geno);
  reorg_geno(n_ind, n_pos, argmax, &Argmax);
  allocate_imatrix(n_pos, n_gen, &traceback);
  allocate_double(n_gen, &gamma);
  allocate_double(n_gen, &tempgamma);

  for(i=0; i<n_ind; i++) { /* i = individual */

    /* begin viterbi algorithm */
    if(n_pos > 1) { /* multiple markers */
      for(v=0; v<n_gen; v++) 
	gamma[v] = initf(v+1) + emitf(Geno[0][i], v+1, error_prob);
    
      for(j=0; j<n_pos-1; j++) {
	for(v=0; v<n_gen; v++) {
	  tempgamma[v] = s = gamma[0] + stepf(1, v+1, rf[j], rf2[j]);
	  traceback[j][v] = 0;
	  
	  for(v2=1; v2<n_gen; v2++) {
	    t = gamma[v2] + stepf(v2+1, v+1, rf[j], rf2[j]);
	    if(t > s || (fabs(t-s) < TOL && unif_rand() < 0.5)) {
	      tempgamma[v] = s = t;
	      traceback[j][v] = v2;
	    }
	  }
	  gamma[v] = tempgamma[v] + emitf(Geno[j+1][i], v+1, error_prob);
	}

      }
    
      /* finish off viterbi and then traceback to get most 
	 likely sequence of genotypes */
      Argmax[n_pos-1][i] = 0;
      s = gamma[0];
      for(v=1; v<n_gen; v++) {
	if(gamma[v] > s || (fabs(gamma[v]-s) < TOL && 
			    unif_rand() < 0.5)) {
	  s = gamma[v];
	  Argmax[n_pos-1][i] = v;
	}
      }
      for(j=n_pos-2; j >= 0; j--) 
	Argmax[j][i] = traceback[j][Argmax[j+1][i]];
    }
    else {  /* for exactly one marker */
      s = initf(1) + emitf(Geno[0][i], 1, error_prob);
      Argmax[0][i] = 0;
      for(v=1; v<n_gen; v++) {
	t = initf(v+1)+emitf(Geno[0][i], v+1, error_prob);
	if(t > s || (fabs(t-s) < TOL && unif_rand() < 0.5)) {
	  s = t;
	  Argmax[0][i] = v;
	}
      }
    }
    
    /* code genotypes as 1, 2, ... */
    for(j=0; j<n_pos; j++) 
      Argmax[j][i]++;
    
  } /* loop over individuals */
  
  
  /* write R's random seed */
  PutRNGstate();
}



/**********************************************************************
 * 
 * calc_errorlod
 *
 * Uses the results of calc_genoprob to calculate a LOD score for 
 * each genotype, indicating whether it is likely to be in error.
 *
 * n_ind, n_mar, n_gen, geno        These are all as in the above funcs
 * error_prob, genoprob 
 *
 * errlod          The output, as a single vector stored by columns, 
 *                 of size n_ind x n_mar
 * 
 * errorlod        Function taking observed genotype, genotype probs,
 *                 and error probability, and returning the error LOD
 *
 **********************************************************************/

void calc_errorlod(int n_ind, int n_mar, int n_gen, int *geno, 
		   double error_prob, double *genoprob, double *errlod, 
		   double errorlod(int, double *, double))
{
  int i, j, k, **Geno;
  double *p, ***Genoprob, **Errlod;

  /* reorganize geno, genoprob and errlod */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  reorg_genoprob(n_ind, n_mar, n_gen, genoprob, &Genoprob);
  reorg_errlod(n_ind, n_mar, errlod, &Errlod);
  allocate_double(n_gen, &p);

  for(i=0; i<n_ind; i++) {
    for(j=0; j<n_mar; j++) {
      for(k=0; k<n_gen; k++) p[k] = Genoprob[k][j][i];
      Errlod[j][i] = errorlod(Geno[j][i], p, error_prob);
    }
  }

}



/**********************************************************************
 * 
 * est_rf
 *
 * Estimate sex-averaged recombination fractions for all pairs of loci
 *
 * This is for f2 and 4way crosses; backcrosses don't need the EM 
 * algorithm, since there is no partially missing data.
 *
 * n_ind        Number of individuals
 *
 * n_mar        Number of markers
 *
 * geno         Matrix of genotype data (n_ind x n_mar), stored as a 
 *              single vector (by columns)
 * 
 * rf           The output: matrix of doubles (n_mar x n_mar), stored
 *              as a single vector (by columns).  The diagonal will 
 *              contain the number of meioses, the lower triangle will
 *              contain the est'd rec fracs, and the upper triangle
 *              will contain the LOD scores (testing rf=0.5)
 *
 * erec         Function returning the expected number of recombination
 *              events given observed marker genotypes
 * 
 * logprec      Function returning the log probability of a pair of 
 *              observed genotypes, given the recombination fraction
 *              (for calculating the LOD score)
 * 
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence of the EM 
 *
 **********************************************************************/

void est_rf(int n_ind, int n_mar, int *geno, double *rf, 
	    double erec(int, int, double), 
	    double logprec(int, int, double), 
	    int maxit, double tol)
{
  int i, j1, j2, s, **Geno, n_mei=0, flag=0;
  double **Rf, next_rf=0.0, cur_rf=0.0;

  /* reorganize geno and rf */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  reorg_errlod(n_mar, n_mar, rf, &Rf);

  for(j1=0; j1<n_mar; j1++) {

    /* count number of meioses */
    for(i=0, n_mei=0; i<n_ind; i++) 
      if(Geno[j1][i] != 0) n_mei += 2;
    Rf[j1][j1] = (double) n_mei;
    
    for(j2=j1+1; j2<n_mar; j2++) {
      
      /* count meioses */
      n_mei = flag = 0;
      for(i=0; i<n_ind; i++) {
	if(Geno[j1][i] != 0 && Geno[j2][i] != 0) {
	  n_mei += 2;
	  /* check if informatve */
	  if(logprec(Geno[j1][i], Geno[j2][i], 0.5) < 0.0) flag = 1;
	}
      }

      if(n_mei != 0 && flag == 1) {

	/* begin EM algorithm; start with cur_rf = 0.5 */
	for(s=0, cur_rf=0.5; s < maxit; s++) {
	  next_rf = 0.0; 
	  for(i=0; i<n_ind; i++) {
	    if(Geno[j1][i] != 0 && Geno[j2][i] != 0) 
	      next_rf += erec(Geno[j1][i], Geno[j2][i], cur_rf);
	  }

	  next_rf /= (double) n_mei;
	  
	  if(fabs(next_rf - cur_rf) < tol) { 
	    flag = 1;
	    break;
	  }
	  cur_rf = next_rf;
	}
	if(!flag) warning("Markers (%d,%d) didn't converge\n", j1+1, j2+1);

	/* calculate LOD score */
	Rf[j1][j2] = next_rf;
	Rf[j2][j1] = 0.0;
	for(i=0; i<n_ind; i++) {
	  if(Geno[j1][i] != 0 && Geno[j2][i] != 0) {
	    Rf[j2][j1] += logprec(Geno[j1][i],Geno[j2][i], next_rf);
	    Rf[j2][j1] -= logprec(Geno[j1][i],Geno[j2][i], 0.5);
	  }
	}
	Rf[j2][j1] /= log(10.0);
	
      }
      else { /* no informative meioses */
	Rf[j1][j2] = NA_REAL;
	Rf[j2][j1] = 0.0;
      }
	  
    } /* end loops over markers */
  }
}

/* end of hmm_main.c */
