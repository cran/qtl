/**********************************************************************
 * 
 * hmm_bci.c
 * 
 * copyright (c) 2006, Karl W Broman, Johns Hopkins University
 *         (Some code adapted from code from Nicola Armstrong)
 *
 * last modified Nov, 2006
 * first written Aug, 2006
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * Contains: est_map_bci, emit_bci, nrec_bci, R_est_map_bci,
 *           mf_stahl, R_mf_stahl, imf_stahl, R_imf_stahl,
 *           imf_stahl_sub
 *
 * These are functions for the HMM under the Stahl model
 * (with chiasmata coming from two mechanisms: one following a 
 * chi-square model and one following a no interference model).
 * m = interference parameter in the chi-square model (m=0 == NI)
 * p = proportion of chiasmata from the NI model (p=1 == NI)
 *
 * This is exclusively for a backcross.
 *
 * Genotype codes:  0, ..., 2(m+1) - 1, with the first (m+1) 
 *                  corresponding to AA and the others to AB
 * Phenotype codes: 0=missing; 1=AA; 2=AB
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
#include "hmm_bci.h"

/**********************************************************************
 * 
 * est_map_bci
 *
 * This function re-estimates the genetic map for a chromosome
 * with the Stahl model, taking m and p known
 *
 * n_ind        Number of individuals
 *
 * n_mar        Number of markers 
 *
 * n_gen        Number of different genotypes (strictly = 2), 
 *              but included for later expansion to intercross
 *
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * d            inter-marker distances in cM
 *              (on exit, contains the new estimates)
 *
 * m            Interference parameter (non-negative integer)
 *
 * p            Proportion of chiasmata from the NI mechanism
 *
 * error_prob   Genotyping error probability
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 * nrecf        Function returning number of recombinations associated
 *              with (g_1, g_2)
 *
 * loglik       Loglik at final estimates of recombination fractions
 *
 * maxit        Maximum number of iterations to perform
 * 
 * tol          Tolerance for determining convergence
 * 
 **********************************************************************/

/* Note: true genotypes coded as 0, 1, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void est_map_bci(int n_ind, int n_mar, int n_gen, int *geno, double *d, 
		 int m, double p, double error_prob, 
		 double emitf(int, int, double, int),
		 void stepf(int, int, double ***, double *, int, double, 
			    int, double), 
		 double nrecf(int, int, int), 
		 double *loglik, int maxit, double tol, int verbose)
{
  int i, j, j2, v, v2, it, flag=0, **Geno, n_states;
  double s, **alpha, **beta, **gamma, *cur_d, *rf;
  double ***tm, *temp;
  double curloglik;
  double initprob;
  
  n_states = 2*(m+1); /* number of states in the HMM */
  initprob = -log((double)n_states);

  /* allocate space for beta and reorganize geno */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  allocate_alpha(n_mar, n_states, &alpha);
  allocate_alpha(n_mar, n_states, &beta);
  allocate_dmatrix(n_states, n_states, &gamma);
  allocate_double(n_mar-1, &cur_d);
  allocate_double(n_mar-1, &rf);

  /* allocate space for the transition matrices */
  /* size n_states x n_states x (n_mar-1) */
  /* tm[state1][state2][interval] */
  tm = (double ***)R_alloc(n_states, sizeof(double **));
  tm[0] = (double **)R_alloc(n_states * n_states, sizeof(double *));
  for(i=1; i<n_states; i++) tm[i] = tm[i-1] + n_states;
  tm[0][0] = (double *)R_alloc(n_states * n_states * (n_mar - 1), 
			       sizeof(double));
  temp = tm[0][0];
  for(i=0; i < n_states; i++) {
    for(j=0; j < n_states; j++) {
      tm[i][j] = temp;
      temp += n_mar-1;
    }
  }

  if(verbose) {
    /* print initial estimates */
    Rprintf("      "); 
    for(j=0; j<n_mar-1; j++) Rprintf("%.3lf ", d[j]);
    Rprintf("\n"); 
  }

  for(j=0; j<n_mar-1; j++) d[j] /= 100.0; /* convert to Morgans */

  /* begin EM algorithm */
  for(it=0; it<maxit; it++) {

    for(j=0; j<n_mar-1; j++) {
      cur_d[j] = d[j];
      rf[j] = 0.0;
    }

    /* calculate the transition matrices */
    step_bci(n_mar, n_states, tm, cur_d, m, p, maxit, tol);

    for(i=0; i<n_ind; i++) { /* i = individual */

      /* initialize alpha and beta */
      for(v=0; v<n_states; v++) {
	alpha[v][0] = initprob + emitf(Geno[0][i], v, error_prob, m);
	beta[v][n_mar-1] = 0.0;
      }

      /* forward-backward equations */
      for(j=1,j2=n_mar-2; j<n_mar; j++, j2--) {
	
	for(v=0; v<n_states; v++) {
	  alpha[v][j] = alpha[0][j-1] + tm[0][v][j-1];
	  
	  beta[v][j2] = beta[0][j2+1] + tm[v][0][j2] +
	    emitf(Geno[j2+1][i], 0, error_prob, m);
	  
	  for(v2=1; v2<n_states; v2++) {
	    alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
				 tm[v2][v][j-1]);
	    beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
				 tm[v][v2][j2] +
				 emitf(Geno[j2+1][i], v2, error_prob, m));
	  }
	  
	  alpha[v][j] += emitf(Geno[j][i], v, error_prob, m);
		 
	}

      }

      for(j=0; j<n_mar-1; j++) {

	/* calculate gamma = log Pr(v1, v2, O) */
	for(v=0, s=0.0; v<n_states; v++) {
	  for(v2=0; v2<n_states; v2++) {
	    gamma[v][v2] = alpha[v][j] + beta[v2][j+1] + 
	      emitf(Geno[j+1][i], v2, error_prob, m) +
	      tm[v][v2][j];

	    if(v==0 && v2==0) s = gamma[v][v2];
	    else s = addlog(s, gamma[v][v2]);
	  }
	}

	for(v=0; v<n_states; v++) {
	  for(v2=0; v2<n_states; v2++) {
	    rf[j] += nrecf(v, v2, m) * exp(gamma[v][v2] - s);
	  }
	}
      }

    } /* loop over individuals */

    /* rescale */
    for(j=0; j<n_mar-1; j++) {
      rf[j] /= (double)n_ind;
      if(rf[j] < tol/100.0) rf[j] = tol/100.0;
      else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;
    }

    /* use map function to convert back to distances */
    for(j=0; j<n_mar-1; j++)
      d[j] = imf_stahl(rf[j], m, p, 1e-10, 1000);

    if(verbose > 1) { /* print some debugging stuff */
      if(verbose == 2) Rprintf("Iteration");
      Rprintf(" %4d ", it+1);
      if(verbose > 2) 
	for(j=0; j<n_mar-1; j++) Rprintf("%.3lf ", d[j]*100);
      Rprintf("\n"); 
    }

    /* check convergence */
    for(j=0, flag=0; j<n_mar-1; j++) {
      if(fabs(d[j] - cur_d[j]) > tol*(cur_d[j]+tol*100.0)) {
	flag = 1; 
	break;
      }
    }

    if(!flag) break;

  } /* end EM algorithm */
  
  if(flag) warning("Didn't converge!\n");

  /* re-calculate transition matrices */
  step_bci(n_mar, n_states, tm, d, m, p, maxit, tol);

  /* calculate log likelihood */
  *loglik = 0.0;
  for(i=0; i<n_ind; i++) { /* i = individual */
    /* initialize alpha */
    for(v=0; v<n_states; v++) 
      alpha[v][0] = initprob + emitf(Geno[0][i], v, error_prob, m);

    /* forward equations */
    for(j=1; j<n_mar; j++) {
      for(v=0; v<n_states; v++) {
	alpha[v][j] = alpha[0][j-1] + tm[0][v][j-1];
	for(v2=1; v2<n_states; v2++) 
	  alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
			       tm[v2][v][j-1]);
	alpha[v][j] += emitf(Geno[j][i], v, error_prob, m);
      }
    }

    curloglik = alpha[0][n_mar-1];
    for(v=1; v<n_states; v++) 
      curloglik = addlog(curloglik, alpha[v][n_mar-1]);
    *loglik += curloglik;
  }

  /* convert distances back to cM */
  for(j=0; j<n_mar-1; j++) d[j] *= 100.0;

  if(verbose) {
    /* print final estimates */
    Rprintf(" %4d ", it+1);
    for(j=0; j<n_mar-1; j++) Rprintf("%.3lf ", d[j]);
    Rprintf("\n");
    
    Rprintf("loglik: %10.4lf\n\n", *loglik);
  }

}



double emit_bci(int obs_gen, int true_gen, double error_prob,
		int m)
{
  if(true_gen < m+1) true_gen = 1;
  else true_gen = 2;

  switch(obs_gen) {
  case 0: return(0.0);
  case 1: case 2:
    if(obs_gen == true_gen) return(log(1.0-error_prob));
    else return(log(error_prob));
  }
  return(0.0); /* shouldn't get here */
}

double nrec_bci(int gen1, int gen2, int m)  
{
  if(gen1 < m+1) gen1=0;
  else gen1=1;
  if(gen2 < m+1) gen2=0;
  else gen2=1;

  if(gen1==gen2) return(0.0);
  else return(1.0);
}

/* R wrapper for est_map_bci */
void R_est_map_bci(int *n_ind, int *n_mar, int *geno, double *d, 
		   int *m, double *p, double *error_prob, 
		   double *loglik, int *maxit, double *tol, int *verbose)
{

  est_map_bci(*n_ind, *n_mar, 2, geno, d, *m, *p,
	      *error_prob, emit_bci, step_bci, nrec_bci,
	      loglik, *maxit, *tol, *verbose);
}


/***********************************************************************
 * R_mf_stahl: wrapper for R
 * 
 * n_d = length of vector d
 * d   = genetic distances (in Morgans)
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 * result = vector of length n_d to contain the results
 **********************************************************************/
void R_mf_stahl(int *n_d, double *d, int *m, double *p, double *result)
{
  int i;
  
  for(i=0; i<*n_d; i++)
    result[i] = mf_stahl(d[i], *m, *p);
}
  
/**********************************************************************
 * mf_stahl: map function for Stahl model
 * 
 * d   = genetic distances (in cM)
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 **********************************************************************/
double mf_stahl(double d, int m, double p)
{
  int i;
  double result, lam1, lam2;
  
  lam1 = (double)(m+1) * d *(1.0-p) * 2.0;
  lam2 = d * p * 2.0;

  result = 0.0;
  for(i=0; i<m+1; i++)
    result += (1.0 - (double)i/(double)(m+1))*dpois((double)i, lam1, 0);

  return(0.5*(1.0 - result * exp(-lam2)));
}

/**********************************************************************
 * R_imf_stahl: wrapper for R
 * 
 * n_r = length of vector r
 * r   = recombination fractions
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 * result = vector of length n_r to contain the results
 * tol = tolerance for convergence
 * maxit = number of interations
 **********************************************************************/
void R_imf_stahl(int *n_r, double *r, int *m, double *p,
		 double *result, double *tol, int *maxit)
{
  int i;
  
  for(i=0; i<*n_r; i++)
    result[i] = imf_stahl(r[i], *m, *p, *tol, *maxit);
}


/**********************************************************************
 * imf_stahl: inverse map function for chi-square model
 * 
 * r   = recombination fraction
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 * tol = tolerance for convergence
 * maxit = number of interations
 **********************************************************************/
double imf_stahl(double r, int m, double p, double tol, int maxit)
{
  double result;
  struct imf_stahl_data info;

  info.r = r;
  info.m = m;
  info.p = p;
  
  result = R_zeroin(r, -log(1.0-2.0*r)/2.0,  /* lower and upper of range */
		    imf_stahl_sub, (void *)(&info), &tol, &maxit);
  return(result);
}
  

/**********************************************************************
 * imf_stahl_sub: utility function for imf_stahl
 * 
 * r   = recombination fraction
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 * tol = tolerance for convergence
 * maxit = number of interations
 **********************************************************************/
double imf_stahl_sub(double d, void *info)
{
  int m;
  double r, p;
  double result;
  struct imf_stahl_data temp;

  temp = *((struct imf_stahl_data *)info);

  result = mf_stahl(d, temp.m, temp.p);
  
  return(result - temp.r);
}

/**********************************************************************
 * step_bci
 * 
 * Calculate transition probabilities (for all intervals) for
 * the Stahl model
 **********************************************************************/
void step_bci(int n_mar, int n_states, double ***tm, double *d, 
	      int m, double p, int maxit, double tol)
{
  int i, v1, v2;
  double *the_distinct_tm;
  double *fms_bci_result;
  double lambda1, lambda2, rfp;

  allocate_double(2*m+1, &fms_bci_result);
  allocate_double(3*m+2, &the_distinct_tm);

  for(i=0; i<n_mar-1; i++) {
    lambda1 = d[i]*(1-p)*(double)(m+1)*2.0;
    lambda2 = d[i]*p*2.0;
    rfp = 0.5*(1.0 - exp(-lambda2));

    fms_bci(lambda1, fms_bci_result, m, tol, maxit);
    distinct_tm_bci(lambda1, the_distinct_tm, m, fms_bci_result);

    for(v1=0; v1<n_states; v1++) {
      for(v2=0; v2<n_states; v2++) {
	tm[v1][v2][i] = alltm_bci(v1, v2, the_distinct_tm, m);
	if(p > 0) 
	  tm[v1][v2][i] = (1.0-rfp)*tm[v1][v2][i] + 
	    rfp*alltm_bci(v1, (v2+m+1) % (2*m+2), the_distinct_tm, m);
	tm[v1][v2][i] = log(tm[v1][v2][i]);
      }
    }
  }
} 



/*****************************************************************************
 * alltm_bci: this function calculates the required transition probability for the
 * backcross case
 ****************************************************************************/

double alltm_bci(int i, int j, double *the_distinct_tm, int m)
{
  int s, tempi, tempj;
  
  if ((i<=m && j<=m) || (i>m && j>m)) {
    s=j-i;
    if (s>=0) {
      return(the_distinct_tm[s]);
    }
    else {
      return(the_distinct_tm[abs(s)+2*m+1]);
    }
  }
  else if (i<=m && j>m) {
    if (j>(i+m)) {
      s=j-i;
      return(the_distinct_tm[s]);
    }
    else /* j <=i+m */ {
      s=j-i-(m+1);
      return(the_distinct_tm[abs(s)+2*m+1]);
    }
  }
  else /* i>m && j<=m */ {
    tempi=i-(m+1);
    tempj=j+(m+1);
    if (tempj>(tempi+m)) {
      s=tempj-tempi;
      return(the_distinct_tm[s]);
    }
    else /* tempj <=tempi+m */ {
      s=tempj-tempi-(m+1);
      return(the_distinct_tm[abs(s)+2*m+1]);
    }
  }
}

/*****************************************************************************
 * fms_bci: this function calculates the sum to infinity part of the
 * transition probabilities for a given lambda_t
 *
 * f should have length 2m+1
 ****************************************************************************/
void fms_bci(double lambda, double *f, int m, double tol, int maxit)
{
  int i,k;
  double diff;

  for (i=0; i<2*m+1; i++) {
    k=1;
    f[i]=0;
    if (i <= m) {
      f[i] = dpois((double)(k*(m+1)+i), lambda, 0);

      for(k=2; k<maxit; k++) {
        diff = dpois((double)(k*(m+1)+i), lambda, 0);
	f[i] += diff;

	if(diff < tol) break;
      }
    }
    else /* i<m */ {
      f[i] += dpois((double)(k*(m+1)+(m-i)), lambda, 0);

      for(k=2; k<maxit; k++) {
	diff = dpois((double)(k*(m+1)+(m-i)), lambda, 0);
	f[i] += diff;

	if(diff < tol) break;
      }  
    }  
    f[i] *= 0.5;
  }  
}


/*****************************************************************************
 * distinct_tm_bci: this function calculates the 3m+2 distinct transition 
 * probabilities for a given lambda_t
 ****************************************************************************/

void distinct_tm_bci(double lambda, double *the_distinct_tm, int m, 
		     double *fms_bci_result)
{
  int i;
  
  for (i=0;i<3*m+2;i++) {
    if (i<=m) 
      the_distinct_tm[i] = fms_bci_result[i]+ dpois((double)i, lambda, 0);

    else 
      the_distinct_tm[i] = fms_bci_result[i-(m+1)];

  }

}

/* end of hmm_bci.c */
