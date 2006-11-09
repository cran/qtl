/**********************************************************************
 * 
 * hmm_bci.h
 * 
 * copyright (c) 2006, Karl W Broman, Johns Hopkins University
 *         (Some code adapted from code from Nicola Armstrong)
 *
 * last modified Aug, 2006
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

/* structure used by imf_stahl and imf_stahl_sub */
struct imf_stahl_data {
  double r;
  int m;
  double p;
};

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
		 double *loglik, int maxit, double tol, int verbose);

double emit_bci(int obs_gen, int true_gen, double error_prob,
		int m);

double nrec_bci(int gen1, int gen2, int m);

/* R wrapper for est_map_bci */
void R_est_map_bci(int *n_ind, int *n_mar, int *geno, double *d, 
		   int *m, double *p, double *error_prob, 
		   double *loglik, int *maxit, double *tol, int *verbose);


/***********************************************************************
 * R_mf_stahl: wrapper for R
 * 
 * n_d = length of vector d
 * d   = genetic distances (in Morgans)
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 * result = vector of length n_d to contain the results
 **********************************************************************/
void R_mf_stahl(int *n_d, double *d, int *m, double *p, double *result);
  
/**********************************************************************
 * mf_stahl: map function for Stahl model
 * 
 * d   = genetic distances (in cM)
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 **********************************************************************/
double mf_stahl(double d, int m, double p);

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
		 double *result, double *tol, int *maxit);

/**********************************************************************
 * imf_stahl: inverse map function for chi-square model
 * 
 * r   = recombination fraction
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 * tol = tolerance for convergence
 * maxit = number of interations
 **********************************************************************/
double imf_stahl(double r, int m, double p, double tol, int maxit);

/**********************************************************************
 * imf_stahl_sub: utility function for imf_stahl
 * 
 * r   = recombination fraction
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 * tol = tolerance for convergence
 * maxit = number of interations
 **********************************************************************/
double imf_stahl_sub(double d, void *info);

/**********************************************************************
 * step_bci
 * 
 * Calculate transition probabilities (for all intervals) for
 * the Stahl model
 **********************************************************************/
void step_bci(int n_mar, int n_states, double ***tm, double *d, 
	      int m, double p, int maxit, double tol);

/*****************************************************************************
 * alltm_bci: this function calculates the required transition probability for the
 * backcross case
 ****************************************************************************/

double alltm_bci(int i, int j, double *the_distinct_tm, int m);

/*****************************************************************************
 * fms_bci: this function calculates the sum to infinity part of the
 * transition probabilities for a given lambda_t
 *
 * f should have length 2m+1
 ****************************************************************************/
void fms_bci(double lambda, double *f, int m, double tol, int maxit);

/*****************************************************************************
 * distinct_tm_bci: this function calculates the 3m+2 distinct transition 
 * probabilities for a given lambda_t
 ****************************************************************************/
void distinct_tm_bci(double lambda, double *the_distinct_tm, int m, 
		     double *fms_bci_result);

/* end of hmm_bci.h */
