/**********************************************************************
 * 
 * hmm_main.h
 *
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 * Aug, 2001; Feb, 2001
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for the main HMM engine
 *
 * Contains: calc_genoprob, sim_geno, est_map, argmax_geno,
 *           calc_errorlod, est_rf
 *  
 **********************************************************************/

#define TOL    1.0e-12

/**********************************************************************
 * 
 * calc_genoprob
 *
 * This function uses the Lander-Green algorithm to calculate the 
 * genotype probabilities at each of marker and (optionally) at points
 * in-between markers, conditional on all marker data for a chromosome.
 * This assumes data on a single chromosome
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
		   double stepf(int, int, double, double));


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
	      double stepf(int, int, double, double));


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
 * loglik       Value of loglik at final estimates of rec fracs.
 *
 * maxit        Maximum number of iterations to perform
 * 
 * tol          Tolerance for determining convergence
 * 
 * sexsp        Indicates whether sex-specific maps should be estimated
 *
 * prnt         Indicates whether to print initial and final rec fracs
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
	     int prnt);


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
 *              find most likely genotypes
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
		   double stepf(int, int, double, double));

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
		   double errorlod(int, double *, double));


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
 *              contain the number of meioses, the upper triangle will
 *              contain the est'd rec fracs, and the lower triangle
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
	    int maxit, double tol);

/* end of hmm_main.h */
