/**********************************************************************
 * 
 * discan.h
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

/**********************************************************************
 * 
 * R_discan_mr
 *
 * Wrapper for call from R; reorganizes genotype and result matrix
 * and calls discan_mr.
 * 
 **********************************************************************/

void R_discan_mr(int *n_ind, int *n_pos, int *n_gen,
		    int *geno, double *pheno, double *result);

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
		 int *maxit, double *tol);

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
 *              the first column contains the RSS, the next set contain
 *              estimated genotype-specific probabilities
 *
 **********************************************************************/

void discan_mr(int n_ind, int n_pos, int n_gen, int **Geno, 
		  double *pheno, double **Result);

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
 * Result       Result matrix of size [n_pos x (n_gen+2)]; upon return, 
 *              the first column contains the log10 likelihood, the 
 *              next set contain estimated genotype-specific means, and 
 *              the last column contains the estimated residual SD
 *
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * work         Workspace of length n_gen
 *
 **********************************************************************/

void discan_im(int n_ind, int n_pos, int n_gen, double ***Genoprob,
	       double *pheno, double **Result, 
	       int maxit, double tol, double **work);

/* end of discan.h */

