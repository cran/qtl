/**********************************************************************
 * 
 * scanone_imp.h
 *
 * copyright (c) 2001-2, Karl W Broman, Johns Hopkins University
 *                 and Hao Wu, The Jackson Laboratory
 *
 * This file is written by Hao Wu (hao@jax.org), 
 * with slight modifications by Karl Broman.
 *
 * last modified Oct, 2002
 * first written Nov, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a genome scan with a 
 * single QTL model by imputation.  
 *
 * Contains: R_scanone_imp, scanone_imp, nullRss, altRss
 *
 **********************************************************************/

/**********************************************************************
 * 
 * R_scanone_imp
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_imp.
 * 
 **********************************************************************/

void R_scanone_imp(int *n_ind, int *n_pos, int *n_gen, int *n_draws, 
		   int *draws, double *addcov, int *n_addcov, 
		   double *intcov, int *n_intcov, double *pheno, 
		   double *weights,
		   double *result, int *trim, int *direct);

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
 * Addcov	Additive covariates matrix, Addcov[mar][ind]
 *
 * n_addcov     Number of additive covariates
 *
 * Intcov	Interacting covariates matrix, Intcov[mar][ind]
 *
 * n_intcov     Number of interacting covariates
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * result       Result vector of length [n_pos]; upon return, contains
 *              the "LPD" (log posterior distribution of QTL location).
 * 
 * trim         If 1, trim off the top and bottom log2(n.draws) LODs 
 *
 * direct       If 1, return log10[mean(10^LOD)]; if 0, return
 *              mean(LOD) + var(LOD)/2
 *
 **********************************************************************/

void scanone_imp(int n_ind, int n_pos, int n_gen, int n_draws, 
		 int ***Draws, double **Addcov, int n_addcov, 
		 double **Intcov, int n_intcov, double *pheno, 
		 double *weights,
		 double *result, int trim, int direct);

/* function to calculate the null model RSS for scanone_imp */
double nullRss(double *pheno, double *weights, int n_ind, 
	       double **Addcov, int n_addcov, 
	       double *dwork, int *iwork);

/* function to calculate the alternative model RSS. 
   This function is called by scanone_imp */
double altRss(double *pheno, double *weights, int n_ind, int n_gen, 
	      int *Draws, double **Addcov, int n_addcov, double **Intcov, 
	      int n_intcov, double *dwork, int *iwork);

/* end of scanone_imp.h */
