/**********************************************************************
 * 
 * scanone_imp.h
 *
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 *                 and Hao Wu, The Jackson Laboratory
 *
 *
 * This file is written by Hao Wu (hao@jax.org), 
 * with slight modifications by Karl Broman.
 *
 * last modified Nov, 2001
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
 * addcov       Additive covariates matrix, addcov[mar][ind]
 *  
 * n_addcov     Number of additive covariates
 *
 * intcov       Interacting covariates matrix, intcov[mar][ind]
 *
 * n_intcov     Number of interacting covariates
 *
 * pheno        Phenotype data, as a vector
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
		 int ***Draws, double **addcov, int n_addcov, 
		 double **intcov, int n_intcov, double *pheno, 
		 double *result, int trim, int direct);

double nullRss(double *pheno, int n_ind, double **addcov, int n_addcov,
	       double *dwork, int *iwork);

double altRss(double *pheno, int n_ind, int n_gen, int *Draws,
	      double **addcov, int n_addcov, double **intcov, int n_intcov,
	      double *dwork, int *iwork);

/* end of scanone_imp.h */
