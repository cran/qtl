/**********************************************************************
 *
 * scantwo_imp.h
 *
 * copyright (c) 2001-2, Karl W Broman, Johns Hopkins University
 *                     and Hao Wu, The Jackson Lab
 *
 * This file was written by Hao Wu with modifications by 
 * Karl Broman.
 *
 * last modified Oct, 2002 
 * first written Nov, 2001 
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a 2-dimensional genome scan 
 * with a 2-QTL model by imputation.
 *
 * Contains: R_scantwo_imp, scantwo_imp, altRss2
 *
 **********************************************************************/

/**********************************************************************
 *
 * R_scantwo_imp
 *
 * Wrapper for call from R; reorganizes genotype prob, additive and 
 * interactive covariates and result matrix. Then calls scantwo_imp.
 *
 **********************************************************************/

void R_scantwo_imp(int *n_ind, int *same_chr, int *n_pos1, int *n_pos2, 
		   int *n_gen1, int *n_gen2, int *n_draws, int *draws1, 
		   int *draws2, double *addcov, int *n_addcov, 
		   double *intcov, int *n_intcov, double *pheno, 
		   double *weights, double *result);

/**********************************************************************
 * 
 * scantwo_imp
 *
 * Performs genotype pair scan using the pseudomarker algorithm 
 * (imputation) method of Sen and Churchill (2001).
 * 
 * n_ind        Number of individuals
 *
 * same_chr     If = 1, work only with Draws1 and do 2-QTL model with
 *              QTLs on the same chromosome.
 *
 * chr2         Chromesome id 2
 *
 * n_pos1       Number of marker positions in chromesome 1
 *
 * n_pos2       Number of marker positions in chromesome 2
 *
 * n_gen1       Number of different genotypes on chr 1
 *
 * n_gen2       Number of different genotypes on chr 2
 *
 * n_draws      Number of impiutations
 *
 * Draws1       Array of genotype imputations in chromesome 1, 
 *              indexed as Draws1[repl][mar][ind]
 * 
 * Draws2       Array of genotype imputations in chromesome 2, 
 *              indexed as Draws2[repl][mar][ind]
 *
 * addcov	Additive covariates matrix, addcov[mar][ind]
 *
 * n_addcov     Number of additive covariates
 *
 * intcov	Interacting covariates matrix, intcov[mar][ind]
 *
 * n_intcov     Number of interacting covariates
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * result       Result vector of length [n_pos1*n_pos2];
 *
 **********************************************************************/

void scantwo_imp(int n_ind, int same_chr, int n_pos1, int n_pos2, 
		 int n_gen1, int n_gen2, int n_draws, int ***Draws1, 
		 int ***Draws2, double **Addcov, int n_addcov, 
		 double **Intcov, int n_intcov, double *pheno, 
		 double *weights, double *result);

void altRss2(double *pheno, double *weights, 
	     int n_ind, int n_gen1, int n_gen2, 
	     int *Draws1, int *Draws2, double **Addcov, int n_addcov, 
	     double **Intcov, int n_intcov, double *lrss_add, 
	     double *lrss_full, double *dwork, int *iwork);

/* end of scantwo_imp.h */
