/**********************************************************************
 * 
 * scanone_hk.h
 *
 * copyright (c) 2001-2, Karl W Broman, Johns Hopkins University
 *
 * last modified Oct, 2002
 * first written Nov, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a genome scan with a 
 * single QTL model by Haley-Knott regression
 *
 * Contains: R_scanone_hk, scanone_hk
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_scanone_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_hk.
 * 
 **********************************************************************/

void R_scanone_hk(int *n_ind, int *n_pos, int *n_gen,
		  double *genoprob, double *addcov, int *n_addcov, 
                  double *intcov, int *n_intcov, double *pheno,
		  double *weights, double *result);

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
 * Addcov       Matrix of additive covariates: Addcov[cov][ind]
 * 
 * n_addcov     Number of columns of Addcov
 *
 * Intcov       Number of interactive covariates: Intcov[cov][ind]
 *
 * n_intcov     Number of columns of Intcov
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * Result       Result matrix of size [n_pos x (n_gen+2)]; upon return, 
 *              the first column contains the LOD, the next set contain
 *              estimated genotype-specific means, and the last column
 *              contains the estimated residual SD
 *
 **********************************************************************/

void scanone_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
                double **Addcov, int n_addcov, double **Intcov, 
		int n_intcov, double *pheno, double *weights,
		double **Result);

/* end of scanone_hk.h */
