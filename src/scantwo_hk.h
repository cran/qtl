/**********************************************************************
 * 
 * scantwo_hk.h
 *
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 *
 * last modified Nov, 2001
 * first written Nov, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a two-dimensional genome scan with  
 * a two-QTL model by Haley-Knott regression
 *
 * Contains: R_scantwo_1chr_hk, scantwo_1chr_hk, 
 *           R_scantwo_2chr_hk, scantwo_2chr_hk
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_scantwo_1chr_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_1chr_hk.
 * 
 **********************************************************************/

void R_scantwo_1chr_hk(int *n_ind, int *n_pos, int *n_gen,
		       double *genoprob, double *pairprob, 
		       double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *result);

/**********************************************************************
 * 
 * scantwo_1chr_hk
 *
 * Performs a 2-dimensional genome scan using the Haley-Knott 
 * regression method (regressing phenotypes on conditional genotype 
 * probabilities) for a two-QTL model with the two QTL residing on
 * the same chromosome.
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *              Indexed as Genoprob[gen][pos][ind]
 *
 * Pairprob     Array of joint genotype probabilities for QTL
 *              pairs; indexed as Pairprob[gen1][gen2][pos1][pos2][ind]
 *              where pos2 > pos1 (for pos2 <= pos1, points to nothing)
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
 * Result       Result matrix of size [n_pos x n_pos]; the upper 
 *              triangle (row < col) contains the joint LODs while 
 *              the lower triangle (row > col) contains the LODs for 
 *              testing epistasis.
 *              Note: indexed as Result[col][row]
 *
 **********************************************************************/

void scantwo_1chr_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		     double *****Pairprob, double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     double **Result);

/**********************************************************************
 * 
 * R_scantwo_2chr_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_2chr_hk.
 * 
 **********************************************************************/

void R_scantwo_2chr_hk(int *n_ind, int *n_pos1, int *n_pos2, 
		       int *n_gen1, int *n_gen2,
		       double *genoprob1, double *genoprob2,
		       double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *result_full,
		       double *result_int);

/**********************************************************************
 * 
 * scantwo_2chr_hk
 *
 * Performs a 2-dimensional genome scan using the Haley-Knott 
 * regression method (regressing phenotypes on conditional genotype 
 * probabilities) for a two-QTL model with the two QTL residing on
 * the different chromosomes.
 * 
 * n_ind        Number of individuals
 *
 * n_pos1       Number of marker positions on first chromosome
 *
 * n_pos2       Number of marker positions on second chromosome
 *
 * n_gen1       Number of different genotypes on first chromosome
 *
 * n_gen2       Number of different genotypes on second chromosome
 *
 * Genoprob1    Array of conditional genotype probs for 1st chr
 *              Indexed as Genoprob[gen][pos][ind]
 *
 * Genoprob2    Array of conditional genotype probs for 2nd chr
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
 * Result_full  Result matrix of size [n_pos1 x n_pos2]
 *              containing the joint LODs
 *              Note: indexed as Result[pos2][pos1]
 *
 * Result_int   Result matrix of size [n_pos1 x n_pos2] 
 *              containing the LODs testing interactions
 *
 **********************************************************************/

void scantwo_2chr_hk(int n_ind, int n_pos1, int n_pos2, int n_gen1, 
		     int n_gen2,
		     double ***Genoprob1, double ***Genoprob2, 
		     double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     double **Result_full, double **Result_int);

/* end of scantwo_hk.h */
