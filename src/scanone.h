/**********************************************************************
 * 
 * scanone.h
 *
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 * Oct, 2001; Aug, 2001; May, 2001
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a genome scan with a 
 * single QTL model
 *
 * Contains: R_scanone_anova, scanone_anova, scanone_anova_perm
 *           R_scanone_im, scanone_im, scanone_im_perm
 *           R_scanone_hk, scanone_hk, scanone_hk_perm
 *           R_scanone_imp, scanone_imp, scanone_imp_perm
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_scanone_anova
 *
 * Wrapper for call from R; reorganizes genotype and result matrix
 * and calls scanone_anova.
 * 
 **********************************************************************/

void R_scanone_anova(int *n_ind, int *n_pos, int *n_gen, int *geno,
		     double *pheno, double *result);

/**********************************************************************
 * 
 * R_scanone_im
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_im.
 * 
 **********************************************************************/

void R_scanone_im(int *n_ind, int *n_pos, int *n_gen, 
		  double *genoprob, double *pheno,
		  double *result, int *std_start, double *start,
		  int *maxit, double *tol);

/**********************************************************************
 * 
 * R_scanone_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_hk.
 * 
 **********************************************************************/

void R_scanone_hk(int *n_ind, int *n_pos, int *n_gen,
		  double *genoprob, double *pheno,
		  double *result);

/**********************************************************************
 * 
 * R_scanone_imp
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_imp.
 * 
 **********************************************************************/

void R_scanone_imp(int *n_ind, int *n_pos, int *n_gen, int *n_draws,
		   int *draws, double *pheno, double *result);

/**********************************************************************
 * 
 * scanone_anova
 *
 * Performs genotype scan using ANOVA, assuming complete genotypes 
 * (reconstructed by the Viterbi algorithm in argmax.geno)
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
 * Result       Result matrix of size [n_pos x (n_gen+2)]; upon return, 
 *              the first column contains the RSS, the next set contain
 *              estimated genotype-specific means, and the last column
 *              contains the estimated residual SD
 *
 **********************************************************************/

void scanone_anova(int n_ind, int n_pos, int n_gen, int **Geno, 
		   double *pheno, double **Result);

/**********************************************************************
 * 
 * scanone_im
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
 * work1        Workspace to store conditional probabilities
 * 
 * work2        Workspace to store genotype-specific RSS 
 *
 **********************************************************************/

void scanone_im(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		double *pheno, double **Result, int std_start, 
		double *start, int maxit, double tol, double **work);

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
 * pheno        Phenotype data, as a vector
 *
 * Result       Result matrix of size [n_pos x (n_gen+2)]; upon return, 
 *              the first column contains the RSS, the next set contain
 *              estimated genotype-specific means, and the last column
 *              contains the estimated residual SD
 *
 **********************************************************************/

void scanone_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		double *pheno, double **Result);


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
 * pheno        Phenotype data, as a vector
 *
 * result       Result vector of length [n_pos]; upon return, contains
 *              the "LPD" (log posterior distribution of QTL location).
 *
 **********************************************************************/

void scanone_imp(int n_ind, int n_pos, int n_gen, int n_draws,
		 int ***Draws, double *pheno, double *result);


/**********************************************************************
 * 
 * scanone_anova_perm
 * 
 * Do permutation test with anova
 * 
 * arguments are like that for R_scanone_anova, except result, 
 * which is a vector of length n_perm (to contain the output)
 *
 **********************************************************************/

void scanone_anova_perm(int *n_ind, int *n_pos, int *n_gen, int *geno,
			double *pheno, int *n_perm, double *result);

void scanone_im_perm(int *n_ind, int *n_pos, int *n_gen, 
		     double *genoprob, double *pheno, 
		     int *n_perm, double *result,
		     int *std_start, double *start,
		     int *maxit, double *tol);

void scanone_hk_perm(int *n_ind, int *n_pos, int *n_gen, 
		     double *genoprob, double *pheno, 
		     int *n_perm, double *result);

/* end of scanone.h */
