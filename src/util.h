/**********************************************************************
 * 
 * util.h
 *
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 * Feb, 2001
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These are utility functions, mostly for the HMM engine.
 *
 * Other functions: addlog, subtrlog, reorg_geno, reorg_genoprob,
 *                  allocate_alpha, reorg_draws, allocate_double,
 *                  sample_int, allocate_imatrix, allocate_dmatrix,
 *                  reorg_errlod, double_permute, random_int
 *
 **********************************************************************/

/**********************************************************************
 * 
 * addlog
 *
 * Calculate addlog(a,b) = log[exp(a) + exp(b)]
 *
 * This makes use of the function log1p(x) = log(1+x) provided
 * in R's math library.
 *
 **********************************************************************/

double addlog(double a, double b);
		       
/**********************************************************************
 * 
 * subtrlog
 *
 * Calculate subtrlog(a,b) = log[exp(a) - exp(b)]
 *
 * This makes use of the function log1p(x) = log(1+x) provided
 * in R's math library.
 *
 **********************************************************************/

double subtrlog(double a, double b);

/**********************************************************************
 * 
 * reorg_geno
 *
 * Reorganize the genotype data so that it is a doubly indexed array
 * rather than a single long vector
 *
 * Afterwards, geno indexed like Geno[mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno);

/**********************************************************************
 * 
 * reorg_genoprob
 *
 * Reorganize the genotype probability data so that it is a triply 
 * indexed array rather than a single long vector
 *
 * Afterwards, genoprob indexed like Genoprob[gen][mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void reorg_genoprob(int n_ind, int n_pos, int n_gen, 
		    double *genoprob, double ****Genoprob);

/**********************************************************************
 * 
 * allocate_alpha
 *
 * Allocate space for alpha and beta matrices
 *
 * Afterwards, indexed like alpha[gen][mar]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void allocate_alpha(int n_pos, int n_gen, double ***alpha);

/**********************************************************************
 * 
 * reorg_draws
 *
 * Reorganize the simulated genotypes so that it is a triply 
 * indexed array rather than a single long vector
 *
 * Afterwards, draws indexed like Draws[repl][mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void reorg_draws(int n_ind, int n_pos, int n_draws, 
		 int *draws, int ****Draws);

/**********************************************************************
 * 
 * allocate_double
 *
 * Allocate space for a vector of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void allocate_double(int n, double **vector);

/**********************************************************************
 * 
 * allocate_dmatrix
 *
 * Allocate space for a matrix of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void allocate_dmatrix(int n_row, int n_col, double ***matrix);

/**********************************************************************
 * 
 * allocate_imatrix
 *
 * Allocate space for a matrix of ints
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void allocate_imatrix(int n_row, int n_col, int ***matrix);

/**********************************************************************
 * 
 * sample_int
 *
 * Make a single draw from (1, ..., n) with probs (p_0, ..., p_(n-1))
 *
 **********************************************************************/

int sample_int(int n, double *p);


/**********************************************************************
 * 
 * reorg_errlod
 *
 * Just like reorg_geno(), only for a matrix of doubles.
 *
 * Afterwards, errlod indexed like Errlod[mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void reorg_errlod(int n_ind, int n_mar, double *errlod, double ***Errlod);

/**********************************************************************
 * 
 * double_permute
 *
 *   This function randomly permutes a vector of doubles
 *   
 * Input:
 * 
 *   array = vector of doubles; on output, it contains a random 
 *           permutation of the input vector
 *
 *   len   = length of the vector
 *
 **********************************************************************/

void double_permute(double *array, long len);

/**********************************************************************
 * 
 * random_int
 *   
 * Generates a random int integer between "low" and "high", inclusive.
 *
 *  Input:
 * 
 *    low
 *
 *    high
 *
 **********************************************************************/

int random_int(int low, int high);

/* end of util.h */
