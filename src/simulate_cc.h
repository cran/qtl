/**********************************************************************
 * 
 * simulate_cc.h
 *
 * copyright (c) 2005, Karl W Broman, Johns Hopkins University
 *
 * last modified Mar, 2005
 * first written Mar, 2005
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for simulating experimental cross data;
 * I start with the simulation of RILs
 *
 * Contains: R_sim_ril, sim_all_ril, sim_ril, meiosis, sim_cc, R_sim_cc
 *  
 **********************************************************************/

/* wrapper for sim_all_ril, to be called from R */
void R_sim_ril(int *n_chr, int *n_mar, int *n_ril, double *map,
	       int *n_str, int *m, int *include_x, int *random_cross,
	       int *cross, int *ril);
	  
/**********************************************************************
 * 
 * sim_all_ril
 * 
 * n_chr   Number of chromosomes
 * n_mar   Number of markers on each chromosome (vector of length n_chr)
 * n_ril   Number of RILs to simulate
 * 
 * map     Vector of marker locations, of length sum(n_mar)
 *         First marker on each chromosome should be at 0.
 *
 * n_str   Number of parental strains (either 2, 4, or 8)
 *
 * m       Interference parameter (0 is no interference)
 *
 * include_x   Whether the last chromosome is the X chromosome
 *
 * random_cross  Indicates whether the order of the strains in the cross
 *               should be randomized.
 *
 * cross   On output, the cross used for each line 
 *         (vector of length n_ril x n_str 
 *
 * ril     On output, the simulated data 
 *         (vector of length sum(n_mar) x n_ril)
 *
 **********************************************************************/
void sim_all_ril(int n_chr, int *n_mar, int n_ril, double *map, 
		 int n_str, int m, int include_x, int random_cross,
		 int *cross, int *ril);

/**********************************************************************
 * 
 * sim_ril: simulate a single RIL to fixation
 *  
 * n_chr   Number of chromosomes
 * n_mar   Number of markers on each chr (vector of length n_chr)
 * tot_mar  sum(n_mar)
 *
 * map     Positions of markers (referred to as
 *         map[mar][chr]; map[0][i] should = 0)
 *
 * n_str   Number of strains (2, 4, or 8)
 * m       Interference parameter (0 is no interference)
 *
 * ril     On exit, vector of genotypes (length tot_mar)
 * 
 * include_x    Indicates whether the last chromosome is the X chromosome
 *
 * random_cross Indicates whether to randomize order of parents in the cross
 *
 * cross    On exit, vector indicating the cross done (length n_str)
 *
 * Par1a...Kid1b  Workspace to contain intermediate genotypes
 *                Dimension n_chr x n_mar, referred to as Par1a[chr][mar]
 *
 **********************************************************************/
void sim_ril(int n_chr, int *n_mar, int tot_mar, double **map, int n_str, 
	     int m, int *ril, int include_x, int random_cross, int *cross,
	     int **Par1a, int **Par1b, int **Par2a, int **Par2b, 
	     int **Kid1a, int **Kid1b, int **Kid2a, int **Kid2b);

/**********************************************************************
 * 
 * meiosis
 *
 * n_mar  Number of markers
 *
 * chr1   Vector of alleles along one chromosome
 * chr2   Vector of alleles along the other chromosome
 *
 * map    Marker locations in cM (need first position to be 0)
 *
 * m      interference parameter (0 corresponds to no interference)
 *
 * product  vector of length n_mar: the chromosome produced
 *
 **********************************************************************/
void meiosis(int n_mar, int *chr1, int *chr2, double *map, int m, int *product);

/**********************************************************************
 * 
 * sim_cc    Use the result of sim_all_ril with n_str=8 plus data on
 *           the SNP genotypes of the 8 parental strains to create 
 *           real SNP data for the Collaborative Cross
 *
 * n_ril     Number of RILs to simulate
 * tot_mar   Total number of markers
 *
 * Parents   SNP data for the 8 parental lines [dim tot_mar x 8]
 * 
 * Geno      On entry, the detailed genotype data; on exit, the 
 *           SNP data written bitwise.
 * 
 * error_prob  Probability of genotyping error
 * missing_prob  Probability a genotype will be missing
 *
 **********************************************************************/
void sim_cc(int n_ril, int tot_mar, int **Parents, int **Geno,
	    double error_prob, double missing_prob);

/* wrapper for calling sim_cc from R */
void R_sim_cc(int *n_ril, int *tot_mar, int *parents, int *geno,
	      double *error_prob, double *missing_prob);

/* end of simulate_cc.h */

