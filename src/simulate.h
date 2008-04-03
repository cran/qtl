/**********************************************************************
 * 
 * simulate.h
 *
 * copyright (c) 2006, Karl W Broman
 *
 * last modified Jul, 2006
 * first written Jul, 2006
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for simulating backcross genotype data
 *
 * Contains: sim_bc_ni, sim_bc, R_sim_bc, R_sim_bc_ni
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_sim_bc_ni   Wrapper for sim_bc_ni
 *
 * geno is empty, of size n_mar * n_ind
 *
 **********************************************************************/

void R_sim_bc_ni(int *n_mar, int *n_ind, double *rf, int *geno);

/**********************************************************************
 * 
 * sim_bc_ni    Simulate backcross under no interference
 *
 * n_mar    Number of markers
 * n_ind    Number of individuals
 * rf       recombination fractions (length n_mar-1)
 * Geno     Matrix of size n_ind x n_mar to contain genotype data
 *
 **********************************************************************/

void sim_bc_ni(int n_mar, int n_ind, double *rf, int **Geno);

/**********************************************************************
 * 
 * R_sim_bc   Wrapper for sim_bc
 *
 * geno is empty, of size n_mar * n_ind
 *
 **********************************************************************/

void R_sim_bc(int *n_mar, int *n_ind, double *pos,
	      int *m, double *p, int *geno);

/**********************************************************************
 * 
 * sim_bc    Simulate backcross under Stahl's interference model
 *
 * n_mar    Number of markers
 * n_ind    Number of individuals
 * pos      Positions of markers (in cM)
 * m        Interference parameter (integer > 0)
 * p        Probability chiasma comes from no interference mechanism
 * Geno     Matrix of size n_ind x n_mar to contain genotype data
 *
 **********************************************************************/

void sim_bc(int n_mar, int n_ind, double *pos, int m, double p, int **Geno);

/* end of simulate.h */
