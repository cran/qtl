/**********************************************************************
 * 
 * hmm_bc.h
 * 
 * copyright (c) 2001-4, Karl W Broman, Johns Hopkins University
 *
 * last modified Nov, 2004
 * first written Feb, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * Contains: init_bc, emit_bc, step_bc, nrec_bc, calc_genoprob_bc,
 *           sim_geno_bc, est_map_bc, argmax_geno_bc, errorlod_bc,
 *           calc_errorlod_bc, est_rf_bc, calc_pairprob_bc
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the backcross.
 *
 * Genotype codes:  0=AA; 1=AB
 * Phenotype codes: 0=missing; 1=AA; 2=AB
 *
 **********************************************************************/

double init_bc(int true_gen);

double emit_bc(int obs_gen, int true_gen, double error_prob);

double step_bc(int gen1, int gen2, double rf, double junk);

double nrec_bc(int gen1, int gen2);

void calc_genoprob_bc(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob);

void sim_geno_bc(int *n_ind, int *n_pos, int *n_draws, int *geno,
		 double *rf, double *error_prob, int *draws);

void est_map_bc(int *n_ind, int *n_mar, int *geno, double *rf, 
		double *error_prob, double *loglik, int *maxit, 
		double *tol, int *verbose);

void argmax_geno_bc(int *n_ind, int *n_pos, int *geno, 
		   double *rf, double *error_prob, int *argmax);

double errorlod_bc(int obs, double *prob, double error_prob);

void calc_errorlod_bc(int *n_ind, int *n_mar, int *geno, 
		      double *error_prob, double *genoprob, 
		      double *errlod);

void est_rf_bc(int *n_ind, int *n_mar, int *geno, double *rf);

void calc_pairprob_bc(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob,
		      double *pairprob);

/* end of hmm_bc.h */
