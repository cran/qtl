/**********************************************************************
 * 
 * hmm_f2.h
 * 
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 * Aug, 2001; Feb, 2001
 * Licensed under the GNU General Public License version 2 (June, 1991)
 * 
 * C functions for the R/qtl package
 *
 * Contains: init_f2, emit_f2, step_f2, init_f2b, emit_f2b, step_f2b,
 *           calc_genoprob_f2, sim_genoprob_f2, est_map_f2, 
 *           argmax_geno_f2, errorlod_f2, calc_errorlod_f2, nrec2_f2,
 *           logprec_f2, est_rf_f2
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the F2 intercross.
 *
 * Genotype codes:  0=AA; 1=AB; 2=BB
 * Phenotype codes: 0=missing; 1=AA; 2=AB; 3=BB; 4=not BB; 5=not AA
 *
 **********************************************************************/

double init_f2(int true_gen);

double emit_f2(int obs_gen, int true_gen, double error_prob);
  
double step_f2(int gen1, int gen2, double rf, double junk);

double init_f2b(int true_gen);

double emit_f2b(int obs_gen, int true_gen, double error_prob);
  
double step_f2b(int gen1, int gen2, double rf, double junk);

double nrec_f2b(int gen1, int gen2);

void calc_genoprob_f2(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob);
  
void sim_geno_f2(int *n_ind, int *n_pos, int *n_draws, int *geno,
		 double *rf, double *error_prob, int *draws);

void est_map_f2(int *n_ind, int *n_mar, int *geno, double *rf, 
		double *error_prob, double *loglik, int *maxit, 
		double *tol, int *prnt);

void argmax_geno_f2(int *n_ind, int *n_pos, int *geno, 
		   double *rf, double *error_prob, int *argmax);

double errorlod_f2(int obs, double *prob, double error_prob);

void calc_errorlod_f2(int *n_ind, int *n_mar, int *geno, 
		      double *error_prob, double *genoprob, 
		      double *errlod);

double nrec2_f2(int obs1, int obs2, double rf);

double logprec_f2(int obs1, int obs2, double rf);

void est_rf_f2(int *n_ind, int *n_mar, int *geno, double *rf, 
	       int *maxit, double *tol);

/* end of hmm_f2.h */
