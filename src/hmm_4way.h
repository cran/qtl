/**********************************************************************
 * 
 * hmm_4way.h
 * 
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 * Aug, 2001; Feb, 2001
 * Licensed under the GNU General Public License version 2 (June, 1991)
 * 
 * C functions for the R/qtl package
 * 
 * Contains: init_4way, emit_4way, step_4way, nrec_4way, nrec_4way1,
 *           nrec_4way2, calc_genoprob_4way, sim_geno_4way, 
 *           est_map_4way, argmax_geno_4way, errorlod_4way, 
 *           calc_errorlod_4way, nrec2_4way, logprec_4way, est_rf_4way
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the "4-way" cross (autosomal data)
 *
 * Genotype codes:  0=AC, 1=AD, 2=BC, 3=BD
 * Phenotype codes: 0=missing, 1=AC, 2=BC, 3=AD, 4=BD,
 *                  5=(AC/AD=A), 6=(BC/BD=B),7=(AC/BC=C),8=(AD/BD=D)
 *                  9=(AC/BD), 10=(AD/BC)
 *
 **********************************************************************/

double init_4way(int true_gen);

double emit_4way(int obs_gen, int true_gen, double error_prob);

double step_4way(int gen1, int gen2, double rf1, double rf2);

double nrec_4way(int gen1, int gen2);

double nrec_4way1(int gen1, int gen2);

double nrec_4way2(int gen1, int gen2);

void calc_genoprob_4way(int *n_ind, int *n_mar, int *geno, 
			double *rf1, double *rf2, double *error_prob, 
			double *genoprob);

void sim_geno_4way(int *n_ind, int *n_pos, int *n_draws, int *geno,
		   double *rf1, double *rf2, double *error_prob, int *draws);

void est_map_4way(int *n_ind, int *n_mar, int *geno, double *rf1, double *rf2,
		  double *error_prob, double *loglik, int *maxit, 
		  double *tol, int *sexsp, int *prnt);

void argmax_geno_4way(int *n_ind, int *n_pos, int *geno, 
		      double *rf1, double *rf2, 
		      double *error_prob, int *argmax);

double errorlod_4way(int obs, double *prob, double error_prob);

void calc_errorlod_4way(int *n_ind, int *n_mar, int *geno, 
		      double *error_prob, double *genoprob, 
		      double *errlod);

double nrec2_4way(int obs1, int obs2, double rf);

double logprec_4way(int obs1, int obs2, double rf);

void est_rf_4way(int *n_ind, int *n_mar, int *geno, double *rf, 
	       int *maxit, double *tol);

/* end of hmm_4way.h */
