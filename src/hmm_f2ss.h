/**********************************************************************
 * 
 * hmm_f2ss.h
 * 
 * copyright (c) 2001, Karl W Broman, Johns Hopkins University
 *
 * last modified Sep, 2001
 * first written Sep, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 * 
 * C functions for the R/qtl package
 *
 * Contains: init_f2ss, emit_f2ss, step_f2ss, nrec_f2ss1,
 *           nrec_f2ss2, est_map_f2ss, argmax_geno_f2ss
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the F2 intercross with 
 * sex-specific maps.
 *
 * Genotype codes:  1=AA; 2=AB; 3=BA; 4=BB
 * Phenotype codes: 0=missing; 1=AA; 2=AB; 3=BB; 4=not BB; 5=not AA
 *
 **********************************************************************/

double init_f2ss(int true_gen);

double emit_f2ss(int obs_gen, int true_gen, double error_prob);
    
double step_f2ss(int gen1, int gen2, double rf1, double rf2);

double nrec_f2ss1(int gen1, int gen2);

double nrec_f2ss2(int gen1, int gen2);

void est_map_f2ss(int *n_ind, int *n_mar, int *geno, double *rf1,
		  double *rf2, double *error_prob, double *loglik,
		  int *maxit, double *tol, int *junk, int *trace);

void argmax_geno_f2ss(int *n_ind, int *n_pos, int *geno, 
		      double *rf1, double *rf2, 
		      double *error_prob, int *argmax);

/* end of hmm_f2ss.h */
