/**********************************************************************
 * 
 * fitqtl_imp.h
 *
 * copyright (c) 2002, Hao Wu, The Jackson Laboratory
 *
 * last modified May, 2002
 * first written Apr, 2002
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a general genome scan by 
 * imputation. 
 *
 * Contains: R_fitqtl_imp, fitqtl_imp, nullRss0, galtRss
 *
 **********************************************************************/

void R_fitqtl_imp(int *n_ind, int *n_qtl, int *n_gen, int *n_draws,
		  int *draws, int *n_cov, double *cov, int *model, 
		  int *n_int, double *pheno,
		  /* return variables */
		  double *lod, int* df);


void fitqtl_imp(int n_ind, int n_qtl, int *n_gen, int n_draws, 
		int ***Draws, double **Cov, int n_cov, 
		int *model, int n_int, double *pheno, double *lod, int* df);


double nullRss0(double *pheno, int n_ind);

double galtRss(double *pheno, int n_ind, int *n_gen, int n_qtl, 
	       int **Draws, double **Cov, int n_cov, int *model, 
	       int n_int, double *dwork, int *iwork, int sizefull);

/* end of fitqtl_imp.h */
