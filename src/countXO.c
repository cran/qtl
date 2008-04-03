/**********************************************************************
 * 
 * countXO.c
 *
 * copyright (c) 2008, Karl W Broman
 *
 * last modified Feb, 2008
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for comparing marker orders by counts of 
 * obligate crossovers
 *
 * Contains: R_countXO_bc, R_countXO_f2, R_countXO_4way, countXO
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "ripple.h"
#include "countXO.h"

/**********************************************************************
 * 
 * countXO
 *
 * This function counts the number of obligate crossovers for each 
 * individual on a chromosome
 *
 * Input:
 *
 * n_ind    = no. individuals
 *
 * n_mar    = no. markers
 *
 * n_gen    = no. possible genotypes
 *
 * geno     = genotype data [n_ind x n_mar]
 *
 * nxo      = the output; the number of obligate crossovers for each
 *            individual
 *
 * countxo  = function to count the number of obligate crossovers in
 *            an interval and to update the current inferred genotype
 *            (specific for backcross, intercross, and four-way cross)
 *
 **********************************************************************/

void countXO(int n_ind, int n_mar, int n_gen, int *geno,
	     int *nxo, int countxo(int *curgen, int nextgen))
{
  int **Geno;
  int j, k, curgen;

  /* reorganize genotype data and marker order matrix */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  
  for(j=0; j<n_ind; j++) { /* loop over individuals */

    R_CheckUserInterrupt(); /* check for ^C */
    nxo[j] = 0;

    /* genotype at first marker */
    curgen = Geno[0][j];
    for(k=1; k<n_mar; k++) /* loop over markers */
      /* count no obligate crossovers and update current genotype */
      nxo[j] += countxo(&curgen, Geno[k][j]);
  }
}


/**********************************************************************
 * 
 * R_countXO_bc
 *
 * Wrapper for call from R for a backcross
 * 
 **********************************************************************/

void R_countXO_bc(int *n_ind, int *n_mar, int *geno, 
		  int *nxo)
{
  countXO(*n_ind, *n_mar, 2, geno, nxo, countxo_bc);
}




/**********************************************************************
 * 
 * R_countXO_f2
 *
 * Wrapper for call from R for an intercross
 * 
 **********************************************************************/

void R_countXO_f2(int *n_ind, int *n_mar, int *geno, 
		  int *nxo)
{
  countXO(*n_ind, *n_mar, 4, geno, nxo, countxo_f2);
}


/**********************************************************************
 * 
 * R_countXO_4way
 *
 * Wrapper for call from R for a four-way cross
 * 
 **********************************************************************/

void R_countXO_4way(int *n_ind, int *n_mar, int *geno, 
		    int *nxo)
{
  countXO(*n_ind, *n_mar, 4, geno, nxo, countxo_4way);
}


/* end of countXO.c */

