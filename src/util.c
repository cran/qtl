/**********************************************************************
 * 
 * util.c
 *
 * copyright (c) 2001-4, Karl W Broman, Johns Hopkins University
 *                       and Hao Wu, The Jackson Laboratory
 *
 * This file written mostly by Karl Broman with some additions
 * from Hao Wu.
 *
 * last modified Aug, 2004
 * first written Feb, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These are utility functions, mostly for the HMM engine.
 *
 * Other functions: addlog, subtrlog, reorg_geno, reorg_genoprob,
 *                  reorg_pairprob, allocate_int
 *                  allocate_alpha, reorg_draws, allocate_double,
 *                  sample_int, allocate_imatrix, allocate_dmatrix
 *                  reorg_errlod, double_permute, random_int
 *                  wtaverage
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "util.h"

#define THRESH 200.0

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

double addlog(double a, double b)
{
  if(b > a + THRESH) return(b);
  else if(a > b + THRESH) return(a);
  else return(a + log1p(exp(b-a)));
}
		       
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

double subtrlog(double a, double b)
{
  if(a > b + THRESH) return(a);
  else return(a + log1p(-exp(b-a)));
}

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

void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno)
{
  int i;

  *Geno = (int **)R_alloc(n_pos, sizeof(int *));

  (*Geno)[0] = geno;
  for(i=1; i< n_pos; i++) 
    (*Geno)[i] = (*Geno)[i-1] + n_ind;

}

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
		    double *genoprob, double ****Genoprob)
{
  int i, j;
  double **a;

  *Genoprob = (double ***)R_alloc(n_gen, sizeof(double **));

  a = (double **)R_alloc(n_pos*n_gen, sizeof(double *));

  (*Genoprob)[0] = a;
  for(i=1; i< n_gen; i++) 
    (*Genoprob)[i] = (*Genoprob)[i-1]+n_pos;
  
  for(i=0; i<n_gen; i++) 
    for(j=0; j<n_pos; j++) 
      (*Genoprob)[i][j] = genoprob + i*n_ind*n_pos + j*n_ind;
}

/**********************************************************************
 * 
 * reorg_pairprob
 *
 * Reorganize the joint genotype probabilities so that they form a 
 * quintuply indexed array rather than a single long vector
 *
 * Afterwards, pairprob indexed like 
 *    Pairprob[gen1][gen2][pos1][pos2][ind] with pos2 > pos1
 * 
 * You *must* refer to cases with pos2 > pos1, as cases with
 * pos2 <= pos1 point off into the ether.
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void reorg_pairprob(int n_ind, int n_pos, int n_gen, 
		    double *pairprob, double ******Pairprob)
{
  int i, j, k, s, n_pairs;
  double ****ptr1, ***ptr2, **ptr3;

  /* note: n_pos must be at least 2 */
  n_pairs = n_pos*(n_pos-1)/2;

  *Pairprob = (double *****)R_alloc(n_gen, sizeof(double ****));

  ptr1 = (double ****)R_alloc(n_gen*n_gen, sizeof(double ***));
  (*Pairprob)[0] = ptr1;
  for(i=1; i<n_gen; i++) 
    (*Pairprob)[i] = ptr1 + i*n_gen;

  ptr2 = (double ***)R_alloc(n_gen*n_gen*n_pos, sizeof(double **));
  for(i=0; i<n_gen; i++)
    for(j=0; j<n_gen; j++) 
      (*Pairprob)[i][j] = ptr2 + (i*n_gen+j)*n_pos;

  ptr3 = (double **)R_alloc(n_gen*n_gen*n_pos*n_pos, sizeof(double **));
  for(i=0; i<n_gen; i++) 
    for(j=0; j<n_gen; j++)
      for(k=0; k<n_pos; k++)
	(*Pairprob)[i][j][k] = ptr3 + ((i*n_gen+j)*n_pos + k)*n_pos;
  
  for(i=0; i<n_gen; i++) 
    for(j=0; j<n_gen; j++)
      for(k=0; k<n_pos-1; k++)
	for(s=(k+1); s<n_pos; s++) 
	  (*Pairprob)[i][j][k][s] = pairprob + (i*n_gen+j)*n_ind*n_pairs +
	    + n_ind*(2*n_pos-1-k)*k/2 + (s-k-1)*n_ind;
}

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

void allocate_alpha(int n_pos, int n_gen, double ***alpha)
{
  int i;

  *alpha = (double **)R_alloc(n_gen, sizeof(double *));

  (*alpha)[0] = (double *)R_alloc(n_gen*n_pos, sizeof(double));

  for(i=1; i< n_gen; i++) 
    (*alpha)[i] = (*alpha)[i-1] + n_pos;
}

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
		 int *draws, int ****Draws)
{
  int i, j;
  int **a;

  *Draws = (int ***)R_alloc(n_draws, sizeof(int **));

  a = (int **)R_alloc(n_pos*n_draws, sizeof(int *));
  (*Draws)[0] = a;
  for(i=1; i<n_draws; i++) 
    (*Draws)[i] = (*Draws)[i-1]+n_pos;
  
  for(i=0; i<n_draws; i++) 
    for(j=0; j<n_pos; j++) 
      (*Draws)[i][j] = draws + (i*n_pos+j)*n_ind;
}

/**********************************************************************
 * 
 * allocate_double
 *
 * Allocate space for a vector of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void allocate_double(int n, double **vector)
{
  *vector = (double *)R_alloc(n, sizeof(double));
}

/**********************************************************************
 * 
 * allocate_int
 *
 * Allocate space for a vector of ints
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void allocate_int(int n, int **vector)
{
  *vector = (int *)R_alloc(n, sizeof(int));
}

/**********************************************************************
 * 
 * allocate_dmatrix
 *
 * Allocate space for a matrix of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void allocate_dmatrix(int n_row, int n_col, double ***matrix)
{
  int i;

  *matrix = (double **)R_alloc(n_row, sizeof(double *));
  
  (*matrix)[0] = (double *)R_alloc(n_col*n_row, sizeof(double));

  for(i=1; i<n_row; i++) 
    (*matrix)[i] = (*matrix)[i-1]+n_col;
}

/**********************************************************************
 * 
 * allocate_imatrix
 *
 * Allocate space for a matrix of ints
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void allocate_imatrix(int n_row, int n_col, int ***matrix)
{
  int i;

  *matrix = (int **)R_alloc(n_row, sizeof(int *));
  
  (*matrix)[0] = (int *)R_alloc(n_col*n_row, sizeof(int));

  for(i=1; i<n_row; i++) 
    (*matrix)[i] = (*matrix)[i-1]+n_col;
}

/**********************************************************************
 * 
 * sample_int
 *
 * Make a single draw from (1, ..., n) with probs (p_0, ..., p_(n-1))
 *
 **********************************************************************/

int sample_int(int n, double *p)
{
  int i;
  double r;

  /* R's random number generator */
  r = unif_rand();

  for(i=0; i<n; i++) {
    if(r < p[i]) return(i+1);
    else r -= p[i];
  }
  return(n); /* this shouldn't happen */
}

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

void reorg_errlod(int n_ind, int n_mar, double *errlod, double ***Errlod)
{
  int i;

  *Errlod = (double **)R_alloc(n_mar, sizeof(double *));

  (*Errlod)[0] = errlod;
  for(i=1; i< n_mar; i++) 
    (*Errlod)[i] = (*Errlod)[i-1] + n_ind;
}

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

void double_permute(double *array, long len)
{
  long i, which;
  double tmp;
  
  for(i=0; i < len; i++) {
    which = random_int(i, len-1);
    tmp = array[which];
    array[which] = array[i];
    array[i] = tmp;
  }
}

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

int random_int(int low, int high)
{
  return((int)(unif_rand()*(double)(high - low + 1)) + low);
}

/**********************************************************************
 * wtaverage
 * calculate the weight average of the LOD scores
 *********************************************************************/

double wtaverage(double *LOD, int n_draws)
{
  int k, idx, nnewLOD;
  double sum, sums, meanLOD, varLOD, *newLOD;

  /* calculate the number of LOD needs to be thrown */
  idx = (int) floor( 0.5*log(n_draws)/log(2) );
  nnewLOD = n_draws-2*idx; /* number of items in newLOD vector */
  /* allocate memory for newLOD */  
  newLOD = (double *)R_alloc( nnewLOD, sizeof(double) );

  /* sort the LOD scores in ascending order */
  R_rsort(LOD, n_draws);

  /* get a new list of LOD scores, throwing the biggest 
     and smallest idx LOD scores. */
  for(k=idx, sum=0.0; k<n_draws-idx; k++) {
    newLOD[k-idx] = LOD[k];
    sum += LOD[k]; /* calculate the sum of newLOD in the same loop */
  }

  /* calculate the mean of newLOD */
  meanLOD = sum / nnewLOD; 
  /* calculate the variance of newLOD */
  for(k=0,sums=0.0; k<nnewLOD; k++) 
    sums += (newLOD[k]-meanLOD) * (newLOD[k]-meanLOD);
  varLOD = sums/(nnewLOD-1);

  /* return the weight average */
  return( meanLOD+0.5*log(10.0)*varLOD );

}

/* end of util.c */
