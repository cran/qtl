/**********************************************************************
 * 
 * simulate.c
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "util.h"
#include "simulate.h"

/* wrapper for sim_all_ril, to be called from R */
void R_sim_ril(int *n_chr, int *n_mar, int *n_ril, double *map,
	       int *n_str, int *m, int *include_x, int *random_cross,
	       int *cross, int *ril)
{
  GetRNGstate();

  sim_all_ril(*n_chr, n_mar, *n_ril, map, *n_str, *m, *include_x, 
	      *random_cross, cross, ril);

  PutRNGstate();
}
	  
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
		 int *cross, int *ril)
{
  int i, j, tot_mar;
  int *par1a, *par1b, *par2a, *par2b;
  int *kid1a, *kid1b, *kid2a, *kid2b;
  int **Par1a, **Par1b, **Par2a, **Par2b;
  int **Kid1a, **Kid1b, **Kid2a, **Kid2b;
  int **Ril, **Cross;
  double **Map;

  /* count total number of markers */
  for(i=0, tot_mar=0; i<n_chr; i++) 
    tot_mar += n_mar[i];

  reorg_geno(tot_mar, n_ril, ril, &Ril);
  reorg_geno(n_str, n_ril, cross, &Cross);

  /* allocate space */
  Map = (double **)R_alloc(n_chr, sizeof(double *));
  Map[0] = map;
  for(i=1; i<n_chr; i++)
    Map[i] = Map[i-1] + n_mar[i-1];

  allocate_int(tot_mar, &par1a);
  allocate_int(tot_mar, &par1b);
  allocate_int(tot_mar, &par2a);
  allocate_int(tot_mar, &par2b);
  allocate_int(tot_mar, &kid1a);
  allocate_int(tot_mar, &kid1b);
  allocate_int(tot_mar, &kid2a);
  allocate_int(tot_mar, &kid2b);
  Par1a = (int **)R_alloc(n_chr, sizeof(int *));
  Par1b = (int **)R_alloc(n_chr, sizeof(int *));
  Par2a = (int **)R_alloc(n_chr, sizeof(int *));
  Par2b = (int **)R_alloc(n_chr, sizeof(int *));
  Kid1a = (int **)R_alloc(n_chr, sizeof(int *));
  Kid1b = (int **)R_alloc(n_chr, sizeof(int *));
  Kid2a = (int **)R_alloc(n_chr, sizeof(int *));
  Kid2b = (int **)R_alloc(n_chr, sizeof(int *));
  Par1a[0] = par1a;
  Par1b[0] = par1b;
  Par2a[0] = par2a;
  Par2b[0] = par2b;
  Kid1a[0] = kid1a;
  Kid1b[0] = kid1b;
  Kid2a[0] = kid2a;
  Kid2b[0] = kid2b;
  for(i=1; i<n_chr; i++) {
    Par1a[i] = Par1a[i-1] + n_mar[i-1];
    Par1b[i] = Par1b[i-1] + n_mar[i-1];
    Par2a[i] = Par2a[i-1] + n_mar[i-1];
    Par2b[i] = Par2b[i-1] + n_mar[i-1];
    Kid1a[i] = Kid1a[i-1] + n_mar[i-1];
    Kid1b[i] = Kid1b[i-1] + n_mar[i-1];
    Kid2a[i] = Kid2a[i-1] + n_mar[i-1];
    Kid2b[i] = Kid2b[i-1] + n_mar[i-1];
  }

  for(i=0; i<n_ril; i++) 
    sim_ril(n_chr, n_mar, tot_mar, Map, n_str, m, Ril[i], include_x,
	    random_cross, Cross[i], 
	    Par1a, Par1b, Par2a, Par2b, Kid1a, Kid1b, Kid2a, Kid2b);
}

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
	     int **Kid1a, int **Kid1b, int **Kid2a, int **Kid2b)
{
  int i, j, k, s, flag;

  for(i=0; i<n_str; i++) cross[i] = i+1;
  if(random_cross) int_permute(cross, n_str);

  if(n_str==2) {
    for(i=0; i<n_chr; i++) {
      for(j=0; j<n_mar[i]; j++) {
	Par1a[i][j] = Par2a[i][j] = 1;
	Par1b[i][j] = Par2b[i][j] = 2;
      }
    }
  }
  else if(n_str==4) {
    for(i=0; i<n_chr; i++) {
      for(j=0; j<n_mar[i]; j++) {
	Par1a[i][j] = 1;
	Par1b[i][j] = 2;
	Par2a[i][j] = 3;
	Par2b[i][j] = 4;
      }
    }
  }
  else { 
    for(i=0; i<n_chr; i++) {
      for(j=0; j<n_mar[i]; j++) {
	Par1a[i][j] = 1;
	Par1b[i][j] = 2;
	Par2a[i][j] = 3;
	Par2b[i][j] = 4;
      }
      meiosis(n_mar[i], Par1a[i], Par1b[i], map[i], m, Kid1a[i]);
      if(include_x && i==n_chr-1) /* X chromosome */
	for(j=0; j<n_mar[i]; j++) Kid1b[i][j] = Par2a[i][j];
      else
	meiosis(n_mar[i], Par2a[i], Par2b[i], map[i], m, Kid1b[i]);

      for(j=0; j<n_mar[i]; j++) {
	Par1a[i][j] = 5;
	Par1b[i][j] = 6;
	Par2a[i][j] = 7;
	Par2b[i][j] = 8;
      }
      meiosis(n_mar[i], Par1a[i], Par1b[i], map[i], m, Kid2a[i]);
      if(include_x && i==n_chr-1) /* X chromosome */
	for(j=0; j<n_mar[i]; j++) Kid2b[i][j] = Kid2a[i][j];
      else
	meiosis(n_mar[i], Par2a[i], Par2b[i], map[i], m, Kid2b[i]);

      for(j=0; j<n_mar[i]; j++) {
	Par1a[i][j] = Kid1a[i][j];
	Par1b[i][j] = Kid1b[i][j];
	Par2a[i][j] = Kid2a[i][j];
	Par2b[i][j] = Kid2b[i][j];
      }
    }
  } 

  while(1) { /* now do inbreeding to fixation */
    for(i=0; i<n_chr; i++) {
      meiosis(n_mar[i], Par1a[i], Par1b[i], map[i], m, Kid1a[i]);
      if(include_x && i==n_chr-1) /* X chromosome */
	for(j=0; j<n_mar[i]; j++) Kid1b[i][j] = Par2a[i][j];
      else
	meiosis(n_mar[i], Par2a[i], Par2b[i], map[i], m, Kid1b[i]);

      meiosis(n_mar[i], Par1a[i], Par1b[i], map[i], m, Kid2a[i]);
      if(include_x && i==n_chr-1) /* X chromosome */
	for(j=0; j<n_mar[i]; j++) Kid2b[i][j] = Kid2a[i][j];
      else
	meiosis(n_mar[i], Par2a[i], Par2b[i], map[i], m, Kid2b[i]);
    }

    flag = 0;
    for(i=0; i<n_chr; i++) {
      for(j=0; j<n_mar[i]; j++) {
	if(!(Kid1a[i][j] == Kid1b[i][j] && Kid1a[i][j] == Kid2a[i][j] &&
	     Kid1a[i][j] == Kid2b[i][j])) {
	  flag=1;
	  break;
	}
      }
    }

    if(!flag) { /* we're done! */
      for(i=0, k=0; i<n_chr; i++) 
	for(j=0; j<n_mar[i]; j++, k++) 
	  ril[k] = cross[Kid1a[i][j]-1];
      return;
    }

    for(i=0; i<n_chr; i++) {
      for(j=0; j<n_mar[i]; j++) {
	Par1a[i][j] = Kid1a[i][j];
	Par1b[i][j] = Kid1b[i][j];
	Par2a[i][j] = Kid2a[i][j];
	Par2b[i][j] = Kid2b[i][j];
      }
    }
  }
  
}

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
void meiosis(int n_mar, int *chr1, int *chr2, double *map, int m, int *product)
{
  double L, *xo;
  int i, j, first, n, cur;

  L = map[n_mar-1];

  if(m > 0) { /* crossover interference */

    /* simulate number of XOs and intermediates */
    n = (int) rpois(L*(double)(m+1)/50.0);

    /* simulate locations */
    allocate_double(n, &xo);
    for(i=0; i<n; i++) 
      xo[i] = L*unif_rand();
    /* sort them */
    R_rsort(xo, n);

    /* which is the first crossover? */
    first = random_int(1,m+1);

    for(i=first, j=0; i<n; i += (m+1), j++) 
      xo[j] = xo[i];
    n = j;
  
    /* thin with probability 1/2 */
    for(i=0, j=0; i<n; i++) {
      if(unif_rand() < 0.5) {
	xo[j] = xo[i]; 
	j++;
      }
    }
    n = j;
  }
  else { /* no crossover interference */
    n = (int) rpois(L/100.0);

    /* simulate locations */
    allocate_double(n, &xo);
    for(i=0; i<n; i++) 
      xo[i] = L*unif_rand();
    /* sort them */
    R_rsort(xo, n);
  }

  /* determine which intervals recombined */
  for(i=1,cur=0; i<n_mar; i++) {
    product[i] = 0;

    if(cur < n && xo[cur] <= map[i]) {
      while(1) {
	product[i] = 1-product[i];
	cur++;
	if(cur >= n || xo[cur] > map[i]) 
	  break;
      }
    }
  }
  
  if(unif_rand() < 0.5) product[0] = 0;
  else product[0] = 1;
  for(i=1; i<n_mar; i++) {
    if(product[i] == 1)
      product[i] = 1-product[i-1];
    else
      product[i] = product[i-1];
  }
  for(i=0; i<n_mar; i++) {
    if(product[i]==0) product[i] = chr1[i];
    else product[i] = chr2[i];
  }

}


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
	    double error_prob, double missing_prob)
{
  int i, j, k, temp;

  for(i=0; i<n_ril; i++) {
    for(j=0; j<tot_mar; j++) {
      temp = Parents[Geno[j][i]-1][j];
      if(unif_rand() < error_prob)  /* switch the SNP genotype */
	temp = 1-temp;

      Geno[j][i] = 0;
      if(unif_rand() > missing_prob) {/* no error; convert to bit string */
	for(k=0; k<8; k++) 
	  if(temp == Parents[k][j]) Geno[j][i] += (1<<k);
      }
    }
  }
}

/* wrapper for calling sim_cc from R */
void R_sim_cc(int *n_ril, int *tot_mar, int *parents, int *geno,
	      double *error_prob, double *missing_prob)
{
  int **Parents, **Geno;

  reorg_geno(*tot_mar, 8, parents, &Parents);
  reorg_geno(*n_ril, *tot_mar, geno, &Geno);

  GetRNGstate();

  sim_cc(*n_ril, *tot_mar, Parents, Geno, *error_prob, *missing_prob);

  PutRNGstate();
}

/* end of simulate.c */

