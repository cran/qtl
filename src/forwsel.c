/**********************************************************************
 * 
 * forwsel.c
 * 
 * copyright (c) 2007, Karl W Broman
 *
 * last modified Mar, 2007
 * first written Jan, 2007
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * This is a simple routine to do forward selection in regression
 * to a fixed number of covariates
 *
 * Contains: R_forwsel, forwsel
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include "forwsel.h"


/* R wrapper */
void R_forwsel(int *n, int *m, double *x, double *y,
	       int *maxsize, int *chosen, double *rss)
{
  double **X;
  int i;

  /* reorganize x matrix */
  X = (double **)R_alloc(*m, sizeof(double *));
  X[0] = x;
  for(i=1; i< *m; i++) X[i] = X[i-1] + *n;

  forwsel(*n, *m, X, y, *maxsize, chosen, rss);
}

/**********************************************************************
 * forwsel 
 * 
 * n = number of individuals
 * m = number of covariates (not including intercept)
 *
 * X = covariate matrix, indexed as X[covariate][individual]
 * y = outcome
 *
 * maxsize = maximum number of covariates
 *
 * chosen = on output, index [0, 1, ..., (m-1)] of chosen covariates
 * rss = on output, rss for those models
 *
 **********************************************************************/
void forwsel(int n, int m, double **X, double *y,
	     int maxsize, int *chosen, double *rss)
{
  double minrss, *work, sxx, sxy, syy; 
  double minsxx=0.0, minsxy=0.0, currss;
  int *ignore, ym;
  int i, j, k;

  /* allocate space */
  work = (double *)R_alloc(m, sizeof(double));
  ignore = (int *)R_alloc(m, sizeof(int));
  
  for(j=0; j<m; j++) {
    ignore[j] = 0;
    work[j] = 0.0;
  }

  /* recenter everything */
  ym = 0.0;
  for(i=0; i<n; i++) {
    ym += y[i];
    for(j=0; j<m; j++) work[j] += X[j][i];
  }

  ym /= (double)n;
  for(j=0; j<m; j++) work[j] /= (double)n;

  minrss = 0.0;
  for(i=0; i<n; i++) {
    y[i] -= ym;
    minrss += (y[i]*y[i]);
    for(j=0; j<m; j++) X[j][i] -= work[j];
  }

  for(k=0; k<maxsize; k++) { 
    chosen[k] = -1;

    syy = minrss;

    /* loop over markers */
    for(j=0; j<m; j++) {
      if(!ignore[j]) {

	/* calculate RSS */
	sxx = sxy = 0.0;
	for(i=0; i<n; i++) {
	  sxx += (X[j][i]*X[j][i]);
	  sxy += (X[j][i]*y[i]);
	}
	
	currss = syy - sxy*sxy/sxx;

	if(currss < minrss) { /* best so far */
	  rss[k] = minrss = currss;
	  minsxx = sxx;
	  minsxy = sxy;
	  chosen[k] = j;
	}
      }

    }

    ignore[chosen[k]] = 1;

    /* recenter y */
    for(i=0; i<n; i++) 
      y[i] -= (X[chosen[k]][i]*minsxy/minsxx);

    /* recenter other x's */
    for(j=0; j<m; j++) {
	
      if(!ignore[j]) {
      
	sxy = 0.0;
	for(i=0; i<n; i++) 
	  sxy += (X[j][i]*X[chosen[k]][i]);
	for(i=0; i<n; i++) 
	  X[j][i] -= (X[chosen[k]][i] * sxy/minsxx);
      }
    }

  }
}

/* end of forwsel.c */
