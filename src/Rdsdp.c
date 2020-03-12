/*
 * R interface to DSDP semidefinite programming library
 *
 * Created: Zhisu Zhu <zhuzhisu@gmail.com>
 * December 08, 2014
 */

#include <R.h>
#include <Rinternals.h>
#include <dsdp5.h>
#include <R_ext/Rdynload.h>

SEXP int_vector_dsdp2R(int, int*);
SEXP double_vector_dsdp2R(int, double*);

int rReadSDPAFile(char*, char*, double**, int*, double**, int*, double**, int*);

SEXP dsdp(SEXP data_filename_p, SEXP options_filename_p)
{
  SEXP X_p, y_p, ret;
  SEXP STATS_p;

  double *y, *X, *STATS; 
  int nvars, nXvars, nstats;  
  char* data_filename;
  char* options_filename;

  data_filename = (char*) CHAR(STRING_ELT(data_filename_p,0));
  options_filename = (char*) CHAR(STRING_ELT(options_filename_p,0));

  int status, info;


  /*
   * Solve the problem
   */
  status = rReadSDPAFile(data_filename, options_filename, &y, &nvars, &X, &nXvars, &STATS, &nstats);
  if (status) Rprintf("Error: reading sdpa file %s, status: %d.\n",data_filename,status);

  /* 
   * Get the results
   */
  
  X_p = PROTECT(double_vector_dsdp2R(nXvars,X));
  DSDPFREE(&X, &info);
  
  
  /* Copy y */
  y_p = PROTECT(double_vector_dsdp2R(nvars, y));
  DSDPFREE(&y, &info);

  /* Copy STATS */
  STATS_p = PROTECT(double_vector_dsdp2R(nstats, STATS));
  DSDPFREE(&STATS, &info);

  ret = PROTECT(allocVector(VECSXP,3));
  SET_VECTOR_ELT(ret,0,X_p);
  SET_VECTOR_ELT(ret,1,y_p);
  SET_VECTOR_ELT(ret,2,STATS_p);

  UNPROTECT(4);
  return ret;
}

SEXP int_vector_dsdp2R(int n, int *y)
{
  SEXP ret;
  int i;
  int *intvec;

  ret = allocVector(INTSXP,n);
  intvec = INTEGER(ret);
  for (i=0; i<n; i++)
    intvec[i] = y[i];
  return ret;
}

SEXP double_vector_dsdp2R(int n, double *y)
{
  SEXP ret;
  int i;
  double *dblvec;

  ret = allocVector(REALSXP,n);
  dblvec = REAL(ret);
  for (i=0; i<n; i++)
    dblvec[i] = y[i];
  return ret;
}

// RegisteringDynamic Symbols
void R_init_Rdsdp(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

