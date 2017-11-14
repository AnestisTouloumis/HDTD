#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _HDTD_centerdatamatrix(SEXP, SEXP);
extern SEXP _HDTD_covmathat_statistics(SEXP, SEXP);
extern SEXP _HDTD_covmathat_statistics_centered(SEXP, SEXP);
extern SEXP _HDTD_covmathat_statistics_trans(SEXP, SEXP);
extern SEXP _HDTD_covmathat_statistics_trans_centered(SEXP, SEXP);
extern SEXP _HDTD_crossprod2cpp(SEXP, SEXP);
extern SEXP _HDTD_crossprodcpp(SEXP);
extern SEXP _HDTD_meanmatts_statistics(SEXP, SEXP);
extern SEXP _HDTD_pmat(SEXP);
extern SEXP _HDTD_projectionmatrix(SEXP);
extern SEXP _HDTD_sampleSigmaR(SEXP, SEXP);
extern SEXP _HDTD_statistics(SEXP, SEXP);
extern SEXP _HDTD_statistics_centered(SEXP, SEXP);
extern SEXP _HDTD_statistics_trans(SEXP, SEXP);
extern SEXP _HDTD_statistics_trans_centered(SEXP, SEXP);
extern SEXP _HDTD_sumdatamatrix(SEXP, SEXP);
extern SEXP _HDTD_tcrossprod2cpp(SEXP, SEXP);
extern SEXP _HDTD_tcrossprodcpp(SEXP);
extern SEXP _HDTD_transposedatamatrix(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_HDTD_centerdatamatrix",                    (DL_FUNC) &_HDTD_centerdatamatrix,                    2},
  {"_HDTD_covmathat_statistics",                (DL_FUNC) &_HDTD_covmathat_statistics,                2},
  {"_HDTD_covmathat_statistics_centered",       (DL_FUNC) &_HDTD_covmathat_statistics_centered,       2},
  {"_HDTD_covmathat_statistics_trans",          (DL_FUNC) &_HDTD_covmathat_statistics_trans,          2},
  {"_HDTD_covmathat_statistics_trans_centered", (DL_FUNC) &_HDTD_covmathat_statistics_trans_centered, 2},
  {"_HDTD_crossprod2cpp",                       (DL_FUNC) &_HDTD_crossprod2cpp,                       2},
  {"_HDTD_crossprodcpp",                        (DL_FUNC) &_HDTD_crossprodcpp,                        1},
  {"_HDTD_meanmatts_statistics",                (DL_FUNC) &_HDTD_meanmatts_statistics,                2},
  {"_HDTD_pmat",                                (DL_FUNC) &_HDTD_pmat,                                1},
  {"_HDTD_projectionmatrix",                    (DL_FUNC) &_HDTD_projectionmatrix,                    1},
  {"_HDTD_sampleSigmaR",                        (DL_FUNC) &_HDTD_sampleSigmaR,                        2},
  {"_HDTD_statistics",                          (DL_FUNC) &_HDTD_statistics,                          2},
  {"_HDTD_statistics_centered",                 (DL_FUNC) &_HDTD_statistics_centered,                 2},
  {"_HDTD_statistics_trans",                    (DL_FUNC) &_HDTD_statistics_trans,                    2},
  {"_HDTD_statistics_trans_centered",           (DL_FUNC) &_HDTD_statistics_trans_centered,           2},
  {"_HDTD_sumdatamatrix",                       (DL_FUNC) &_HDTD_sumdatamatrix,                       2},
  {"_HDTD_tcrossprod2cpp",                      (DL_FUNC) &_HDTD_tcrossprod2cpp,                      2},
  {"_HDTD_tcrossprodcpp",                       (DL_FUNC) &_HDTD_tcrossprodcpp,                       1},
  {"_HDTD_transposedatamatrix",                 (DL_FUNC) &_HDTD_transposedatamatrix,                 2},
  {NULL, NULL, 0}
};

void R_init_HDTD(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}