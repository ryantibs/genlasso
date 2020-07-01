#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern void C_givens(double a, double b, double *c, double *s);
extern void C_rowrot(double *A, int i1, int i2, int m, int n, int j1, int j2, double c, double s);
extern void C_colrot(double *A, int j1, int j2, int m, int n, int i1, int i2, double c, double s);
extern void C_downdate1(double *Q1, double *R, int *j0p, int *mp, int *np);
extern void C_update1(double *Q2, double *w, int *mp, int *kp);
extern void C_downdate2(double *Q, double *R, int *mp, int *np);
extern void C_update2(double *y, double *D, double *r, int *mp, int *np, int *qp);
extern void C_maketri1(double *y, double *A, double *R, int *mp, int *np, int *kp);
extern void C_maketri2(double *y, double *A, double *R, int *mp, int *np, int *kp);
extern void C_maketri3(double *y, double *A, double *R, int *m1p, int *m2p, int *np, int *qp);
extern void C_maketri4(double *y, double *A, double *Q, double *R, int *m1p, int *m2p, int *np, int *qp, int *kp);

static const R_CMethodDef R_CDef[] = {
   CALLDEF(C_colrot, 9),
   CALLDEF(C_downdate1, 5),
   CALLDEF(C_downdate2, 4),
   CALLDEF(C_givens, 4),
   CALLDEF(C_maketri1, 6),
   CALLDEF(C_maketri2, 6),
   CALLDEF(C_maketri3, 7),
   CALLDEF(C_maketri4, 9),
   CALLDEF(C_rowrot, 9),
   CALLDEF(C_update1, 4),
   CALLDEF(C_update2, 6),
   {NULL, NULL, 0}
};

void R_init_genlasso(DllInfo *dll)
{
    R_registerRoutines(dll, R_CDef, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
