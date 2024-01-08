/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include <math.h>
#include <string.h>

#include "lib_poisson1D.h"

void eig_poisson1D(double *eigval, int *la)
{
  memset(eigval, 0, sizeof(double) * *la);
  for (int k = 1; k <= *la; ++k)
  {
    eigval[k - 1] = 2 - (2 * cos(k * M_PI / (*la + 1)));
  }
}

double eigmax_poisson1D(int *la)
{
  return (2 - (2 * cos(*la * M_PI / (*la + 1))));
}

double eigmin_poisson1D(int *la){
  return (2 - (2 * cos(1 * M_PI / (*la + 1))));
}

double richardson_alpha_opt(int *la){
  return (2 / (eigmax_poisson1D(la) + eigmin_poisson1D(la)));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  // RHS : right hand side
  // tol : tolérance, epsilon
  // maxit : nombre maximal d'itération
  // resvec : stocke le résidu
  // nbites : nombre d'itération
  
  /*

  r = norm(b - Ax) / norm(b)
  -> calcule norme de b
  -> dgbmv sur b-Ax puis norme de nouveau b 
  RESVEC = RSH - AB * X
  while RESVEC[k+1] > tol && maxit > 0 
  Xk+1] = X[k] + alpha_rich(RHS - AB * X[k]) 
  RESVEC
  
  --maxit
  nbr itération ++
  */

  double* Y = (double *) calloc(*la, sizeof(double));
  const double norm_b = (1 / cblas_dnrm2(*la, RHS, 1)); // division de norme de b

  // copy de b dans y
  cblas_dcopy(*la, RHS, 1.0, Y, 1.0);
  // y = y - Ax
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1.0, 1.0, Y,  1.0);

  // calcul residu
  resvec[*nbite] = cblas_dnrm2(*la, Y, 1) * norm_b;
  while (resvec[*nbite] > *tol && *maxit > *nbite)
  {
    // calcule de x(k+1) -> X = X + alphaY
    cblas_daxpy(*la, *alpha_rich, Y, 1.0, X, 1.0);

    // résidu
    cblas_dcopy(*la, RHS, 1.0, Y, 1.0);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1.0, 1.0, Y, 1.0);
    ++*nbite;
    resvec[*nbite] = cblas_dnrm2(*la, Y, 1) * norm_b;
  }
  free(Y);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv)
{
  // M = D 
  memset(MB, 0, sizeof(double) * *la * *lab);
  for (int indice = *kv + 1; indice < (*lab) * (*la); indice += *lab)
  {
    MB[indice] = AB[indice];
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv)
{
  // M = D - E
  memset(MB, 0, sizeof(double) * *la * *lab);
  for (int indice = *kv + 1; indice < (*lab) * (*la); indice += *lab)
  {
    MB[indice] = AB[indice];
    MB[indice + 1] = AB[indice + 1];
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  double norm_b = (1.0 / cblas_dnrm2(*la, RHS, 1));
  int NRHS = 1;
  int ku_MB = *ku - 1; // jsp pk faut le mettre
  int *ipiv = (int *)calloc(*la, sizeof(int));
  int info = 0;

  // Calcul du résidu initial
  double *Y = (double *)calloc(*la, sizeof(double));
  cblas_dcopy(*la, RHS, 1, Y, 1);
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, Y, 1);
  resvec[*nbite] = cblas_dnrm2(*la, Y, 1) * norm_b;

  // Résolution itérative
  dgbtrf_(la, la, kl, &ku_MB, MB, lab, ipiv, &info);
  while (resvec[*nbite] > *tol && *nbite < *maxit)
  {
    // factorisation LU + résolution
    dgbtrs_("N", la, kl, &ku_MB, &NRHS, MB, lab, ipiv, Y, la, &info, 1);

    // Mise à jour de X : (X(k+1) -> X = X + Y)
    cblas_daxpy(*la, 1.0, Y, 1, X, 1);

    // Calcul du résidu
    cblas_dcopy(*la, RHS, 1, Y, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, Y, 1);
    ++*nbite;
    resvec[*nbite] = cblas_dnrm2(*la, Y, 1) * norm_b;
  }

  free(Y);
  free(ipiv);
}